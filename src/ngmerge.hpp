#pragma once

#include "pear.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

struct NgmergeConfig {
    int    min_overlap  = 20;    // -m
    double max_mismatch = 0.10;  // -p
    bool   prefer_short = false; // -s
    int    phred_base   = 33;    // -q
    int    max_qual     = 40;    // -u
};

class NgmergeMerger {
public:
    explicit NgmergeMerger(const NgmergeConfig& cfg = NgmergeConfig())
        : cfg_(cfg) {}

    // Loads qual_profile.txt: two (max_qual + 1) square matrices under
    // "#match" and "#mismatch" headers. Returns false and leaves the
    // fastq-join rule in place on any problem, so a missing file degrades
    // rather than aborts.
    bool load_qual_profile(const std::string& path) {
        std::ifstream fin(path);

        if (!fin) {
            fprintf(stderr,
                    "Warning: cannot open NGmerge quality profile: %s\n",
                    path.c_str());
            return false;
        }

        const int n = cfg_.max_qual + 1;

        std::vector<std::vector<int>> cur, match, mism;
        int which = -1; // 0 = match, 1 = mismatch

        std::string line;

        auto flush = [&]() {
            if (which == 0)      match = cur;
            else if (which == 1) mism  = cur;
            cur.clear();
        };

        while (std::getline(fin, line)) {
            while (!line.empty() &&
                   (line.back() == '\n' || line.back() == '\r'))
                line.pop_back();

            if (line.empty())
                continue;

            if (line[0] == '#') {
                flush();

                std::string tag = line.substr(1);
                for (char& c : tag) c = (char)tolower((unsigned char)c);

                if (tag.rfind("mismatch", 0) == 0)   which = 1;
                else if (tag.rfind("match", 0) == 0) which = 0;
                else                                 which = -1;

                continue;
            }

            if (which < 0)
                continue;

            std::vector<int> row;
            std::string tok;
            std::stringstream ss(line);

            while (std::getline(ss, tok, ',')) {
                std::stringstream s2(tok);
                std::string t2;
                while (s2 >> t2)
                    row.push_back(atoi(t2.c_str()));
            }

            if (!row.empty())
                cur.push_back(std::move(row));
        }

        flush();

        // Some copies carry a row of column labels and a leading label on
        // each row; drop them when that is what makes the matrix one too
        // large in each direction.
        auto trim = [&](std::vector<std::vector<int>>& m) {
            if ((int)m.size() == n + 1)
                m.erase(m.begin());

            for (auto& r : m)
                if ((int)r.size() == n + 1)
                    r.erase(r.begin());
        };

        trim(match);
        trim(mism);

        auto square = [&](const std::vector<std::vector<int>>& m) {
            if ((int)m.size() != n)
                return false;
            for (const auto& r : m)
                if ((int)r.size() != n)
                    return false;
            return true;
        };

        if (!square(match) || !square(mism)) {
            fprintf(stderr,
                    "Warning: NGmerge quality profile %s is not two %dx%d\n"
                    "         matrices; the fastq-join rule stays in use.\n",
                    path.c_str(), n, n);
            return false;
        }

        match_ = std::move(match);
        mism_ = std::move(mism);
        have_profile_ = true;
        return true;
    }

    PearMergeResult merge(
        const PearFastqRecord& r1,
        const PearFastqRecord& r2) const
    {
        PearMergeResult result;

        if (r1.seq.empty() || r2.seq.empty() ||
            r1.seq.size() != r1.qual.size() ||
            r2.seq.size() != r2.qual.size())
        {
            return result;
        }

        PearFastqRecord rc;
        pear_reverse_complement(r2, rc);

        const int n1 = (int)r1.seq.size();
        const int n2 = (int)rc.seq.size();

        std::vector<int> q1(n1), q2(n2);

        for (int i = 0; i < n1; ++i)
            q1[i] = qual_of(r1.qual[i]);

        for (int i = 0; i < n2; ++i)
            q2[i] = qual_of(rc.qual[i]);

        // findPos(): offsets are walked downward, and the best is only
        // replaced on a strictly smaller fraction, so an equal fraction
        // keeps the earlier, longer alignment unless -s asked otherwise.
        const int fail_pos = n1 - cfg_.min_overlap + 1;

        int pos = fail_pos;
        double best = NOTMATCH;

        for (int i = n1 - cfg_.min_overlap; i > -1 && n1 - i <= n2; --i) {
            const double res = compare(r1.seq, rc.seq, i, n1 - i);

            if (res < best || (res == best && cfg_.prefer_short)) {
                best = res;
                pos = i;
            }

            if (res == 0.0 && !cfg_.prefer_short)
                break;                 // exact match; nothing can beat it
        }

        if (pos == fail_pos)
            return result;             // stitch failure

        // createSeq()
        const int len = n2 + pos;

        std::string seq = r1.seq;
        std::string qual = r1.qual;

        seq.resize(len);
        qual.resize(len);

        for (int i = pos; i < len; ++i) {
            const int k = i - pos;

            if (i >= n1) {             // 2nd read only
                seq[i] = rc.seq[k];
                qual[i] = rc.qual[k];
                continue;
            }

            const char a = r1.seq[i];
            const char b = rc.seq[k];

            if (b == 'N')              // 2nd read uncalled: keep R1
                continue;

            if (a == 'N') {            // 1st read uncalled: take R2
                seq[i] = b;
                qual[i] = rc.qual[k];
                continue;
            }

            const int qa = q1[i];
            const int qb = q2[k];

            if (a != b) {
                // Equal quality is settled by position: the first half of
                // the stitched read belongs to R1's 5' end, the second to
                // R2's.
                const bool take_r2 =
                    (qa < qb) || (qa == qb && i >= len / 2.0);

                seq[i] = take_r2 ? b : a;

                qual[i] = to_ascii(mismatch_score(qa, qb));

                result.mismatches.push_back(
                    PearAlternativeBase{
                        i,
                        take_r2 ? b : a,
                        take_r2 ? a : b,
                        take_r2 ? qb : qa,
                        take_r2 ? qa : qb});
            }
            else {
                qual[i] = to_ascii(match_score(qa, qb));
            }
        }

        result.assembled = true;
        result.seq = std::move(seq);
        result.qual = std::move(qual);
        result.overlap = std::min(n1 - pos, n2);
        result.r2_start = pos;

        return result;
    }

private:
    static constexpr double NOTMATCH = 1.0;

    // compare(): the mismatch fraction of one candidate alignment, or
    // NOTMATCH once it cannot stay under the allowance. An N shortens the
    // divisor, and the alignment is abandoned when that takes the length
    // below the minimum overlap.
    double compare(
        const std::string& s1,
        const std::string& s2,
        int offset,
        int length) const
    {
        int mis = 0;
        int len = length;
        double allow = len * cfg_.max_mismatch;

        for (int i = 0; i < length; ++i) {
            const char a = s1[offset + i];
            const char b = s2[i];

            if (a == 'N' || b == 'N') {
                if (--len < cfg_.min_overlap ||
                    mis > len * cfg_.max_mismatch)
                    return NOTMATCH;

                allow = len * cfg_.max_mismatch;
            }
            else if (a != b && ++mis > allow) {
                return NOTMATCH;
            }
        }

        return (len > 0) ? (double)mis / (double)len : NOTMATCH;
    }

    int qual_of(char c) const {
        int q = (int)(unsigned char)c - cfg_.phred_base;
        return std::max(0, std::min(q, cfg_.max_qual));
    }

    char to_ascii(int q) const {
        q = std::max(0, std::min(q, cfg_.max_qual));
        return (char)(q + cfg_.phred_base);
    }

    // Both profiles are indexed [R1's score][R2's score], as createSeq()
    // does; which base won does not enter into it.
    int match_score(int qa, int qb) const {
        if (have_profile_)
            return match_[qa][qb];

        return std::max(qa, qb);          // fastq-join (-y)
    }

    int mismatch_score(int qa, int qb) const {
        if (have_profile_)
            return mism_[qa][qb];

        return std::abs(qa - qb);         // fastq-join (-y)
    }

    NgmergeConfig cfg_;

    bool have_profile_ = false;
    std::vector<std::vector<int>> match_;
    std::vector<std::vector<int>> mism_;
};