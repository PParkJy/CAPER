#pragma once

#include "pear.hpp"

#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>

struct FlashConfig {
    int    min_overlap   = 10;    // -m
    int    max_overlap   = 65;    // -M
    double max_mismatch  = 0.25;  // -x
    bool   cap_mismatch  = false; // --cap-mismatch-quals
    int    phred_base    = 33;    // -p
    int    max_qual      = 41;
};

class FlashMerger {
public:
    explicit FlashMerger(const FlashConfig& cfg = FlashConfig())
        : cfg_(cfg) {}

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

        // pair_align(): the second read may not be contained entirely in
        // the first, and the overlap may not fall below min_overlap.
        const int start = std::max(0, n1 - n2);
        const int end   = n1 - cfg_.min_overlap + 1;

        double best_density = cfg_.max_mismatch + 1.0;
        double best_qual_score = 0.0;
        int best_position = -1;
        int best_overlap = 0;

        for (int i = start; i < end; ++i) {
            int overlap_len = n1 - i;

            int mismatches = 0;
            int uncalled = 0;
            long qual_total = 0;

            for (int k = 0; k < overlap_len; ++k) {
                const char a = r1.seq[i + k];
                const char b = rc.seq[k];

                if (a == 'N' || b == 'N') {
                    ++uncalled;
                    continue;
                }

                if (a != b) {
                    ++mismatches;
                    qual_total += std::min(q1[i + k], q2[k]);
                }
            }

            overlap_len -= uncalled;

            if (overlap_len < cfg_.min_overlap)
                continue;

            // The divisor is capped; the counts above are not.
            const double score_len =
                (double)std::min(overlap_len, cfg_.max_overlap);

            const double density = (double)mismatches / score_len;
            const double qual_score = (double)qual_total / score_len;

            if (density <= best_density &&
                (density < best_density || qual_score < best_qual_score))
            {
                best_density = density;
                best_qual_score = qual_score;
                best_position = i;
                best_overlap = overlap_len;
            }
        }

        if (best_position < 0 || best_density > cfg_.max_mismatch)
            return result;

        // generate_combined_read()
        const int pos = best_position;
        const int overlap_len = n1 - pos;
        const int remaining = n2 - overlap_len;
        const int len = n1 + remaining;

        std::string seq;
        std::string qual;

        seq.reserve(len);
        qual.reserve(len);

        for (int i = 0; i < pos; ++i) {           // read 1 before the overlap
            seq.push_back(r1.seq[i]);
            qual.push_back(r1.qual[i]);
        }

        for (int k = 0; k < overlap_len; ++k) {
            const int i = pos + k;

            const char a = r1.seq[i];
            const char b = rc.seq[k];
            const int qa = q1[i];
            const int qb = q2[k];

            if (a == b) {
                seq.push_back(a);
                qual.push_back(to_ascii(std::max(qa, qb)));
                continue;
            }

            qual.push_back(
                to_ascii(cfg_.cap_mismatch
                             ? std::min(std::min(qa, qb), 2)
                             : std::max(std::abs(qa - qb), 2)));

            char sel, alt;
            int qs, qal;

            if (qa > qb) {
                sel = a; alt = b; qs = qa; qal = qb;
            }
            else if (qa < qb) {
                sel = b; alt = a; qs = qb; qal = qa;
            }
            else {
                // Equal quality: read 2's base, unless it is uncalled.
                if (b == 'N') { sel = a; alt = b; qs = qa; qal = qb; }
                else          { sel = b; alt = a; qs = qb; qal = qa; }
            }

            seq.push_back(sel);

            // An N carries no alternative worth trying, so only a
            // disagreement between two called bases is recorded.
            if (a != 'N' && b != 'N') {
                result.mismatches.push_back(
                    PearAlternativeBase{i, sel, alt, qs, qal});
            }
        }

        for (int k = overlap_len; k < n2; ++k) {  // read 2 after the overlap
            seq.push_back(rc.seq[k]);
            qual.push_back(rc.qual[k]);
        }

        result.assembled = true;
        result.seq = std::move(seq);
        result.qual = std::move(qual);
        result.overlap = best_overlap;
        result.r2_start = pos;

        return result;
    }

private:
    int qual_of(char c) const {
        int q = (int)(unsigned char)c - cfg_.phred_base;
        return std::max(0, std::min(q, cfg_.max_qual));
    }

    char to_ascii(int q) const {
        q = std::max(0, std::min(q, cfg_.max_qual));
        return (char)(q + cfg_.phred_base);
    }

    FlashConfig cfg_;
};