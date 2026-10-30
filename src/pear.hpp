#ifndef PEAR_HPP
#define PEAR_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "pvalue.h"


struct PearFastqRecord {
    std::string seq;
    std::string qual;
};

using PearFastqDB =
    std::unordered_map<std::string, PearFastqRecord>;

struct PearEmpFreq {
    uint64_t freqa = 0;
    uint64_t freqc = 0;
    uint64_t freqg = 0;
    uint64_t freqt = 0;
    uint64_t freqn = 0;

    double total = 0.0;
    double pa = 0.25;
    double pc = 0.25;
    double pg = 0.25;
    double pt = 0.25;
    double q  = 0.25;
};

struct PearConfig {
    double p_value       = 0.01;
    int min_overlap      = 10;
    int min_asm_len      = 50;
    int max_asm_len      = 0;   // 0 = unlimited
    double max_uncalled  = 1.0;
    int score_method     = 2;
    int phred_base       = 33;
    bool nbase           = false;
    int test_method      = 1;   // 1: test at min_overlap, 2: observed overlap

    // PEAR's -c / --cap option. 0 disables capping.
    // PEAR's default is 40.
    int cap              = 40;
};

struct PearAlternativeBase {
    int position = -1;           // position in merged sequence
    char selected_base = 'N';    // base PEAR kept
    char alternative_base = 'N'; // base PEAR discarded
    int q_selected = 0;          // internal PEAR Q+1 of selected_base
    int q_alternative = 0;       // internal PEAR Q+1 of alternative_base

    // PEAR always keeps the higher-Q base (R2 on a tie), so this is >= 0.
    // 0 means the two reads disagreed at identical quality.
    int q_gap() const {
        return q_selected - q_alternative;
    }

    bool equal_q() const {
        return q_selected == q_alternative;
    }
};

struct PearMergeResult {
    bool assembled = false;
    std::string seq;
    std::string qual;

    int overlap = 0;
    int r2_start = 0; // start coordinate of reverse-complemented R2 relative to R1
    double score = 0.0;
    double oes = 0.0;

    std::vector<PearAlternativeBase> mismatches;
};


inline const signed char* pear_map_nt_table() {
    static const signed char map_nt[256] = {
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 15, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 15,
        -1,  1, 14,  2, 13, -1, -1,  4, 11, -1, -1, 12, -1,  3, 15, 15,
        -1, -1,  5,  6,  8,  8,  7,  9, 15, 10, -1, -1, -1, -1, -1, -1,
        -1,  1, 14,  2, 13, -1, -1,  4, 11, -1, -1, 12, -1,  3, 15, 15,
        -1, -1,  5,  6,  8,  8,  7,  9, 15, 10, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
    };

    return map_nt;
}

inline const char* pear_revmap_nt_table() {
    static const char revmap_nt[16] = {
        '\0', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
        'T',  'W', 'Y', 'H', 'K', 'D', 'B', 'N'
    };

    return revmap_nt;
}

inline const int* pear_cpl_nt_table() {
    static const int cpl_nt[16] = {
        0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15
    };

    return cpl_nt;
}

// ASCII base -> 4-bit mask.
// PEAR aborts on illegal characters; here they are treated as N (15).
inline int pear_base_mask(char c) {
    const signed char m =
        pear_map_nt_table()[static_cast<unsigned char>(c)];

    return (m >= 0) ? static_cast<int>(m) : 15;
}

// 4-bit mask -> ASCII base.
inline char pear_mask_to_base(int mask) {
    if (mask <= 0 || mask > 15)
        return 'N';

    return pear_revmap_nt_table()[mask];
}

inline bool pear_degenerate_mask(int x) {
    return x != 1 && x != 2 && x != 4 && x != 8;
}

inline std::string pear_trim_line(const std::string& s) {
    size_t start = 0;
    size_t end = s.size();

    while (start < end &&
           (s[start] == '\n' || s[start] == '\r' ||
            s[start] == ' '  || s[start] == '\t')) {
        ++start;
    }

    while (end > start &&
           (s[end - 1] == '\n' || s[end - 1] == '\r' ||
            s[end - 1] == ' '  || s[end - 1] == '\t')) {
        --end;
    }

    return s.substr(start, end - start);
}

inline std::string pear_normalize_header(const std::string& s) {
    std::string t = pear_trim_line(s);
    size_t pos = t.find_first_of(" \t");
    if (pos != std::string::npos)
        t.resize(pos);
    return t;
}

inline PearFastqDB pear_load_fastq(const char* path) {
    PearFastqDB db;

    std::ifstream fin(path);
    if (!fin) {
        std::fprintf(stderr, "Error: Cannot open FASTQ: %s\n", path);
        return db;
    }

    std::string header, seq, plus, qual;

    while (std::getline(fin, header)) {
        if (!std::getline(fin, seq) ||
            !std::getline(fin, plus) ||
            !std::getline(fin, qual)) {
            std::fprintf(
                stderr,
                "Warning: incomplete FASTQ record in %s\n",
                path);
            break;
        }

        header = pear_normalize_header(header);
        seq  = pear_trim_line(seq);
        qual = pear_trim_line(qual);

        if (header.empty())
            continue;

        if (seq.size() != qual.size()) {
            std::fprintf(
                stderr,
                "Warning: sequence/quality length mismatch: %s\n",
                header.c_str());
            continue;
        }

        db[header] = {seq, qual};
    }

    return db;
}

inline char pear_complement(char c) {
    const signed char m =
        pear_map_nt_table()[static_cast<unsigned char>(c)];

    if (m < 0)
        return 'N';

    return pear_mask_to_base(
        pear_cpl_nt_table()[static_cast<int>(m)]);
}

inline void pear_reverse_complement(
    const PearFastqRecord& input,
    PearFastqRecord& output)
{
    const size_t n = input.seq.size();

    output.seq.resize(n);
    output.qual.resize(n);

    for (size_t i = 0; i < n; ++i) {
        output.seq[i] =
            pear_complement(input.seq[n - 1 - i]);

        output.qual[i] =
            input.qual[n - 1 - i];
    }
}


struct PearPairedLoadResult {
    PearFastqDB r1;
    PearFastqDB r2;
    PearEmpFreq ef;

    size_t total_pairs = 0;
    size_t selected_pairs = 0;
    size_t header_mismatch = 0;
};

inline void pear_count_sequence(
    const std::string& seq,
    PearEmpFreq& ef)
{
    for (char c : seq) {
        switch (c) {
            case 'A': case 'a': ++ef.freqa; break;
            case 'C': case 'c': ++ef.freqc; break;
            case 'G': case 'g': ++ef.freqg; break;
            case 'T': case 't': ++ef.freqt; break;
            default:            ++ef.freqn; break;
        }
    }
}

inline PearPairedLoadResult pear_load_paired_selected(
    const char* r1_path,
    const char* r2_path,
    const std::unordered_set<std::string>& wanted_headers)
{
    PearPairedLoadResult out;

    std::ifstream f1(r1_path);
    std::ifstream f2(r2_path);

    if (!f1) {
        std::fprintf(stderr, "Error: Cannot open R1 FASTQ: %s\n", r1_path);
        return out;
    }

    if (!f2) {
        std::fprintf(stderr, "Error: Cannot open R2 FASTQ: %s\n", r2_path);
        return out;
    }

    out.r1.reserve(wanted_headers.size() * 2);
    out.r2.reserve(wanted_headers.size() * 2);

    std::string h1, s1, p1, q1;
    std::string h2, s2, p2, q2;

    while (std::getline(f1, h1) &&
           std::getline(f2, h2))
    {
        if (!std::getline(f1, s1) ||
            !std::getline(f1, p1) ||
            !std::getline(f1, q1) ||
            !std::getline(f2, s2) ||
            !std::getline(f2, p2) ||
            !std::getline(f2, q2)) {
            std::fprintf(
                stderr,
                "Warning: incomplete paired FASTQ record\n");
            break;
        }

        ++out.total_pairs;

        h1 = pear_normalize_header(h1);
        h2 = pear_normalize_header(h2);

        s1 = pear_trim_line(s1);
        q1 = pear_trim_line(q1);
        s2 = pear_trim_line(s2);
        q2 = pear_trim_line(q2);

        // PEAR empirical frequency calculation uses both input files.
        pear_count_sequence(s1, out.ef);
        pear_count_sequence(s2, out.ef);

        if (h1 != h2)
            ++out.header_mismatch;

        // sampled_idx header는 R1의 normalized header를 기준으로 lookup.
        if (!wanted_headers.count(h1))
            continue;

        if (s1.size() != q1.size() ||
            s2.size() != q2.size()) {
            continue;
        }

        out.r1[h1] = {s1, q1};
        out.r2[h1] = {s2, q2};
        ++out.selected_pairs;
    }

    out.ef.total =
        static_cast<double>(
            out.ef.freqa +
            out.ef.freqc +
            out.ef.freqg +
            out.ef.freqt);

    if (out.ef.total > 0.0) {
        out.ef.pa = out.ef.freqa / out.ef.total;
        out.ef.pc = out.ef.freqc / out.ef.total;
        out.ef.pg = out.ef.freqg / out.ef.total;
        out.ef.pt = out.ef.freqt / out.ef.total;

        out.ef.q =
            out.ef.pa * out.ef.pa +
            out.ef.pc * out.ef.pc +
            out.ef.pg * out.ef.pg +
            out.ef.pt * out.ef.pt;
    }

    return out;
}

inline PearEmpFreq pear_compute_empirical_freqs(
    const PearFastqDB& r1,
    const PearFastqDB& r2)
{
    PearEmpFreq ef;

    auto count_db = [&ef](const PearFastqDB& db) {
        for (const auto& kv : db) {
            for (char c : kv.second.seq) {
                switch (c) {
                    case 'A': case 'a': ++ef.freqa; break;
                    case 'C': case 'c': ++ef.freqc; break;
                    case 'G': case 'g': ++ef.freqg; break;
                    case 'T': case 't': ++ef.freqt; break;
                    default:            ++ef.freqn; break;
                }
            }
        }
    };

    count_db(r1);
    count_db(r2);

    ef.total =
        static_cast<double>(
            ef.freqa + ef.freqc + ef.freqg + ef.freqt);

    if (ef.total > 0.0) {
        ef.pa = ef.freqa / ef.total;
        ef.pc = ef.freqc / ef.total;
        ef.pg = ef.freqg / ef.total;
        ef.pt = ef.freqt / ef.total;

        ef.q =
            ef.pa * ef.pa +
            ef.pc * ef.pc +
            ef.pg * ef.pg +
            ef.pt * ef.pt;
    }

    return ef;
}

class PearMerger {
public:
    PearMerger(
        const PearEmpFreq& empirical_freqs,
        const PearConfig& config = PearConfig())
        : ef_(empirical_freqs),
          cfg_(config)
    {
        init_scores();
    }

    PearMergeResult merge(
        const PearFastqRecord& forward_ascii,
        const PearFastqRecord& reverse_ascii) const
    {
        PearMergeResult result;

        if (forward_ascii.seq.empty() ||
            reverse_ascii.seq.empty() ||
            forward_ascii.seq.size() != forward_ascii.qual.size() ||
            reverse_ascii.seq.size() != reverse_ascii.qual.size()) {
            return result;
        }

        PearFastqRecord reverse_rc;
        pear_reverse_complement(reverse_ascii, reverse_rc);

        const int n1 = static_cast<int>(forward_ascii.seq.size());
        const int n2 = static_cast<int>(reverse_rc.seq.size());

        // PEAR internally stores quality values as Q+1.
        std::vector<int> q1(n1);
        std::vector<int> q2(n2);

        for (int i = 0; i < n1; ++i) {
            q1[i] =
                static_cast<unsigned char>(forward_ascii.qual[i])
                - cfg_.phred_base + 1;

            q1[i] = std::clamp(q1[i], 1, 63);
        }

        for (int i = 0; i < n2; ++i) {
            q2[i] =
                static_cast<unsigned char>(reverse_rc.qual[i])
                - cfg_.phred_base + 1;

            q2[i] = std::clamp(q2[i], 1, 63);
        }

        bool found = false;
        double best_score = 0.0;
        double best_oes = 0.0;
        int best_start = 0;
        int best_overlap = 0;

        // r2_start is the position where reverse-complemented R2 starts
        // relative to R1. Scanning all relative offsets covers the same
        // overlap geometries handled separately by PEAR's equal/forward-
        // longer/reverse-longer functions.
        for (int r2_start = -(n2 - 1);
             r2_start <= n1 - 1;
             ++r2_start)
        {
            const int left =
                std::max(0, r2_start);

            const int right =
                std::min(n1, r2_start + n2);

            const int overlap =
                right - left;

            if (overlap <= 0)
                continue;

            double score = 0.0;
            double oes = 0.0;

            for (int p = left; p < right; ++p) {
                const int r1_idx = p;
                const int r2_idx = p - r2_start;

                scoring_ef_nm(
                    pear_base_mask(forward_ascii.seq[r1_idx]),
                    pear_base_mask(reverse_rc.seq[r2_idx]),
                    q1[r1_idx],
                    q2[r2_idx],
                    cfg_.score_method,
                    score,
                    oes);
            }

            // Original PEAR selects a new candidate only for score > best_score.
            if (score > best_score) {
                found = true;
                best_score = score;
                best_oes = oes;
                best_start = r2_start;
                best_overlap = overlap;
            }
        }

        if (!found)
            return result;

        const int test_overlap =
            (cfg_.test_method == 1)
                ? cfg_.min_overlap
                : best_overlap;

        if (!stat_test2(
                cfg_.p_value,
                best_oes,
                test_overlap,
                ef_.q)) {
            return result;
        }

        // PEAR's assembly geometry always keeps coordinate 0 at the
        // beginning of R1 and trims any R2-only prefix before R1.
        // If R2 ends before R1, PEAR also trims the trailing R1 region.
        const int asm_len =
            best_start + n2;

        if (asm_len <= 0)
            return result;

        if (best_overlap < cfg_.min_overlap)
            return result;

        if (cfg_.min_asm_len > 0 &&
            asm_len < cfg_.min_asm_len) {
            return result;
        }

        if (cfg_.max_asm_len > 0 &&
            asm_len > cfg_.max_asm_len) {
            return result;
        }

        std::string merged_seq(
            static_cast<size_t>(asm_len),
            'N');

        std::vector<int> merged_q(
            static_cast<size_t>(asm_len),
            1);

        int uncalled = 0;

        for (int pos = 0; pos < asm_len; ++pos) {
            const bool has_r1 =
                pos >= 0 && pos < n1;

            const int r2_idx =
                pos - best_start;

            const bool has_r2 =
                r2_idx >= 0 && r2_idx < n2;

            if (has_r1 && has_r2) {
                const char a =
                    forward_ascii.seq[pos];

                const char b =
                    reverse_rc.seq[r2_idx];

                const int ma =
                    pear_base_mask(a);

                const int mb =
                    pear_base_mask(b);

                const int qa =
                    q1[pos];

                const int qb =
                    q2[r2_idx];

                if (pear_degenerate_mask(ma) &&
                    pear_degenerate_mask(mb))
                {
                    // PEAR: ai->data = x | y, qscore = min(qx, qy).
                    // The union of the two IUPAC masks is kept rather
                    // than collapsing unconditionally to 'N'.
                    merged_seq[pos] =
                        pear_mask_to_base(ma | mb);

                    merged_q[pos] = std::min(qa, qb);
                    ++uncalled;
                }
                else if (pear_degenerate_mask(ma))
                {
                    merged_seq[pos] = b;
                    merged_q[pos] = qb;
                }
                else if (pear_degenerate_mask(mb))
                {
                    merged_seq[pos] = a;
                    merged_q[pos] = qa;
                }
                else if (ma == mb)
                {
                    // PEAR: qscore = qx + qy - 1, with no cap at this
                    // stage. Capping happens once at output time.
                    merged_seq[pos] = a;
                    merged_q[pos] = qa + qb - 1;
                }
                else
                {
                    if (cfg_.nbase) {
                        merged_seq[pos] = 'N';
                        merged_q[pos] = std::max(qa, qb);
                        ++uncalled;
                    }
                    else if (qa > qb) {
                        merged_seq[pos] = a;
                        merged_q[pos] = qa;

                        // R1 won on quality; R2 is the alternative.
                        result.mismatches.push_back(
                            PearAlternativeBase{
                                pos, a, b, qa, qb
                            });
                    }
                    else {
                        // Same tie-breaking as PEAR:
                        // qleft <= qright -> choose right.
                        merged_seq[pos] = b;
                        merged_q[pos] = qb;

                        // R2 won on quality, or the two tied.
                        // R1 is the alternative either way.
                        // This position is inside the overlap because
                        // both has_r1 and has_r2 are true here.
                        result.mismatches.push_back(
                            PearAlternativeBase{
                                pos, b, a, qb, qa
                            });
                    }
                }
            }
            else if (has_r1) {
                merged_seq[pos] =
                    forward_ascii.seq[pos];

                merged_q[pos] =
                    q1[pos];

                if (pear_degenerate_mask(
                        pear_base_mask(merged_seq[pos]))) {
                    ++uncalled;
                }
            }
            else if (has_r2) {
                merged_seq[pos] =
                    reverse_rc.seq[r2_idx];

                merged_q[pos] =
                    q2[r2_idx];

                if (pear_degenerate_mask(
                        pear_base_mask(merged_seq[pos]))) {
                    ++uncalled;
                }
            }
        }

        const double uncalled_ratio =
            static_cast<double>(uncalled) /
            static_cast<double>(asm_len);

        if (uncalled_ratio > cfg_.max_uncalled)
            return result;

        // PEAR's add_base_qscores():
        //   if (cap) { cap = cap + 1; if (q > cap) q = cap; }
        //   ascii = q + phred_base - 1;
        std::string merged_qual(
            static_cast<size_t>(asm_len),
            '!');

        const int cap_internal =
            (cfg_.cap > 0) ? cfg_.cap + 1 : 0;

        for (int i = 0; i < asm_len; ++i) {
            int q = merged_q[i];

            if (cap_internal > 0 && q > cap_internal)
                q = cap_internal;

            int ascii =
                q + cfg_.phred_base - 1;

            // Defensive only. With the default cap this never binds;
            // it prevents a non-printable byte when capping is disabled.
            ascii = std::clamp(ascii, 33, 126);

            merged_qual[i] =
                static_cast<char>(ascii);
        }

        result.assembled = true;
        result.seq = std::move(merged_seq);
        result.qual = std::move(merged_qual);
        result.overlap = best_overlap;
        result.r2_start = best_start;
        result.score = best_score;
        result.oes = best_oes;

        return result;
    }

private:
    PearEmpFreq ef_;
    PearConfig cfg_;

    // PEAR quality/base-frequency score table.
    std::array<std::array<double, 4096>, 10> penalty_{};

    void init_scores() {
        const double pa2 = ef_.pa * ef_.pa;
        const double pc2 = ef_.pc * ef_.pc;
        const double pg2 = ef_.pg * ef_.pg;
        const double pt2 = ef_.pt * ef_.pt;

        const double pacg2 = pa2 + pc2 + pg2;
        const double pact2 = pa2 + pc2 + pt2;
        const double pagt2 = pa2 + pg2 + pt2;
        const double pcgt2 = pc2 + pg2 + pt2;

        const double pacg = ef_.pa + ef_.pc + ef_.pg;
        const double pact = ef_.pa + ef_.pc + ef_.pt;
        const double pagt = ef_.pa + ef_.pg + ef_.pt;
        const double pcgt = ef_.pc + ef_.pg + ef_.pt;

        const double p2acg = pacg * pacg;
        const double p2act = pact * pact;
        const double p2agt = pagt * pagt;
        const double p2cgt = pcgt * pcgt;

        const double pac2 =
            (ef_.pa + ef_.pc) * (ef_.pa + ef_.pc);

        const double pag2 =
            (ef_.pa + ef_.pg) * (ef_.pa + ef_.pg);

        const double pat2 =
            (ef_.pa + ef_.pt) * (ef_.pa + ef_.pt);

        const double pcg2 =
            (ef_.pc + ef_.pg) * (ef_.pc + ef_.pg);

        const double pct2 =
            (ef_.pc + ef_.pt) * (ef_.pc + ef_.pt);

        const double pgt2 =
            (ef_.pg + ef_.pt) * (ef_.pg + ef_.pt);

        for (int i = 1; i < 64; ++i) {
            for (int j = 1; j < 64; ++j) {
                const double ex =
                    std::pow(10.0, -(i - 1) / 10.0);

                const double ey =
                    std::pow(10.0, -(j - 1) / 10.0);

                const int k =
                    (i << 6) | j;

                penalty_[0][k] =
                    (1 - ex) * (1 - ey) +
                    (ex * ey) * pcgt2 / p2cgt;

                penalty_[1][k] =
                    (1 - ex) * (1 - ey) +
                    (ex * ey) * pagt2 / p2agt;

                penalty_[2][k] =
                    (1 - ex) * (1 - ey) +
                    (ex * ey) * pact2 / p2act;

                penalty_[3][k] =
                    (1 - ex) * (1 - ey) +
                    (ex * ey) * pacg2 / p2acg;

                penalty_[4][k] =
                    1
                    - (1 - ey) * ex * (ef_.pc / pcgt)
                    - (1 - ex) * ey * (ef_.pa / pagt)
                    - ex * ey * (pg2 + pt2) / pgt2;

                penalty_[5][k] =
                    1
                    - (1 - ey) * ex * (ef_.pg / pcgt)
                    - (1 - ex) * ey * (ef_.pa / pact)
                    - ex * ey * (pc2 + pt2) / pct2;

                penalty_[6][k] =
                    1
                    - (1 - ey) * ex * (ef_.pt / pcgt)
                    - (1 - ex) * ey * (ef_.pa / pacg)
                    - ex * ey * (pc2 + pg2) / pcg2;

                penalty_[7][k] =
                    1
                    - (1 - ey) * ex * (ef_.pg / pagt)
                    - (1 - ex) * ey * (ef_.pc / pact)
                    - ex * ey * (pa2 + pt2) / pat2;

                penalty_[8][k] =
                    1
                    - (1 - ey) * ex * (ef_.pt / pagt)
                    - (1 - ex) * ey * (ef_.pc / pacg)
                    - ex * ey * (pa2 + pg2) / pag2;

                penalty_[9][k] =
                    1
                    - (1 - ey) * ex * (ef_.pt / pact)
                    - (1 - ex) * ey * (ef_.pg / pacg)
                    - ex * ey * (pa2 + pc2) / pac2;
            }
        }
    }

    static int translate_pair(int a, int b) {
        static constexpr int translate[16] = {
            -1, 0, 1, 4, 2, 5, 7, -1,
             3, 6, 8, -1, 9, -1, -1, -1
        };

        return translate[a | b];
    }

    void scoring_ef_nm(
        int dleft,
        int dright,
        int qleft,
        int qright,
        int score_method,
        double& score,
        double& oes) const
    {
        if (pear_degenerate_mask(dleft) ||
            pear_degenerate_mask(dright))
        {
            switch (score_method) {
                case 1:
                    score += ef_.q - (1 - ef_.q);
                    oes = score;
                    break;

                case 2: {
                    const double tmp = 1 - ef_.q;
                    oes += ef_.q - tmp;
                    score -= tmp;
                    break;
                }

                case 3:
                    score -= 1;
                    oes += ef_.q - (1 - ef_.q);
                    break;
            }

            return;
        }

        if (dleft == dright) {
            const int idx =
                translate_pair(dleft, dright);

            if (idx < 0)
                return;

            const int k =
                (qright << 6) | qleft;

            switch (score_method) {
                case 1:
                    score +=
                        penalty_[idx][k]
                        - (1 - penalty_[idx][k]);
                    oes = score;
                    break;

                case 2:
                    score += penalty_[idx][k];
                    oes +=
                        penalty_[idx][k]
                        - (1 - penalty_[idx][k]);
                    break;

                case 3:
                    oes +=
                        penalty_[idx][k]
                        - (1 - penalty_[idx][k]);
                    score += 1;
                    break;
            }

            return;
        }

        if (dleft > dright) {
            std::swap(dleft, dright);
            std::swap(qleft, qright);
        }

        const int idx =
            translate_pair(dleft, dright);

        if (idx < 0)
            return;

        const int k =
            (qleft << 6) | qright;

        switch (score_method) {
            case 1:
                score =
                    score - 2 * penalty_[idx][k] + 1;
                oes = score;
                break;

            case 2:
                score -= penalty_[idx][k];
                oes =
                    oes - 2 * penalty_[idx][k] + 1;
                break;

            case 3:
                score -= 1;
                oes =
                    oes - 2 * penalty_[idx][k] + 1;
                break;
        }
    }

    static bool stat_test2(
        double pval,
        double oes,
        int min_overlap,
        double q)
    {
        if (min_overlap > 99)
            min_overlap = 99;

        if (min_overlap < 0)
            min_overlap = 0;

        if (q > 0.49)
            q = 0.49;

        if (q < 0.0)
            q = 0.0;

        if (pval == 1.0)
            return true;

        double (*table)[100] = nullptr;

        if (pval == 0.01)
            table = precomp2_01;
        else if (pval == 0.05)
            table = precomp2_05;
        else if (pval == 0.001)
            table = precomp2_001;
        else
            table = precomp2_0001;

        const double cutoff =
            table[min_overlap][static_cast<int>(q * 100)];

        return oes > cutoff;
    }
};

struct PearAltMask {
    int mask = 0;        // bit i set => flip plan.positions[i]
    int flip_count = 0;
    int q_gap_sum = 0;   // sum of q_gap() over flipped positions
    int q_sel_sum = 0;   // sum of q_selected over flipped positions
};

struct PearAltPlan {
    // Positions actually searched, already ranked.
    std::vector<PearAlternativeBase> positions;

    // All 2^positions.size() - 1 non-empty combinations, ranked.
    std::vector<PearAltMask> masks;

    int total_mismatches = 0;    // every recorded mismatch
    int eligible_mismatches = 0; // those within max_q_gap

    // total_mismatches > max_positions. Such a read is not searched.
    bool overflow = false;
};


inline PearAltPlan pear_build_alt_plan(
    const std::vector<PearAlternativeBase>& mismatches,
    int max_positions,
    int max_q_gap = -1)
{
    PearAltPlan plan;

    plan.total_mismatches =
        static_cast<int>(mismatches.size());

    if (mismatches.empty() || max_positions <= 0)
        return plan;

    // Stage 1: reject reads carrying too many mismatches.
    if (plan.total_mismatches > max_positions) {
        plan.overflow = true;
        return plan;
    }

    // Stage 2: keep the positions within the quality bound.
    if (max_q_gap < 0) {
        plan.positions = mismatches;
    }
    else {
        for (const PearAlternativeBase& m : mismatches)
            if (m.q_gap() <= max_q_gap)
                plan.positions.push_back(m);

        if (plan.positions.empty())
            return plan;
    }

    plan.eligible_mismatches =
        static_cast<int>(plan.positions.size());

    std::stable_sort(
        plan.positions.begin(),
        plan.positions.end(),
        [](const PearAlternativeBase& lhs,
           const PearAlternativeBase& rhs)
        {
            const int lg = lhs.q_gap();
            const int rg = rhs.q_gap();

            if (lg != rg)
                return lg < rg;

            if (lhs.q_selected != rhs.q_selected)
                return lhs.q_selected < rhs.q_selected;

            return lhs.position < rhs.position;
        });

    const int n =
        static_cast<int>(plan.positions.size());

    const int max_mask = 1 << n;

    plan.masks.reserve(
        static_cast<size_t>(max_mask - 1));

    for (int mask = 1; mask < max_mask; ++mask) {
        PearAltMask m;
        m.mask = mask;

        for (int i = 0; i < n; ++i) {
            if ((mask & (1 << i)) == 0)
                continue;

            ++m.flip_count;
            m.q_gap_sum += plan.positions[i].q_gap();
            m.q_sel_sum += plan.positions[i].q_selected;
        }

        plan.masks.push_back(m);
    }

    std::sort(
        plan.masks.begin(),
        plan.masks.end(),
        [](const PearAltMask& lhs, const PearAltMask& rhs)
        {
            if (lhs.q_gap_sum != rhs.q_gap_sum)
                return lhs.q_gap_sum < rhs.q_gap_sum;

            if (lhs.flip_count != rhs.flip_count)
                return lhs.flip_count < rhs.flip_count;

            if (lhs.q_sel_sum != rhs.q_sel_sum)
                return lhs.q_sel_sum < rhs.q_sel_sum;

            return lhs.mask < rhs.mask;
        });

    return plan;
}

#endif // PEAR_HPP