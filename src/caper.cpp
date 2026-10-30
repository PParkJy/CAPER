#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <filesystem>
#include <atomic>
#include <mutex>
#include <thread>
#include <memory>

#include "pear.hpp"
#include "ngmerge.hpp"
#include "flash.hpp"

static int LT_L          = 256;
static constexpr int LT_SEED_BIT   = 32;
static constexpr int RS_PARITY_NUM = 2;
static int LT_K          = 16050;
static int STRAND_WIDTH  = LT_SEED_BIT + 256 + 8 * RS_PARITY_NUM;
static int MAX_STRANDS   = 80000;
static int RS_CODE_BYTES = 38;
static int DNA_LENGTH    = 152;

static constexpr int MAX_CODE_BYTES = 64;

static const char* const DATA_DIR = "../data";

static const char* const R1_FASTQ_NAME = "OP2-2_S1_L001_R1_RMother.fastq";
static const char* const R2_FASTQ_NAME = "OP2-2_S1_L001_R2_RMother.fastq";

static const char* const DEFAULT_IMAGE_FILE =
    "../data/image_restart_pad.txt";

static const char* const DEFAULT_SEED_TABLE_FILE =
    "../data/all_seed_to_linear_table_noconst.txt";

static const char* g_r1_path_opt    = nullptr;
static const char* g_r2_path_opt    = nullptr;
static const char* g_seed_table_opt = nullptr;
static const char* g_image_opt      = nullptr;
static const char* g_base_dir_opt   = nullptr;

static int         g_tool_threads = 5;
static int         g_msa_jobs     = 1;
static std::string g_tmp_dir      = "/tmp";
static int         g_merge_cache  = 1;
static int         g_merge_jobs   = 1;

static std::vector<char> tmp_template(const char* stem) {
    std::string t = g_tmp_dir + "/" + stem + "_XXXXXX";
    return std::vector<char>(t.begin(), t.end() + 1);
}



static constexpr int MAX_ALT_CAP           = 20;
static constexpr int DEFAULT_MAX_ALT       = 5;
static int           g_max_alt             = DEFAULT_MAX_ALT;

static constexpr int DEFAULT_MAX_FLIPS = 5;
static int  g_max_flips     = DEFAULT_MAX_FLIPS;
static bool g_max_flips_set = false;

static inline int alt_flip_budget() {
    return std::min(g_max_flips, g_max_alt);
}

static constexpr int DEFAULT_MAX_QSUM = 20;
static int g_max_qsum = DEFAULT_MAX_QSUM;

static int g_i1_mm_order = 1;

static constexpr int DEFAULT_STARCODE_DISTANCE = 7;

static constexpr int DEFAULT_MULTILEN = 5;

static const char* const RESULT_DIR = "../result";

static const char* const TOOL_DIR = "./sequence_analysis";

static const char* const STARCODE_CANDIDATES[] = {
    "./sequence_analysis/starcode",
    "./sequence_analysis/starcode/starcode",
    nullptr
};

static const char* const CDHIT_CANDIDATES[] = {
    "./sequence_analysis/cd-hit-est",
    "./sequence_analysis/cdhit/cd-hit-est",
    "./sequence_analysis/cd-hit/cd-hit-est",
    nullptr
};

static const char* const MSA_CANDIDATES[] = {
    "./sequence_analysis/muscle",
    "./sequence_analysis/muscle-linux-x86.v5.3",
    "./sequence_analysis/muscle/muscle",
    "./sequence_analysis/muscle/muscle-linux-x86.v5.3",
    nullptr
};

static const char* find_tool(const char* const* candidates) {
    for (const char* const* c = candidates; *c; ++c) {
        std::error_code ec;

        if (std::filesystem::is_regular_file(*c, ec) &&
            access(*c, X_OK) == 0)
            return *c;
    }

    return candidates[0];
}

static const char* g_starcode_bin = nullptr;
static const char* g_cdhit_bin    = nullptr;

enum ClusterTool {
    CLUSTER_TOOL_STARCODE = 0,
    CLUSTER_TOOL_CDHIT    = 1
};

static int    g_cluster_tool = CLUSTER_TOOL_STARCODE;
static double g_cdhit_id     = 0.94;

enum AltMode {
    ALT_OFF     = 0,
    ALT_EQUAL_Q = 1,
    ALT_ALL     = 2
};

enum ClusterMode {
    CLUSTER_EXACT    = 0,
    CLUSTER_STARCODE = 1
};

enum MergerMode {
    MERGER_PEAR     = 0,
    MERGER_EXTERNAL = 1,
    MERGER_NGMERGE  = 2,
    MERGER_FLASH    = 3
};

static int g_flash_max_overlap = 150;

static const char* const DEFAULT_NGMERGE_PROFILE = "qual_profile.txt";

static const char* g_ngmerge_profile = DEFAULT_NGMERGE_PROFILE;
static bool        g_ngmerge_profile_set = false;

static double g_ngmerge_mismatch = 0.40;
static double g_pear_pvalue      = 0.01;
static double g_flash_mismatch   = 0.40;

static const char* g_tag = nullptr;

static bool alt_rank_less(const PearAlternativeBase& lhs,
                          const PearAlternativeBase& rhs)
{
    const int lg = lhs.q_gap();
    const int rg = rhs.q_gap();

    if (lg != rg)
        return lg < rg;

    if (lhs.q_selected != rhs.q_selected)
        return lhs.q_selected < rhs.q_selected;

    return lhs.position < rhs.position;
}

static const char* const KNOWN_MERGERS[] = {
    "flash", "ngmerge", "bbmerge", "fastp", "vsearch", "usearch"
};

static bool is_known_merger(const char* name) {
    for (const char* k : KNOWN_MERGERS)
        if (strcmp(name, k) == 0)
            return true;
    return false;
}

struct Dataentry {
    unsigned long seed;
    std::vector<int> indices;
};

int acgt_to_bits(const char* buf, int len, int* out) {
    int idx = 0;

    for (int i = 0; i < len; i++) {
        switch (buf[i]) {
            case 'A':
                out[idx++] = 0;
                out[idx++] = 0;
                break;
            case 'C':
                out[idx++] = 0;
                out[idx++] = 1;
                break;
            case 'G':
                out[idx++] = 1;
                out[idx++] = 0;
                break;
            case 'T':
                out[idx++] = 1;
                out[idx++] = 1;
                break;
            default:
                return -1;
        }
    }

    return idx;
}

bool dna_to_rs_bytes(
    const std::string& seq,
    uint8_t out[])
{
    if ((int)seq.size() != DNA_LENGTH)
        return false;

    for (int i = 0; i < RS_CODE_BYTES; i++) {
        uint8_t value = 0;

        for (int j = 0; j < 4; j++) {
            uint8_t two_bits;

            switch (seq[i * 4 + j]) {
                case 'A': two_bits = 0; break;
                case 'C': two_bits = 1; break;
                case 'G': two_bits = 2; break;
                case 'T': two_bits = 3; break;
                default:
                    return false;
            }

            value =
                static_cast<uint8_t>(
                    (value << 2) | two_bits);
        }

        out[i] = value;
    }

    return true;
}

std::vector<std::string>
scan_fastq_headers(const char* path)
{
    std::vector<std::string> headers;

    std::ifstream fin(path);

    if (!fin) {
        fprintf(
            stderr,
            "Error: Cannot open FASTQ: %s\n",
            path);

        return headers;
    }

    std::string line;
    long long line_no = 0;

    while (std::getline(fin, line)) {
        if (line_no % 4 == 0) {
            std::string h =
                pear_normalize_header(line);

            if (!h.empty())
                headers.push_back(h);
        }

        ++line_no;
    }

    return headers;
}

PearFastqDB load_merged_selected(
    const char* path,
    const std::unordered_set<std::string>& wanted,
    size_t& total_records,
    size_t& kept_records)
{
    PearFastqDB db;
    total_records = 0;
    kept_records = 0;

    std::ifstream fin(path);

    if (!fin) {
        fprintf(
            stderr,
            "Error: Cannot open merged FASTQ: %s\n",
            path);

        return db;
    }

    db.reserve(wanted.size() * 2);

    std::string header, seq, plus, qual;

    while (std::getline(fin, header)) {
        if (!std::getline(fin, seq) ||
            !std::getline(fin, plus) ||
            !std::getline(fin, qual))
        {
            fprintf(
                stderr,
                "Warning: incomplete record in %s\n",
                path);
            break;
        }

        ++total_records;

        header = pear_normalize_header(header);
        seq  = pear_trim_line(seq);
        qual = pear_trim_line(qual);

        if (header.empty() || !wanted.count(header))
            continue;

        if (seq.size() != qual.size())
            continue;

        db[header] = {seq, qual};
        ++kept_records;
    }

    return db;
}

std::vector<std::string>
sample_headers(
    const std::vector<std::string>& pool,
    size_t sample_n,
    unsigned long long seed)
{
    const size_t n = pool.size();

    std::vector<std::string> out;

    if (n == 0)
        return out;

    if (sample_n > n)
        sample_n = n;

    std::mt19937_64 rng(seed);

    std::unordered_map<size_t, size_t> swapped;
    swapped.reserve(sample_n * 2);

    auto value_at =
        [&](size_t i) -> size_t {
            auto it = swapped.find(i);
            return (it == swapped.end()) ? i : it->second;
        };

    out.reserve(sample_n);

    for (size_t i = 0; i < sample_n; ++i) {
        std::uniform_int_distribution<size_t>
            dist(i, n - 1);

        const size_t j = dist(rng);

        const size_t vj = value_at(j);
        const size_t vi = value_at(i);

        out.push_back(pool[vj]);

        swapped[j] = vi;
    }

    return out;
}

class GF256 {
public:
    GF256() {
        init();
    }

    uint8_t mul(uint8_t a, uint8_t b) const {
        if (a == 0 || b == 0)
            return 0;

        return exp_table[log_table[a] + log_table[b]];
    }

    uint8_t alpha_pow(int exponent) const {
        exponent %= 255;
        if (exponent < 0)
            exponent += 255;

        return exp_table[exponent];
    }

private:
    uint8_t exp_table[512]{};
    uint8_t log_table[256]{};

    void init() {
        unsigned int x = 1;

        for (int i = 0; i < 255; i++) {
            exp_table[i] = static_cast<uint8_t>(x);
            log_table[x] = static_cast<uint8_t>(i);

            x <<= 1;

            if (x & 0x100)
                x ^= 0x11D;
        }

        for (int i = 255; i < 512; i++)
            exp_table[i] = exp_table[i - 255];

        log_table[0] = 0;
    }
};

bool rs_detect_syndrome(
    const uint8_t codeword[],
    const GF256& gf)
{
    for (int root = 1; root <= RS_PARITY_NUM; root++) {
        const uint8_t x = gf.alpha_pow(root);

        uint8_t syndrome = 0;

        for (int i = 0; i < RS_CODE_BYTES; i++) {
            syndrome =
                static_cast<uint8_t>(
                    gf.mul(syndrome, x) ^ codeword[i]);
        }

        if (syndrome != 0)
            return false;
    }

    return true;
}

struct CachedPearMerge {
    int status = 0;
    PearMergeResult merged;
};

struct StarCluster {
    std::string canonical;
    int cluster_size = 0;
    std::vector<std::string> members;
};

static int dna_base_index(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1;
    }
}

static char dna_index_base(int i) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    return (i >= 0 && i < 4) ? BASES[i] : 'A';
}

struct StarVoteProfile {
    int observations = 0;
    std::map<int, std::array<double, 4>> delta_by_pos;
};

static bool star_q_position_eligible(
    const PearAlternativeBase& m,
    int alt_mode)
{
    if (alt_mode == ALT_OFF)
        return false;

    if (alt_mode == ALT_EQUAL_Q)
        return m.equal_q();

    return true;
}

static std::vector<PearAlternativeBase> star_q_positions(
    const std::vector<PearAlternativeBase>& mismatches,
    int alt_mode)
{
    std::vector<PearAlternativeBase> eligible;

    if (alt_mode == ALT_OFF)
        return eligible;

    for (const PearAlternativeBase& m : mismatches) {
        if (star_q_position_eligible(m, alt_mode))
            eligible.push_back(m);
    }

    if ((int)eligible.size() > g_max_alt) {
        std::stable_sort(
            eligible.begin(),
            eligible.end(),
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

        eligible.resize(g_max_alt);
    }

    return eligible;
}

static void star_split_weight_from_gap(
    int ,
    double& selected_weight,
    double& alternative_weight)
{
    selected_weight = 0.5;
    alternative_weight = 0.5;
}

static void add_star_vote_profile(
    std::map<std::string, StarVoteProfile>& profiles,
    const std::string& output_seq,
    const std::vector<PearAlternativeBase>& mismatches,
    int alt_mode)
{
    if (output_seq.empty())
        return;

    StarVoteProfile& profile = profiles[output_seq];
    ++profile.observations;

    if (alt_mode == ALT_OFF)
        return;

    const std::vector<PearAlternativeBase> eligible =
        star_q_positions(mismatches, alt_mode);

    for (const PearAlternativeBase& m : eligible) {
        if (m.position < 0 ||
            m.position >= (int)output_seq.size())
            continue;

        const int current_b = dna_base_index(output_seq[m.position]);
        const int selected_b = dna_base_index(m.selected_base);
        const int alternative_b = dna_base_index(m.alternative_base);

        if (current_b < 0 || selected_b < 0 || alternative_b < 0)
            continue;

        double w_selected = 1.0;
        double w_alternative = 0.0;
        star_split_weight_from_gap(
            m.q_gap(),
            w_selected,
            w_alternative);

        auto& delta = profile.delta_by_pos[m.position];

        delta[current_b] -= 1.0;
        delta[selected_b] += w_selected;
        delta[alternative_b] += w_alternative;
    }
}

static bool run_starcode_clusters(
    const std::map<std::string, int>& counts,
    std::vector<StarCluster>& out_clusters,
    int distance = DEFAULT_STARCODE_DISTANCE)
{
    out_clusters.clear();

    if (counts.empty())
        return true;

    std::vector<char> tmp_name_buf = tmp_template("dna_starcode");
    char* tmp_name = tmp_name_buf.data();
    const int fd = mkstemp(tmp_name);

    if (fd < 0) {
        perror("mkstemp");
        return false;
    }

    FILE* fp = fdopen(fd, "w");
    if (!fp) {
        perror("fdopen");
        close(fd);
        unlink(tmp_name);
        return false;
    }

    for (const auto& kv : counts) {
        if (kv.second > 0)
            fprintf(fp, "%s\t%d\n", kv.first.c_str(), kv.second);
    }
    fclose(fp);

    std::string cmd =
        std::string(g_starcode_bin) + " -i " + tmp_name +
        " -d " + std::to_string(distance) +
        " -s --print-clusters -t " + std::to_string(g_tool_threads) +
        " -q 2>/dev/null";

    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        perror("popen");
        unlink(tmp_name);
        return false;
    }

    char* line = nullptr;
    size_t cap = 0;

    while (::getline(&line, &cap, pipe) != -1) {
        std::string row(line);
        while (!row.empty() && (row.back() == '\n' || row.back() == '\r'))
            row.pop_back();

        const size_t t1 = row.find('\t');
        if (t1 == std::string::npos)
            continue;

        const size_t t2 = row.find('\t', t1 + 1);

        StarCluster c;
        c.canonical = row.substr(0, t1);

        const std::string count_str =
            (t2 == std::string::npos)
                ? row.substr(t1 + 1)
                : row.substr(t1 + 1, t2 - t1 - 1);

        c.cluster_size = atoi(count_str.c_str());

        if (t2 != std::string::npos) {
            std::string members = row.substr(t2 + 1);
            std::stringstream ss(members);
            std::string seq;

            while (std::getline(ss, seq, ',')) {
                if (!seq.empty())
                    c.members.push_back(seq);
            }
        }

        if (c.members.empty())
            c.members.push_back(c.canonical);

        out_clusters.push_back(std::move(c));
    }

    if (line)
        free(line);

    const int status = pclose(pipe);
    unlink(tmp_name);

    if (status != 0) {
        fprintf(stderr,
                "Warning: starcode failed (status=%d).\n"
                "         Check the binary path: %s\n",
                status, g_starcode_bin);
        out_clusters.clear();
        return false;
    }

    return true;
}

static bool run_cdhit_clusters(
    const std::map<std::string, int>& counts,
    std::vector<StarCluster>& out_clusters,
    double identity)
{
    out_clusters.clear();

    if (counts.empty())
        return true;

    auto word_length = [](double id) {
        if (id >= 0.95) return 10;
        if (id >= 0.90) return 8;
        if (id >= 0.88) return 7;
        if (id >= 0.85) return 6;
        if (id >= 0.80) return 5;
        return 4;
    };

    std::vector<char> in_name_buf  = tmp_template("dna_cdhit_in");
    std::vector<char> out_name_buf = tmp_template("dna_cdhit_out");
    char* in_name  = in_name_buf.data();
    char* out_name = out_name_buf.data();

    const int in_fd = mkstemp(in_name);

    if (in_fd < 0) {
        perror("mkstemp");
        return false;
    }

    const int out_fd = mkstemp(out_name);

    if (out_fd < 0) {
        perror("mkstemp");
        close(in_fd);
        unlink(in_name);
        return false;
    }

    close(out_fd);

    FILE* fp = fdopen(in_fd, "w");

    if (!fp) {
        perror("fdopen");
        close(in_fd);
        unlink(in_name);
        unlink(out_name);
        return false;
    }

    std::vector<const std::string*> seqs;
    seqs.reserve(counts.size());

    for (const auto& kv : counts) {
        if (kv.second <= 0)
            continue;

        fprintf(fp, ">s%zu\n%s\n", seqs.size(), kv.first.c_str());
        seqs.push_back(&kv.first);
    }

    fclose(fp);

    if (seqs.empty()) {
        unlink(in_name);
        unlink(out_name);
        return true;
    }

    char id_buf[32];
    snprintf(id_buf, sizeof(id_buf), "%.4f", identity);

    const std::string cmd =
        std::string(g_cdhit_bin) +
        " -i " + in_name +
        " -o " + out_name +
        " -c " + id_buf +
        " -n " + std::to_string(word_length(identity)) +
        " -r 0 -g 1 -d 0 -M 0 -T " + std::to_string(g_tool_threads) +
        " > /dev/null 2>&1";

    const int rc = system(cmd.c_str());

    const std::string clstr = std::string(out_name) + ".clstr";

    if (rc != 0) {
        fprintf(stderr,
                "Warning: cd-hit-est exited %d.\n"
                "         Check the binary path: %s\n",
                rc, g_cdhit_bin);

        unlink(in_name);
        unlink(out_name);
        unlink(clstr.c_str());
        return false;
    }

    std::ifstream fin(clstr);

    if (!fin) {
        fprintf(stderr,
                "Warning: cd-hit-est produced no cluster file: %s\n",
                clstr.c_str());

        unlink(in_name);
        unlink(out_name);
        return false;
    }

    StarCluster cur;
    bool open_cluster = false;

    auto close_cluster = [&]() {
        if (!open_cluster)
            return;

        if (!cur.members.empty()) {
            if (cur.canonical.empty())
                cur.canonical = cur.members.front();

            out_clusters.push_back(cur);
        }

        cur = StarCluster();
        open_cluster = false;
    };

    std::string line;

    while (std::getline(fin, line)) {
        while (!line.empty() &&
               (line.back() == '\n' || line.back() == '\r'))
            line.pop_back();

        if (line.empty())
            continue;

        if (line[0] == '>') {
            close_cluster();
            open_cluster = true;
            continue;
        }

        const size_t a = line.find(">s");

        if (a == std::string::npos)
            continue;

        const size_t b = line.find("...", a);

        const std::string idx_s =
            line.substr(a + 2,
                        (b == std::string::npos)
                            ? std::string::npos
                            : b - (a + 2));

        const long idx = atol(idx_s.c_str());

        if (idx < 0 || idx >= (long)seqs.size())
            continue;

        const std::string& seq = *seqs[(size_t)idx];

        auto cit = counts.find(seq);

        cur.members.push_back(seq);
        cur.cluster_size += (cit == counts.end()) ? 0 : cit->second;

        if (line.find('*') != std::string::npos)
            cur.canonical = seq;
    }

    close_cluster();

    unlink(in_name);
    unlink(out_name);
    unlink(clstr.c_str());

    return true;
}

static bool run_clusters(
    const std::map<std::string, int>& counts,
    std::vector<StarCluster>& out_clusters,
    int distance)
{
    if (g_cluster_tool == CLUSTER_TOOL_CDHIT)
        return run_cdhit_clusters(counts, out_clusters, g_cdhit_id);

    return run_starcode_clusters(counts, out_clusters, distance);
}

struct MsaConsensus {
    std::string plain;
    std::string weighted;
    std::vector<int>    weak_pos;
    std::vector<char>   weak_runner_up;
    std::vector<double> weak_margin;
};

static bool msa_run_rows(
    const std::vector<std::string>& seqs,
    const char* muscle_bin,
    int threads,
    std::vector<std::string>& rows)
{
    rows.clear();

    std::vector<char> in_buf  = tmp_template("dna_msa_in");
    std::vector<char> out_buf = tmp_template("dna_msa_out");
    char* in_name  = in_buf.data();
    char* out_name = out_buf.data();

    const int in_fd = mkstemp(in_name);
    if (in_fd < 0)
        return false;

    const int out_fd = mkstemp(out_name);
    if (out_fd < 0) {
        close(in_fd);
        unlink(in_name);
        return false;
    }
    close(out_fd);

    FILE* fin = fdopen(in_fd, "w");
    if (!fin) {
        close(in_fd);
        unlink(in_name);
        unlink(out_name);
        return false;
    }

    for (size_t i = 0; i < seqs.size(); ++i)
        fprintf(fin, ">s%zu\n%s\n", i, seqs[i].c_str());

    fclose(fin);

    const std::string cmd =
        std::string(muscle_bin) +
        " -align " + in_name +
        " -output " + out_name +
        " -threads " + std::to_string(threads) +
        " > /dev/null 2>&1";

    const int rc = system(cmd.c_str());

    if (rc != 0) {
        static std::atomic<bool> warned{false};

        if (!warned.exchange(true))
            fprintf(stderr,
                    "Warning: MUSCLE exited %d. Every mixed-length\n"
                    "         cluster will produce nothing. Try the\n"
                    "         command by hand:\n           %s\n",
                    rc, cmd.c_str());
    }

    std::map<std::string, std::string> aligned;

    if (rc == 0) {
        std::ifstream fout(out_name);
        std::string line, id, acc;

        while (std::getline(fout, line)) {
            while (!line.empty() &&
                   (line.back() == '\n' || line.back() == '\r'))
                line.pop_back();

            if (!line.empty() && line[0] == '>') {
                if (!id.empty())
                    aligned[id] = acc;

                id = line.substr(1);

                const size_t cut = id.find_first_of(" \t");
                if (cut != std::string::npos)
                    id.resize(cut);

                acc.clear();
            }
            else {
                acc += line;
            }
        }

        if (!id.empty())
            aligned[id] = acc;
    }

    unlink(in_name);
    unlink(out_name);

    if (aligned.size() != seqs.size())
        return false;

    size_t width = 0;

    for (size_t i = 0; i < seqs.size(); ++i) {
        auto it = aligned.find("s" + std::to_string(i));

        if (it == aligned.end() || it->second.empty())
            return false;

        if (width == 0)
            width = it->second.size();
        else if (it->second.size() != width)
            return false;

        rows.push_back(it->second);
    }

    return true;
}

struct MsaCached { bool ok; std::vector<std::string> rows; };
static std::unordered_map<std::string, MsaCached> g_msa_cache;
static std::mutex g_msa_cache_mx;

static std::string msa_key(const std::vector<std::string>& seqs) {
    std::string k;
    for (const auto& x : seqs) { k += x; k += '|'; }
    return k;
}

static bool msa_rows_cached(
    const std::vector<std::string>& seqs,
    const char* muscle_bin,
    std::vector<std::string>& rows)
{
    if (g_msa_jobs > 1) {
        std::lock_guard<std::mutex> lk(g_msa_cache_mx);
        auto it = g_msa_cache.find(msa_key(seqs));

        if (it != g_msa_cache.end()) {
            rows = it->second.rows;
            return it->second.ok;
        }
    }

    return msa_run_rows(seqs, muscle_bin, g_tool_threads, rows);
}

static std::vector<std::string> msa_members(
    const StarCluster& cluster,
    const std::map<std::string, int>& counts)
{
    std::vector<std::string> v;

    for (const std::string& m : cluster.members) {
        auto it = counts.find(m);

        if (it != counts.end() && it->second > 0)
            v.push_back(m);
    }

    return v;
}

static void msa_prefetch(
    const std::vector<StarCluster>& clusters,
    const std::map<std::string, int>& counts,
    const char* muscle_bin)
{
    if (g_msa_jobs <= 1)
        return;

    {
        std::lock_guard<std::mutex> lk(g_msa_cache_mx);
        g_msa_cache.clear();
    }

    std::vector<std::vector<std::string>> work;

    for (const StarCluster& sc : clusters) {
        std::vector<std::string> mem = msa_members(sc, counts);
        if (mem.size() >= 2)
            work.push_back(std::move(mem));
    }

    if (work.empty())
        return;

    std::atomic<size_t> next{0};
    const int n_threads = std::min<int>(g_msa_jobs, (int)work.size());
    std::vector<std::thread> pool;

    for (int t = 0; t < n_threads; ++t) {
        pool.emplace_back([&]() {
            for (size_t i = next++; i < work.size(); i = next++) {
                MsaCached c;

                c.ok = msa_run_rows(work[i], muscle_bin, 1, c.rows);

                std::lock_guard<std::mutex> lk(g_msa_cache_mx);
                g_msa_cache.emplace(msa_key(work[i]), std::move(c));
            }
        });
    }

    for (auto& th : pool)
        th.join();
}

static MsaConsensus msa_consensus(
    const StarCluster& cluster,
    const std::map<std::string, int>& counts,
    const std::map<std::string, StarVoteProfile>& profiles,
    const char* muscle_bin,
    bool want_weighted,
    bool allow_len_fix)
{
    MsaConsensus out;

    std::vector<std::pair<std::string, int>> mem;

    for (const std::string& m : cluster.members) {
        auto it = counts.find(m);

        if (it != counts.end() && it->second > 0)
            mem.push_back({m, it->second});
    }

    if (mem.empty())
        return out;

    if (mem.size() == 1) {
        out.plain = mem[0].first;
        out.weighted = mem[0].first;
        return out;
    }

    std::vector<std::string> rows;

    {
        std::vector<std::string> seqs;
        seqs.reserve(mem.size());
        for (const auto& m : mem) seqs.push_back(m.first);

        if (!msa_rows_cached(seqs, muscle_bin, rows))
            return out;
    }

    const size_t width = rows[0].size();

    std::vector<std::array<double, 5>> hard(width, {0, 0, 0, 0, 0});
    std::vector<std::array<double, 5>> full(width, {0, 0, 0, 0, 0});

    for (size_t i = 0; i < rows.size(); ++i) {
        std::vector<int> col_of_pos;
        col_of_pos.reserve(width);

        for (size_t c = 0; c < width; ++c) {
            const char ch = rows[i][c];

            if (ch != '-')
                col_of_pos.push_back((int)c);

            const int b = (ch == '-') ? 4 : dna_base_index(ch);

            if (b >= 0) {
                hard[c][b] += mem[i].second;
                full[c][b] += mem[i].second;
            }
        }

        if (!want_weighted)
            continue;

        auto pit = profiles.find(mem[i].first);

        if (pit == profiles.end())
            continue;

        for (const auto& kv : pit->second.delta_by_pos) {
            const int pos = kv.first;

            if (pos < 0 || pos >= (int)col_of_pos.size())
                continue;

            const int c = col_of_pos[pos];

            for (int b = 0; b < 4; ++b)
                full[c][b] += kv.second[b];
        }
    }

    const uint64_t tie_seed =
        (uint64_t)std::hash<std::string>{}(cluster.canonical);

    auto tie_pick =
        [tie_seed](size_t col, int n) -> int
    {
        uint64_t x =
            tie_seed ^ (0x9E3779B97F4A7C15ULL * (uint64_t)(col + 1));

        x ^= x >> 33;
        x *= 0xFF51AFD7ED558CCDULL;
        x ^= x >> 33;

        return (int)(x % (uint64_t)n);
    };

    struct WeakCol {
        int pos;
        char runner_up;
        double margin;
    };

    auto build =
        [&](const std::vector<std::array<double, 5>>& t,
            std::vector<WeakCol>* weak)
    {
        constexpr double EPS = 1e-12;

        std::vector<char> keep(width, 0);
        std::vector<int> best_base(width, 0);

        std::vector<double> gap_margin(width, 0.0);

        size_t kept = 0;

        for (size_t c = 0; c < width; ++c) {
            int tied[4];
            int n_tied = 0;
            double top = -1.0;

            for (int b = 0; b < 4; ++b) {
                if (t[c][b] > top + EPS) {
                    top = t[c][b];
                    n_tied = 0;
                    tied[n_tied++] = b;
                }
                else if (std::fabs(t[c][b] - top) <= EPS) {
                    tied[n_tied++] = b;
                }
            }

            const int bb =
                (n_tied <= 1) ? tied[0] : tied[tie_pick(c, n_tied)];

            best_base[c] = bb;
            gap_margin[c] = t[c][bb] - t[c][4];

            if (t[c][bb] <= EPS && t[c][4] <= EPS)
                continue;

            if (t[c][bb] + EPS >= t[c][4]) {
                keep[c] = 1;
                ++kept;
            }
        }

        const int delta = DNA_LENGTH - (int)kept;

        if (allow_len_fix && delta != 0) {
            std::vector<size_t> cand;

            for (size_t c = 0; c < width; ++c) {
                if (delta < 0 && keep[c])
                    cand.push_back(c);
                else if (delta > 0 && !keep[c] &&
                         t[c][best_base[c]] > EPS)
                    cand.push_back(c);
            }

            std::sort(
                cand.begin(),
                cand.end(),
                [&](size_t a, size_t b) {
                    const double ma =
                        (delta < 0) ? gap_margin[a] : -gap_margin[a];
                    const double mb =
                        (delta < 0) ? gap_margin[b] : -gap_margin[b];

                    if (ma != mb)
                        return ma < mb;

                    return a < b;
                });

            if ((int)cand.size() >= std::abs(delta)) {
                for (int i = 0; i < std::abs(delta); ++i)
                    keep[cand[i]] = (delta > 0) ? 1 : 0;
            }
        }

        std::string seq;
        seq.reserve(width);

        for (size_t c = 0; c < width; ++c) {
            if (!keep[c])
                continue;

            const int best = best_base[c];

            const int out_pos = (int)seq.size();

            seq.push_back(dna_index_base(best));

            if (!weak)
                continue;

            for (int b = 0; b < 4; ++b) {
                if (b == best || t[c][b] <= EPS)
                    continue;

                weak->push_back(
                    WeakCol{out_pos,
                            dna_index_base(b),
                            t[c][best] - t[c][b]});
            }
        }

        return seq;
    };

    std::vector<WeakCol> weak;

    out.plain = build(hard, nullptr);

    out.weighted =
        want_weighted
            ? build(full, &weak)
            : build(hard, &weak);

    if (!want_weighted)
        out.weighted = out.plain;

    std::sort(
        weak.begin(),
        weak.end(),
        [](const WeakCol& a, const WeakCol& b) {
            if (a.margin != b.margin)
                return a.margin < b.margin;
            return a.pos < b.pos;
        });

    if ((int)weak.size() > g_max_alt)
        weak.resize(g_max_alt);

    for (const WeakCol& w : weak) {
        out.weak_pos.push_back(w.pos);
        out.weak_runner_up.push_back(w.runner_up);
        out.weak_margin.push_back(w.margin);
    }

    return out;
}

struct StarConsensus {
    std::string seq;
    std::vector<int>    weak_pos;
    std::vector<char>   weak_runner_up;
    std::vector<double> weak_margin;
};

struct AltCombo {
    std::vector<int> idx;
    double score = 0.0;
};

static void alt_combos_recurse(
    int n,
    int max_flips,
    int start,
    std::vector<int>& cur,
    std::vector<AltCombo>& out)
{
    for (int i = start; i < n; ++i) {
        cur.push_back(i);
        out.push_back(AltCombo{cur, 0.0});

        if ((int)cur.size() < max_flips)
            alt_combos_recurse(n, max_flips, i + 1, cur, out);

        cur.pop_back();
    }
}

static std::vector<AltCombo> alt_combos(
    int n,
    int max_flips,
    const std::vector<double>& score,
    bool include_empty = false)
{
    std::vector<AltCombo> out;

    if (include_empty)
        out.push_back(AltCombo{{}, 0.0});

    if (n < 1 || max_flips < 1)
        return out;

    std::vector<int> cur;
    alt_combos_recurse(n, std::min(max_flips, n), 0, cur, out);

    for (AltCombo& c : out) {
        double s = 0.0;

        for (int i : c.idx)
            s += score[(size_t)i];

        c.score = s;
    }

    std::stable_sort(
        out.begin(),
        out.end(),
        [](const AltCombo& a, const AltCombo& b) {
            if (a.idx.size() != b.idx.size())
                return a.idx.size() < b.idx.size();

            if (a.score != b.score)
                return a.score < b.score;

            return a.idx < b.idx;
        });

    return out;
}

static bool star_repair_consensus(
    const StarConsensus& cons,
    const GF256& gf,
    std::string& out_seq,
    int& out_flips,
    bool try_unchanged = false)
{
    const int n = (int)cons.weak_pos.size();

    if (n < 1 && !try_unchanged)
        return false;

    std::vector<double> score(cons.weak_margin.begin(),
                              cons.weak_margin.end());

    const std::vector<AltCombo> combos =
        alt_combos(n, alt_flip_budget(), score, try_unchanged);

    for (const AltCombo& cb : combos) {
        bool clash = false;

        for (size_t a = 0; a + 1 < cb.idx.size() && !clash; ++a) {
            for (size_t b = a + 1; b < cb.idx.size(); ++b) {
                if (cons.weak_pos[cb.idx[a]] == cons.weak_pos[cb.idx[b]]) {
                    clash = true;
                    break;
                }
            }
        }

        if (clash)
            continue;

        std::string trial = cons.seq;

        for (int i : cb.idx)
            trial[cons.weak_pos[i]] = cons.weak_runner_up[i];

        uint8_t cw[MAX_CODE_BYTES];
        if (!dna_to_rs_bytes(trial, cw))
            continue;

        if (rs_detect_syndrome(cw, gf)) {
            out_seq = std::move(trial);
            out_flips = (int)cb.idx.size();
            return true;
        }
    }

    return false;
}

struct DecodeResult {
    int    result = 0;
    int    num_clusters = 0;
    int    trial_number = 0;
    int    LT_real_N = 0;
    int    LT_decoding_observance = 0;

    int    inferred_count = 0;
    double base_error = 0.0;
};

static DecodeResult run_decode(
    const PearFastqDB& r1_db,
    const PearFastqDB& r2_db,
    const PearFastqDB* merged_db,
    std::unordered_map<std::string, CachedPearMerge>& merge_cache,
    const PearMerger& pear_merger,
    const NgmergeMerger* ng_merger,
    const FlashMerger* fl_merger,
    const std::vector<std::string>& sampled_headers,
    int read_limit,
    int alt_mode,
    int cluster_mode,
    int starcode_enabled,
    int multilen,
    int ml_weighted,
    int star_dist,
    const char* msa_bin,
    const GF256& gf,
    const std::vector<Dataentry>& seed_data,
    const std::unordered_map<unsigned long, int>& seed_map,
    const int* Binary_Data_input,
    int* LT_parity_flat,
    int* Inferred_bitstream,
    int* Inferred_idx,
    std::vector<std::vector<int>>& node_to_strands,
    std::vector<std::vector<int>>& strand_nodes,
    int* Current_bitstream_flat)
{
    std::map<std::string, int> read_count;

    std::unordered_map<std::string, int> i1_mm;

    std::map<std::string, int> recovered_count;

    std::unordered_map<std::string, std::pair<int, int>> recovered_meta;

    std::map<std::string, int> final_rs_fail_count;

    std::map<std::string, int> all_reads_count;

    std::map<std::string, int> multilen_count;

    const bool want_weighted_vote =
        (multilen > 0) ? (ml_weighted != 0) : (alt_mode != ALT_OFF);

    std::map<std::string, StarVoteProfile> star_vote_profiles;

    const int sampled_total =
        static_cast<int>(sampled_headers.size());

    const int selected_count =
        std::min(read_limit, sampled_total);

    auto make_entry =
        [&](const std::string& header,
            const PearMerger& pm,
            const NgmergeMerger* ngm,
            const FlashMerger* flm) -> CachedPearMerge
    {
        CachedPearMerge entry;

        if (merged_db) {
            auto mit = merged_db->find(header);

            if (mit == merged_db->end()) {
                entry.status = 1;
            }
            else {
                entry.merged.assembled = true;
                entry.merged.seq = mit->second.seq;
                entry.merged.qual = mit->second.qual;
                entry.status = 2;
            }
        }
        else {
            auto r1_it = r1_db.find(header);
            auto r2_it = r2_db.find(header);

            if (r1_it == r1_db.end() ||
                r2_it == r2_db.end())
            {
                entry.status = 0;
            }
            else {
                if (ngm)
                    entry.merged = ngm->merge(r1_it->second, r2_it->second);
                else if (flm)
                    entry.merged = flm->merge(r1_it->second, r2_it->second);
                else
                    entry.merged = pm.merge(r1_it->second, r2_it->second);

                entry.status = entry.merged.assembled ? 2 : 1;
            }
        }

        return entry;
    };

    if (g_merge_jobs > 1) {
        std::vector<const std::string*> todo;
        std::unordered_set<std::string> queued;

        for (int i = 0; i < selected_count; ++i) {
            const std::string& h = sampled_headers[i];

            if (!merge_cache.count(h) && queued.insert(h).second)
                todo.push_back(&h);
        }

        if (!todo.empty()) {
            std::vector<CachedPearMerge> built(todo.size());
            const int n_threads =
                std::min<int>(g_merge_jobs, (int)todo.size());

            std::atomic<size_t> next{0};
            std::vector<std::thread> pool;

            for (int t = 0; t < n_threads; ++t) {
                pool.emplace_back([&]() {
                    PearMerger pm_local = pear_merger;
                    std::unique_ptr<NgmergeMerger> ng_local;
                    std::unique_ptr<FlashMerger> fl_local;

                    if (ng_merger) ng_local.reset(new NgmergeMerger(*ng_merger));
                    if (fl_merger) fl_local.reset(new FlashMerger(*fl_merger));

                    for (size_t k = next++; k < todo.size(); k = next++)
                        built[k] = make_entry(*todo[k], pm_local,
                                              ng_local.get(), fl_local.get());
                });
            }

            for (auto& th : pool)
                th.join();

            for (size_t k = 0; k < todo.size(); ++k)
                merge_cache.emplace(*todo[k], std::move(built[k]));
        }
    }

    for (int i = 0; i < selected_count; ++i) {
        const std::string& header = sampled_headers[i];

        CachedPearMerge* cached = nullptr;

        auto cache_it = merge_cache.find(header);

        if (cache_it != merge_cache.end()) {
            cached = &cache_it->second;
        }
        else {
            auto inserted =
                merge_cache.emplace(
                    header,
                    make_entry(header, pear_merger, ng_merger, fl_merger));

            cached = &inserted.first->second;
        }

        if (cached->status != 2)
            continue;

        const std::string& seq = cached->merged.seq;

        if ((int)seq.size() != DNA_LENGTH) {
            if (multilen > 0 &&
                std::abs((int)seq.size() - DNA_LENGTH) <= multilen)
            {
                multilen_count[seq]++;

                add_star_vote_profile(
                    star_vote_profiles, seq, cached->merged.mismatches,
                    alt_mode);
            }

            continue;
        }

        uint8_t codeword[MAX_CODE_BYTES];

        if (!dna_to_rs_bytes(seq, codeword))
            continue;

        const bool rs_ok = rs_detect_syndrome(codeword, gf);

        if (cluster_mode == CLUSTER_STARCODE) {
            all_reads_count[seq]++;

            add_star_vote_profile(
                star_vote_profiles, seq, cached->merged.mismatches,
                alt_mode);

            continue;
        }

        auto keep_rs_fail = [&]() {
            if (starcode_enabled) {
                final_rs_fail_count[seq]++;
                add_star_vote_profile(
                    star_vote_profiles, seq, cached->merged.mismatches,
                    alt_mode);
            }
        };

        if (!rs_ok) {
            if (alt_mode == ALT_OFF) {
                keep_rs_fail();
                continue;
            }

            const std::vector<PearAlternativeBase>& all_mm =
                cached->merged.mismatches;

            if (all_mm.empty()) {
                keep_rs_fail();
                continue;
            }

            std::vector<PearAlternativeBase> eligible;

            if (alt_mode == ALT_EQUAL_Q) {
                for (const PearAlternativeBase& m : all_mm)
                    if (m.q_gap() <= 0)
                        eligible.push_back(m);
            }
            else {
                eligible = all_mm;
            }

            if (eligible.empty()) {
                keep_rs_fail();
                continue;
            }

            if ((int)eligible.size() > g_max_alt) {
                std::stable_sort(
                    eligible.begin(),
                    eligible.end(),
                    alt_rank_less);

                eligible.resize(g_max_alt);
            }

            const int n_alt = (int)eligible.size();

            if (n_alt < 1) {
                keep_rs_fail();
                continue;
            }

            std::vector<double> alt_score(n_alt);

            for (int a = 0; a < n_alt; ++a)
                alt_score[a] = (double)eligible[a].q_gap();

            const std::vector<AltCombo> combos =
                alt_combos(n_alt, alt_flip_budget(), alt_score);

            bool recovered = false;
            int winning_flips = 0;
            int winning_qsum = 0;
            std::string recovered_seq;

            for (const AltCombo& cb : combos) {
                int cand_qsum = 0;
                for (int a : cb.idx)
                    cand_qsum += eligible[(size_t)a].q_gap();

                if (g_max_qsum >= 0 && cand_qsum > g_max_qsum)
                    continue;

                std::string candidate_seq = seq;
                bool valid_candidate = true;

                for (int a : cb.idx) {
                    const PearAlternativeBase& alt = eligible[(size_t)a];

                    if (alt.position < 0 ||
                        alt.position >= (int)candidate_seq.size())
                    {
                        valid_candidate = false;
                        break;
                    }

                    if (candidate_seq[alt.position] != alt.selected_base) {
                        valid_candidate = false;
                        break;
                    }

                    candidate_seq[alt.position] = alt.alternative_base;
                }

                if (!valid_candidate)
                    continue;

                uint8_t recovered_codeword[MAX_CODE_BYTES];

                if (!dna_to_rs_bytes(candidate_seq, recovered_codeword))
                    continue;

                if (rs_detect_syndrome(recovered_codeword, gf)) {
                    recovered = true;
                    winning_flips = (int)cb.idx.size();
                    winning_qsum = cand_qsum;
                    recovered_seq = std::move(candidate_seq);
                    break;
                }
            }

            if (recovered) {
                recovered_count[recovered_seq]++;

                const std::pair<int, int> m{winning_flips, winning_qsum};
                auto mit = recovered_meta.find(recovered_seq);

                if (mit == recovered_meta.end())
                    recovered_meta.emplace(recovered_seq, m);
                else if (m < mit->second)
                    mit->second = m;

                if (starcode_enabled) {
                    add_star_vote_profile(
                        star_vote_profiles, recovered_seq, all_mm,
                        alt_mode);
                }
            }
            else {
                keep_rs_fail();
            }

            continue;
        }

        read_count[seq]++;
        {
            const int nmm = (int)cached->merged.mismatches.size();
            auto ins = i1_mm.emplace(seq, nmm);
            if (!ins.second && nmm < ins.first->second)
                ins.first->second = nmm;
        }

        if (starcode_enabled) {
            add_star_vote_profile(
                star_vote_profiles, seq, cached->merged.mismatches,
                alt_mode);
        }
    }

    std::vector<std::pair<std::string, int>> rec_orphans;

    for (const auto& kv : recovered_count) {
        auto it = read_count.find(kv.first);

        if (it != read_count.end())
            it->second += kv.second;
        else
            rec_orphans.push_back(kv);
    }

    if (!rec_orphans.empty()) {
        const int seed_nt = LT_SEED_BIT / 2;

        std::unordered_set<std::string> seed_i1;

        for (const auto& kv : read_count) {
            const std::string& q = kv.first;

            if ((int)q.size() >= seed_nt)
                seed_i1.insert(q.substr(0, seed_nt));
        }

        std::vector<std::pair<std::string, int>> kept;
        kept.reserve(rec_orphans.size());

        for (const auto& kv : rec_orphans) {
            const std::string& o = kv.first;

            const bool conflict =
                (int)o.size() >= seed_nt &&
                seed_i1.count(o.substr(0, seed_nt));

            if (!conflict)
                kept.push_back(kv);
        }

        rec_orphans.swap(kept);
    }

    auto by_count_then_seq =
        [](const auto& a, const auto& b) {
            if (a.second != b.second)
                return a.second > b.second;

            return a.first < b.first;
        };

    auto cluster_consensus =
        [&](const StarCluster& sc,
            const std::map<std::string, int>& counts,
            std::string& out) -> int
    {
        const MsaConsensus mc =
            msa_consensus(
                sc, counts, star_vote_profiles,
                msa_bin, want_weighted_vote,
                alt_mode != ALT_OFF);

        if ((int)mc.plain.size() == DNA_LENGTH) {
            uint8_t cw[MAX_CODE_BYTES];

            if (dna_to_rs_bytes(mc.plain, cw) &&
                rs_detect_syndrome(cw, gf))
            {
                out = mc.plain;
                return 1;
            }
        }

        const bool try_weighted = (mc.weighted != mc.plain);

        if (alt_mode != ALT_OFF &&
            (int)mc.weighted.size() == DNA_LENGTH &&
            (!mc.weak_pos.empty() || try_weighted))
        {
            StarConsensus cons;
            cons.seq = mc.weighted;
            cons.weak_pos = mc.weak_pos;
            cons.weak_runner_up = mc.weak_runner_up;
            cons.weak_margin = mc.weak_margin;

            int flips = 0;

            if (star_repair_consensus(cons, gf, out, flips, try_weighted))
                return (flips == 0) ? 2 : 3;
        }

        return 0;
    };

    std::vector<std::pair<std::string, int>> clusters;

    std::vector<signed char> cluster_tier;

    int tier1_n = 0;
    int tier2_n = 0;
    int tier3_n = 0;

    if (cluster_mode == CLUSTER_STARCODE) {
        for (const auto& kv : multilen_count)
            all_reads_count[kv.first] += kv.second;

        std::vector<StarCluster> sclusters;

        if (!run_clusters(all_reads_count, sclusters, star_dist)) {
            fprintf(stderr,
                    "Error: primary clustering failed; this round\n"
                    "       produces no strand.\n");
        }
        else {
            std::vector<std::pair<std::string, int>> tier_cons[4];

            msa_prefetch(sclusters, all_reads_count, msa_bin);

            for (const StarCluster& sc : sclusters) {
                std::string cons;
                const int t = cluster_consensus(sc, all_reads_count, cons);

                if (t > 0)
                    tier_cons[t].push_back({cons, sc.cluster_size});
            }

            for (int t = 1; t <= 3; ++t)
                std::sort(tier_cons[t].begin(), tier_cons[t].end(),
                          by_count_then_seq);

            tier1_n = (int)tier_cons[1].size();
            tier2_n = (int)tier_cons[2].size();
            tier3_n = (int)tier_cons[3].size();

            clusters = std::move(tier_cons[1]);
            clusters.insert(clusters.end(), tier_cons[2].begin(), tier_cons[2].end());
            clusters.insert(clusters.end(), tier_cons[3].begin(), tier_cons[3].end());
        }
    }
    else {
        std::vector<std::pair<std::string, int>> star_rescued;

        if (starcode_enabled) {
            std::map<std::string, int> star_counts;

            for (const auto& kv : final_rs_fail_count)
                star_counts[kv.first] += kv.second;

            for (const auto& kv : multilen_count)
                star_counts[kv.first] += kv.second;

            std::vector<std::pair<std::string, int>> rescue_cons[4];
            std::vector<StarCluster> star_clusters;

            if (!run_clusters(star_counts, star_clusters, star_dist)) {
                fprintf(stderr,
                        "Warning: cluster rescue skipped; continuing with "
                        "the original decoding path.\n");
            }
            else {
                msa_prefetch(star_clusters, star_counts, msa_bin);

                for (const StarCluster& sc : star_clusters) {
                    if (sc.members.size() < 2)
                        continue;

                    std::string cons;
                    const int t = cluster_consensus(sc, star_counts, cons);

                    if (t > 0)
                        rescue_cons[t].push_back({cons, sc.cluster_size});
                }
            }

            for (int t = 1; t <= 3; ++t)
                std::sort(rescue_cons[t].begin(), rescue_cons[t].end(),
                          by_count_then_seq);

            star_rescued = std::move(rescue_cons[1]);
            star_rescued.insert(star_rescued.end(),
                                rescue_cons[2].begin(), rescue_cons[2].end());
            star_rescued.insert(star_rescued.end(),
                                rescue_cons[3].begin(), rescue_cons[3].end());
        }

        clusters.assign(read_count.begin(), read_count.end());

        tier1_n = (int)clusters.size();

        if (g_i1_mm_order) {
            auto mm_of = [&](const std::string& q) {
                auto it = i1_mm.find(q);
                return (it != i1_mm.end()) ? it->second : INT_MAX;
            };

            std::sort(
                clusters.begin(),
                clusters.end(),
                [&](const std::pair<std::string, int>& a,
                    const std::pair<std::string, int>& b)
                {
                    if (a.second != b.second)
                        return a.second > b.second;
                    if (a.second == 1) {
                        const int ma = mm_of(a.first), mb = mm_of(b.first);
                        if (ma != mb)
                            return ma < mb;
                    }
                    return a.first < b.first;
                });
        }
        else {
            std::sort(clusters.begin(), clusters.end(), by_count_then_seq);
        }

        cluster_tier.assign(clusters.size(), 1);

        std::unordered_set<std::string> already_present;
        already_present.reserve(clusters.size() * 2 + rec_orphans.size() +
                                star_rescued.size() + 1);

        for (const auto& kv : clusters)
            already_present.insert(kv.first);

        for (const auto& kv : rec_orphans)
            already_present.insert(kv.first);

        if (starcode_enabled) {
            for (const auto& kv : star_rescued) {
                if (already_present.insert(kv.first).second) {
                    clusters.push_back(kv);
                    cluster_tier.push_back(3);
                    ++tier3_n;
                }
            }
        }

        auto meta_of = [&](const std::string& q) {
            auto it = recovered_meta.find(q);
            return (it != recovered_meta.end())
                ? it->second
                : std::pair<int, int>{MAX_ALT_CAP + 1, 1 << 20};
        };

        std::sort(
            rec_orphans.begin(),
            rec_orphans.end(),
            [&](const std::pair<std::string, int>& a,
                const std::pair<std::string, int>& b)
            {
                const auto ma = meta_of(a.first);
                const auto mb = meta_of(b.first);

                if (ma != mb)
                    return ma < mb;

                return by_count_then_seq(a, b);
            });

        for (const auto& kv : rec_orphans) {
            clusters.push_back(kv);
            cluster_tier.push_back(2);
            ++tier2_n;
        }
    }

    if (tier3_n > 0) {
        const int seed_nt = LT_SEED_BIT / 2;
        const bool have_tier = (cluster_tier.size() == clusters.size());

        auto tier_of = [&](size_t i) -> int {
            return have_tier ? cluster_tier[i]
                 : ((int)i < tier1_n) ? 1
                 : ((int)i < tier1_n + tier2_n) ? 2 : 3;
        };

        std::unordered_set<std::string> i1_seq, i1_seed;

        for (size_t i = 0; i < clusters.size(); ++i) {
            if (tier_of(i) != 1) continue;
            const std::string& q = clusters[i].first;
            i1_seq.insert(q);
            if ((int)q.size() >= seed_nt) i1_seed.insert(q.substr(0, seed_nt));
        }

        std::vector<std::pair<std::string, int>> kept;
        std::vector<signed char> kept_tier;
        kept.reserve(clusters.size());
        kept_tier.reserve(clusters.size());

        int i3_conflict = 0;

        for (size_t i = 0; i < clusters.size(); ++i) {
            const int t = tier_of(i);
            const std::string& o = clusters[i].first;

            const bool conflict =
                t == 3 && !i1_seq.count(o) &&
                (int)o.size() >= seed_nt &&
                i1_seed.count(o.substr(0, seed_nt));

            if (conflict) {
                ++i3_conflict;
            }
            else {
                kept.push_back(clusters[i]);
                kept_tier.push_back((signed char)t);
            }
        }

        if (i3_conflict > 0) {
            clusters.swap(kept);
            if (have_tier) cluster_tier.swap(kept_tier);
            tier3_n -= i3_conflict;
        }
    }

    if (!clusters.empty()) {
        const int seed_nt = LT_SEED_BIT / 2;
        const bool have_tier = (cluster_tier.size() == clusters.size());

        auto seed_of = [&](const std::string& q, unsigned long& g) {
            if ((int)q.size() < seed_nt) return false;
            g = 0;
            for (int i = 0; i < seed_nt; ++i) {
                unsigned long v;
                switch (q[i]) {
                    case 'A': v = 0; break;
                    case 'C': v = 1; break;
                    case 'G': v = 2; break;
                    case 'T': v = 3; break;
                    default:  return false;
                }
                g = (g << 2) | v;
            }
            return true;
        };

        std::vector<std::pair<std::string, int>> kept;
        std::vector<signed char> kept_tier;
        kept.reserve(clusters.size());
        kept_tier.reserve(clusters.size());

        int sc_drop[4] = {0, 0, 0, 0};

        for (size_t i = 0; i < clusters.size(); ++i) {
            const int t = have_tier
                ? cluster_tier[i]
                : ((int)i < tier1_n) ? 1
                : ((int)i < tier1_n + tier2_n) ? 2 : 3;

            unsigned long g = 0;
            const bool ok = seed_of(clusters[i].first, g) && seed_map.count(g);

            if (ok) {
                kept.push_back(clusters[i]);
                kept_tier.push_back((signed char)t);
            }
            else {
                ++sc_drop[t >= 1 && t <= 3 ? t : 3];
            }
        }

        if (sc_drop[1] + sc_drop[2] + sc_drop[3] > 0) {
            clusters.swap(kept);
            if (have_tier) cluster_tier.swap(kept_tier);
            tier1_n -= sc_drop[1];
            tier2_n -= sc_drop[2];
            tier3_n -= sc_drop[3];
        }
    }

    DecodeResult dr{};

    dr.num_clusters = (int)clusters.size();

    if (clusters.empty()) {
        dr.result = 0;
        dr.inferred_count = 0;
        return dr;
    }

    std::string acgt_concat;
    acgt_concat.reserve((size_t)clusters.size() * DNA_LENGTH);

    for (const auto& cluster : clusters)
        acgt_concat += cluster.first;

    std::vector<int> raw_bits(acgt_concat.size() * 2);

    int total_bits =
        acgt_to_bits(
            acgt_concat.c_str(),
            (int)acgt_concat.size(),
            raw_bits.data());

    int number_of_strands = total_bits / STRAND_WIDTH;

    if (number_of_strands > MAX_STRANDS)
        number_of_strands = MAX_STRANDS;

    memcpy(
        LT_parity_flat,
        raw_bits.data(),
        (long long)number_of_strands *
            STRAND_WIDTH *
            sizeof(int));

    int LT_decoding_observance = 0;
    int LT_real_N              = 0;
    int trial_number           = 0;
    int inferred_count         = 0;

    std::queue<std::pair<int, int>> degree1_queue;

    auto process_degree1 =
        [&](int node_idx, int* seq)
        {
            if (Inferred_idx[node_idx] != 0)
                return;

            memcpy(
                Inferred_bitstream + (long long)node_idx * LT_L,
                seq,
                LT_L * sizeof(int));

            Inferred_idx[node_idx] = 1;
            inferred_count++;

            for (int strand_idx : node_to_strands[node_idx]) {
                auto& nodes = strand_nodes[strand_idx];

                for (int k = 0; k < (int)nodes.size(); k++) {
                    if (nodes[k] != node_idx)
                        continue;

                    int* current =
                        Current_bitstream_flat +
                        (long long)strand_idx * LT_L;

                    int* inferred =
                        Inferred_bitstream +
                        (long long)node_idx * LT_L;

                    for (int b = 0; b < LT_L; b++)
                        current[b] =
                            (current[b] + inferred[b]) & 1;

                    nodes[k] = nodes.back();
                    nodes.pop_back();
                    break;
                }

                if (nodes.size() == 1)
                    degree1_queue.push(
                        {nodes[0], strand_idx});
            }
        };

    std::vector<const std::vector<int>*> strand_list(number_of_strands, nullptr);

    for (int s = 0; s < number_of_strands; ++s) {
        const int* cand = LT_parity_flat + (long long)s * STRAND_WIDTH;
        unsigned long g = 0;

        for (int i = 0; i < LT_SEED_BIT; i++)
            if (cand[i] == 1)
                g |= (1UL << (31 - i));

        auto it = seed_map.find(g);

        if (it != seed_map.end())
            strand_list[s] = &seed_data[it->second].indices;
    }

    memset(Inferred_bitstream, 0,
           (long long)LT_K * LT_L * sizeof(int));
    memset(Inferred_idx, 0, LT_K * sizeof(int));

    for (int i = 0; i < LT_K; i++)
        node_to_strands[i].clear();

    for (int i = 0; i < number_of_strands; i++)
        strand_nodes[i].clear();

    while (trial_number < number_of_strands &&
           inferred_count != LT_K)
    {
        while (!degree1_queue.empty() &&
               inferred_count != LT_K)
        {
            auto [node_idx, strand_idx] = degree1_queue.front();
            degree1_queue.pop();

            process_degree1(
                node_idx,
                Current_bitstream_flat + (long long)strand_idx * LT_L);
        }

        if (inferred_count == LT_K)
            break;

        int* binary_candidate =
            LT_parity_flat + (long long)trial_number * STRAND_WIDTH;

        const int cur_strand = trial_number;
        trial_number++;

        const std::vector<int>* nlp = strand_list[(size_t)cur_strand];

        if (!nlp)
            continue;

        const std::vector<int>& node_list = *nlp;

        bool has_uninferred = false;

        for (int node : node_list) {
            if (Inferred_idx[node] == 0) {
                has_uninferred = true;
                break;
            }
        }

        if (!has_uninferred)
            continue;

        LT_decoding_observance++;

        std::vector<int> sub_seq(LT_L);

        memcpy(sub_seq.data(),
               binary_candidate + LT_SEED_BIT,
               LT_L * sizeof(int));

        std::vector<int> remaining;

        for (int node : node_list) {
            if (Inferred_idx[node] == 1) {
                int* inferred =
                    Inferred_bitstream + (long long)node * LT_L;

                for (int b = 0; b < LT_L; b++)
                    sub_seq[b] = (sub_seq[b] + inferred[b]) & 1;
            }
            else {
                remaining.push_back(node);
            }
        }

        if (remaining.size() == 1) {
            process_degree1(remaining[0], sub_seq.data());
        }
        else {
            if (LT_real_N >= MAX_STRANDS) {
                fprintf(stderr,
                        "Warning: LT_real_N reached MAX_STRANDS=%d\n",
                        MAX_STRANDS);
                break;
            }

            int strand_idx = LT_real_N;

            memcpy(Current_bitstream_flat + (long long)strand_idx * LT_L,
                   sub_seq.data(),
                   LT_L * sizeof(int));

            strand_nodes[strand_idx] = remaining;

            for (int node : remaining)
                node_to_strands[node].push_back(strand_idx);

            LT_real_N++;
        }
    }

    while (!degree1_queue.empty() &&
           inferred_count != LT_K)
    {
        auto [node_idx, strand_idx] = degree1_queue.front();
        degree1_queue.pop();

        process_degree1(
            node_idx,
            Current_bitstream_flat + (long long)strand_idx * LT_L);
    }

    if (inferred_count != LT_K) {
        dr.result = 0;
    }
    else {
        bool equal =
            memcmp(
                Binary_Data_input,
                Inferred_bitstream,
                (long long)LT_K *
                    LT_L *
                    sizeof(int)) == 0;

        dr.result = equal ? 1 : 2;
    }

    int bit_errors = 0;

    for (int i = 0; i < LT_K * LT_L; i++) {
        if (Binary_Data_input[i] !=
            Inferred_bitstream[i]) {
            bit_errors++;
        }
    }

    dr.trial_number           = trial_number;
    dr.LT_real_N              = LT_real_N;
    dr.LT_decoding_observance = LT_decoding_observance;
    dr.inferred_count         = inferred_count;

    dr.base_error =
        (double)bit_errors / (LT_K * LT_L);

    return dr;
}

static int decode_dna(
    int start_index,
    int end_index,
    int rsn,
    int alt_mode,
    int cluster_mode,
    int starcode_enabled,
    int merger_mode,
    const char* merger_name,
    int multilen,
    int ml_weighted,
    int star_dist,
    const char* msa_bin)
{
    clock_t tStart = clock();

    const char* image_path =
        g_image_opt ? g_image_opt : DEFAULT_IMAGE_FILE;

    const char* seed_table_path =
        g_seed_table_opt ? g_seed_table_opt : DEFAULT_SEED_TABLE_FILE;

    char base_dir[512];

    if (g_base_dir_opt)
        snprintf(base_dir, sizeof(base_dir), "%s", g_base_dir_opt);
    else
        snprintf(base_dir, sizeof(base_dir), "%s", DATA_DIR);

    FILE* fp_img = fopen(image_path, "r");

    if (!fp_img) {
        fprintf(
            stderr,
            "Error: Cannot open image: %s\n",
            image_path);
        return 1;
    }

    const int img_total = LT_K * LT_L;

    int* Binary_Data_input =
        (int*)malloc(
            (long long)img_total * sizeof(int));

    if (!Binary_Data_input) {
        fprintf(stderr, "Error: Cannot allocate image buffer\n");
        fclose(fp_img);
        return 1;
    }

    int img_loaded = 0;

    while (img_loaded < img_total &&
           fscanf(
               fp_img,
               "%1d",
               &Binary_Data_input[img_loaded]) == 1)
    {
        img_loaded++;
    }

    fclose(fp_img);

    if (img_loaded != img_total) {
        fprintf(
            stderr,
            "Error: image bits = %d, expected = %d\n",
            img_loaded,
            img_total);

        free(Binary_Data_input);
        return 1;
    }

    std::vector<Dataentry> seed_data;
    std::unordered_map<unsigned long, int> seed_map;

    {
        std::ifstream fin(seed_table_path);

        if (!fin) {
            fprintf(
                stderr,
                "Error: Cannot open seed table: %s\n",
                seed_table_path);

            free(Binary_Data_input);
            return 1;
        }

        std::string sline;

        while (std::getline(fin, sline)) {
            if (sline.empty() || sline[0] == '#')
                continue;

            std::istringstream is(sline);
            Dataentry e;
            unsigned long seed_v;

            if (!(is >> seed_v))
                continue;

            e.seed = seed_v;
            long value;

            while (is >> value)
                if (value > 0)
                    e.indices.push_back((int)(value - 1));

            seed_map[e.seed] = (int)seed_data.size();
            seed_data.push_back(std::move(e));
        }
    }

    printf(
        "Seed table loaded: %d entries\n",
        (int)seed_data.size());

    char r1_path[640];
    char r2_path[640];

    if (g_r1_path_opt)
        snprintf(r1_path, sizeof(r1_path), "%s", g_r1_path_opt);
    else
        snprintf(
            r1_path,
            sizeof(r1_path),
            "%s/%s",
            base_dir, R1_FASTQ_NAME);

    if (g_r2_path_opt)
        snprintf(r2_path, sizeof(r2_path), "%s", g_r2_path_opt);
    else
        snprintf(
            r2_path,
            sizeof(r2_path),
            "%s/%s",
            base_dir, R2_FASTQ_NAME);

    printf("image      : %s\n", image_path);
    printf("seed_table : %s\n", seed_table_path);
    printf("base_dir   : %s\n", base_dir);
    printf("R1         : %s\n", r1_path);
    printf("R2         : %s\n\n", r2_path);

    std::vector<std::string> header_pool =
        scan_fastq_headers(r1_path);

    if (header_pool.empty()) {
        fprintf(
            stderr,
            "Error: No header was read from R1: %s\n",
            r1_path);

        free(Binary_Data_input);
        return 1;
    }

    printf(
        "R1 reads in FASTQ  : %zu\n",
        header_pool.size());

    size_t cap = (size_t)rsn;

    if (cap > header_pool.size()) {
        fprintf(
            stderr,
            "Warning: RSN (%d) > reads in R1 (%zu); using all %zu reads.\n",
            rsn,
            header_pool.size(),
            header_pool.size());

        cap = header_pool.size();
    }

    const int read_limit = (int)cap;

    const int total_frames =
        end_index - start_index + 1;

    std::vector<std::vector<std::string>>
        sampled_by_idx(total_frames);

    std::unordered_set<std::string>
        wanted_headers;

    for (int recursive = start_index;
         recursive <= end_index;
         ++recursive)
    {
        const int col = recursive - start_index;

        sampled_by_idx[col] =
            sample_headers(
                header_pool,
                cap,
                (unsigned long long)recursive);

        for (const auto& h : sampled_by_idx[col])
            wanted_headers.insert(h);
    }

    header_pool.clear();
    header_pool.shrink_to_fit();

    printf(
        "read_limit         : %d\n",
        read_limit);

    printf(
        "Unique sampled headers needed: %zu\n",
        wanted_headers.size());

    PearFastqDB merged_db;
    PearFastqDB* merged_ptr = nullptr;

    if (merger_mode == MERGER_EXTERNAL) {
        char merger_path[720];

        snprintf(
            merger_path,
            sizeof(merger_path),
            "%s/merged_%s.fastq",
            base_dir,
            merger_name);

        printf("merged FASTQ : %s\n", merger_path);

        size_t merged_total = 0;
        size_t merged_kept = 0;

        merged_db =
            load_merged_selected(
                merger_path,
                wanted_headers,
                merged_total,
                merged_kept);

        printf(
            "Merged records scanned : %zu\n"
            "Merged records kept    : %zu / %zu sampled\n",
            merged_total,
            merged_kept,
            wanted_headers.size());

        if (merged_db.empty()) {
            fprintf(
                stderr,
                "\nError: no sampled header matched the merged FASTQ.\n"
                "  %zu records were scanned. Keys are normalized with\n"
                "  pear_normalize_header() (trimmed, cut at the first\n"
                "  whitespace, leading '@' kept). Compare:  head -1 %s\n",
                merged_total,
                merger_path);

            free(Binary_Data_input);
            return 1;
        }

        merged_ptr = &merged_db;
    }

    PearPairedLoadResult paired;

    if (merger_mode != MERGER_EXTERNAL) {
        paired =
            pear_load_paired_selected(
                r1_path,
                r2_path,
                wanted_headers);

        printf(
            "FASTQ pairs scanned : %zu\n",
            paired.total_pairs);

        printf(
            "Selected pairs kept : %zu\n",
            paired.selected_pairs);

        printf(
            "Header mismatches   : %zu\n",
            paired.header_mismatch);

        if (paired.r1.empty() ||
            paired.r2.empty())
        {
            fprintf(
                stderr,
                "\nError: No selected R1/R2 pairs were loaded.\n");

            if (paired.total_pairs == 0) {
                fprintf(
                    stderr,
                    "  The FASTQ itself yielded no pair. Check that both\n"
                    "  files exist and are readable:\n"
                    "    %s\n"
                    "    %s\n",
                    r1_path,
                    r2_path);
            }
            else {
                fprintf(
                    stderr,
                    "  %zu pairs were scanned but none matched the sampled\n"
                    "  headers: R1 and R2 are not in the same order, or the\n"
                    "  two files belong to different runs.\n",
                    paired.total_pairs);
            }

            free(Binary_Data_input);
            return 1;
        }

        if (paired.selected_pairs * 2 < wanted_headers.size()) {
            fprintf(
                stderr,
                "Warning: only %zu of %zu sampled headers matched.\n"
                "         header_mismatch=%zu suggests R1/R2 ordering.\n",
                paired.selected_pairs,
                wanted_headers.size(),
                paired.header_mismatch);
        }

        printf("\n");
    }

    PearConfig pear_cfg;
    pear_cfg.p_value      = g_pear_pvalue;
    pear_cfg.min_overlap  = 10;
    pear_cfg.min_asm_len  = 50;
    pear_cfg.max_asm_len  = 0;
    pear_cfg.max_uncalled = 1.0;
    pear_cfg.score_method = 2;
    pear_cfg.phred_base   = 33;
    pear_cfg.nbase        = false;
    pear_cfg.test_method  = 1;
    pear_cfg.cap          = 40;

    PearMerger pear_merger(
        paired.ef,
        pear_cfg);

    NgmergeConfig ng_cfg;
    ng_cfg.max_mismatch = g_ngmerge_mismatch;

    NgmergeMerger ng_merger(ng_cfg);

    if (merger_mode == MERGER_NGMERGE &&
        g_ngmerge_profile && g_ngmerge_profile[0])
    {
        if (ng_merger.load_qual_profile(g_ngmerge_profile))
            printf("NGmerge quality profile: %s\n", g_ngmerge_profile);
        else
            printf("NGmerge quality profile: fastq-join fallback\n");
    }

    const NgmergeMerger* ng_ptr =
        (merger_mode == MERGER_NGMERGE) ? &ng_merger : nullptr;

    FlashConfig fl_cfg;
    fl_cfg.max_overlap = g_flash_max_overlap;
    fl_cfg.max_mismatch = g_flash_mismatch;

    FlashMerger fl_merger(fl_cfg);

    const FlashMerger* fl_ptr =
        (merger_mode == MERGER_FLASH) ? &fl_merger : nullptr;

    std::unordered_map<std::string, CachedPearMerge>
        merge_cache;

    merge_cache.reserve(
        wanted_headers.size() * 2);

    wanted_headers.clear();
    wanted_headers.rehash(0);

    GF256 gf;

    int* LT_parity_flat =
        (int*)malloc(
            (long long)MAX_STRANDS *
            STRAND_WIDTH *
            sizeof(int));

    int* Inferred_bitstream =
        (int*)malloc(
            (long long)LT_K *
            LT_L *
            sizeof(int));

    int* Inferred_idx =
        (int*)calloc(
            LT_K,
            sizeof(int));

    int* Current_bitstream_flat =
        (int*)malloc(
            (long long)MAX_STRANDS *
            LT_L *
            sizeof(int));

    if (!LT_parity_flat ||
        !Inferred_bitstream ||
        !Inferred_idx ||
        !Current_bitstream_flat)
    {
        fprintf(
            stderr,
            "Error: Cannot allocate LT decoder buffers\n");

        free(Binary_Data_input);
        free(LT_parity_flat);
        free(Inferred_bitstream);
        free(Inferred_idx);
        free(Current_bitstream_flat);

        return 1;
    }

    std::vector<std::vector<int>>
        node_to_strands(LT_K);

    std::vector<std::vector<int>>
        strand_nodes(MAX_STRANDS);

    char suffix[256];

    if (alt_mode == ALT_OFF)
        snprintf(suffix, sizeof(suffix), "org");
    else if (alt_mode == ALT_EQUAL_Q)
        snprintf(suffix, sizeof(suffix), "eqq");
    else
        snprintf(suffix, sizeof(suffix), "allmm");

    auto add = [&](const char* part) {
        strncat(suffix, part, sizeof(suffix) - strlen(suffix) - 1);
    };

    char part[64];

    if (cluster_mode == CLUSTER_STARCODE)
        add("_sc");

    if (starcode_enabled)
        add("_rescue");

    if (merger_mode != MERGER_PEAR) {
        add("_");
        add(merger_name);
    }

    if (g_cluster_tool == CLUSTER_TOOL_STARCODE &&
        (starcode_enabled || cluster_mode == CLUSTER_STARCODE))
    {
        snprintf(part, sizeof(part), "_d%d", star_dist);
        add(part);
    }

    if (multilen > 0) {
        snprintf(part, sizeof(part), "_ml%d%s", multilen,
                 ml_weighted ? "w" : "");
        add(part);
    }

    if (g_cluster_tool == CLUSTER_TOOL_CDHIT) {
        snprintf(part, sizeof(part), "_cdhit%d",
                 (int)(g_cdhit_id * 1000.0 + 0.5));
        add(part);
    }

    if (merger_mode == MERGER_PEAR) {
        snprintf(part, sizeof(part), "_p%d",
                 (int)(g_pear_pvalue * 1000.0 + 0.5));
        add(part);
    }

    if (merger_mode == MERGER_NGMERGE) {
        snprintf(part, sizeof(part), "_p%d",
                 (int)(g_ngmerge_mismatch * 1000.0 + 0.5));
        add(part);
    }

    if (merger_mode == MERGER_FLASH) {
        snprintf(part, sizeof(part), "_M%d", g_flash_max_overlap);
        add(part);
        snprintf(part, sizeof(part), "_x%d",
                 (int)(g_flash_mismatch * 1000.0 + 0.5));
        add(part);
    }

    if (alt_flip_budget() != g_max_alt) {
        snprintf(part, sizeof(part), "_f%d", alt_flip_budget());
        add(part);
    }

    if (g_max_alt != DEFAULT_MAX_ALT) {
        snprintf(part, sizeof(part), "_k%d", g_max_alt);
        add(part);
    }

    if (alt_mode != ALT_OFF && g_max_qsum != DEFAULT_MAX_QSUM) {
        if (g_max_qsum < 0)
            snprintf(part, sizeof(part), "_qsoff");
        else
            snprintf(part, sizeof(part), "_qs%d", g_max_qsum);
        add(part);
    }

    if (!g_i1_mm_order)
        add("_i1mm0");

    if (g_tag && g_tag[0]) {
        add("_");
        add(g_tag);
    }

    const char* out_dir = RESULT_DIR;

    {
        std::error_code ec;
        std::filesystem::create_directories(out_dir, ec);

        if (ec)
            fprintf(stderr, "Warning: Cannot create %s: %s\n",
                    out_dir, ec.message().c_str());
    }

    char out_path[1024];
    char summary_path[1024];

    snprintf(out_path, sizeof(out_path), "%s/result_%d_%d_%s.txt",
             out_dir, start_index, end_index, suffix);

    snprintf(summary_path, sizeof(summary_path), "%s/summary_%d_%d_%s.txt",
             out_dir, start_index, end_index, suffix);

    FILE* fp_out = fopen(out_path, "wt");
    FILE* fp_summary = fopen(summary_path, "wt");

    if (!fp_out)
        fprintf(stderr, "Warning: Cannot open result file: %s\n", out_path);
    else {
        fprintf(fp_out,
                "# idx  result  clusters  trial  LT_real_N  observance  "
                "inferred  base_error  read_limit\n");
        fflush(fp_out);
    }

    if (!fp_summary)
        fprintf(stderr, "Warning: Cannot open summary file: %s\n",
                summary_path);
    else {
        fprintf(fp_summary,
                "# read_limit  success_count  total\n");
        fflush(fp_summary);
    }

    std::vector<bool>
        success_flag(total_frames, false);

    int success_count = 0;

    printf("\nDECODING START =======================\n");
    printf("suffix=_%s, read_limit=%d\n\n", suffix, read_limit);

    for (int recursive = start_index;
         recursive <= end_index;
         recursive++)
    {
        const int col = recursive - start_index;

        const std::vector<std::string>& sampled_headers =
            sampled_by_idx[col];

        if (sampled_headers.empty()) {
            printf("[%d] no header was sampled\n", recursive);

            if (fp_out) {
                fprintf(fp_out,
                        "%d  -1  -1  -1  -1  -1  -1  -1.0  %d\n",
                        recursive, read_limit);
                fflush(fp_out);
            }

            continue;
        }

        if (!g_merge_cache)
            merge_cache.clear();

        DecodeResult dr =
            run_decode(
                paired.r1,
                paired.r2,
                merged_ptr,
                merge_cache,
                pear_merger,
                ng_ptr,
                fl_ptr,
                sampled_headers,
                read_limit,
                alt_mode,
                cluster_mode,
                starcode_enabled,
                multilen,
                ml_weighted,
                star_dist,
                msa_bin,
                gf,
                seed_data,
                seed_map,
                Binary_Data_input,
                LT_parity_flat,
                Inferred_bitstream,
                Inferred_idx,
                node_to_strands,
                strand_nodes,
                Current_bitstream_flat);

        if (fp_out) {
            fprintf(fp_out,
                    "%d  %d  %d  %d  %d  %d  %d  %.6f  %d\n",
                    recursive, dr.result, dr.num_clusters,
                    dr.trial_number, dr.LT_real_N,
                    dr.LT_decoding_observance, dr.inferred_count,
                    dr.base_error, read_limit);
            fflush(fp_out);
        }

        if (dr.result == 1) {
            success_flag[col] = true;
            success_count++;
        }

        printf("[%d] %s (result=%d)\n",
               recursive, dr.result == 1 ? "SUCCESS" : "FAIL", dr.result);

        fflush(stdout);
    }

    if (fp_summary) {
        fprintf(fp_summary, "%d  %d  %d\n",
                read_limit, success_count, total_frames);
        fclose(fp_summary);
    }

    if (fp_out) {
        fclose(fp_out);
        printf("\nResult     : %s\n", out_path);
    }

    if (fp_summary)
        printf("Summary    : %s\n", summary_path);

    printf(
        "Success: %d / %d at read_limit=%d\n",
        success_count,
        total_frames,
        read_limit);

    printf(
        "\nTotal time: %.2f minutes\n",
        (double)(clock() - tStart) /
            CLOCKS_PER_SEC /
            60.0);

    free(Binary_Data_input);
    free(LT_parity_flat);
    free(Inferred_bitstream);
    free(Inferred_idx);
    free(Current_bitstream_flat);

    return 0;
}

static void print_usage(const char* prog) {
    printf(
        "CAPER: DNA Fountain LT decoder with RS-guided read recovery\n\n"
        "Usage: %s <start_index> <end_index> <RSN> [alt_mode]\n"
        "          [options]\n\n",
        prog);

    printf(
        "Input : ../data/%s\n"
        "        ../data/%s\n"
        "Output: ../result/\n\n"
        "Positional arguments:\n"
        "  start_index  first idx (trial). The idx is the RNG seed of its\n"
        "               read draw, so rerunning an idx reproduces it.\n"
        "  end_index    last idx; every idx in [start, end] is run.\n"
        "  RSN          number of reads (read_limit) drawn and decoded per\n"
        "               idx. Each idx is decoded once at this read_limit.\n"
        "  alt_mode     which R1/R2 mismatch positions may be flipped when a\n"
        "               read fails RS (and how the consensus vote is weighted):\n"
        "                 0 = off, merged sequence used as-is  (Baseline)\n"
        "                 1 = EQ, equal-quality ties only      (Proposed EQ)\n"
        "                 2 = DQ, every mismatch               (Proposed DQ)\n"
        "               default 1\n\n",
        R1_FASTQ_NAME, R2_FASTQ_NAME);

    printf(
        "Merging:\n"
        "  --merger NAME    pear (default), ngmerge or flash, run in process.\n"
        "                   Any other name reads {base_dir}/merged_{NAME}.fastq\n"
        "                   (no mismatch list, so nothing can be flipped).\n"
        "  --pear-pvalue F  PEAR -p (default %.2f)\n"
        "  --ngmerge-mismatch F   NGmerge -p (default %.2f)\n"
        "  --ngmerge-profile PATH NGmerge quality profile (merged quality\n"
        "                   string only; default ./%s)\n"
        "  --flash-mismatch F     FLASH -x (default %.2f)\n"
        "  --flash-max-overlap N  FLASH -M (default %d)\n\n",
        g_pear_pvalue, g_ngmerge_mismatch, DEFAULT_NGMERGE_PROFILE,
        g_flash_mismatch, g_flash_max_overlap);

    printf(
        "Read-level search (I2, alt_mode 1 or 2):\n"
        "  --max-alt K      K: candidate positions per read, the K cheapest\n"
        "                   eligible mismatches (default %d, max %d)\n"
        "  --max-flips M    M: positions one candidate may change at once\n"
        "                   (default %d, at most K)\n"
        "  --max-qsum T     Q_max: skip candidates whose flipped positions\n"
        "                   have a q_gap sum above T (default %d, -1 = off)\n\n",
        DEFAULT_MAX_ALT, MAX_ALT_CAP, DEFAULT_MAX_FLIPS, DEFAULT_MAX_QSUM);

    printf(
        "Clustering (I3):\n"
        "  --starcode N     cluster-level rescue of the reads that are left\n"
        "                   after RS / read-level search: 0 off (Naive,\n"
        "                   default), 1 on (NTL)\n"
        "  --cluster N      0 = exact match, RS per read (default)\n"
        "                   1 = cluster every read first, RS on the consensus\n"
        "  --cluster-tool NAME   starcode (default) or cdhit\n"
        "  --star-dist D    Starcode edit distance, 1..8 (default %d)\n"
        "  --cdhit-id F     CD-HIT identity threshold, 0.80..1.0 (default %.2f)\n"
        "  --multilen L     reads of length 152 +- L also enter clustering,\n"
        "                   0..8 (default %d; 0 = 152 only). Ignored without\n"
        "                   a clustering stage.\n"
        "  --ml-consensus NAME   majority (default) or weighted. weighted\n"
        "                   needs alt_mode 1 or 2 and --multilen > 0.\n"
        "  --starcode-bin PATH   Starcode (default %s/starcode)\n"
        "  --cdhit-bin PATH      cd-hit-est (default %s/cd-hit-est)\n"
        "  --msa PATH            MUSCLE v5 (default %s/muscle)\n"
        "                   A tool may also sit in its own subfolder of\n"
        "                   sequence_analysis/ (e.g. starcode/starcode,\n"
        "                   cdhit/cd-hit-est, muscle/muscle-linux-x86.v5.3).\n\n",
        DEFAULT_STARCODE_DISTANCE, g_cdhit_id, DEFAULT_MULTILEN,
        TOOL_DIR, TOOL_DIR, TOOL_DIR);

    printf(
        "Decoder order:\n"
        "  --i1-mm-order N  1 = single-read I1 sequences ordered by fewest\n"
        "                   R1/R2 mismatches (default); 0 = by sequence\n\n");

    printf(
        "Paths (defaults are under ../data/):\n"
        "  --image PATH       original bits    (default %s)\n"
        "  --seed-table PATH  '<seed> <chunk> <chunk> ...', 1-based\n"
        "                     (default %s)\n"
        "  --base-dir PATH    folder of the FASTQs (default %s)\n"
        "  --r1 PATH          R1 FASTQ (default {base_dir}/%s)\n"
        "  --r2 PATH          R2 FASTQ (default {base_dir}/%s)\n"
        "  --tmp-dir PATH     temporary files (default /tmp; /dev/shm keeps\n"
        "                     them in memory)\n"
        "  --tag NAME         append _NAME to the output filenames\n\n",
        DEFAULT_IMAGE_FILE, DEFAULT_SEED_TABLE_FILE, DATA_DIR,
        R1_FASTQ_NAME, R2_FASTQ_NAME);

    printf(
        "Speed (no effect on results):\n"
        "  --merge-jobs N   merge each round's reads on N threads (default 1)\n"
        "  --msa-jobs N     MUSCLE processes in parallel (default 1)\n"
        "  --tool-threads N threads per starcode / cd-hit-est / MUSCLE call\n"
        "                   (default 5)\n"
        "  --merge-cache N  1 = reuse merges across idx (default),\n"
        "                   0 = every decode merges its own reads\n\n");

    printf(
        "Output (in %s/, created if missing):\n"
        "  result_<start>_<end>_<suffix>.txt   one line per idx:\n"
        "      idx result clusters trial LT_real_N observance inferred\n"
        "      base_error read_limit\n"
        "      result: 1 success, 0 inference failure, 2 decoded with bit\n"
        "      errors\n"
        "  summary_<start>_<end>_<suffix>.txt  one line:\n"
        "      read_limit success_count total\n\n",
        RESULT_DIR);

    printf("Examples (idx 1-100, PEAR):\n");
    printf("  # Baseline-Naive\n");
    printf("  %s 1 100 80200 0 --merger pear"
           " --starcode 0\n\n", prog);
    printf("  # Proposed-Naive (EQ)\n");
    printf("  %s 1 100 80200 1 --merger pear"
           " --starcode 0 \\\n"
           "      --max-alt 5 --max-flips 5 --max-qsum 20\n\n", prog);
    printf("  # Proposed-Naive (DQ)\n");
    printf("  %s 1 100 80200 2 --merger pear"
           " --starcode 0 \\\n"
           "      --max-alt 5 --max-flips 5 --max-qsum 20\n\n", prog);
    printf("  # Baseline-NTL\n");
    printf("  %s 1 100 78000 0 --merger pear"
           " --starcode 1 \\\n"
           "      --cluster-tool starcode --star-dist 7 --multilen 5\n\n", prog);
    printf("  # Proposed-NTL (EQ)\n");
    printf("  %s 1 100 78000 1 --merger pear"
           " --starcode 1 \\\n"
           "      --cluster-tool starcode --star-dist 7 --multilen 5"
           " --ml-consensus weighted \\\n"
           "      --max-alt 5 --max-flips 5 --max-qsum 20\n\n", prog);
    printf("  # Proposed-NTL (DQ)\n");
    printf("  %s 1 100 78000 2 --merger pear"
           " --starcode 1 \\\n"
           "      --cluster-tool starcode --star-dist 7 --multilen 5"
           " --ml-consensus weighted \\\n"
           "      --max-alt 5 --max-flips 5 --max-qsum 20\n\n", prog);
    printf(
        "  NTL with CD-HIT instead of Starcode: replace\n"
        "      --cluster-tool starcode --star-dist 7\n"
        "  with\n"
        "      --cluster-tool cdhit --cdhit-id 0.94\n"
        "  With --multilen 0, leave out --ml-consensus.\n"
        "  Faster: add --merge-jobs 5 --msa-jobs 3 --tool-threads 3"
        " --tmp-dir /dev/shm\n");
}

int main(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--help") == 0 ||
            strcmp(argv[i], "-h") == 0)
        {
            print_usage(argv[0]);
            return 0;
        }
    }

    if (argc < 4) {
        print_usage(argv[0]);
        return 1;
    }

    int start_index  = atoi(argv[1]);
    int end_index    = atoi(argv[2]);
    int rsn          = atoi(argv[3]);

    if (end_index < start_index) {
        fprintf(stderr,
                "Error: end_index (%d) < start_index (%d)\n",
                end_index, start_index);
        return 1;
    }

    if (rsn <= 0) {
        fprintf(stderr, "Error: RSN must be positive\n");
        return 1;
    }

    int alt_mode = ALT_EQUAL_Q;

    int argi = 4;

    if (argi < argc &&
        strncmp(argv[argi], "--", 2) != 0)
    {
        alt_mode = atoi(argv[argi]);

        if (alt_mode < ALT_OFF || alt_mode > ALT_ALL) {
            fprintf(stderr,
                    "Error: alt_mode must be 0, 1 or 2 (got \"%s\")\n",
                    argv[argi]);
            return 1;
        }

        argi++;
    }

    int cluster_mode     = CLUSTER_EXACT;
    int starcode_enabled = 0;
    int merger_mode      = MERGER_PEAR;
    const char* merger_name = "pear";
    int multilen    = DEFAULT_MULTILEN;
    int ml_weighted = 0;
    int star_dist   = DEFAULT_STARCODE_DISTANCE;
    const char* msa_bin = nullptr;

    auto opt = [&](const char* name) {
        return strcmp(argv[argi], name) == 0 && argi + 1 < argc;
    };

    auto bad01 = [&](int v, const char* name) {
        if (v != 0 && v != 1) {
            fprintf(stderr, "Error: %s must be 0 or 1\n", name);
            return true;
        }
        return false;
    };

    while (argi < argc) {
        const char* val = (argi + 1 < argc) ? argv[argi + 1] : nullptr;

        if (opt("--max-alt")) {
            g_max_alt = atoi(val);

            if (g_max_alt < 1 || g_max_alt > MAX_ALT_CAP) {
                fprintf(stderr,
                        "Error: --max-alt must be between 1 and %d\n",
                        MAX_ALT_CAP);
                return 1;
            }
        }
        else if (opt("--max-flips")) {
            g_max_flips = atoi(val);
            g_max_flips_set = true;

            if (g_max_flips < 1 || g_max_flips > MAX_ALT_CAP) {
                fprintf(stderr,
                        "Error: --max-flips must be between 1 and %d\n",
                        MAX_ALT_CAP);
                return 1;
            }
        }
        else if (opt("--max-qsum")) {
            g_max_qsum = atoi(val);

            if (g_max_qsum < -1) {
                fprintf(stderr,
                        "Error: --max-qsum must be >= 0 (or -1 = off)\n");
                return 1;
            }
        }
        else if (opt("--cluster")) {
            cluster_mode = atoi(val);

            if (cluster_mode != CLUSTER_EXACT &&
                cluster_mode != CLUSTER_STARCODE)
            {
                fprintf(stderr, "Error: --cluster must be 0 or 1\n");
                return 1;
            }
        }
        else if (opt("--starcode")) {
            starcode_enabled = atoi(val);
            if (bad01(starcode_enabled, "--starcode")) return 1;
        }
        else if (opt("--merger")) {
            merger_name = val;

            if (strcmp(merger_name, "pear") == 0) {
                merger_mode = MERGER_PEAR;
            }
            else if (strcmp(merger_name, "ngmerge") == 0) {
                merger_mode = MERGER_NGMERGE;
            }
            else if (strcmp(merger_name, "flash") == 0) {
                merger_mode = MERGER_FLASH;
            }
            else {
                merger_mode = MERGER_EXTERNAL;

                if (!is_known_merger(merger_name)) {
                    fprintf(stderr,
                            "Note: \"%s\" is not a known merger name; its FASTQ\n"
                            "      must exist at {base_dir}/merged_%s.fastq.\n",
                            merger_name, merger_name);
                }
            }
        }
        else if (opt("--multilen")) {
            multilen = atoi(val);

            if (multilen < 0 || multilen > 8) {
                fprintf(stderr, "Error: --multilen must be between 0 and 8\n");
                return 1;
            }
        }
        else if (opt("--ml-consensus")) {
            if (strcmp(val, "majority") == 0) {
                ml_weighted = 0;
            }
            else if (strcmp(val, "weighted") == 0) {
                ml_weighted = 1;
            }
            else {
                fprintf(stderr,
                        "Error: --ml-consensus must be majority or weighted\n");
                return 1;
            }
        }
        else if (opt("--cluster-tool")) {
            if (strcmp(val, "starcode") == 0) {
                g_cluster_tool = CLUSTER_TOOL_STARCODE;
            }
            else if (strcmp(val, "cdhit") == 0 ||
                     strcmp(val, "cd-hit") == 0)
            {
                g_cluster_tool = CLUSTER_TOOL_CDHIT;
            }
            else {
                fprintf(stderr,
                        "Error: --cluster-tool must be starcode or cdhit\n");
                return 1;
            }
        }
        else if (opt("--star-dist")) {
            star_dist = atoi(val);

            if (star_dist < 1 || star_dist > 8) {
                fprintf(stderr, "Error: --star-dist must be between 1 and 8\n");
                return 1;
            }
        }
        else if (opt("--cdhit-id")) {
            g_cdhit_id = atof(val);

            if (g_cdhit_id < 0.80 || g_cdhit_id > 1.0) {
                fprintf(stderr,
                        "Error: --cdhit-id must be between 0.80 and 1.0\n");
                return 1;
            }
        }
        else if (opt("--starcode-bin")) {
            g_starcode_bin = val;
        }
        else if (opt("--cdhit-bin")) {
            g_cdhit_bin = val;
        }
        else if (opt("--msa")) {
            msa_bin = val;
        }
        else if (opt("--pear-pvalue")) {
            g_pear_pvalue = atof(val);

            if (g_pear_pvalue <= 0.0 || g_pear_pvalue > 1.0) {
                fprintf(stderr, "Error: --pear-pvalue must be in (0, 1]\n");
                return 1;
            }
        }
        else if (opt("--ngmerge-mismatch")) {
            g_ngmerge_mismatch = atof(val);

            if (g_ngmerge_mismatch < 0.0 || g_ngmerge_mismatch >= 1.0) {
                fprintf(stderr,
                        "Error: --ngmerge-mismatch must be in [0, 1)\n");
                return 1;
            }
        }
        else if (opt("--ngmerge-profile")) {
            g_ngmerge_profile = val;
            g_ngmerge_profile_set = true;
        }
        else if (opt("--flash-mismatch")) {
            g_flash_mismatch = atof(val);

            if (g_flash_mismatch < 0.0 || g_flash_mismatch >= 1.0) {
                fprintf(stderr, "Error: --flash-mismatch must be in [0, 1)\n");
                return 1;
            }
        }
        else if (opt("--flash-max-overlap")) {
            g_flash_max_overlap = atoi(val);

            if (g_flash_max_overlap < 1) {
                fprintf(stderr, "Error: --flash-max-overlap must be positive\n");
                return 1;
            }
        }
        else if (opt("--i1-mm-order")) {
            g_i1_mm_order = atoi(val);
            if (bad01(g_i1_mm_order, "--i1-mm-order")) return 1;
        }
        else if (opt("--tag")) {
            g_tag = val;

            for (const char* c = g_tag; *c; ++c) {
                if (!(isalnum((unsigned char)*c) ||
                      *c == '-' || *c == '_' || *c == '.'))
                {
                    fprintf(stderr,
                            "Error: --tag may contain only letters, digits,\n"
                            "       '-', '_' and '.' (got \"%s\")\n", g_tag);
                    return 1;
                }
            }
        }
        else if (opt("--r1")) {
            g_r1_path_opt = val;
        }
        else if (opt("--r2")) {
            g_r2_path_opt = val;
        }
        else if (opt("--seed-table")) {
            g_seed_table_opt = val;
        }
        else if (opt("--image")) {
            g_image_opt = val;
        }
        else if (opt("--base-dir")) {
            g_base_dir_opt = val;
        }
        else if (opt("--tmp-dir")) {
            g_tmp_dir = val;
        }
        else if (opt("--tool-threads")) {
            g_tool_threads = std::max(1, atoi(val));
        }
        else if (opt("--msa-jobs")) {
            g_msa_jobs = std::max(1, atoi(val));
        }
        else if (opt("--merge-jobs")) {
            g_merge_jobs = std::max(1, atoi(val));
        }
        else if (opt("--merge-cache")) {
            g_merge_cache = atoi(val);
            if (bad01(g_merge_cache, "--merge-cache")) return 1;
        }
        else {
            fprintf(stderr,
                    "Error: unknown or incomplete option \"%s\"\n"
                    "       (see %s --help)\n",
                    argv[argi], argv[0]);
            return 1;
        }

        argi += 2;
    }

    if (!g_starcode_bin) g_starcode_bin = find_tool(STARCODE_CANDIDATES);
    if (!g_cdhit_bin)    g_cdhit_bin    = find_tool(CDHIT_CANDIDATES);
    if (!msa_bin)        msa_bin        = find_tool(MSA_CANDIDATES);

    if (g_ngmerge_profile_set && merger_mode != MERGER_NGMERGE) {
        fprintf(stderr,
                "Error: --ngmerge-profile applies to --merger ngmerge only.\n");
        return 1;
    }

    if (cluster_mode == CLUSTER_EXACT && !starcode_enabled)
        multilen = 0;

    if (ml_weighted && multilen <= 0) {
        fprintf(stderr,
                "Error: --ml-consensus applies to --multilen > 0 only.\n");
        return 1;
    }

    if (g_max_flips_set && g_max_flips > g_max_alt) {
        fprintf(stderr,
                "Error: --max-flips (%d) cannot exceed --max-alt (%d).\n",
                g_max_flips, g_max_alt);
        return 1;
    }

    if (ml_weighted && alt_mode == ALT_OFF) {
        fprintf(stderr,
                "Error: --ml-consensus weighted needs alt_mode 1 or 2.\n");
        return 1;
    }

    if (cluster_mode == CLUSTER_STARCODE && starcode_enabled) {
        fprintf(stderr,
                "Error: --starcode 1 applies to --cluster 0 only.\n");
        return 1;
    }

    const char* mode_name =
        (alt_mode == ALT_OFF)     ? "OFF" :
        (alt_mode == ALT_EQUAL_Q) ? "EQ" : "DQ";

    printf("start=%d, end=%d, frames=%d, RSN=%d, alt_mode=%s",
           start_index, end_index, end_index - start_index + 1,
           rsn, mode_name);

    printf(", merger=%s", merger_name);

    printf(", cluster=%s",
           (cluster_mode == CLUSTER_STARCODE) ? "starcode" : "exact");

    printf(", rescue=%s", starcode_enabled ? "on" : "off");

    if (starcode_enabled || cluster_mode == CLUSTER_STARCODE) {
        if (g_cluster_tool == CLUSTER_TOOL_CDHIT)
            printf(", cluster_tool=cdhit(%.4f)", g_cdhit_id);
        else
            printf(", cluster_tool=starcode(d=%d)", star_dist);

        printf(", multilen=%d", multilen);

        if (multilen > 0)
            printf(", ml_consensus=%s", ml_weighted ? "weighted" : "majority");
    }

    if (alt_mode != ALT_OFF)
        printf(", max_alt=%d, max_flips=%d, max_qsum=%d",
               g_max_alt, alt_flip_budget(), g_max_qsum);

    printf("\n");

    return decode_dna(
        start_index,
        end_index,
        rsn,
        alt_mode,
        cluster_mode,
        starcode_enabled,
        merger_mode,
        merger_name,
        multilen,
        ml_weighted,
        star_dist,
        msa_bin);
}