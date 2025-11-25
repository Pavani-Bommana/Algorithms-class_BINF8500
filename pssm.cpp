/*

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <string>
#include <iomanip>

using namespace std;

const vector<char> BASES = {'A', 'C', 'G', 'T'};
const double DEFAULT_PSEUDOCOUNT = 1.0;

// Reverse complement
string reverse_complement(const string& seq) {
    map<char, char> comp = {{'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}};
    string rev;
    for (int i = seq.size() - 1; i >= 0; --i)
        rev += comp.count(seq[i]) ? comp[seq[i]] : 'N';
    return rev;
}

// Read motif sequences line by line
vector<string> read_motifs(const string& filename) {
    ifstream file(filename);
    vector<string> motifs;
    string line;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;
        motifs.push_back(line);
    }

    if (!motifs.empty()) {
        size_t motif_len = motifs[0].size();
        for (const string& m : motifs) {
            if (m.size() != motif_len) {
                cerr << "Error: Motif lengths are inconsistent.\n";
                exit(1);
            }
        }
    }

    return motifs;
}

// Read DNA sequence from FASTA
string read_dna(const string& filename) {
    ifstream file(filename);
    string line, dna;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;
        dna += line;
    }
    return dna;
}

// Background probabilities from DNA
map<char, double> compute_background(const string& dna) {
    map<char, int> counts;
    int total = 0;
    for (char c : dna) {
        c = toupper(c);
        if (find(BASES.begin(), BASES.end(), c) != BASES.end()) {
            counts[c]++;
            total++;
        }
    }
    map<char, double> bg;
    for (char b : BASES)
        bg[b] = (double)counts[b] / total;
    return bg;
}

// Frequency matrix (returns normalized frequencies and raw counts)
pair<vector<map<char, double>>, vector<map<char, int>>> build_frequency_matrix(const vector<string>& motifs, double pseudocount) {
    int L = motifs[0].size();
    vector<map<char, int>> raw_counts(L);
    vector<map<char, double>> freq(L);

    for (int i = 0; i < L; ++i)
        for (char b : BASES) {
            raw_counts[i][b] = 0;
            freq[i][b] = 0.0;
        }

    for (const string& motif : motifs)
        for (int i = 0; i < L; ++i)
            raw_counts[i][toupper(motif[i])]++;

    for (int i = 0; i < L; ++i) {
        double total = 0;
        for (char b : BASES) {
            freq[i][b] = raw_counts[i][b] + pseudocount;
            total += freq[i][b];
        }
        for (char b : BASES)
            freq[i][b] /= total;
    }

    return {freq, raw_counts};
}

// PSSM: log-ratio scores
vector<map<char, double>> build_pssm(const vector<map<char, double>>& freq, const map<char, double>& bg) {
    int L = freq.size();
    vector<map<char, double>> pssm(L);
    for (int i = 0; i < L; ++i)
        for (char b : BASES)
            pssm[i][b] = log2(freq[i].at(b) / bg.at(b));
    return pssm;
}

// Score a sequence using PSSM
double score_sequence(const string& seq, const vector<map<char, double>>& pssm) {
    double score = 0.0;
    for (int i = 0; i < seq.size(); ++i) {
        char b = toupper(seq[i]);
        if (pssm[i].count(b))
            score += pssm[i].at(b);
        else
            return -INFINITY;
    }
    return score;
}

// Percentile-based cutoff
double percentile_cutoff(const vector<double>& scores, double percentile) {
    vector<double> sorted = scores;
    sort(sorted.begin(), sorted.end());
    int index = (int)(percentile * sorted.size());
    return sorted[index];
}

// Scan strand for motif hits
void scan_strand(const string& dna, const vector<map<char, double>>& pssm, double cutoff, bool reverse) {
    int L = pssm.size();
    string seq = reverse ? reverse_complement(dna) : dna;

    for (int i = 0; i <= seq.size() - L; ++i) {
        string window = seq.substr(i, L);
        double score = score_sequence(window, pssm);

        if (score >= cutoff) {
            if (reverse) {
                int rev_start = dna.size() - i - L;
                int rev_end = dna.size() - i - 1;
                cout << rev_start + 1 << "\t" << rev_end + 1 << "\t-\t"
                     << window << "\t\t" << fixed << setprecision(3) << score << endl;
            } else {
                int start = i;
                int end = start + L - 1;
                cout << start + 1 << "\t" << end + 1 << "\t+\t"
                     << window << "\t\t" << fixed << setprecision(3) << score << endl;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: ./pssm <motif_file> <dna_file> [score_cutoff]" << endl;
        return 1;
    }

    string motif_file = argv[1];
    string dna_file = argv[2];
    double cutoff = -INFINITY;

    vector<string> motifs = read_motifs(motif_file);
    string dna = read_dna(dna_file);

    map<char, double> bg = compute_background(dna);
    auto [freq, raw_counts] = build_frequency_matrix(motifs, DEFAULT_PSEUDOCOUNT);
    auto pssm = build_pssm(freq, bg);

    // Frequency matrix (raw counts)
    cout << "Frequency matrix:\n";
    cout << "pos.\tA\tC\tG\tT\n";
    for (int i = 0; i < raw_counts.size(); ++i) {
        cout << i + 1;
        for (char b : BASES)
            cout << "\t" << raw_counts[i][b];
        cout << endl;
    }

    // Probability matrix
    cout << "\nProbability matrix:\n";
    cout << "pos.\tA\tC\tG\tT\n";
    for (int i = 0; i < freq.size(); ++i) {
        cout << i + 1;
        for (char b : BASES)
            cout << "\t" << fixed << setprecision(3) << freq[i][b];
        cout << endl;
    }

    // PSSM matrix
    cout << "\nPSSM:\n";
    cout << "pos.\tA\tC\tG\tT\n";
    for (int i = 0; i < pssm.size(); ++i) {
        cout << i + 1;
        for (char b : BASES)
            cout << "\t" << fixed << setprecision(3) << pssm[i][b];
        cout << endl;
    }

    // Training motif scores
    vector<double> motif_scores;
    cout << "\nTraining set motif scores:\n";
    for (const string& m : motifs) {
        double s = score_sequence(m, pssm);
        motif_scores.push_back(s);
        cout << m << "\t\t" << fixed << setprecision(3) << s << endl;
    }

    // Score cutoff
    if (argc >= 4)
        cutoff = stod(argv[3]);
    else
        cutoff = percentile_cutoff(motif_scores, 0.1); // 10th percentile

    cout << "\nMatches with score " << fixed << setprecision(3) << cutoff
         << " or higher found in " << dna_file << " (length " << dna.size() << " bp):\n\n";
    cout << "Start\tEnd\tStrand\tSequence\tScore\n";

    scan_strand(dna, pssm, cutoff, false);
    scan_strand(dna, pssm, cutoff, true);

    return 0;
}
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <string>
#include <iomanip>

using namespace std;

const vector<char> BASES = {'A', 'C', 'G', 'T'};
const double PSEUDOCOUNT = 0.25;

// Reverse complement
string reverse_complement(const string& seq) {
    map<char, char> comp = {{'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}};
    string rev;
    for (int i = seq.size() - 1; i >= 0; --i)
        rev += comp.count(seq[i]) ? comp[seq[i]] : 'N';
    return rev;
}

// Read motif sequences line by line
vector<string> read_motifs(const string& filename) {
    ifstream file(filename);
    vector<string> motifs;
    string line;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;
        motifs.push_back(line);
    }

    if (!motifs.empty()) {
        size_t motif_len = motifs[0].size();
        for (const string& m : motifs) {
            if (m.size() != motif_len) {
                cerr << "Error: Motif lengths are inconsistent.\n";
                exit(1);
            }
        }
    }

    return motifs;
}

// Read DNA sequence from FASTA
string read_dna(const string& filename) {
    ifstream file(filename);
    string line, dna;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;
        dna += line;
    }
    return dna;
}

// Background probabilities from DNA
map<char, double> compute_background(const string& dna) {
    map<char, int> counts;
    int total = 0;
    for (char c : dna) {
        c = toupper(c);
        if (find(BASES.begin(), BASES.end(), c) != BASES.end()) {
            counts[c]++;
            total++;
        }
    }
    map<char, double> bg;
    for (char b : BASES)
        bg[b] = (double)counts[b] / total;
    return bg;
}

// Frequency matrix (raw counts and normalized probabilities with pseudocounts)
pair<vector<map<char, double>>, vector<map<char, int>>> build_frequency_matrix(const vector<string>& motifs) {
    int L = motifs[0].size();
    vector<map<char, int>> raw_counts(L);
    vector<map<char, double>> freq(L);

    for (int i = 0; i < L; ++i)
        for (char b : BASES) {
            raw_counts[i][b] = 0;
            freq[i][b] = 0.0;
        }

    for (const string& motif : motifs)
        for (int i = 0; i < L; ++i)
            raw_counts[i][toupper(motif[i])]++;

    for (int i = 0; i < L; ++i) {
        double total = 0;
        for (char b : BASES) {
            freq[i][b] = raw_counts[i][b] + PSEUDOCOUNT;
            total += freq[i][b];
        }
        for (char b : BASES)
            freq[i][b] /= total;
    }

    return {freq, raw_counts};
}

// PSSM: log-ratio scores
vector<map<char, double>> build_pssm(const vector<map<char, double>>& freq, const map<char, double>& bg) {
    int L = freq.size();
    vector<map<char, double>> pssm(L);
    for (int i = 0; i < L; ++i)
        for (char b : BASES)
            pssm[i][b] = freq[i].at(b) > 0 ? log2(freq[i].at(b) / bg.at(b)) : -INFINITY;
    return pssm;
}

// Score a sequence using PSSM
double score_sequence(const string& seq, const vector<map<char, double>>& pssm) {
    double score = 0.0;
    for (int i = 0; i < seq.size(); ++i) {
        char b = toupper(seq[i]);
        if (pssm[i].count(b))
            score += pssm[i].at(b);
        else
            return -INFINITY;
    }
    return score;
}

// Percentile-based cutoff
double percentile_cutoff(const vector<double>& scores, double percentile) {
    vector<double> sorted = scores;
    sort(sorted.begin(), sorted.end());
    int index = (int)(percentile * sorted.size());
    return sorted[index];
}

// Scan strand for motif hits
void scan_strand(const string& dna, const vector<map<char, double>>& pssm, double cutoff, bool reverse) {
    int L = pssm.size();
    string seq = reverse ? reverse_complement(dna) : dna;

    for (int i = 0; i <= seq.size() - L; ++i) {
        string window = seq.substr(i, L);
        double score = score_sequence(window, pssm);

        if (score >= cutoff) {
            if (reverse) {
                int rev_start = dna.size() - i - L;
                int rev_end = dna.size() - i - 1;
                cout << rev_start + 1 << "\t" << rev_end + 1 << "\t-\t"
                     << window << "\t\t" << fixed << setprecision(3) << score << endl;
            } else {
                int start = i;
                int end = start + L - 1;
                cout << start + 1 << "\t" << end + 1 << "\t+\t"
                     << window << "\t\t" << fixed << setprecision(3) << score << endl;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: ./pssm <motif_file> <dna_file> [score_cutoff]" << endl;
        return 1;
    }

    string motif_file = argv[1];
    string dna_file = argv[2];
    double cutoff = -INFINITY;

    vector<string> motifs = read_motifs(motif_file);
    string dna = read_dna(dna_file);

    map<char, double> bg = compute_background(dna);
    auto [freq, raw_counts] = build_frequency_matrix(motifs);
    auto pssm = build_pssm(freq, bg);

    // Frequency matrix (raw counts)
    cout << "Frequency matrix:\n";
    cout << "pos.\tA\tC\tG\tT\n";
    for (int i = 0; i < raw_counts.size(); ++i) {
        cout << i + 1;
        for (char b : BASES)
            cout << "\t" << raw_counts[i][b];
        cout << endl;
    }

    // Probability matrix (with pseudocounts)
    cout << "\nProbability matrix (with pseudocount = " << PSEUDOCOUNT << "):\n";
    cout << "pos.\tA\tC\tG\tT\n";
    for (int i = 0; i < freq.size(); ++i) {
        cout << i + 1;
        for (char b : BASES)
            cout << "\t" << fixed << setprecision(3) << freq[i][b];
        cout << endl;
    }

    // PSSM matrix
    cout << "\nPSSM:\n";
    cout << "pos.\tA\tC\tG\tT\n";
    for (int i = 0; i < pssm.size(); ++i) {
        cout << i + 1;
        for (char b : BASES)
            cout << "\t" << fixed << setprecision(3) << pssm[i][b];
        cout << endl;
    }

    // Training motif scores
    vector<double> motif_scores;
    cout << "\nTraining set motif scores:\n";
    for (const string& m : motifs) {
        double s = score_sequence(m, pssm);
        motif_scores.push_back(s);
        cout << m << "\t\t" << fixed << setprecision(3) << s << endl;
    }

    // Score cutoff
    if (argc >= 4)
        cutoff = stod(argv[3]);
    else
        cutoff = percentile_cutoff(motif_scores, 0.5); // 10th percentile

    cout << "\nMatches with score " << fixed << setprecision(3) << cutoff
         << " or higher found in " << dna_file << " (length " << dna.size() << " bp):\n\n";
    cout << "Start\tEnd\tStrand\tSequence\tScore\n";

scan_strand(dna, pssm, cutoff, false);
scan_strand(dna, pssm, cutoff, true);

return 0;
}