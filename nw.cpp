/*#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

// Function to read sequence from a FASTA file
string read_fasta(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }

    string line, sequence;
    while (getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') continue; // Skip header
        sequence += line;
    }
    return sequence;
}

// Function to perform Needleman-Wunsch alignment
void needleman_wunsch(const string& seq1, const string& seq2, int match, int mismatch, int gap) {
    size_t m = seq1.length();
    size_t n = seq2.length();

    // Create score and direction matrices
    vector<vector<int>> score(m + 1, vector<int>(n + 1));
    vector<vector<char>> traceback(m + 1, vector<char>(n + 1));

    // Initialize matrices
    for (size_t i = 0; i <= m; ++i) {
        score[i][0] = i * gap;
        traceback[i][0] = 'U'; // Up
    }
    for (size_t j = 0; j <= n; ++j) {
        score[0][j] = j * gap;
        traceback[0][j] = 'L'; // Left
    }
    traceback[0][0] = '0';

    // Fill matrices
    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            int score_diag = score[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch);
            int score_up = score[i - 1][j] + gap;
            int score_left = score[i][j - 1] + gap;

            score[i][j] = max({score_diag, score_up, score_left});

            if (score[i][j] == score_diag)
                traceback[i][j] = 'D';
            else if (score[i][j] == score_up)
                traceback[i][j] = 'U';
            else
                traceback[i][j] = 'L';
        }
    }

    // Traceback
    string aligned_seq1, aligned_seq2;
    size_t i = m, j = n;

    while (i > 0 || j > 0) {
        if (traceback[i][j] == 'D') {
            aligned_seq1 += seq1[i - 1];
            aligned_seq2 += seq2[j - 1];
            --i; --j;
        } else if (traceback[i][j] == 'U') {
            aligned_seq1 += seq1[i - 1];
            aligned_seq2 += '-';
            --i;
        } else { // 'L'
            aligned_seq1 += '-';
            aligned_seq2 += seq2[j - 1];
            --j;
        }
    }

    reverse(aligned_seq1.begin(), aligned_seq1.end());
    reverse(aligned_seq2.begin(), aligned_seq2.end());

    // Output
    cout << "\nAlignment Result:\n";
    cout << "Sequence 1: " << aligned_seq1 << endl;
    cout << "Sequence 2: " << aligned_seq2 << endl;
    cout << "Alignment Score: " << score[m][n] << endl;
}

// Main function
int main(int argc, char* argv[]) {
    if (argc != 6) {
        cerr << "Usage: " << argv[0] << " <file1.fasta> <file2.fasta> <match> <mismatch> <gap>\n";
        return 1;
    }

    string file1 = argv[1];
    string file2 = argv[2];
    int match = stoi(argv[3]);
    int mismatch = stoi(argv[4]);
    int gap = stoi(argv[5]);

    string seq1 = read_fasta(file1);
    string seq2 = read_fasta(file2);

    needleman_wunsch(seq1, seq2, match, mismatch, gap);

    return 0;
}*/




#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

// Function to read sequence from a FASTA file
string read_fasta(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }

    string line, sequence;
    while (getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') continue; // Skip header
        sequence += line;
    }
    return sequence;
}

// Custom function to get max of three integers
int max3(int a, int b, int c) {
    if (a >= b && a >= c) return a;
    else if (b >= a && b >= c) return b;
    else return c;
}
void reverse_string(string& s) {
    size_t i = 0, j = s.length() - 1;
    while (i < j) {
        swap(s[i], s[j]);
        ++i;
        --j;
    }
}

// Function to perform Needleman-Wunsch alignment
void needleman_wunsch(const string& seq1, const string& seq2, int match, int mismatch, int gap) {
    size_t m = seq1.length();
    size_t n = seq2.length();

    // Create score and direction matrices
    vector<vector<int>> score(m + 1, vector<int>(n + 1));
    vector<vector<char>> traceback(m + 1, vector<char>(n + 1));

    // Initialize matrices
    for (size_t i = 0; i <= m; ++i) {
        score[i][0] = i * gap;
        traceback[i][0] = 'U'; // Up
    }
    for (size_t j = 0; j <= n; ++j) {
        score[0][j] = j * gap;
        traceback[0][j] = 'L'; // Left
    }
    traceback[0][0] = '0';

    // Fill matrices
    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            int score_diag = score[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch);
            int score_up = score[i - 1][j] + gap;
            int score_left = score[i][j - 1] + gap;

            score[i][j] = max3(score_diag, score_up, score_left);

            if (score[i][j] == score_diag)
                traceback[i][j] = 'D';
            else if (score[i][j] == score_up)
                traceback[i][j] = 'U';
            else
                traceback[i][j] = 'L';
        }
    }

    // Traceback
    string aligned_seq1, aligned_seq2;
    size_t i = m, j = n;

    while (i > 0 || j > 0) {
        if (traceback[i][j] == 'D') {
            aligned_seq1 += seq1[i - 1];
            aligned_seq2 += seq2[j - 1];
            --i; --j;
        } else if (traceback[i][j] == 'U') {
            aligned_seq1 += seq1[i - 1];
            aligned_seq2 += '-';
            --i;
        } else { // 'L'
            aligned_seq1 += '-';
            aligned_seq2 += seq2[j - 1];
            --j;
        }
    }

	reverse_string(aligned_seq1);
	reverse_string(aligned_seq2); 
    // Output
    cout << "\nAlignment Result:\n";
    cout << "Sequence 1: " << aligned_seq1 << endl;
    cout << "Sequence 2: " << aligned_seq2 << endl;
    cout << "Alignment Score: " << score[m][n] << endl;
}

// Main function
int main(int argc, char* argv[]) {
    if (argc != 6) {
        cerr << "Usage: " << argv[0] << " <file1.fasta> <file2.fasta> <match> <mismatch> <gap>\n";
        return 1;
    }

    string file1 = argv[1];
    string file2 = argv[2];
    int match = stoi(argv[3]);
    int mismatch = stoi(argv[4]);
    int gap = stoi(argv[5]);

    string seq1 = read_fasta(file1);
    string seq2 = read_fasta(file2);

    needleman_wunsch(seq1, seq2, match, mismatch, gap);

    return 0;
}
