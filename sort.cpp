
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>

using namespace std;

vector<vector<string>> read_fastq(const string& filename);
int median_of_three(const vector<vector<string>>& arr, int low, int high);
int hoare_partition(vector<vector<string>>& arr, int low, int high);
void quicksort(vector<vector<string>>& arr, int low, int high);
void write_fastq(const vector<vector<string>>& reads, ostream& out);

int main(int argc, char **argv) {
    if (argc < 2) {
        cout << "Use as:  " << argv[0] << " <input.fastq>\n";
        cout << "Example: " << argv[0] << " reads.fastq > sorted_reads.fastq\n";
        return 0;
    }

    ifstream testFile(argv[1]);
    if (!testFile.is_open()) {
        cout << "Cannot open file \"" << argv[1] << "\"\n";
        return 1;
    }
    testFile.close();

    string input_file = argv[1];

    
    vector<vector<string>> reads = read_fastq(input_file);

    quicksort(reads, 0, reads.size() - 1);


    write_fastq(reads, cout);

    return 0;
}


vector<vector<string>> read_fastq(const string& filename) {
    vector<vector<string>> reads;
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error: Cannot open input file " << filename << endl;
        exit(1);
    }

    string line;
    vector<string> block;
    while (getline(infile, line)) {
        block.push_back(line);
        if (block.size() == 4) {
            reads.push_back(block);
            block.clear();
        }
    }
    infile.close();
    return reads;
}

int median_of_three(const vector<vector<string>>& arr, int low, int high) {
    int mid = (low + high) / 2;
    string a = arr[low][1];
    string b = arr[mid][1];
    string c = arr[high][1];

    if (a < b) {
        if (b < c) return mid;
        else if (a < c) return high;
        else return low;
    } else {
        if (a < c) return low;
        else if (b < c) return high;
        else return mid;
    }
}

int hoare_partition(vector<vector<string>>& arr, int low, int high) {
    int pivot_index = median_of_three(arr, low, high);
    string pivot = arr[pivot_index][1];
    int i = low - 1;
    int j = high + 1;

    while (true) {
        do { i++; } while (arr[i][1] < pivot);
        do { j--; } while (arr[j][1] > pivot);
        if (i >= j) return j;
        swap(arr[i], arr[j]);
    }
}


void quicksort(vector<vector<string>>& arr, int low, int high) {
    if (low < high) {
        int p = hoare_partition(arr, low, high);
        quicksort(arr, low, p);
        quicksort(arr, p + 1, high);
    }
}

void write_fastq(const vector<vector<string>>& reads, ostream& out) {
    for (const auto& read : reads) {
        for (const auto& line : read) {
            out << line << "\n";
        }
    }
}

/*

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

struct Fastq {
    string header;
    string seq;
    string sign;
    string score;
};

vector<Fastq> read_fastq(const string& filename);
void Sort3(Fastq* X, Fastq* Y, Fastq* Z);
int hoare_partition(vector<Fastq>& arr, int low, int high);
void quicksort(vector<Fastq>& arr, int low, int high);
void write_fastq(const vector<Fastq>& reads, ostream& out);

int main(int argc, char **argv) {
    if (argc < 2) {
        cout << "Use as:  " << argv[0] << " <input.fastq>\n";
        cout << "Example: " << argv[0] << " reads.fastq > sorted_reads.fastq\n";
        return 0;
    }

    ifstream testFile(argv[1]);
    if (!testFile.is_open()) {
        cout << "Cannot open file \"" << argv[1] << "\"\n";
        return 1;
    }
    testFile.close();

    string input_file = argv[1];

    vector<Fastq> reads = read_fastq(input_file);

    quicksort(reads, 0, reads.size() - 1);

    write_fastq(reads, cout);

    return 0;
}


vector<Fastq> read_fastq(const string& filename) {
    vector<Fastq> reads;
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error: Cannot open input file " << filename << endl;
        exit(1);
    }

    string h, s, si, sc;
    while (getline(infile, h) && getline(infile, s) &&
           getline(infile, si) && getline(infile, sc)) {
        reads.push_back({h, s, si, sc});
    }

    infile.close();
    return reads;
}

void Sort3(Fastq* X, Fastq* Y, Fastq* Z) {
    if (Y->seq > Z->seq) swap(*Y, *Z);
    if (X->seq > Y->seq) {
        swap(*X, *Y);
        if (Y->seq > Z->seq) swap(*Y, *Z);
    }
}


int hoare_partition(vector<Fastq>& arr, int low, int high) {
    int mid = (low + high) / 2;

    Sort3(&arr[low], &arr[mid], &arr[high]);
    swap(arr[low], arr[mid]); // move median to low

    string pivot = arr[low].seq;
    int i = low - 1;
    int j = high + 1;

    while (true) {
        do { i++; } while (arr[i].seq < pivot);
        do { j--; } while (arr[j].seq > pivot);
        if (i >= j) return j;
        swap(arr[i], arr[j]);
    }
}

// --- QuickSort function ---
void quicksort(vector<Fastq>& arr, int low, int high) {
    if (low < high) {
        int p = hoare_partition(arr, low, high);
        quicksort(arr, low, p);
        quicksort(arr, p + 1, high);
    }
}

// --- Write sorted FASTQ to output ---
void write_fastq(const vector<Fastq>& reads, ostream& out) {
    for (const auto& read : reads) {
        out << read.header << '\n';
        out << read.seq    << '\n';
        out << read.sign   << '\n';
        out << read.score  << '\n';
    }
}
*/