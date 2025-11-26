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
int median_of_three(const vector<Fastq>& arr, int low, int high);
int hoare_partition(vector<Fastq>& arr, int low, int high);
void quicksort(vector<Fastq>& arr, int low, int high);
void write_fastq(const vector<Fastq>& reads, ostream& out);

int main(int argc, char **argv) {
    if (argc < 2) {
        cout << "Use as:  " << argv[0] << " <input.fastq>\n";
        cout << "Example: " << argv[0] << " reads.fastq > sorted_reads.fastq\n";
        return 0;
    }

/*    ifstream testFile(argv[1]);
    if (!testFile.is_open()) {
        cout << "Cannot open file \"" << argv[1] << "\"\n";
        return 1;
    }
    testFile.close();*/
//JM: you are testing it inside read_fastq, so there is no need to do it also here. It requires opening additional stream just for the test.

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


int median_of_three(const vector<Fastq>& arr, int low, int high) {
    int mid = (low + high) / 2;
    string a = arr[low].seq;
    string b = arr[mid].seq;
    string c = arr[high].seq;

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

int hoare_partition(vector<Fastq>& arr, int low, int high) {
    int pivot_index = median_of_three(arr, low, high);
    string pivot = arr[pivot_index].seq;
    int i = low - 1;
    int j = high + 1;

    while (true) {
        do { i++; } while (arr[i].seq < pivot);
        do { j--; } while (arr[j].seq > pivot);
        if (i >= j) return j;
        swap(arr[i], arr[j]);
    }
}


void quicksort(vector<Fastq>& arr, int low, int high) {
    if (low < high) {
        int p = hoare_partition(arr, low, high);
        quicksort(arr, low, p);
        quicksort(arr, p + 1, high);
    }
}


void write_fastq(const vector<Fastq>& reads, ostream& out) {
    for (const auto& read : reads) {
        out << read.header << '\n';
        out << read.seq    << '\n';
        out << read.sign   << '\n';
        out << read.score  << '\n';
    }
}

