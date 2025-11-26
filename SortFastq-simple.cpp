/*

Very straightforward quicksort implementation
no median-of-three

(using code from an unnamed student)

*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

// Define FastqRecord structure
struct FastqRecord {
	string header;
	string sequence;
	string plus;
	string quality;
};

// Lomuto partition
int partition(vector<FastqRecord>& records, int min, int max){

	string pivotSequence = records[max].sequence;

	int i = min;
	for (int j = min; j < max; ++j){
		if (records[j].sequence < pivotSequence){
			FastqRecord temp = records[i];
			records[i] = records[j];
			records[j] = temp;
			i++;
		}
	}
	FastqRecord temp = records[i];
	records[i]=records[max];
	records[max]=temp;
	return (i);

}

// REcursive quicksort
void quickSort (vector<FastqRecord>& records, int min, int max){
	if (min < max){
		int pivot = partition(records, min, max);
		quickSort(records, min, pivot-1);
		quickSort(records, pivot+1, max);
		
	}
}

/*************************************************************/

// Main function
int main (int argc, char * argv[]){
	if (argc != 2) {
		cout << "Usage: " << argv[0] << "<input.fastq>\n";
		return 1;
	}

	ifstream inFile(argv[1]);
	if (!inFile.is_open()){
		cout << "Error: Can not open file \" " << argv[1] << "\"\n";
		return 1;
	}

	// Read input into a vector
	vector<FastqRecord> allRecords;
	FastqRecord tempRecord;

	while (getline(inFile, tempRecord.header)) {
		if (getline(inFile, tempRecord.sequence))
			if (getline(inFile, tempRecord.plus))
				if (getline(inFile, tempRecord.quality))
					allRecords.push_back(tempRecord);
	}

	inFile.close();

	// Sort the records
	quickSort(allRecords, 0, allRecords.size()-1);

	// Write the sorted records to the output
	for (int i=0; i < allRecords.size(); ++i) {
		cout << allRecords[i].header << '\n';
		cout << allRecords[i].sequence << '\n';
		cout << allRecords[i].plus << '\n';
		cout << allRecords[i].quality << '\n';
	}

	return 0;

}

