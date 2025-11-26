/*

Very straightforward quicksort implementation (no median-of-three)

Updated to sort an arrray of pointers (which are small) instead of the vector (or array) of actual records
As a result, the records (which can be large and swapping them may require moving large amount of data in memory)
never move. Instead, the sorting takes place in an array of pointers, which take only 8 bytes each and swapping 
them can be mmuch faster.

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


/*************************************************************/

// Lomuto partition
//int partition(vector<FastqRecord>& records, int min, int max){
//NEW: same change as in quickSort below:
int partition(FastqRecord **records, int min, int max){

//	string pivotSequence = records[max].sequence;
//NEW: records[max] is now a pointer to the structure, not the structure itself, and when the structure is
//     itself referenced by a pointer, one has to use operator -> instead of . to refer to its members:
	string pivotSequence = records[max]->sequence;

	int i = min;
	for (int j = min; j < max; ++j){
//		if (records[j].sequence < pivotSequence){
//NEW: change . to ->
		if (records[j]->sequence < pivotSequence){
//			FastqRecord temp = records[i];
//NEW: since records[i] is a pointer to FastqRecord, temp has to be declared as a pointer, too
			FastqRecord *temp = records[i];
			records[i] = records[j];
			records[j] = temp;
			i++;
		}
	}
//	FastqRecord temp = records[i];
//NEW: same here:
	FastqRecord *temp = records[i];
	records[i]=records[max];
	records[max]=temp;
	return (i);

}

/*************************************************************/

// Recursive quicksort
//void quickSort (vector<FastqRecord>& records, int min, int max){
//NEW: the first argument is now array of pointers to FastqRecord, passed by a pointer.
//     So, records is a pointer to a pointer to FastqRecord
void quickSort (FastqRecord **records, int min, int max){
//NEW: records has a different meaning but otherwise the function remains the same
	if (min < max){
		int pivot = partition(records, min, max);
		quickSort(records, min, pivot-1);
		quickSort(records, pivot+1, max);
	}
}

/*************************************************************/
/*************************************************************/
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
	
//NEW: set up an array of pointers to each record:
		FastqRecord **recptr=new FastqRecord*[allRecords.size()];  //array recptr with elements of type "pointer to FastqRecord"
		for (int i=0;i<allRecords.size();++i)  recptr[i]=&allRecords[i];
// recptr[i] now points to where the record allRecords[i] is stored.

	// Sort the records
//	quickSort(allRecords, 0, allRecords.size()-1);
//NEW: instead of the vector of whole records, we'll be sorting the arrray of pointers:
	quickSort(recptr, 0, allRecords.size()-1);

	// Write the sorted records to the output
	for (int i=0; i < allRecords.size(); ++i) {
//NEW: print the output using the pointers, so that it comes out in teh proper sorted order:
//		cout << allRecords[i].header << '\n';
//		cout << allRecords[i].sequence << '\n';
//		cout << allRecords[i].plus << '\n';
//		cout << allRecords[i].quality << '\n';
		cout << recptr[i]->header << '\n';
		cout << recptr[i]->sequence << '\n';
		cout << recptr[i]->plus << '\n';
		cout << recptr[i]->quality << '\n';
// when referencing members of a structure that is itself referenced by a pointer, replace operator . with ->
	}

	return 0;

}

