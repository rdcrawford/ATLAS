#include <string>
#include <vector>
using namespace std;

// -------------------------------------------------------------------------- //
// Author  Ryan D. Crawford
// Date    02/17/2018
// Purpose This file defines a Barcode class functions used to asses and define
//         attributes: sequence, indicies of reads containing that read in a
//         respective fastq file, and the number of reads aligning to that
//         Barcode.
// -------------------------------------------------------------------------- //

#ifndef _BARCODE_H_
#define _BARCODE_H_
class Barcode
{
  public:
    // Default Ctor
    Barcode():numMappedReads(0)
    { ; }

    // Input: barcode Sequence
    // Description: Value Ctor to create a barcode class object from the input
    //              input parameters.
    Barcode(string inBarcodeSeq);

    // Return Number of reads aligning to that barcode.
    int getReadCount();

    // Desctipion Getter function for the sequence of a specific barcode
    //            specific barcode.
    string getBarcodeSeq();

    //Input  barcode to match
    //Return true if querry is a match for this barcode
    bool matchBarcodes(string querryBarcode, bool const &IsReverseStrand);

    // Adds one to the count of reads mapping to a barcode
    void incrementReadCount();

  private:
    string barcodeSeq;     // DNA sequence of the barcode
    int    numMappedReads; // Count of the reads mapping to the barcode

    // Function to calculate the reverse compliment of querry sequence
    string getRevComp();
};
#endif
