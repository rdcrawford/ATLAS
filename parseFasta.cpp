#include <fstream>
#include <iostream>
#include <cstdlib>
#include "BarcodeClass.h"

// -------------------------------------------------------------------------- //
// Author  Ryan D. Crawford
// Date    03/20/2018
// Purpose Parses a fasta file and dynamically allocates a vector of barcode
//         class objects each containing a unique sequence contained in a line
//         in the fasta file.
// -------------------------------------------------------------------------- //

// Input: test file in fasta fromat
// Return: vector containing the sequence in the fasta file
vector <Barcode> *parseFasta(ifstream &fastaFile)
{
  string           line;      // Store the sequence of the even lines
  int              count = 0; // Counter for even vs. odd lines
  vector <Barcode> *barcodes; // Vector to store all of the seqs

  // dynamically allocate memeory for the barcode objects
  barcodes = new std::vector<Barcode>;

  // Iterate through the lines of the file and add each line containing a
  // sequence to the vector
  while(getline(fastaFile, line))
  {
    count ++;
    if (count % 2 == 0)
    {
      barcodes->push_back(Barcode(line));
    }
  }
  return barcodes;
}
