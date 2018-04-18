#include "BarcodeClass.h"
#include <string>
#include <string.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdlib>
using namespace std;

// -------------------------------------------------------------------------- //
// Author  Ryan D. Crawford
// Date    02/17/2018
// Purpose This file defines a Barcode class functions used to asses and define
//         attributes: sequence, indicies of reads containing that read in a
//         respective fastq file, and the number of reads aligning to that
//         Barcode.
// -------------------------------------------------------------------------- //

// Input       1) Barcode Sequence
//             2) basePairCount
// Description value Ctor to create a Barcode class object from the input
//             input parameters.
Barcode::Barcode(string inBarcodeSeq):barcodeSeq(inBarcodeSeq),
         numMappedReads(0)
{ ; }

// Desctipion Getter function for the sequence of a specific Barcode
//            specific Barcode.
string Barcode::getBarcodeSeq()
{
  return barcodeSeq;
}

// Input      The number of the mate pair
// Return     Number of reads aligning to that barcode.
// Desctipion Getter function for the number of reads associated with a
//            specific barcode.
int Barcode::getReadCount()
{
  return numMappedReads;
}

// Adds one to the count of reads mapping to a barcode
void Barcode::incrementReadCount()
{
  numMappedReads ++;
}

// Input  barcode to match
// Return true if querry is a match for this barcode
bool Barcode::matchBarcodes(string querry, bool const &IsReverseStrand)
{
  int seqPos;
  //querry.erase(std::remove(querry.begin(), querry.end(), '-'), querry.end());
  if (IsReverseStrand)
  {
    seqPos = querry.find(barcodeSeq);
  }
  else
  {
    string seqRevComp = getRevComp();
    seqPos = querry.find(seqRevComp);
  }

  if (seqPos == -1)
  {
    return false;
  }
  else
  {
    return true;
  }
}

// Function to calculate the reverse compliment of querry sequence
string Barcode::getRevComp()
{
  string revCompSeq = barcodeSeq;
  reverse(revCompSeq.begin(), revCompSeq.end());
  for(int i = 0; i < revCompSeq.length(); i++)
  {
    if(revCompSeq[i] == 'A')
    {
      revCompSeq[i] = 'T';
    }
    else if(revCompSeq[i] == 'T')
    {
      revCompSeq[i] = 'A';
    }
    else if(revCompSeq[i] == 'C')
    {
      revCompSeq[i] = 'G';
    }
    else if(revCompSeq[i] == 'G')
    {
      revCompSeq[i] = 'C';
    }
  }
  return revCompSeq;
}
