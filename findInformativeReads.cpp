#include "api/BamAlignment.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "BarcodeClass.h"
#include <string>
#include "constants.h"
#include "readData.h"
using namespace BamTools;

// -------------------------------------------------------------------------- //
// Author  Ryan D. Crawford
// Date    03/20/2018
// Purpose Step 4 of ATLAS. Creates a dynamically allocated vector with all of
//         the informativeReads.
// -------------------------------------------------------------------------- //

// Input       1) Name of the bamFile
//             2) Name of the bam index file
//             3) Name of the chromosome to be searched
//             4-5) Start and end positions for the barcode within the chr
//             6) Vector containing the barcode sequences
// Return      A vector containing the indicies of the mate pairs for reads
//             containing barcode sequences
// Description Uses the bamtools api to find reads mapping to a specific region
//             specified by the arguments and returns the index of the mate pair
//             corresponding to that read
std::vector <readData> *findInformativeReads(std::string bamFile,
                        std::string baiFile, std::string chr, int startPos,
                        int endPos, vector <Barcode> *barcodes)
{
  // Variable Initializations:
  BamReader    reader; // open the BAM reader
  BamWriter    writer; // Open the BAM writer
  BamAlignment algn; // sore the alignment
  bool         isOpen; // indicates wheter the bam file was opened
  int          refID; // Locus of interest for the alignment
  int          mateID = 0; // ID of the mate pair for the current alignment
  readData     currentRead; // struct to store data on the current read
  int          counter = 0;

  // stores the index of informative mates
  std::vector <readData> *informativeReads;
  informativeReads = new std::vector <readData>;

  // Open the bam file
  isOpen = reader.Open(bamFile);

  // If the BAM file was not opened properly print an error and exit
  if (!isOpen)
  {
    std::cout << "  Unable to open BAM file: " << bamFile << std::endl;
    exit(1);
  }

  // Define the region of interest
  refID = 0;//reader.GetReferenceID(chr);

  // BamTools::BamRegion region(refID, startPos, refID, endPos);
  // reader.SetRegion(region);

  std::string alignedBases; // string to store the sequence of alignment

  // iterate through alignments in this region,
  // ignoring alignments with a Map quality below the minimium below some cutoff
  while (reader.GetNextAlignment(algn))
  {
    // Update to the bases aligned in this read
    alignedBases  = algn.QueryBases;
    bool isBcAlgn = false; // indicates that a read aligns to the barcode locus

    if (algn.RefID == refID)
    {
      if (algn.Position <= startPos && algn.Position + algn.Length >= endPos)
      {
        isBcAlgn = true;
      }
      else if(algn.GetEndPosition() - algn.Length <= startPos &&
              algn.GetEndPosition() >= endPos)
      {
        isBcAlgn = true;
      }
    }

    if (isBcAlgn)
    {
      int  bcNum    = 0;     // Counter for the barcode index
      bool isMatch  = false; // Indicates matching barcode in the read

      // iterate through the barcodes until the barcode is found or all the
      // barcodes have been found
      while(!isMatch && bcNum < barcodes->size())
      {
        // Look for a matching barcode
        isMatch = (*barcodes)[bcNum].matchBarcodes(alignedBases,
                                                   algn.IsReverseStrand());

        if (!isMatch) bcNum ++; // increment the counts
      }

      if (isMatch) // if the read is a match, mark it for further analysis
      {
        (*barcodes)[bcNum].incrementReadCount(); // increment the count
        // Dudate the current read data
        currentRead.bcNum = bcNum;

        // Append read index to the barcode's vector for that mate pair
        if (algn.IsFirstMate())
        {
          currentRead.readIndex = counter + 1; //For first mate: next read
        }
        else
        {
          currentRead.readIndex = counter - 1; //For second mate: previous read
        }
        // Add current read data to the vector
        informativeReads->push_back(currentRead);

        // Get the index of the mate pair
        mateID ++;
      }
    }
    counter ++;
  }

  // close the reader
  reader.Close();
  return informativeReads;
}
