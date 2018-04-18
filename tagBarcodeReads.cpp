#include <fstream>
#include <algorithm>
#include <string>
#include "api/BamAlgorithms.h"
#include "api/BamAlignment.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "UnmappedRead.h"
#include "constants.h"
#include "readData.h"
using namespace BamTools;

// -------------------------------------------------------------------------- //
// Author  Ryan D. Crawford
// Date    04/4/2018
// Purpose Step 7 of ATLAS. This function for ATLAS uses the bamtools C++ API to
//         add flags to the BAM file containing the processed reads to indicate
//         which barcode the restriction fragment is paired to. The BAM file is
//         then sorted by flag.
// -------------------------------------------------------------------------- //

// Input       1) Vector containting the vectro of structs with information on
//                the informative reads
//             2) Bam file containing the alignemts for the informative reads
// Description This function takes the data on the informative reads and adds
//             a tag with the number of the barcode that the read belongs to
void tagBarcodeReads(vector <readData> *informativeReads, string bamFile)
{
  BamReader    reader;       //Reader fo the bam file
  BamWriter    writer;       //Writer for the alignments with tags added
  BamAlignment algn;         //Alignment in bam file
  int          readNum = -1; //Counter for which read is being processed
  int          currentBcNum; //Barcode contained in the opposite mate
  bool         isOpen;       //Indicates a file was opended properly
  string       prevAlgnName; //Name of the previous alignment

  //Attempt to open our BamMultiReader
  isOpen = reader.Open(bamFile);
  if (!isOpen)
  {
    cout << "Could not open input BAM files" << endl;
    exit(1);
  }

  // Retrieve metadata from BAM files, these are required by BamWriter
  const SamHeader HEADER     = reader.GetHeader();
  const RefVector REFERENCES = reader.GetReferenceData();

  //Create a new vector for the dynamically allovated alignments
  std::vector <BamAlignment> *algnVec;
  algnVec = new std::vector <BamAlignment>;

  //Iterate through all alignments, only keeping ones with high map quality
  while (reader.GetNextAlignmentCore(algn))
  {
    //If this is not a restriction fragment from the previously processed read,
    if (prevAlgnName != algn.Name)
    {
      readNum ++; // Increment the alignment number

      //Change the name of the variable to the that of the alignment that was
      //just processed
      prevAlgnName = algn.Name;
    }

    currentBcNum = (*informativeReads)[readNum].bcNum; //Find which barcode

    //Add a tag to the bam file so it can be sorted later
    bool isUpdated = algn.AddTag(TAG, TAG_TYPE, currentBcNum);

    //Only keep alignments with high map quality
    if (algn.MapQuality >= MIN_MAP_QUALITY)
    {
      algnVec->push_back(algn);
    }
  }

  //Sort by the barcode tag defined for the barcode num
  std::sort(algnVec->begin(), algnVec->end(), Algorithms::Sort::ByTag<int>(TAG));

  //Attempt to open our BamWriter
  writer.Open(bamFile, HEADER, REFERENCES);

  //For each alignement in the vector write to the output file
  for (int i = 0; i < algnVec->size(); i ++)
  {
    writer.SaveAlignment((*algnVec)[i]); //Save the current alignment
  }

  delete algnVec; //Deallocate the vector

  //Close the reader and writer
  reader.Close();
  writer.Close();
}
