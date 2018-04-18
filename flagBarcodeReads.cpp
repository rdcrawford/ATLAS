#include <fstream>
#include <algorithm>
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
// Purpose This function for ATLAS uses the bamtools C++ API to add flags to
//         the BAM file containing the processed reads to indicate which barcode
//         the restriction fragment is paired to. The BAM file is then sorted
//         by flag.
// -------------------------------------------------------------------------- //

// Input
// Return
// Description
void tagBarcodeReads(vector <readData> *informativeReads, string bamFile)
{
  string       outputFilename = bamFile;
  BamReader    reader;
  BamWriter    writer;
  BamAlignment algn;
  int          readNum = 0;
  bool         isOpen;
  string       prevAlgnName;
  string       tag  = "XX";
  string       type = "i";

  // attempt to open our BamMultiReader
  isOpen = reader.Open(bamFile);
  if (!isOpen)
  {
    cout << "Could not open input BAM files" << endl;
    exit(1);
  }

  // Retrieve metadata from BAM files, these are required by BamWriter
  const SamHeader HEADER     = reader.GetHeader();
  const RefVector REFERENCES = reader.GetReferenceData();

  // Create a new vector for the dynamically allovated alignments
  std::vector<BamAlignment> *algnVec;
  algnVec = new std::vector <BamAlignment>;

  // iterate through all alignments, only keeping ones with high map quality
  while (reader.GetNextAlignmentCore(algn))
  {
    //Add a tag to the bam file so it can be sorted later
    bool isUpdated = algn.AddTag(tag, type, readNum);

    if (algn.MapQuality >= 10) //Only keep alignments with high map quality
    {
      algnVec->push_back(algn);
    }

    //If this is not a restriction fragment from the previously processed read,
    //
    if (prevAlgnName != algn.Name)
    {
      readNum ++; // Increment the alignment number

      //Change the name of the variable to the that of the alignment that was
      //just processed
      prevAlgnName = algn.Name;
    }
  }

  // Sort by the barcode tag defined for the barcode num
  std::sort(algnVec->begin(), algnVec->end(), Algorithms::Sort::ByTag<int>(tag));

  int currenctBcTag = 0; //Barcode currently being written to the BAM file
  int readTag       = 0; //Take of the current read being iterated over
  readNum           = 0; //Number of the read currently being iterated over


  ostringstream sortStream;                 // Open the string stream
  string        vcCallString;                // String to store the command
  vcCallStream.str("");                       // Clear the string stream
  // Append arguments to the string stream
  vcCallStream << "bwa mem -U 0 -M " << reference << " " << fqFileRead1 << " "
               << fqFileRead2 << " | samtools view -bS - > " << bamFile;
  bwaMemCommand = bwaMemString.str(); // Convert to a c string

  ostringstream vcCallStream;                 // Open the string stream
  string        vcCallString;                // String to store the command
  vcCallStream.str("");                       // Clear the string stream
  // Append arguments to the string stream
  vcCallStream << "bwa mem -U 0 -M " << reference << " " << fqFileRead1 << " "
               << fqFileRead2 << " | samtools view -bS - > " << bamFile;
  bwaMemCommand = bwaMemString.str(); // Convert to a c string

  ostringstream bwaMemString;                 // Open the string stream
  string        bwaMemCommand;                // String to store the command
  bwaMemString.str("");                       // Clear the string stream
  bamFile = "firstPass.bam";

  // Append arguments to the string stream
  bwaMemString << "bwa mem -U 0 -M " << reference << " " << fqFileRead1 << " "
               << fqFileRead2 << " | samtools view -bS - > " << bamFile;
  bwaMemCommand = bwaMemString.str(); // Convert to a c string

  //Attempt to open our BamWriter
  writer.Open(tempBam, HEADER, REFERENCES);

  //For each alignement in the vector write to the output file
  for(int i = 0; i < algnVec->size(); i ++)
  {
    writer.SaveAlignment((*algnVec)[i]); //Save the current alignment
  }

  delete algnVec; //Deallocate the vector

  // Close the reader and writer
  reader.Close();
  writer.Close();
}
