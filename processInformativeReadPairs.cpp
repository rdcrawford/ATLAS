#include <fstream>
#include <string>
#include "api/BamAlignment.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "UnmappedRead.h"
#include "constants.h"
#include "readData.h"
using namespace BamTools;

// -------------------------------------------------------------------------- //
// Author  Ryan D. Crawford
// Date    04/3/2018
// Purpose Step five of ATLAS. Digests the mate of informative read pairs
// -------------------------------------------------------------------------- //

// Input       1) Bamfile
//             2) Sequence of the restriction site used to make the library
//             3) Vector containting the vectro of structs with information on
//                the informative reads
// Description This function takes the vector of reads that contains the
//             information on all of the informative reads and if reads are
//             unmapped or duplicately mapped the function digests them and
//             writes a new fastq file with the digrest products
void processInformativeReads(std::string &bamFile, std::string &restrSite,
     vector <readData> *informativeReads)
{
  // Variable Initializations:
  BamReader    reader;          //open the BAM reader
  BamWriter    writer;          //Open the BAM writer
  BamAlignment algn;            //Store the alignment
  bool         isOpen;          //Indicates wheter the bam file was opened
  int          count = 0;       //Counter for the read number
  int          infoReadNum = 0; //Stores counter value if a read is infromative
  ofstream     newFastqFile;    //Fastq file being written to
  string       saString;        //Stores the informationn on split reads

  isOpen = reader.Open(bamFile); // Open the bam file

  if (!isOpen) // If the BAM file wasn't opened properly print an error and exit
  {
    std::cout << "  Unable to open BAM file: " << bamFile << std::endl;
    exit(1);
  }

  //Open the output file stream
  newFastqFile.open("processedReads.fq");

  std::string alignedBases; // string to store the sequence of alignment

  //Iterate through alignments in this region,
  //ignoring alignments with a Map quality below the minimium below some cutoff
  while (reader.GetNextAlignment(algn) &&
         infoReadNum < informativeReads->size())
  {
    if (count == (*informativeReads)[infoReadNum].readIndex)
    {
      infoReadNum ++;

      saString = ""; //Empty the split read string

      //Get the SA string from the read
      algn.GetTag(SA_TAG, saString);

      //Create a new unmapped read object to store the attributes of the fastq
      //Sequence
      UnmappedRead currentRead(algn.Name, algn.QueryBases, algn.Qualities);

      if (saString.length() > 1 || !algn.IsMapped()) //If the read is split...
      {
        //Digest the read, and write the fragments
        currentRead.digestSeq(restrSite, newFastqFile);
      }
      else
      {
        //Write the whole sequence of the read
        currentRead.writeComponentValues(newFastqFile);
      }
    }
    count ++;
  }
  // close the reader
  reader.Close();
  newFastqFile.close();
}
