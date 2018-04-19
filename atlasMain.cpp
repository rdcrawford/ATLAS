#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <string.h>
#include <vector>
#include "constants.h"
#include "BarcodeClass.h"
#include "readData.h"
#include "UnmappedRead.h"
using namespace std;

// -------------------------------------------------------------------------- //
// Author  Ryan D. Crawford
// Date    02/14/2018
// Purpose This Program aligns TLA sequencing data using the Bamtools API and
//         finds informative read pairs where one read contaings a barcode
//         sequence. The opposite mate, if unmapped of duplicately mapped is
//         then digested and realigned to the reference. Finally, alignements
//         from the processed reads are then flagged with the barcode number
//         contained in the opposite read.
// -------------------------------------------------------------------------- //

// ---- Global Function Prototypes ------------------------------------------ //

// Input   test file in fasta fromat
// Return  vector containing the sequence in the fasta file
vector <Barcode> *parseFasta(ifstream &fastaFile);

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
vector <readData> *findInformativeReads(string bamFile, string baiFile,
                   string chr,  int startPos, int endPos,
                   vector <Barcode> *barcodes);

// Input       1) Bamfile
//             2) Sequence of the restriction site used to make the library
//             3) Vector containting the vectro of structs with information on
//                the informative reads
// Description This function takes the vector of reads that contains the
//             information on all of the informative reads and if reads are
//             unmapped or duplicately mapped the function digests them and
//             writes a new fastq file with the digrest products
void processInformativeReads(string &bamFile, string &restrSite,
     vector <readData> *informativeReads);

// Input       1) Vector containting the vectro of structs with information on
//                the informative reads
//             2) Bam file containing the alignemts for the informative reads
// Description This function takes the data on the informative reads and adds
//             a tag with the number of the barcode that the read belongs to
void tagBarcodeReads(vector <readData> *informativeReads, string bamFile);

// ---- Main Function ------------------------------------------------------- //

int main(int argc, char *argv[])
{
  // Variable declarations:
  string   reference;       // Reference genome (BWA indexed)
  string   fqFileRead1;     // Name of the BAM file for the reads
  string   fqFileRead2;     // Name of the Bam index file fot the reads
  string   barcodeFile;     // Name of the barcode file
  string   restrSite;       // String corresponding the ristriction site
  string   chrName;         // Name of the "chromosome" containing the bc
  string   bamFile;         // Name of the BAM file
  string   baiFile;         // Name of the BAM index file
  bool     areInvalidArgs;  // Indicates valid input parameters
  int      bcStart;         // Start of the locus containing the bc
  int      bcEnd;           // End of the locus containing the bc
  ifstream matePairOneFile; // First mate pair from the paired end data
  ifstream matePairTwoFile; // Second mate pair from the paired end data
  ifstream barcodeSeqs;     // Fasta formated Barcode sequence
  istringstream argStream;  // input sream for the command line arguments
  ostringstream oss;        // Out string stream for command line arguments

  // Parse the command line arguments
  for (int j = 1; j < argc; j++)
    oss << " " << argv[j];

  argStream.str(oss.str()); //Convert stream stream to c string

  //Asign args to variables
  argStream >> reference;
  argStream >> fqFileRead1;
  argStream >> fqFileRead2;
  argStream >> barcodeFile;
  argStream >> restrSite;
  argStream >> chrName;
  argStream >> bcStart;
  argStream >> bcEnd;

  //If there is an error with the command line agruments, print and error
  if (argc != COMMAND_ARGS_LEN) //Incorrect number of args...
  {
    cout << endl << endl <<
    "    Incorrect number of command line arguments!" << endl;
    areInvalidArgs = false;
  }
  else if(!argStream) //String stream fail state...
  {
    cout << endl << endl <<
    "    Invalid format for parameters!" << endl;
    areInvalidArgs = false;
  }
  else
  {
    areInvalidArgs = true;
  }

  if(!areInvalidArgs)
  {
    cout << endl << "  Usage:" << endl << endl
    << "    " << argv[0]
    << " <reference> <matePairOneFile.fq> <matePairTwoFile.fq> <Barcodeseq>"
    << " <Restriction Site> <Anchor chromosome> <Bacrode start> <Barcode end>"
    << endl                                       << endl
    << "      1) Reference (fasta format)"        << endl
    << "      2) Mate pair 1 (fastq format)"      << endl
    << "      3) Mate pair 2 (fastq format)"      << endl
    << "      4) Barcode Sequence (fasta format)" << endl
    << "      5) Restriction site sequence"       << endl
    << "      6) Anchor chromosome             "  << endl
    << "      7) Barcode start position        "  << endl
    << "      8) Barcode end position          "  << endl
    << endl;
    exit(2);
  }
  else //Output the parameters that were input
  {
    cout                                                            << endl
    << "    The input paramaters are:  "             << endl        << endl
    << "      1) Reference (fasta format):         " << reference   << endl
    << "      2) Mate Pair 1 (fastq format):       " << fqFileRead1 << endl
    << "      3) Mate pair 2 (fastq format):       " << fqFileRead2 << endl
    << "      4) Barcode sequences (fasta format): " << barcodeFile << endl
    << "      5) Restriction site sequence:        " << restrSite   << endl
    << "      6) Anchor chromosome                 " << chrName     << endl
    << "      7) Barcode start position:           " << bcStart     << endl
    << "      8) Barcode end position:             " << bcEnd       << endl
    << endl;
  }

  // Open the text files
  matePairOneFile.open(fqFileRead1.c_str());
  matePairTwoFile.open(fqFileRead1.c_str());
  barcodeSeqs.open(barcodeFile.c_str());

  if (barcodeSeqs.fail()) //If the barcode file failed to open, print an error
  {
    cout << "  Exiting: Unable to open input " << argv[4] << "!" << endl
    << endl;
    exit(1);
  }

  ostringstream indexString;                // Open the string stream
  string        indexCommand;               // String for the index command
  indexString.str("");                      // Clear the string stream
  indexString << "bwa index " << reference; // Fill the string stream
  indexCommand = indexString.str();         // Convert to a c string

  ostringstream bwaMemString;               // Open the string stream
  string        bwaMemCommand;              // String to store the command
  bwaMemString.str("");                     // Clear the string stream
  bamFile = "processedReads.bam";           // Name of the bam file to write 

  // Append arguments to the string stream
  bwaMemString << "bwa mem -U 0 -M " << reference << " " << fqFileRead1 << " "
               << fqFileRead2 << " | samtools view -bS - > " << bamFile;
  bwaMemCommand = bwaMemString.str(); // Convert to a c string

  cout << endl << "    Step 1 of 7: Indexing reference genome" << endl << endl;
  system(indexCommand.c_str());  // Index the reerence genome

  cout << endl << "    Step 2 of 7: Aligning reads with BWA mem"<< endl << endl;
  system(bwaMemCommand.c_str()); // Align reads with BWA and write to BAM file

  cout << endl << "   Step 3 of 7: Parsing barcode sequences" << endl << endl;
  vector <Barcode> *barcodes;         // Initiate a vector for barcodes
  barcodes = parseFasta(barcodeSeqs); // Generate a vector of barcode objects
  barcodeSeqs.close();                // close the file stream for the barcodes

  // Create a pointer to the vector containing the informative reads
  vector <readData> *informativeReads;

  cout << "   Step 4 of 7: Finding informative reads"    << endl << endl
       << "     - Searching for " << barcodes->size() << " barcodes"
       << endl << endl;

  // Create a dynamically allocated vector with all of the informativeReads
  informativeReads = findInformativeReads(bamFile, baiFile, chrName, bcStart,
                     bcEnd, barcodes);

  cout << "     - Identified " << informativeReads->size()
       << " informative reads" << endl << endl;

  //Write out the data for all of the barcodes to a text file
  ofstream barcodeData;
  barcodeData.open("barcodeData.tsv");
  for (int i = 0; i < barcodes->size(); i ++)
  {
    barcodeData << i                              << "\t"
                << (*barcodes)[i].getBarcodeSeq() << "\t"
                << (*barcodes)[i].getReadCount()  << endl;
  }
  barcodeData.close();

  // Deallocate the barcodes on the heap
  delete barcodes;

  // Digests the mate of informative read pairs
  cout << "   Step 5 of 7: Processing informative mate pairs" << endl << endl;
  processInformativeReads(bamFile, restrSite, informativeReads);

  cout << "   Step 6 of 7: Aligning processed reads" << endl << endl;
  string        realignBamFile = "processedReads.bam";
  ostringstream bwaMemReAlgnStr;    // Open the string stream
  string        bwaMemReAlgnStrCmd; // String to store the command
  bwaMemReAlgnStr.str("");          // Clear the string stream

  // Append arguments to the string stream
  bwaMemReAlgnStr << "bwa mem " << reference <<
  " processedReads.fq  | samtools view -bS - > " << realignBamFile;
  bwaMemReAlgnStrCmd = bwaMemReAlgnStr.str();
  cout << bwaMemReAlgnStrCmd << endl;

  // Align reads with BWA and write to BAM file
  system(bwaMemReAlgnStrCmd.c_str());

  // Add flags indicating which barcode each read is paired with
  cout << endl << endl << "   Step 7 of 7: Adding tags to BAM file and sorting"
       << " reads by tag"<< endl
       << endl;

  //Add tags to the bam file indicating which barcode is contained in the
  //opposite read
  tagBarcodeReads(informativeReads, realignBamFile);

  delete informativeReads; //Deallocate memeory used for the informativeReads

  return 0;
}

// -------------------------------------------------------------------------- //
