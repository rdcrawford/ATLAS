#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include "SortedListClass.h"
using namespace std;

// -------------------------------------------------------------------------- //
// Author:  Ryan D. Crawford
// Date:    02/17/2018
// Purpose: This file defines an "UnmappedRead" class object to store the
//          attributes of an entry in a fastq file.
// -------------------------------------------------------------------------- //

#ifndef _UNMAPPEDREAD_H_
#define _UNMAPPEDREAD_H_
class UnmappedRead
{
  public:
    // Default Ctor
    UnmappedRead()
    { ; }

    // Input: 1:3) attributes of a fastq entry
    // Description: Value Ctor
    UnmappedRead(string inName, string inSeq, string inQualScores);

    // Input: Output stream to write reads to
    // Write the values of the digested sequence to the output
    void writeComponentValues(ofstream &newFastqFile);

    // Input: String corresponding to a restriction site
    // Return: Vector of digested fastq files
    // Destription:
    bool digestSeq(string restrSite, ofstream &outFile);

  private:
    string name;                   // The name of a sequence
    string seq;                    // Raw sequence read
    string qualScores;             // The quality scores for the fq read
    SortedListClass restrSiteLoci; // list of restriction site positions

    // Calculate the reverse compliment of string
    string getRevComp(string restrSite);

    // Input: 1) A sequence opposite a mate pair containing a Barcode.
    //        2) Restriction site sequence
    //        3) Count coressponding to the number of restriction sites within
    //           the read.
    //        4) Location of restriction sites within the read
    // Return: True if a restriction site was found.
    // Description: For a given read all of the restiction sites within the
    //              read are found. The count of restricion sites are then
    //              updated.
    bool findResSites(string restrSite);

    // Getter function to return the number of restriction sites
    int getResSiteCount();
};
#endif
