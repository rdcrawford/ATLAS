#include <string>
#include <algorithm>
#include <vector>
#include <list>
#include <fstream>
#include "UnmappedRead.h"
#include "SortedListClass.h"
using namespace std;

// -------------------------------------------------------------------------- //
// Author:  Ryan D. Crawford
// Date:    02/17/2018
// Purpose: This file defines an "UnmappedRead" Seq class object to store the
//          attributes of an entry in a fastq file.
// -------------------------------------------------------------------------- //

// value Ctor
UnmappedRead::UnmappedRead(string inName, string inSeq, string inQualScores)
{
  name       = inName;       // The name of a sequence
  seq        = inSeq;        // Raw sequence read
  qualScores = inQualScores; // The quality scores for the fq read
}

// Input: 1) A fastq read which is the opposite mate pair of a read containing a
//           Barcode
//        2) A restriction site to use to digest the fastq Sequence
// Return: True if a restriction site was found in the read.
// Description: The sequences between the identified restricion sites are output
//              into a vector of strings.
bool UnmappedRead::digestSeq(string restrSite, ofstream &newFastqFile)
{
  // get the length of the restriction site and read
  int    restrSiteLen  = restrSite.length();
  int    readLength    = seq.length();
  int    startPos;
  int    endPos;
  int    restrSiteCount = restrSiteLoci.getNumElems();
  string restrFragSeq;
  string restrFragQuals;
  string restrFragName;

  // Find the restriction sites in the current sequence
  findResSites(restrSite);

  if(restrSiteLoci.getNumElems() == 0)
  {
    writeComponentValues(newFastqFile);
    return false;
  }
  else if (restrSiteLoci.getNumElems() == 1 && restrSiteLoci[0] == 0)
  {
    writeComponentValues(newFastqFile);
  }
  else if (restrSiteLoci.getNumElems() == 1 && restrSiteLoci[0] == seq.length() -4)
  {
    writeComponentValues(newFastqFile);
  }
  else
  {
    for (int i = 0; i <= restrSiteLoci.getNumElems(); i++)
    {
      if (i == 0)
      {
        startPos = 0;
        endPos   = restrSiteLoci[i] - startPos + restrSite.length();
      }
      else if (i == restrSiteLoci.getNumElems())
      {
        startPos = restrSiteLoci[i - 1];
        endPos   = seq.length() - startPos + restrSite.length();
      }
      else
      {
        startPos = restrSiteLoci[i - 1];
        endPos   = restrSiteLoci[i] - startPos + restrSite.length();
      }

      // Write the restriction fragments to the output file
      newFastqFile << '@'
                   << name                                << endl
                   << seq.substr(startPos, endPos)        << endl
                   << '+'                                 << endl
                   << qualScores.substr(startPos, endPos) << endl;
    }
  }
  return true;
}

// Input: 1) A sequence opposite a mate pair containing a Barcode.
//        2) Restriction site sequence
//        3) Count of restriction sites within the read
//        4) Location of restriction sites within the read
// Return: True if a restriction site was found.
// Description:
bool UnmappedRead::findResSites(string restrSite)
{
  int resSiteLocus = -1; // start the value of locus at 0 to enter the loop
  // if a restriction site is found add it to the vector and add one to the
  // count of restriction sites.
  do
  {
    resSiteLocus = resSiteLocus + 1;
    resSiteLocus = seq.find(restrSite, resSiteLocus);
    if (resSiteLocus != -1)
    {
      restrSiteLoci.insertValue(resSiteLocus);
    }
  } while (resSiteLocus != -1);

  restrSite = getRevComp(restrSite);
  resSiteLocus = -1; // start the value of locus at 0 to enter the loop

  // if a restriction site is found add it to the vector and add one to the
  // count of restriction sites.
  do
  {
    resSiteLocus = resSiteLocus + 1;
    resSiteLocus = seq.find(restrSite, resSiteLocus);
    if (resSiteLocus != -1)
    {
      restrSiteLoci.insertValue(resSiteLocus);
    }
  } while (resSiteLocus != -1);

  if (restrSiteLoci.getNumElems() == 0)
  {
    return false;
  }
  else
  {
    return true;
  }
}

// Write the values of the digested sequence to the output
void UnmappedRead::writeComponentValues(ofstream &newFastqFile)
{
  //outFile
    newFastqFile << '@'
                 << name       << endl
                 << seq        << endl
                 << "+"        << endl
                 << qualScores << endl;
}

// Getter function to return the number of restriction sites
int UnmappedRead::getResSiteCount()
{
  return restrSiteLoci.getNumElems();
}

// Calculate the reverse compliment of string
string UnmappedRead::getRevComp(string restrSite)
{
  string revCompSeq = restrSite;
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
