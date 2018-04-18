#include <string>
using namespace std;

// ---- Constant Declarations and Initializations --------------------------- //

#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

// Number of command line arguments required for this program
const int COMMAND_ARGS_LEN = 9;

// Threshold for the quality of the mapping of the reads
const int MIN_MAP_QUALITY = 0;

const string TAG      = "XB";  //Bag under which the bacrode number is added
const string TAG_TYPE = "i";   //Type of tag used by "tag"

const string SA_TAG = "SA"; //String for split reads

#endif
