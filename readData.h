// -------------------------------------------------------------------------- //
// Author  Ryan D. Crawford
// Date    02/14/2018
// Purpose Define a struct to store data on each informative read
// -------------------------------------------------------------------------- //

  struct readData
  {
    int readIndex; // Indicates the index of the read in the fastq file
    int bcNum;     // Indicates the barces conatined in the read pair
  };
