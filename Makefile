CC=g++
IDIR =/usr/local/include/bamtools/
export LD_LIBRARY_PATH=/usr/local/lib:LD_LIBRARY_PATH
LINKER=bamtools
BAMTOOLS_LIB_DIR=/usr/local/lib/
CPPFLAGS=-I$(IDIR)
CXXFLAGS=-Wl,-rpath,/usr/local/lib/
DEPS = BarcodeClass.h UnmappedRead.h constants.h readData.h \
	SortedListClass.h LinkedNodeClass.h
OBJ = atlasMain.o processInformativeReadPairs.o BarcodeClass.o UnmappedRead.o \
  SortedListClass.o findInformativeReads.o parseFasta.o LinkedNodeClass.o \
	tagBarcodeReads.o

%.o: %.c $(DEPS)
	$(CC) -std=gnu++0x -c -o $@ $< $(CPPFLAGS) $(CXXFLAGS) -L$(LD_LIBRARY_PATH) \
	-lbamtools
ATLAS: $(OBJ)
	$(CC) -std=gnu++0x -o $@ $^ $(CPPFLAGS) $(CXXFLAGS) -L$(LD_LIBRARY_PATH) \
	-lbamtools
clean:
	rm -f $(OBJ) ATLAS
