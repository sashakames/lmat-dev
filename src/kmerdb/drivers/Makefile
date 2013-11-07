include ../../../Makefile.inc

ifeq ($(LOCATION_FLAG), USE_GENOME_LOCATIONS=0)
  ALL = sortKmerNodes_no_locations dumpKmerNodes_no_locations printMetaData verifyKmerDBFile
else 
  ALL = sortKmerNodes dumpKmerNodesWithLocations fixOffsetError printMetaData
endif


DRIVERS = sortKmerNodes_no_locations dumpKmerNodes_no_locations sortKmerNodes dumpKmerNodesWithLocations fixOffsetError 


all : $(ALL)

verifyKmerDBFile : ../lib/libmetag.a verifyKmerDBFile.cpp
	$(CXX) $(CXXFLAGS) verifyKmerDBFile.cpp $(LIBS) -o $@

printMetaData : ../lib/libmetag.a printMetaData.cpp
	$(CXX) $(CXXFLAGS) printMetaData.cpp $(LIBS) -o $@

sortKmerNodes : ../lib/libmetag.a sortKmerNodes.cpp
	$(CXX) $(CXXFLAGS) sortKmerNodes.cpp $(LIBS) -o $@

sortKmerNodes_no_locations : ../lib/libmetag.a sortKmerNodes.cpp
	$(CXX) $(CXXFLAGS) sortKmerNodes.cpp $(LIBS) -o $@

dumpKmerNodesWithLocations : ../lib/libmetag.a dumpKmerNodes.cpp
	$(CXX) $(CXXFLAGS) dumpKmerNodes.cpp $(LIBS) -o $@

dumpKmerNodes_no_locations: ../lib/libmetag.a dumpKmerNodes.cpp
	$(CXX) $(CXXFLAGS) dumpKmerNodes.cpp $(LIBS) -o $@

clean: 
	rm -f *.o $(ALL) $(DRIVERS) *.mc
