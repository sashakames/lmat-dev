include ../Makefile.inc


.PHONY: all
all: $(APPS)

.PHONY: clean
clean: 
	$(RM) $(OBJECTS) $(APPS)

all : convertTaxHistoTo16Bits convertTaxHisto_ascii2bin convertTaxHisto_bin2ascii countTaxHistoTidFrequency dumpTaxHistoOutput dumpTaxHistoOutput_new_fmt kmerPrefixCounter kmerPrefixCounterSort mergeKmerPrefixOutputs readProfiler run_tax_histo_MPI sortKmerNodes tax_histo_new_fmt verifyKmerEntries verifyTaxHistoEntries


convertTaxHistoToAscii : $(METAG_LIB) convertTaxHistoToAscii.cpp
	mpiCC $(CXXFLAGS) $@.cpp $(LIBS) -o $@

getAllTaxNodes : $(METAG_LIB) getAllTaxNodes.cpp
	mpiCC $(CXXFLAGS) $@.cpp $(LIBS) -o $@

countTaxHistoTidFrequency : $(METAG_LIB) countTaxHistoTidFrequency.cpp
	mpiCC $(CXXFLAGS) $@.cpp $(LIBS) -o $@

verifyKmerEntries : $(METAG_LIB) verifyKmerEntries.cpp
	mpiCC $(CXXFLAGS) $@.cpp $(LIBS) -o $@


convertTaxHistoTo16Bits : $(METAG_LIB) convertTaxHistoTo16Bits.cpp
	mpiCC $(CXXFLAGS) $@.cpp $(LIBS) -o $@

xx : $(METAG_LIB) xx.cpp
	mpiCC $(CXXFLAGS) $@.cpp $(LIBS) -o $@

test : $(METAG_LIB) test.cpp
	mpiCC $(CXXFLAGS) $@.cpp $(LIBS) -o $@

readProfiler : $(METAG_LIB) readProfiler.cpp
	mpiCC $(CXXFLAGS) $@.cpp $(LIBS) -o $@

testEncoder : $(METAG_LIB) testEncoder.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

testPulledAgainstGenomes : $(METAG_LIB) testPulledAgainstGenomes.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

newKmer : $(METAG_LIB) newKmer.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

find : $(METAG_LIB) find.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

countHumanInTaxHisto : $(METAG_LIB) countHumanInTaxHisto.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

countHumanInKmerDBFiles : $(METAG_LIB) countHumanInKmerDBFiles.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

doit : $(METAG_LIB) doit.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@


dumpBinaryKmers : $(METAG_LIB) dumpBinaryKmers.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

hunt : $(METAG_LIB) hunt.cpp
hunt : $(METAG_LIB) hunt.cpp

hunt : $(METAG_LIB) hunt.cpp
hunt : $(METAG_LIB) hunt.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

partitionTaxHistoBinFile : $(METAG_LIB) partitionTaxHistoBinFile.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

run_tax_histo_MPI : $(METAG_LIB) run_tax_histo_MPI.cpp
	mpiCC $(CXXFLAGS) $@.cpp $(LIBS) -o $@

tax_histo_omp : $(METAG_LIB) tax_histo_omp.cpp
	mpiCC $(CXXFLAGS) $@.cpp $(LIBS) -o $@

makeSingleTidDB : $(METAG_LIB) makeSingleTidDB.cpp
	mpiCC $(CXXFLAGS) $@.cpp $(LIBS) -o $@

tax_histo_fast_limited : $(METAG_LIB) tax_histo_fast_limited.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

countTaxHistoEntries : $(METAG_LIB) countTaxHistoEntries.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

makeMergeList : $(METAG_LIB) makeMergeList.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

testTaxHistoForDuplicates : $(METAG_LIB) testTaxHistoForDuplicates.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

dumpKmers : $(METAG_LIB) dumpKmers.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

dumpTaxHistoOutput : $(METAG_LIB) dumpTaxHistoOutput.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

convertTaxHisto_bin2ascii : $(METAG_LIB) convertTaxHisto_bin2ascii.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

kmerPrefixCounter : $(METAG_LIB) kmerPrefixCounter.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

kmerPrefixCounter_MPI : $(METAG_LIB) kmerPrefixCounter_MPI.cpp
	mpiCC $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

kmerPrefixCounterSort : $(METAG_LIB) kmerPrefixCounterSort.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

kmerPrefixCounterSort_V2 : $(METAG_LIB) kmerPrefixCounterSort_V2.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

kmerPrefixCounterSort_withLocations : $(METAG_LIB) kmerPrefixCounterSort_withLocations.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

countKmers_tax_histo : $(METAG_LIB) countKmers_tax_histo.cpp
	mpiCC $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

tax_histo_new_fmt : $(METAG_LIB) tax_histo_new_fmt.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

dumpTaxHistoOutput_new_fmt : $(METAG_LIB) dumpTaxHistoOutput_new_fmt.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

dumpJellylistOutput : $(METAG_LIB) dumpJellylistOutput.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

convertTaxHisto_ascii2bin : $(METAG_LIB) convertTaxHisto_ascii2bin.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

