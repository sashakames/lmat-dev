include ../../../../Makefile.inc
#include /usr//local/tools/tau-2.17.1/lib/Makefile.tau-pdt-openmp-opari

DRIVERS =  test_KmerDB test_KmerDB_dump test_KmerDB_merge_dump test_error dump_ascii_kmers getKmersNaive test_offsets_2 test_offsets test_taxonomy test_taxonomy_2 test_taxonomy_3 test_taxonomy_4

all: $(DRIVERS)

testme : 
	-@echo "current CXXXFLAGS: $(CXXFLAGS)"
	-@echo
	-@echo "current INCLUDE: $(INCLUDE)"
	-@echo
	-@echo "current LIBS: $(LIBS)"

test_JellyReader : $(METAG_LIB)  test_JellyReader.cpp
	g++ -g -O3 -Wall -fopenmp -I/g/g10/hysom/perm-je-0.9.1/include -I. -I/g/g10/hysom/metag_repo_6_feb_2012/dave_h/src/kmerdb -I/g/g10/hysom/metag_repo_6_feb_2012/dave_h/include  -DWITH_PJMALLOC=0 -DUSE_GENOME_LOCATIONS=0    test_JellyReader.cpp /g/g10/hysom/metag_repo_6_feb_2012/third_party/jellyfish-1.1.2_modified_statistical/jellyfish/JellyReader.o /g/g10/hysom/metag_repo_6_feb_2012/dave_h/src/kmerdb/lib/libmetag.a    -lpthread  -o test_JellyReader /g/g10/hysom/metag_repo_6_feb_2012/third_party/jellyfish-1.1.2_modified_statistical/jellyfish/.libs/libjellyfish.a /g/g10/hysom/metag_repo_6_feb_2012/third_party/jellyfish-1.1.2_modified_statistical/jellyfish/kmer.cc

test_TaxTree_threaded : $(METAG_LIB)  test_TaxTree_threaded.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@


test_TaxTree : $(METAG_LIB)  test_TaxTree.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

test_TaxTableAndTaxHisto : $(METAG_LIB)  test_TaxTableAndTaxHisto.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

test_TaxTable : $(METAG_LIB)  test_TaxTable.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

test_GenomeIdToTaxId : $(METAG_LIB)  test_GenomeIdToTaxId.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

test_lca : $(METAG_LIB)  test_lca.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

test_taxonomy : $(METAG_LIB)  test_taxonomy.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

test_taxonomy_2 : $(METAG_LIB)  test_taxonomy_2.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

test_taxonomy_3 : $(METAG_LIB)  test_taxonomy_3.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

test_taxonomy_4 : $(METAG_LIB)  test_taxonomy_4.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

test_KmerDB_merge_dump : $(METAG_LIB)  test_KmerDB_merge_dump.cpp
	$(CXX) $(CXXFLAGS) test_KmerDB_merge_dump.cpp $(LIBS) -o $@

test_offsets : $(METAG_LIB)  test_offsets.cpp
	$(CXX) $(CXXFLAGS) test_offsets.cpp $(LIBS) -o $@

test_offsets_2 : $(METAG_LIB)  test_offsets_2.cpp
	$(CXX) $(CXXFLAGS) test_offsets_2.cpp $(LIBS) -o $@

test_KmerDB : $(METAG_LIB)  test_KmerDB.cpp
	$(CXX) $(CXXFLAGS) test_KmerDB.cpp $(LIBS) -o $@

test_KmerDB_dump : $(METAG_LIB)  test_KmerDB_dump.cpp
	$(CXX) $(CXXFLAGS) test_KmerDB_dump.cpp $(LIBS) -o $@

test_query_time : $(METAG_LIB)  test_query_time.cpp
	$(CXX) $(CXXFLAGS) test_query_time.cpp  $(LIBS) -o $@

test_error : $(METAG_LIB)  test_error.cpp
	$(CXX) $(CXXFLAGS) test_error.cpp  $(LIBS) -o $@

getKmersNaive : $(METAG_LIB)  getKmersNaive.cpp
	$(CXX) $(CXXFLAGS) getKmersNaive.cpp  $(LIBS) -o $@

dump_ascii_kmers : $(METAG_LIB)  dump_ascii_kmers.cpp
	$(CXX) $(CXXFLAGS) dump_ascii_kmers.cpp  $(LIBS) -o $@

clean: 
	rm -f *.o $(DRIVERS) *.mc *.out *TVD* tmp
