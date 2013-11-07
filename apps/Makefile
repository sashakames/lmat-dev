include ../Makefile.inc

ifeq ($(USE_OMP), 0)
READ_SEARCH_SRC=read_search
else
READ_SEARCH_SRC=read_search-tt
endif

ifeq ($(NO_WRITE), 1)
NWFLAG=-D_NO_WRITE
endif


APPS = read_search read_label tax_histo_fast_limited 

OBJECTS = $(addsuffix .o,$(APPS))
DATA = $(METAG_DIR)/dev/src/kmerdb/examples/tests/data

.PHONY: all
all: $(APPS)

.PHONY: clean
clean: 
	$(RM) $(OBJECTS) $(APPS)

.PHONY: test
test:
	./read_search -t 4 -d $(DATA)/test.fa.int.bin.no_loc.12 -i $(DATA)/reads.fa -o $(DATA)


read_search: $(METAG_LIB) $(READ_SEARCH_SRC).cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(READ_SEARCH_SRC).cpp $(LIBS) $(NWFLAG) -o $@

force_search: $(METAG_LIB) force_search.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS)   -o $@

intersect_with_lib: $(METAG_LIB) intersect_with_lib.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

kmerSearch : $(METAG_LIB) kmerSearch.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

kmer-search-th: $(METAG_LIB) kmer-search-th.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

gene_label: $(METAG_LIB) gene_label.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -DUSE_SORTED_DB=$(USE_SORTED_DB) -o $@

gene_label_v2: $(METAG_LIB) gene_label_v2.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -DUSE_SORTED_DB=$(USE_SORTED_DB) -o $@

content_summ: $(METAG_LIB) content_summ.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -DUSE_SORTED_DB=$(USE_SORTED_DB) -o $@


rand_read_label_v2: $(METAG_LIB) rand_read_label_v2.cpp
	$(CXX) $(CXXFLAGS) -std=c++0x $(LDFLAGS) $@.cpp $(LIBS) -DUSE_SORTED_DB=$(USE_SORTED_DB) -o $@

read_label: $(METAG_LIB) read_label.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -DUSE_SORTED_DB=$(USE_SORTED_DB) -o $@


rand_read_label: $(METAG_LIB) rand_read_label.cpp
	$(CXX) $(CXXFLAGS) -std=c++0x $(LDFLAGS) $@.cpp $(LIBS) -DUSE_SORTED_DB=$(USE_SORTED_DB) -o $@

read_label_v4: $(METAG_LIB) read_label_v4.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -DUSE_SORTED_DB=$(USE_SORTED_DB) -o $@

read_label_v3: $(METAG_LIB) read_label_v3.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -DUSE_SORTED_DB=$(USE_SORTED_DB) -o $@


read_label_v2: $(METAG_LIB) read_label_v2.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -DUSE_SORTED_DB=$(USE_SORTED_DB) -o $@

read_label_lo: $(METAG_LIB) read_label_lo.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS)  -o $@

read_lookup: $(METAG_LIB) read_lookup.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS)  -o $@

read_label_mine : $(METAG_LIB) read_label_mine.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@


read_label_sparse: $(METAG_LIB) read_label_sparse.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

read_label_sparse_static: $(METAG_LIB) read_label_sparse_static.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

read_lookup_sdb: $(METAG_LIB) read_lookup_sdb.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

node_trimmer : $(METAG_LIB) node_trimmer.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

pull_kmer_ids : $(METAG_LIB) pull_kmer_ids.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

get_kmer_data: $(METAG_LIB) get_kmer_data.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

buildTaxHistoFromDBTable: $(METAG_LIB) buildTaxHistoFromDBTable.cpp
	echo "$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@"
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

getKmersFromTaxHisto : $(METAG_LIB) getKmersFromTaxHisto.cpp
	echo "$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@"
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@


verifyJellyfishTaxHisto : $(METAG_LIB) verifyJellyfishTaxHisto.cpp
	echo "$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@"
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

dumpKmersFromDBTable : $(METAG_LIB) dumpKmersFromDBTable.cpp
	echo "$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@"
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

jellyfish_verifier : $(METAG_LIB) jellyfish_verifier.cpp
	echo "$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@"
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

get_kmer_data_sdb: $(METAG_LIB) get_kmer_data_sdb.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@


fasta.o: fasta.c
	gcc -c $(CXXFLAGS) fasta.c 
