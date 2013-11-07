include ../../Makefile.inc

APPS = make_db_table  read_db make_db
OBJECTS = $(addsuffix .o,$(APPS))
DATA = $(METAG_DIR)/dev/src/kmerdb/examples/tests/data

.PHONY: all
all: $(APPS)

.PHONY: clean
clean:
	$(RM) $(OBJECTS) $(APPS) *.dump

.PHONY: test
test:
	./make_db -i $(DATA)/test.fa.int.bin.no_loc.12 -o $(DATA)/test.fa.mmap
	./read_db -i $(DATA)/test.fa.mmap

#$(APPS): $(METAG_LIB) $(OBJECTS)

db_fix: $(METAG_LIB) db_fix.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

make_db: $(METAG_LIB) make_db.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

make_db_table: $(METAG_LIB) make_db_table.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

make_db_sorted: $(METAG_LIB) make_db_sorted.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

make_db_new: $(METAG_LIB) make_db_new.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

read_db: $(METAG_LIB) read_db.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

read_db_new: $(METAG_LIB) read_db_new.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

amake_db: $(METAG_LIB) amake_db.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

aread_db: $(METAG_LIB) aread_db.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@

fasta_to_db: $(METAG_LIB) fasta_to_db.cpp
	$(CXX)	$(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) ../fasta.o -o $@
