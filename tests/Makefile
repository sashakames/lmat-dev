include ../Makefile.inc

test_KmerDB : $(METAG_LIB) test_KmerDB.cpp
#test_KmerDB : ../lib/libmetag.so test_KmerDB.cpp
	$(CXX) $(LDFLAGS) test_KmerDB.cpp $(LIBS) -o test_KmerDB

test_KmerDB_2 : $(METAG_LIB) test_KmerDB_2.cpp
#test_KmerDB : ../lib/libmetag.so test_KmerDB.cpp
	$(CXX) $(LDFLAGS) test_KmerDB_2.cpp $(LIBS) -o test_KmerDB_2

test_KmerNode : ../lib/libmetag.so test_KmerNode.cpp
	$(CXX) $(LDFLAGS) test_KmerNode.cpp $(LIBS) -o test_KmerNode

test_system : $(METAG_LIB) test_system.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) test_system.cpp $(LIBS) -o test_system

test_persistent_system : $(METAG_LIB) test_persistent_system.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) test_persistent_system.cpp $(LIBS) -o test_persistent_system

test_promotion : $(METAG_LIB) test_promotion.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) test_promotion.cpp $(LIBS) -o test_promotion

test_find_children: $(METAG_LIB) test_find_children.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) test_find_children.cpp $(LIBS) -o test_find_children

test_read_search: $(METAG_LIB) test_read_search.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) test_read_search.cpp $(LIBS) -o test_read_search

gen_loc_histo: $(METAG_LIB) gen_loc_histo.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) gen_loc_histo.cpp $(LIBS) -o gen_loc_histo

clean: 
	rm -f *.o
