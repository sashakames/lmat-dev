include ../../Makefile.inc

.PHONY: all
all: $(APPS)

.PHONY: clean
clean: 
	$(RM) $(OBJECTS) $(APPS)

.PHONY: test

test_for_duplicate_tax_histo_IDs: $(METAG_LIB) test_for_duplicate_tax_histo_IDs.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@


testTaxHisto: $(METAG_LIB) testTaxHisto.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@
