include ../Makefile.inc


APPS = countTaxidFrequency

OBJECTS = $(addsuffix .o,$(APPS))
DATA = $(METAG_DIR)/dev/src/kmerdb/examples/tests/data

.PHONY: all
all: $(APPS)

.PHONY: clean
clean: 
	$(RM) $(OBJECTS) $(APPS)

countTaxidFrequency: $(METAG_LIB) countTaxidFrequency.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@


