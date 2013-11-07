include ../../Makefile.inc

.PHONY: all
all: $(APPS)

.PHONY: clean
clean: 
	$(RM) $(OBJECTS) $(APPS)

doit :

countCanonicalKmers : $(METAG_LIB) countCanonicalKmers.cpp doit
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

countCanonicalKmers_scott : $(METAG_LIB) countCanonicalKmers_scott.cpp doit
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

countCanonicalKmers_scott2 : $(METAG_LIB) countCanonicalKmers_scott2.cpp doit
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

dumpBinaryKmers : $(METAG_LIB) dumpBinaryKmers.cpp  doit
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

test_Encoder : $(METAG_LIB) test_Encoder.cpp doit
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

test_Encoder2 : $(METAG_LIB) test_Encoder2.cpp doit
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@

genCanonicalKmers : $(METAG_LIB) genCanonicalKmers.cpp doit
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@
