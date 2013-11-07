include ../../../../Makefile.inc
#include /usr//local/tools/tau-2.17.1/lib/Makefile.tau-pdt-openmp-opari

DRIVERS =  ex1 ex2

all: $(DRIVERS)

#test_query_time

ex1 : $(METAG_LIB)  ex1.cpp
	$(CXX) $(CXXFLAGS) ex1.cpp $(LIBS) -o $@

ex2 : $(METAG_LIB)  ex2.cpp
	$(CXX) $(CXXFLAGS) ex2.cpp $(LIBS) -o $@

ex3 : $(METAG_LIB)  ex3.cpp
	$(CXX) $(CXXFLAGS) ex3.cpp $(LIBS) -o $@

ex4 : $(METAG_LIB)  ex4.cpp
	$(CXX) $(CXXFLAGS) ex4.cpp $(LIBS) -o $@

ex5 : $(METAG_LIB)  ex5.cpp
	$(CXX) $(CXXFLAGS) ex5.cpp $(LIBS) -o $@

ex6 : $(METAG_LIB)  ex6.cpp
	$(CXX) $(CXXFLAGS) $@.cpp $(LIBS) -o $@


clean: 
	rm -f *.o $(DRIVERS)
