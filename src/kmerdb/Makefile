include ../../Makefile.inc

ifeq ($(USE_SORTED_DB),1)
INDEXDB=SortedDb
else
INDEXDB=TaxTable
endif


CXXHDR = $(INDEXDB).hpp KmerNode.hpp Utils.hpp KmerFileMetaData.hpp \
         TaxTree.hpp TaxNode.hpp 


CXXSRCS = $(INDEXDB).cpp Utils.cpp KmerFileMetaData.cpp TaxTree.cpp TaxNode.cpp 


CXXOBJS = $(INDEXDB).o Utils.o KmerFileMetaData.o TaxTree.o TaxNode.o 

lib: clean $(CXXHDR) $(CXXSRCS) $(CXXOBJS)
	ar -cvq libmetag.a $(CXXOBJS)
	mv libmetag.a lib


clean : 
	rm -f *.o tmp *.mc a.out *TVD*
	rm -f lib/libmetag.a
