CXX = g++
INCLUDE = -I$(PJPATH)/include -I. -I$(METAG_DIR)/src/kmerdb -I$(METAG_DIR)/include 


CXXFLAGS = -O3 -g  -Wall -fopenmp $(INCLUDE) -D$(PJOPT)  \
            -DIDX_CONFIG=$(SDBIDX) -DTID_SIZE=$(TAXID_SIZE)  -DDBTID_T=$(TID_T) -Wno-deprecated \
	    -DUSE_SORTED_DB=$(USE_SORTED_DB)

METAG_LIB = $(METAG_DIR)/src/kmerdb/lib/libmetag.a

#note: libgzstream.a *must* appear before libz.a
LIBS = $(METAG_LIB)  $(PJMALLOC_LIB) -lpthread \
       $(METAG_DIR)/lib/libgzstream.a -ldl \
       $(METAG_DIR)/lib/libz.a 

LDENV= LD_RUN_PATH=/lib:/usr/lib:$(HOME)/lib:/usr/local/lib:$(PJPATH)/lib
RUNENV= LD_LIBRARY_PATH=/lib:/usr/lib:$(HOME)/lib:/usr/local/lib:$(PJPATH)/lib

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $< -o $@


#--------------------------------------------------------------------
# macros for conditional compilation
#--------------------------------------------------------------------

#compile and link with PJMALLOC?
ifeq ($(BOOST), 1)
 PJOPT = USE_BOOST=1 -DWITH_PJMALLOC=0
else
ifeq ($(WITH_PJMALLOC), 1)
 PJOPT = WITH_PJMALLOC=1
 PJMALLOC_LIB =  -L$(PJPATH)/lib -ljemalloc
else 
 PJOPT = WITH_PJMALLOC=0
endif
endif

ifndef USE_SORTED_DB
USE_SORTED_DB=1	
endif


#location of PJMALLOC lib
ifndef PJPATH
PJPATH= $(HOME)
endif

ifndef SDBIDX
SDBIDX=2027
endif

#default METAG_DIR
ifndef METAG_DIR
METAG_DIR=$(HOME)/lmat-git
endif

ifndef TAXID_SIZE
TID_T = uint32_t
TAXID_SIZE=32
else
ifeq ($(TAXID_SIZE), 16)
TID_T = uint16_t
else
TID_T = uint32_t
endif
endif

#--------------------------------------------------------------------
# usage message
#--------------------------------------------------------------------
help :
	-@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	-@echo "important environmental variables you can set on the cmd line:"
	-@echo "METAG_DIR=xxx"
	-@echo "  pass on make cmd line, or: setenv METAG_DIR xxx"
	-@echo "WITH_PJMALLOC=0|1"
	-@echo "PJPATH=..."
	-@echo "USE_SORTED_DB=0|1"
	-@echo "TAXID_SIZE=32|16"
	-@echo
	-@echo "current CXXXFLAGS: $(CXXFLAGS)"
	-@echo
	-@echo "current INCLUDE: $(INCLUDE)"
	-@echo "current LIBS: $(LIBS)"
	-@echo
	-@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	-@echo
	-@echo
