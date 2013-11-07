CXX = g++

INCLUDE = -I$(PJPATH)/include -I. -I$(METAG_DIR)/src/kmerdb -I$(METAG_DIR)/include 



CXXFLAGS =  -g -O3 -Wall -fopenmp $(INCLUDE) -D$(PJOPT)  \
            -DIDX_CONFIG=$(SDBIDX) -DTID_SIZE=$(TAXID_SIZE) \
            -DDBTID_T=$(TID_T) -Wno-deprecated -DUSE_SORTED_DB=$(USE_SORTED_DB)

#LDFLAGS = -static


METAG_LIB = $(METAG_DIR)/src/kmerdb/lib/libmetag.a
LIBS = $(METAG_LIB)  $(PJMALLOC_LIB) -lpthread -ldl \
       $(METAG_DIR)/lib/libgzstream.a \
       $(METAG_DIR)/lib/libz.a

LDENV= LD_RUN_PATH=/lib:/usr/lib:$(HOME)/lib:/usr/local/lib:$(PJPATH)/lib
RUNENV= LD_LIBRARY_PATH=/lib:/usr/lib:$(HOME)/lib:/usr/local/lib:$(PJPATH)/lib

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $< -o $@


#--------------------------------------------------------------------
# macros for conditional compilation
#--------------------------------------------------------------------


#compile and link with PJMALLOC?
ifndef WITH_PJMALLOC
 PJOPT = WITH_PJMALLOC=1
 PJMALLOC_LIB =  $(PJPATH)/lib/libjemalloc.a
else
 PJOPT = WITH_PJMALLOC=$(WITH_PJMALLOC)
 ifeq ($(WITH_PJMALLOC), 1)
 PJMALLOC_LIB =  $(PJPATH)/lib/libjemalloc.a  
 endif 
endif

#location of PJMALLOC lib
ifndef PJPATH
PJPATH= $(METAG_DIR)/
endif

ifeq ($(USE_SORTED_DB),0)
USE_SORTED_DB=0
else
USE_SORTED_DB=1 
endif

ifndef SDBIDX
SDBIDX=2027
endif

ifndef TAXID_SIZE
TID_T = uint16_t
TAXID_SIZE=16
else
ifeq ($(TAXID_SIZE), 32)
TID_T = uint32_t
else
TID_T = uint16_t
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
	-@echo "WITH_GENOME_LOCATIONS=0|1"
	-@echo "WITH_PJMALLOC=0|1"
	-@echo "PJPATH=..."
	-@echo "PJMALLOC_LIB =...\n";
	-@echo "USE_OMP=0|1 - 1 is default"
	-@echo 
	-@echo "current CXXXFLAGS: $(CXXFLAGS)"
	-@echo 
	-@echo "current INCLUDE: $(INCLUDE)"
	-@echo 
	-@echo "current LIBS: $(LIBS)"
	-@echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	-@echo 
	-@echo 
