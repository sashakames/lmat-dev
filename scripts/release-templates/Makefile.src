ifndef PJPATH
PJPATH=$(PWD)
endif

ifndef HAVE_MPI
HAVE_MPI=0
endif


all:
	make -C third-party/ all METAG_DIR=$(PWD)
	make -C src/kmerdb lib METAG_DIR=$(PWD) PJPATH=$(PJPATH) HAVE_MPI=$(HAVE_MPI)
	make -C src/ all METAG_DIR=$(PWD) PJPATH=$(PJPATH)  HAVE_MPI=$(HAVE_MPI)

clean :

	make -C src/kmerdb clean METAG_DIR=$(PWD)
	make -C src/ clean METAG_DIR=$(PWD)
