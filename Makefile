include Makefile.inc

.PHONY: all
all:
	make -C src/kmerdb clean
	make -C src/kmerdb lib
	make -C apps all
	make -C apps/make_db all

.PHONY: clean
clean:
	make -C tests clean
	make -C apps clean
	make -C apps/make_db clean
	make -C src/kmerdb clean
	make -C src/kmerdb/examples/tests/data clean
	make -C src/kmerdb/examples/tests clean
	make -C src/kmerdb/examples/tutorials clean
	make -C src/kmerdb/drivers clean

.PHONY: atest
atest: # all
	make -C apps test
	make -C apps/make_db test

.PHONY: test
test:
	./run_tests.py
