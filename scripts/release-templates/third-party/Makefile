all :
	tar xvfz perm-je-0.9.7.tgz
	cd perm-je-0.9.7; ./configure --with-jemalloc-prefix=pj ; make; cp lib/*.a ../../lib; cp -r include/jemalloc ../../include;  cd ..
	tar xvfz gzstream.tgz
	cd gzstream ; make; cp *.a ../../lib ; cp gzstream.h ../../include ; cd ..
	tar xvfz zlib-1.2.7.tar.gz
	cd zlib-1.2.7; ./configure ; make ; cp *.a ../../lib; cd ..
	tar xvfz parallel-20130822.tar.gz
	cd parallel-20130822; ./configure ; make ; cp src/parallel ../../bin
