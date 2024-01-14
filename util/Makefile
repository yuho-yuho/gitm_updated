INSTALLFILE = \
	DATAREAD/srcIndices/Makefile.DEPEND \
	EMPIRICAL/srcEE/Makefile.DEPEND \
	EMPIRICAL/srcIE/Makefile.DEPEND \
	EMPIRICAL/srcUA/Makefile.RULES \
	NOMPI/src/Makefile.RULES \
	CRASH/src/Makefile.DEPEND \
	CRASH/src/Makefile.RULES

install:
	touch ${INSTALLFILE}
	@(if [ -d HYPRE ]; then cd HYPRE;  make install; fi);

clean:
	touch ${INSTALLFILE}
	cd NOMPI/src;                 make clean
	cd TIMING/src;                make clean
	cd TIMING/srcEmpty;           make clean
	cd TIMING/doc;                make clean
	cd DATAREAD/srcIndices;       make clean
	cd DATAREAD/srcMagnetogram;   make clean
	cd DATAREAD/srcDemt;          make clean
	cd EMPIRICAL/srcEE;           make clean
	cd EMPIRICAL/srcIE;           make clean
	cd EMPIRICAL/srcGM;           make clean
	cd EMPIRICAL/srcUA;           make clean
	cd CRASH/src;                 make clean
	cd CRASH/doc/Tex;	      make clean

distclean:
	touch ${INSTALLFILE}
	cd NOMPI/src;                 make distclean
	cd TIMING/src;                make distclean
	cd TIMING/srcEmpty;           make distclean
	cd TIMING/doc;                make distclean
	cd DATAREAD/srcIndices;       make distclean
	cd DATAREAD/srcMagnetogram;   make distclean
	cd DATAREAD/srcDemt;          make distclean
	cd EMPIRICAL/srcEE;           make distclean
	cd EMPIRICAL/srcIE;           make distclean
	cd EMPIRICAL/srcGM;           make distclean
	cd EMPIRICAL/srcUA;           make distclean
	cd CRASH/src;                 make distclean
	cd CRASH/doc/Tex;             make distclean
	@(if [ -d HYPRE ]; then cd HYPRE;  make distclean; fi);
	rm -f *~
	rm -f ${INSTALLFILE}

dist: distclean
	/bin/tar -cf util.tar .
