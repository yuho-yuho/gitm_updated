# Makefile for the EFVM module in GITM
# Created: Qingyu Zhu, 08/20/2020
# ----------------------------------------------------------------------------

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libEFVM.a

include Makefile.DEPEND

OBJECTS = \
	EFVM_Interface.o \
	EFVM_Initialization.o \
	EFVM_util.o\
	EFVM_Ualibrary.o \
	EFVM.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl -I${SHAREDIR} ${OBJECTS}

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

distclean: clean
	rm -f Makefile.DEPEND
