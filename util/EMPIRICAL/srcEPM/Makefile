# Makefile for the SAM module in GITM
# Created: Qingyu Zhu, 10/16/2019
# ----------------------------------------------------------------------------

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libEPM.a

include Makefile.DEPEND

OBJECTS = \
	ModSPHARM.o \
	EPM_Interface.o \
	EPM_Initialization.o \
	EPM_util.o\
	EPM_Ualibrary.o \
	epm.o

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
