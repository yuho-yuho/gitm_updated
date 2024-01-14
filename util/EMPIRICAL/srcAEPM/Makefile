# Makefile for the SAM module in GITM
# Created: Qingyu Zhu, 10/16/2019
# ----------------------------------------------------------------------------

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libAEPM.a

include Makefile.DEPEND

OBJECTS = \
	AEPM_Interface.o \
	AEPM_Initialization.o \
	AEPM_util_base.o\
	AEPM_util.o\
	AEPM_Ualibrary.o \
	aepm.o

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
