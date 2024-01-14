######################################################################
# Make 3D electricdynamics module from NCAR
# The module was originally developed by A. Maute and A. Richmond
# QZ 10/29/2016
######################################################################
include ../../Makefile.def
include ../../Makefile.conf

MY_LIB = ${LIBDIR}/libEDYNAMICS.a

include Makefile.DEPEND

OBJECTS = \
	params.o\
	calc_efield.o\
	areas.o\
	remap_wt.o\
	grid.o\
	fieldline_p.o\
	fieldline_s.o\
	calc_B.o\
	fieldline_r.o\
	qd_grid.o\
	calc_Iqd.o\
	apex.o\
	get_apex.o\
	calc_I.o\
	calc_J.o\
	calc_Jpg.o\
	calc_S1approx.o\
	calc_MN_S1S2.o\
	calc_je_s1s2.o\
	calc_S.o\
	readin.o\
	calc_coef.o\
	ModFAC.o\
	edynamicsout.o

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

