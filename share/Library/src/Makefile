#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES

# These modules provide the "infrastructure" for SWMF and components.
# A component does not have to use the whole library.
# The minimum reqguirement is ModIoUnit.o

# This rule is used to regenerate ModMpiInterface.f90 from a modified SWMF_MPI_routines.dat
MPIINTERFACE:
	${SCRIPTDIR}/CreateModMpi.pl

MODULES = \
	${BLAS}\
	ModBlasLapack.o\
	ModConst.o\
	ModCoordTransform.o\
	ModExactRS.o\
	ModFreq.o\
	ModInitialState.o\
	ModInterpolate.o\
        ModInterGen.o\
        ModInterWrapper.o\
	ModIoUnit.o\
	ModKind.o\
	ModLinearAdvection.o\
	ModLinearSolver.o\
	ModLookupTable.o\
	ModMpiInterfaces.o\
	ModMpi.o\
	ModMpiOrig.o\
        ModNumConst.o\
	ModPlanetConst.o\
	ModPlotFile.o\
	ModProcessVarName.o\
	ModReadParam.o\
	ModSort.o\
	ModTimeConvert.o\
	ModTriangulate.o\
	ModUtilities.o\
	CON_axes.o\
	CON_geopack.o\
	CON_line_extract.o\
	CON_planet.o\
	CON_planet_field.o\
	ModHdf5Utils.o\
	CON_ray_trace.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${MODULES}

MY_LIB = ${LIBDIR}/libSHARE.a

LIB:	DEPEND
	${MAKE} ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES}


distclean: clean
	rm -f Makefile.DEPEND


