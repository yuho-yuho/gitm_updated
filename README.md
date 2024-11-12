# gitm_on_tacc
Default version of the Global Ionosphere/Thermosphere Model used by the Upper Atmosphere Group at University of Texas at Arlington.

GITM has been developed in fortran-90. Original code and copyright by University of Michigan. Please refer to the paper by Ridley, Deng & Tóth (2006) at https://doi.org/10.1016/j.jastp.2006.01.008. 

Other useful links:

TACC LS6 User Guide: https://docs.tacc.utexas.edu/hpc/lonestar6/

## HPC Environments & Dependencies:

1. GITM runs on TACC machine. 
2. GITM needs MPI to work.

## Quick Start:

1\. Clone the repository on your TACC Home directory

```shell
git clone https://github.com/yuho-yuho/gitm_default.git
```

2\. Go into the folder (You can change to whatever you want)

```shell
cd gitm_default
```

3\. Configure the Fortran compiler with ifort (By default)

```shell
./Config.pl -install -compiler=ifortmpif90 -earth
```

4\. Add the following four lines to the end of the /GITM/Makefile.def file:

```shell
EDYNAMICSDIR = ${UADIR}/util/EDYNAMICS
EMPIRICALAEPMDIR    = ${UADIR}/util/EMPIRICAL/srcAEPM
EMPIRICALEPMDIR    = ${UADIR}/util/EMPIRICAL/srcEPM
EMPIRICALEFVMDIR    = ${UADIR}/util/EMPIRICAL/srcEFVM
```

5\. Compile your GITM codes

```shell
gmake
```

6\. Create your run directory

```shell
make rundir
```

7\. Link input data for your rundir

```shell
cd run/DataIn
ln -s ~/GITM/srcData/* .
```

8\. Apply idev with 4 nodes & 144 mpi tasks

```shell
cd run
idev -m 10 -N 4 -n 144
```

9\. Run your GITM on idev (I use ibrun instead of mpirun)

```shell
ibrun ./GITM.exe
```
