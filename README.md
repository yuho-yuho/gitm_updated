# gitm_updated
Updated version of the Global Ionosphere/Thermosphere Model used by the Upper Atmosphere Group at University of Texas at Arlington.

For basic info about running GITM on TACC, please refer to the repository at /gitm_default/. 

## Main Updates:

1. Updated low-latitude dynamo: NCAR Electrodynamo added 
2. Updated high-latitude forcing specifications, e.g., AMIE
3. Updated interhemispheric Asymmetry, e.g., IMF By effect
4. Coupled with new empirical models for high-latitude forcing
5. Modefied auroral precipitation profile and energy types

## Quick Start:

1\. Clone the repository on your TACC Home directory

```shell
git clone https://github.com/yuho-yuho/gitm_updated.git
```

2\. Go into the folder (You can change to whatever you want)

```shell
cd gitm_updated
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

Go to this Google Drive folder: https://drive.google.com/drive/folders/1xLMNRsTMPyqVM-TcP0lm4lwFOURQARow?usp=sharing
Directly copy/replace Makefile, netcdf.mod, ModGITM.f90 and output_common.f90 to your /gitm_updated/src/ folder

5\. Compile your GITM codes

```shell
gmake
cd src
```

Open the ModSize.f90 file and modify the 'nLons' and 'nlats' from x to 12:

```shell
integer, parameter :: nLons = 12
integer, parameter :: nLats = 12
```

Then return to the previous directory:

```shell
cd ..
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

Open the UAM.in file and modify the 'lons' and 'lats' from 2 to 12:

```shell
#GRID
12              lons
12              lats
```

9\. Run your GITM on idev (I use ibrun instead of mpirun)

```shell
ibrun ./GITM.exe
```

Alternatively, you can run GITM with batch. First, make a batch script with the following steps:

```shell
emacs standard_run
```

Copy & Paste the following content to your "standard_run" script:

```shell
#!/bin/bash
#SBATCH -J hongyu            # job name
#SBATCH -o gitm.o%j        # output and error file name (%j expands to jobID)
#SBATCH -N 2                # number of nodes requested
#SBATCH -n 144               # total number of mpi tasks requested
#SBATCH -p normal      # queue (partition) -- normal, development, etc.
#SBATCH -t 00:10:00         # run time (hh:mm:ss) - 1.5 hours

module load netcdf/4.6.2
# run the executable named a.out                                                                    
ibrun ./GITM.exe

./pGITM
```

Then, submit your GITM simulation as a batch job by using:

```shell
sbatch standard_run
```
