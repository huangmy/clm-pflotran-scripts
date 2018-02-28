
We provide detailed notes on running the coupled CLM-PFLOTRAN model
on NERSC's Edison supercomputer.


### Install PETSc

```
module rm PrgEnv-intel
module rm PrgEnv-cray 
module rm PrgEnv-gnu
module rm cray-hdf5
module rm cray-hdf5-parallel

module load PrgEnv-intel
module swap intel intel/16.0.3.210
module load cray-hdf5-parallel/1.10.0.3
module load cmake/3.8.2

setenv BASE_DIR <directory-of-choice>
cd $BASE_DIR
setenv PETSC_HASH 03c0fad465
git clone https://bitbucket.org/petsc/petsc petsc_$PETSC_HASH
cd petsc_$PETSC_HASH
git checkout $PETSC_HASH

setenv PETSC_DIR $BASE_DIR/petsc_$PETSC_HASH
setenv PETSC_ARCH cori_intel_O

./config/configure.py         \
--PETSC_ARCH=$PETSC_ARCH      \
--with-cc=cc                  \
--with-cxx=CC                 \
--with-fc=ftn                 \
--CFLAGS='-fast -no-ipo'      \
--CXXFLAGS='-fast -no-ipo'    \
--FFLAGS='-fast -no-ipo'      \
--with-shared-libraries=0     \
--with-debugging=0            \
--with-clanguage=c            \
--with-x=0                    \
--download-parmetis=1         \
--download-metis=1            \
--with-hdf5=1                 \
--with-hdf5-dir=$HDF5_DIR     \
--with-c2html=0               \
--with-clib-autodetect=0      \
--with-fortranlib-autodetect=0 \
--with-cxxlib-autodetect=0    \
--LIBS=-lstdc++


make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all

```


### Download code repository


#### PFLOTRAN code
```
cd $BASE_DIR
git clone https://bitbucket.org/pflotran/pflotran
setenv PFLOTRAN_SRC_DIR $BASE_DIR/pflotran
cd $PFLOTRAN_SRC_DIR

```

#### CLM code
```
cd $BASE_DIR
git clone git@github.com:CLM-PFLOTRAN/clm-pflotran.git 
setenv CLM_SRC_DIR $BASE_DIR/clm-pflotran
cd $CLM_SRC_DIR
```

### Download data repository
```
cd $BASE_DIR
mkdir cases
git clone git@bitbucket.org:pnnl_sbr_sfa/notes-for-gmd-2017-35.git
setenv CASE_DIR $BASE_DIR/cases
cd $BASE_DIR/notes-for-gmd-2017-35
setenv INPUTDATA_DIR ${PWD}/datasets
```

#### Download data from NCAR repo, instruction for registration can be found at http://www.cesm.ucar.edu/models/cesm1.2

```
svn export https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/atm/cam/chem/trop_mozart/emis/megan21_emis_factors_c20120313.nc  \
${INPUTDATA_DIR}/cesm-inputdata/atm/cam/chem/trop_mozart/emis/megan21_emis_factors_c20120313.nc

svn export https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_monthly_1850_mean_1.9x2.5_c090421.nc \
${INPUTDATA_DIR}/cesm-inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_monthly_1850_mean_1.9x2.5_c090421.nc

svn export https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/lai_streams/MODISPFTLAI_0.5x0.5_c140711.nc \
${INPUTDATA_DIR}/cesm-inputdata/lnd/clm2/lai_streams/MODISPFTLAI_0.5x0.5_c140711.nc

svn export https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/paramdata/clm_params.c140423.nc \
${INPUTDATA_DIR}/cesm-inputdata/lnd/clm2/paramdata/clm_params.c140423.nc

svn export https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc \
${INPUTDATA_DIR}/cesm-inputdata/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc

svn export https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c090915.nc \
${INPUTDATA_DIR}/cesm-inputdata/lnd/clm2/snicardata/snicar_optics_5bnd_c090915.nc

```

### Create a CLM-PFLOTRAN case
```
./create_case.sh                       \
-clm_source_dir      $CLM_SRC_DIR      \
-pflotran_source_dir $PFLOTRAN_SRC_DIR \
-inputdata_dir       $INPUTDATA_DIR    \
-case_dir            $CASE_DIR         \
-petsc_dir           $PETSC_DIR        \
-petsc_arch          $PETSC_ARCH       \
-case_name           10m_case          \
-resolution          10m
```

### Reference
Bisht, G., Huang, M., Zhou, T., Chen, X., Dai, H., Hammond, G., Riley, W., Downs, J., Liu, Y., and Zachara, J.:
Coupling a three-dimensional subsurface flow and transport model with a land surface model to simulate stream-aquifer-land interactions (CP v1.0), Geosci. Model Dev., 10(12):4539-4562.  doi:10.5194/gmd-10-4539-2017.

### Declaimer
CLM4.5 is an open-source software released as part of the Community Earth System Model (CESM) version 1.2 (http://www.cesm.ucar.edu/models/cesm1.2). The version of CLM4.5 used in CP v1.0 is a branch from the CLM developer's repository. Its functionality is scientifically consistent with descriptions in Oleson et al. [2013] with source codes refactored for a modular code design. Additional minor code modifications were added by the authors to support coupling with PFLOTRAN.  Permission from the CESM Land Model Working Group has been obtained to release this CLM4.5 development branch but the National Center for Atmospheric Research cannot provide technical support for this version of the code CP v1.0. PFLOTRAN is an open-source software distributed under the terms of the GNU Lesser General Public License as published by the Free Software Foundation either version 2.1 of the License, or any later version. The CP v1.0 has two separate, open-source repositories for CLM4.5 and PFLOTRAN at:
https://bitbucket.org/CLM-PFLOTRAN/clm-pflotran
https://bitbucket.org/pflotran/pflotran
The README guide for the CP v1.0 and dataset used in this study are available from the open-source repository https://github.com/huangmy/clm-pflotran-scripts.git.
