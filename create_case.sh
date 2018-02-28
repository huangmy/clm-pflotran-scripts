#!/bin/sh

BASE_DIR=
CLM_SRC_DIR=
PFLOTRAN_SRC_DIR=
INPUTDATA_DIR=
CASE_DIR=

CESM_COMPSET=I1850CLM45

usage ()
{
 echo "usage: $0 [options] "
 echo     "OPTIONS: "
 echo     "        -clm_source_dir                  <dirpath>     Full path to the directory containing CLM source code (required)"
 echo     "        -pflotran_source_dir             <dirpath>     Full path to the directory containing PFLOTRAN source code (required)"
 echo     "        -inputdata_dir                   <dirpath>     Full path to the directory containing inputdata (required)"
 echo     "        -case_dir                        <dirpath>     Full path to the directory to put the new case (optional)"
 echo     "        -case_name                       <name>        Specifies the case name (optional)"
 echo     "        -petsc_dir                       <dirpath>     Path to PETSc directory (optional)"
 echo     "        -petsc_arch                      <name>        The PETSC_ARCH variable used during configuration and build of PETSc(optional)"
 echo     "        -resolution                      [2m|10m|20m]  The PETSC_ARCH variable used during configuration and build of PETSc(optional)"
}


# Get command line arguments
while [ $# -gt 0 ]
do
  case "$1" in
    -clm_source_dir                       ) CLM_SRC_DIR="$2"; shift;;
    -pflotran_source_dir                  ) PFLOTRAN_SRC_DIR="$2"; shift;;
    -inputdata_dir                        ) INPUTDATA_DIR="$2"; shift;;
    -case_dir                             ) CASE_DIR="$2"; shift;;
    -case_name                            ) CASE_NAME="$2"; shift;;
    -petsc_dir                            ) PETSC_DIR_LOC="$2"; shift;;
    -petsc_arch                           ) PETSC_ARCH_LOC="$2"; shift;;
    -resolution                           ) RES_NAME="$2"; shift;;
    -*) usage
      exit 1;;
    *)  break;;	# terminate while loop
  esac
  shift
done

#
# Check case-related options
#
case "$CLM_SRC_DIR" in
  "")
    echo "Invalid value for --clm_source_dir."
    usage
    exit 1
    ;;
   *)
    if [ ! -d "$CLM_SRC_DIR" ]; then
      echo "Bailing out as the following clm source directory does not exist: " ${CLM_SRC_DIR}
      usage
      exit 1
    fi
esac

case "$PFLOTRAN_SRC_DIR" in
  "")
    echo "Invalid value for --pflotran_source_dir."
    usage
    exit 1
    ;;
   *)
    if [ ! -d "$CLM_SRC_DIR" ]; then
      echo "Bailing out as the following PFLOTRAN source directory does not exist: " ${CLM_SRC_DIR}
      usage
      exit 1
    fi
esac

case "$INPUTDATA_DIR" in
  "")
    echo "Invalid value for --inputdata_dir."
    usage
    exit 1
    ;;
   *)
    if [ ! -d "$INPUTDATA_DIR" ]; then
      echo "Bailing out as the following inputdata directory does not exist: " ${INPUTDATA_DIR}
      usage
      exit 1
    fi
esac


case "$CASE_DIR" in
  ".")
    ;;
   *)
    if [ ! -d "$CASE_DIR" ]; then
      echo "Bailing out as the following case directory does not exist: " ${CASE_DIR}
      usage
      exit 1
    fi
esac

case "$RES_NAME" in
  "2m")
    PFLOTRAN_INPUTDATA_DIR=${INPUTDATA_DIR}/pflotran_2m
    CLM_USRDAT_NAME=2mx2m_300A
    DOMAINFILE_CYYYYMMDD=c20150407
    SURFFILE_CYYYYMMDD=c20150915
    PERM_FILE=PriorFields_Plume_CLMPF_2mx2mxhalfm_par1to2.h5
    RIVER_FILE=clmpf_400x400x31_2mx2mxhalfRes_material_mapped2.h5
    CLM2PF_SUBSURFACE=clm2pf_map_200x200x62PFLOTRAN_200x200x10CLM_c150226.meshmap
    CLM2PF_SURFACE=clm2pf_map_surface_200x200x62PFLOTRAN_200x200x10CLM_c150226.meshmap
    PF2CLM_SUBSURFACE=pf2clm_map_200x200x62PFLOTRAN_200x200x10CLM_c150226.meshmap
    ;;
  "10m")
    PFLOTRAN_INPUTDATA_DIR=${INPUTDATA_DIR}/pflotran_10m
    CLM_USRDAT_NAME=10mx10m_300A
    DOMAINFILE_CYYYYMMDD=c20150407
    SURFFILE_CYYYYMMDD=c20150915
    PERM_FILE=PriorFields_Plume_CLMPF_10mx10mxhalfm_par1to2.h5
    RIVER_FILE=clmpf_400x400x31_10mx10mxhalfRes_material_mapped2.h5
    CLM2PF_SUBSURFACE=clm2pf_map_40x40x62PFLOTRAN_40x40x10CLM_c150226.meshmap
    CLM2PF_SURFACE=clm2pf_map_surface_40x40x62PFLOTRAN_40x40x10CLM_c150226.meshmap
    PF2CLM_SUBSURFACE=pf2clm_map_40x40x62PFLOTRAN_40x40x10CLM_c150226.meshmap
    ;;
  "20m")
    PFLOTRAN_INPUTDATA_DIR=${INPUTDATA_DIR}/pflotran_20m
    CLM_USRDAT_NAME=20mx20m_300A
    DOMAINFILE_CYYYYMMDD=c20150407
    SURFFILE_CYYYYMMDD=c20150915
    PERM_FILE=PriorFields_Plume_CLMPF_20mx20mxhalfm_par1to2.h5
    RIVER_FILE=clmpf_400x400x31_20mx20mxhalfRes_material_mapped2.h5
    CLM2PF_SUBSURFACE=clm2pf_map_20x20x62PFLOTRAN_20x20x10CLM_c150220.meshmap
    CLM2PF_SURFACE=clm2pf_map_surface_20x20x62PFLOTRAN_20x20x10CLM_c150220.meshmap
    PF2CLM_SUBSURFACE=pf2clm_map_20x20x62PFLOTRAN_20x20x10CLM_c150220.meshmap
    ;;
  *)
    echo "Bailing out as the -resolution specified is invalid."
    usage
    exit 1
    ;;
esac

case "$CASE_NAME" in
  "")
    CASE_NAME=${CLM_USRDAT_NAME}_clmpfv-${CESM_COMPSET}-`date "+%Y-%m-%d"`
    ;;
  *)
    ;;
esac

case "$PETSC_DIR_LOC" in
  "")
    PETSC_DIR_LOC=`printenv PETSC_DIR`
    if [ "$PETSC_DIR_LOC" == "" ]; then
      echo "Could not find PETSc installation. "
      echo "PETSc installation directory is not specified via --petsc_dir and PETSC_DIR environmental variable is undefined.\n"
      usage
      exit 1
    fi
    ;;
  *)
esac

if [ ! -d "$PETSC_DIR_LOC" ]; then
  echo "Bailing out as the following PETSc directory does not exist: " ${PETSC_DIR_LOC} "\n"
  usage
  exit 1
fi

case "$PETSC_ARCH_LOC" in
  "")
    PETSC_ARCH_LOC=`printenv PETSC_ARCH`
    if [ "$PETSC_ARCH" == "" ]; then
      echo "Could not find PETSC_ARCH installation. "
      echo "PETSC_ARCH is not specified via --petsc_arch and PETSC_ARCH environmental variable is undefined.\n"
      usage
      exit 1
    fi
    ;;
  *)
esac

if [ ! -d "$PETSC_DIR_LOC/$PETSC_ARCH_LOC" ]; then
  echo "Bailing out as the following PETSc installation directory does not exist: " ${PETSC_DIR_LOC}/${PETSC_ARCH_LOC} "\n"
  usage
  exit 1
fi

export CESM_INPUTDATA_DIR=${INPUTDATA_DIR}/cesm-inputdata


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create soft links for CESM inputdata
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

mkdir -p ${CESM_INPUTDATA_DIR}/atm/datm7/${CLM_USRDAT_NAME}/CLM1PT_data
rm -rf ${INPUTDATA_DIR}/cesm_inputdata/atm/datm7/${CLM_USRDAT_NAME}/CLM1PT_data/*.nc

ls -l ${INPUTDATA_DIR}/CLM1PT_data/1x1_300A/*.nc | awk '{ print $9}' | awk -F'.' '{print $3}' | \
awk -v INPUTDATA_DIR=${INPUTDATA_DIR} -v CLM_USRDAT_NAME=${CLM_USRDAT_NAME} \
'{ system( "ln -s " INPUTDATA_DIR "/CLM1PT_data/1x1_300A/clmforc.300A." $1 ".nc " INPUTDATA_DIR"/cesm-inputdata/atm/datm7/" CLM_USRDAT_NAME "/CLM1PT_data/"$1".nc") }'

mkdir -p ${CESM_INPUTDATA_DIR}/share/domains/domain.clm
rm -rf ${CESM_INPUTDATA_DIR}/share/domains/domain.clm/domain.lnd.${CLM_USRDAT_NAME}_navy.nc
rm -rf ${CESM_INPUTDATA_DIR}/share/domains/domain.clm/domain.lnd.1x1_300A_navy.nc

ln -s ${INPUTDATA_DIR}/${CLM_USRDAT_NAME}/domain.lnd.${CLM_USRDAT_NAME}_ugrid_${DOMAINFILE_CYYYYMMDD}.nc ${CESM_INPUTDATA_DIR}/share/domains/domain.clm/domain.lnd.${CLM_USRDAT_NAME}_navy.nc
ln -s ${INPUTDATA_DIR}/CLM1PT_data/1x1_300A/domain.lnd.1x1_300A_ugrid_c20150107.nc ${CESM_INPUTDATA_DIR}/share/domains/domain.clm/domain.lnd.1x1_300A_navy.nc

mkdir -p ${CESM_INPUTDATA_DIR}/lnd/clm2/surfdata_map/
rm -rf ${CESM_INPUTDATA_DIR}/lnd/clm2/surfdata_map/surfdata_${CLM_USRDAT_NAME}_simyr1850.nc
ln -s ${INPUTDATA_DIR}/${CLM_USRDAT_NAME}/surfdata.${CLM_USRDAT_NAME}_ugrid_${SURFFILE_CYYYYMMDD}.nc ${CESM_INPUTDATA_DIR}/lnd/clm2/surfdata_map/surfdata_${CLM_USRDAT_NAME}_simyr1850.nc

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Do the PFLOTRAN stuff
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cd ${PFLOTRAN_SRC_DIR}/src/clm-pflotran
make clean
./remove_linked_files.sh
./link_files.sh
make libpflotran.a

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Now do the CLM stuff
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cd ${CLM_SRC_DIR}/scripts

# Creating case with command :
./create_newcase -case ${CASE_DIR}/${CASE_NAME} -res CLM_USRDAT -compset ${CESM_COMPSET} -mach corip1 -compiler intel

# Configuring case :
cd ${CASE_DIR}/${CASE_NAME}

# Create Macros
cp -f ${INPUTDATA_DIR}/scripts/shell/Macros.corip1         ${CASE_DIR}/${CASE_NAME}/Macros
perl -w -i -p -e "s@PETSC-DIR@${PETSC_DIR}@"               ${CASE_DIR}/${CASE_NAME}/Macros
perl -w -i -p -e "s@PETSC-ARCH@${PETSC_ARCH}@"             ${CASE_DIR}/${CASE_NAME}/Macros
perl -w -i -p -e "s@PFLOTRAN-SRC-DIR@${PFLOTRAN_SRC_DIR}@" ${CASE_DIR}/${CASE_NAME}/Macros

# Modifying : env_mach_pes.xml
if [ $RES_NAME == "2m" ]; then
./xmlchange  -file env_mach_pes.xml -id NTASKS_LND -val 1440
./xmlchange  -file env_mach_pes.xml -id NTASKS_CPL -val 240
else
./xmlchange  -file env_mach_pes.xml -id NTASKS_LND -val 40
./xmlchange  -file env_mach_pes.xml -id NTASKS_CPL -val 1
fi
 
# Modifying : env_build.xml

# Modifying : env_run.xml
./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_END   -val 2015
./xmlchange -file env_run.xml -id DATM_MODE             -val CLM1PT
./xmlchange -file env_run.xml -id STOP_N                -val 2557
./xmlchange -file env_run.xml -id REST_N                -val 73
./xmlchange -file env_run.xml -id RUN_STARTDATE         -val 0001-01-01
./xmlchange -file env_run.xml -id STOP_OPTION           -val ndays
./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_START -val 2009
./xmlchange -file env_run.xml -id DATM_CLMNCEP_YR_ALIGN -val 1
./xmlchange -file env_run.xml -id DIN_LOC_ROOT          -val ${CESM_INPUTDATA_DIR}
./xmlchange -file env_run.xml -id DIN_LOC_ROOT_CLMFORC  -val "\$DIN_LOC_ROOT/atm/datm7"
./xmlchange -file env_run.xml -id RUNDIR                -val ${SCRATCH}/${CASE_NAME}/run
./xmlchange -file env_run.xml -id CLM_USRDAT_NAME       -val ${CLM_USRDAT_NAME}
./xmlchange -file env_run.xml -id RUN_STARTDATE         -val "0001-01-01"

./cesm_setup

# Modify run script

# Modify user_nl_clm
cat >> user_nl_clm << EOF
&clm_inparm
  hist_mfilt = 1
  use_pflotran = .true.
  use_clm_soils = .false.
  hist_nhtfrq = -24
/
&clm_pflotran_inparm
  pflotran_prefix = '${CASE_NAME}'
/
EOF

# Modify user_nl_datm
cat >> user_nl_datm << EOF
&shr_strdata_nml taxmode = 'cycle, extend'
EOF

# Modify datm streams
cp ${CASE_DIR}/${CASE_NAME}/CaseDocs/datm.streams.txt.CLM1PT.CLM_USRDAT ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLM1PT.CLM_USRDAT
chmod +rw ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLM1PT.CLM_USRDAT
perl -w -i -p -e "s@domain.lnd.${CLM_USRDAT_NAME}_navy.nc@domain.lnd.1x1_300A_navy.nc@" ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLM1PT.CLM_USRDAT
sed -i '/ZBOT/d' ${CASE_DIR}/${CASE_NAME}/user_datm.streams.txt.CLM1PT.CLM_USRDAT


# Build the case
./${CASE_NAME}.build

# Copy PFLOTRAN related files in the run directory
ln -s ${PFLOTRAN_INPUTDATA_DIR}/BC_UK1_Plume2009-2015.h5                 ${SCRATCH}/${CASE_NAME}/run
ln -s ${PFLOTRAN_INPUTDATA_DIR}/BC_UK1_Plume2009-2015_Initial_Head.h5    ${SCRATCH}/${CASE_NAME}/run
ln -s ${PFLOTRAN_INPUTDATA_DIR}/DatumH_River_2009-2015.txt               ${SCRATCH}/${CASE_NAME}/run
ln -s ${PFLOTRAN_INPUTDATA_DIR}/Gradients_River_2009-2015.txt            ${SCRATCH}/${CASE_NAME}/run
ln -s ${PFLOTRAN_INPUTDATA_DIR}/$PERM_FILE                               ${SCRATCH}/${CASE_NAME}/run
ln -s ${PFLOTRAN_INPUTDATA_DIR}/$RIVER_FILE                              ${SCRATCH}/${CASE_NAME}/run
ln -s ${INPUTDATA_DIR}/scripts/matlab/mapping_files/${CLM2PF_SUBSURFACE} ${SCRATCH}/${CASE_NAME}/run
ln -s ${INPUTDATA_DIR}/scripts/matlab/mapping_files/${CLM2PF_SURFACE}    ${SCRATCH}/${CASE_NAME}/run
ln -s ${INPUTDATA_DIR}/scripts/matlab/mapping_files/${PF2CLM_SUBSURFACE} ${SCRATCH}/${CASE_NAME}/run

#obs inland condition
cp ${PFLOTRAN_INPUTDATA_DIR}/pflotran_inputdeck_for_clmpflotran_run_obs2009-2015.in ${SCRATCH}/${CASE_NAME}/run/${CASE_NAME}.in

# Running case :
./${CASE_NAME}.submit


