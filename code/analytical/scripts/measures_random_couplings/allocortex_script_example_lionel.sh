#!/bin/bash

# base name of this script (without extension)
scriptname=$(basename $0 .sh)

echo -e "\n*** Running batch script '"$scriptname"' ***\n"

# path to Matlab code directory
codedir=$LOCALREPO/matlab/ncomp/PsyMEG

# current directory (the directory this script is run in - log files will go here)
currdir=$(pwd -P)

# data root directory (environmental variable: where to find PsyMEG directory)
# DATADIR=/shared

datadir=SDT
ds=3
regmode=OLS
maxmo=20

# run multiple concurrent Matlab sessions
for drugnum in 1 2 3 4; do

	# cores (note: zero offset!) - adjust per machine
	case $drugnum in
		1 ) cores="0-5"   ;;
		2 ) cores="6-11"  ;;
		3 ) cores="12-17" ;;
		4 ) cores="18-23" ;;
	esac

	# Matlab invocation
	runmatlab="nohup nice taskset -c $cores matlab -nojvm -nodisplay"

    # the log file
    logfile=$currdir/$scriptname\_drug\_$drugnum.log

    # Matlab commands
    matcmds="\
		datadir = '$datadir';\
		drugnum =  $drugnum;\
		ds      =  $ds;\
		regmode = '$regmode';\
		maxmo   =  $maxmo;\
		plotm   =  [];\
		verb    =  false;\
		VARMOs;\
		quit"

    # run Matlab
    cd $codedir && $runmatlab -r "$matcmds" > $logfile < /dev/null 2>&1 &

done
