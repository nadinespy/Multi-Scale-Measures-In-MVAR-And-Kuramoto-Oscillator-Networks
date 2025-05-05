#!/bin/bash

# Name of this Bash script
scriptname=$0

echo -e $(tput bold)"\n*** Batch script: $scriptname\n"$(tput sgr0)

# Simulation name (Matlab script)
simname=mec_mvar_rand_cross_coup_rmi

# Base directory for code,etc
basedir="/its/home/ns508"

# Directory to run Matlab in
rundir=$basedir"/code/multivar_autoreg/rand_cross_coup_rmi"

# Define measures as an array so that different cores are used for each measure in parallel
# Comment out as appropriate ("phiid_based_measures_ccs" missing currently).
#
# If you like you can even run concurrent instances of the *same* measure, possibly with
# different parameter scripts (see below).
#
# Example:

measures=(
#    	"dd_ce_co_info"
#	"phiid_based_measures_mmi"
#	"phiid_based_measures_ccs"
#	"integrated_info_measures"
	"integrated_info_measures"
#	"multi_info"
)

# Parameter (MATLAB) scripts - must correspond to measures array above
parmscripts=(
#    	"params_dd_ce_co_info_mvar_rand_cross_coup_rmi"
#	"params_phiid_mmi_mvar_rand_cross_coup_rmi"
#	"params_phiid_ccs_mvar_rand_cross_coup_rmi"
#	"params_integrated_info_mvar_rand_cross_coup_rmi"
	"params_control_mvar_rand_cross_coup_rmi"
#	"params_multi_info_mvar_rand_cross_coup_rmi"
)

# Number of measures
num_measures=${#measures[@]}

# Total effective number of cores (cores x threads-per-core x sockets)
tot_cores=$(nproc --all)

# Cores per measure
cores_per_measure=$(( $tot_cores / $num_measures ))

echo -e "Number of measures : $num_measures"
echo -e "Cores available    : $tot_cores"
echo -e "Cores per measure  : $cores_per_measure\n"

# Function to return cores specification for 'taskset -c'
get_cores() {
	# $1 is number of cores per measure
	# $2 is measure index
	local core_start=$(( $1 * $2 ))
	local core_end=$(( $core_start + $1 - 1 ))
	echo $core_start-$core_end
}

# Loop over the measures and run MATLAB on a different set of cores for each measure
for (( i = 0; i < $num_measures; i++ )); do

	# The measure
	measure=${measures[i]}

	# The parameter (MATLAB) script
	parmscript=${parmscripts[i]}

	# MATLAB command string
	matcmds="maxNumCompThreads($cores_per_measure); basedir = '$basedir'; measure = '$measure'; $parmscript; $simname; quit"

	# Cores to run on
	cores=$(get_cores $cores_per_measure $i)

	# for control measures computed via integrated info java toolbox
	logfile=$simname\_$measure\_"control".log
	
	# Run MATLAB with appropriate command on a specific set of cores
	#
	# 'nohup' prevents the process terminating if the invoking shell exits (e.g., on logout from terminal)
	#
	# 'nice' says to reduce the process priority (be nice to other users!)
	#
	# 'tasket -c $cores" restricts the set of cores that the process may run on
	#
	# '> $logfile 2>&1' says that there is no input, andboth standard output and error output go to the logfile
	#
	# The '&' at the end puts the process in the background

	cd $rundir && nohup nice taskset -c $cores matlab -nodisplay -r "$matcmds" > $logfile < /dev/null 2>&1 &

	echo -e "Process $(($i+1)) running on cores $cores"
	echo -e "MATLAB commands : $matcmds"
	echo -e "Logfile         : $logfile\n"
done

echo -e $(tput bold)"*** Script completed: $num_measures MATLAB processes running in the background (you can log out now)\n"$(tput sgr0)
