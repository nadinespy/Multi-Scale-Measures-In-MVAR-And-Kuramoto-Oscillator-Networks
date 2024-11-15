#!/bin/bash

# Base name of this script (without extension)
scriptname=$(basename $0 .sh)

echo -e "\n*** Running batch script '$scriptname' ***\n"

# define measures as an array so that different cores are used for each measure in parallel ("phiid_based_measures_ccs" missing currently)
measures=(
    "DD_CE_CO_INFO"
    #"phiid_based_measures_mmi"
    #"integrated_info_measures"
    #"multi_info"
)

# define the number of cores to use (equal to number of measures)
num_cores_per_measure=24
num_cores=$((${#measures[@]} * num_cores_per_measure))
#num_cores=${#measures[@]}

# define integrated info measures to be calculated when 'measures' is equal to 'integrated_info_measures'
calcNames="IntegratedInformation,IntegratedInteraction,DecoderIntegration,CausalDensity,IntegratedSynergy,AverageCorrelation,TimeDelayedMutualInfo"

# define other variables
n_nodes=8
m_dim=1
dim_reduction="pca,grassmanian"
time_lag_for_model=1
n_samples_couplings=50
n_samples_noise_corrs=50
seed=1
n_rmi=100
n_norms=100
spectral_radius=0.9
filename_table="results_dd_ce_co_info_table"
filename="results_dd_ce_co_info"

# define script name
scriptname="bash_measures_random_couplings_neocortex.ssh"

# function to run MATLAB for a specific measure on a core
run_matlab() {
    local measure="$1"
    local core_index="$2"

    # Adjust the core indices for each measure to ensure they are assigned to different cores
    local cores_assigned=$((core_index * num_cores_per_measure))

    # create the command string for the specific measure
    cmd="global MYSCRIPT_ARGS; MYSCRIPT_ARGS={'$n_nodes','$m_dim','$dim_reduction','$time_lag_for_model','$n_samples_couplings','$n_samples_noise_corrs','$seed','$measure','$calcNames','$n_rmi','$n_norms','$spectral_radius','$filename_table','$filename'}; measures_random_couplings_neocortex;"

    # define log file name (one for each measure)
    logfile="$scriptname_$measure.log"

    # run MATLAB on a specific core with the generated command string; '&' at the end runs Matlab in the background
    nohup taskset -c $cores_assigned-$((cores_assigned + num_cores_per_measure - 1)) matlab -nodisplay -r "rehash; $cmd; exit" > $logfile 2>&1 &
    #nohup taskset -c $core_index matlab -nodisplay -r "rehash; $cmd; exit" > $logfile 2>&1 &
}

# loop over the measures and run MATLAB for each measure on a different core
for ((i = 0; i < ${#measures[@]}; i++)); do
    core_start=$((i * num_cores_per_measure))
    measure="${measures[i]}"
    run_matlab "$measure" "$core_start"
done

#for ((i = 0; i < num_cores; i++)); do
#    measure="${measures[i]}"
#    run_matlab "$measure" "$i"
#done

echo "script is running in the background, you can log out now"

