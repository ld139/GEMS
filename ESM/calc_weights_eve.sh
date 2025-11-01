#!/bin/bash 

# source ../zero_shot_config.sh
# conda activate  SaProt
# EVE example
export DMS_id="tyrh"
export DMS_MSA_data_folder="/mnt/g/software/GEMS/MSA"
# export DMS_reference_file_path_subs="./DMS_substitutions.csv"
export DMS_MSA_weights_folder="/mnt/g/software/GEMS/MSA_weights"
python EVE/calc_weights.py \
    --MSA_data_folder ${DMS_MSA_data_folder} \
    --DMS_id "${DMS_id}" \
    --MSA_weights_location ${DMS_MSA_weights_folder} \
    --num_cpus -1 \
    --calc_method evcouplings \
    --threshold_focus_cols_frac_gaps 1 \
    --skip_existing
    #--overwrite

