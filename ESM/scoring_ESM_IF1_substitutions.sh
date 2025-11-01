#!/bin/bash 

# source ../zero_shot_config.sh
# source activate proteingym_env

## Regression weights are at: https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t33_650M_UR50S-contact-regression.pt
#https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t33_650M_UR50S-contact-regression.pt


# pdbid=tyrh
export DMS_id="tyrh"
export model_checkpoint='/home/luod/.cache/torch/hub/checkpoints/esm_if1_gvp4_t16_142M_UR50.pt'
export DMS_output_score_folder="/mnt/g/software/GEMS/DATASET/${DMS_id}/ESM-IF1"



python compute_fitness_esm_if1.py \
    --model_location ${model_checkpoint} \
    --structure_folder /mnt/g/software/GEMS/DATASET/${DMS_id} \
    --DMS_id ${DMS_id} \
    --DMS_data_folder /mnt/g/software/GEMS/DATASET/${DMS_id} \
    --output_scores_folder ${DMS_output_score_folder}