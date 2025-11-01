#!/bin/bash 

# source ../zero_shot_config.sh
# source activate proteingym_env
# conda activate
# MSA transformer checkpoint 
# pdbid=tyrh
DMS_id="tyrh"
# export DMS_reference_file_path_subs="./DMS_substitutions.csv"
DMS_data_folder_subs="/mnt/g/software/GEMS/DATASET/${DMS_id}"
dms_output_folder="/mnt/g/software/GEMS/DATASET/${DMS_id}/MSA_Transformer/"
model_checkpoint="/home/luod/.cache/torch/hub/checkpoints/esm_msa1b_t12_100M_UR50S.pt"

DMS_MSA_weights_for_MSA_Transformer_folder="/mnt/g/software/GEMS/MSA_weights"
DMS_MSA_data_folder="/mnt/g/software/GEMS/MSA"

# export dms_output_folder="${DMS_output_score_folder_subs}/MSA_Transformer/"
scoring_strategy=masked-marginals # MSA transformer only supports "masked-marginals"
model_type=MSA_transformer
scoring_window="optimal"
random_seeds="1 2 3 4 5"
# export DMS_MSA_weights_for_MSA_Transformer_folder="${DMS_MSA_weights_folder}/DMS_msa_weights_for_MSA_Transformer" # Use weights recomputed post MSA filtering used in MSA Transformer

python compute_fitness.py \
    --model-location ${model_checkpoint} \
    --model_type ${model_type} \
    --dms_id ${DMS_id} \
    --dms-input ${DMS_data_folder_subs} \
    --dms-output ${dms_output_folder} \
    --scoring-strategy ${scoring_strategy} \
    --scoring-window ${scoring_window} \
    --msa-path ${DMS_MSA_data_folder} \
    --msa-weights-folder ${DMS_MSA_weights_for_MSA_Transformer_folder} \
    --seeds ${random_seeds}
