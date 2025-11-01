from ast import main
import subprocess
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='runing SaProt, ESM-IF1 and MSA Transformer')
    parser.add_argument("--data_dir", default= "DATASET", help='input directory')
    parser.add_argument('--DMS_id',type=str,required=True, help='id of DMS in DMS reference file')
    # parser.add_argument('-wt', '--wildtype', required=False, help='wildtype fasta file')
    parser.add_argument('--offset',type=int,help='offset of mutation sites',default=1)
    parser.add_argument('--sites', required=False, 
                        help='mutation sites file OR comma-separated list of sites (e.g., "10,20,30")')
    parser.add_argument('--all_sites', action='store_true', 
                        help='perform saturation mutagenesis on all sites in the sequence')
    parser.add_argument('--combinatorial', type=str, required=False,
                        help='comma-separated list of sites for combinatorial saturation mutagenesis (includes single and multiple mutations)')
    # parser.add_argument('--output', required=False, help='output directory')

    parser.add_argument("--ckpt_path_saprot", type=str, default="./SaProt/weights/PLMs/SaProt_650M_AF2", help="The path of the SaProt checkpoint")
    parser.add_argument("--ckpt_path_esmif1", type=str, default="/home/luod/.cache/torch/hub/checkpoints/esm_if1_gvp4_t16_142M_UR50.pt", help="The path of the ESM-IF1 checkpoint")
    parser.add_argument("--ckpt_path_msatransformer", type=str, default="/home/luod/.cache/torch/hub/checkpoints/esm_msa1b_t12_100M_UR50S.pt", help="The path of the MSA Transformer checkpoint")
    parser.add_argument('--esmfold_env', type=str, default='/home/luod/miniconda3/envs/esmfold/bin/python', help='path to the python executable in the esmfold conda environment')
    parser.add_argument('--saprot_env', type=str, default='/home/luod/miniconda3/envs/SaProt/bin/python', help='path to the python executable in the saprot conda environment')

    return parser.parse_args()


def main():
    args = parse_args()
    DMS_id = args.DMS_id
    data_dir = args.data_dir
    wt_fasta = os.path.join(data_dir, DMS_id, f'{DMS_id}.fasta')
    # pdb = os.path.join(data_dir, DMS_id, f'{DMS_id}.pdb')
  
    offset = args.offset
    sites = args.sites
    all_sites = args.all_sites
    combinatorial = args.combinatorial
    output = "." + os.sep + data_dir
    saprot_env = args.saprot_env
    esmfold_env = args.esmfold_env
    msa_path ="." + os.sep + "MSA" 
    MSA_weights = "." +  os.sep + "MSA_weights" 
    # run mutation generation
    cmd = f"{saprot_env} saturation.py -wt {wt_fasta} -offset {offset} -o {data_dir}/{DMS_id}"
    if all_sites:
        cmd += " --all_sites "
    elif sites:
        cmd += f" --sites {sites} "
    elif combinatorial:
        cmd += f" --combinatorial {combinatorial} "
    print("Generating mutation combinations...")
    subprocess.run(cmd, shell=True, check=True)

    # cal msa weights for msa trans

    cmd = f"{saprot_env} ./EVE/calc_weights.py --MSA_data_folder .{msa_path} --DMS_id {DMS_id} --MSA_weights_location .{MSA_weights} --num_cpus -1 --calc_method evcouplings --threshold_focus_cols_frac_gaps 1 --skip_existing"
    subprocess.run(cmd, shell=True, check=True, cwd="./ESM")

    # run SaProt
    # 如果
    if not os.path.exists(os.path.join(output, DMS_id, "SaProt")):
        cmd = f"{saprot_env} ./SaProt/compute_fitness.py --id {DMS_id} --dataset_dir {output} --ckpt_path {args.ckpt_path_saprot}"
        print("Running SaProt...")
        subprocess.run(cmd, shell=True, check=True)

    # run ESM-IF1
    if not os.path.exists(os.path.join(output, DMS_id, "ESM-IF1")):
        cmd = f"{esmfold_env} ./ESM/compute_fitness_esm_if1.py --DMS_data_folder {data_dir}/{DMS_id} --structure_folder {data_dir}/{DMS_id} --DMS_id {DMS_id} --model_location {args.ckpt_path_esmif1} --output_scores_folder {output}/{DMS_id}/ESM-IF1 "
        print("Running ESM-IF1...")
        subprocess.run(cmd, shell=True, check=True)

    # run MSA Transformer
    if not os.path.exists(os.path.join(output, DMS_id, "MSA_Transformer")):
        cmd = f"{esmfold_env} ./ESM/compute_fitness.py --model-location {args.ckpt_path_msatransformer} --model_type MSA_transformer --dms_id {DMS_id} --dms-input {data_dir}/{DMS_id} --dms-output {output}/{DMS_id}/MSA_Transformer --scoring-strategy masked-marginals --scoring-window optimal --msa-path {msa_path} --msa-weights-folder {MSA_weights} --seeds 1 2 3 4 5"
        print("Running MSA Transformer...")
        subprocess.run(cmd, shell=True, check=True)

    # run GEMME
    if not os.path.exists(os.path.join(output, DMS_id, "GEMME")):
        cmd = f"sh run_gemme.sh {DMS_id}"
        print("Running GEMME...")
        subprocess.run(cmd, shell=True, check=True)

    # run ranking
    cmd = f"python ranking.py --data_path {data_dir} --id {DMS_id} --offset {offset-1}" # offset -1
    print("Ranking ensemble models...")
    subprocess.run(cmd, shell=True, check=True)
    
if __name__ == "__main__":
    main()
