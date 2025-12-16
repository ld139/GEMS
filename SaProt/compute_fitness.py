from model.saprot.saprot_foldseek_mutation_model import SaprotFoldseekMutationModel
from typing import List, Optional, Dict, Tuple
from utils.foldseek_util import get_struc_seq
from Bio import PDB
import pandas as pd
import numpy as np
import argparse
import os
def spearmanr(x, y):

    x = np.array(x)
    y = np.array(y)

    if x.shape != y.shape:
        raise ValueError("Input arrays must have the same length.")
    
 
    rank_x = np.argsort(np.argsort(x))
    rank_y = np.argsort(np.argsort(y))
    

    d = rank_x - rank_y
    d_squared = d ** 2
    

    n = len(x)
    rho = 1 - (6 * np.sum(d_squared)) / (n * (n**2 - 1))
    
    return rho
    
def extract_plddt_from_pdb(pdb_file, chain_id):
    
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

  
    chain = structure[0][chain_id]

    plddt_scores = {}
    for residue in chain:
        for atom in residue:
            if atom.get_name() == "N":  
                plddt = atom.get_bfactor()  

                residue_id = residue.get_resname() + str(residue.get_id()[1])
                plddt_scores[residue_id] = plddt

    return plddt_scores


def get_str_seq_from_pdb(pdb_path: str, chains: Optional[List[str]] = None, plddt_mask: bool = False, binary_path: str = "./bin/foldseek") -> Tuple[str, str, str]:
    
    parsed_seqs = get_struc_seq(binary_path, pdb_path, chains, plddt_mask)[chains[0]]
    seq, foldseek_seq, combined_seq = parsed_seqs

    
    plddt_per_residue = {}
    for chain in chains:
        plddt_per_residue[chain] = extract_plddt_from_pdb(pdb_path, chain)

    
    combined_seq_list = list(combined_seq)  

    seq_index = 0
    for chain in chains:
        plddt_scores = plddt_per_residue[chain]
        for i, (residue, plddt) in enumerate(plddt_scores.items()):
            if plddt < 70:
                
                if combined_seq_list[seq_index].islower():
                    combined_seq_list[seq_index] = "#" 
            seq_index += 1

    
    modified_combined_seq = ''.join(combined_seq_list)

    return seq, foldseek_seq, modified_combined_seq

def get_mut_value(foldseek_path: str = "./bin/foldseek", config_path: str = "./weights/PLMs/SaProt_650M_AF2", load_pretrained: bool = True, device: str = "cuda",
                combined_seq : str = None, mut_list: List[str] = None):

    config = {
        "foldseek_path": foldseek_path,
        "config_path": config_path,
        "load_pretrained": load_pretrained
    }
    model = SaprotFoldseekMutationModel(**config)
    tokenizer = model.tokenizer
    model.to(device)
    model.eval()
    mut_values = []
    for mut in mut_list:
        mut_info = mut
        # print(mut_info)
        mut_value = model.predict_mut(combined_seq, mut_info)
        # print(mut_value)
        mut_values.append(mut_value)
    return mut_values

def get_mute_from_csv(csv_file: str, pdb_path: str, chains: Optional[List[str]] = None, plddt_mask: bool = False,
                        foldseek_path: str = "./SaProt/bin/foldseek", config_path: str = "./weights/PLMs/SaProt_650M_AF2",
                         load_pretrained: bool = True, device: str = "cuda"):
    pdb_name = pdb_path.split("/")[-1].split(".")[0]
    ckpt_name = config_path.split("/")[-1]
    
    df = pd.read_csv(csv_file)
    # print(df.head())
    mut_list = df["mutant"].tolist()
    dm_score = df["DMS_score"].tolist()
    print("Processing", pdb_name, "using", ckpt_name, ",", "mute_numbers:", len(mut_list))

    assert len(mut_list) == len(dm_score)
    seq, foldseek_seq, combined_seq = get_str_seq_from_pdb(pdb_path, chains, plddt_mask, foldseek_path)
    mut_value = get_mut_value(foldseek_path, config_path, load_pretrained, device, combined_seq, mut_list)
    spearman = spearmanr(mut_value, dm_score)
    # print(mut_value[0])
    # print(len(mut_value))
    print("spearman for ", pdb_name, "using", ckpt_name, "is", spearman)
    df[ckpt_name] = mut_value

    output_dir = os.path.dirname(csv_file)
    csv_file = os.path.basename(csv_file)
    output_dir = os.path.join(output_dir, "SaProt")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    df.to_csv(os.path.join(output_dir, csv_file), index=False)
    return mut_value
    

def main(args):
    dataset_dir = args.dataset_dir
    CKPT_PATH = args.ckpt_path


    csv_file_name = os.path.join(dataset_dir, args.id, f"{args.id}.csv")
    pdb_file = os.path.join(dataset_dir, args.id, f"{args.id}.pdb")
    get_mute_from_csv(csv_file_name, pdb_file, chains=["A"], config_path=CKPT_PATH)

    
def parse_args():
    parser = argparse.ArgumentParser(description="Compute fitness score for each mutant in the dataset")
    parser.add_argument("--id", type=str,required=True, help='id of DMS in DMS reference file')
    parser.add_argument("--dataset_dir", type=str, default="../data/mutant_example/ProteinGym_v1/DATASET_Joint", help="The directory of the dataset")
    parser.add_argument("--ckpt_path", type=str, default="./weights/PLMs/SaProt_650M_AF2", help="The path of the checkpoint")
    return parser.parse_args()
        

if __name__ == "__main__":
    args = parse_args()
    main(args)
    print("Done!")
    
    
