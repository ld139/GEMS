import pandas as pd
import numpy as np
import os
import re


def combine_normalized_lists(nested_lists, weights=None):

    arr = np.array(nested_lists)
    means = arr.mean(axis=1, keepdims=True)
    stds = arr.std(axis=1, keepdims=True)
    stds[stds == 0] = 1e-9  # 
    normalized = (arr - means) / stds
    
    weights = np.array(weights if weights else [1/len(arr)]*len(arr))
    weights /= weights.sum()  # 
    
    return np.dot(weights, normalized).tolist()

def read_fitness_from_txt(file_path):

    df = pd.read_csv(file_path, sep='\s+', header=None, usecols=[1], skiprows=1)
    # print(df.head())
    return df[1].astype(float).tolist()

def get_fitness_csv(csv_file, model_name):

    col_map = {
        "ESM-IF1": "esmif1_ll",
        "MSA_Transformer": "esm_msa1b_t12_100M_UR50S_ensemble",
        "SaProt": "SaProt_650M_AF2"
    }
    return pd.read_csv(csv_file)[col_map[model_name]].astype(float).tolist()


MUT_PATTERN = re.compile(r"([A-Z])(\d+)([A-Z])")

def process_mutations(mutant_str, offset):

    unified_str = mutant_str.replace(',', ':')
    mutations = unified_str.split(':')
    
    processed = []
    for mut in mutations:
        match = MUT_PATTERN.match(mut)
        if match:
            wt = match.group(1)
            site = str(int(match.group(2)) + offset)
            mut_aa = match.group(3)
            processed.append(f"{wt}{site}{mut_aa}")
    return ":".join(processed)

def extract_sites(mutant_str):

    return [MUT_PATTERN.match(mut).group(2) for mut in mutant_str.split(':')]

def rank_ensemble_models(data_path, id, offset, weights=None):
  
    esm_if1 = get_fitness_csv(os.path.join(data_path, id, "ESM-IF1", f"{id}.csv"), "ESM-IF1")
    msa_trans = get_fitness_csv(os.path.join(data_path, id, "MSA_Transformer", f"{id}.csv"), "MSA_Transformer")
    sa_prot = get_fitness_csv(os.path.join(data_path, id, "SaProt", f"{id}.csv"), "SaProt")
    
    gemme_dir = os.path.join(data_path, id, "GEMME")
    gemme_file = os.listdir(gemme_dir)[0]
    gemme = read_fitness_from_txt(os.path.join(gemme_dir, gemme_file))

  
    combined = combine_normalized_lists([esm_if1, msa_trans, sa_prot, gemme], weights or [0.25]*4)


    raw_mutants = pd.read_csv(os.path.join(data_path, id, f"{id}.csv"))["mutant"]
    processed_mutants = raw_mutants.apply(lambda x: process_mutations(x, offset))
    

    rank_df = pd.DataFrame({
        'mutant': processed_mutants,
        'fitness': combined,
        'sites': processed_mutants.apply(extract_sites)
    }).sort_values('fitness', ascending=False)
    

    rank_df[['mutant', 'fitness']].to_csv(os.path.join(data_path, id, f"{id}_rank.csv"), index=False)


    exploded_df = rank_df.explode('sites')
    for site, group in exploded_df.groupby('sites'):
        group[['mutant', 'fitness']].drop_duplicates() \
                                   .sort_values('fitness', ascending=False) \
                                   .to_csv(os.path.join(data_path, id, f"{id}_{site}_rank.csv"), index=False)
    
    return rank_df


def main(args):
    data_path = args.data_path
    id = args.id
    offset = args.offset
    rank_ensemble_models(data_path, id, offset)
def arg_parser():
    import argparse
    parser = argparse.ArgumentParser(description="Rank ensemble models")
    parser.add_argument("--data_path", type=str, help="Path to the data directory")
    parser.add_argument("--id", type=str, help="ID of the dataset")
    parser.add_argument("--offset", type=int, default=0, help="Offset of the mutation site")
    return parser

if __name__ == "__main__":
    args = arg_parser().parse_args()
    main(args)
    # print("Done!")
    # data_path = "./DATASET"
    # id = "D7M15"
    # offset = 8#0 offset -1
    # rank_ensemble_models(data_path, id, offset)

    
    