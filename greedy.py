import itertools
import os
import csv
import argparse
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import heapq
import re

def arg_parser():
    parser = argparse.ArgumentParser(description='Generate saturation mutagenesis combinations')
    parser.add_argument('--DMS_id', type=str, required=True, help='id of DMS in DMS reference file')
    parser.add_argument('--data_dir', type=str, default="DATASET", help='input directory')
    parser.add_argument('-offset', type=int, help='offset of mutation sites', default=1)
    parser.add_argument('--subset_size', type=int, default=20, help='size of subset to consider at each iteration')
    parser.add_argument('--weight', type=float, default=0.55, help='weight for combining scores')
    parser.add_argument('--top_n', type=int, default=1000, help='number of top mutants to select')

    return parser.parse_args()

import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import heapq
import re

# Define functions (same as provided, with validation)
def parse_mutation_string(mut_str):
    mutations = {}
    parts = mut_str.split(':')
    for part in parts:
        match = re.match(r'([A-Z])(\d+)([A-Z])', part)
        if match:
            orig_aa, position, mut_aa = match.groups()
            mutations[int(position)] = mut_aa
    return mutations

def get_all_positions(total_library):
    all_positions = set()
    for mut_str in total_library:
        mutations = parse_mutation_string(mut_str)
        all_positions.update(mutations.keys())
    return sorted(all_positions)

def convert_to_position_string(mut_str, positions, wild_type_sequence, offset):
    mutations = parse_mutation_string(mut_str)
    position_string = []
    
    for pos in positions:
        if pos in mutations:
            position_string.append(mutations[pos])
        else:
            idx = pos - offset
            if 0 <= idx < len(wild_type_sequence):
                position_string.append(wild_type_sequence[idx])
            else:
                position_string.append('.')
    
    result = ''.join(position_string)
    return result

def calculate_subset_diversity_optimized(position_strings):
    if not position_strings:
        return 0.0
    
    M = len(position_strings[0])
    total_count = len(position_strings)
    
    if total_count == 0:
        raise ValueError("Empty position_strings list")
    
    diversity = 0.0
    for i in range(M):
        site_i = [seq[i] for seq in position_strings]
        aa_counts = {}
        for aa in site_i:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
        
        H_p_i = 0.0
        for aa, count in aa_counts.items():
            p_i = count / total_count
            if p_i > 0:
                H_p_i -= p_i * np.log2(p_i)
        
        diversity += H_p_i
    
    return diversity

def normalize_value(value, min_val, max_val):
    if not all(isinstance(x, (int, float)) for x in [value, min_val, max_val]):
        raise ValueError(f"Non-numeric inputs to normalize_value: value={value}, min_val={min_val}, max_val={max_val}")
    if max_val == min_val:
        return 0.0
    return (value - min_val) / (max_val - min_val)

def evaluate_candidate(seq, selected_subset, remaining_scores, positions, wild_type_sequence, offset):
    temp_subset = selected_subset + [seq]
    scores = [remaining_scores.get(s, 0.0) for s in temp_subset]
    if not all(isinstance(score, (int, float)) for score in scores):
        raise ValueError(f"Non-numeric scores found in temp_subset: {scores}")
    temp_activity = np.mean(scores)
    
    position_strings = [
        convert_to_position_string(s, positions, wild_type_sequence, offset) 
        for s in temp_subset
    ]
    temp_diversity = calculate_subset_diversity_optimized(position_strings)
    
    return temp_activity, temp_diversity

def get_max_diversity(position_strings):
    if not position_strings:
        return 0.0
    
    M = len(position_strings[0])
    max_entropy = 0.0
    for i in range(M):
        site_amino_acids = set(seq[i] for seq in position_strings)
        num_amino_acids = len(site_amino_acids)
        if num_amino_acids <= 1:
            continue
        p = 1.0 / num_amino_acids
        max_entropy_i = -num_amino_acids * (p * np.log2(p))
        max_entropy += max_entropy_i
    
    return max_entropy

def select_subset_optimized(total_library, activity_scores, wild_type_sequence, offset, 
                           subset_size=500, w=0.55, top_n=30000):
    if not total_library:
        return []
    
    positions = get_all_positions(total_library)
    if not positions:
        raise ValueError("无法从突变字符串中解析到位点位置")
    
    position_strings = [
        convert_to_position_string(mut_str, positions, wild_type_sequence, offset) 
        for mut_str in total_library
    ]
    lengths = [len(s) for s in position_strings]
    if len(set(lengths)) > 1:
        raise ValueError(f"Inconsistent position string lengths: {lengths}")
    
    if len(total_library) != len(activity_scores):
        raise ValueError("total_library 和 activity_scores 的长度不匹配")
    
    if isinstance(activity_scores, list):
        remaining_scores = {mut_str: score for mut_str, score in zip(total_library, activity_scores)}
    elif isinstance(activity_scores, dict):
        remaining_scores = activity_scores.copy()
    else:
        raise ValueError("activity_scores 必须是列表或字典")
    
    for mutant, score in remaining_scores.items():
        if not isinstance(score, (int, float)) or pd.isna(score):
            raise ValueError(f"Non-numeric or NaN fitness score found for mutant {mutant}: {score}")
    
    if len(remaining_scores) > top_n:
        top_sequences = heapq.nlargest(top_n, remaining_scores.items(), key=lambda x: x[1])
        remaining_library = [mut_str for mut_str, _ in top_sequences]
        remaining_scores = {mut_str: score for mut_str, score in top_sequences}
    else:
        remaining_library = total_library.copy()
    
    selected_subset = []
    all_activities = list(remaining_scores.values())
    if not all_activities:
        raise ValueError("no activity scores available to determine min and max")
    min_activity, max_activity = min(all_activities), max(all_activities)
    
    max_diversity = get_max_diversity(position_strings)
    min_diversity = 0.0
    
    for _ in range(subset_size):
        best_score = -float('inf')
        best_seq = None
        
        results = Parallel(n_jobs=-1)(
            delayed(evaluate_candidate)(
                seq, selected_subset, remaining_scores, positions, wild_type_sequence, offset
            ) 
            for seq in remaining_library
        )
        
        for seq, (temp_activity, temp_diversity) in zip(remaining_library, results):
            if seq not in remaining_scores:
                continue
            norm_activity = normalize_value(temp_activity, min_activity, max_activity)
            norm_diversity = normalize_value(temp_diversity, min_diversity, max_diversity)
            if not isinstance(norm_activity, (int, float)) or not isinstance(norm_diversity, (int, float)):
                raise ValueError(f"Non-numeric normalized values: norm_activity={norm_activity}, norm_diversity={norm_diversity}")
            score = w * norm_activity + (1 - w) * norm_diversity
            
            if score > best_score:
                best_score = score
                best_seq = seq
        
        if best_seq:
            selected_subset.append(best_seq)
            remaining_library.remove(best_seq)
        else:
            break
    
    return selected_subset

def main():
    args = arg_parser()
    DMS_id = args.DMS_id
    # wt_fasta = args.wildtype
    offset = args.offset
    input_csv = args.data_dir + os.sep + DMS_id + os.sep + f'{DMS_id}_rank.csv'
    output_dir = args.data_dir + os.sep + DMS_id + os.sep + "greedy_" + f'{DMS_id}.csv'
    subset_size = args.subset_size
    # print(subset_size)
    # print(args.weight)
    weight = args.weight
    top_n = args.top_n
    wt_fasta = args.data_dir + os.sep + DMS_id + os.sep + f'{DMS_id}.fasta'
    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir)
    
    with open(wt_fasta, 'r') as f:
        lines = f.readlines()
        wild_type_sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    
    df = pd.read_csv(input_csv, dtype={'mutant': str, 'fitness': float})
    total_library = df['mutant'].tolist()
    activity_scores = dict(zip(df['mutant'], df['fitness']))
    
    selected_subset = select_subset_optimized(
        total_library, activity_scores, wild_type_sequence, offset,
        subset_size=subset_size, w=weight, top_n=top_n
    )
    
    result_df = df[df['mutant'].isin(selected_subset)].copy()
    result_df = result_df.set_index('mutant').loc[selected_subset].reset_index()
    
    output_csv = os.path.join(output_dir)
    result_df.to_csv(output_csv, index=False)
    print(f"Selected subset saved to {output_csv}")

if __name__ == '__main__':
    main()
