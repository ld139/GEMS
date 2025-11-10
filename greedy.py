import os
import argparse
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import heapq
import re
import math
from collections import Counter

def arg_parser():
    parser = argparse.ArgumentParser(description='Greedy selection of mutants balancing activity and diversity.')
    parser.add_argument('--DMS_id', type=str, required=True, help='ID of DMS used to locate input files')
    parser.add_argument('--data_dir', type=str, default="DATASET", help='Root directory containing DMS datasets')
    parser.add_argument('--offset', type=int, default=1, help='Position offset (mutation numbering vs Python indexing)')
    parser.add_argument('--subset_size', type=int, default=50, help='Number of mutants to greedily select')
    parser.add_argument('--weight', type=float, default=0.55, help='Weight w: final score = w*activity + (1-w)*diversity')
    parser.add_argument('--top_n', type=int, default=1000, help='Pre-filter: keep only top-N by fitness before selection')
    parser.add_argument('--mode', choices=['single', 'multi'], default='single',
                        help='Mutation mode: single (specialized schemes) or multi (original entropy logic)')
    parser.add_argument('--single_strategy', choices=['coverage', 'entropy', 'hybrid'], default='coverage',
                        help='Strategy for single-mutant diversity: coverage | entropy | hybrid')
    parser.add_argument('--new_pos_reward', type=float, default=1.0,
                        help='Reward for selecting a mutant from a previously unseen position (single coverage strategy)')
    parser.add_argument('--repeat_pos_decay', type=float, default=0.3,
                        help='Base reward factor for already covered position; actual reward = repeat_pos_decay/(1+count)')
    parser.add_argument('--entropy_weight', type=float, default=0.5,
                        help='Weight of entropy component inside hybrid single strategy (0..1)')
    return parser.parse_args()

# ---------------- Parsing utilities ----------------

def parse_mutation_string(mut_str):
    """
    Parse mutation string like 'A42G' or 'A42G:Y50F' into {position: mutated_amino_acid}.
    """
    mutations = {}
    parts = mut_str.split(':')
    for part in parts:
        match = re.match(r'([A-Z])(\d+)([A-Z])', part)
        if match:
            _, position, mut_aa = match.groups()
            mutations[int(position)] = mut_aa
    return mutations

def is_single_mutant(mut_str):
    """
    Return True if mutation string represents exactly one mutated position.
    """
    return len(parse_mutation_string(mut_str)) == 1

def get_all_positions(total_library):
    """
    Collect all unique mutation positions across the library.
    """
    all_positions = set()
    for mut_str in total_library:
        all_positions.update(parse_mutation_string(mut_str).keys())
    return sorted(all_positions)

# ---------------- Position string conversion ----------------

def convert_to_position_string(mut_str, positions, wild_type_sequence, offset):
    """
    For a given mutant, return a string with length = len(positions).
    Mutated position: mutated AA; otherwise wild-type AA; '.' if out of bounds.
    """
    mutations = parse_mutation_string(mut_str)
    chars = []
    for pos in positions:
        if pos in mutations:
            chars.append(mutations[pos])
        else:
            idx = pos - offset
            if 0 <= idx < len(wild_type_sequence):
                chars.append(wild_type_sequence[idx])
            else:
                chars.append('.')
    return ''.join(chars)

# ---------------- Multi-mutation entropy diversity ----------------

def calculate_subset_diversity_multi(position_strings):
    """
    Shannon entropy per site including wild-type residues; sum over sites.
    """
    if not position_strings:
        return 0.0
    M = len(position_strings[0])
    total_count = len(position_strings)
    diversity = 0.0
    for i in range(M):
        site_i = [seq[i] for seq in position_strings]
        aa_counts = {}
        for aa in site_i:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
        H = 0.0
        for count in aa_counts.values():
            p = count / total_count
            if p > 0:
                H -= p * math.log2(p)
        diversity += H
    return diversity

def get_max_diversity_multi(position_strings):
    """
    Theoretical max entropy sum (uniform over observed residues per site).
    """
    if not position_strings:
        return 0.0
    M = len(position_strings[0])
    max_entropy = 0.0
    for i in range(M):
        site_residues = set(seq[i] for seq in position_strings)
        k = len(site_residues)
        if k <= 1:
            continue
        p = 1.0 / k
        max_entropy += -k * (p * math.log2(p))
    return max_entropy

# ---------------- Single-mutation entropy diversity (existing scheme) ----------------

def calculate_subset_diversity_single_entropy(position_strings, wild_type_chars):
    """
    Entropy-only single scheme:
    For each site, compute entropy over mutated residues only (exclude wild-type),
    multiply by coverage ratio, sum across sites.
    """
    if not position_strings:
        return 0.0
    total = len(position_strings)
    M = len(position_strings[0])
    diversity = 0.0
    for i in range(M):
        wt = wild_type_chars[i]
        mutated = [seq[i] for seq in position_strings if seq[i] != wt]
        if not mutated:
            continue
        counts = Counter(mutated)
        mutated_total = len(mutated)
        H = 0.0
        for c in counts.values():
            p = c / mutated_total
            H -= p * math.log2(p)
        coverage = mutated_total / total
        diversity += coverage * H
    return diversity

def get_max_diversity_single_entropy(position_strings, wild_type_chars):
    """
    Max theoretical entropy sum (coverage=1, uniform mutated residues).
    """
    if not position_strings:
        return 0.0
    M = len(position_strings[0])
    max_div = 0.0
    for i in range(M):
        wt = wild_type_chars[i]
        aas = set(s[i] for s in position_strings if s[i] != wt)
        k = len(aas)
        if k <= 1:
            continue
        p = 1.0 / k
        max_div += -k * (p * math.log2(p))
    return max_div

# ---------------- Single-mutation coverage gain ----------------

def compute_coverage_gain_for_single(mut_str, selected_subset, offset, new_pos_reward, repeat_pos_decay):
    """
    Reward for a candidate based on position coverage:
    - New position: new_pos_reward
    - Existing position: repeat_pos_decay / (1 + count_at_position)
    """
    mutations = parse_mutation_string(mut_str)
    # Single mutant assumption
    pos = list(mutations.keys())[0]

    # Count occurrences in selected subset
    counts = {}
    for s in selected_subset:
        p = list(parse_mutation_string(s).keys())[0]
        counts[p] = counts.get(p, 0) + 1

    if pos not in counts:
        return new_pos_reward
    else:
        return repeat_pos_decay / (1 + counts[pos])

def calculate_entropy_component_for_subset(selected_subset, positions, wild_type_sequence, offset, wild_type_chars):
    """
    Compute entropy-based diversity on current subset (used in hybrid).
    """
    position_strings = [
        convert_to_position_string(s, positions, wild_type_sequence, offset)
        for s in selected_subset
    ]
    return calculate_subset_diversity_single_entropy(position_strings, wild_type_chars)

# ---------------- Utility ----------------

def normalize_value(value, min_val, max_val):
    """
    Normalize to [0,1]; return 0 if max==min.
    """
    if max_val == min_val:
        return 0.0
    return (value - min_val) / (max_val - min_val)

def evaluate_activity_only(seq, selected_subset, remaining_scores):
    """
    Mean activity for temporary subset (selected + seq).
    """
    temp_subset = selected_subset + [seq]
    scores = [remaining_scores[s] for s in temp_subset]
    return float(np.mean(scores))

# ---------------- Main selection routine ----------------

def select_subset_optimized(total_library,
                            activity_scores,
                            wild_type_sequence,
                            offset,
                            subset_size=500,
                            w=0.55,
                            top_n=1000,
                            mode='single',
                            single_strategy='coverage',
                            new_pos_reward=1.0,
                            repeat_pos_decay=0.3,
                            entropy_weight=0.5):
    """
    Greedy selection with activity/diversity tradeoff.
    Single-mode strategies:
      coverage: diversity = coverage_gain
      entropy: diversity = mutated-only entropy * coverage
      hybrid: diversity = normalized combination of coverage_gain and entropy
    """
    if not total_library:
        return []

    positions = get_all_positions(total_library)
    if not positions:
        raise ValueError("No mutation positions parsed.")

    position_strings_full = [
        convert_to_position_string(m, positions, wild_type_sequence, offset)
        for m in total_library
    ]
    lengths = {len(s) for s in position_strings_full}
    if len(lengths) != 1:
        raise ValueError(f"Inconsistent position string lengths: {lengths}")

    # Build scores mapping
    if isinstance(activity_scores, list):
        if len(activity_scores) != len(total_library):
            raise ValueError("Mismatch between total_library and activity_scores length.")
        remaining_scores = {m: s for m, s in zip(total_library, activity_scores)}
    elif isinstance(activity_scores, dict):
        remaining_scores = activity_scores.copy()
    else:
        raise ValueError("activity_scores must be list or dict.")

    for mutant, score in remaining_scores.items():
        if not isinstance(score, (int, float)) or pd.isna(score):
            raise ValueError(f"Invalid fitness score for {mutant}: {score}")

    # Pre-filter by top_n
    if len(remaining_scores) > top_n:
        top_sequences = heapq.nlargest(top_n, remaining_scores.items(), key=lambda x: x[1])
        remaining_library = [m for m, _ in top_sequences]
        remaining_scores = {m: s for m, s in top_sequences}
    else:
        remaining_library = total_library.copy()

    selected_subset = []
    all_activities = list(remaining_scores.values())
    min_activity, max_activity = min(all_activities), max(all_activities)

    # Wild-type chars aligned to positions
    wild_type_chars = []
    for pos in positions:
        idx = pos - offset
        if 0 <= idx < len(wild_type_sequence):
            wild_type_chars.append(wild_type_sequence[idx])
        else:
            wild_type_chars.append('.')

    # Diversity bounds per mode
    if mode == 'multi':
        max_diversity = get_max_diversity_multi(position_strings_full)
        min_diversity = 0.0
    else:
        if single_strategy == 'entropy':
            max_diversity = get_max_diversity_single_entropy(position_strings_full, wild_type_chars)
            min_diversity = 0.0
        elif single_strategy == 'coverage':
            # coverage_gain raw max is new_pos_reward; single-threat normalization uses that
            max_diversity = new_pos_reward
            min_diversity = 0.0
        else:  # hybrid
            # we will dynamically recompute combined raw; need provisional max
            # approximate hybrid max: (1 - entropy_weight)*new_pos_reward + entropy_weight*max_entropy
            max_entropy = get_max_diversity_single_entropy(position_strings_full, wild_type_chars)
            max_diversity = (1 - entropy_weight) * new_pos_reward + entropy_weight * max_entropy
            min_diversity = 0.0

    for _ in range(subset_size):
        if not remaining_library:
            break
        best_seq = None
        best_score = -float('inf')

        if mode == 'multi':
            # Evaluate all candidates in parallel
            activities_divs = Parallel(n_jobs=-1)(
                delayed(_eval_multi)(
                    seq, selected_subset, remaining_scores, positions, wild_type_sequence, offset
                )
                for seq in remaining_library
            )
            for seq, (temp_activity, temp_diversity) in zip(remaining_library, activities_divs):
                norm_activity = normalize_value(temp_activity, min_activity, max_activity)
                norm_diversity = normalize_value(temp_diversity, min_diversity, max_diversity)
                score = w * norm_activity + (1 - w) * norm_diversity
                if score > best_score:
                    best_score = score
                    best_seq = seq

        else:
            # Single-mode evaluation: parallel only for activity
            activities = Parallel(n_jobs=-1)(
                delayed(evaluate_activity_only)(seq, selected_subset, remaining_scores)
                for seq in remaining_library
            )

            # Precompute entropy component if needed (entropy / hybrid)
            if single_strategy in ('entropy', 'hybrid'):
                # We need entropy for each candidate subset -> recompute per candidate (costly).
                # Optimization: compute current entropy once; then for candidate, update by addition.
                # Simplicity: full recompute for correctness (library size moderate).
                current_entropy = calculate_entropy_component_for_subset(
                    selected_subset, positions, wild_type_sequence, offset, wild_type_chars
                )
            else:
                current_entropy = 0.0  # unused

            for seq, temp_activity in zip(remaining_library, activities):
                # Coverage gain
                coverage_gain = compute_coverage_gain_for_single(
                    seq, selected_subset, offset, new_pos_reward, repeat_pos_decay
                )

                if single_strategy == 'coverage':
                    temp_div_raw = coverage_gain
                elif single_strategy == 'entropy':
                    # Build temp subset and compute entropy afresh
                    temp_subset = selected_subset + [seq]
                    temp_entropy = calculate_entropy_component_for_subset(
                        temp_subset, positions, wild_type_sequence, offset, wild_type_chars
                    )
                    temp_div_raw = temp_entropy
                else:  # hybrid
                    temp_subset = selected_subset + [seq]
                    temp_entropy = calculate_entropy_component_for_subset(
                        temp_subset, positions, wild_type_sequence, offset, wild_type_chars
                    )
                    # Combined raw diversity before normalization
                    temp_div_raw = (1 - entropy_weight) * coverage_gain + entropy_weight * temp_entropy

                norm_activity = normalize_value(temp_activity, min_activity, max_activity)
                norm_diversity = normalize_value(temp_div_raw, min_diversity, max_diversity)
                score = w * norm_activity + (1 - w) * norm_diversity

                if score > best_score:
                    best_score = score
                    best_seq = seq

        if best_seq is None:
            break
        selected_subset.append(best_seq)
        remaining_library.remove(best_seq)

    return selected_subset

def _eval_multi(seq, selected_subset, remaining_scores, positions, wild_type_sequence, offset):
    """
    Helper for multi-mode parallel evaluation.
    """
    temp_subset = selected_subset + [seq]
    scores = [remaining_scores[s] for s in temp_subset]
    temp_activity = float(np.mean(scores))
    position_strings = [
        convert_to_position_string(s, positions, wild_type_sequence, offset)
        for s in temp_subset
    ]
    temp_diversity = calculate_subset_diversity_multi(position_strings)
    return temp_activity, temp_diversity

# ---------------- Main entry ----------------

def main():
    args = arg_parser()
    DMS_id = args.DMS_id
    data_dir = args.data_dir
    offset = args.offset
    subset_size = args.subset_size
    weight = args.weight
    top_n = args.top_n
    mode = args.mode
    single_strategy = args.single_strategy
    new_pos_reward = args.new_pos_reward
    repeat_pos_decay = args.repeat_pos_decay
    entropy_weight = args.entropy_weight

    fasta_path = os.path.join(data_dir, DMS_id, f'{DMS_id}.fasta')
    input_csv = os.path.join(data_dir, DMS_id, f'{DMS_id}_rank.csv')
    output_csv = os.path.join(data_dir, DMS_id, f'greedy_{DMS_id}.csv')

    # Read wild-type sequence
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
        wild_type_sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))

    # Read mutant fitness data
    df = pd.read_csv(input_csv, dtype={'mutant': str, 'fitness': float})
    total_library = df['mutant'].tolist()
    activity_scores = dict(zip(df['mutant'], df['fitness']))

    all_single = all(is_single_mutant(m) for m in total_library)
    any_multi = any(not is_single_mutant(m) for m in total_library)

    if mode == 'single' and any_multi:
        print("WARNING: mode='single' but multi-position mutants detected. Single strategies may be suboptimal.")
    if mode == 'multi' and all_single:
        print("INFO: mode='multi' but all mutants are single; diversity may be weakly discriminative.")

    selected_subset = select_subset_optimized(
        total_library=total_library,
        activity_scores=activity_scores,
        wild_type_sequence=wild_type_sequence,
        offset=offset,
        subset_size=subset_size,
        w=weight,
        top_n=top_n,
        mode=mode,
        single_strategy=single_strategy,
        new_pos_reward=new_pos_reward,
        repeat_pos_decay=repeat_pos_decay,
        entropy_weight=entropy_weight
    )

    # Preserve order of selection
    result_df = df[df['mutant'].isin(selected_subset)].copy()
    result_df = result_df.set_index('mutant').loc[selected_subset].reset_index()
    result_df.to_csv(output_csv, index=False)

    print(f"Selected subset saved to {output_csv}")
    print(f"Total selected: {len(selected_subset)}")
    print(f"Mode: {mode} | Single strategy: {single_strategy}")
    if mode == 'single' and single_strategy == 'coverage':
        # Simple coverage stats
        pos_counts = {}
        for m in selected_subset:
            p = list(parse_mutation_string(m).keys())[0]
            pos_counts[p] = pos_counts.get(p, 0) + 1
        unique_positions = len(pos_counts)
        print(f"Unique positions covered: {unique_positions}")
        print(f"Position frequency (first 20): {list(pos_counts.items())[:20]}")

if __name__ == '__main__':
    main()