import itertools
import os
import csv
import argparse

def get_Saturation_Mutagenesis(wt_fasta, mute_sites, offset):
    """
    生成单点饱和突变列表
    参数:
        wt_fasta (str): 野生型氨基酸序列文件路径
        mute_sites (list): 突变点位置列表（1-based 索引）
        offset (int): 偏移量，将突变点位置调整为 0-based
    返回:
        combinations (list): 单点突变组合列表 (A3R)
        combinations_str (list): 偏移后的单点突变组合 (A4R)
        mutated_seqs (list): 突变后的序列列表
    """
    # 读取野生型氨基酸序列
    with open(wt_fasta, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                wt_seq = line.strip()
                break

    # 定义 20 种氨基酸
    amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 
                   'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    combinations = []       # 突变组合，如 A3R
    combinations_str = []   # 偏移后的突变组合，如 A4R
    mutated_seqs = []       # 突变后的序列列表
    mute_sites = sorted(mute_sites)  # 突变位点排序
    
    # 遍历突变位点，逐一生成单点突变
    for site in mute_sites:
        site_0_based = site - offset  # 调整为 0-based 索引
        if site_0_based < 0 or site_0_based >= len(wt_seq):
            raise ValueError(f"突变位点 {site} 超出序列长度范围")
        
        wt_residue = wt_seq[site_0_based]  # 获取野生型残基

        for aa in amino_acids:
            if aa != wt_residue:  # 避免自我突变
                # 生成突变组合
                mutation = f"{wt_residue}{site-offset+1}{aa}"
                mutation_str = f"{wt_residue}{site+1}{aa}"
                combinations.append(mutation)
                combinations_str.append(mutation_str)
                
                # 生成突变后的序列
                mutated_seq = list(wt_seq)
                mutated_seq[site_0_based] = aa
                mutated_seqs.append(''.join(mutated_seq))

    return combinations, combinations_str, mutated_seqs

def get_Combinatorial_Mutagenesis(wt_fasta, mute_sites, offset):
    """
    生成组合饱和突变（包括单突变和多突变）
    
    参数:
        wt_fasta (str): 野生型氨基酸序列文件路径
        mute_sites (list): 突变点位置列表（1-based 索引）
        offset (int): 偏移量，将突变点位置调整为 0-based
    
    返回:
        combinations (list): 组合突变列表（每个元素是单个突变字符串列表）
        combinations_str (list): 偏移后的组合突变字符串（用冒号分隔）
        mutated_seqs (list): 突变后的序列列表
    """
    # 读取野生型氨基酸序列
    with open(wt_fasta, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                wt_seq = line.strip()
                break

    # 定义 20 种氨基酸
    amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 
                   'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    combinations = []       # 突变组合列表（每个元素是单个突变字符串列表）
    combinations_str = []   # 偏移后的组合突变字符串（用冒号分隔）
    mutated_seqs = []       # 突变后的序列列表
    mute_sites = sorted(mute_sites)  # 突变位点排序
    
    # 准备位点信息
    site_info = []
    for site in mute_sites:
        site_0_based = site - offset  # 调整为 0-based 索引
        if site_0_based < 0 or site_0_based >= len(wt_seq):
            raise ValueError(f"突变位点 {site} 超出序列长度范围")
        
        wt_residue = wt_seq[site_0_based]  # 获取野生型残基
        possible_aas = [aa for aa in amino_acids if aa != wt_residue]
        site_info.append((site, site_0_based, wt_residue, possible_aas))
    
    # 生成所有可能的组合（包括单突变和多突变）
    # 使用所有非空子集（从1个位点到所有位点）
    for subset_size in range(1, len(site_info) + 1):
        for site_subset in itertools.combinations(site_info, subset_size):
            # 获取该子集的所有可能氨基酸组合
            mutations_per_site = [site[3] for site in site_subset]
            
            # 使用笛卡尔积生成所有组合
            for combo in itertools.product(*mutations_per_site):
                mutated_seq = list(wt_seq)  # 转换为列表以便修改
                mutation_list = []          # 存储单个突变字符串
                mutation_str_list = []      # 存储偏移后的单个突变字符串
                
                # 应用所有突变
                for (site, site_0_based, wt_residue, _), new_aa in zip(site_subset, combo):
                    mutated_seq[site_0_based] = new_aa
                    mutation = f"{wt_residue}{site-offset+1}{new_aa}"
                    mutation_str = f"{wt_residue}{site+1}{new_aa}"
                    mutation_list.append(mutation)
                    mutation_str_list.append(mutation_str)
                
                # 转换为字符串
                mutated_seq_str = ''.join(mutated_seq)
                combo_str = ':'.join(mutation_list)
                combo_str_offset = ':'.join(mutation_str_list)
                
                # 保存结果
                combinations.append(mutation_list)
                combinations_str.append(combo_str_offset)
                mutated_seqs.append(mutated_seq_str)
    
    return combinations, combinations_str, mutated_seqs

def save_mutation_combinations_to_file(wt_id, combinations, combinations_str, mutated_seqs, output_dir):
    """
    将所有突变组合保存到文件
    参数:
        combinations (list): 突变组合列表
        output_dir (str): 输出目录
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 保存突变组合信息到文本文件
    with open(os.path.join(output_dir, f'{wt_id}.txt'), 'w') as f:
        for combo in combinations:
            if isinstance(combo, list):  # 组合突变
                f.write(':'.join(combo) + '\n')
            else:  # 单点突变
                f.write(combo + '\n')
    
    # 保存到CSV文件（包含偏移）
    with open(os.path.join(output_dir, f'{wt_id}_str.csv'), 'w') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["mutant", "mutated_sequence", "DMS_score"])
        for combo, mutated_seq in zip(combinations_str, mutated_seqs):
            writer.writerow([combo, mutated_seq, '0.0'])
    
    # 保存到CSV文件（不含偏移）
    with open(os.path.join(output_dir, f'{wt_id}.csv'), 'w') as f:
        writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["mutant", "mutated_sequence", "DMS_score"])
        for combo, mutated_seq in zip(combinations, mutated_seqs):
            if isinstance(combo, list):  # 组合突变
                writer.writerow([':'.join(combo), mutated_seq, '0.0'])
            else:  # 单点突变
                writer.writerow([combo, mutated_seq, '0.0'])
    
    # 生成统计信息
    num_mutants = len(mutated_seqs)
    with open(os.path.join(output_dir, "summary.txt"), 'w') as f:
        f.write(f"Total mutants generated: {num_mutants}\n")
        if num_mutants > 0:
            if isinstance(combinations[0], list):  # 组合突变
                # 计算最大突变点数
                max_mutations = max(len(c) for c in combinations)
                min_mutations = min(len(c) for c in combinations)
                f.write(f"Combinatorial mutagenesis with {min_mutations} to {max_mutations} mutations\n")
                f.write(f"Sites: {', '.join(str(site) for site in combinations[0])}\n")
            else:  # 单点突变
                f.write(f"Single-site mutagenesis at {len(set(c.split('_')[0][1:-1] for c in combinations))} sites\n")

def read_fasta(file_path):
    """读取FASTA文件并返回序列字符串"""
    sequence = ""
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue  # 跳过头部行
            sequence += line.strip()
    return sequence

def arg_parser():
    parser = argparse.ArgumentParser(description='Generate saturation mutagenesis combinations')
    parser.add_argument('-wt', '--wildtype', required=True, help='wildtype fasta file')
    parser.add_argument('-sites', '--sites', required=False, 
                        help='mutation sites file OR comma-separated list of sites (e.g., "10,20,30")')
    parser.add_argument('-offset', '--offset', required=True, type=int, help='offset of mutation sites')
    parser.add_argument('-o', '--output', required=True, help='output directory')
    parser.add_argument('--all_sites', action='store_true', 
                        help='perform saturation mutagenesis on all sites in the sequence')
    parser.add_argument('-comb','--combinatorial', type=str, required=False,
                        help='comma-separated list of sites for combinatorial saturation mutagenesis (includes single and multiple mutations)')
    return parser.parse_args()

def main():
    args = arg_parser()
    
    # 读取野生型序列
    wt_seq = read_fasta(args.wildtype)
    wt_id = os.path.basename(args.wildtype).split(".fasta")[0]
    
    # 处理突变位点
    if args.combinatorial:
        # 组合突变模式（包括单突变和多突变）
        mute_sites = [int(site) for site in args.combinatorial.split(',')]
        combinations, combinations_str, mutated_seqs = get_Combinatorial_Mutagenesis(
            args.wildtype, mute_sites, int(args.offset)
        )
    elif args.all_sites:
        # 全序列突变：生成1到序列长度的位点列表
        mute_sites = list(range(1, len(wt_seq) + 1))
        combinations, combinations_str, mutated_seqs = get_Saturation_Mutagenesis(
            args.wildtype, mute_sites, int(args.offset)
        )
    elif args.sites:
        # 直接输入位点（文件或逗号分隔列表）
        if os.path.isfile(args.sites):
            with open(args.sites, 'r') as f:
                mute_sites = [int(line.strip()) for line in f]
        else:
            mute_sites = [int(site) for site in args.sites.split(',')]
        combinations, combinations_str, mutated_seqs = get_Saturation_Mutagenesis(
            args.wildtype, mute_sites, int(args.offset)
        )
    else:
        raise ValueError("Must provide mutation sites (--sites), use --all_sites, or specify combinatorial sites (--combinatorial)")
    
    # 保存结果
    save_mutation_combinations_to_file(
        wt_id, combinations, combinations_str, mutated_seqs, args.output
    )

if __name__ == '__main__':
    main()