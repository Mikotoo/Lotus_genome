from Bio import AlignIO
from collections import Counter
import argparse
import os

try:
    from Levenshtein import distance as levenshtein_distance
except ImportError:
    def levenshtein_distance(s1, s2):
        """纯Python实现的Levenshtein距离计算（适用于短序列）"""
        if len(s1) < len(s2):
            return levenshtein_distance(s2, s1)
        if len(s2) == 0:
            return len(s1)
        
        previous_row = range(len(s2) + 1)
        for i, c1 in enumerate(s1):
            current_row = [i + 1]
            for j, c2 in enumerate(s2):
                insertions = previous_row[j + 1] + 1
                deletions = current_row[j] + 1
                substitutions = previous_row[j] + (c1 != c2)
                current_row.append(min(insertions, deletions, substitutions))
            previous_row = current_row
        return previous_row[-1]

def calculate_consensus(alignment, threshold=0.5, include_gap=False, ambiguity=False):
    """增强型共有序列计算"""
    consensus = []
    valid_bases = set('ACGTUacgtu-') if include_gap else set('ACGTUacgtu')
    iupac_code = {
        frozenset(['A']): 'A', frozenset(['C']): 'C', frozenset(['G']): 'G',
        frozenset(['T']): 'T', frozenset(['U']): 'U', frozenset(['A','G']): 'R',
        frozenset(['C','T']): 'Y', frozenset(['G','C']): 'S', frozenset(['A','T']): 'W',
        frozenset(['G','T']): 'K', frozenset(['A','C']): 'M', frozenset(['C','G','T']): 'B',
        frozenset(['A','G','T']): 'D', frozenset(['A','C','T']): 'H', frozenset(['A','C','G']): 'V',
        frozenset(['A','C','G','T']): 'N'
    }

    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx].upper()
        filtered = [b for b in column if b in valid_bases] or ['N']  # 处理全无效列
        
        counts = Counter(filtered)
        total = len(filtered)
        
        candidates = []
        current_sum = 0.0
        for base, count in counts.most_common():
            freq = count / total
            if current_sum + freq > threshold:
                if current_sum == 0:  # 单个碱基超过阈值
                    candidates = [base]
                else:  # 多个碱基组合超过阈值
                    candidates.append(base)
                break
            candidates.append(base)
            current_sum += freq

        if ambiguity:
            key = frozenset(c.upper() for c in candidates)
            consensus.append(iupac_code.get(key, 'N'))
        else:
            consensus.append(candidates[0] if candidates else 'N')

    return ''.join(consensus)

def main():
    parser = argparse.ArgumentParser(
        description='Generate consensus sequence and calculate distances',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True, help='Input aligned FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output consensus FASTA file')
    parser.add_argument('-d', '--distances', required=True, help='Output TSV file for distances')
    parser.add_argument('-t', '--threshold', type=float, default=0.5,
                        help='Consensus frequency threshold')
    parser.add_argument('--include-gap', action='store_true',
                        help='Include gaps in consensus calculation')
    parser.add_argument('--ambiguity', action='store_true',
                        help='Use IUPAC ambiguity codes')
    parser.add_argument('--name', help='Custom consensus name')

    args = parser.parse_args()

    # 读取比对文件并生成共识序列
    try:
        alignment = AlignIO.read(args.input, "fasta")
    except Exception as e:
        raise SystemExit(f"Error reading alignment: {str(e)}")

    consensus_name = args.name if args.name else f"Consensus_{os.path.basename(args.input)}"
    consensus_seq = calculate_consensus(
        alignment,
        threshold=args.threshold,
        include_gap=args.include_gap,
        ambiguity=args.ambiguity
    )

    # 写入共识序列

    consensus_clean = consensus_seq.replace('-', '').upper()

    with open(args.output, 'w') as f:
        f.write(f">{consensus_name}\n")
        for i in range(0, len(consensus_clean), 80):
            f.write(consensus_clean[i:i+80] + "\n")

    # 计算并保存距离数据
    with open(args.distances, 'w') as f_dist:
        f_dist.write("SequenceID\tEffectiveLength\tDistance\tNormalizedDistance\n")
        
        for record in alignment:
            # 去除gap并转为大写
            seq_clean = str(record.seq).upper().replace('-', '')
            
            # 计算实际参数
            effective_length = len(seq_clean)
            distance = levenshtein_distance(consensus_clean, seq_clean)
            normalized = distance / effective_length if effective_length > 0 else 0
            
            f_dist.write(f"{record.id}\t{effective_length}\t{distance}\t{normalized:.4f}\n")

if __name__ == "__main__":
    main()
