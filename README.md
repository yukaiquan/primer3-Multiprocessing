# 1 SSR引物批量设计
##### 🙈: yukaiquan
##### 📧: 1962568272@qq.com
## 1.1 通过misa结果制作bed文件以及剔除一些连续多段的SSR
那么写一个python脚本吧
```python
# yukaiquan
# email:1962568272@qq.com
import pandas as pd
import sys


input_file = sys.argv[1]
output_file = sys.argv[2]
output_bed = sys.argv[3]
misa_df = pd.read_csv(input_file, sep='\t')
misa_df = misa_df[misa_df['SSR type'].isin(['p2', 'p3', 'p4', 'p5', 'p6'])]
# output
misa_df.to_csv(output_file, sep='\t', index=False)

misa_df = misa_df[['ID', 'SSR nr.']]

misa_df['name'] = misa_df['ID'].str.replace(
    ':', '_').str.replace('-', '_') + '_' + misa_df['SSR nr.'].astype(str)

misa_df['chr'] = misa_df['ID'].str.split(':').str.get(0)
misa_df['start'] = misa_df['ID'].str.split(
    ':').str.get(1).str.split('-').str.get(0)
misa_df['end'] = misa_df['ID'].str.split(
    ':').str.get(1).str.split('-').str.get(1)

misa_df = misa_df[['chr', 'start', 'end', 'name']]
# output
misa_df.to_csv(output_bed,
               sep='\t', index=False, header=False)

```
运行一下`python misa_result_2_bed.py Avena_eriantha.fasta.gz_second.tsv Avena_eriantha.fasta.gz_second_vp1c.misa Avena_eriantha.fasta.gz_second_vp1c.bed`
## 1.2 根据bed文件制作fasta文件
```shell
bedtools getfasta -fi ACD_sativa_sfs_genome.fasta -bed /mnt/e/oatdatabase/05_SSR/01_SFS/SFSchr.fasta_second_vp1c.fasta.bed -nameOnly -fo SFSchr.fasta_second_vp1c.fasta
```
## 1.2 利用fasta进行引物设计
因为数据非常庞大，如果利用primer3单进程设计引物肯定是不可行的，我们需要进行改造让primer3支持多进程。于是写了一个多进程的py程序。
```python
# yukaiquan
# email:1962568272@qq.com


import primer3
import pandas as pd
import sys
import argparse
import os
import time
from progressbar import ProgressBar
from multiprocessing import Pool
from functools import partial

# primer_condition


def main():
    start = time.time()
    args = check_options(get_options())
    print(args.input + ' ' + str(args.num) + ' ' + str(args.min_size) + ' ' + str(args.max_size) + ' ' +
          str(args.min_tm) + ' ' + str(args.max_tm) + ' ' + str(args.min_gc) + ' ' + str(args.max_gc) + ' ' + str(args.thread) + ' ' + args.output)
    global_args = gloabal_args(int(args.num), int(args.min_size), int(args.max_size), int(
        args.min_gc), int(args.max_gc), int(args.min_tm), int(args.max_tm))
    # primer3_py 进行引物设计，输出结果文件####primer3模块
    # read fasta
    pbar = ProgressBar()
    pbar.start()
    with open(args.input, 'r') as f:
        lines = f.readlines()
        index, seq = readfasta(lines)
    num_treads = int(args.thread)
    # 分割数据
    n = len(index) // num_treads
    print(len(index))
    indexs = [index[int(i):int(i + n)] for i in range(0, len(index), n)]
    seqs = [seq[int(i):int(i + n)] for i in range(0, len(seq), n)]
    pool = Pool(processes=num_treads)
    primer_result = pd.DataFrame()
    partial_func = partial(
        primer3_design, global_args=global_args)
    # primer_dfs = pd.DataFrame()
    for i, j in zip(indexs, seqs):
        zip_args = list(zip(i, j))
        primer_df = pool.starmap_async(
            partial_func, zip_args)
        # primer_df = pd.concat([primer_dfs, primer_df])
        # print(primer_df)
        # print(primer_df.get())
        for k in primer_df.get():
            primer_result = pd.concat([primer_result, k])
        # primer_result = primer_result.append(k)
    pool.close()
    pool.join()
    # primer_result = primer3_design(index, seq, global_args)
    primer_result.to_csv(args.output, sep="\t")
    end = time.time()
    print("time: " + str(end - start))
    pbar.finish()
    print("finish")


def primer3_design(index, seq, global_args):
    # build table
    index = index.split("\n")[0]
    seq = seq.split("\n")
    # primer_df = pd.DataFrame()
    # primer finder, dic -> datafrme
    # for i in range(len(index)):
    seq_args = {
        'SEQUENCE_ID': str(index),
        'SEQUENCE_TEMPLATE': str(seq),
        'SEQUENCE_INCLUDED_REGION': [0, len(str(seq))-1],
    }
    GeneID = str(index)

    primer3_result = primer3.bindings.designPrimers(seq_args, global_args)
    # print(primer3_result)
    # change dic
    primer3_result_table_dict = {}
    for j in range(primer3_result["PRIMER_PAIR_NUM_RETURNED"]):
        primer_id = str(j)
        for key in primer3_result:
            if primer_id in key:
                # 要将每个信息中的数字和下划线去掉
                info_tag = key.replace("_" + primer_id, "")
                # 就是把不同的引物对结果归到一起
                try:
                    primer3_result_table_dict[info_tag]
                except:
                    primer3_result_table_dict[info_tag] = []
                finally:
                    primer3_result_table_dict[info_tag].append(
                        primer3_result[key])

    df_index = []

    # append dataframe
    for m in range(primer3_result["PRIMER_PAIR_NUM_RETURNED"]):
        df_index.append(GeneID + "_" + str(m + 1))
    primer3_result_df = pd.DataFrame(
        primer3_result_table_dict, index=df_index)
    # primer_df = pd.concat([primer_df, primer3_result_df])
    # primer_df = primer_df.append(primer3_result_df)
    return primer3_result_df


def readfasta(lines):
    seq = []
    index = []
    seqplast = ""
    numlines = 0
    for i in lines:
        if ">" in i:
            index.append(i.replace("\n", "").replace(">", ""))
            seq.append(seqplast.replace("\n", ""))
            seqplast = ""
            numlines += 1
        else:
            seqplast = seqplast + i.replace("\n", "")
            numlines += 1
        if numlines == len(lines):
            seq.append(seqplast.replace("\n", ""))
    seq = seq[1:]
    return index, seq


def gloabal_args(return_num, min_size, max_size, min_gc, max_gc, min_tm, max_tm):
    global_args = {
        'PRIMER_NUM_RETURN': return_num,  # number of primer return
        ### 引物 TM 值设定: ###
        # TM 的计算方法， 1 表示使用 the SantaLucia parameters (Proc Natl Acad Sci 95:1460-65)
        'PRIMER_TM_FORMULA': 1,
        'PRIMER_MIN_TM': min_tm,  # 55
        'PRIMER_OPT_TM': int((max_tm + min_tm)/2),
        'PRIMER_MAX_TM': max_tm,  # 65
        'PRIMER_PAIR_MAX_DIFF_TM': 2.0,         # 两个引物之间的 TM 值最多相差 5 摄氏度
        'PRIMER_WT_TM_LT': 0,
        'PRIMER_WT_TM_GT': 0,
        'PRIMER_PAIR_WT_DIFF_TM': 0.0,
        'PRIMER_MIN_SIZE': min_size,  # 18
        'PRIMER_OPT_SIZE': int((max_size + min_size)/2),
        'PRIMER_MAX_SIZE': max_size,  # 22
        'PRIMER_WT_SIZE_LT': 0.0,
        'PRIMER_WT_SIZE_GT': 0.0,
        # GC含量控制
        'PRIMER_MIN_GC': min_gc,  # minimum GC content 40
        'PRIMER_MAX_GC': max_gc,  # maximum GC content 65
        # 引物热力学计算#
        # 开启热力学计算，开启后，根据热力学的值 TH 值来计算罚分。 TH 的罚分方法为：罚分系数 * (1 / (引物TM - 4 - TH值))。
        # 该方法好处是，TH 值越大，罚分力度越重。
        'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,  # alignment of oligo
        # 引物自身进行反向互补，罚分系数计算为 9 / ( 1/(60-4-45) - 1/(60-4-0) ) = 123.2
        # 这样计算罚分系数的原则是，若 TM 为最佳 60 摄氏度的时候，所允许的最大罚分值减去最小罚分值为 9，和 3' 端碱基的稳定性的罚分额度一致。
        'PRIMER_MAX_SELF_ANY': 8,
        'PRIMER_WT_SELF_ANY': 0.0,
        'PRIMER_MAX_SELF_ANY_TH': 45.0,
        'PRIMER_WT_SELF_ANY_TH': 123.2,
        # 引物自身进行 3' 端反向互补形成引物二聚体，罚分系数计算为 9 / ( 1/(60-4-35) - 1/(60-4-0) ) = 302.4
        'PRIMER_MAX_SELF_END': 3,
        'PRIMER_WT_SELF_END': 0.0,
        'PRIMER_MAX_SELF_END_TH': 35.0,
        'PRIMER_WT_SELF_END_TH': 302.4,
        # left primer 和 right primer 序列的反向互补
        'PRIMER_PAIR_MAX_COMPL_ANY': 8,
        'PRIMER_WT_PAIR_COMPL_ANY': 0.0,
        'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45.0,
        'PRIMER_WT_PAIR_COMPL_ANY_TH': 123.2,
        # left primer 和 right primer 进行 3' 端反向互补形成引物二聚体
        'PRIMER_PAIR_MAX_COMPL_END': 3,
        'PRIMER_WT_PAIR_COMPL_END': 0.0,
        'PRIMER_PAIR_MAX_COMPL_END_TH': 35.0,
        'PRIMER_WT_PAIR_COMPL_END_TH': 302.4,
        # 发夹结构，罚分系数计算为 9 / ( 1/(60-4-24) - 1/(60-4-0) ) = 672
        'PRIMER_MAX_HAIRPIN_TH': 24.0,
        'PRIMER_WT_HAIRPIN_TH': 672.0,
        # 3' 端碱基的稳定性
        'PRIMER_MAX_END_STABILITY': 9,
        'PRIMER_WT_END_STABILITY': 1,
        ### 碱基序列设定： ###
        # 模板序列中包含小写字符不影响引物设计
        'PRIMER_LOWERCASE_MASKING': 0,
        # 引物序列中不能包含单核苷酸连续长度超过 4 bp
        'PRIMER_MAX_POLY_X': 4,
        # 引物中允许的 N 的数目
        'PRIMER_MAX_NS_ACCEPTED': 0,
        # 每个 N 的罚分
        'PRIMER_WT_NUM_NS': 0.0,
        # 引物中 3' 端 5bp 碱基中允许的最大 Gs 或 Cs 的数目
        'PRIMER_MAX_END_GC': 4,
        # 引物中 3' 端碱基中不能出现连续的 Gs 和 Cs 序列
        'PRIMER_GC_CLAMP': 0,
        # 是否允许有诸如 N A R Y 等类型的碱基。必须在设定PRIMER_MAX_NS_ACCEPTED 不为 0 后方有效。
        'PRIMER_LIBERAL_BASE': 1,
        # 如果设置为 1，则 C 能与 S 完美匹配，任意碱基能和 N 完美匹配。
        'PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS': 0,
        ### 碱基质量设定： ###
        # primer 序列允许最小的碱基质量
        'PRIMER_MIN_QUALITY': 0,
        'PRIMER_MIN_END_QUALITY': 0,
        'PRIMER_QUALITY_RANGE_MIN': 0,
        'PRIMER_QUALITY_RANGE_MAX': 100,
        'PRIMER_WT_SEQ_QUAL': 0.0,
        'PRIMER_WT_END_QUAL': 0.0,
        ### PCR 反应体系设定: ###
        # 单价盐离子浓度(mM)
        'PRIMER_SALT_MONOVALENT': 50.0,
        # PRIMER_SALT_CORRECTIONS=1 means use the salt correction in SantaLucia et al 1998
        'PRIMER_SALT_CORRECTIONS': 1,
        # 二价镁离子浓度(mM)
        'PRIMER_SALT_DIVALENT': 1.5,
        # 总dNTP浓度(mM)
        'PRIMER_DNTP_CONC': 0.6,
        # DNA产物浓度(mM)
        'PRIMER_DNA_CONC': 50.0,
        ### PCR 产物的设定： ###
        'PRIMER_PRODUCT_MIN_TM': -1000000.0,
        'PRIMER_PRODUCT_OPT_TM': 0.0,
        'PRIMER_PRODUCT_MAX_TM': 1000000.0,
        'PRIMER_PAIR_WT_PRODUCT_TM_LT': 0.0,
        'PRIMER_PAIR_WT_PRODUCT_TM_GT': 0.0,
        'PRIMER_PRODUCT_OPT_SIZE': 0,
        'PRIMER_PAIR_WT_PRODUCT_SIZE_LT': 0.0,
        'PRIMER_PAIR_WT_PRODUCT_SIZE_GT': 0.0,
        # 罚分因子
        # left primer和right primer的罚分之和。此和乘以此系数，再加其它罚分作为最终罚分。
        'PRIMER_PAIR_WT_PR_PENALTY': 1.0,
        # 将 internal oligo 的罚分乘以此系数，加入到引物的最终罚分中。
        'PRIMER_PAIR_WT_IO_PENALTY': 0.0,
        ### internal oligo 的设定：###
        # TM 值
        'PRIMER_INTERNAL_MIN_TM': 57.0,
        'PRIMER_INTERNAL_OPT_TM': 60.0,
        'PRIMER_INTERNAL_MAX_TM': 63.0,
        'PRIMER_INTERNAL_WT_TM_LT': 1.0,
        'PRIMER_INTERNAL_WT_TM_GT': 1.0,
        # 长度
        'PRIMER_INTERNAL_MIN_SIZE': 18,
        'PRIMER_INTERNAL_OPT_SIZE': 20,
        'PRIMER_INTERNAL_MAX_SIZE': 27,
        'PRIMER_INTERNAL_WT_SIZE_LT': 1.0,
        'PRIMER_INTERNAL_WT_SIZE_GT': 1.0,
        # GC 含量
        'PRIMER_INTERNAL_MIN_GC': 20.0,
        'PRIMER_INTERNAL_MAX_GC': 80.0,
        'PRIMER_INTERNAL_OPT_GC_PERCENT': 50.0,
        'PRIMER_INTERNAL_WT_GC_PERCENT_LT': 0.0,
        'PRIMER_INTERNAL_WT_GC_PERCENT_GT': 0.0,
        # 热力学
        'PRIMER_INTERNAL_MAX_SELF_ANY': 12.00,
        'PRIMER_INTERNAL_WT_SELF_ANY': 0.0,
        'PRIMER_INTERNAL_MAX_SELF_END': 12.00,
        'PRIMER_INTERNAL_WT_SELF_END': 0.0,
        # 碱基序列
        'PRIMER_INTERNAL_MAX_POLY_X': 5,
        'PRIMER_INTERNAL_MAX_NS_ACCEPTED': 0,
        # 碱基质量
        'PRIMER_INTERNAL_MIN_QUALITY': 0,
        'PRIMER_INTERNAL_WT_END_QUAL': 0.0,
        'PRIMER_INTERNAL_WT_SEQ_QUAL': 0.0,
        # 非引物数据库
        # PRIMER_INTERNAL_MAX_LIBRARY_MISHYB=12.00
        # PRIMER_INTERNAL_WT_LIBRARY_MISHYB=0.0
        # PCR 反应体系设定: ###
        # 单价盐离子浓度(mM)
        'PRIMER_INTERNAL_SALT_MONOVALENT': 50.0,
        'PRIMER_INTERNAL_SALT_DIVALENT': 1.5,
        # 总dNTP浓度(mM)
        'PRIMER_INTERNAL_DNTP_CONC': 0.0,
        # DNA产物浓度(mM)
        'PRIMER_INTERNAL_DNA_CONC': 50.0,
        ### PCR 产物的设定： ###
        # 'PRIMER_PRODUCT_SIZE_RANGE': [140, 160],

    }
    return global_args


def check_options(parser):
    args = parser.parse_args()
    if not os.path.exists(args.input):
        parser.print_help()
        sys.exit(1)
    if not args.output:
        parser.print_help()
        sys.exit(1)
    return args


def get_options():
    parser = argparse.ArgumentParser(description="Primer3")
    parser.add_argument(
        "-i", "--input", help="fasta input file (seq name must uniq)", required=True)
    parser.add_argument(
        "-n", "--num", help="number of primers to be designed", default=1
    )
    parser.add_argument(
        "-min_size", "--min_size", help="minimal size of primers", default=18
    )
    parser.add_argument(
        "-max_size", "--max_size", help="maximal size of primers", default=22
    )
    parser.add_argument(
        "-min_tm", "--min_tm", help="minimal melting temperature of primers", default=55
    )
    parser.add_argument(
        "-max_tm", "--max_tm", help="maximal melting temperature of primers", default=65
    )
    parser.add_argument(
        "-min_gc", "--min_gc", help="minimal GC content of primers", default=40
    )
    parser.add_argument(
        "-max_gc", "--max_gc", help="maximal GC content of primers", default=70
    )
    parser.add_argument(
        "-t", "--thread", help="number of threads", default=1
    )
    parser.add_argument("-o", "--output", help="output file", required=True)
    return parser


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nCTRL-C detected! Exiting...")
        sys.exit(0)

```
那么运行一下吧----
`python primer3_design.py -i SFSchr.fasta_second.fasta_vp1c_uniq.fasta -t 16 -o SFSchr.fasta_second.fasta_vp1c_uniq.out`

