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

那么运行一下吧----
`python primer3_design.py -i SFSchr.fasta_second.fasta_vp1c_uniq.fasta -t 16 -o SFSchr.fasta_second.fasta_vp1c_uniq.out`

