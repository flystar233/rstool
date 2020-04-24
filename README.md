# rstool使用说明
## 功能介绍
rstool 是使用 python 编写的重测序流程软件,当前版本为 v1.6。主要功能如下：	
1. index	（参考基因组建立索引）
2. data_filter	(下机数据过滤)
3. bwa	（bwa比对）
4. realign	（比对过滤，GATK 局部重新比对）
5. snpcall	(GATK + freebayes + bcftools 变异检测)
8. cnvcall (cnvcaller 变异检测)
7. vcf_marge	（3种变异检测结果整合）

## 功能内容
1. 使用 cutadapt 和 sickle 进行接头处理和数据质控	
**cutadapt** 是处理数据接头最常使用的软件。**sickle** 是数据质控中较为准确的软件，内部使用滑动窗口式过滤方法，适应于如今测序数据。
2. 使用 bcftools 和 freebayes 同时进行 snp calling	
除了 **GATK** 之外，最常用的变异检测软件就是这两个，我们使用 2 个软件同时 snp calling 再取 snp 位点交集，之后还可以与 **GATK** 软件结果取 snp位点交集。这样处理的结果很大可能会一刀切减少 snp 位点数据，可以根据实际情况和测序方法的差别，改变过滤条件。位点数据应在十万到百万级别。
3. vcf 文件的聚合
使用 bcftools 取位点交集。
4. **cnvcaller** 是研究群体 cnv 较为方便的软件。使用去 PCR 重复的 bam 文件进行 cnv 检测。

## 特点
1. 详细的命令说明
2. 子命令方法，各个命令间基本独立
3. 自动识别染色体号
4. 生成流程脚本
5. python源码加二进制文件（无需python和流程所需扩展包）
6. 可扩展的功能

## 使用方法
一：数据过滤
```bash
./rstool_v1.6 data_filter -l raw_data.list -a CACTCGACTAGCATCA
```
参数解释：	
-l 下机数据文件表，raw_data.list 格式如下：
```
testa	/zfssz3/NASCT_BACKUP/MS_PMO2017/testa_1.fq.gz	/zfssz3/NASCT_BACKUP/MS_PMO2017/testa_2.fq.gz
testb	/zfssz3/NASCT_BACKUP/MS_PMO2017/testb_1.fq.gz	/zfssz3/NASCT_BACKUP/MS_PMO2017/testb_2.fq.gz
testc	/zfssz3/NASCT_BACKUP/MS_PMO2017/testc_1.fq.gz	/zfssz3/NASCT_BACKUP/MS_PMO2017/testc_2.fq.gz
```	
-a 3'接头序列.
当不提供接头序列时，则不使用 cutadapt 去除接头。

二：建立参考基因组索引	
```bash
rstool_v1.6 index -r genome.fa
```
参数解释：	
-r 参考基因组文件	

三：bwa 比对
```bash
rstool_v1.6 bwa -r 02.Index/genome.fa -l 01.Data_filter/clean_list.txt
```
参数解释：	
-r 在建立参考基因组索引时会默认建立源参考基因组的软链，在比对时需要已经建立好的文件，所以 -r 请使用软链文件。	
-l 数据过滤中默认生成clean data 的文件列表，请使用过滤数据的文件列表。	

如果想加快比对速度,建议使用bwa-mem2.可在软件列表（input.toml）进行替换路径。

四：局部重新比对
```bash
rstool_v1.6 realign -d 03.Bwa/chrbam -rd genome_cut
```
参数解释：	
-d bwa 比对后，分割染色体比对文件的路径	
-rd 参考基因组分割染色体路径	

五：SNP 变异检测
```bash
rstool_v1.6 snpcall -r 02.Index/genome_cut -l 04.Realign -i F
```
参数解释：	
-r 参考基因组分割染色体后路径（一般变异结果要求按染色体分割）	
-l 所有染色体 bam文件列表所在路径
-d 单染色体变异检测主目录	
-i 使用 beagle 进行基因组填补	

六：变异整合
```bash
rstool_v1.6 vcf_marge -a freebayes.vcf -b bcftools.vcf
```
参数解释：	
-a 由 freebayes 检测的变异文件	
-b 由 bcftools 检测的变异文件

七：CNV 变异检测
```bash
rstool_v1.6 snpcall -r 02.Index/genome.fa  -l 04.Realign/dup.bam.list
```
参数解释：	
-r 参考基因组文件
-l 经过去 PCR 重复的 bam 文件列表。（cnvcaller 不推荐对 bam 文件进行过滤）

## NOTE：
1. 流程使用中，各个文件、目录尽量使用绝对路径。
2. 现在 GATK 变异检测可以在一个自动处理染色体问题，只需输入参考基因组按染色体分割文件路径即可。
3. 在变异整合时，同时只能处理 2 个文件，如果还想整合 GATK 的结果，按脚本替换成 GATK 文件即可。
4. 变异过滤条件较为严格，建议样本数大于 100 使用此流程。
5. 流程软件使用集群稳定版本，如需替换可在源码 software path 下替换。
