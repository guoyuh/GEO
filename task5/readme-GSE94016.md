# 分析分析 GSE94016 数据集 来做练习

首先花费3个小时阅读原文查看作者收据设计及分析流程：

### Whole-genome expression profile

To assess whole-genome expression, we individually extracted total RNAs from liver tumours of orthotopic xenograft mice using TRIZOL Reagent (Life technologies) and checked for a RIN number to inspect RNA integrity by an Agilent Bioanalyzer 2100 (Agilent technologies). The qualified total RNAs were further purified by RNeasy micro kit (QIAGEN) and RNase-Free DNase Set (QIAGEN). Then, total RNAs were amplified, labelled and purified by using GeneChip 3′IVT Express Kit to obtain biotin labelled cRNA (Affymetrix). Array hybridisation and wash were performed with constant rotation on the PrimeView Human Gene Expression Assay (Affymetrix). Microarrays were scanned by GeneChip Scanner 3000 (Affymetrix) and Command Console Software 3.1 (Affymetrix) with default settings. 

Raw data were normalized by robust multiarray analysis **(RMA)** algorithm, Gene Spring Software 11.0 (Agilent technologies). 

All microarray data were deposited in the Gene Expression Omnibus database (<http://www.ncbi.nlm.nih.gov/geo/>) Accession Number GSE94016.

To identify DEGs, we compared gene expression intensities among samples at two different time points using **Welch’s *t*-test with two-tailed *P*-value < 0.05**, which was adjusted by FDR for multiple testing.

通常我们默认作者对其芯片处理是准确的，所以直接下载其表达矩阵即可。

作者鉴定到的差异基因数量有点多，是13,247 differentially expressed genes (DEGs) by multiple comparisons with false discovery rate (FDR) adjustment (*P* < 0.05; Supplementary Data [1](https://www.nature.com/articles/s41467-018-03024-2#MOESM4)).  

聚类分析发现4个时间点的样本(W2, W3, W4, W5)可以粗浅的分成2组：

- one primarily including all samples at the second week after orthotopic implantation, two samples at the third week and one sample at the fifth week; 
- the other including all samples at the fourth and fifth weeks (except W5-a) and three samples at the third week