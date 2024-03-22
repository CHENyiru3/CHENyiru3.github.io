---
title: '[Study markdown][2] GO analysis and GSEA practice'
date: 2024-3-22
permalink: /posts/2024/3/22/OGandGESA
tags:
  - gene
  - biomedical
  - my_study_note
---
GOa and GSEA? How to choose....

Welcome to Chen's learning blog! This is a space where I share my thoughts, experiences, and insights on mostly bioinformatics, many biomedical, some computor science! Whether you're a fellow enthusiast or simply curious, I invite you to join me on this journey as we explore the captivating world of gene, protein, RNA... and codes.

So, sit back, grab a beverage, and let's dive into the fascinating realm of gene ontology and gene set enrichment analysis together!

# Interaction between proteins
https://string-db.org
This website can search and generate interaction graph of the interested proteins, here we can input the proteins we selected:
![](https://i.imgur.com/FYtrC3i.png)

## Result: interaction graph of these proteins
![](https://i.imgur.com/dkCsJcp.png)
- There are tow sections:
	1. node: present proteins, each node represents all the proteins produced by a single, protein-coding gene locus
	2. Edges present protein-protein associations
	![](https://i.imgur.com/dv3qf5H.png)
- Setting can be changed
	- ![](https://i.imgur.com/6vilPzL.png)
	- here you can change some settings to update the interaction graph, because some relationship might be lack of evidences
- Stats
	- ![](https://i.imgur.com/OizfLwn.png)
	- PPI enrichment p-values: evaluate a series of protein/gene in interaction network is connected better/not better than random situation
		- When STRING does comparison to random networks it uses the whole genome as the default statistical background
		- ![](https://i.imgur.com/9SrL8FF.png)
		- But what if the number of calculated interaction edges equals expected number?
		- Here claud said:![](https://i.imgur.com/8AogTTz.png)
- gene ontology
	- ![](https://i.imgur.com/pGI9vcR.png)


# Gene Ontology: Unveiling the Layers of Biological Function

Gene Ontology (GO) is a powerful framework that provides a standardized and comprehensive system for describing the attributes of genes and gene products across different organisms. This ontology encompasses three main domains: molecular function, cellular component, and biological process. Let's explore each of these domains in detail.

GO like a database, having 3 terms:
- Molecular function
- Cellular component
- biologicla process

## Molecular Function

The molecular function domain of Gene Ontology focuses on the molecular-level activities carried out by gene products, such as proteins and RNAs. These activities include binding, catalysis, transport, and signaling, among others. This domain describes the biochemical actions that a gene product performs, independent of the context in which it operates.

For example, the molecular function of the enzyme "DNA polymerase" would be "DNA-directed DNA polymerase activity," which captures its role in catalyzing the synthesis of DNA. Similarly, the molecular function of a transcription factor might be "DNA-binding transcription factor activity," reflecting its ability to bind to DNA and regulate gene expression.

## Cellular Component

The cellular component domain of Gene Ontology describes the locations within a cell where gene products are found and perform their functions. This includes organelles, membranes, protein complexes, and other cellular structures. Understanding the cellular component of a gene product provides insights into the specific context in which it operates.

For instance, the cellular component of a mitochondrial enzyme would be "mitochondrion," indicating its localization within this organelle. Likewise, a cell surface receptor would have the cellular component "plasma membrane," reflecting its positioning on the outer membrane of the cell.

## Biological Process

The biological process domain of Gene Ontology encompasses the broader, higher-level functions that gene products contribute to within an organism. These processes can range from fundamental cellular activities, such as cell division and metabolism, to complex physiological processes, like immune response and nervous system development.

For example, the biological process of "DNA repair" would capture the collective actions of various gene products involved in maintaining the integrity of the genetic material. Similarly, the biological process of "photosynthesis" would describe the integrated set of molecular and cellular activities that enable plants to convert light energy into chemical energy.

By considering these three domains of Gene Ontology, researchers can gain a comprehensive understanding of the roles and functions of genes and gene products, both at the molecular level and within the broader context of cellular and organismal biology. This knowledge is crucial for unraveling the complexities of biological systems and advancing our understanding of the mechanisms underlying health and disease.

![](https://i.imgur.com/Y9KjwT6.png)

![](https://i.imgur.com/IcJrffC.png)

## Do GO analysis
### Hit list
First we need a hit list. This list contains the genes of interest.

>假设我们做了一个基因表达实验,比较了正常细胞和肿瘤细胞的基因表达情况。我们发现有100个基因在肿瘤细胞中表达显著上调(比较差异显著性P值<0.05)。这100个基因就构成了我们的"hit list"

### hypergeometric test
>The hypergeometric test uses the hypergeometric distribution to measure the statistical significance of having drawn a sample consisting of a specific number of k successes (out of n total draws) from a population of size N containing K successes

- 举个例子：设总共有29个人，其中11个吸烟者，18个非吸烟者，现从中随机抽取16个样本（在此实验中对应着肺癌病人），有10个是吸烟者，这样的事件是否显著？这里面就需要用到超几何分布

### GO analysis using hypergeometric test
- searching the GO information about the hit list, and comparing with the background genes 
- In other words, this simply measures what is the probability of certain genes belonging to a specific GO class in your list to appear compared to the probability of picking genes randomly from the background gene set and getting enrichment of that same specific gene ontology class.

### GO enrichment analysis using DAVID
I followed this link: https://zhuanlan.zhihu.com/p/560010219 

![](https://i.imgur.com/Sw1lzwe.png)
Input the list of your proteins/genes, DAVID will give you a analysis result in txt file.
But we need to do some visualization.
```R
# 加载所需的库
library(ggplot2)
library(forestplot)
library(pheatmap)
library(igraph)

data <- read.table("result.txt",sep = "\t",header = TRUE)
head(data)

ggplot(data, aes(x = Term, y = Count)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Term", y = "Count") +
  ggtitle("DAVID Analysis Results")

head(data)


```
The enrichment result is like this:
![](https://i.imgur.com/en8tOMg.png)

# GSEA
- GSEA: gene set enrichment analysis 是一种用于分析基因表达数据的统计方法。对比的是感兴趣的基因和某个特定的已知的gene set之前的关系。
- 分析时候的输入：
	- 已知功能的基因集
	- 表达矩阵，也可以是排序好的列表
		- 给定一个基因表达数据集,我们可以根据在两种状态下的表达差异,为每个基因计算一个相关性排序统计量,并排列形成一个有序基因列表L

假设我们研究肺癌与正常肺组织的基因表达差异。我们有10个肺癌样本和10个正常肺组织样本的基因表达谱数据。对于每一个基因,我们可以计算其在肺癌组和正常组的平均表达值,然后对两组平均值进行统计学检验(如t检验或Wilcoxon秩和检验),得到一个相关性评分,用来衡量该基因在两组之间的差异表达程度。

例如:  
基因A在肺癌组平均表达值为10,在正常组平均表达值为2,p值=0.001  
基因B在肺癌组平均表达值为8,在正常组平均表达值为6,p值=0.05  
基因C在肺癌组平均表达值为3,在正常组平均表达值为3,p值=1

我们可以使用-log10(p值)作为相关性评分的统计量,对所有基因按此值从大到小排序,形成一个有序基因列表L:

L = [基因A(-log10(0.001)=3), 基因B(-log10(0.05)=1.3), 基因C(-log10(1)=0), ...]

这样基因A因为在两组间差异最大,就排在L的最前面。基因B次之,基因C无差异排最后。

通过这个排序过程,我们可以看到与肺癌表型关联最密切的基因排在L的前列。接下来就可以检验我们感兴趣的基因集(如某代谢通路基因)在L中是否显著富集在顶端,如果是,就说明该基因集可能参与了肺癌的发生发展过程。

这就是GSEA分析的第一步,构建有序基因列表L的基本思路。希望这个例子能更清楚地解释相关性排序的过程

- 计算出enrichment score。ES值越大，标识改基因集在两个状态之间的差异越显著从排序列表的第一个基因开始，计算累计统计值，当遇到一个落在gene set里面的基因，则增加统计值，当遇到不在gene set里面的，就降低统计值
	- `正值`ES表示基因集在列表的顶部富集，`负值`ES表示基因集在列表的底部富集