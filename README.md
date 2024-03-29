# *DRaWR* - Discriminative Random Walk with Restart
Charles Blatti [blatti@illinois.edu] and Saurabh Sinha  
KnowEnG BD2K Center of Excellence  
University of Illinois Urbana-Champaign  

## Table of Contents
1. [Motivation](#motivation)
2. [Installation](#installation)
3. [Tutorial](#tutorial)
    1. [Creating Gene Sets](#creating-gene-sets)
    2. [Building Heterogeneous Networks](#building-heterogeneous-networks)
    3. [Example DRaWR Runs](#example-drawr-runs)
4. [DRaWR Resources](#drawr-resources)
    1. [Function Parameters](#function-parameters)
    2. [Input File Formats](#input-file-formats)
    3. [Output File Formats](#output-file-formats)

## Motivation

Analysis of co-expressed gene sets typically involves testing for enrichment of different “properties” such as biological processes, pathways, transcription factor binding, etc., one property at a time. This approach ignores any known relationships among the properties or genes themselves. Previous work has sought to exploit these relationships by building biological networks that combine multiple types of gene-gene or gene-property relationships, and performing network analysis to identify other genes and properties most relevant to a given gene set. However, these existing network-based method often collapse information about individual property annotations to create simplified, homogeneous networks.

We present DRaWR, a network-based method for ranking genes or properties related to a given gene set. Such related genes or properties are identified from among the nodes of a large, heterogeneous network of biological information. Our method involves a random walk with restarts, performed on an initial network with multiple node and edge types, preserving more of the original, specific property information than current methods that operate on homogeneous networks. In this first stage of our algorithm, we find the properties that are the most relevant to the given gene set and extract a subnetwork of the original network, comprising only the relevant properties. We then rerank genes by their similarity to the given gene set, based on a second random walk with restarts, performed on the above subnetwork.

![Method Overview](images/DRaWR_method.small.png)

[Return to TOC](#table-of-contents)

## Installation

### Local copy of DRaWR repository

If you wish to use the sample files necessary to complete this tutorial or the datasets from the paper, first clone this repository from github:
```
git clone https://github.com/cblatti3/DRaWR.git
```

To use the 5 species networks from the paper, you must first unzip their contents:
```
gunzip networks/5ins_cdhmw.names.edge.gz
gunzip networks/5sp_adhiw.names.edge.gz
```

### DRaWR R Package Installation

The [DRaWR package](https://cran.r-project.org/web/packages/DRaWR/index.html) has been uploaded to the Comprehensive R Archive Network [CRAN](https://cran.r-project.org/index.html).  To install DRaWR, enter the R environment:
```
R
```

and use the standard R package installation method:
```
install.packages("DRaWR")
```
You will be prompted to select a local CRAN mirror to install from.  The DRaWR package and its dependencies will be automatically downloaded and installed.  Details on alternative methods for installing CRAN packages are available [here](https://cran.r-project.org/doc/manuals/r-release/R-admin.html). 

Every time you open R and wish to use the installed DRaWR library, load its functions with
```
library(DRaWR)
```

[Return to TOC](#table-of-contents)

## Tutorial

This section of the README is meant to walk a user through a process of using DRaWR to find related genes and features/properties that relate to a gene set (or sets) of interest.  In these examples, we will examine the sets of genes that are expressed in 12 specific spatial-temporal domains of the developing *Drosophila* embryo.   

### Creating Gene Sets

The first step is to create a file for each gene set that lists the members of that set.  This query [gene set file](#gene-set-file) format should list one gene name on each row. For example, the [Fly Brain Primordium Gene Set Query File](gene_sets/dmel/5268_brain_primordium.names.txt) contains:
```
FBgn0008636
FBgn0040534
FBgn0040918
FBgn0041105
...
```

If the genes have distinct weights (e.g. different expression fold changes), a second column (separated by a tab) can be added.  Node weights must be positive with larger values meaning stronger evidence.  For example:
```
FBgn0041156 5.2
FBgn0041186 1.3
FBgn0010105 2.6
```

We also must create a file that specifies the universe of genes we wish to rank.  Our query gene sets should all be a subset of the universe set.  Typically, it will contain all of the genes on the species of interest. This is done with the same [gene set file](#gene-set-file) format above.  
For example, the *Drosophila* universe file, [gene_sets/universes/dmel.ids.uni.txt](gene_sets/universes/dmel.ids.uni.txt) starts:
```
FBgn0000032 1
FBgn0000289 1
FBgn0008636 1
FBgn0040505 1
...
```

Finally, we need a way to refer to all 12 query gene set files we wish to examine with DRaWR.  For this, we create a [gene set list](#gene-set-list).  This file must contain a line for each gene set with the location of that gene set file relative to the gene set list.  In the DRaWR repository, the 12 *Drosophila* gene set list is located at [setlists/dmel_possets12.list.txt](setlists/dmel_possets12.list.txt), so it contents must list the location of the gene set files relative to the setlists directory:
```
../gene_sets/dmel/2270_procephalic_ectoderm_anlage_in_statu_nascendi.names.txt
../gene_sets/dmel/4227_ventral_ectoderm_primordium.names.txt
../gene_sets/dmel/5155_ventral_epidermis_primordium.names.txt
../gene_sets/dmel/6069_embryonic_dorsal_epidermis.names.txt
```

### Building Heterogeneous Networks

We must also assemble a heterogeneous network before starting our DRaWR analysis.  This network will contain both gene nodes and property (feature) nodes as well as many different types of relationships connecting them.  To represent the network, we will use a [network edge file](#network-edge-file) format.  Each line in this file will represent an edge in our network.  The network edge file is expected to be a tab separated file with 4 columns [name of node 1, name of node 2, weight of edge/relationship, type of edge/relationship]. The edge weights must be positive with larger values meaning stronger relationships. 

For example, from the heterogeneous *Drosophila* network, [networks/dmel_cdhmw.names.edge](networks/dmel_cdhmw.names.edge), there are weighted gene-gene (GG) edges representing protein sequence similarity.  The edge type is 'homol' and the larger weights indicate greater similarity.  
```
...
FBgn0040765 FBgn0000289 6.14    homol
FBgn0041105 FBgn0000289 2.32    homol
FBgn0042205 FBgn0000289 8.2     homol
FBgn0043364 FBgn0000289 7.65    homol
...
```

To represent edges between a property node and a gene node, the property node must be listed as node1. The second stage of DRaWR will select the most relevant property nodes based on the query genes. In the heterogeneous *Drosophila* network, [networks/dmel_cdhmw.names.edge](networks/dmel_cdhmw.names.edge), there are three types of property-gene (PG) edges, each connecting multiple property nodes to their related genes.
```
...
mt_fkh_u5_gc        FBgn0040717 3.17    motif_u5
mt_grh_u5_gc        FBgn0040717 2.67    motif_u5
...
Pdom_Rad17          FBgn0032244 5.26    pfam_domain
Pdom_Rep_fac_C      FBgn0032244 8.2     pfam_domain
...
chip_CAD_Bchip_s5   FBgn0033062 2.59    chip_binding
chip_EVE_Mseq_s14   FBgn0033062 2.66    chip_binding
...
```

DRaWR does not currently support edges between the property nodes.  We would recommend to limit the number of gene-gene (GG) edge types in order to get the greatest improvements with the second stage of DRaWR.  Although edges from the network edge file are initially directional, by default, the network will be converted to an undirected, normalized by edge type adjacency matrix.  The maximum edge weights will be preserved if there are redundant edges.

### Example DRaWR Runs

Now that we have our gene sets and our network prepared, we are ready to run DRaWR.  Once in R, we must first load the DRaWR functions:
```
library(DRaWR)
```

To run our 12 gene set queries with the heterogeneous *Drosophila* network with the default settings, the command requires that we specify the setlist ("possetfile"), the gene universe ("unifile"), the location of the heterogeneous network ("networkfile"), the property-gene edge types ("property_types"), and a location for the outputs ("outdir").  We also will decide how many property nodes we wish to keep for each PG edge type ("st2keep") and the restart probability that we will use for our random walks ("restarts").  More details are available below for all of the DRaWR [function parameters](#function-parameters).
```
DRaWR(possetfile = "setlists/dmel_possets12.list.txt",
    unifile = "gene_sets/universes/dmel.ids.uni.txt",
    networkfile = "networks/dmel_cdhmw.names.edge",
    property_types = c("motif_u5", "pfam_domain", "chip_binding"),
    outdir = "sample_results/",
    st2keep = 50, 
    restarts = c(0.3))
```

First, DRaWR will make the network undirected, normalized, and remove redundant edges.  Then for each of our 12 query gene sets, it will read in the query set, run the 'baseline', 'stage 1', and 'stage 2' random walks with restart (RWR).  For each random walk, it will calculate an Area Under the Receiver Operating Characteristics Curve (AUROC) using left out genes from the gene set as well as record the number of iterations it took for the RWR to converge.  

DRaWR will by default produce two files in output directory.  The name of these files will contain a concatenation of the arguments used to produce the run.  The first file is the converged probability distribution from the 'baseline' RWR.  This [prediction file](#prediction-file) ends with the suffix '.base'.  It contains a row for each node in the network and its four columns indicate the node name, the node type, the node's presence in the gene universe, and the probability of a 'baseline' walker being at that node.  

For example, we have the output from above, [sample_results/dmel.ids.uni.dmel_cdhmw.names.undir.weight.type.50.1e-04.0.3.base](sample_results/dmel.ids.uni.dmel_cdhmw.names.undir.weight.type.50.1e-04.0.3.base):
```
node            type    universe    baseline
FBgn0000008     -1      1           5.63872111193178e-05
FBgn0000014     -1      1           8.07988313120129e-05
FBgn0000015     -1      1           5.93092943412521e-05
FBgn0000017     -1      1           0.000112563218659276
FBgn0000018     -1      1           3.0293883832415e-05
FBgn0000022     -1      1           5.32024158248407e-05
...
```

The second output of DRaWR will be the performance of each RWR stage in ranking the left out genes from the query gene set.  This [prediction file](#prediction-file) ending in '.stats' contains summary information from each RWR stage, fold, RWR restart parameter, and query gene set completed by the DRaWR() call as well as the size of gene set, number of iterations, and AUROC metrics. 

From the example above, we output, [sample_results/dmel.ids.uni.dmel_cdhmw.names.undir.weight.type.50.1e-04.50.1.0.3.dmel_possets12.list.stats](sample_results/dmel.ids.uni.dmel_cdhmw.names.undir.weight.type.50.1e-04.50.1.0.3.dmel_possets12.list.stats):
```
network             direct  weight  normalize   uni             restart maxiters    thresh  st2keep posset                                                      nfolds  iter    stage       aucval  rwr_iters   ntrain
dmel_cdhmw.names    undir   weight  type        dmel.ids.uni    0.3     50          1e-04   50      2270_procephalic_ectoderm_anlage_in_statu_nascendi.names    1       1       baseline    0.667   0           13609
dmel_cdhmw.names    undir   weight  type        dmel.ids.uni    0.3     50          1e-04   50      2270_procephalic_ectoderm_anlage_in_statu_nascendi.names    1       1       stage1      1       17          222
dmel_cdhmw.names    undir   weight  type        dmel.ids.uni    0.3     50          1e-04   50      2270_procephalic_ectoderm_anlage_in_statu_nascendi.names    1       1       diff        1       0           222
dmel_cdhmw.names    undir   weight  type        dmel.ids.uni    0.3     50          1e-04   50      2270_procephalic_ectoderm_anlage_in_statu_nascendi.names    1       1       stage2      1       22          222
...
```

The default behavior of the DRaWR method is to do no cross-validation.  The performance values reported above are the 'training' AUROCs.  For a better estimate of performance, DRaWR should be run with cross validation.  To do this, the user just needs to specify the number of folds ('nfolds'). For example, 10-fold cross-validation:
```
DRaWR(possetfile = "setlists/dmel_possets12.list.txt",
    unifile = "gene_sets/universes/dmel.ids.uni.txt",
    networkfile = "networks/dmel_cdhmw.names.edge",
    property_types = c("motif_u5", "pfam_domain", "chip_binding"),
    outdir = "sample_results/",
    st2keep = 50, 
    restarts = c(0.3),
    nfolds = 10)
```

The default behavior of DRaWR is to not print out the converged probability distributions for the RWRs of each stage.  However, these values are necessary if one wants to return the top ranked genes, or the most related property nodes.  To override the default behavior and produce the [prediction file](#prediction-file), you need to specify the 'writepreds' option:

```
DRaWR(possetfile = "setlists/dmel_possets12.list.txt",
    unifile = "gene_sets/universes/dmel.ids.uni.txt",
    networkfile = "networks/dmel_cdhmw.names.edge",
    property_types = c("motif_u5", "pfam_domain", "chip_binding"),
    outdir = "sample_results/",
    st2keep = 50, 
    restarts = c(0.3),
    writepreds = 1)
```

This produces a file ending in '.rwr' for each query gene set containing information about each node in the network after running DRaWR.  The 'baseline', 'stage1', and 'stage2' columns contain the converged RWR probabilities from each stage of DRaWR.  

Using our example above, we produce for the brain primordium query set, [sample_results/dmel.ids.uni.dmel_cdhmw.names.undir.weight.type.50.1e-04.50.1.0.3.5268_brain_primordium.names.1.rwr](sample_results/dmel.ids.uni.dmel_cdhmw.names.undir.weight.type.50.1e-04.50.1.0.3.5268_brain_primordium.names.1.rwr)
```
node                    type            universe    baseline                train   test    stage1                  diff                    keep    stage2
FBgn0000008             -1              1           5.63872111193178e-05    0       0       2.9902923974918e-05     -2.64842871443998e-05   -1      9.09485166256704e-06
FBgn0000014             -1              1           8.07988313120129e-05    0       0       0.000118572506546682    3.77736752346694e-05    -1      0.000167470108206141
FBgn0000015             -1              1           5.93092943412521e-05    0       0       9.72564314976961e-05    3.7947137156444e-05     -1      0.000141947317124086
... 
chip_CAD_Bseq_s5        chip_binding    0           0.000714960664334332    0       0       0.00109212186511493     0.000377161200780602    1       0.00163031185472102
chip_CAD_Mchip_s5_9     chip_binding    0           0.000593923399913552    0       0       0.000847676120446827    0.000253752720533275    1       0.00136124484478671
chip_CHINMO_Mchip_s5_14 chip_binding    0           0.000380019396277925    0       0       0.000300343904418371    -7.96754918595537e-05   0       0
...
```

To find genes/nodes most related to the query gene set, we can sorting this table by 'stage2' in descending order.  We expect that that gene nodes in the query gene set 'train'=1 or 'test'=1 to show up near the top.  To find the ranking of the important property nodes at the end of stage 1, we can sort the 'diff' column in descending order.  The 'keep' column indicates which property nodes were the most related and were selected for the subnetwork used in stage two. 

There are several other [function parameters](#function-parameters) that can be modified to run DRaWR and more details for these parameters are found below.

The repository contains all 92 *Drosophila* embryonic developement gene expression sets from the paper, as well as the social species aggression related gene sets.  They are available in the gene_sets and the setlists folders.  Multi-species insect and aggression heterogeneous networks are also provided in the networks directory.  These data files are available to the user to recreate the results reported in the paper.

[Return to TOC](#table-of-contents)

## DRaWR Resources

### Function Parameters

#### Relating to input files
| Parameter      | Type    | Default  | Description                                                                                                                                                                                                                                                                             |
|----------------|---------|----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| possetfile     | string  | required | location of file containing relative location of gene sets to test.                                                                                                                                                                                                                     |
| unifile        | string  | required | location of file listing gene universe.                                                                                                                                                                                                                                                 |
| networkfile    | string  | required | location of file containing network contents.                                                                                                                                                                                                                                           |

#### Relating to output files
| Parameter      | Type    | Default  | Description                                                                                                                                                                                                                                                                             |
|----------------|---------|----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| outdir         | string  | required | prefix of location of file to write performance results (optionally prediction results).                                                                                                                                                                                                |
| writepreds     | boolean | FALSE    | write predictions out to a file.                                                                                                                                                                                                                                                        |
| nfolds         | int     | 1        | number of folds for cross validation, Default is 1, no cross-validation.                                                                                                                                                                                                                |

#### Relating to network processing
| Parameter      | Type    | Default  | Description                                                                                                                                                                                                                                                                               |
|----------------|---------|----------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| undirected     | bool    | TRUE     | boolean to make network undirected.                                                                                                                                                                                                                                                       |
| unweighted     | bool    | FALSE    | boolean to make network unweighted.                                                                                                                                                                                                                                                       |
| normalize      | string  | "type"   | "type" or "none" normalization method.                                                                                                                                                                                                                                                    |
| property_types | vector  | required | vector containing *ALL* names of property-gene (PG) edge types. May contain PG edge types that are not in the given network, but should not contain gene-gene (GG) edge types. Default is c("allen_brain_atlas", "chip_binding", "gene_ontology", "motif_u5", "pfam_domain", "T1", "T2"). |
| st2keep        | int     | 1        | number of property nodes to keep in second stage for each property type. To skip the second stage of DRaWR just set st2keep to 0.                                                                                                                                                         |

#### Relating to random walks
| Parameter      | Type    | Default  | Description                                                                                                                                                                                                                                                                             |
|----------------|---------|----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| restarts       | vector  | c(0.7)   | vector of restart values to test. Must be between 0 and 1.                                                                                                                                                                                                                              |
| maxiters       | int     | 50       | maximum number of allowable iterations.                                                                                                                                                                                                                                                 |
| thresh         | float   | 0.0001   | threshold for L1 norm convergence.                                                                                                                                                                                                                                                      |


### Input File Formats
#### Network Edge File
The network edge file is expected to be a 4 column, tab separated file. Each row represents an edge in the network and contains the values for:

| Columns       | Type      | Description                                                                                                       |
|---------------|-----------|-------------------------------------------------------------------------------------------------------------------|
| node1         | string    | String name of node1.  If PG edge, must be the property node name.                                                |
| node2         | string    | String name of node2, always a gene node.                                                                         |
| edge_weight   | float     | Positive value for weight of edge where larger numbers mean stronger relationships.                               |
| edge_type     | string    | String name for the type of the edge. PG edge types should be listed in the 'property_types' argument of DRaWR(). |

The order of the nodes typically does not matter for gene-gene edges. By default the network edge file will be converted to an undirected, normalized by edge type adjacency matrix with the maximum edge weight taken in the case of redundancies.

Example network edge file:
```
G1  G6  0.76    typeGG
G2  G1  0.41    typeGG
G2  G7  0.73    typeGG
G4  G5  0.89    typeGG
P1  G2  0.57    typePG.1
P1  G6  0.30    typePG.1
P2  G1  0.07    typePG.1
P2  G4  0.89    typePG.1
P3  G3  0.66    typePG.1
P3  G6  0.80    typePG.1
P4  G4  0.24    typePG.2
P4  G1  0.74    typePG.2
P5  G5  0.20    typePG.2
P5  G2  0.95    typePG.2
P5  G7  0.92    typePG.2
```

#### Gene Set File
The query gene set files or gene universe files should have one node name listed on each row.  If the gene nodes have distinct weights, a second column (separated by a tab) can be added.  Node weights must be positive with larger values meaning stronger evidence.

| Columns   | Type      | Description                                                                                               |
|-----------|-----------|-----------------------------------------------------------------------------------------------------------|
| node      | string    | String name of gene.  Must match network file.                                                            |
| weight    | string    | Optional. Positive value of weight of node where larger numbers mean stronger evidence of gene in set.    |

Example gene set file:
```
G3  2.5
G4  1.3
G5  4.0
```

The gene universe file must contain the appropriate universe of genes contain all genes in the query gene sets. Typically, this will be all of the genes on the species of interest.

#### Gene Set List
This file facilitates the examination of multiple gene sets with the identical settings of DRaWR.  The file must contain the location of a gene set file relative to itself on each line

Example gene set list:
```
../gene_sets/dmel/2270_procephalic_ectoderm_anlage_in_statu_nascendi.names.txt
../gene_sets/dmel/4227_ventral_ectoderm_primordium.names.txt
../gene_sets/dmel/5155_ventral_epidermis_primordium.names.txt
../gene_sets/dmel/6069_embryonic_dorsal_epidermis.names.txt
```

### Output File Formats
#### Statistics File
This file ending in '.stats' is produced to analyze the performance results. Each row reports the size of the restart set ("ntrain"), iterations for the random walk to converge ("rwr_iters"), and AUROC ("aucval") for the corresponding stage ("stage1") of DRaWR, fold ("iter") of cross validation, RWR restart parameter ("restart") setting, and query gene set ("posset").

| Columns   | Type      | Description                                                                                       |
|-----------|-----------|---------------------------------------------------------------------------------------------------|
| network   | string    | Name of heterogeneous network                                                                     |
| direct    | string    | Whether network was left directed ('direct') or made undirected (''undir')                        |
| weight    | string    | If the network was left weighted ('weight') or force to be unweighted ('unweight')                |
| normalize | string    | If the network was normalized by edge type ('type') or not ('none')                               |
| uni       | string    | Name of the gene universe file                                                                    |
| restart   | float     | Random walk restart parameter                                                                     |
| maxiters  | int       | Random walk maxiters parameter                                                                    |
| thresh    | float     | Random walk convergence threshold parameter                                                        |
| st2keep   | int       | Number of property nodes to keep in second stage for each property-gene (PG) edge type            |
| posset    | string    | Name of query gene set file                                                                       |
| nfolds    | int       | Total number of folds specified in the run                                                        |
| iter      | int       | Number of the fold of the current result                                                          |
| stage     | string    | Stage of DRaWR of the current result. Can be 'baseline', 'stage1', 'diff', or 'stage2'            |
| aucval    | float     | AUROC from ranking the left out genes with the RWR probability distribution of the current stage  |
| rwr_iters | int       | Number of iterations for the RWR to converge                                                      |
| ntrain    | int       | Number of genes in the RWR restart set                                                            |

#### Prediction File
This file ending in '.rwr' contains information about each node in the network after running DRaWR.  This file is only produced if writepreds = TRUE.  There is one row for each node in the network, and the columns are:

| Columns   | Type      | Description                                                                                   |
|-----------|-----------|-----------------------------------------------------------------------------------------------|
| node      | string    | The name of that node                                                                         |
| type      | string    | The type of that node. -1 if node is a gene node                                              |
| universe  | int       | 1 if node in gene universe set. 0 otherwise                                                   |
| baseline  | float     | probability of being in the node from converged baseline RWR                                  |
| train     | int       | 1 if node in the gene query set training set, 0 otherwise                                     |
| test      | int       | 1 if node in the gene query set test set, 0 otherwise                                         |
| stage1    | float     | probability of being in the node from converged stage1 RWR                                    |
| diff      | float     | difference between "stage1" and "baseline"                                                    |
| keep      | int       | 1 if the property node was kept for the stage2 RWR, 0 otherwise. -1 if node is a gene node    |
| stage2    | float     | probability of being in the node from converged stage2 RWR                                    |


The file ending in '.base' contains the first 4 columns of the standard prediction file and is used in the acceleration of the convergence of the RWRs.

[Return to TOC](#table-of-contents)