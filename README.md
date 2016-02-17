# *DRaWR* - Discriminative Random Walk with Restart
Charles Blatti [blatti@illinois.edu] and Saurabh Sinha  
KnowEnG BD2K Center of Excellence  
University of Illinois Urbana-Champaign  

## Table of Contents
1. [Motivation](#motivation)
2. [Installation](#installation)
3. [Tutorial](#tutorial)
    1. [Gene Sets](#gene-sets)
    2. [Heterogeneous Network](#heterogeneous-network)
    3. [Example Runs](#example-runs)
4. [Resources](#resources)
    1. [Function Parameters](#function-parameters)
    2. [Input Formats](#input-formats)
    3. [Output Formats](#output-formats)

## Motivation

Analysis of co-expressed gene sets typically involves testing for enrichment of different “properties” such as biological processes, pathways, transcription factor binding, etc., one property at a time. This approach ignores any known relationships among the properties or genes themselves. Previous work has sought to exploit these relationships by building biological networks that combine multiple types of gene-gene or gene-property relationships, and performing network analysis to identify other genes and properties most relevant to a given gene set. However, these existing network-based method often collapse information about individual property annotations to create simplified, homogeneous networks.

We present DRaWR, a network-based method for ranking genes or properties related to a given gene set. Such related genes or properties are identified from among the nodes of a large, heterogeneous network of biological information. Our method involves a random walk with restarts, performed on an initial network with multiple node and edge types, preserving more of the original, specific property information than current methods that operate on homogeneous networks. In this first stage of our algorithm, we find the properties that are the most relevant to the given gene set and extract a subnetwork of the original network, comprising only the relevant properties. We then rerank genes by their similarity to the given gene set, based on a second random walk with restarts, performed on the above subnetwork.

![Method Overview](images/DRaWR_method.small.png)

[Return to TOC](#table-of-contents)

## Installation

Clone this repository
```
git clone https://github.com/cblatti3/DRaWR.git
```

Enter the repository and install packages into a local library
```
cd ./DRaWR/
mkdir library
R CMD INSTALL -l ./library packages/lattice_0.20-33.tar.gz\
    packages/Matrix_1.2-2.tar.gz \
    packages/KernSmooth_2.23-15.tar.gz \
    packages/bitops_1.0-6.tar.gz \
    packages/caTools_1.17.1.tar.gz \
    packages/gtools_3.5.0.tar.gz \
    packages/gdata_2.17.0.tar.gz \
    packages/gplots_2.17.0.tar.gz \
    packages/ROCR_1.0-7.tar.gz
R CMD INSTALL -l ./library packages/DRaWR/
```

To use the 5 species networks
```
gunzip networks/5ins_cdhmw.names.edge.gz
gunzip networks/5sp_adhiw.names.edge.gz
```

## Tutorial

### Gene Sets
### Heterogeneous Networks
### Example Runs
Load libraries
```
.libPaths("library")
library(DRaWR)
```

Simple example no cross validation
```
DRaWR(possetfile = "packages/DRaWR/data/sample_inputs/test.setlist",
    unifile = "packages/DRaWR/data/sample_inputs/test.uni",
    networkfile = "packages/DRaWR/data/sample_inputs/test.edge",
    outdir = "packages/DRaWR/data/sample_outs/results_",
    restarts = c(.3), nfolds = 1, st2keep = 1,
    undirected = TRUE, unweighted = FALSE, normalize = "type",
    maxiters = 50, thresh = 0.000001,
    property_types = c("T1", "T2"), writepreds = 1)
```

12 fly gene sets with no cross validation
```
DRaWR(possetfile = "setlists/dmel_possets12.list.txt",
    unifile = "gene_sets/universes/dmel.ids.uni.txt",
    networkfile = "networks/dmel_cdhmw.names.edge",
    outdir = "sample_results/",
    restarts = c(.3), nfolds = 1, st2keep = 50,
    undirected = TRUE, unweighted = FALSE, normalize = "type",
    maxiters = 50, thresh = 0.0001,
    property_types = c("chip_binding", "motif_u5", "pfam_domain"), writepreds = 1)
```

12 fly gene sets with 4fold cross validation
```
DRaWR(possetfile = "setlists/dmel_possets12.list.txt",
    unifile = "gene_sets/universes/dmel.ids.uni.txt",
    networkfile = "networks/dmel_cdhmw.names.edge",
    outdir = "sample_results/",
    restarts = c(.3), nfolds = 4, st2keep = 50,
    undirected = TRUE, unweighted = FALSE, normalize = "type",
    maxiters = 50, thresh = 0.0001,
    property_types = c("chip_binding", "motif_u5", "pfam_domain"), writepreds = 0)
```

Combined aggression set on 5sp aggression network without cross validation
```
DRaWR(possetfile = "setlists/3sps.fdr.1.list.txt",
    unifile = "gene_sets/universes/3sps_ids.uni.txt",
    networkfile = "networks/5sp_adhiw.names.edge",
    outdir = "sample_results/",
    restarts = c(.3), nfolds = 1, st2keep = 50,
    undirected = TRUE, unweighted = FALSE, normalize = "type",
    maxiters = 50, thresh = 0.0001,
    property_types = c("allen_brain_atlas", "gene_ontology","pfam_domain"), writepreds = 0)
```

### From command line
Simple example no cross validation
```
Rscript -e ".libPaths(\"library\"); library(\"DRaWR\"); \
    DRaWR(possetfile = \"packages/DRaWR/data/sample_inputs/test.setlist\", \
    unifile = \"packages/DRaWR/data/sample_inputs/test.uni\", \
    networkfile = \"packages/DRaWR/data/sample_inputs/test.edge\", \
    outdir = \"packages/DRaWR/data/sample_outs/cli_results_\", \
    restarts = c(.3), nfolds = 1, st2keep = 1, \
    undirected = TRUE, unweighted = FALSE, normalize = \"type\", \
    maxiters = 50, thresh = 0.000001, \
    property_types = c(\"T1\", \"T2\"), writepreds = 1) "
```

## Resources

### Function Parameters
#### relating to input files
* *possetfile (string)*: location of file containing relative location of gene sets to test.
* *unifile (string)*: location of file listing gene universe.
* *networkfile (string)*: location of file containing network contents.

#### relating to output files
* *outdir (string)*: prefix of location of file to write performance results (optionally prediction results).
* *writepreds (boolean)*: write predictions out to a file. Default is FALSE.
* *nfolds (int)*: number of folds for cross validation, Default is 1, no cross-validation.

#### relating to network processing
* *undirected (bool)*: boolean to make network undirected. Default is TRUE.
* *unweighted (bool)*: boolean to make network unweighted. Default is FALSE.
* *normalize (string)*: "type" or "none" normalization method. Default is 'type'.
* *property_types (vector)*: vector containing *ALL* names of property->gene (PG) edge_types. May contain PG edge_types that are not in the given network, but should not contain gene<->gene edge_types. Default is c("allen_brain_atlas", "chip_binding", "gene_ontology", "motif_u5", "pfam_domain", "T1", "T2").
* *st2keep (int)*: number of property nodes to keep in second stage for each property type. To skip the second stage of DRaWR can set st2keep to 0. Default is 1.

#### relating to random walks
* *restarts (vector)*: vector of restart values to test. Default is c(0.7).
* *maxiters (int)*: maximum number of allowable iterations. Default is 50.
* *thresh (float)*: threshold for L1 norm convergence. Default is 0.0001.

### Input File Formats
#### Network Edge File
The network edge file is expected to be a 4 column, tab separated file. Each row represents and edge in the network and contains the values for:
*col1: node1 (string)
*col2: node2 (string)
*col3: edge_weight (float)
*col4: edge_type (string)

For edges between a property node and a gene node, the property node must be listed as node 1.  The order of the nodes does not matter for gene-gene edges.  The edge weights must be positive with large values meaning stronger relationships.  The second stage of DRaWR is intended to limit the property nodes from many different property->gene (PG) edge_types.  The number of gene<->gene (GG) edge_types should be limited.  By default the network edge file will be converted to an undirected, normalized by edge_type adjacency matrix with the edge weights taken as the maximum if there is repetition.

Example network edge file:
```
G1	G6	0.76	typeGG
G2	G1	0.41	typeGG
G2	G7	0.73	typeGG
G4	G5	0.89	typeGG
P1	G2	0.57	typePG.1
P1	G6	0.30	typePG.1
P2	G1	0.07	typePG.1
P2	G4	0.89	typePG.1
P3	G3	0.66	typePG.1
P3	G6	0.80	typePG.1
P4	G4	0.24	typePG.2
P4	G1	0.74	typePG.2
P5	G5	0.20	typePG.2
P5	G2	0.95	typePG.2
P5	G7	0.92	typePG.2
```

#### Gene Set File / Gene Universe File
The query gene set files or gene universe files should one node name listed on each row.  If the gene nodes have distinct weights, a second column (separated by a tab) can be added.  Node weights must be positive with larger values meaning stronger evidence.

Example gene set file:
```
G3  2.5
G4  1.3
G5  4.0
```

The gene universe file must contain the appropriate universe of genes that the genes in the query gene sets may be chosen from. Typically, this will be all of the genes on the species of interest.

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
#### File of results statistics
This file ending in '.stats' is produced to analyzing the performance results.

For every:
* stage ("stage1") of DRaWR
* fold ("iter") of cross validation
* RWR restart parameter ("restart") setting
* query gene set ("posset")

this file outputs on every line the settings of the DRaWR run as well as the:
* AUROC evaluation metric ("aucval")
* number of RWR interations ("rwr_iters")
* number of training examples ("ntrain")

#### File of predictions
This file ending in '.rwr' contains information about each node in the network after running DRaWR.  This file is only produced if writepreds = TRUE.  There is one row for each node in the network, and the columns are:
* "node": the name of that node
* "type": the type of that node, -1 if node is a gene node
* "universe": 1 if node in gene universe set, 0 otherwise
* "baseline": probability of being in the node from converged baseline RWR
* "train": 1 if node in the gene query set training set, 0 otherwise
* "test": 1 if node in the gene query set test set, 0 otherwise
* "stage1": probability of being in the node from converged stage1 RWR
* "diff": difference between "stage1" and "baseline"
* "keep": 1 if the property node was kept for the stage two RWR, 0 otherwise. -1 if node is a gene node
* "stage2": probability of being in the node from converged stage2 RWR

The file ending in '.base' contains the first 4 columns of the '.rwr' file and is used to accelerate the convergence of the RWRs.

