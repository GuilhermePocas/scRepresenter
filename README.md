
# scBLOOM: Improving Cell Type Classification by Bridging Single-cell Expression and Biological Knowledge

## Overview

![scBLOOM Overview](docs/methodology.png)

**scBLOOM** (single-cell Biologically-enhanced Learning from Ontologies and Omics-based Models) is a computational workflow geared towards cell type classification, with a particular focus on accuratly predicting rare cell types. It achieves this by bridging large-scale, expression-only Foundation Models, known for their accurate predictions, with domain-specific biological knowledge drawn from large graphs of relevant information.

The workflow employs two methods for this purpose, the biological knowledge integration method [scNET, 2025](https://github.com/madilabcode/scNET), and the Foundation Model [scGPT, 2024](https://github.com/bowang-lab/scGPT).


##  Installation

### **1. Setting up**
To clone the repository, use the following command:
```
git clone https://github.com/GuilhermePocas/scBLOOM.git
```

The required packages can be installed in a docker enviroment, with two Dockerfiles available depending on the availability of the GPU.

If there is a GPU available, which is recommended, then run the following commands:

```
docker build -f Dockerfile.gpu -t env .
docker run -it --gpus all --rm -v $(pwd)/output:/app/output env bash
```

If there is no GPU available, run the following command:
```
docker build -f Dockerfile.cpu -t env .
docker run -it --rm -v $(pwd)/output:/app/output env bash
```


##  API

### scBLOOM

To run the scBLOOM pipeline, first load a scRNAseq dataset using Scanpy, and then choose the appropriate scGPT model checkpoint from https://github.com/bowang-lab/scGPT/tree/main#pretrained-scGPT-checkpoints, and place it in ```./src/scgpt/checkpoints```. Then, run the following function:

```
scnet_embs, scgpt_embs, avg_embs, conq_embs, labels = run_scBLOOM(model_name, data, network, dir, scnet_epochs, scgpt_epochs)
```

with the following args:

- **model_name(str)**: the name of the current run.

- **data(AnnData)**: a scRNAseq AnnData object.

- **network(str)**: the file path of the gene similarity network.

- **dir(str)**: the output directory where the results and embeddings will be saved.

- **scnet_epochs(int)**: the number of steps when training scNET.

- **scgpt_epochs(int)**: the number of steps when training scGPT.

The resulting output objects are:

- **scnet_embs**: the embeddings from the scNET model.

- **scgpt_embs**: the embeddings from the scGPT model.

- **avg_embs**: the scBLOOM embeddings with the average aggregation strategy.

- **conq_embs**: the scBLOOM embeddings with the concatenation aggregation strategy.

- **labels**: the cell type labels, in the same order as the embeddings.

An example script is provided in ```/main.py```, using the pbmc3k dataset from Scanpy, and the human scGPT model checkpoint. To run the example, start the docker enviroment as previously explained, and run the following command:

```
python main.py --model_name test --network PPI.csv --scnet_epochs 800 --scgpt_epochs 60
```

### Networks

scBLOOM utilizes a gene similarity graph to run the scNET model. Four pre-made graphs can already be found in ``` ./src/networks ```, constructed from GO annotations, Protein-Protein Interactions and Transcription Factor proteins. 

To build custom gene similarity graphs (barring the PPI graph), use the methods found in ``` ./src/networks ```. 

For the GO-based methods, a GO annotations file is required, found in the GO archive https://release.geneontology.org/. The human file was used, ````goa_human.gaf```, from version 2025-02-06.
For the Transcription Factors method, the full interaction table of the TFLink database was used, found in https://tflink.net/download/. The small and large-scale full interaction table was used for the human organism. Both of these files need to be placed in ``` ./src/networks/knowledge sources ``` to run the corresponding methods, which can be done with the following command:

``` 
python create_METHOD_graph.py --top_p 100 --var_genes 2000
```

Whith the following args:

- **top_p(int)**: the number of pairs to consider for each gene, higher number might generate really large graphs.

- **var_genes(int)**: the number of highly variable genes to be selected.


