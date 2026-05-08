
# scRepresenter: a unified framework for computing and benchmarking biologically informed cellular representations in single-cell transcriptional data

## Overview

![scRepresenter Overview](docs/methodology.png)

**scRepresenter** is a computational workflow geared towards cell type classification, with a particular focus on accuratly predicting rare cell types. It achieves this by bridging large-scale, expression-only Foundation Models, known for their accurate predictions, with domain-specific biological knowledge drawn from large graphs of relevant information.

The workflow employs two methods for this purpose, the biological knowledge integration method [scNET, 2025](https://github.com/madilabcode/scNET), and the Foundation Model [scGPT, 2024](https://github.com/bowang-lab/scGPT).


##  Installation

### **1. Setting up**
To clone the repository, use the following command:
```
git clone https://github.com/GuilhermePocas/scRepresenter.git
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


##  Usage

### scRepresenter

To run the scRepresenter pipeline, first load a scRNAseq dataset using Scanpy. The pipeline automatically downloads the Human scGPT checkpoint, in order to use another organ checkpoint from https://github.com/bowang-lab/scGPT/tree/main#pretrained-scGPT-checkpoints download it and place it in ```./src/scgpt/checkpoints```. Then, run the following function:

```
scnet_embs, scgpt_embs, avg_embs, conq_embs, labels = run_scRepresenter(model_name, data, network, dir, scnet_epochs, scgpt_epochs)
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

- **avg_embs**: the scRepresenter embeddings with the average aggregation strategy.

- **conq_embs**: the scRepresenter embeddings with the concatenation aggregation strategy.

- **labels**: the cell type labels, in the same order as the embeddings.

An example script is provided in ```/main.py```, using the pbmc3k dataset from Scanpy, and the human scGPT model checkpoint. To run the example, start the docker enviroment as previously explained, and run the following command:

```
python main.py --model_name test --network PPI.csv --scnet_epochs 800 --scgpt_epochs 60
```

with the following args:

- **model_name(str)**: the name of the current run.

- **network**: the name of the gene similarity network to be used in this run. They can be found ```./src/networks``` 

- **scnet_epochs**: the number of training steps scNET will run for. If its 0 this model is skipped.

- **scgpt_epochs**: the number of training steps scGPT will run for. If its 0 this model is skipped.

### Networks

scRepresenter utilizes a gene similarity graph to run the scNET model. Four pre-made graphs can already be found in ``` ./src/networks ```, constructed from GO annotations, Protein-Protein Interactions and Transcription Factor proteins. 

To build custom gene similarity graphs (barring the PPI graph), use the methods found in ``` ./src/networks ```. 

For the GO-based methods, a GO annotations file is required, found in the GO archive https://release.geneontology.org/. The human file was used, ````goa_human.gaf```, from version 2025-02-06.
For the Transcription Factors method, the full interaction table of the TFLink database was used, found in https://tflink.net/download/. The small and large-scale full interaction table was used for the human organism. Both of these files need to be placed in ``` ./src/networks/knowledge sources ``` to run the corresponding methods, which can be done with the following command:

``` 
python create_METHOD_graph.py --top_p 100 --var_genes 2000
```

Whith the following args:

- **top_p(int)**: the number of pairs to consider for each gene, higher number might generate really large graphs.

- **var_genes(int)**: the number of highly variable genes to be selected.

### Embedding analysis

A shiny application is also provided to better analyse the embedding objects produced by scRepresenter, allowing for the easy viewing of several insights and the creation of different types of plots and graphs. The application can be run locally from a Docker container, starting with building the docker image from the provided Dockerfile in the ```/app``` directory using the following command:

```
docker build -t shiny-app .
```

Then, by simply running the docker container:

```
docker run -p 3838:3838 shiny-app
```

The application can then acessed trough ```http://localhost:3838/``` using your internet browser. Any H5 embedding object produced by the pipeline can be uploaded here, and the specific embeddings described in the scRepresenter paper can be found at https://drive.google.com/file/d/19NVlnrTCTDOW9tHsVgjekZEyGVTTS1wJ/view?usp=sharing.