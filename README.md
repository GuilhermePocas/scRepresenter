

# scRepresenter: a unified framework for computing and benchmarking biologically informed cellular representations in single-cell transcriptional data

![scRepresenter Overview](docs/methodology.png)

In single-cell RNA-sequencing (scRNA-seq) analysis, the quality of cellular embeddings is crucial for accurate downstream biological interpretation. This issue is especially significant in neural tissues, which display a wide diversity of cell types, and in disease contexts characterized by transient cellular states. Foundation models, when applied individually, effectively capture global transcriptomic structures across large-scale datasets. In contrast, biological knowledge–guided approaches use struc-tured priors to enhance interpretability and local coherence. However, despite the complementary ad-vantages of these methods, there is currently no unified framework for their systematic comparison and integration.

To address this gap, we introduce scRepresenter—an open-source, systematic workflow for computing, integrating, and validating multiple cellular embedding strategies within standard single-cell analysis pipelines. scRepresenter supports four categories of embedding representations: (1) expres-sion-based, (2) knowledge-guided, (3) foundation model–derived, and (4) hybrid approaches that com-bine global transcriptomic context with structured biological priors. All strategies are evaluated under harmonized preprocessing and benchmarking conditions to ensure fair and systematic comparison. We benchmarked scRepresenter using data from a neurodegenerative disease model based on human brain organoids—multicellular, patient-specific, lab-grown cortical tissues. Our results show that hybrid embeddings consistently achieve higher macro F1-scores and better structural preservation, especially for rare and transcriptionally similar cell populations. Furthermore, unsupervised analyses revealed en-hanced recovery of reference cell-type structures across a range of clustering resolutions. By enabling transparent comparison and principled integration of complementary embedding paradigms, scRepre-senter provides a robust and systematic approach for interpreting complex single-cell datasets, particu-larly in challenging biological contexts. 

<!-- **scRepresenter** is a computational workflow geared towards cell type classification, with a particular focus on accuratly predicting rare cell types. It achieves this by bridging large-scale, expression-only Foundation Models, known for their accurate predictions, with domain-specific biological knowledge drawn from large graphs of relevant information.

The workflow employs two methods for this purpose, the biological knowledge integration method [scNET, 2025](https://github.com/madilabcode/scNET), and the Foundation Model [scGPT, 2024](https://github.com/bowang-lab/scGPT). -->


##  Installation

To clone the repository, use the following command:
```
git clone https://github.com/GuilhermePocas/scRepresenter.git
```

A docker enviroment is recommended to avoid versioning issues. For this, we provide two Dockerfiles depending on the availability of GPU capabilities, which is highly recommended. In this case, run the following commands:

```
docker build -f Dockerfile.gpu -t env .
docker run -it --gpus all --rm -v $(pwd)/output:/app/output env bash
```

If there is no available GPU, then use the CPU-based dockerfile trough the following command:
```
docker build -f Dockerfile.cpu -t env .
docker run -it --rm -v $(pwd)/output:/app/output env bash
```

Alternatively, you can install the packages directly in your computer. This requires Python <= 3.9 and should be done in a separate python enviroment.
```
python -m venv .venv
source venv/bin/activate
pip install -r requirements.txt
```

##  Usage

### Data preprocessing

In order to get the best performance, preprocessing the scRNA dataset is highly recommended. See the Preprocessing script for an example.

### Networks

scRepresenter utilizes a gene similarity graph to run the scNET model. Four pre-made graphs can already be found in ``` ./src/networks ``` , constructed from GO annotations, Protein-Protein Interactions and Transcription Factor proteins. 

To build custom gene similarity graphs (barring the PPI graph), use the methods found in ``` ./src/networks ```. 

For the GO-based methods, a GO annotations file is required, found in the GO archive https://release.geneontology.org/. The human file was used, ```goa_human.gaf```, from version 2025-02-06.
For the Transcription Factors method, the full interaction table of the TFLink database was used, found in https://tflink.net/download/. The small and large-scale full interaction table was used for the human organism. Both of these files need to be placed in ``` ./src/networks/knowledge sources ``` to run the corresponding methods, which can be done with the following command:

``` 
python create_METHOD_graph.py --top_p 100 --var_genes 2000
```

Whith the following args:

- **top_p(int)**: the number of pairs to consider for each gene, higher number might generate really large graphs.

- **var_genes(int)**: the number of highly variable genes to be selected.

For detailed steps on how to create each of the similarity graphs, see the Gene Networks script.

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

For a detailed example see the Training script, where scRepresented is trained on the PBMC3k dataset from Scanpy.

<!--An example script is provided in ```/main.py```, using the pbmc3k dataset from Scanpy, and the human scGPT model checkpoint. To run the example, start the docker enviroment as previously explained, and run the following command:

```
python main.py --model_name test --network PPI.csv --scnet_epochs 800 --scgpt_epochs 60
```

with the following args:

- **model_name(str)**: the name of the current run.

- **network**: the name of the gene similarity network to be used in this run. They can be found ```./src/networks``` 

- **scnet_epochs**: the number of training steps scNET will run for. If its 0 this model is skipped.

- **scgpt_epochs**: the number of training steps scGPT will run for. If its 0 this model is skipped.-->

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