# GNorm2
***
GNorm2 is a gene name recognition and normalization tool with optimized functions and customizable configuration to the user preferences. The GNorm2 integrates multiple deep learning-based methods and achieves state-of-the-art performance. GNorm2 is freely available to download for stand-alone usage. [Download GNorm2 here](https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/tmTools/download/GNorm2/GNorm2.tar.gz)

## Content
- [Dependency package](#package)
- [Introduction of folders](#intro)
- [Running GNorm2](#pipeline)

## Dependency package
<a name="package"></a>
The codes have been tested using Python3.8/3.9 on CentOS and uses the following main dependencies on a CPU and GPU:
- [TensorFlow 2.3.0](https://www.tensorflow.org/)
- [Transformer 4.18.0](https://huggingface.co/docs/transformers/installation)
- [stanza 1.4.0](stanfordnlp.github.io/stanza/)

To install all dependencies automatically using the command:

    $ pip install -r requirements.txt


## Introduction of folders
<a name="intro"></a>

- src_python
	- GeneNER: the codes for gene recognition
	- SpeAss: the codes for species assignment
- src_Java
	- GNormPluslib : the codes for gene normalization and species recogntion
- GeneNER_SpeAss_run.py: the script for runing pipeline
- GNormPlus.jar: the upgraded GNormPlus tools for gene normalization
- gnorm_trained_models:pre-trianed models and trained NER/SA models
	- bioformer-cased-v1.0: the original bioformer model
	- BiomedNLP-PubMedBERT-base-uncased-abstract: the original pubmedbert model
	- geneNER
		- GeneNER-Bioformer/PubmedBERT-Allset.h5: the Gene NER models trained by all datasets
		- GeneNER-Bioformer/PubmedBERT-Trainset.h5: the Gene NER models trained by the training set only
	- SpeAss
		- SpeAss-Bioformer/PubmedBERT-SG-Allset.h5: the Species Assignment models trained by all datasets
		- SpeAss-Bioformer/PubmedBERT-SG-Trainset.h5: the Species Assignment models trained by the trianing set only
	- stanza
		- downloaded stanza library for offline usage
- vocab: label files for the machine learning models of GeneNER and SpeAss
- Dictionary: The dictionary folder contains all required files for gene normalization
- CRF: CRF++ library (called by GNormPlus.sh)
- Library: Ab3P library
- tmp/tmp_GNR/tmp_SA/tmp_SR folders: temp folder
- input/output folders: input and output folders. BioC (abstract or full text) and PubTator (abstract only) formats are both avaliable.
- GNorm2.sh: the script to run GNorm2
- setup.GN.txt/setup.SR.txt/setup.txt the setup files for GNorm2.

## Running GNorm2
<a name="pipeline"></a>
Please firstly download [GNorm2](https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/tmTools/download/GNorm2/GNorm2.tar.gz) to your local.
Below are the well-trained models (i.e., PubmedBERT/Bioformer) for Gene NER and Species Assignment.
Models for Gene NER:
- gnorm_trained_models/geneNER/GeneNER-PubmedBERT.h5
- gnorm_trained_models/geneNER/GeneNER-Bioformer.h5
Models for Species Assignment:
- gnorm_trained_models/SpeAss/SpeAss-PubmedBERT.h5
- gnorm_trained_models/SpeAss/SpeAss-Bioformer.h5

The parameters of the input/output folders:

- INPUT, default="input"
- OUTPUT, default="output"

[BioC-XML](bioc.sourceforge.net) or [PubTator](https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/tmTools/Format.html) formats are both avaliabel to GNorm2.

1. Run GNorm2

Run Example:

    $ ./GNorm2.sh input output

## Acknowledgments
This research was supported by the Intramural Research Program of the National Library of Medicine (NLM), National Institutes of Health.
