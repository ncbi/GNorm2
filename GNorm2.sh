#!/bin/sh

#
# ./GNorm2.sh input output
#
INPUT=input
OUTPUT=output

# SR
java -Xmx60G -Xms30G -jar GNormPlus.jar ${INPUT} tmp_SR setup.SR.txt

# GNR+SA (CUDA_VISIBLE_DEVICES=3 )
python GeneNER_SpeAss_run.py -i tmp_SR -r tmp_GNR -a tmp_SA -n gnorm_trained_models/geneNER/GeneNER-Bioformer.h5 -s gnorm_trained_models/SpeAss/SpeAss-Bioformer.h5

# GN
java -Xmx60G -Xms30G -jar GNormPlus.jar tmp_SA ${OUTPUT} setup.GN.txt

rm tmp_SR/*
rm tmp_GNR/*
rm tmp_SA/*
