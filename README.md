# Target Prediction

This is a ligand-based target prediction model developed by the ChEMBL team. They trained the model using pairs of small molecules and their protein targets, and produced a predictor for two cut-offs: 1uM and 10uM. The model predicts the main target of a small molecule with an accuray of 69%. These predictors are available since ChEMBL_18 and have been subsequently updated since. You can read more about it in the ChEMBL [blogpost](https://chembl.github.io/ligand-based-target-predictions-in/). Here we use the ChEMBL_25 version (2019)

## Summary

* Predicts **protein targets** of small molecules
* Takes **compound structures** as input
* Trained with **experimental** bioactivity data from ChEMBL
* Based on a dataset of **>1000** molecules
* Results **validated in-silico**
* Processed data can be downloaded [here](https://github.com/chembl/target_predictions)

## Specifications

* Input: SMILES string (also accepts an InChIKey string or a molecule name string, and converts them to SMILES) 
* Endpoint: probability that a certain protein (identified by ChEMBL Id) is the target of the small molecule
* Results interpretation: 0 to 1

## History

1. Model was downloaded on 06/05/2021 from [ChEMBL](http://ftp.ebi.ac.uk/pub/databases/chembl/target_predictions/NB/) following the direct link.
4. Model was incorporated to Ersilia on 06/05/2021.
