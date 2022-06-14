# ChEMBL multi-task descriptor

Protein target prediction based on ChEMBL data

| Description | Input  | Output  | Training Data | Experimental Validation |
| ------- | --- | --- | --- | --- |
| A ligand-based prediction model by the ChEMBL team. The model predicts the main target of a small molecule with an accuracy of 69%  | SMILES | Protein target | ChEMBL_18 and 1244 targets | No |

## Source code
This model has been published by George Papadatos. Ligand-based target predictions. *ChEMBL* https://chembl.github.io/ligand-based-target-predictions-in/# (2014)

Code: https://github.com/chembl/target_predictions

## Extended description
This model is a ligand-based target prediction model that learns what features of ligands have mutual relations with activity against a certain target and assign a score to each of these features. 

### Summary
- Trained using pairs of compounds and protein targets
- Predict targets for small molecule compounds
- Licensed using Apache 2.0 License

### Specification
- Input: SMILES compound
- Output: protein targets

## History
- Model was downloaded on September 14, 2021
- Model was incorporated on September 15, 2021

