# ChEMBL Multi-task descriptor

Protein target prediction based on ChEMBL data

| Description | Input  | Output  | Training Data | Experimental Validation |
| ------- | --- | --- | --- | --- |
| Predicts targets of small molecule compounds | SMILES | Protein target | ChEMBL_18 and 1244 targets | No |

## Source code
This model has been published by George, Ligand-based target predictions in ChEMBL April 03, 2014 https://chembl.github.io/ligand-based-target-predictions-in/#

Code: https://github.com/chembl/target_predictions

## Extended description
This model is a ligand-based target prediction model that learns what substructural features (encoded as fingerprints) of ligands correlate with activity against a certain target and assign a score to each of these features. It also comes up with a sorted list of likely targets with the highest scores.

### Summary
- Trained using ligand information only
- Predict targets for small molecule drugs
- Licensed using Apache 2.0 License

### Specification
- Input: SMILES compound
- Output: protein targets

## History
- Model was downloaded on September 14, 2021
- Model was incorporated on September 15, 2021

