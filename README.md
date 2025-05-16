# Multi-target prediction based on ChEMBL data

This is a ligand-based target prediction model developed by the ChEMBL team. They trained the model using pairs of small molecules and their protein targets, and produced a multitask predictor. The thresholds of activity where determined by protein families (kinases: <= 30nM,  GPCRs: <= 100nM, Nuclear Receptors: <= 100nM, Ion Channels: <= 10μM, Non-IDG Family Targets: <= 1μM). Here we provide the model trained on ChEMBL\_28, which showed an accuracy of 85%.

This model was incorporated on 2021-05-07.

## Information
### Identifiers
- **Ersilia Identifier:** `eos1vms`
- **Slug:** `chembl-multitask-descriptor`

### Domain
- **Task:** `Representation`
- **Subtask:** `Featurization`
- **Biomedical Area:** `Any`
- **Target Organism:** `Not Applicable`
- **Tags:** `Bioactivity profile`, `Target identification`, `ChEMBL`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `1`
- **Output Consistency:** `Fixed`
- **Interpretation:** Probability of having the protein (identified by ChEMBL ID), as target

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| chembl1075104 | float | high | Predicted probability of binding to protein CHEMBL1075104 |
| chembl1075110 | float | high | Predicted probability of binding to protein CHEMBL1075110 |
| chembl1075126 | float | high | Predicted probability of binding to protein CHEMBL1075126 |
| chembl1075138 | float | high | Predicted probability of binding to protein CHEMBL1075138 |
| chembl1075145 | float | high | Predicted probability of binding to protein CHEMBL1075145 |
| chembl1075189 | float | high | Predicted probability of binding to protein CHEMBL1075189 |
| chembl1075232 | float | high | Predicted probability of binding to protein CHEMBL1075232 |
| chembl1075317 | float | high | Predicted probability of binding to protein CHEMBL1075317 |
| chembl1163101 | float | high | Predicted probability of binding to protein CHEMBL1163101 |
| chembl1163125 | float | high | Predicted probability of binding to protein CHEMBL1163125 |

_10 of 616 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos1vms](https://hub.docker.com/r/ersiliaos/eos1vms)
- **Docker Architecture:** `AMD64`, `ARM64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos1vms.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos1vms.zip)

### Resource Consumption
- **Model Size (Mb):** `18`
- **Environment Size (Mb):** `590`


### References
- **Source Code**: [https://github.com/chembl/chembl_multitask_model/](https://github.com/chembl/chembl_multitask_model/)
- **Publication**: [http://chembl.blogspot.com/2019/05/multi-task-neural-network-on-chembl.html](http://chembl.blogspot.com/2019/05/multi-task-neural-network-on-chembl.html)
- **Publication Type:** `Other`
- **Publication Year:** `2019`
- **Ersilia Contributor:** [miquelduranfrigola](https://github.com/miquelduranfrigola)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [MIT](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos1vms
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos1vms
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
