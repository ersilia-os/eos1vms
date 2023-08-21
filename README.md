# Multi-target prediction based on ChEMBL data

This is a ligand-based target prediction model developed by the ChEMBL team. They trained the model using pairs of small molecules and their protein targets, and produced a multitask predictor. The thresholds of activity where determined by protein families (kinases: <= 30nM,  GPCRs: <= 100nM, Nuclear Receptors: <= 100nM, Ion Channels: <= 10μM, Non-IDG Family Targets: <= 1μM). Here we provide the model trained on ChEMBL\_28, which showed an accuracy of 85%.

## Identifiers

* EOS model ID: `eos1vms`
* Slug: `chembl-multitask-descriptor`

## Characteristics

* Input: `Compound`
* Input Shape: `Single`
* Task: `Classification`
* Output: `Probability`
* Output Type: `Float`
* Output Shape: `List`
* Interpretation: Probability of having the protein (identified by ChEMBL ID), as target

## References

* [Publication](http://chembl.blogspot.com/2019/05/multi-task-neural-network-on-chembl.html)
* [Source Code](https://github.com/chembl/chembl_multitask_model/)
* Ersilia contributor: [miquelduranfrigola](https://github.com/miquelduranfrigola)

## Ersilia model URLs
* [GitHub](https://github.com/ersilia-os/eos1vms)
* [AWS S3](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos1vms.zip)
* [DockerHub](https://hub.docker.com/r/ersiliaos/eos1vms) (AMD64, ARM64)

## Citation

If you use this model, please cite the [original authors](http://chembl.blogspot.com/2019/05/multi-task-neural-network-on-chembl.html) of the model and the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff).

## License

This package is licensed under a GPL-3.0 license. The model contained within this package is licensed under a None license.

Notice: Ersilia grants access to these models 'as is' provided by the original authors, please refer to the original code repository and/or publication if you use the model in your research.

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission!