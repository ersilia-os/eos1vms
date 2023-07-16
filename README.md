# Multi-target prediction based on ChEMBL data

This is a ligand-based target prediction model developed by the ChEMBL team. They trained the model using pairs of small molecules and their protein targets, and produced a predictor for two cut-offs: 1uM and 10uM. The model predicts the main target of a small molecule with an accuracy of 69%. These predictors are available since ChEMBL\_18 and have been subsequently updated since. You can read more about it in the ChEMBL blogpost. Here we use the ChEMBL\_25 version (2019).

## Identifiers

* EOS model ID: `eos1vms`
* Slug: `chembl-multitask-descriptor`

## Characteristics

* Input: `Compound`
* Input Shape: `Single`
* Task: `Classification`
* Output: `Text`
* Output Type: `String`
* Output Shape: `Flexible List`
* Interpretation: Main protein target. Thresholds at 10 uM and 1 uM

## References

* [Publication](https://chembl.github.io/ligand-based-target-predictions-in/)
* [Source Code](https://github.com/chembl/target_predictions)
* Ersilia contributor: [miquelduranfrigola](https://github.com/miquelduranfrigola)

## Ersilia model URLs
* [GitHub](https://github.com/ersilia-os/eos1vms)
* [AWS S3](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos1vms.zip)
* [DockerHub](https://hub.docker.com/r/ersiliaos/eos1vms) (AMD64, ARM64)

## Citation

If you use this model, please cite the [original authors](https://chembl.github.io/ligand-based-target-predictions-in/) of the model and the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff).

## License

This package is licensed under a GPL-3.0 license. The model contained within this package is licensed under a Apache-2.0 license.

Notice: Ersilia grants access to these models 'as is' provided by the original authors, please refer to the original code repository and/or publication if you use the model in your research.

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission!