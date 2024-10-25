# Mode-of-inheritance prediciton.

We developed an ensemble machine learning model for variant-level mode-of-inheritance prediction (MOI-Pred). Additionally, we aggregated variant-level mode-of-inheritance prediction from two independent ensemble approached to develop a consensus predictor of mode-of-inheritance (ConMOI).

MOI-Pred (Mode Of Inheritance-Predictor) is a three-way variant-level mode of inheritance prediction tool that specifically tackles recessive identification of missense variants. MOI-Pred uses a random forest model to predict variants pathogenic for autosomal recessive (AR) disease, pathogenic for autosomal dominant (AD) disease, or benign. The model learns from 78 functional, evolutionary and population frequency annotations.
We provide pre-computed MOI-Pred predictions for 71M missense variants in the human genome (hg38), please use: https://doi.org/10.5281/zenodo.5565246.

ConMOI is a consensus approach that combines predictions from MOI-Pred with two previously published methods that produce predictions incorporating mode of inheritance, MAPPIN (PMID 28977528) and MAVERICK (PMID 37443090). We provide pre-computed ConMOI predictions for 53M missense variants in the human genome (hg38), please use: https://doi.org/10.5281/zenodo.5565246.

More details about our approach to mode-of-inheritance prediciton, including the methods, features and performance assessment are described in our preprint:
Ben O. Petrazzini, Daniel J. Balick, Iain S. Forrest, Judy Cho, Ghislain Rocheleau, Daniel M. Jordan*, Ron Do*. Prediction of recessive inheritance for missense variants in human disease. _medRxiv_ 2021.10.25.21265472; doi: https://doi.org/10.1101/2021.10.25.21265472.
