# Fairness between Afrian Americans (AA) and white Americans (WA) in RSFC-based behavioral prediction using the ABCD dataset

## Reference

Jingwei Li, Danilo Bzdok, Jianzhong Chen, Angela Tam, Leon Qi Rong Ooi, Avram J. Holmes, Tian Ge, Kaustubh R. Patil, Mbemba Jabbi, Simon B. Eickhoff, B.T. Thomas Yeo*, Sarah Genon*, _**Cross-ethnicity/race generalization failure of behavioral prediction from resting-state functional connectivity**_, under review.

## Background

Algorithmic biases that favor majority populations pose a key challenge to the application of machine learning for precision medicine. Here, we assessed such bias in prediction models of behavioral phenotypes from brain functional magnetic resonance imaging. We examined 
the prediction bias using two independent datasets (pre-adolescent versus adult) of mixed ethnic/racial composition. When predictive models were trained on data dominated by white Americans (WA), out-of-sample prediction errors were generally higher in African Americans (AA) than for WA. This bias towards WA corresponds to more WA-like brain-behavioral association patterns learned by models. When models were trained on AA only, compared to training only on WA or an equal number of AA and WA participants, AA prediction accuracy improved but stayed below that for WA. Overall, the results point to the need for caution and further research regarding the application of current brain-behavior prediction models in minority population.

## Usage

First, this repository relies on multiple utility functions in the Computational Brain Imaging Group repository (CBIG; https://github.com/ThomasYeoLab/CBIG), e.g. kernel ridge regression package. Please follow the configuration instructions of CBIG repository before you use the current repository. Also, make sure you have the HCP csv files prepared on your devices.

After that, this repository should be used as the following steps:

1. Run `HCP_addpath` when everytime you open a new matlab session, to add all subfolders of the current repository into your matlab paths.

2. Follow the README in `preparation` folder to collect all the necessary phenotypes and resting-state functional connectivity matrices.

3. Follow the README in `match_split` folder to find matched AA and WA pairs and split subjects into training versus test folds.

4. Follow the README in `KRR` folder to perform kernel ridge regression.

5. Follow the README in `AAvsWA` folder to calculate prediction accuracies of matched AA and WA for the KRR models trained on whole populations.

6. Follow the README in `Interpretation` folder to calculate model-learned brain-behavioral associations and true brain-behavioral associations.

7. Use the scripts in `plot` folder to create whisker plots, scatter plots, etc.
