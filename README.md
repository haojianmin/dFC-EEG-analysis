## dFC-EEG-analysis

This repository contains the MATLAB code for dynamic functional connectivity (dFC) feature extraction and statistical analysis used in the paper.

## System Requirements

- MATLAB R2021b or later
- Statistics and Machine Learning Toolbox
- Bioinformatics Toolbox
- EEGLAB

## Usage

Run the following commands in MATLAB:
matlab

cfg = config_analysis();

run_dFC_preprocessing(cfg);

run_dFC_feature_analysis(cfg);

The first script computes dynamic wPLI matrices for each subject and saves them
to derivatives/. The second script extracts dFC features, performs ANOVA with
FDR correction, and generates the bar plots shown in the paper.
Configuration

Data Availability

The raw EEG data are not publicly available due to patient privacy restrictions.
De-identified derived data supporting the findings of this study are available
from the corresponding author upon reasonable request.

Code Availability

The custom MATLAB code used to compute dynamic wPLI features and perform
statistical analyses is available at Zenodo with the DOI:
10.5281/zenodo.XXXXXXX (to be replaced after archiving).

The repository contains all scripts, a README with execution instructions,
and a configuration file. Third-party functions (eegfilt, mafdr, pdist)
are publicly available through EEGLAB, the MATLAB Bioinformatics Toolbox,
and the MATLAB Statistics and Machine Learning Toolbox, respectively.
License

If you use this code in your research, please cite this paper.



