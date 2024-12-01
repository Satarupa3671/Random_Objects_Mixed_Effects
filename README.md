# Geodesic Mixed Effects Models for Repeatedly Observed/Longitudinal Random Objects

This repository contains the implementation for the paper ["Geodesic Mixed Effects Model for Repeatedly Observed/Longitudinal Random Objects"](https://arxiv.org/pdf/2307.05726) in R.

## Summary of the paper

Mixed effect modeling for longitudinal data is challenging when the observed data are random objects, which are complex data taking values in a general metric space without either global linear or local linear (Riemannian) structure. In such settings the classical additive error model and distributional assumptions are unattainable. Due to the rapid advancement of technology, longitudinal data containing complex random objects, such as covariance matrices, data on Riemannian manifolds, and probability distributions are becoming more common. Addressing this challenge, we develop a mixed-effects regression for data in geodesic spaces, where the underlying mean response trajectories are geodesics in the metric space and the deviations of the observations from the model are quantified by perturbation maps or transports. A key finding is that the geodesic trajectories assumption for the case of random objects is a natural extension of the linearity assumption in the standard Euclidean scenario to the case of general geodesic metric spaces. Geodesics can be recovered from noisy observations by exploiting a connection between the geodesic path and the path obtained by global Fréchet regression for random objects. The effect of baseline Euclidean covariates on the geodesic paths is modeled by another Fréchet regression step. We study the asymptotic convergence of the proposed estimates and provide illustrations through simulations and real-data applications.

## File Overview

-   `code/` implementation for the paper [""Geodesic Mixed Effects Model for Repeatedly Observed/Longitudinal Random Objects"](https://arxiv.org/pdf/2307.05726) in R.
    -   `data_adni/` Folder containing preprocessed data, analysis code, and outputs (visualizations and tables) generated in Section 5.1 and continued on Section S.4.3 in the Supplement.

        -   `data/` For the ADNI data (adni.loni.usc.edu), one has to create an account at adni.loni.usc.edu and, after logging in, need to follow the steps at <https://adni.loni.usc.edu/data-samples/access-data/>. The rs-fMRI data was downloaded, accessed, and preprocessed by the UC Davis group led by Dr. Hans Mueller in 2020. The raw fMRI image data is above 100 GB and cannot be shared without ADNI permission. An outline of the preprocessing steps from the raw fMRI images using the necessary tools and software is now attached in the folder code/data_adni/data/adni_preprocessing.pdf. We used this preprocessed data in the application and supplied it here in the folder code/data_adni/data as .Rda files for Cognitive Normal (CN) and Mild Cognitive Impairment (MCI) groups of subjects. The relevant covariates are also stored as separate .Rda files in the same location. This is described as follows.
            -   resp_final_06262022_AAL_CN.Rda– mxm correlation matrix data obtained from the fMRI time series following the AAL parcellation for Cognitive Normal (CN) subjects.
            -   resp_final_06262022_AAL_MCI.Rda – m\*m correlation matrix data obtained from the fMRI time series following the AAL parcellation for Mild Cognitive Impairment (MCI) subjects.
            -   covariate_final_06262022_AAL_CN.Rda – data on age, sex, and total cognitive score (C-Score) that are used as baseline covariates for Cognitive Normal (CN) subjects.
            -   covariate_final_06262022_AAL_MCI.Rda – data on age, sex, and total cognitive score (C-Score) that are used as baseline covariates for Mild Cognitive Impairment (MCI) subjects.
        -   `code/` The analysis and visulations included in the paper is generated in the file analysis.R
        -   `output/` The figures and tables included in the paper are saved in this folder.

    -   `data_mort/` Folder containing preprocessed data, analysis code, and outputs (visualizations and tables) generated in Section S.4.2 in the Supplement.

        -   `data/` For the mortality data, the life tables for males and females are available on the Human Mortality Database (www.mortality.org) according to the yearly age group for ages 0 to 110. We selected life tables for females from the year 1990-2019 from the period life table. The data on the relevant covariates for the same countries over the same years is obtained from the World Bank Database at <https://data.worldbank.org>. We selected 28 countries for which data is available on both the Human Mortality Database and the World Bank Database for the 30 years 1990-2019. The pre-processing code from the life table data (downloaded from the above two sources and supplied as .csv files in the folder code/data_mort/data/ ) is given in code/data_mort/code/preprocess.R. In this pre-processing step, the remaining life distribution and its density are obtained from the available lifetable data that correspond to histograms with a bin width of one year by adding a smoothing step, for which we used the R package “frechet” with bandwidth 2 years. This is described as follows.

        1.  Raw Data – Downloaded from World Bank Database at <https://data.worldbank.org>
            -   Fertility/API_SP.DYN.TFRT.IN_DS2_EN_csv_v2_3404027.csv – Data on Fertility rate, total (births per woman)
            -   GDP/API_NY.GDP.PCAP.PP.CD_DS2_en_csv_v2_3401652.csv – Data on GDP per capita based on purchasing power parity
            -   Pop_Growth/API_SP.POP.GROW_DS2_en_csv_v2_3404396.csv – Data on Population Growth (annual percentage)
            -   Unemployment/API_SL.UEM.TOTL.ZS_DS2_en_csv_v2_3401869.csv – Data on Unemployment, total (% of total labor force) (modeled ILO estimate)
        2.  Raw Data – Downloaded from Human Mortality Database at www.mortality.org
            -   data/fltper_1x1 – The folder contains 50 .txt files for the data on Period life tables for the female population indexed by calendar year for each country. Each file has the life table data on Deaths,3 exposure-to-risk, death rates, and life tables are given in similar formats of age and time 1x1 (by age and year).
        3.  Processed Data – Using the code provided in code/data_mort/code/preprocess.R
            -   lt_subset.csv – data on the life tables for females in 28 countries over 30 calendar years from the human mortality database.
            -   covariate_1990.csv – data on the relevant covariates for 28 countries measured in the year 1990 as the baseline, obtained from the World Bank database.

        -   `code/` The preprocessing steps from the raw data is done in the file preprocess.R. The analysis and visulations included in the paper is generated in the file analysis.R
        -   `output/` The figures and tables included in the paper are saved in this folder.

    -   `simulation/` Folder containing preprocessed data, analysis code, and outputs (visualizations and tables) generated in Section 4 and continued on Sections S.4.1, and S.4.4 in the Supplement.

        -   `data/` For the simulation studies, the data used is generated using the mechanisms described in the respective sections in the manuscripts and supplementary materials. For example, lines 7-76 on the file code/simulation/code/analysis_density.R implements the steps involved in the four different data generation mechanisms for univariate distributional data in Section 4 of the main manuscript; lines 14-72 describe the same for data on the surface of the sphere discussed in Section S.4.1. Furthermore, for the comparison with classical Euclidean approaches, we use the dataset on the study of the influence on Menarche on changes in body fat accretion from Fritzmaurice et al. (2012), which is provided in the file code/simulation/data/fat.dta. This is described as follows:
            -   code/simulation/data/fta.dta -- dataset on the study of the influence of Menarche on changes in body fat accretion from Fritzmaurice et al. (2012).
            -   code/simulation/data/MISE_density_boxplot_comp_dense_nA_settingB.Rda – Mean Integrated Squared Error (Wasserstein Distance) between the observed and estimated distributional objects for dense design (50 observations per subject along the geodesic); A = 50, 400, 1000; B = 1, 2, 3, 4. (So in all, there are 12 files for each combination of A and B).
            -   code/simulation/data/MISE_density_boxplot_comp_sparse_nA_settingB.Rda – Mean Integrated Squared Error (Wasserstein Distance) between the observed and estimated distributional objects for sparse design (50 observations per subject along the geodesic); A = 50, 400, 1000; B = 1, 2, 3, 4. (So in all, there are 12 files for each combination of A and B).
            -   code/simulation/data/subj_sp_over_alpha_Ak_Setting4.Rda – the true, observed (perturbed with parameter alpha), and estimated distributional object responses as densities for a randomly selected simulation sample generated under setting IV with a sparse design where each subject has 2 to 5 repeated measurements, comparing varying perturbation levels 0.01,0.1,0.3, corresponding to the indices k = 1, 2, 3. (So in all, there are 4 files for each of the three alpha levels).
            -   code/simulation/data/subj_specific_dens_zz_settingB.Rda – The time-dynamic effect of the baseline covariate for distributional objects represented as densities for a randomly selected simulation sample. The true and estimated densities at time points t = 0, 0.5, and 1 along the geodesics for simulation settings I-IV are saved as .Rda. Data were generated under a sparse design, where each subject has 2 to 5 repeated measurements and where response distributions were perturbed with a fixed small perturbation level α = 0.1. The estimated/predicted densities are computed for the 10%, 50%, and 90% quantile levels of the baseline covariate Z; B = 1, 2, 3, 4. (So in all, there are 4 files for each of the four settings).
        -   `code/` The relevant files are as follows: \* For univariate distributions in Section 4, see analysis_density.R, \* For data on the surface of a sphere in Section S.4.1, see analysis_spehere.R \* For comparison with Euclidean responses data in Section S.4.4, see analysis_Eucl.R
        -   `output/` The figures and tables included in the paper are saved in this folder.

## Citation

Please cite our paper ["Geodesic Mixed Effects Model for Repeatedly Observed/Longitudinal Random Objects"](https://arxiv.org/pdf/2307.05726).

```         
@article{bhattacharjee2023geodesic,
  title={Geodesic mixed effects models for repeatedly observed/longitudinal random objects},
  author={Bhattacharjee, Satarupa and M{\"u}ller, Hans-Georg},
  journal={arXiv preprint arXiv:2307.05726},
  year={2023}
}
```
