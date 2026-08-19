<p align="center"> <img src="assets/COMPASS_logo_dark.png" alt="COMPASS Logo" width="400" height="400"> </p>

This software was developed in the framework of the cams2\_35\_bis (ECWMF) by [Athanasios Tsikerdekis (KNMI)](mailto:thanos.tsikerdekis@knmi.nl) to compare experiments of IFS-COMPO. Creates mean field maps for the specified period and timeseries per aerosol species for two experiments. Experiments can either be using the AER or the HAMM7 scheme.

## 📜 Attribution

[COMPosition Analysis Software for Simulations (COMPASS)](https://github.com/atsikerdekis/COMPASS)

## 🛠 Environment

1. Install Miniforge3:  
   🔗 https://github.com/conda-forge/miniforge

2. Create the environment using mamba:  
   `mamba env create -f environment/COMPASS.yml`

3. Activate the environment:  
   `mamba activate COMPASS`

## ▶️ Running

1. Change parameters under the EXPIREMENTS settings in the config.R file

2. Set runtype = "download", submit and wait until all file have been downloaded:
`Rscript 00.submit.R`

3. Set runtype = "plot" and submit again to create figures
`Rscript 00.submit.R`
