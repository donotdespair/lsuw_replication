# lsuw_replication
Replication files for Lütkepohl, Shang, Uzeda, Woźniak (2024) Partial Identification of Heteroskedastic Structural VARs: Theory and Bayesian Inference, DOI: https://doi.org/10.48550/arXiv.2404.11057


## Replications of results from Section 7. Empirical application: identification of tax shocks

1. Set working directory to `/estimations/`. Install the **bsvars** package. Scripts from the `/estimations/source_code/` are necessary to run the estimations.
2. Files `tax*.R` are the **R** scripts and files `tax*.slurm` contain **slurm** files that run the estimations on the HPC server. The results will be saved in the `/estimations/results/` folder.
3. Run the normalisation of the estimation results by executing the following **R** scripts `tax*_our.R`.
4. Obtain the particular figures and tables by running scripts:
```
DATA.R
HV_ce.R
HV_MR.R
HV_PM.R
HV_PMall.R
HV.R
IRF_NORM.R
IRF_PM.R
IRF.R
OMEGA_MR.R
OMEGA_PM.R
OMEGA.R
PAR.R
SDDR.R
SHOCK.R
tax_get_sddr_table.R
tax_get_sddr.R
```
5. Obtain the reproductions of the results from other papers by running scripts from folder `/estimations/BP_reproductions/`.

## Replications of results from Section 6. A Monte Carlo study

1. Set working directory to `/simulations/`. Install the **bsvars** package. Scripts from the `/simulations/src_cpp/` are necessary to run the simulations.
2. Run the `/simulations/dgps/dgps260.R` and `/simulations/dgps/dgps780.R` scripts to generate the data for the Monte Carlo simulations. The data will be saved in the `/simulations/dgps/` folder and are necessary to run subsequent files and steps of the simulation.
3. Files 
```
bsvar_ce.R
bsvar_msh.R
bsvar_sv.R
bsvar_sv_chan.R
bsvars_ce.R
bsvars_msh.R
bsvars_sv_chan.R
bsvars_sv.R
```
contain **R** scripts that run the simulations and save the results in the `/simulations/results/` folder.
4. Run the simulations by executing the following `*.slurm` files on your HPC server. These files contain a total of 120 days of computations, so parallel computing is necessary.
```
bsvar_ce_ga.slurm
bsvar_ce_ms.slurm
bsvar_ce_sv.slurm
bsvar_msh_ga.slurm
bsvar_msh_ms.slurm
bsvar_msh_sv.slurm
bsvar_sv_chan_ga.slurm
bsvar_sv_chan_ms.slurm
bsvar_sv_chan_sv.slurm
bsvar_sv_ga.slurm
bsvar_sv_ms.slurm
bsvar_sv_sv.slurm

bsvars_ce_ga.slurm
bsvars_ce_ms.slurm
bsvars_ce_sv.slurm
bsvars_msh_ga.slurm
bsvars_msh_ms.slurm
bsvars_msh_sv.slurm
bsvars_sv_chan_ga.slurm
bsvars_sv_chan_ms.slurm
bsvars_sv_chan_sv.slurm
bsvars_sv_ga.slurm
bsvars_sv_ms.slurm
bsvars_sv_sv.slurm
```
The results will be saved in the `/simulations/results/` folder.
5. Run files from the `/simulations/` folder on your HPC server to collect the results from the simulations. The results will be saved in the `/simulations/results/` folder.
```
bsvar_sddr_wrap.R
bsvars_sddr_wrap.R
bsvar_mse_wrap.R
bsvar_spartan_wrap.R
bsvar_sv_wrap.R
```
    
