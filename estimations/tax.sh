# setup the password as an environmental variable
source ~/.bash-profile

# just to check the progress of simulations on spartan
sshpass -e ssh twozniak@spartan.hpc.unimelb.edu.au
squeue -u twozniak
exit

# Upload files
sshpass -e scp JOE_RR_sim/spartan/bsvar*.* twozniak@spartan.hpc.unimelb.edu.au:/data/projects/punim0093/Carnaval-gulasz/JOE_RR_sim/
sshpass -e scp JOE_RR_sim/spartan/dgps/*.csv twozniak@spartan.hpc.unimelb.edu.au:/data/projects/punim0093/Carnaval-gulasz/JOE_RR_sim/dgps/
sshpass -e scp JOE_RR_sim/spartan/dgps/*.rda twozniak@spartan.hpc.unimelb.edu.au:/data/projects/punim0093/Carnaval-gulasz/JOE_RR_sim/dgps/
sshpass -e scp JOE_RR_sim/spartan/src_cpp/*.cpp twozniak@spartan.hpc.unimelb.edu.au:/data/projects/punim0093/Carnaval-gulasz/JOE_RR_sim/src_cpp/

sshpass -e scp spartan/*_ce.* twozniak@spartan.hpc.unimelb.edu.au:/data/projects/punim0093/Carnaval-gulasz/JOE_RR_sim/
sshpass -e scp spartan/MR2006.RData twozniak@spartan.hpc.unimelb.edu.au:/data/projects/punim0093/Carnaval-gulasz/JOE_RR_sim/
sshpass -e scp spartan/source_code/bsvar_sv_ce.cpp twozniak@spartan.hpc.unimelb.edu.au:/data/projects/punim0093/Carnaval-gulasz/JOE_RR_sim/src_cpp/

# Download files
sshpass -e scp twozniak@spartan.hpc.unimelb.edu.au:/data/projects/punim0093/Carnaval-gulasz/JOE_RR_sim/results/sv_780_0001.rda JOE_RR_sim/spartan/results/
scp twozniak@spartan.hpc.unimelb.edu.au:/data/projects/punim0093/Carnaval-gulasz/JOE_RR_sim/results/*.rda JOE_RR_sim/spartan/results/

scp twozniak@spartan.hpc.unimelb.edu.au:/data/projects/punim0093/Carnaval-gulasz/JOE_RR_sim/results/tax23_ce.rda spartan/results/
scp twozniak@spartan.hpc.unimelb.edu.au:/data/projects/punim0093/Carnaval-gulasz/JOE_RR_sim/results/tax06_ce.rda spartan/results/
scp twozniak@spartan.hpc.unimelb.edu.au:/data/projects/punim0093/Carnaval-gulasz/JOE_RR_sim/results/taxMR_ce.rda spartan/results/

for i in $(seq -w 1 130); do
  sshpass -e scp twozniak@spartan.hpc.unimelb.edu.au:/data/projects/punim0093/Carnaval-gulasz/JOE_RR_sim/results/msh_78*$i.rda JOE_RR_sim/spartan/results/
done

say I am finished

# working with svar_betel on spartan
#################################################
sshpass -e ssh twozniak@spartan.hpc.unimelb.edu.au
cd /data/projects/punim0093/Carnaval-gulasz/JOE_RR_sim/

sbatch tax23_ce.slurm
sbatch tax06_ce.slurm
sbatch taxMR_ce.slurm

sbatch bsvar_sv_sv.slurm
sbatch bsvar_sv_ga.slurm
sbatch bsvar_sv_ms.slurm
sbatch bsvar_sv_chan_sv.slurm
sbatch bsvar_sv_chan_ga.slurm 
sbatch bsvar_sv_chan_ms.slurm 
sbatch bsvar_msh_sv.slurm 
sbatch bsvar_msh_ga.slurm 
sbatch bsvar_msh_ms.slurm 
sbatch bsvar_ce_sv.slurm
sbatch bsvar_ce_ga.slurm
sbatch bsvar_ce_ms.slurm
squeue -u twozniak

sbatch bsvars_sv_sv.slurm
sbatch bsvars_sv_ga.slurm
sbatch bsvars_sv_ms.slurm
sbatch bsvars_sv_chan_sv.slurm
sbatch bsvars_sv_chan_ga.slurm 
sbatch bsvars_sv_chan_ms.slurm 
sbatch bsvars_msh_sv.slurm 
sbatch bsvars_msh_ga.slurm 
sbatch bsvars_msh_ms.slurm 
sbatch bsvars_ce_sv.slurm
sbatch bsvars_ce_ga.slurm
sbatch bsvars_ce_ms.slurm

sbatch bsvar_spartan_wrap.slurm

squeue -u twozniak

rm *.out
cat *9500.out
cat tax23_ce.R
cat *spartan_wrap.R
tail -200 src_cpp/bsvar_sv_ce.cpp
ls
ls results/
# rm results/*.rda
exit
#

# install bsvars package on spartan
#################################################
sshpass -e ssh twozniak@spartan.hpc.unimelb.edu.au
sinteractive
module load R/4.2.1
R
devtools::install_git("https://github.com/bsvars/bsvars.git")
bsvars::bsvar
q("no")
exit

