
# bsvars package installed from the developer's repo using:
# devtools::install_git("https://github.com/bsvars/bsvars.git")

files_in_results  = list.files("results/")
files_sddr        = files_in_results[grepl("_", files_in_results)]
I                 = length(files_sddr)

for (i in 1:I) {
  load(paste0("results/", files_sddr[i]))
  save(sddr, file = paste0("results/sddr_", files_sddr[i]))
}

