# Run everything from the root directory.

################################################################################################
################                       Generate networks                        ################
################################################################################################

# generate network topology replicates
# from analysis_scripts folder
Rscript scripts/RGenNets/generateERNetTopos.R
Rscript scripts/RGenNets/generateBANetTopos.R
Rscript scripts/RGenNets/generateWSNetTopos.R

# Simulations on the cluster

################################################################################################
################                       Data preparation 1                       ################
################################################################################################

###### calculate network metrics
# in /results
bash ../scripts/sh_scripts/getNetMetrics.sh 20210609_ER_40node_0.05dens 2000;

###### summarize GPF data from evolved populations
# in /results
# summarize the results of all replicates per net sample
bash ../scripts/sh_scripts/calculateSelStrength.sh 20210609_ER_40node_0.05dens 2000;

###### combined network metrics and summarized GPF
# in /results
Rscript ../scripts/analysis/combineNetAndEvol.R 20210609_ER_40node_0.05dens ER 0.05 2000;



################################################################################################
################                       Data preparation 2                       ################
################################################################################################

## Prepare data
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_generate_data/0_prepareData.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_generate_data/0_prepareData_networkProperties.Rmd", "html_document")'

########################################     end     ###########################################
