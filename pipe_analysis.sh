# Run everything from the root directory.

################################################################################################
################                           Main results                         ################
################################################################################################
# Each script creates a subfolder with the same name in the /plots folder, e.g.
# running mainResults_modelValidation.Rmd creates a subfolder called mainResults_modelValidation
# in the plots folder. The subfolders will hold the outputs of each script (figures, tables, etc.)

Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_modelValidation.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_evolPhenotypes.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_evolGenotypes.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_networkProperties.Rmd", "html_document")'



################################################################################################
################                              SI                                ################
################################################################################################

### Realization robustness
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/S_1.Rmd", "html_document")'
## Check colinearity of network metrics because they are explanatory variables
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/S_checkColinearityNetMetrics.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/S_FDR.Rmd", "html_document")'
### SI - Other topos
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_modelValidation_topo_BA.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_evolPhenotypes_topo_BA.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_evolGenotypes_topo_BA.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_modelValidation_topo_WS.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_evolPhenotypes_topo_WS.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_evolGenotypes_topo_WS.Rmd", "html_document")'
### SI - Filtering
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_modelValidation_filtered_1.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_evolPhenotypes_filtered_1.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_evolGenotypes_filtered_1.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_networkProperties_filtered_1.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_modelValidation_filtered_2.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_evolPhenotypes_filtered_2.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_evolGenotypes_filtered_2.Rmd", "html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("scripts_analysis/mainResults_networkProperties_filtered_2.Rmd", "html_document")'

########################################     end     ###########################################
