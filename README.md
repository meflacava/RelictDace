# RelictDace

Scripts and data for relict dace project (samples collected 2019-2022, RADseq data generated in 2022)


Sequence data have been deposited with NCBI into the Sequence Read Archive under BioProject PRJNA933690.


## Workflow documents:

BioinformaticsWorkflow.Rmd = initial bioinformatics processing and filtering steps

Plot_PCAngsd.Rmd = Run ANGSD for SNP discovery and filtering; produce PCA plots

Plot_NGSadmix.Rmd = Run NGSadmix on beagle file output from ANGSD; produce admixture plots

Plot_ThetaDiversity.Rmd = generate SFS and calculate diversity stats; produce boxplot of theta pi

Plot_Fst.Rmd = calculate pairwise Fst values for a given list of populations (for which you have already calculated SFS); produce Fst tile plot

RelictDace_ContractFigs.R = figures for contract report



## Folder contents:

Bamlists = bamlists and population lists

Metadata = sampling metadata, library prep data (e.g., barcodes), spatial data for mapping

Results = results files for plotting

Scripts = scripts referenced in R Markdown workflows
