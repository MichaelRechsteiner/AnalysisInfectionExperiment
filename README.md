# AnalysisInfectionExperiment
R scripts and data to reproduce figures and analyses from the manuscript  "Rechsteiner M. et al. Host resistance and pathogen aggressiveness vary along an elevation gradient in a wild plant-pathogen interaction"\n

Script scripts/Figures.R will produce unedited figures 1C, 1D, 2, 3 and 4, as well as Supplementary Figures 1 and 2.\n
Script scripts/StatisticalAnalysis.R will produce unedited Table 1 and 2, as well as Supplementary Tables 1-14.\n

Execute scripts in your terminal from the repository root directory to swiftly reproduce tables and figures:\n
Rscript scripts/StatisticalAnalyses.R --outdir "$WhereYouWantTheTables"\n
Rscript scripts/Figures.R --outdir "$WhereYouWantTheFigures"\n

The file sessionInfo.txt contains a snapshot of the R evironment where the analyses were carried out.