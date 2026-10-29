# AnalysisInfectionExperiment
R scripts and data to reproduce figures and analyses from the manuscript  "Rechsteiner M. et al. Host resistance and pathogen aggressiveness vary along an elevation gradient in a wild plant-pathogen interaction

Script Figures.R will produce unedited figures 1C, 1D, 2, 3 and 4, as well as Supplementary Figures 1 and 2.
Script StatisticalAnalysis.R will produce unedited Table 1 and 2, as well as Supplementary Tables 1-14.

Data files:
PLLAAbundanceData.csv
MildewAbundanceData.csv
InfectionSeverityData.csv
SporulationTimeData.csv
ChasmotheciaData.csv

Execute scripts in your terminal from the directory where the Data files and Scripts are stored to swiftly reproduce tables and figures:
Rscript StatisticalAnalyses.R --outdir "$WhereYouWantTheTables"
Rscript Figures.R --outdir "$WhereYouWantTheFigures"

