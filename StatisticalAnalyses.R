#### Script to reproduce the statistical analyses performed in the manuscript "Patterns of local adaptation reverse along an elevation gradient in a wild plant-pathogen system" ####
#### required libraries/packages ####
require(tidyverse)
require(lme4)
require(ordinal)
require(glmmTMB)
require(report)
require(RVAideMemoire)
require(car)
require(emmeans)
require(argparser)

#### Set up argparser for command line to read in data and to set the output directory for tables####
## Insert the whole file paths under default in this script or add the file paths after the corresponding arguments
args <- arg_parser(description="Parser to analyse data for the manuscript 'Patterns of local adaptation reverse along elevation in a wild plant-pathogen interaction")

args <- add_argument(args,"--InfectionSeverityFile", default="InfectionSeverityData.csv", help="Infection severity file",type="character")
args <- add_argument(args,"--SporulationFile", default="SporulationTimeData.csv", help="Sporulation time file",type="character")
args <- add_argument(args,"--ChasmotheciaFile", default="ChasmotheciaData.csv", help="Chasmothecia file",type="character")
args <- add_argument(args,"--outdir", default="", help="Set output directory for result tables",type="character")
args <- add_argument(args,"--PlantAbundanceFile", default="PLLAAbundanceData.csv", help="Plant abundance file",type="character")
args <- add_argument(args,"--PathogenAbundanceFile", default="MildewAbundanceData.csv", help="Pathogen abundance file",type="character")

# Parse the command line arguments
pargs <- parse_args(args)

InfectionSeverityFile <- pargs[["InfectionSeverityFile"]]
SporulationFile <- pargs[["SporulationFile"]]
ChasmotheciaFile <- pargs[["ChasmotheciaFile"]]
outdir <- pargs[["outdir"]]
PlantAbundanceFile <- pargs[["PlantAbundanceFile"]]
PathogenAbundanceFile <- pargs[["PathogenAbundanceFile"]]
# read in the data files
## Infection status/Infection severity
Infection_Severity_Analysis <- read.csv(InfectionSeverityFile,stringsAsFactors = T)
## Sporulation time
Days_to_spore_Analysis <- read.csv(SporulationFile,stringsAsFactors = T)
## Presence of chasmothecia
Chasmothecia_Analysis <- read.csv(ChasmotheciaFile,stringsAsFactors = T)
##
PlantAbundance_Analysis <- read.csv(PlantAbundanceFile,stringsAsFactors = T)
##
PathogenAbundance_Analysis <- read.csv(PathogenAbundanceFile,stringsAsFactors = T)

JoinedAbundance <- merge(PlantAbundance_Analysis, PathogenAbundance_Analysis, by = c("Season","Population","Elevation","Site"))

#### STATISTICAL ANALYSES ####
#### Binomial regression: Infection status 0/1 across the two experimental blocks ####
Infection_Status_Analysis <- Infection_Severity_Analysis %>% 
  # Assign an infection status (0/1) for each Host-Pathogen pair based on infection severity reached in both experimental blocks.
  group_by(Plant.genotype,Pathogen.isolate,Host.population,Pathogen.population, Geography) %>%
  summarise(Bevan.mean = mean(Bevan.Score)) %>% 
  mutate(Infection = case_when(Bevan.mean >= 1 ~ 1,
                               Bevan.mean < 1 ~ 0)) %>%
  mutate(Host.population = fct_relevel(Host.population, c("Low","Mid-low","Mid-high","High")),
         Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High")))

# Binomial regression model of the infection status
Infection_Status_Model <- glmmTMB(data=Infection_Status_Analysis, Infection ~ Host.population*Pathogen.population + (1|Plant.genotype)+(1|Pathogen.isolate), family = binomial) # glmmTMB() was used instead of glmer() because of model convergence problems
Infection_Status_Model_Anova <- Anova(Infection_Status_Model)
Infection_Status_Model_report_table <-  report_table(Infection_Status_Model)

write.csv(Infection_Status_Model_Anova, paste0(outdir,"Infection_Status_AnovaTable.csv"))
write.csv(Infection_Status_Model_report_table, paste0(outdir,"Infection_Status_Model_report_table.csv"))

#### Ordinal regression: Infection severity of successful infections of both experimental blocks ####
Infection_Severity_for_model <- Infection_Severity_Analysis %>% 
# Prepare variables for analysis
  mutate(Bevan.Cat=as.factor(Bevan.Score),
         Block=as.factor(Block)) %>%
  filter(Bevan.Score >0) %>%
  mutate(Host.population = fct_relevel(Host.population, c("Low","Mid-low","Mid-high","High")),
         Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High")))

# Ordinal regression model of the ordered infection severity categories
Infection_Severity_Model <- clmm(data= Infection_Severity_for_model, Bevan.Cat ~ Block*Host.population*Pathogen.population + (1|Plant.genotype)+ (1|Pathogen.isolate))
Infection_Severity_Model_Anova <- RVAideMemoire::Anova.clmm(Infection_Severity_Model)
Infection_Severity_Model_emmeans_Pathogen <- emmeans(Infection_Severity_Model, pairwise ~ Pathogen.population)
Infection_Severity_Model_emmeans_Host <- emmeans(Infection_Severity_Model, pairwise ~ Host.population)
Infection_Severity_Model_emmeans_PathogenHost <- emmeans(Infection_Severity_Model, pairwise ~ Host.population:Pathogen.population)
Infection_Severity_Model_emmeans_HostBlock <- emmeans(Infection_Severity_Model, pairwise ~ Host.population:Block)

write.csv(Infection_Severity_Model_Anova, paste0(outdir,"Infection_Severity_AnovaTable.csv"))
write.csv(summary(Infection_Severity_Model_emmeans_Pathogen$contrasts, infer = T), paste0(outdir,"Infection_Severity_Model_emmeans_contrasts_Pathogen.csv"))
write.csv(summary(Infection_Severity_Model_emmeans_PathogenHost$contrasts, infer = T), paste0(outdir,"Infection_Severity_Model_emmeans_contrasts_PathogenHost.csv"))
write.csv(summary(Infection_Severity_Model_emmeans_Host$contrasts, infer = T), paste0(outdir,"Infection_Severity_Model_emmeans_contrasts_Host.csv"))
write.csv(summary(Infection_Severity_Model_emmeans_Pathogen$emmeans, infer = T), paste0(outdir,"Infection_Severity_Model_emmeans_Pathogen.csv"))
write.csv(summary(Infection_Severity_Model_emmeans_PathogenHost$emmeans, infer = T), paste0(outdir,"Infection_Severity_Model_emmeans_PathogenHost.csv"))
write.csv(summary(Infection_Severity_Model_emmeans_Host$emmeans, infer = T), paste0(outdir,"Infection_Severity_Model_emmeans_Host.csv"))
write.csv(summary(Infection_Severity_Model_emmeans_HostBlock$emmeans, infer = T), paste0(outdir,"Infection_Severity_Model_emmeans_HostBlock.csv"))
write.csv(summary(Infection_Severity_Model_emmeans_HostBlock$contrasts, infer = T), paste0(outdir,"Infection_Severity_Model_emmeans_contrasts_HostBlock.csv"))

# # Function to calculate the probabilities and odds from the model estimated log odds
# log_odds_to_prob_odds <- function(log_odds) {
#   # Ensure input is numeric
#   log_odds <- as.numeric(log_odds)
#   
#   # Calculate probabilities and odds
#   prob <- 1 / (1 + exp(-log_odds))
#   odds <- exp(log_odds)
#   return(c(prob, odds))
# }
# 
# # Extra section with code to calculate specific probabilities for different "cuts" of categories.
# # In the following section the difference in the probability of and infection being in category 4 or higher (only 5 in this case) 
# # between the low and high host populations is calculated
# emm <- emmeans(Infection_Severity_Model, ~ cut | Host.population, mode="cum.prob")
# probs <- as.data.frame(emm)
# Low_probs <- subset(probs, Host.population == "Low" & cut == "4|5")
# High_probs <- subset(probs, Host.population == "High" & cut == "4|5")
# prob_4_or_higher_pop1 <- 1 - Low_probs$cumprob
# prob_4_or_higher_pop4 <- 1 - High_probs$cumprob
# prob_diff <- prob_4_or_higher_pop1 - prob_4_or_higher_pop4
# prob_diff
# 
# # Here the difference in the probability of and infection being in category 4 or higher 
# # between low-elevation hosts infected with low-elevation pathogens and high-elevation host infected with high-elevation pathogens is calculated
# emm <- emmeans(Infection_Severity_Model, ~ cut | Host.population:Pathogen.population, mode="cum.prob")
# probs <- as.data.frame(emm)
# Low_probs <- subset(probs, Host.population == "Low" & Pathogen.population == "Low" & cut == "4|5")
# High_probs <- subset(probs, Host.population == "High" &Pathogen.population == "High"& cut == "4|5")
# prob_4_or_higher_pop1 <- 1 - Low_probs$cumprob
# prob_4_or_higher_pop4 <- 1 - High_probs$cumprob
# prob_diff <- prob_4_or_higher_pop1 - prob_4_or_higher_pop4
# prob_diff

#### Linear regression: Sporulation times of successful infections of both blocks ####
Days_to_spore_for_model <- Days_to_spore_Analysis %>%
  # Prepare variables for analysis
  mutate(Block=as.factor(Block))  %>%
  mutate(Host.population = fct_relevel(Host.population, c("Low","Mid-low","Mid-high","High")),
         Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High")))

Days_to_spore_model <- lmer(data = Days_to_spore_for_model, Sporulation.time ~ Host.population*Pathogen.population*Block + (1|Plant.genotype)+(1|Pathogen.isolate))
Days_to_spore_model_Anova <- Anova(Days_to_spore_model)
Days_to_spore_model_emmeans <- emmeans(Days_to_spore_model, pairwise ~ Host.population)
Days_to_spore_model_report_table <- report_table(Days_to_spore_model)

write.csv(Days_to_spore_model_Anova , paste0(outdir,"Days_to_spore_model_Anova.csv"))
write.csv(summary(Days_to_spore_model_emmeans$contrasts, infer = T) , paste0(outdir,"Days_to_spore_model_contrasts.csv"))
write.csv(summary(Days_to_spore_model_emmeans$emmeans, infer = T) , paste0(outdir,"Days_to_spore_model_emmeans.csv"))
write.csv(Days_to_spore_model_report_table, paste0(outdir,"Days_to_spore_model_report_table.csv"))

#### Binomial regression: Presence of chasmothecia at the day 15 of the first experimental block  ####
Chasmothecia_Analysis<- Chasmothecia_Analysis %>%
  mutate(Host.population = fct_relevel(Host.population, c("Low","Mid-low","Mid-high","High")),
         Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High")))
Chasmothecia_Model <- glmmTMB(data = Chasmothecia_Analysis, Chasmothecia ~ Host.population*Pathogen.population + (1|Plant.genotype) + (1|Pathogen.isolate), family = binomial)# glmmTMB() was used instead of glmer() because of model convergence problems
Chasmothecia_Model_Anova <- Anova(Chasmothecia_Model)
Chasmothecia_Model_report_table <- report_table(Chasmothecia_Model)

write.csv(Chasmothecia_Model_Anova, paste0(outdir,"Chasmothecia_Model_Anova.csv"))
write.csv(Chasmothecia_Model_report_table, paste0(outdir,"Chasmothecia_Model_report_table.csv"))
#### Linear regression: Field observational data ####
# lmer_fit_plantago <- lmer(Plantago.lanceolata~Elevation + (1|Season), data=PlantAbundance_Analysis)
# #### Fix this, two sites have the same elevation.... causes problems when merging
# lmer_fit_pathogen <- lmer(Infection.aggregates~Plantago.lanceolata+Elevation + (1|Season), data=JoinedAbundance)
# 
# 
# lmer_fit_pathogen2 <- JoinedAbundance %>% filter(Infected.plants < 900) %>% lmer(Infected.plants~Elevation + (1|Season), data=.)
# lmer_fit_pathogen3 <- JoinedAbundance %>% glm(Infected.plants~Elevation, data=., family = poisson)
# lmer_fit_pathogen4 <- JoinedAbundance %>% glmmTMB(Infected.plants ~ Plantago.lanceolata+Elevation + (1|Season), data=., family = nbinom2)
# 
# report_table_plantago <- report_table(lmer_fit_plantago)
# report_table_pathogen <- report_table(lmer_fit_pathogen)
# write.csv(report_table_plantago, paste0(outdir,"PlantagoAbundance_Model_report_table.csv"))
# write.csv(report_table_pathogen, paste0(outdir,"PathogenAbundance_Model_report_table.csv"))

# mean_count <- mean(JoinedAbundance$Infected.plants)
# variance_count <- var(JoinedAbundance$Infected.plants)
# 
# library(DHARMa)
# modres<-simulateResiduals(lmer_fit_pathogen4 )
# plot(modres)
# 
# testDispersion(modres)
# Anova(lmer_fit_pathogen4 )
# summary(lmer_fit_pathogen4)