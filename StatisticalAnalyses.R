#### Script to reproduce the statistical analyses performed in the manuscript "Rechsteiner M. et al. Host resistance and pathogen aggressiveness vary along an elevation gradient in a wild plant-pathogen system" ####
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
require(DHARMa)
require(stringr)

#### Set up argparser for command line to read in data and to set the output directory for tables####
## Insert the whole file paths under default in this script or add the file paths after the corresponding arguments
args <- arg_parser(description="Parser to analyse data for the manuscript 'Host resistance and pathogen aggressiveness vary along an elevation gradient in a wild plant-pathogen system'")
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
JoinedAbundance$Population <- factor(JoinedAbundance$Population, levels = c("High","Mid-high","Mid-low","Low"))
JoinedAbundance$Season <- factor(JoinedAbundance$Season, levels = c("2019","2020","2021"))

#### STATISTICAL ANALYSES ####
#### Linear and binomial regression: Field observational data, Supplementary Table 1-5  ####
## Linear model of Host Density
JoinedAbundance$HostDensity <- JoinedAbundance$Plantago.lanceolata/36
JoinedAbundance$Hosts_per_2500m2 <- round(JoinedAbundance$HostDensity*2500)
JoinedAbundance$Healthy.plants <- JoinedAbundance$Hosts_per_2500m2 - JoinedAbundance$Infected.plants

lmm_HostDensity <- lmer(HostDensity ~ Population*Season+ (1 | Site), data = JoinedAbundance)

Anova_hostdensity <- car::Anova(lmm_HostDensity) # Season and Season:Population 

Emmeans_Population <- emmeans(lmm_HostDensity, ~ Population)
Emmeans_Season <- emmeans(lmm_HostDensity, pairwise ~ Season)$contrasts
Emmeans_Season_within_Pops <- emmeans(lmm_HostDensity, pairwise ~ Season | Population)$contrasts

write.csv(Anova_hostdensity, paste0(outdir,"/Supplementary_Table1.csv")) # Supplementary Table 1
write.csv(Emmeans_Season , paste0(outdir,"/Supplementary_Table2.csv"))# Supplementary Table 2
write.csv(Emmeans_Season_within_Pops , paste0(outdir,"/Supplementary_Table3.csv"))# Supplemetary Table 3

## Grouped Binomial Model of infeciton prevalence
grouped_binomial_prevalence <- JoinedAbundance %>% glmmTMB(cbind(Infected.plants, Healthy.plants) ~ Population+Season+ (1 | Site), data = ., family = binomial)

GroupedBinom_anova <- car::Anova(grouped_binomial_prevalence)
ContrastsGroupedBinom_Population <- emmeans(grouped_binomial_prevalence, pairwise ~ Population, type = "response")$contrast
ContrastsGroupedBinom_Season <- emmeans(grouped_binomial_prevalence, pairwise ~ Season, type = "response")$contrast

write.csv(GroupedBinom_anova, paste0(outdir,"/Supplementary_Table4.csv"))# Supplemetary Table 4
write.csv(ContrastsGroupedBinom_Population, paste0(outdir,"/Supplementary_Table5.1.csv"))# Supplemetary Table 5
write.csv(ContrastsGroupedBinom_Season, paste0(outdir,"/Supplementary_Table5.2.csv"))# Supplemetary Table 5

#### Binomial regression: Infection status 0/1 across the two experimental blocks, Table 1 and Supplementary Table 6 ####
Infection_Status_Analysis <- Infection_Severity_Analysis %>% 
  # Assign an infection status (0/1) for each Host-Pathogen pair based on infection severity reached in both experimental blocks.
  group_by(Plant.genotype,Pathogen.isolate,Host.population,Pathogen.population, Geography) %>%
  summarise(Bevan.mean = mean(Bevan.Score)) %>% 
  mutate(Infection = case_when(Bevan.mean >= 1 ~ 1,
                               Bevan.mean < 1 ~ 0)) %>%
  mutate(Host.population = fct_relevel(Host.population, c("Low","Mid-low","Mid-high","High")),
         Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High")))
Infection_Status_Analysis$PlantSite <- substr(Infection_Status_Analysis $Plant.genotype, 1, 2)
Infection_Status_Analysis$PlantSite <- as.factor(Infection_Status_Analysis $PlantSite)

# Binomial regression model of the infection status
Infection_Status_Model <- glmmTMB(data=Infection_Status_Analysis, Infection ~ Host.population*Pathogen.population + (1|Host.population:PlantSite:Plant.genotype)+(1|Pathogen.population:Pathogen.isolate), family = binomial) # glmmTMB() was used instead of glmer() because of model convergence problems

Infection_Status_Model_Anova <- Anova(Infection_Status_Model)
Infection_Status_Model_report_table <-  report_table(Infection_Status_Model)

write.csv(Infection_Status_Model_Anova, paste0(outdir,"/Table1.csv")) # Table 1
write.csv(Infection_Status_Model_report_table, paste0(outdir,"/Supplementary_Table6.csv")) # Supplementary Table 6

#### Ordinal regression: Infection severity of successful infections of both experimental blocks, Table 2, Supplementary Table 7-10  ####
Infection_Severity_for_model <- Infection_Severity_Analysis %>% 
# Prepare variables for analysis
  mutate(Bevan.Cat=as.factor(Bevan.Score),
         Block=as.factor(Block)) %>%
  filter(Bevan.Score > 0) %>%
  mutate(Host.population = fct_relevel(Host.population, c("Low","Mid-low","Mid-high","High")),
         Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High"))) %>%
  mutate(Bevan.final = case_when(
    Bevan.Cat == "0" ~ "0",
    Bevan.Cat == "1" ~ "1",
    Bevan.Cat == "1.5" ~ "2",
    Bevan.Cat == "2.5" ~ "3",
    Bevan.Cat == "3" ~ "4",
    Bevan.Cat == "4" ~ "5"
  )) %>%
  mutate(Bevan.final = factor(Bevan.final, levels = c("0","1","2","3","4","5")))

Infection_Severity_for_model$PlantSite <- substr(Infection_Severity_for_model$Plant.genotype, 1, 2)
Infection_Severity_for_model$PlantSite <- as.factor(Infection_Severity_for_model$PlantSite)

# Ordinal regression model of the ordered infection severity categories
Infection_Severity_Model <- clmm(data= Infection_Severity_for_model, Bevan.final ~ Block*Host.population*Pathogen.population + (1|Host.population:PlantSite:Plant.genotype)+ (1|Pathogen.population:Pathogen.isolate))

# Run anovas for the infection severity models
Infection_Severity_Model_Anova <- RVAideMemoire::Anova.clmm(Infection_Severity_Model)

# Function to transform the log of the odds estimate into an odds ratio
calculate_OR <- function(emmeans_df) {
  emmeans_df$OddsRatio    <- exp(emmeans_df$estimate)
  emmeans_df$InvOddsRatio <- 1 / emmeans_df$OddsRatio
  emmeans_df$OR_lowerCI   <- exp(emmeans_df$estimate - 1.96 * emmeans_df$SE)
  emmeans_df$OR_upperCI   <- exp(emmeans_df$estimate + 1.96 * emmeans_df$SE)
  return(emmeans_df)
}

Infection_Severity_Model_contrasts_Host <- as.data.frame(emmeans(Infection_Severity_Model, pairwise ~ Host.population,type = "response")$contrast)
Infection_Severity_Model_contrasts_Host <- calculate_OR(Infection_Severity_Model_contrasts_Host) 

Infection_Severity_Model_contrasts_Pathogen <- as.data.frame(emmeans(Infection_Severity_Model, pairwise ~ Pathogen.population, type = "response")$contrast)
Infection_Severity_Model_contrasts_Pathogen <- calculate_OR(Infection_Severity_Model_contrasts_Pathogen)

Infection_Severity_Model_contrasts_PathogenHost <- as.data.frame(emmeans(Infection_Severity_Model, pairwise ~ Host.population|Pathogen.population, type = "response")$contrast)
Infection_Severity_Model_contrasts_PathogenHost <- calculate_OR(Infection_Severity_Model_contrasts_PathogenHost)

Infection_Severity_Model_contrasts_HostPathogen <- as.data.frame(emmeans(Infection_Severity_Model, pairwise ~ Pathogen.population|Host.population, type = "response")$contrast)
Infection_Severity_Model_contrasts_HostPathogen <- calculate_OR(Infection_Severity_Model_contrasts_HostPathogen)

Infection_Severity_Model_contrasts_HostBlock <- as.data.frame(emmeans(Infection_Severity_Model, pairwise ~ Block | Host.population, type = "response")$contrast)
Infection_Severity_Model_contrasts_HostBlock <- calculate_OR(Infection_Severity_Model_contrasts_HostBlock)

Infection_Severity_Model_emmeans_PathogenHost <- as.data.frame(emmeans(Infection_Severity_Model, pairwise ~ Pathogen.population|Host.population, type = "response",mode = "mean.class")$emmeans)

write.csv(Infection_Severity_Model_Anova, paste0(outdir,"/Table2.csv")) # Table 2
write.csv(Infection_Severity_Model_contrasts_Pathogen, paste0(outdir,"/Supplementary_Table7.1.csv")) # Supplementary Table 7
write.csv(Infection_Severity_Model_contrasts_PathogenHost , paste0(outdir,"/Supplementary_Table8.1.csv"))# Supplementary Table 8
write.csv(Infection_Severity_Model_contrasts_Host, paste0(outdir,"/Supplementary_Table7.2.csv")) # Supplemetary Table 7
write.csv(Infection_Severity_Model_contrasts_HostPathogen, paste0(outdir,"/Supplementary_Table8.2.csv")) # Supplementary Table 8
write.csv(Infection_Severity_Model_emmeans_PathogenHost, paste0(outdir,"/Supplementary_Table9.csv"))# Supplementary Table 9
write.csv(Infection_Severity_Model_contrasts_HostBlock, paste0(outdir,"/Supplementary_Table10.csv"))# Supplementary Table 10

## EXTRA SECTION 
# Function to calculate the probabilities and odds from the model estimated log odds
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
# prob_4_or_higher_pop1 <- 1 - Low_probs$cum.prob
# prob_4_or_higher_pop4 <- 1 - High_probs$cum.prob
# prob_diff <- prob_4_or_higher_pop1 - prob_4_or_higher_pop4
# prob_diff
# 
# # Here the difference in the probability of an infection being in category 4 or higher 
# # between low-elevation hosts infected with low-elevation pathogens and high-elevation host infected with high-elevation pathogens is calculated
# emm <- emmeans(Infection_Severity_Model, ~ cut | Host.population:Pathogen.population, mode="cum.prob")
# probs <- as.data.frame(emm)
# Low_probs <- subset(probs, Host.population == "Low" & Pathogen.population == "Low" & cut == "4|5")
# High_probs <- subset(probs, Host.population == "High" &Pathogen.population == "High"& cut == "4|5")
# prob_4_or_higher_pop1 <- 1 - Low_probs$cumprob
# prob_4_or_higher_pop4 <- 1 - High_probs$cumprob
# prob_diff <- prob_4_or_higher_pop1 - prob_4_or_higher_pop4
# prob_diff

#### Linear regression: Sporulation times of successful infections of both blocks, Supplementary Table 11 and 12 ####
Days_to_spore_for_model <- Days_to_spore_Analysis %>%
  # Prepare variables for analysis
  mutate(Block=as.factor(Block),
         Sporulation.timeLog = log(Sporulation.time))  %>%
  mutate(Host.population = fct_relevel(Host.population, c("Low","Mid-low","Mid-high","High")),
         Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High")))
Days_to_spore_for_model$PlantSite <- substr(Days_to_spore_for_model$Plant.genotype, 1, 2)
Days_to_spore_for_model$PlantSite <- as.factor(Days_to_spore_for_model$PlantSite)

Days_to_spore_model <- lmer(data = Days_to_spore_for_model, Sporulation.time ~ Host.population*Pathogen.population*Block + (1|Host.population:PlantSite:Plant.genotype)+(1|Pathogen.population:Pathogen.isolate))

Days_to_spore_model_Anova <- Anova(Days_to_spore_model)
Days_to_spore_model_emmeans <- emmeans(Days_to_spore_model, pairwise ~ Host.population)

write.csv(Days_to_spore_model_Anova , paste0(outdir,"/Supplementary_Table11.csv")) # Supplementary Table 11
write.csv(summary(Days_to_spore_model_emmeans$contrasts, infer = T) , paste0(outdir,"/Supplementary_Table12.csv")) # Supplementary Table 12
#### Binomial regression: Presence of chasmothecia at the day 15 of the first experimental block, Supplementary Table 13  ####
Chasmothecia_Analysis<- Chasmothecia_Analysis %>%
  mutate(Host.population = fct_relevel(Host.population, c("Low","Mid-low","Mid-high","High")),
         Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High")))
Chasmothecia_Analysis$PlantSite <- substr(Chasmothecia_Analysis$Plant.genotype, 1, 2)
Chasmothecia_Analysis$PlantSite <- as.factor(Chasmothecia_Analysis$PlantSite)


Chasmothecia_Model <- glmmTMB(data = Chasmothecia_Analysis, Chasmothecia ~ Host.population*Pathogen.population + (1|Host.population:PlantSite:Plant.genotype)+(1|Pathogen.population:Pathogen.isolate), family = binomial)# glmmTMB() was used instead of glmer() because of model convergence problems

Chasmothecia_Model_Anova <- car::Anova(Chasmothecia_Model)

write.csv(Chasmothecia_Model_Anova, paste0(outdir,"/Supplementary_Table13.csv")) # Supplementary Table 13

#### Calculate Pseudo R2 values to compare model predictors of the ordinal regression model####
# Check the relative explained variance

null <- clmm(Bevan.Cat ~ 1 + (1|Plant.genotype)+ (1|Pathogen.isolate), data = Infection_Severity_for_model, Hess = TRUE)

# M0: block only (nuisance backbone)
M0 <- clmm(Bevan.Cat ~ Block + (1|Plant.genotype)+ (1|Pathogen.isolate), data = Infection_Severity_for_model, Hess = TRUE)

# M1: add host + pathogen main effects
M1 <- clmm(Bevan.Cat ~ Block + Host.population + Pathogen.population + (1|Plant.genotype)+ (1|Pathogen.isolate),
           data = Infection_Severity_for_model, Hess = TRUE)

# M2: add block's own interactions with host/pathogen (nuisance, not your focus,
#     but needed before you can cleanly isolate Host:Pathogen)
M2 <- clmm(Bevan.Cat ~ Block + Host.population + Pathogen.population +
             Block:Host.population + Block:Pathogen.population + (1|Plant.genotype)+ (1|Pathogen.isolate),
           data = Infection_Severity_for_model, Hess = TRUE)

# M3: add the interaction of interest
M3 <- clmm(Bevan.Cat ~ Block + Host.population + Pathogen.population +
                  Block:Host.population + Block:Pathogen.population+
                  Host.population:Pathogen.population + (1|Plant.genotype)+ (1|Pathogen.isolate),
                data = Infection_Severity_for_model, Hess = TRUE)

# M4 = full model: add the three-way to complete Block*Host.population*Pathogen.population
M4 <- clmm(Bevan.Cat ~ Block*Host.population*Pathogen.population + (1|Plant.genotype)+ (1|Pathogen.isolate),
           data = Infection_Severity_for_model, Hess = TRUE)

R2 <- function(model, null) 1 - as.numeric(logLik(model)) / as.numeric(logLik(null))

R2_block            <- R2(M0, null)
R2_mains            <- R2(M1, null) - R2(M0, null)   # Host + Pathogen mains
R2_block_interax    <- R2(M2, null) - R2(M1, null)   # Block:Host + Block:Pathogen 
R2_HostPath_interax <- R2(M3, null) - R2(M2, null)   # <-- the comparison you want
R2_threeway         <- R2(M4, null) - R2(M3, null)   # Block:Host:Pathogen

R2_full <- R2(M4, null)  # should equal sum of all five increments above


anova(M0, M1)   # mains: significant?
anova(M2, M3)   # Host:Pathogen interaction, controlling for block-by-predictor interactions: significant?

## Combined pseudo-R2 / LRT summary table for the CLMM model sequence
## Run this AFTER fitting null, M0, M1, M2, M3, M4 as in your existing code

# Helper: pseudo-R2 (McFadden 1974) relative to a null model
R2 <- function(model, null) 1 - as.numeric(logLik(model)) / as.numeric(logLik(null))

# Helper: number of estimated fixed-effect parameters (excludes thresholds/random effects)
n_fixef <- function(model) length(coef(model)) - length(model$alpha)
# NOTE: for clmm objects, coef() returns thresholds (alpha) + fixed effects (beta).
# model$alpha holds threshold estimates, so subtracting its length isolates the
# fixed-effect parameter count. Double-check this matches str(coef(M4)) for your model.

# ---- Build the ordered list of models ----
models <- list(Null = null, M0 = M0, M1 = M1, M2 = M2, M3 = M3, M4 = M4)

term_labels <- c(
  Null = "\u2014",
  M0   = "Block",
  M1   = "Host + Pathogen (mains)",
  M2   = "Block-by-population interactions",
  M3   = "Host:Pathogen interaction",
  M4   = "Three-way interaction (Block:Host:Pathogen)"
)

# ---- Extract per-model stats ----
model_names <- names(models)
df_vec      <- sapply(models, n_fixef)
logLik_vec  <- sapply(models, function(m) as.numeric(logLik(m)))
R2_vec      <- sapply(models, function(m) R2(m, null))

# ---- Pairwise LRTs between successive nested models ----
chisq_vec <- c(NA, rep(NA, length(models) - 1))
ddf_vec   <- c(NA, rep(NA, length(models) - 1))
p_vec     <- c(NA, rep(NA, length(models) - 1))

for (i in 2:length(models)) {
  lrt <- anova(models[[i - 1]], models[[i]])
  # anova.clmm returns a data.frame; the second row holds the comparison stats
  chisq_vec[i] <- lrt$LR.stat[2]
  ddf_vec[i]   <- lrt$df[2]
  p_vec[i]     <- lrt$`Pr(>Chisq)`[2]
}

# ---- Delta pseudo-R2 relative to the preceding model ----
dR2_vec <- c(NA, diff(R2_vec))

# ---- Assemble the table ----
summary_table <- data.frame(
  Model        = model_names,
  Term_added   = term_labels[model_names],
  df           = df_vec,
  logLik       = round(logLik_vec, 3),
  PseudoR2     = round(R2_vec, 4),
  Delta_PseudoR2 = round(dR2_vec, 4),
  ChiSq        = round(chisq_vec, 3),
  Delta_df     = ddf_vec,
  p_value      = signif(p_vec, 3),
  row.names    = NULL
)

print(summary_table)

# ---- Export ----
write.csv(summary_table, paste0(outdir, "/Supplementary_Table14.csv"), row.names = FALSE)
