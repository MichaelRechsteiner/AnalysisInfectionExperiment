#### Script to reproduce the figures of the manuscript "Patterns of local adaptation reverse along an elevation gradient in a wild plant-pathogen system" ####
#### required libraries/packages ####
library(ggplot2)
library(ggpubr)
library(tidyverse)
require(argparser)
require(lme4)
require("ggmosaic")
#### Set up argparser for command line to read in data and to set the output directory for tables####
## Insert the whole file paths under default in this script or add the file paths after the corresponding arguments
args <- arg_parser(description="Parser to analyse data for the manuscript 'Patterns of local adaptation reverse along elevation in a wild plant-pathogen interaction'")
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
PlantAbundanceFile <- pargs[["PlantAbundanceFile"]]
PathogenAbundanceFile <- pargs[["PathogenAbundanceFile"]]

outdir <- pargs[["outdir"]]

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

#### Figures ####
#### Plant and pathogen abundance, Figure 1C and 1D ####
PlantAbundance_Analysis$Population <- factor(PlantAbundance_Analysis$Population, levels = c("Low","Mid-low","Mid-high","High"))
PlantAbundance_Analysis$Season <- as.factor(PlantAbundance_Analysis$Season)

PathogenAbundance_Analysis$Population <- factor(PathogenAbundance_Analysis$Population, levels = c("Low","Mid-low","Mid-high","High"))
PathogenAbundance_Analysis$Season <- as.factor(PathogenAbundance_Analysis$Season)
PathogenAbundance_Analysis$LNInfected.plants <- log(PathogenAbundance_Analysis$Infected.plants)

JoinedAbundance <- merge(PlantAbundance_Analysis, PathogenAbundance_Analysis, by = c("Season","Population","Elevation","Site"))

## Run regression models to add to the data plots
# create mixed linear model
# lmer_fit_plantago <- lmer(Plantago.lanceolata~Elevation + (1|Season), data=PlantAbundance_Analysis)
# summary(lmer_fit_plantago)
# lmer_fit_pathogen <- lmer(Infection.aggregates~Plantago.lanceolata+Elevation + (1|Season), data=JoinedAbundance)
# summary(lmer_fit_pathogen)
# lmer_fit_pathogen2 <- lmer(Infection.aggregates~Elevation + (1|Season), data=JoinedAbundance)
# summary(lmer_fit_pathogen)

# add predictions of the model in the existing data frame 
# PlantAbundance_Analysis <- PlantAbundance_Analysis %>% 
#   mutate(fit.m = predict(lmer_fit_plantago, re.form = NA), # marginal fit (only fixed effect)
#          fit.c = predict(lmer_fit_plantago, re.form = NULL), # conditional fir (fixed and random effect)
#          resid = resid(lmer_fit_plantago))
# 
# PathogenAbundance_Analysis <- JoinedAbundance %>% 
#   mutate(fit.m = predict(lmer_fit_pathogen, re.form = NA),
#          fit.c = predict(lmer_fit_pathogen, re.form = NULL),
#          resid = resid(lmer_fit_pathogen))

#### Plots

PlantAbundance_plot <- PlantAbundance_Analysis %>%
  ggplot(aes(x = Elevation, y = Plantago.lanceolata)) +geom_point(aes(shape = Season, color = Population), size = 3)+
  # geom_line(aes(y = fit.c, linetype = Season), size = 2)+
  scale_color_manual(values = c("#E11919","#EFA1A1","#8EA6CF","#184CA1")) +
  scale_linetype_manual(values = c("dotted", "longdash","solid"))+
  ylab("Number of plants") +
  xlab("Meters above sea level")+
  # ggtitle("Plantago lanceolata density\non Calanda from 2019-2021")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), text = element_text(size = 20),legend.key.width=unit(5,"cm")) +
  guides(color = guide_legend(override.aes = list(size = 7.5))) +
  guides(shape = guide_legend(override.aes = list(size = 7.5)))
pdf(paste0(outdir,"/PlantAbundance.pdf"),width=10,height=7.5)
print(PlantAbundance_plot)
dev.off()

PathogenAbundance_plot <- PathogenAbundance_Analysis %>%
  ggplot(aes(x = Elevation, y = LNInfected.plants)) +geom_point(aes(shape = Season, color = Population), size = 3)+
  # geom_line(aes(y = fit.c, linetype = Season), size = 2)+
  scale_color_manual(values = c("#E11919","#EFA1A1","#8EA6CF","#184CA1")) +
  scale_linetype_manual(values = c("dotted", "longdash","solid"))+
  ylab("ln(Number of infected plants)") +
  xlab("Meters above sea level")+
  # ggtitle("Infected plants\non Calanda from 2019-2021")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), text = element_text(size = 20),legend.key.width=unit(5,"cm")) +
  guides(color = guide_legend(override.aes = list(size = 7.5))) +
  guides(shape = guide_legend(override.aes = list(size = 7.5)))
pdf(paste0(outdir,"/PathogenAbundance.pdf"),width=10,height=7.5)
print(PathogenAbundance_plot)
dev.off()

#### Infection Matrix, Figure 2 ####
Infection_Severity_Analysis$Plant.genotype <- factor(Infection_Severity_Analysis$Plant.genotype, levels = c("I63","I42","I65","I59","I33","I55","I48","I35","I32","I11","I66","I58","I12","I13","I47",
                                                                                "A54","A51","A68","A52","A45","A611","A37","A44","A43","A65","A11","A15","A16","A33","A35",
                                                                                "N33","N42","N11","N59","N44","N58","N35","N34","N32","N43","N54","N511","N21","N510","N22",
                                                                                "O32","U21","U25","U11","U19","U15","O43","U110","O31","O312","U27","O41","U13","U23","O311"))

Infection_Severity_Analysis$Pathogen.isolate <- factor(Infection_Severity_Analysis$Pathogen.isolate, levels = c("U1M1","O4M2","O4M1","N5M1","N4M1","A5M2","A6M2","A5M3","A1M2","I5M2","I5M1"))
Infection_Severity_Analysis$Pathogen.isolate <- factor(Infection_Severity_Analysis$Pathogen.isolate, levels = c("I5M1","I5M2","A1M2","A5M3","A6M2","A5M2","N4M1","N5M1","O4M1","O4M2","U1M1"))

#### Mean Bevan scores day 12
MeanBevanBlock12Day12PI <- Infection_Severity_Analysis %>% group_by(Plant.genotype,Pathogen.isolate,Host.population,Pathogen.population,Geography) %>% summarise(Bevan.mean=mean(Bevan.Score)) %>%
  ggplot( aes(Plant.genotype,Pathogen.isolate, fill = Bevan.mean)) +
  theme_void() +
  geom_tile(color = "grey15",
            lwd = 0.7,
            linetype = 1) +
  scale_fill_gradient(low = "white", high = "grey15") +
  labs(fill="Infection score")+
  theme(axis.text.x = element_text(angle = 45, hjust=1), axis.text.y= element_text(hjust=1),text = element_text(size = 10),plot.margin = unit(c(1,1,1,1), "cm"))

pdf(paste0(outdir,"/MeanBevanBlock12Day12PI.pdf"),width=15,height=5)
print(MeanBevanBlock12Day12PI)
dev.off()

#### Infection outcome, Figure 3 ####
Infection_Status_plot <- Infection_Severity_Analysis %>% 
  # Assign an infection status (0/1) for each Host-Pathogen pair based on infection severity reached in both experimental blocks.
  group_by(Plant.genotype,Pathogen.isolate,Host.population,Pathogen.population, Geography) %>%
  summarise(Bevan.mean = mean(Bevan.Score)) %>% 
  mutate(Infection = case_when(Bevan.mean >= 1 ~ 1,
                               Bevan.mean < 1 ~ 0)) %>%
  mutate(Host.population = fct_relevel(Host.population, c("Low","Mid-low","Mid-high","High")),
         Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High")))

pdf(paste0(outdir,"/Infection.status.pdf"),width=10,height=5)
Infection_Status_plot %>%
  ggplot(mapping=aes(x=Pathogen.population:Host.population, y = Infection, fill = Pathogen.population:Host.population), width = 4) + 
  scale_fill_manual(values=c("grey25",rep("grey70",4),"grey25",rep("grey70",4),"grey25",rep("grey70",4),"grey25"))+
  stat_summary(fun.data=mean_sdl, geom="bar", width = 0.6)+
  ylim(0,1)+
  ylab("Proportion infected hosts")+
  # xlab("")+
  theme_classic()+
  theme(legend.position="none") + 
  ggtitle("Proportion of \nsuccessfully infected hosts")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),text = element_text(size=20),legend.position="none",axis.text.x = element_text(angle = 30,hjust=1))+
  stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=0.3)
dev.off()
#### Infection Severity, Figure 4 ####
Infection_Severity_plot <- Infection_Severity_Analysis %>% 
  # Prepare variables for analysis
  mutate(Bevan.Cat=as.factor(Bevan.Score),
         Block=as.factor(Block)) %>%
  filter(Bevan.Score >0) %>%
  mutate(Host.population = fct_relevel(Host.population, c("Low","Mid-low","Mid-high","High")),
         Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High")))

pdf(paste0(outdir,"/Infection.severity.pdf"),width=10,height=5)
ggplot(Infection_Severity_plot, aes(fill=Bevan.Cat, x = Host.population:Pathogen.population)) + 
  geom_bar(position="fill", width = 0.8)+ theme_classic()+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.5))+scale_fill_manual(values=c("#C6C6C6","#ABACAC","#737373","#565656","#333334"))
dev.off()

Infection_Severity_plot <- Infection_Severity_Analysis %>% 
  # Prepare variables for analysis
  mutate(Bevan.Cat=as.factor(Bevan.Score),
         Block=as.factor(Block)) %>%
  filter(Bevan.Score >0) %>%
  mutate(Host.population = fct_relevel(Host.population, c("Low","Mid-low","Mid-high","High")),
         Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High"))) %>%
  group_by(Host.population, Pathogen.population, Bevan.Cat) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(Host.population, Pathogen.population) %>%
  mutate(Frequency = Count / sum(Count)) %>%
  ungroup()

# Create the mosaic plot
Infection_Severity_plot<- droplevels(Infection_Severity_plot)
Colors1 = c("#D2D2D2" ,"#BDBEBE","#929292","#7D7D7D","#5C5C5C")
Colors2 = c("#C6C6C6","#ABACAC","#737373","#565656","#333334")

pdf(paste0(outdir,"/Infection.severity.matrix.pdf"),width=10,height=5)
Infection_Severity_plot %>% ggplot() +
  geom_mosaic(aes(weight = Frequency, x = product(Host.population, Pathogen.population), fill = Bevan.Cat), alpha = 1, width = 0.5) +
  labs(title = "Mosaic Plot of Infection Severity Frequencies",
       x = "Host and Pathogen Populations",
       y = "Host") +
  theme_classic()+scale_fill_manual(values= Colors2)
dev.off()
#### Days per sporulation, Supp. Figure 1 ####
filtered_data <- Days_to_spore_Analysis %>%
  filter(Host.population == "Low") %>%
  mutate(Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High")))
mean_data <- filtered_data %>%
  group_by(Pathogen.population) %>%
  summarise(Sporulation.time.mean = mean(Sporulation.time))
overall_mean <- mean(mean_data$Sporulation.time.mean)
plot <- ggplot(data = filtered_data, aes(x = Pathogen.population, y =Sporulation.time, fill = Geography)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "linerange", size = 2, color = c("grey70", "grey70", "grey70", "grey25")) +
  stat_summary(fun.y = mean, geom = "point", pch = 19, size = 7, color = c("grey70", "grey70", "grey70", "grey25")) +
  ylab("Sporulation time") +
  xlab("Pathogen population") +
  theme_minimal() +
  theme(legend.position = "none") + 
  ggtitle("Low elevation\n hosts") +  scale_y_continuous(breaks = 5:12, limits = c(5, 12)) +
  # scale_x_discrete(
  #   labels = c("Im Bofel","Arella", "Nesselboden", "Unter Alp")
  # )+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), text = element_text(size = 20), legend.position = "none",axis.text.x = element_text(angle = 30,hjust=1))

I_plot <- plot + 
  geom_hline(yintercept = overall_mean, line,  linetype = "dashed",color = "red", size = 1)


filtered_data <- Days_to_spore_Analysis %>%
  filter(Host.population == "Mid-low") %>%
  mutate(Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High")))
mean_data <- filtered_data %>%
  group_by(Pathogen.population) %>%
  summarise(Sporulation.time.mean = mean(Sporulation.time))
overall_mean <- mean(mean_data$Sporulation.time.mean)
plot <- ggplot(data = filtered_data, aes(x = Pathogen.population, y =Sporulation.time, fill = Geography)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "linerange", size = 2, color = c("grey70", "grey70", "grey70", "grey25")) +
  stat_summary(fun.y = mean, geom = "point", pch = 19, size = 7, color = c("grey70", "grey70", "grey70", "grey25")) +
  ylab("") +
  xlab("Pathogen population") +
  theme_minimal() +
  theme(legend.position = "none") + 
  ggtitle("Mid-low \nelevation hosts") +  scale_y_continuous(breaks = 5:12, limits = c(5, 12)) +
  # scale_x_discrete(
  #   labels = c("Im Bofel","Arella", "Nesselboden", "Unter Alp")
  # )+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), text = element_text(size = 20), legend.position = "none",axis.text.y = element_blank(),axis.text.x = element_text(angle = 30,hjust=1))

A_plot <- plot + 
  geom_hline(yintercept = overall_mean, line,  linetype = "dashed",color = "red", size = 1)


filtered_data <- Days_to_spore_Analysis %>%
  filter(Host.population == "Mid-high") %>%
  mutate(Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High")))
mean_data <- filtered_data %>%
  group_by(Pathogen.population) %>%
  summarise(Sporulation.time.mean = mean(Sporulation.time))
overall_mean <- mean(mean_data$Sporulation.time.mean)
plot <- ggplot(data = filtered_data, aes(x = Pathogen.population, y =Sporulation.time, fill = Geography)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "linerange", size = 2, color = c("grey70", "grey70", "grey70", "grey25")) +
  stat_summary(fun.y = mean, geom = "point", pch = 19, size = 7, color = c("grey70", "grey70", "grey70", "grey25")) +
  ylab("") +
  xlab("Pathogen population") +
  theme_minimal() +
  theme(legend.position = "none") + 
  ggtitle("Mid-high\nelevation hosts") +  scale_y_continuous(breaks = 5:12, limits = c(5, 12)) +
  # scale_x_discrete(
  #   labels = c("Im Bofel","Arella", "Nesselboden", "Unter Alp")
  # )+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), text = element_text(size = 20), legend.position = "none",axis.text.y = element_blank(),axis.text.x = element_text(angle = 30,hjust=1))

N_plot <- plot + 
  geom_hline(yintercept = overall_mean, line,  linetype = "dashed",color = "red", size = 1)


filtered_data <- Days_to_spore_Analysis %>%
  filter(Host.population == "High") %>%
  mutate(Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High")))
mean_data <- filtered_data %>%
  group_by(Pathogen.population) %>%
  summarise(Sporulation.time.mean = mean(Sporulation.time))
overall_mean <- mean(mean_data$Sporulation.time.mean)
plot <- ggplot(data = filtered_data, aes(x = Pathogen.population, y =Sporulation.time, fill = Geography)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "linerange", size = 2, color = c("grey70", "grey70", "grey70", "grey25")) +
  stat_summary(fun.y = mean, geom = "point", pch = 19, size = 7, color = c("grey70", "grey70", "grey70", "grey25")) +
  ylab("") +
  xlab("Pathogen population") +
  theme_minimal() +
  theme(legend.position = "none") + 
  ggtitle("High elevation\n hosts") +  scale_y_continuous(breaks = 5:12, limits = c(5, 12)) +
  # scale_x_discrete(
  #   labels = c("Im Bofel","Arella", "Nesselboden", "Unter Alp")
  # )+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), text = element_text(size = 20), legend.position = "none",axis.text.y = element_blank(),axis.text.x = element_text(angle = 30,hjust=1))

OU_plot <- plot + 
  geom_hline(yintercept = overall_mean, line,  linetype = "dashed",color = "red", size = 1)

pdf(paste0(outdir,"/Sporulation.time.pdf"),width=15,height=5)
ggarrange(I_plot,A_plot,N_plot,OU_plot, ncol = 4, nrow = 1,labels = c("", "", "", ""))
dev.off()

#### Chasmothecia, Supp. Figure 2 ####
Chasmothecia_plot <- Chasmothecia_Analysis %>%
  mutate(Host.population = fct_relevel(Host.population, c("Low","Mid-low","Mid-high","High")),
         Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High")))

pdf(paste0(outdir,"/Chasmothecia.pdf"),width=15,height=5)
Chasmothecia_plot %>%
  ggplot(mapping=aes(x=Host.population:Pathogen.population, y = Chasmothecia, fill = Host.population:Pathogen.population), width = 4) + 
  scale_fill_manual(values=c("grey25",rep("grey70",4),"grey25",rep("grey70",4),"grey25",rep("grey70",4),"grey25"))+
  stat_summary(fun.data=mean_sdl, geom="bar", width = 0.6)+
  ylim(0,1)+
  ylab("Proportion")+
  xlab("")+
  theme_classic()+
  theme(legend.position="none") + 
  ggtitle("Proportion of infections \n developing chasmothecia")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),text = element_text(size=25),legend.position="none",axis.text.x = element_text(angle = 30,hjust=1))+
  stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=0.3)
dev.off()

pdf(paste0(outdir,"/Chasmothecia.pdf"),width=10,height=5)
Chasmothecia_plot %>%
  ggplot(mapping=aes(x=Host.population, y = Chasmothecia, fill = Pathogen.population), width = 2) + 
  scale_fill_manual(values=c(rep("grey70",4)))+
  stat_summary(fun.data=mean_sdl, geom="bar", width = 0.6,position = position_dodge(width = 0.7))+
  ylim(0,1)+
  ylab("Proportion of infections\n developing chasmothecia")+
  # xlab("")+
  theme_classic()+
  theme(legend.position="none") + 
  ggtitle("Presence of chasmothecia")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),text = element_text(size=20),legend.position="none",axis.text.x = element_text(angle = 30,hjust=1))+
  stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=0.3,position = position_dodge(width = 0.7))
dev.off()

#### Figures were further edited using Adobe Illustrator