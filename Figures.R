#### Script to reproduce the figures of the manuscript  "Rechsteiner M. et al. Host resistance and pathogen aggressiveness vary along an elevation gradient in a wild plant-pathogen system" ####
#### required libraries/packages ####
require(ggplot2)
require(ggpubr)
require(tidyverse)
require(argparser)
require(lme4)
require(ggmosaic)
require(Hmisc)
require(ggpubr)
#### Set up argparser for command line to read in data and to set the output directory for tables####
## Insert the whole file paths under default in this script, or set the working directory accordingly, or lastly, add the file paths after the corresponding arguments when executing this R script from the command line.

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

JoinedAbundance <- merge(PlantAbundance_Analysis, PathogenAbundance_Analysis, by = c("Season","Population","Elevation","Site"))
JoinedAbundance$Population <- factor(JoinedAbundance$Population, levels = c("High","Mid-high","Mid-low","Low"))
JoinedAbundance$Season <- factor(JoinedAbundance$Season, levels = c("2019","2020","2021"))

#### Figures ####
#### Host density and pathogen prevalence, Figure 1C and 1D ####
JoinedAbundance$HostDensity <- JoinedAbundance$Plantago.lanceolata/36
JoinedAbundance$Hosts_per_2500m2 <- round(JoinedAbundance$HostDensity*2500)
JoinedAbundance$Prevalence<- (JoinedAbundance$Infected.plants/JoinedAbundance$Hosts_per_2500m2)*100


summary_stats <- JoinedAbundance %>%
  group_by(Population, Season) %>%
  summarise(
    mean_average_Plan   = mean(HostDensity, na.rm = TRUE),
    se_average_Plan     = sd(HostDensity, na.rm = TRUE) / sqrt(n()),
    mean_elevation = mean(Elevation, na.rm = TRUE),
    n_sites        = n(),
    .groups = "drop"
  )

elev_range <- diff(range(PlantAbundance_Analysis$Elevation, na.rm = TRUE))
step <- elev_range * 0.04

summary_stats <- summary_stats %>%
  mutate(
    season_offset = case_when(
      Season == levels(Season)[1] ~ -step,
      Season == levels(Season)[2] ~ 0,
      Season == levels(Season)[3] ~ step,
      TRUE ~ 0
    ),
    Elevation_jit = mean_elevation + season_offset
  )

PlantAbundance_plot <- JoinedAbundance %>%
  ggplot(aes(x = Elevation, y = HostDensity)) +
  geom_point(aes(shape = Season, color = Population), size = 2.5, stroke = 1, alpha = 0.4)+
  scale_shape_manual(values = c(0, 5, 2))  +
  #scale_color_manual(values = c("#E11919","#EFA1A1","#8EA6CF","#184CA1")) +
  scale_color_manual(values = c("#184CA1","#8EA6CF","#EFA1A1","#E11919")) +
  scale_linetype_manual(values = c("dotted", "longdash","solid")) +
  ylab("Individuals per m2") +
  xlab("Elevation [m.a.s.l.]") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", colour = alpha("black", 1)),
    text = element_text(size = 20, colour = alpha("black", 1)),
    axis.text.x = element_text(colour = alpha("black", 1)),
    axis.text.y = element_text(colour = alpha("black", 1)),
    legend.key.width = unit(5,"cm"),
    legend.position = "none"
  ) +
  guides(color = guide_legend(override.aes = list(size = 7.5))) +
  guides(shape = guide_legend(override.aes = list(size = 7.5)))

PlantAbundance_plot <- PlantAbundance_plot +
  geom_errorbar(
    data = summary_stats,
    aes(x = Elevation_jit, y = mean_average_Plan,
        ymin = mean_average_Plan - se_average_Plan, ymax = mean_average_Plan + se_average_Plan,
        color = Population),
    width = 0, linewidth = 1.5, inherit.aes = FALSE
  ) +
  geom_point(
    data = summary_stats,
    aes(x = Elevation_jit, y = mean_average_Plan, shape = Season, color = Population),
    size = 6, stroke = 2, inherit.aes = FALSE
  )

### Prevalence plots

summary_stats <- JoinedAbundance %>%
  group_by(Population, Season) %>%
  summarise(
    mean_Prevalence   = mean(Prevalence, na.rm = TRUE),
    se_Prevalence     = sd(Prevalence, na.rm = TRUE) / sqrt(n()),
    mean_elevation = mean(Elevation, na.rm = TRUE),
    n_sites        = n(),
    .groups = "drop"
  )

elev_range <- diff(range(JoinedAbundance$Elevation, na.rm = TRUE))
step <- elev_range * 0.04

summary_stats <- summary_stats %>%
  mutate(
    season_offset = case_when(
      Season == levels(Season)[1] ~ -step,
      Season == levels(Season)[2] ~ 0,
      Season == levels(Season)[3] ~ step,
      TRUE ~ 0
    ),
    Elevation_jit = mean_elevation + season_offset
  )

Prevalence_plot_woOutl <- JoinedAbundance %>%
  ggplot(aes(x = Elevation, y = Prevalence)) +
  geom_point(aes(shape = Season, color = Population), size = 2.5, stroke = 1, alpha = 0.4)+
  scale_shape_manual(values = c(0, 5, 2))  +
  # scale_color_manual(values = c("#E11919","#EFA1A1","#8EA6CF","#184CA1")) +
  scale_color_manual(values = c("#184CA1","#8EA6CF","#EFA1A1","#E11919")) +
  scale_linetype_manual(values = c("dotted", "longdash","solid")) +
  ylab("Percentage of diseased plants") +
  xlab("Elevation [m.a.s.l.]") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", colour = alpha("black", 1)),
    text = element_text(size = 20, colour = alpha("black", 1)),
    axis.text.x = element_text(colour = alpha("black", 1)),
    axis.text.y = element_text(colour = alpha("black", 1)),
    legend.key.width = unit(5,"cm"),
    legend.position = "none"
  ) +
  ylim(0,2.5)+
  guides(color = guide_legend(override.aes = list(size = 7.5))) +
  guides(shape = guide_legend(override.aes = list(size = 7.5)))

Prevalence_plot_woOutl <- Prevalence_plot_woOutl +
  geom_errorbar(
    data = summary_stats,
    aes(x = Elevation_jit, y = mean_Prevalence,
        ymin = mean_Prevalence - se_Prevalence, ymax = mean_Prevalence + se_Prevalence,
        color = Population),
    width = 0, linewidth = 1.5, inherit.aes = FALSE
  ) +
  geom_point(
    data = summary_stats,
    aes(x = Elevation_jit, y = mean_Prevalence, shape = Season, color = Population),
    size = 6, stroke = 2, inherit.aes = FALSE
  )

Prevalence_plot_Outl <- JoinedAbundance %>%
  ggplot(aes(x = Elevation, y = Prevalence)) +
  geom_point(aes(shape = Season, color = Population), size = 2.5, stroke = 1, alpha = 0.4)+
  scale_shape_manual(values = c(0, 5, 2))  +
  #scale_color_manual(values = c("#E11919","#EFA1A1","#8EA6CF","#184CA1")) +
  scale_color_manual(values = c("#184CA1","#8EA6CF","#EFA1A1","#E11919")) +
  scale_linetype_manual(values = c("dotted", "longdash","solid")) +
  ylab("Percentage of diseased plants") +
  xlab("Elevation [m.a.s.l.]") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", colour = alpha("black", 1)),
    text = element_text(size = 20, colour = alpha("black", 1)),
    axis.text.x = element_text(colour = alpha("black", 1)),
    axis.text.y = element_text(colour = alpha("black", 1)),
    legend.key.width = unit(5,"cm"),
    legend.position = "none"
  ) +
  ylim(3.5,6)+
  guides(color = guide_legend(override.aes = list(size = 7.5))) +
  guides(shape = guide_legend(override.aes = list(size = 7.5)))

Prevalence_plot_Outl <- Prevalence_plot_Outl +
  geom_errorbar(
    data = summary_stats,
    aes(x = Elevation_jit, y = mean_Prevalence,
        ymin = mean_Prevalence - se_Prevalence, ymax = mean_Prevalence + se_Prevalence,
        color = Population),
    width = 0, linewidth = 1.5, inherit.aes = FALSE
  ) +
  geom_point(
    data = summary_stats,
    aes(x = Elevation_jit, y = mean_Prevalence, shape = Season, color = Population),
    size = 6, stroke = 2, inherit.aes = FALSE
  )

pdf(paste0(outdir,"/Figure1D.1.pdf"),width=5,height=5)
print(Prevalence_plot_Outl)
dev.off()

pdf(paste0(outdir,"/Figure1D.2.pdf"),width=5,height=5)
print(Prevalence_plot_woOutl)
dev.off()

pdf(paste0(outdir,"/Figure1C.pdf"),width=5,height=5)
print(PlantAbundance_plot)
dev.off()

#### Infection Matrix, Figure 2 ####
Infection_Severity_Analysis$Plant.genotype <- factor(Infection_Severity_Analysis$Plant.genotype, levels = c("I63","I42","I65","I59","I33","I55","I48","I35","I32","I11","I66","I58","I12","I13","I47",
                                                                                "A54","A51","A68","A52","A45","A611","A37","A44","A43","A65","A11","A15","A16","A33","A35",
                                                                                "N33","N42","N11","N59","N44","N58","N35","N34","N32","N43","N54","N511","N21","N510","N22",
                                                                                "O32","U21","U25","U11","U19","U15","O43","U110","O31","O312","U27","O41","U13","U23","O311"))

Infection_Severity_Analysis$Pathogen.isolate <- factor(Infection_Severity_Analysis$Pathogen.isolate, levels = c("U1M1","O4M2","O4M1","N5M1","N4M1","A5M2","A6M2","A5M3","A1M2","I5M2","I5M1"))
Infection_Severity_Analysis$Pathogen.isolate <- factor(Infection_Severity_Analysis$Pathogen.isolate, levels = c("I5M1","I5M2","A1M2","A5M3","A6M2","A5M2","N4M1","N5M1","O4M1","O4M2","U1M1"))

MeanBevanBlock12Day12PI <- Infection_Severity_Analysis %>% group_by(Plant.genotype,Pathogen.isolate,Host.population,Pathogen.population,Geography) %>% summarise(Bevan.mean=mean(Bevan.Score)) %>%
  ggplot( aes(Plant.genotype,Pathogen.isolate, fill = Bevan.mean)) +
  theme_void() +
  geom_tile(color = "grey15",
            lwd = 0.7,
            linetype = 1) +
  scale_fill_gradient(low = "white", high = "grey15") +
  labs(fill="Infection score")+
  theme(axis.text.x = element_text(angle = 45, hjust=1), axis.text.y= element_text(hjust=1),text = element_text(size = 10),plot.margin = unit(c(1,1,1,1), "cm"))

pdf(paste0(outdir,"/Figure2.pdf"),width=15,height=5)
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

pdf(paste0(outdir,"/Figure3.pdf"),width=10,height=5)
Infection_Status_plot %>%
  ggplot(mapping=aes(x=Pathogen.population:Host.population, y = Infection, fill = Pathogen.population:Host.population), width = 4) + 
  scale_fill_manual(values=c(rep("grey70",16)))+
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
         Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High"))) %>%
  mutate(Bevan.final = case_when(
    Bevan.Cat == "0" ~ "0",
    Bevan.Cat == "1" ~ "1",
    Bevan.Cat == "1.5" ~ "2",
    Bevan.Cat == "2.5" ~ "3",
    Bevan.Cat == "3" ~ "4",
    Bevan.Cat == "4" ~ "5"
  )) %>%
  #mutate(Bevan.final = factor(Bevan.final, levels = c("0","1","2","3","4","5")))
  mutate(Bevan.final = factor(Bevan.final, levels = c("5","4","3","2","1","0")))

pdf(paste0(outdir,"/Figure4.pdf"),width=10,height=5)
ggplot(Infection_Severity_plot, aes(fill=Bevan.final, x = Host.population:Pathogen.population)) + 
  geom_bar(position="fill", width = 0.8)+ theme_classic()+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.5))+scale_fill_manual(values = c("#333334","#565656","#737373","#ABACAC","#C6C6C6"))
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
  stat_summary(fun.data = mean_cl_boot, geom = "linerange", size = 2, color = c("grey70", "grey70", "grey70", "grey70")) +
  stat_summary(fun.y = mean, geom = "point", pch = 19, size = 7, color = c("grey70", "grey70", "grey70", "grey70")) +
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
  stat_summary(fun.data = mean_cl_boot, geom = "linerange", size = 2, color = c("grey70", "grey70", "grey70", "grey70")) +
  stat_summary(fun.y = mean, geom = "point", pch = 19, size = 7, color = c("grey70", "grey70", "grey70", "grey70")) +
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
  stat_summary(fun.data = mean_cl_boot, geom = "linerange", size = 2, color = c("grey70", "grey70", "grey70", "grey70")) +
  stat_summary(fun.y = mean, geom = "point", pch = 19, size = 7, color = c("grey70", "grey70", "grey70", "grey70")) +
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
  stat_summary(fun.data = mean_cl_boot, geom = "linerange", size = 2, color = c("grey70", "grey70", "grey70", "grey70")) +
  stat_summary(fun.y = mean, geom = "point", pch = 19, size = 7, color = c("grey70", "grey70", "grey70", "grey70")) +
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

pdf(paste0(outdir,"/SupplementaryFigure1.pdf"),width=15,height=5)
ggarrange(I_plot,A_plot,N_plot,OU_plot, ncol = 4, nrow = 1,labels = c("", "", "", ""))
dev.off()

#### Chasmothecia, Supp. Figure 2 ####
Chasmothecia_plot <- Chasmothecia_Analysis %>%
  mutate(Host.population = fct_relevel(Host.population, c("Low","Mid-low","Mid-high","High")),
         Pathogen.population = fct_relevel(Pathogen.population, c("Low","Mid-low","Mid-high","High")))

pdf(paste0(outdir,"/SupplementaryFigure2.pdf"),width=10,height=5)
Chasmothecia_plot %>%
  ggplot(mapping=aes(x=Host.population, y = Chasmothecia, fill = Pathogen.population)) + 
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

#### Figures were further edited using Inkscape/Adobe Illustrator ####