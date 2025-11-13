##############################################################################################
# A HIGH THROUGHPUT SEQUENCING SCREEN TO IDENTIFY APOMIXIS IN PLANTS - results visualization #
##############################################################################################

library(ggplot2)
library(tidyverse)
library(rmarkdown)
library(dplyr)
library(ggnewscale)

working_directory="/Users/Goeckeritz/Desktop/Desktop - Charity’s MacBook Pro/PRFB_apomixis_project/malus/Genes_and_Development_Screening_Pipeline_Paper/"
setwd(working_directory)

################################################# 6X depth embryos (the pilot study) ##############################################

#read in depth file for pilot study (up to 96 samples on a sequencing lane)
raw_cov_16moms <- read.table("embryo_raw_coverage.txt", header=FALSE) #sorry... sometimes I use coverage and depth interchangeably..!
colnames(raw_cov_16moms) <- c("embryo","raw_cov")

#read in comparison results for pilot study
results_16moms <- read.table("all_results_combined_max2_alleles.tsv", header=TRUE)

#join the coverage information for each embryo with final results, make mother and embryo a factor.
all_results_16moms <- merge(results_16moms, raw_cov_16moms, by="embryo")
all_results_16moms$embryo=as.factor(all_results_16moms$embryo)
all_results_16moms$mother=as.factor(all_results_16moms$mother)

#A simple ANOVA - is there a statistically significant association between mother genotype and raw coverage, which we will define as the sequencing quality?
str(all_results_16moms)
all_results_16moms_no_sims <- all_results_16moms[!grepl("sim", all_results_16moms$embryo), ]
all_results_16moms_sims <- all_results_16moms[grepl("sim", all_results_16moms$embryo), ]
  
library(lme4)
library(car)
model = lm(raw_cov ~ mother, data=all_results_16moms_no_sims) #since mother and lane are confounded factors, I figured it didn't make much sense to test the associations together. 
#But for embryos sequenced on only one lane (6X target depth samples), 
#there was no association of coverage and what lane the embryo was sequenced on. 
#Makes sense, since all sequencing lanes gave me about the same amount of data.
plot(residuals(model))
hist(residuals(model)) #mostly normal, but some skew
Anova(model, type="III") #not unexpected - didn't bother to remove potential outliers from this analysis since ANOVA is robust against deviances from normality.  

#let's plot the coverage of each embryo, but color them by the mother tree, so we can see the bias in different mother plants. (Supplementary Figure 1)
a <- ggplot(all_results_16moms_no_sims, aes(x=embryo, y=raw_cov, fill=mother)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Mspectabilis_588917"="hotpink2","Masiatica_594099"="goldenrod","Mhupehensis_633818"="chocolate4","Msargentii_589400"="forestgreen","Madstringens_588898"="pink","Mbaccata_588907"="cyan3","Mbaccata_613807"="firebrick","Mbrevipes_589170"="lightblue","Mfloribunda_589827"="lavender","Mhalliana_590029"="coral2","Mhupehensis_589756"="steelblue","MIndianSummer_589733"="grey50",
                             "MMaryPotter_588983"="bisque","Mplatycarpa_589415"="hotpink4","MProfusion_589449"="orange","Mtoringoides_588930"="gold"))+
  labs(title="Average depth for 6X embryos", x = "Embryo (grouped by maternal parent)", y="Average read depth per site", fill="Maternal Parent") +
  scale_y_continuous(limits=c(0,25), breaks = seq(0, 25, by = 1), expand=c(0,0))+
  theme_bw()+
  theme(plot.title=element_text(hjust = 0.5, size=50, face="bold"), panel.spacing = unit(1, "lines"), strip.text = element_text(size = 30, face="bold"),
        strip.background = element_rect(fill = "white"), plot.margin = margin(0.5, 0.2, 0.2, 0.2, "in")) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(hjust = 1, face="bold", size=30), axis.title.x=element_text(face="bold", size=40), axis.title.y=element_text(face="bold", size=40), legend.text=element_text(size=30),
        legend.title=element_text(size=40, face="bold"), legend.key.size = unit(1.5, "cm"))
a

ggsave("embryo_coverage_grouped_by_16moms_pilot.png", plot=a, width = 28, height = 24 , units = "in", bg="white")

################## Generate Figure 1... #####################
#Figure 1A
mean(all_results_16moms_no_sims$raw_cov)
median(all_results_16moms_no_sims$raw_cov)
sd(all_results_16moms_no_sims$raw_cov)

b <- ggplot(all_results_16moms_no_sims, aes(x = raw_cov)) +
  geom_histogram(binwidth = 0.25, fill = "#5DADE2", color="white", alpha=0.85) +
  labs(
    title = "Embryo Depth Distribution (6X embryos)",
    x = "Embryo Depth", y = "Number of Embryos"
  ) +
  scale_x_continuous(breaks = seq(0, 24, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 15, 1), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 24), ylim = c(0, 15))+
  geom_vline(xintercept=6.1, linetype="dashed", color = "firebrick3", size=1.5)+
  geom_vline(xintercept=5.9, linetype="dashed", color = "gold", size=1.5)+
  theme_bw()+
  theme(plot.title=element_text(hjust = 0.5, size=40, face="bold"), panel.spacing = unit(1, "lines"), strip.text = element_text(size = 30, face="bold"),
        strip.background = element_rect(fill = "white"), plot.margin = margin(0.2, 0.2, 0.2, 0.2, "in")) +
  theme(axis.text.y = element_text(hjust = 1, face="bold", size=25), axis.text.x = element_text(hjust = 1, face="bold", size=25), axis.title.x=element_text(face="bold", size=40), axis.title.y=element_text(face="bold", size=40), legend.text=element_text(size=30),
        legend.key.size = unit(1.5, "cm"), axis.ticks=element_line(linewidth=1), axis.ticks.length = unit(0.1, "cm"))
b

ggsave("embryo_distribution_6X_depth.png", plot=b, width = 24, height = 12 , units = "in", bg="white")

range(all_results_16moms_no_sims$raw_cov)

#3X embryos - just the true embryos (Figure 1C)
raw_cov_3X_moms <- read.table("true_3X_embryos_only_raw_cov.txt", header=FALSE) #sorry... sometimes I use coverage and depth interchangeably..!
colnames(raw_cov_3X_moms) <- c("embryo","raw_cov")
median(raw_cov_3X_moms$raw_cov)
mean(raw_cov_3X_moms$raw_cov)
sd(raw_cov_3X_moms$raw_cov)
range(raw_cov_3X_moms$raw_cov)

bb <- ggplot(raw_cov_3X_moms, aes(x = raw_cov)) +
  geom_histogram(binwidth = 0.25, fill = "#8E236B", color="white", alpha=0.85) +
  labs(
    title = "Embryo Depth Distribution (3X embryos)",
    x = "Embryo Depth", y = "Number of Embryos"
  ) +
  scale_x_continuous(breaks = seq(0, 14, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 60, 5), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 14), ylim = c(0, 60))+
  geom_vline(xintercept=2.9, linetype="dashed", color = "firebrick3", size=1.5)+
  geom_vline(xintercept=2.5, linetype="dashed", color = "gold", size=1.5)+
  theme_bw()+
  theme(plot.title=element_text(hjust = 0.5, size=40, face="bold"), panel.spacing = unit(1, "lines"), strip.text = element_text(size = 30, face="bold"),
        strip.background = element_rect(fill = "white"), plot.margin = margin(0.2, 0.2, 0.2, 0.2, "in")) +
  theme(axis.text.y = element_text(hjust = 1, face="bold", size=25), axis.text.x = element_text(hjust = 1, face="bold", size=25), axis.title.x=element_text(face="bold", size=40), axis.title.y=element_text(face="bold", size=40), legend.text=element_text(size=30),
        legend.key.size = unit(1.5, "cm"), axis.ticks=element_line(linewidth=1), axis.ticks.length = unit(0.1, "cm"))

bb

ggsave("embryo_distribution_3X_depth.png", plot=bb, width = 24, height = 12 , units = "in", bg="white")

################## plotting the relationship between raw average depth and number of biallelelic sites compared (Figure 1B)

qq <- ggplot(all_results_16moms_no_sims, aes(x=raw_cov, y=log(total_sites), color=mother)) +
  geom_point(shape=17, size=3.5)+
  scale_color_manual(name= "Maternal Parent", values=c("Mspectabilis_588917"="hotpink2","Masiatica_594099"="goldenrod","Mhupehensis_633818"="chocolate4","Msargentii_589400"="forestgreen",
                                                       "Madstringens_588898"="pink","Mbaccata_588907"="cyan3","Mbaccata_613807"="firebrick","Mbrevipes_589170"="lightblue","Mfloribunda_589827"="lavender","Mhalliana_590029"="coral2","Mhupehensis_589756"="steelblue","MIndianSummer_589733"="grey50",
                                                       "MMaryPotter_588983"="bisque","Mplatycarpa_589415"="hotpink4","MProfusion_589449"="orange","Mtoringoides_588930"="gold"))+
  labs(title="Sequencing Depth vs Biallelic Sites Compared", x = "Embryo Depth", y="Number of Biallelic Sites (Log Scale)") +
  scale_y_continuous(limits=c(0,13), breaks = seq(0, 13, by = 1), expand = c(0, 0))+
  scale_x_continuous(limits=c(0,25), breaks = seq(0, 25, by = 1), expand = c(0, 0))+
  #scale_x_continuous(limits=c(60,100))+
  geom_vline(xintercept=3,linetype="dashed",color="black",size=1.5)+
  geom_vline(xintercept=6,linetype="dashed",color="black",size=1.5)+
  geom_hline(yintercept=6.9,linetype="dashed",color="firebrick",size=1.5)+
  geom_hline(yintercept=4.6,linetype="dashed",color="firebrick",size=1.5)+
  geom_hline(yintercept=9.9,linetype="dashed",color="firebrick",size=1.5)+
  geom_text(aes(x=8, label="1000 sites", y=7.2), color="black", size=10)+
  geom_text(aes(x=8, label="100 sites", y=4.9), color="black", size=10)+
  geom_text(aes(x=8.2, label="20000 sites", y=10.2), color="black", size=10)+
  theme_bw()+
  theme(plot.title=element_text(hjust = 0.5),
        strip.background = element_rect(fill = "white"), plot.margin = margin(0.5, 0.2, 0.2, 0.2, "in")) +
  theme(axis.text.x = element_text(face="bold", size=25), axis.text.y = element_text(hjust = 1, face="bold", size=25), axis.title.x=element_text(face="bold", size=30), axis.title.y=element_text(face="bold", size=30), plot.title=element_text(size=35, face = "bold", hjust=0.5), 
        legend.text=element_text(size=25), legend.title=element_text(size=30, face="bold", hjust=0.13), legend.key.size = unit(1.2, "cm")) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 5)))
qq

ggsave("raw_depth_vs_number_of_sites_compared.png", plot=qq, width = 24, height = 12 , units = "in", bg="white")

###################just to see untransformed data #######################
qqq <- ggplot(all_results_16moms_no_sims, aes(x=raw_cov, y=total_sites, color=mother)) +
  geom_point(shape=17, size=4)+
  scale_color_manual(values=c("Mspectabilis_588917"="hotpink2","Masiatica_594099"="goldenrod","Mhupehensis_633818"="chocolate4","Msargentii_589400"="forestgreen",
                              "Madstringens_588898"="pink","Mbaccata_588907"="cyan3","Mbaccata_613807"="firebrick","Mbrevipes_589170"="lightblue","Mfloribunda_589827"="lavender","Mhalliana_590029"="coral2","Mhupehensis_589756"="steelblue","MIndianSummer_589733"="grey50",
                              "MMaryPotter_588983"="bisque","Mplatycarpa_589415"="hotpink4","MProfusion_589449"="orange","Mtoringoides_588930"="gold"))+
  ggtitle("Relationship Between Sequencing Depth and Number of Biallelic Sites Compared")+
  xlab("Average depth of embryo")+
  ylab("Total number of biallelic sites compared (natural log scale)")+
  scale_y_continuous(limits=c(0,275000), breaks = seq(0, 275000, by = 10000), expand = c(0, 0))+
  scale_x_continuous(limits=c(0,25), breaks = seq(0, 25, by = 1), expand = c(0, 0))+
  #scale_x_continuous(limits=c(60,100))+
  theme_bw()+
  theme(plot.title=element_text(hjust = 0.5),
        strip.background = element_rect(fill = "white"), plot.margin = margin(0.5, 0.05, 0.05, 0.05, "in")) +
  theme(axis.text.x = element_text(face="bold", size=25), axis.text.y = element_text(hjust = 1, face="bold", size=25), axis.title.x=element_text(face="bold", size=30), axis.title.y=element_text(face="bold", size=30), plot.title=element_text(size=35, face = "bold", hjust=0.5), 
        legend.text=element_text(size=30), legend.title=element_text(size=30, face="bold", hjust=0.13), legend.key.size = unit(1.2, "cm")) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 5)))

qqq

######################## Get some metrics on the controls, like average number of sites compared. ###############################

tmp <- all_results_16moms_no_sims %>%
  filter(mother %in% c("Mspectabilis_588917", "Masiatica_594099", "Mhupehensis_633818", "Msargentii_589400")) %>% 
  group_by(mother) %>%
  summarise(
    depth_mean_value = mean(raw_cov, na.rm = TRUE),
    sites_avg_value  = mean(total_sites, na.rm = TRUE),
    min_embryo = embryo[which.min(raw_cov)],
    median_embryo = embryo[median(raw_cov)],
    min_raw_cov = min(raw_cov, na.rm = TRUE),
    min_total_sites = total_sites[which.min(raw_cov)],
    max_embryo = embryo[which.max(raw_cov)],
    max_raw_cov = max(raw_cov, na.rm = TRUE),
    max_total_sites = total_sites[which.max(raw_cov)],
    max_total_percent_matches = total_match_percentage[which.max(total_match_percentage)],
    .groups = "drop"
  )

tmp

############################# Plot the distribution of two embryos, highest and lowest depth, for Mspectabilis_588917 and Mhupehensis_633818. For Mspectabilis,
#that would be e21 (low, 0.3X depth) and e5 (high, 14.2X depth) and for Mhupehensis, that would be e6 (low, 1.6X depth) and e22 (high, 12.9X depth)
#just read in a different file and changed the labels and such to plot the others. 

Mspectabilis_588917_e21 <- read.table("Mspectabilis_588917_e21_only_cleaned.table", header=TRUE)
names(Mspectabilis_588917_e21) <- gsub("\\.GT", "", names(Mspectabilis_588917_e21))

#report minimum and maximum distances 
distances <- Mspectabilis_588917_e21 %>%
  arrange(CHROM, POS) %>%
  group_by(CHROM) %>%
  summarise(
    min_dist = min(diff(POS)),   # smallest gap
    max_dist = max(diff(POS)),    # largest gap
    median_dist = median(diff(POS)),
    avg_dist = mean(diff(POS))
  )


q <- ggplot(Mspectabilis_588917_e21, aes(x=POS, y=CHROM)) +
  geom_point(color="hotpink2", size=5) +
  scale_y_continuous(
    breaks = seq(min(Mspectabilis_588917_e21$CHROM),
                 max(Mspectabilis_588917_e21$CHROM),
                 by = 1)
  ) +
  scale_x_continuous(
    breaks = seq(0, 55000000, by = 1e6),
    labels = function(x) x/1e6
  ) +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none"   # 🔹 removes legend
  ) +
  annotate("segment", 
           x = 1, 
           xend = c(33734251,35818875,35674313,30921448,46022136,
                    34967166,35401034,31655921,36128915,
                    42917823,41459566,33526175,44586254,
                    31725111,55653390,39306085,34066592), 
           y = 1:17, 
           yend = 1:17,
           colour = "black") +
  ggtitle("Distribution of Sites Compared to Maternal Genotype") +
  labs(subtitle="M. spectabilis 588917 vs Mspectabilis 588917 e21 (n = 128 sites)", y="Chromosome", x="Position (Mb)")+
  theme(plot.title=element_text(hjust = 0.5, size=50, face="bold"), plot.subtitle=element_text(hjust=0.5, size=30), panel.spacing = unit(1, "lines"), strip.text = element_text(size = 30, face="bold"),
        strip.background = element_rect(fill = "white"), plot.margin = margin(0.05, 0.05, 0.05, 0.05, "in")) +
  theme(axis.text.y = element_text(hjust = 1, face="bold", size=30), axis.text.x = element_text(hjust = 1, face="bold", size=30), axis.title.x=element_text(face="bold", size=40), axis.title.y=element_text(face="bold", size=40), legend.text=element_text(size=30),
        legend.key.size = unit(1.5, "cm"))
q

ggsave("Mspectabilis_58917_e21_sites_distribution.png", plot=q, width = 40, height = 20 , units = "in", bg="white")


################################## Mhup633818
Mhupehensis_633818_e20 <- read.table("Mhupehensis_633818_e20_only_cleaned.table", header=TRUE)
names(Mhupehensis_633818_e20) <- gsub("\\.GT", "", names(Mhupehensis_633818_e20))

q <- ggplot(Mhupehensis_633818_e20, aes(x=POS, y=CHROM)) +
  geom_point(color="chocolate4", size=5) +
  scale_y_continuous(
    breaks = seq(min(Mhupehensis_633818_e20$CHROM),
                 max(Mhupehensis_633818_e20$CHROM),
                 by = 1)
  ) +
  scale_x_continuous(
    breaks = seq(0, max(Mhupehensis_633818_e20$POS), by = 1e6),
    labels = function(x) x/1e6
  ) +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none"   # 🔹 removes legend
  ) +
  annotate("segment", 
           x = 1, 
           xend = c(33734251,35818875,35674313,30921448,46022136,
                    34967166,35401034,31655921,36128915,
                    42917823,41459566,33526175,44586254,
                    31725111,55653390,39306085,34066592), 
           y = 1:17, 
           yend = 1:17,
           colour = "black") +
  ggtitle("Distribution of Sites Compared to Maternal Genotype") +
  labs(subtitle="M. hupehensis 633818 vs M. hupehensis 633818 e20 (n = 80,528 sites)", y="Chromosome", x="Position (Mb)")+
  theme(plot.title=element_text(hjust = 0.5, size=50, face="bold"), plot.subtitle=element_text(hjust=0.5, size=30), panel.spacing = unit(1, "lines"), strip.text = element_text(size = 30, face="bold"),
        strip.background = element_rect(fill = "white"), plot.margin = margin(0.05, 0.05, 0.05, 0.05, "in")) +
  theme(axis.text.y = element_text(hjust = 1, face="bold", size=30), axis.text.x = element_text(hjust = 1, face="bold", size=30), axis.title.x=element_text(face="bold", size=40), axis.title.y=element_text(face="bold", size=40), legend.text=element_text(size=30),
        legend.key.size = unit(1.5, "cm"))
q

ggsave("Mhupehensis_633818_e20_sites_distribution.png", plot=q, width = 40, height = 20 , units = "in", bg="white")

############### what's the range of % similarities of putative clonal embryos for M. hupehensis 633818 and M. sargentii 589400?
#to be conservative, we are going to exclude the three in hupehensis that deviate slightly from the main cluster. Those may be selfs, even if they were sequenced at slightly lower depths.

Mhupehensis_633818 <- all_results_16moms_no_sims %>%
  dplyr::filter(mother %in% c("Mhupehensis_633818")) %>%
  dplyr::filter(!embryo %in% c("Mhupehensis_633818_e21")) #just excludes the one that's probably a BIII
range(Mhupehensis_633818$het_match_percentage)
range(Mhupehensis_633818$total_match_percentage)

Mhupehensis_633818_2 <- all_results_16moms_no_sims %>%
  dplyr::filter(mother %in% c("Mhupehensis_633818")) %>%
  dplyr::filter(!embryo %in% c("Mhupehensis_633818_e21", "Mhupehensis_633818_e6","Mhupehensis_633818_e9","Mhupehensis_633818_e2"))

plot(Mhupehensis_633818_2$total_match_percentage, Mhupehensis_633818_2$raw_cov) #clear correlation
boxplot(Mhupehensis_633818_2$het_match_percentage) #hey I'm pretty good at eyeballing these outliers lol

Msargentii_589400 <- all_results_16moms_no_sims %>%
  dplyr::filter(mother %in% c("Msargentii_589400")) %>%
  dplyr::filter(!embryo %in% c("Msargentii_589400_e1", "Msargentii_589400_e15", "Msargentii_589400_e18", "Msargentii_589400_e2",
                               "Msargentii_589400_e4","Msargentii_589400_e5","Msargentii_589400_e6", "Msargentii_589400_e7", "Msargentii_589400_e8")) 
range(Msargentii_589400$total_match_percentage)
plot(Msargentii_589400$total_match_percentage, Msargentii_589400$raw_cov) #seems there is a clear correlation... but the percent diff between like 3X and 12X is only ~1.5%

######################## combine the "true clonal embryos" datasets from the two moms and find a more mathematically sound range... using the IQR 1.5 rule. 
library(stats)
clonal_embryos2 <- rbind(Mhupehensis_633818_2, Msargentii_589400) #excludes the 3 het outliers from Mhupehensis
boxplot(clonal_embryos2$het_match_percentage) # <-- well now there's another outlier now. Ha. but we'll leave him. can't keep removing them forever. 
boxplot(clonal_embryos2$hom_match_percentage) # no obvious outliers
hist(clonal_embryos2$het_match_percentage)
hist(clonal_embryos2$hom_match_percentage) #little odd, but could just be the automatic bin setting
quantile(clonal_embryos2$het_match_percentage, probs = c(0.025, 0.975))
range(clonal_embryos2$het_match_percentage)
quantile(clonal_embryos2$hom_match_percentage, probs = c(0.025, 0.975))
range(clonal_embryos2$hom_match_percentage)
range(clonal_embryos2$total_match_percentage)


################################### Plot the results for the controls #####################################

#set the order of moms presented
all_results_16moms_no_sims$mother <- factor(all_results_16moms_no_sims$mother, levels=c("Mspectabilis_588917","Masiatica_594099","Madstringens_588898","Mbaccata_588907","Mbaccata_613807","Mbrevipes_589170","Mfloribunda_589827","Mhalliana_590029", 
                                                                  "Mhupehensis_633818","Msargentii_589400","Mhupehensis_589756","MIndianSummer_589733","MMaryPotter_588983" ,"Mplatycarpa_589415","MProfusion_589449","Mtoringoides_588930"))

#set a filter on which embryos to consider. 
controls <- all_results_16moms_no_sims %>%
  dplyr::filter(mother %in% c("Mspectabilis_588917","Masiatica_594099","Mhupehensis_633818","Msargentii_589400")) %>%
  dplyr::filter(het_sites_total > 50) %>%
  dplyr::filter(hom_sites_total > 50)

ctl <- ggplot(controls, aes(x=het_match_percentage, y=hom_match_percentage, color=mother)) +
  geom_point(shape=17, size=5)+
  scale_color_manual(values=c("Mspectabilis_588917"="hotpink2","Masiatica_594099"="goldenrod","Mhupehensis_633818"="chocolate4","Msargentii_589400"="forestgreen"))+
  ggtitle("Genetic Similarity")+
  #labs(subtitle="(No Error Adjustment)")+
  xlab("Mom Het sites (% matching sites)")+
  ylab("Mom Hom sites (% matching sites)")+
  scale_y_continuous(limits=c(20,100), breaks = seq(20, 100, by = 5))+
  scale_x_continuous(limits=c(40,100), breaks = seq(40, 100, by = 5))+
  facet_wrap(~mother, nrow=1,ncol=4)+
  coord_fixed()+
  #scale_x_continuous(limits=c(60,100))+
  geom_vline(xintercept=92.5,linetype="dashed","black",2)+
  geom_hline(yintercept=95,linetype="dashed","black",2)+
  theme_bw()+
  theme(plot.title=element_text(hjust = 0.5), plot.subtitle=element_text(hjust=0.5, size=20), panel.spacing = unit(1, "lines"), strip.text = element_text(size = 20, face="bold"),
        strip.background = element_rect(fill = "white"), plot.margin = margin(0.5, 0.05, 0.05, 0.05, "in")) +
  theme(axis.text.x = element_text(face="bold", size=20), axis.text.y = element_text(hjust = 1, face="bold", size=20), axis.title.x=element_text(face="bold", size=30), axis.title.y=element_text(face="bold", size=30), plot.title=element_text(size=30, face = "bold", hjust=0.5), legend.position = "none")

ctl

ggsave("control_plots_6X.png", plot=ctl, width = 28, height = 24 , units = "in", bg="white")


############################# Read in 3X results ##################

#read in comparison results
results_3Xembryos <- read.table("all_results_combined_3X.tsv", header=TRUE)
raw_cov_all_3X_sims_and_all_incl_spikeins <- read.table("all_3X_data_individuals_and_moms_raw_cov.txt", header=FALSE) 
colnames(raw_cov_all_3X_sims_and_all_incl_spikeins) <- c("embryo","raw_cov")

#join the coverage information for each embryo with final results, make mother and embryo a factor. 861 embryos total - 58 drop out when we require at least 50 sites of each category for comparison.
all_results_3Xembryos <- merge(results_3Xembryos, raw_cov_all_3X_sims_and_all_incl_spikeins, by="embryo")
all_results_3Xembryos$embryo=as.factor(all_results_3Xembryos$embryo)
all_results_3Xembryos$mother=as.factor(all_results_3Xembryos$mother)

all_results_3Xembryos_true_no_sims <- all_results_3Xembryos %>%
  dplyr::filter(str_detect(embryo, "_e")) %>%
  dplyr::filter(!str_detect(embryo, "sim")) %>%
  dplyr::filter(het_sites_total > 50) %>%
  dplyr::filter(hom_sites_total > 50) %>%
  dplyr::filter(!embryo %in% c("Mzhaojiaoensis_633816_e1", "Mmicromalus_594092_e2")) #dropped since there were two notes during extraction about possible contamination

#803/861 left for analysis after the 50 site and contamination filters. Not too bad. 

#################### let's look at the spike-ins##############################
all_spikeins <- all_results_3Xembryos %>%
  dplyr::filter(str_detect(embryo, "_[abc]")) %>%
  dplyr::filter(het_sites_total > 50) %>%
  dplyr::filter(hom_sites_total > 50)

#site filter drops the MHenryFDupont and Mioensis spikeins + one from Mtoringo 633814.
hist(all_spikeins$total_match_percentage)
boxplot(all_spikeins$het_match_percentage) #LOL those platycarpas... what the
all_spikeins$embryo[which(all_spikeins$total_match_percentage < 90)] #what's wrong with these little monsters wtf

all_results_3Xembryos_spikeins_no4X_platycarpas <- all_results_3Xembryos %>%
  dplyr::filter(str_detect(embryo, "_[abc]")) %>%
  dplyr::filter(!mother %in% c("Mplatycarpa_588752", "Mplatycarpa_588847", "Mplatycarpa_589356", "Mplatycarpa_589415")) %>%
  dplyr::filter(het_sites_total > 50) %>%
  dplyr::filter(hom_sites_total > 50)

###### Supplementary Figure 4, showing the spike ins for tetraploid M platycarpas don't match as well to their maternal libraries, as all other spikeins from other species do. 
cc <- ggplot(all_spikeins, aes(x = total_match_percentage)) +
  geom_histogram(binwidth = 0.1, fill = "gray75", color = "black") +
  ggtitle("Distribution of Genetic Similarities") +
  labs(subtitle="(Spike-ins compared to full maternal library)")+
  xlab("% Similarities")+
  ylab("Number of spike-ins")+
  scale_x_continuous(limits=c(80,100), breaks = seq(80, 100, by = 1), expand = c(0, 0))+
  scale_y_continuous(limits=c(0,15), breaks = seq(0, 15, by = 1), expand = c(0, 0))+
  theme_bw()+
  theme(plot.title=element_text(hjust = 0.5, size=50, face="bold"), plot.subtitle=element_text(hjust=0.5, size=30), panel.spacing = unit(1, "lines"), strip.text = element_text(size = 30, face="bold"),
        strip.background = element_rect(fill = "white"), plot.margin = margin(0.05, 0.05, 0.05, 0.05, "in")) +
  theme(axis.text.y = element_text(hjust = 1, face="bold", size=30), axis.text.x = element_text(hjust = 1, face="bold", size=30), axis.title.x=element_text(face="bold", size=40), axis.title.y=element_text(face="bold", size=40), legend.text=element_text(size=30),
        legend.key.size = unit(1.5, "cm"))

cc

ggsave("spikein_distributions.png", plot=cc, width = 28, height = 24 , units = "in", bg="white")


range(all_results_3Xembryos_spikeins_no4X_platycarpas$het_match_percentage)
quantile(all_results_3Xembryos_spikeins_no4X_platycarpas$het_match_percentage, probs =  c(0.025, 0.975)) 

range(all_results_3Xembryos_spikeins_no4X_platycarpas$hom_match_percentage)
range(all_results_3Xembryos_spikeins_no4X_platycarpas$total_match_percentage)
boxplot(all_results_3Xembryos_spikeins_no4X_platycarpas$het_match_percentage)  
boxplot(all_results_3Xembryos_spikeins_no4X_platycarpas$hom_match_percentage) 

# given what we've seen with the spikein ranges and apomixis control ranges, we will go with 92.5% for the het percentage.
# 95% should be a solid cutoff for hom.

#so freaking weird.. but the spikeins suggest we may need to adjust the quadrants for the tetraploid Mplatycarpas.
tetraploid_platycarpas_spikeins <- all_results_3Xembryos %>%
  dplyr::filter(str_detect(embryo, "_[abc]")) %>%
  dplyr::filter(mother %in% c("Mplatycarpa_588752", "Mplatycarpa_588847", "Mplatycarpa_589356", "Mplatycarpa_589415")) %>%
  dplyr::filter(het_sites_total > 50) %>%
  dplyr::filter(hom_sites_total > 50)

plot(tetraploid_platycarpas_spikeins$raw_cov, tetraploid_platycarpas_spikeins$total_match_percentage) #there's that correlation again. 
#it seems like increasing the coverage of the spikeins for the tetraploid M. platycarpas helps alleviate this
#noise issue a little bit. and it looks like the 6X 2023 embryos for M. platycarpa 589415 are shifted upward on the 
#y axis a few percentage points. That means if I sequenced the embryos to maybe 10, 12, 15X... perhaps some of this noise would disappear.
#but why are these the only genotypes with this noise? Where is it coming from? Random mistakes in the PCR-based libraries, maybe?
#That might make a bit of sense. You'd have to sequence deeper to get rid of some of the PCR amplified noise maybe. Making PCR free libraries
#might be able to answer this question someday. 

boxplot(tetraploid_platycarpas_spikeins$total_match_percentage)
boxplot(tetraploid_platycarpas_spikeins$het_match_percentage)
boxplot(tetraploid_platycarpas_spikeins$hom_match_percentage)
range(tetraploid_platycarpas_spikeins$het_match_percentage)
range(tetraploid_platycarpas_spikeins$hom_match_percentage)
range(tetraploid_platycarpas_spikeins$total_match_percentage)

quantile(tetraploid_platycarpas_spikeins$het_match_percentage, probs =  c(0.025, 0.975))
quantile(tetraploid_platycarpas_spikeins$hom_match_percentage, probs =  c(0.025, 0.975))
##following the convention that the range may actually be a little lower for lower coverage, we'll do 85 het, 82 hom

######### PLOTTING TIME - Figure 3 #####################################

#create a separate dataset to define a separate set of quadrants for the tetraploid Mplatycarpas
tetraploid_platycarpas <- all_results_3Xembryos_true_no_sims %>%
  dplyr::filter(mother %in% c("Mplatycarpa_588752", "Mplatycarpa_588847", "Mplatycarpa_589356", "Mplatycarpa_589415"))

#To combine the 3X and 6X results... we need to append a suffix to Mplatycarpa_589415 in one of the years. 2023 seed (the 6X embryos) will get a _2023 added to the name.
all_results_16moms_no_sims$embryo <- gsub("Mplatycarpa_589415", "Mplatycarpa_589415_2023", all_results_16moms_no_sims$embryo)
all_results_16moms_no_sims$mother <- gsub("Mplatycarpa_589415", "Mplatycarpa_589415_2023", all_results_16moms_no_sims$mother)

#then add this  to our tetraploid platycarpas dataframe so we can adjust the boundaries on the plot.
Mplatycarpa_589415_2023 <- all_results_16moms_no_sims %>%
  dplyr::filter(mother %in% c("Mplatycarpa_589415_2023")) #subset
tetraploid_platycarpas <- rbind(tetraploid_platycarpas, Mplatycarpa_589415_2023)

##combine the two big dataframes for plotting - but drop the controls. 
ALL <- rbind(all_results_16moms_no_sims, all_results_3Xembryos_true_no_sims) %>%
  dplyr::filter(!mother %in% c("Mspectabilis_588917", "Masiatica_594099", "Mhupehensis_633818", "Msargentii_589400"))

#order levels for plotting. 
ALL$mother <- factor(ALL$mother, 
                            levels=c("Madstringens_588898","Mbaccata_588907","Mbaccata_613807","Mbrevipes_589170","Mfloribunda_589827", 
                                     "Mhupehensis_589756","MIndianSummer_589733","MMaryPotter_588983","MProfusion_589449","Mtoringoides_588930",
                                     "MMcClintockGrimes_589124","MGoldgelb_589458","Mkansuensis_588944",
                                     "Mrobusta_589383","Mrobusta_588825",
                                     "Mtransitoria_633806","MHenryFDupont_589732","Mxzumi_589840","Mhalliana_590029","Mhalliana_589013",
                                     "MRobinson_589455" ,"Mzhaojiaoensis_633816","Mprunifolia_588914","Mmicromalus_588976", "Mmicromalus_594092",
                                     "Msikkimensis_589599","Msikkimensis_589750", "Msikkimensis_613912",
                                     "MRalphShay_589734", "Mhartwiggii_588757", "Mmagdeburgensis_588959", "Mioensis_590008",
                                     "Mplatycarpa_588752", "Mplatycarpa_588847", "Mplatycarpa_589356", "Mplatycarpa_589415_2023", "Mplatycarpa_589415", "Mplatycarpa_589198",
                                     "Mcoronaria_590014", "Mcoronaria_590000", "Mcoronaria_589344", "Mcoronaria_589988", "Mcoronaria_589977",
                                     "Mrockii_589279", "MYellowAutumnCrabapple_588922", "Mtoringo_633814", "Mtoringo_613945", "Mtoringo_613932", "Mtoringo_590101",
                                     "Msargentii_589372", "Mhybrid_589161", "Mhybrid_589421"))

everybody_else <- ALL %>%
  dplyr::filter(!mother %in% c("Mplatycarpa_588752", "Mplatycarpa_588847", "Mplatycarpa_589356", "Mplatycarpa_589415", "Mplatycarpa_589415_2023"))

da_big_one <- ggplot(ALL, aes(x=het_match_percentage, y=hom_match_percentage, color=mother)) +
  geom_point(shape=17, size=4)+
  scale_color_manual(values=c("Madstringens_588898"="pink","Mbaccata_588907"="cyan3","Mbaccata_613807"="firebrick","Mbrevipes_589170"="lightblue","Mfloribunda_589827"="lavender","Mhupehensis_589756"="steelblue","MIndianSummer_589733"="grey50",
                              "MMaryPotter_588983"="bisque","MProfusion_589449"="orange","Mtoringoides_588930"="gold", "MMcClintockGrimes_589124"="azure2","MGoldgelb_589458"="darkolivegreen1","Mkansuensis_588944"="tomato3",
                              "Mrobusta_589383"="#990033", "Mrobusta_588825"="#CC00CC", 
                              "Mtransitoria_633806"="#CC0000","MHenryFDupont_589732"="#3393FF","Mxzumi_589840"="#CCFFCC","Mhalliana_590029"="coral2","Mhalliana_589013"="#FF728B",
                              "MRobinson_589455"="#660033","Mzhaojiaoensis_633816"="#99CCFF","Mprunifolia_588914"="#FFCCCC","Mmicromalus_588976"="#FFFF00", "Mmicromalus_594092"="#CC3333",
                              "Msikkimensis_589599"="#006600","Msikkimensis_589750"="#99FF99", "Msikkimensis_613912"="#CC6666", 
                              "MRalphShay_589734"="#FF3399", "Mhartwiggii_588757"="#CC6600", "Mmagdeburgensis_588959"="#00FFFF", "Mioensis_590008"="#6699CC",
                              "Mplatycarpa_588752"="#000033", "Mplatycarpa_588847"="#0033CC", "Mplatycarpa_589356"="#0066FF", "Mplatycarpa_589415_2023"="hotpink4", "Mplatycarpa_589415"="hotpink4", "Mplatycarpa_589198"="deeppink3",
                              "Mcoronaria_590014"="#FFCC66","Mcoronaria_590000"="#CC9FF9","Mcoronaria_589344"="#FFCCFF", "Mcoronaria_589988"="#669996", "Mcoronaria_589977"="#330000",
                              "Mrockii_589279"="#FF99CC", "MYellowAutumnCrabapple_588922"="#FFCC33", "Mtoringo_633814"="#66CCFF", "Mtoringo_613945"="#FF9966", "Mtoringo_613932"="#009933", "Mtoringo_590101"="#CC0033", 
                              "Msargentii_589372"="#00CC33", "Mhybrid_589161"="#FFCC99", "Mhybrid_589421"="#FF0099"))+
  ggtitle("Genetic Similarity")+
  #labs(subtitle="(No Error Adjustment)")+
  xlab("Mom Het sites (% matching sites)")+
  ylab("Mom Hom sites (% matching sites)")+
  scale_y_continuous(limits=c(5,100), breaks = seq(5, 100, by = 5))+
  scale_x_continuous(limits=c(25,100), breaks = seq(25, 100, by = 5))+
  facet_wrap(~mother, nrow=5,ncol=11)+
  #stat_smooth(method="lm")+
  coord_fixed()+
  #scale_x_continuous(limits=c(60,100))+
  geom_vline(data=everybody_else, aes(xintercept=92.5), linetype="dashed",color="black", size=0.8)+
  geom_hline(data=everybody_else, aes(yintercept=95), linetype="dashed",color="black", size=0.8)+
  geom_vline(data=tetraploid_platycarpas, aes(xintercept=85), linetype="dashed", color="hotpink3", size=0.8)+
  geom_hline(data=tetraploid_platycarpas, aes(yintercept=82), linetype="dashed", color="hotpink3", size=0.8)+
  theme_bw()+
  theme(plot.title=element_text(hjust = 0.5), plot.subtitle=element_text(hjust=0.5, size=20), panel.spacing = unit(1, "lines"), strip.text = element_text(size = 12, face="bold"),
        strip.background = element_rect(fill = "white"), plot.margin = margin(0.5, 0.05, 0.05, 0.05, "in")) +
  theme(axis.text.x = element_text(face="bold", size=10), axis.text.y = element_text(hjust = 1, face="bold", size=15), axis.title.x=element_text(face="bold", size=30), axis.title.y=element_text(face="bold", size=30), plot.title=element_text(size=30, face = "bold", hjust=0.5), legend.position = "none")

da_big_one #there is one Mbaccata 588907 embryo not on the plot because it has an extremely low het similarity. It is presumably a selfed embryo in Q1. 

ggsave("ALL_no_ctls.png", plot=da_big_one, width = 40, height = 24 , units = "in", bg="white")


############################################ Creating the massive data frame for supplementary table 2 ######################################
#all_results_16moms = 6X results, including sims
#raw_cov_all_3X_sims_and_all_incl_spikeins = 3X results, including sims
str(all_results_16moms)
str(all_results_3Xembryos)

#Should add seed year column to both before combining. 
all_results_16moms <- all_results_16moms %>%
  mutate(seed_year = case_when(
    grepl("_sim", embryo) ~ "simulated",
    grepl("_[abcde]$", embryo) ~ "spikein",
    TRUE ~ "2023"
  ))

all_results_3Xembryos <- all_results_3Xembryos %>%
  mutate(seed_year = case_when(
    grepl("_sim", embryo) ~ "simulated",
    grepl("_[abcde]$", embryo) ~ "spikein",
    TRUE ~ "2024"
  ))

#ready to be concatenated
SUPPTABLE2 <- rbind(all_results_16moms, all_results_3Xembryos)
#add quadrant information. Special boundaries required for tetraploid M. platycarpas. 
#"Mplatycarpa_588752", "Mplatycarpa_588847", "Mplatycarpa_589356", "Mplatycarpa_589415"

SUPPTABLE2 <- SUPPTABLE2 %>%
  mutate(
    quadrant = case_when(
      # Exclude simulated and spikein samples first
      grepl("_sim", embryo) | grepl("_[abcde]$", embryo) ~ NA_character_,
      
      # Too few sites rule (applies if not simulated/spikein)
      het_sites_total < 50 | hom_sites_total < 50 ~ "too few sites to consider",
      
      # Special mother samples (custom thresholds)
      mother %in% c("Mplatycarpa_588752", "Mplatycarpa_588847", 
                    "Mplatycarpa_589356", "Mplatycarpa_589415") &
        het_match_percentage < 85 & hom_match_percentage > 82 ~ "Q1",
      
      mother %in% c("Mplatycarpa_588752", "Mplatycarpa_588847", 
                    "Mplatycarpa_589356", "Mplatycarpa_589415") &
        het_match_percentage > 85 & hom_match_percentage > 82 ~ "Q2",
      
      mother %in% c("Mplatycarpa_588752", "Mplatycarpa_588847", 
                    "Mplatycarpa_589356", "Mplatycarpa_589415") &
        het_match_percentage < 85 & hom_match_percentage < 82 ~ "Q3",
      
      mother %in% c("Mplatycarpa_588752", "Mplatycarpa_588847", 
                    "Mplatycarpa_589356", "Mplatycarpa_589415") &
        het_match_percentage > 85 & hom_match_percentage < 82 ~ "Q4",
      
      # Regular samples (standard thresholds)
      het_match_percentage < 92.5 & hom_match_percentage > 95 ~ "Q1",
      het_match_percentage > 92.5 & hom_match_percentage > 95 ~ "Q2",
      het_match_percentage < 92.5 & hom_match_percentage < 95 ~ "Q3",
      het_match_percentage > 92.5 & hom_match_percentage < 95 ~ "Q4",
      
      # fallback
      TRUE ~ NA_character_
    )
  )
SUPPTABLE2$seed_year <- as.factor(SUPPTABLE2$seed_year)
table(SUPPTABLE2$seed_year)
#write this to a big fat table...
write.csv(SUPPTABLE2, file="Supplementary_Table2.csv", row.names=FALSE)

########## checking association between depth and variant calling error rate. 

ALL_sims <- SUPPTABLE2[grepl("sim", SUPPTABLE2$embryo), ] %>%
  dplyr::filter(het_sites_total > 50) %>%
  dplyr::filter(hom_sites_total > 50)
  
plot(ALL_sims$raw_cov, ALL_sims$total_match_percentage)
range(ALL_sims$total_match_percentage)
ALL_sims$raw_cov_factor <- as.factor(ALL_sims$raw_cov)
sim_model <- lm(total_match_percentage ~ raw_cov_factor, data=ALL_sims)
summary(sim_model)
Anova(sim_model, type="III")

############################################### Lastly, beautify Vasek's example flow spectra (Figure 2B) ###########################################
Mhupehensis_633818 <- read.table(file="Vasek_example_spectra/hupehensis 633818 AP 23 7 CX100_023.txt", header=TRUE, sep="\t")
colnames(Mhupehensis_633818) <- c("Event","Value")
Msargentii_589400 <- read.table(file="Vasek_example_spectra/sargentii 589400 AP 24 7 CX100_049.txt", header=TRUE, sep="\t")
colnames(Msargentii_589400) <- c("Event","Value")
Mspectabilis_588917 <- read.table(file="Vasek_example_spectra/spectabilis 588917 SE 24 1 CX100.txt", header=TRUE, sep="\t")
colnames(Mspectabilis_588917) <- c("Event","Value")
Masiatica_594099 <- read.table(file="Vasek_example_spectra/asiatica 594099 SE 23 2 CX100_120.txt", header=TRUE, sep="\t")
colnames(Masiatica_594099) <- c("Event","Value")

combined_flow <- bind_rows(
  Mspectabilis_588917 = Mspectabilis_588917,
  Masiatica_594099 = Masiatica_594099,
  Mhupehensis_633818 = Mhupehensis_633818,
  Msargentii_589400 = Msargentii_589400,
  .id = "Mom"
)

combined_flow$Mom <- factor(combined_flow$Mom, levels=c("Mspectabilis_588917", "Masiatica_594099", "Mhupehensis_633818", "Msargentii_589400"))

flow <- ggplot(combined_flow, aes(x = Value, fill = Mom)) +
  geom_histogram(binwidth = 1) +
  scale_fill_manual(values=c("Mspectabilis_588917"="hotpink2","Masiatica_594099"="goldenrod","Mhupehensis_633818"="chocolate4","Msargentii_589400"="forestgreen"))+
  ggtitle("Flow Cytometry Spectra") +
  xlab("Intensity Value")+
  ylab("Counts")+
  facet_wrap(~Mom, nrow=1, ncol=4)+
  scale_x_continuous(limits=c(75,1100), breaks = seq(100, 1100, by = 100), expand = c(0, 0))+
  scale_y_continuous(limits=c(0,125), breaks = seq(0, 120, by = 10), expand = c(0, 0))+
  theme_bw()+
  theme(plot.title=element_text(hjust = 0.5, size=30, face="bold"), panel.spacing = unit(1, "lines"), strip.text = element_text(size = 20, face="bold"),
        strip.background = element_rect(fill = "white"), plot.margin = margin(0.5, 0.12, 0.05, 0.05, "in"), axis.text.y = element_text(hjust = 1, face="bold", size=20), axis.text.x = element_text(hjust = 1, face="bold", size=20, angle=45, vjust=1), axis.title.x=element_text(face="bold", size=30), axis.title.y=element_text(face="bold", size=30), legend.text=element_text(size=30),
        legend.position = "none", axis.ticks=element_line(linewidth=1), axis.ticks.length = unit(0.3, "cm"))
flow
ggsave("flow_plots.png", plot=flow, width = 28, height = 8 , units = "in", bg="white")


