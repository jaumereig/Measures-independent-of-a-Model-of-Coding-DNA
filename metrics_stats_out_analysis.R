library("ggplot2")
library("dplyr")
library("cowplot")
library("gridExtra")
#############
#  FOURIER  #
#############
fourier.stats <- read.table("~/Escritorio/Internship_CRG/BishopBook/Measures_Jaume/gffcompare_out/Fourier_stats_out.header", header = TRUE)
fourier.stats$Peak_threshold <- fourier.stats$Peak_threshold/100
head(fourier.stats)

Wlen <- as.factor(fourier.stats$WLen)
Peak_threshold <- as.factor(fourier.stats$Peak_threshold)
Peak_width <- as.factor(fourier.stats$Peak_width)
fourier.stats$mean.sp.sn <- (fourier.stats$SNn + fourier.stats$SPn)/2

f1 <- ggplot(fourier.stats, aes(x=SNn, y=SPn))+
  geom_point(stat="identity", aes(color = as.factor(Peak_threshold)))+
  theme_bw() +
  ggtitle("Fourier Specificity~Sensitivity by Peak Height") +
  theme(plot.title = element_text(hjust = 0.5)) + facet_grid(~Peak_width)
f1
f2 <- ggplot(fourier.stats, aes(x=mean.sp.sn))+
  geom_boxplot(aes(color = as.factor(Peak_threshold)))+
  theme_bw() +
  ggtitle("Fourier Specificity~Sensitivity by Peak Height") +
  theme(plot.title = element_text(hjust = 0.5)) + facet_grid(~WLen)
f2
f3 <- ggplot(fourier.stats, aes(x=SNn, y=SPn))+
  geom_point(stat="identity", aes(color = as.factor(Peak_width)))+
  theme_bw() +
  ggtitle("Fourier Specificity~Sensitivity by Peak Width") +
  theme(plot.title = element_text(hjust = 0.5))
f3
f4 <- ggplot(fourier.stats, aes(x=SNn, y=SPn))+
  geom_point(stat="identity", aes(color = as.factor(Wlen)))+
  theme_bw() +
  ggtitle("Fourier Specificity~Sensitivity by Window Length")+
  theme(plot.title = element_text(hjust = 0.5))
f4


#############
# ASYMMETRY #
#############
asym.stats <- read.table("~/Escritorio/Internship_CRG/BishopBook/Measures_Jaume/gffcompare_out/Assymetry_stats_out.header", header = TRUE)
asym.stats$Peak_threshold <- asym.stats$Peak_threshold/100
head(asym.stats)
Wlen <- as.factor(asym.stats$WLen)
asym.stats$mean.sp.sn <- (asym.stats$SNn + asym.stats$SPn)/2

a1 <- ggplot(asym.stats, aes(x=SNn, y=SPn))+
  geom_point(stat="identity", aes(color = as.factor(Peak_threshold)), size = 2)+
  theme_bw() +
  ggtitle("Position Asymmetry Specificity~Sensitivity by Peak Height") +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
a1
a2 <- ggplot(asym.stats, aes(x=SNn, y=SPn))+
  geom_point(stat="identity", aes(color = as.factor(Peak_width)), size = 2)+
  theme_bw() +
  ggtitle("Position Asymmetry Specificity~Sensitivity by Peak Width") +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
a2
a3 <- ggplot(asym.stats, aes(x=SNn, y=SPn))+
  geom_point(stat="identity", aes(color = as.factor(Wlen)), size = 2)+
  theme_bw() +
  ggtitle("Position Asymmetry Specificity~Sensitivity by Window Length")+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
a3
grid.arrange(a1, a2, a3, ncol=2, nrow = 2)
#############
#    AMI    #
#############
ami.stats <- read.table("~/Escritorio/Internship_CRG/BishopBook/Measures_Jaume/gffcompare_out/AMI_stats_out.header", header = TRUE)
ami.stats$Peak_threshold <- ami.stats$Peak_threshold/100
ami.stats$mean.sp.sn <- (ami.stats$SNn + ami.stats$SPn)/2

head(ami.stats)

ami.stats.filter <- read.table("~/Escritorio/Internship_CRG/BishopBook/Measures_Jaume/gffcompare_out/AMI_stats_filter_out.header", header = TRUE)
ami.stats.filter$Peak_threshold <- ami.stats.filter$Peak_threshold/100
ami.stats.filter$mean.sp.sn <- (ami.stats.filter$SNn + ami.stats.filter$SPn)/2

head(ami.stats.filter)

ami.stats.filter <- read.table("~/Escritorio/Internship_CRG/BishopBook/Measures_Jaume/gffcompare_out/AMI_stats_Test_out.header", header = TRUE)
ami.stats.filter$Peak_threshold <- ami.stats.filter$Peak_threshold/100
ami.stats.filter$mean.sp.sn <- (ami.stats.filter$SNn + ami.stats.filter$SPn)/2

head(ami.stats.filter)

#########################
#   COMPARING MEASURES  #
#########################
measures <- data.frame(Program = c("Asymmetry Position", "Average Mutual Info", "Fourier", "Matching CDS"),
                       Specificity = c(32.1, 14.9, 35.8, 50.5), Sensitivity = c(1.9, 18.8, 2.6, 39.4),
                       Average = c(17, 16.85, 19.20, 44.95))

library(reshape2)

#ggplot(data = measures, aes(x = Specificity, y = Sensitivity, fill = Program)) + geom_bar(stat = 'identity', position = 'dodge')

long.measures <- melt(measures, id.vars = c("Program"))
head(long.measures)

p1 <- ggplot(data = long.measures, aes(x = Program, y = value, label = value)) + 
  geom_bar(stat = 'identity', position = 'dodge', aes(fill = variable)) +
  ylab("Accuracy") +
  theme_bw() +
  ylim(0,100) +
  ggtitle("Gffcompare accuracy at base level by different methods in human chromosome 22", ) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
p1

measures <- data.frame(Program = c("Matching CDS", "Matching ORF"),
                       Specificity = c(50.5, 65.6), Sensitivity = c(39.4, 28.5),
                       Average = c(44.95, 47.05))
long.measures <- melt(measures, id.vars = c("Program"))
head(long.measures)

p2 <- ggplot(data = long.measures, aes(x = Program, y = value, label = value)) + 
  geom_bar(stat = 'identity', position = 'dodge', aes(fill = variable)) +
  ylab("Accuracy") +
  theme_bw() +
  ylim(0,100) +
  ggtitle("Gffcompare accuracy at base level by different methods in human chromosome 22", ) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
p2

# 
df <- read.table("~/Escritorio/Internship_CRG/BishopBook/Measures_Jaume/protein_matches/matches_orf.txt", header = TRUE)
df$num_match <- log2(x = df$num_match)
long.measures <- melt(df, id.vars = c("match_score"))
p3 <- ggplot(data = long.measures, aes(x = match_score, y = value, label = value)) + 
  geom_bar(stat = 'identity', position = 'dodge', aes(fill = variable)) +
  ylab("Accuracy") +
  theme_bw() +
  ylim(0,100) +
  ggtitle("Gffcompare accuracy at base level for different protein match scores in human chromosome 22", ) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
p3
