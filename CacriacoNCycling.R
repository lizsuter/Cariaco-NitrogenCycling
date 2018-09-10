#========================================================#
# Project for Analyzing Xander Results #
# Metagenome and metatranscriptome results for Cariaco
# N Cycling Paper (focusing on nirS and nirK)
#========================================================#
# Script started Jul 23rd 2018
# For notes on using Xander to analyze these data, see lab notebook
# "CariacoNcyclingNotes.md"

rm(list = ls())

# Load packages every time
#library(stats)
library(tidyverse)
library(vegan)
library(ggthemes)

# Make a variable for sample name, which is in file name for all results
samples_metaG <- c("D2a143A", "D2a237B", "D2b143A", "D2b237B", "D3a103A", "D3a234B", "D3b103A",	"D3b234B", 
"D2a143B", "D2a247A", "D2b143B", "D2b247A", "D3a103B", "D3a295A", "D3b103B", "D3b295A", 
"D2a200A", "D2a247B", "D2b200A", "D2b247B", "D3a198A", "D3a295B", "D3b198A", "D3b295B",
"D2a200B", "D2a267A", "D2b200B", "D2b267A", "D3a198B", "D3a314A", "D3b198B", "D3b314A",
"D2a237A", "D2a267B", "D2b237A", "D2b267B", "D3a234A", "D3a314B", "D3b234A", "D3b314B")

samples_metaT <- c("R2a103A", "R2a198A", "R2a234A", "R2a295A", "R2a314A", "R2b198A", "R2b234A", "R2b295A", "R2b314A", "R3b103A",
                   "R2a103B", "R2a198B", "R2a234B", "R2a295B", "R2a314B", "R2b198B", "R2b234B", "R2b295B", "R2b314B", "R3b103B",
                   "R2a148A", "R2a148B", "R2a200A", "R2a200B", "R2a237A", "R2a237B", "R2a247A", "R2a247B", "R2a267A", "R2a267B", 
                   "R2b148B", "R2b148B", "R2b200A", "R2b200B", "R2b237A", "R2b237B", "R2b247A", "R2b247B", "R2b267A", "R2b267B")

# Write empty "total_reads" matrix for each of 3 genes. This will be filled in:
total_reads_metaG <- tibble("samples_metaG" = samples_metaG,
                      "recA" = (1:length(samples_metaG))*0, #empty array for now
                      "nirS" = (1:length(samples_metaG))*0, 
                      "nirK" = (1:length(samples_metaG))*0)

total_reads_metaT <- tibble("samples_metaT" = samples_metaG,
                            "recA" = (1:length(samples_metaG))*0, #empty array for now
                            "nirS" = (1:length(samples_metaG))*0, 
                            "nirK" = (1:length(samples_metaG))*0)


# pull out # of reads from coverage.txt file for each gene:

# First recA. Start with metaG samples
# for recA this is straightforward because recA is in every samples_metaG
recA_reads_metaG <- tibble("samples_metaG" = (1:length(samples_metaG))*0, #empty column for now,
                      "recA" = (1:length(samples_metaG))*0) #empty column for now)
for(i in 1:length(samples_metaG)){
  recA_reads_metaG[i,1] <- samples_metaG[i]
  temp <- read_delim(paste("Xander_results/",samples_metaG[i],"_recA_45_coverage.txt", sep = ""),delim = "  \t", col_names = F)
  recA_reads_metaG[i,2] <- temp[1,3]
}

# Now recA from metaT samples
recAsamples_metaT <- c("R2a148B", "R2a198B", "R2a200B", "R2a234B", "R2a237B", "R2a247A", "R2a247B", 
                       "R2a267A", "R2a267B", "R2a295B", "R2a314B", "R2b198A", "R2b198B", "R2b234B", 
                       "R2b267A", "R2b295B",  "R2b314B", "R3b103A")
recA_reads_metaT <- tibble("samples_metaT" = (1:length(recAsamples_metaT))*0, #empty column for now,
                           "recA" = (1:length(recAsamples_metaT))*0) #empty column for now)
for(i in 1:length(recAsamples_metaT)){
  recA_reads_metaT[i,1] <- recAsamples_metaT[i]
  temp <- read_delim(paste("Xander_results/",recAsamples_metaT[i],"_recA_45_coverage.txt", sep = ""),delim = "  \t", col_names = F)
  recA_reads_metaT[i,2] <- temp[1,3]
}


# NirS
# for nirS, rewrite the "sample" vector with only the samples that had nirS 

# First metaG samples
nirSsamples_metaG <- c("D2a237B", "D2b237B", "D3a234B", "D3b234B", 
            "D2a247A", "D2b247A", "D3a295A", "D3b295A", 
            "D2a247B", "D2b200A", "D2b247B", "D3a198A", "D3a295B", "D3b198A", "D3b295B",
            "D2a267A", "D2b200B", "D2b267A", "D3a314A", "D3b314A",
            "D2a237A", "D2a267B", "D2b237A", "D2b267B", "D3a234A", "D3a314B", "D3b234A", "D3b314B")


nirS_reads_metaG <- tibble("samples_metaG" = (1:length(nirSsamples_metaG))*0, #empty column for now,
                     "nirS" = (1:length(nirSsamples_metaG))*0) #empty column for now)
for(i in 1:length(nirSsamples_metaG)){
  nirS_reads_metaG[i,1] <- nirSsamples_metaG[i]
  temp <- read_delim(paste("Xander_results/",nirSsamples_metaG[i],"_nirS_45_coverage.txt", sep = ""),delim = "  \t", col_names = F)
  nirS_reads_metaG[i,2] <- temp[1,3]
}


# And same for metaT samples. Modify the "Samples" vector first
nirSsamples_metaT <- c("R2a200B", "R2a237B", "R2a247B", 
                       "R2a314B","R2b314B", "R2b237B","R2b314A")
nirS_reads_metaT <- tibble("samples_metaT" = (1:length(nirSsamples_metaT))*0, #empty column for now,
                           "nirS" = (1:length(nirSsamples_metaT))*0) #empty column for now)
for(i in 1:length(nirSsamples_metaT)){
  nirS_reads_metaT[i,1] <- nirSsamples_metaT[i]
  temp <- read_delim(paste("Xander_results/",nirSsamples_metaT[i],"_nirS_45_coverage.txt", sep = ""),delim = "  \t", col_names = F)
  nirS_reads_metaT[i,2] <- temp[1,3]
}



# NirK
# rewrite the "sample" vector with only the samples that had nirK

nirKsamples_metaG <- c("D2a143B","D2a237B","D2b143B","D2b200B","D2b237B","D3a198A",
            "D3a198B","D3a234A","D3a234B", "D3b198B","D3b234A","D3b234B")

nirK_reads_metaG <- tibble("samples_metaG" = (1:length(nirKsamples_metaG))*0, #empty column for now,
                     "nirK" = (1:length(nirKsamples_metaG))*0) #empty column for now)
for(i in 1:length(nirKsamples_metaG)){
  nirK_reads_metaG[i,1] <- nirKsamples_metaG[i]
  temp <- read_delim(paste("Xander_results/",nirKsamples_metaG[i],"_nirK_45_coverage.txt", sep = ""),delim = "  \t", col_names = F)
  nirK_reads_metaG[i,2] <- temp[1,3]
}

# And metaT for nirK- only one sample had it
nirKsamples_metaT <- c("R2a314B")
nirK_reads_metaT <- tibble("samples_metaT" = (1:length(nirKsamples_metaT))*0, #empty column for now,
                           "nirK" = (1:length(nirKsamples_metaT))*0) #empty column for now)
for(i in 1:length(nirKsamples_metaT)){
  nirK_reads_metaT[i,1] <- nirKsamples_metaT[i]
  temp <- read_delim(paste("Xander_results/",nirKsamples_metaT[i],"_nirK_45_coverage.txt", sep = ""),delim = "  \t", col_names = F)
  nirK_reads_metaT[i,2] <- temp[1,3]
}


#================
# Match reads for all 3 genes into 1 matrix

# write empty "total_reads" matrix for each of 3 genes and fill with # of reads
total_reads_metaG <- tibble("samples_metaG" = samples_metaG)
total_reads_metaG <- left_join(total_reads_metaG, nirK_reads_metaG, by = "samples_metaG")
total_reads_metaG <- left_join(total_reads_metaG, nirS_reads_metaG, by = "samples_metaG")
total_reads_metaG <- left_join(total_reads_metaG, recA_reads_metaG, by = "samples_metaG")

total_reads_metaT <- tibble("samples_metaT" = samples_metaT)
total_reads_metaT <- left_join(total_reads_metaT, nirK_reads_metaT, by = "samples_metaT")
total_reads_metaT <- left_join(total_reads_metaT, nirS_reads_metaT, by = "samples_metaT")
total_reads_metaT <- left_join(total_reads_metaT, recA_reads_metaT, by = "samples_metaT")


# replace NA with 0 and remove unecessary variables 
total_reads_metaG[is.na(total_reads_metaG)] <- 0
total_reads_metaT[is.na(total_reads_metaT)] <- 0
rm(nirK_reads_metaG, nirS_reads_metaG, recA_reads_metaG, nirK_reads_metaT, nirS_reads_metaT, recA_reads_metaT, temp)

# make matrix for nirSrecA and nirKrecA ratios
read_ratios_metaG <- tibble("samples_metaG" = samples_metaG, #empty column for now,
                      "nirSrecA_DNA" = (1:length(samples_metaG))*0, 
                      "nirKrecA_DNA" = (1:length(samples_metaG))*0)
read_ratios_metaT <- tibble("samples_metaT" = samples_metaT, #empty column for now,
                            "nirSrecA_RNA" = (1:length(samples_metaT))*0, 
                            "nirKrecA_RNA" = (1:length(samples_metaT))*0)


# Fill in with calculation of ratio (#for tibbles, the entries are characters so you have 
# to tell it "as.numeric" in order to do a normal binary operation)
for(i in 1:length(samples_metaG)){
  if (total_reads_metaG$nirS[i]>0){
    read_ratios_metaG$`nirSrecA_DNA`[i] <- as.numeric(total_reads_metaG$nirS[i]) / as.numeric(total_reads_metaG$recA[i])
  }  else {
    read_ratios_metaG$`nirSrecA_DNA`[i] <- 0
  }
}
for(i in 1:length(samples_metaG)){
  if (total_reads_metaG$nirK[i]>0){
    read_ratios_metaG$`nirKrecA_DNA`[i] <- as.numeric(total_reads_metaG$nirK[i]) / as.numeric(total_reads_metaG$recA[i])
  }  else {
    read_ratios_metaG$`nirKrecA_DNA`[i] <- 0
  }
}

for(i in 1:length(samples_metaT)){
  if (total_reads_metaT$nirS[i]>0){
    read_ratios_metaT$`nirSrecA_RNA`[i] <- as.numeric(total_reads_metaT$nirS[i]) / as.numeric(total_reads_metaT$recA[i])
  }  else {
    read_ratios_metaT$`nirSrecA_RNA`[i] <- 0
  }
}
for(i in 1:length(samples_metaT)){
  if (total_reads_metaT$nirK[i]>0){
    read_ratios_metaT$`nirKrecA_RNA`[i] <- as.numeric(total_reads_metaT$nirK[i]) / as.numeric(total_reads_metaT$recA[i])
  }  else {
    read_ratios_metaT$`nirKrecA_RNA`[i] <- 0
  }
}

# Make into data frame
read_ratios_metaG <- as.data.frame(read_rread_ratios_metaGatios)
rownames(read_ratios_metaG) <- read_ratios_metaG[,1]

read_ratios_metaT <- as.data.frame(read_ratios_metaT)
rownames(read_ratios_metaT) <- read_ratios_metaT[,1]


##### Commented out the following section where I calculated the average ratio of nirS:recA
# and nirK:recA for duplicate samples. I decided not to do that so I could plot the replicates.
# Also had to modify metadata sheet so that it would attach to individual samples, not averaged
# sample.


# # Average duplicate samples and put in new matrix
# samples_metaG <- c("D2a143A", "D2a237B", "D2b143A", "D2b237B", "D3a103A", "D3a234B", "D3b103A",	"D3b234B", 
#             "D2a143B", "D2a247A", "D2b143B", "D2b247A", "D3a103B", "D3a295A", "D3b103B", "D3b295A", 
#             "D2a200A", "D2a247B", "D2b200A", "D2b247B", "D3a198A", "D3a295B", "D3b198A", "D3b295B",
#             "D2a200B", "D2a267A", "D2b200B", "D2b267A", "D3a198B", "D3a314A", "D3b198B", "D3b314A",
#             "D2a237A", "D2a267B", "D2b237A", "D2b267B", "D3a234A", "D3a314B", "D3b234A", "D3b314B")
# samples_metaG_labels <- c("D143A", "D237B", "D103A", "D234B", 
#             "D143B", "D247A", "D103B", "D295A", 
#             "D200A", "D247B", "D198A", "D295B", 
#             "D200B", "D267A", "D198B", "D314A",
#             "D237A", "D267B", "D234A", "D314B")
# 
# # create empty matrix
# mean_read_ratios <- tibble("samples_metaG" = samples_metaG_labels,
#                     "nirSrecA" = (1:length(samples_metaG_labels))*NA,
#                     "nirKrecA" = (1:length(samples_metaG_labels))*NA) #empty column for now)
# 
# # and fill with the mean of the ratios of the reads to recA
# read_ratios <- as.data.frame(read_ratios)# need to convert to data frame in order to call rows by name
# rownames(read_ratios) <- read_ratios[,1]
# 
# mean_read_ratios <- as.data.frame(mean_read_ratios)
# rownames(mean_read_ratios) <- mean_read_ratios[,1]
# 
# 
# # First fill in the averages of nirSrecA. Work from "samples_metaG" variable above that gave all samples 
# # where nirS was detected
# # "D2a237B", "D2b237B", "D3a234B", "D3b234B", 
# # "D2a247A", "D2b247A", "D3a295A", "D3b295A", 
# # "D2a247B", "D2b200A", "D2b247B", "D3a198A", "D3a295B", "D3b198A", "D3b295B",
# # "D2a267A", "D2b200B", "D2b267A", "D3a314A", "D3b314A",
# # "D2a237A", "D2a267B", "D2b237A", "D2b267B", "D3a234A", "D3a314B", "D3b234A", "D3b314B"
# mean_read_ratios["D237B", "nirSrecA"] = sum(read_ratios[grepl(glob2rx("D*237B"), read_ratios$samples_metaG), "nirSrecA"])/2
# mean_read_ratios["D234B", "nirSrecA"] = sum(read_ratios[grepl(glob2rx("D*234B"), read_ratios$samples_metaG), "nirSrecA"])/2
# mean_read_ratios["D247A", "nirSrecA"] = sum(read_ratios[grepl(glob2rx("D*247A"), read_ratios$samples_metaG), "nirSrecA"])/2
# mean_read_ratios["D295A", "nirSrecA"] = sum(read_ratios[grepl(glob2rx("D*295A"), read_ratios$samples_metaG), "nirSrecA"])/2
# mean_read_ratios["D247B", "nirSrecA"] = sum(read_ratios[grepl(glob2rx("D*247B"), read_ratios$samples_metaG), "nirSrecA"])/2
# mean_read_ratios["D200A", "nirSrecA"] = sum(read_ratios[grepl(glob2rx("D*200A"), read_ratios$samples_metaG), "nirSrecA"])/2
# mean_read_ratios["D198A", "nirSrecA"] = sum(read_ratios[grepl(glob2rx("D*198A"), read_ratios$samples_metaG), "nirSrecA"])/2
# mean_read_ratios["D295B", "nirSrecA"] = sum(read_ratios[grepl(glob2rx("D*295B"), read_ratios$samples_metaG), "nirSrecA"])/2
# mean_read_ratios["D267A", "nirSrecA"] = sum(read_ratios[grepl(glob2rx("D*267A"), read_ratios$samples_metaG), "nirSrecA"])/2
# mean_read_ratios["D200B", "nirSrecA"] = sum(read_ratios[grepl(glob2rx("D*200B"), read_ratios$samples_metaG), "nirSrecA"])/2
# mean_read_ratios["D314A", "nirSrecA"] = sum(read_ratios[grepl(glob2rx("D*314A"), read_ratios$samples_metaG), "nirSrecA"])/2
# mean_read_ratios["D237A", "nirSrecA"] = sum(read_ratios[grepl(glob2rx("D*237A"), read_ratios$samples_metaG), "nirSrecA"])/2
# mean_read_ratios["D267B", "nirSrecA"] = sum(read_ratios[grepl(glob2rx("D*267B"), read_ratios$samples_metaG), "nirSrecA"])/2
# mean_read_ratios["D234A", "nirSrecA"] = sum(read_ratios[grepl(glob2rx("D*234A"), read_ratios$samples_metaG), "nirSrecA"])/2
# mean_read_ratios["D314B", "nirSrecA"] = sum(read_ratios[grepl(glob2rx("D*314B"), read_ratios$samples_metaG), "nirSrecA"])/2
# 
# 
# # Now nirKrecA
# # "D2a237B", "D2b237B", "D3a234B", "D3b234B", "D2a143B", "D2b143B",
# # "D3a198A", "D2b200B", "D3a198B", "D3b198B", "D3a234A", "D3b234A"
# mean_read_ratios["D237B", "nirKrecA"] = sum(read_ratios[grepl(glob2rx("D*237B"), read_ratios$samples_metaG), "nirKrecA"])/2
# mean_read_ratios["D234B", "nirKrecA"] = sum(read_ratios[grepl(glob2rx("D*234B"), read_ratios$samples_metaG), "nirKrecA"])/2
# mean_read_ratios["D143B", "nirKrecA"] = sum(read_ratios[grepl(glob2rx("D*143B"), read_ratios$samples_metaG), "nirKrecA"])/2
# mean_read_ratios["D198A", "nirKrecA"] = sum(read_ratios[grepl(glob2rx("D*198A"), read_ratios$samples_metaG), "nirKrecA"])/2
# mean_read_ratios["D200B", "nirKrecA"] = sum(read_ratios[grepl(glob2rx("D*200B"), read_ratios$samples_metaG), "nirKrecA"])/2
# mean_read_ratios["D198B", "nirKrecA"] = sum(read_ratios[grepl(glob2rx("D*198B"), read_ratios$samples_metaG), "nirKrecA"])/2
# mean_read_ratios["D234A", "nirKrecA"] = sum(read_ratios[grepl(glob2rx("D*234A"), read_ratios$samples_metaG), "nirKrecA"])/2
# 
# # Rest of matrix should be zero
# mean_read_ratios[is.na(mean_read_ratios)] <- 0

# Import metadata
meta <- read_delim("Metadata.txt",delim = "\t", col_names = T)

# Link read_ratios to metadata
read_ratios_metaG <- left_join(read_ratios_metaG, meta, by = "samples_metaG")
read_ratios_metaT <- left_join(read_ratios_metaT, meta, by = "samples_metaT")



# Plot by depth


#nirS:RecA
read_ratios_metaG$SizeFraction_f = factor(read_ratios_metaG$SizeFraction, levels = c('PA', 'FL'))
read_ratios_metaT$SizeFraction_f = factor(read_ratios_metaT$SizeFraction, levels = c('PA', 'FL'))

# nirSrecA <- ggplot(read_ratios, aes(y = nirSrecA, x = Depth, group = Depth))+
#   geom_boxplot(width = 10)+
#   scale_y_continuous(limits = c(0, 1.5))+
#   scale_x_reverse()+
#   ylab("nirS:recA")+
#   xlab("Depth (m)")+
#   coord_flip()+
#   facet_wrap(Season~SizeFraction_f)+
#   theme_bw()+
#   theme(axis.text.x = element_text(size=10),
#         axis.text.y = element_text(size=10),
#         axis.title.x= element_text(size=12),
#         axis.title.y= element_text(size=12),
#         strip.text = element_text(size = 14))
# ggsave(nirSrecA, file = "figures/nirSrecA.depthprofile.eps", units = "in", width = 4, height = 6) 


#mod

# Plot PA and FL on same plot                             
nirSrecA_DNA <- ggplot(read_ratios_metaG, aes(x = Depth, y = nirSrecA_DNA, group = Depth))+
  geom_boxplot(data = subset(read_ratios_metaG, SizeFraction == "PA"), width = 10, color = "black", show.legend = TRUE)+
  geom_boxplot(data = subset(read_ratios_metaG, SizeFraction == "FL"), width = 10, color = "grey", show.legend = TRUE)+
  scale_y_continuous(limits = c(0, 1.5))+
  scale_x_reverse()+
  ylab("nirS:recA")+
  xlab("Depth (m)")+
  coord_flip()+
  facet_wrap(~Season)+
  theme_bw()+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x= element_text(size=12),
        axis.title.y= element_text(size=12),
        strip.text = element_text(size = 14))
nirSrecA_DNA
ggsave(nirSrecA_DNA, file = "figures/nirSrecA_DNA.depthprofile.eps", units = "mm", width = 80, height = 80) 


nirSrecA_RNA <- ggplot(read_ratios_metaT, aes(x = Depth, y = nirSrecA_RNA, group = Depth))+
  geom_boxplot(data = subset(read_ratios_metaT, SizeFraction == "PA"), width = 10, color = "black", show.legend = TRUE)+
  geom_boxplot(data = subset(read_ratios_metaT, SizeFraction == "FL"), width = 10, color = "grey", show.legend = TRUE)+
  scale_y_continuous(limits = c(0, 10))+
  scale_x_reverse()+
  ylab("nirS:recA")+
  xlab("Depth (m)")+
  coord_flip()+
  facet_wrap(~Season)+
  theme_bw()+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x= element_text(size=12),
        axis.title.y= element_text(size=12),
        strip.text = element_text(size = 14))
nirSrecA_RNA
ggsave(nirSrecA, file = "figures/nirSrecA_RNA.depthprofile.eps", units = "mm", width = 80, height = 80) 
# NOTE that 2 rows of data did not plot here. The ratio of nirS:recA was infinite because there were no recA
# contigs. But nirS was detected. Those samples were R2b314A (AnoxMay2PA) and R2b237B (SubOxNov2FL)







#nirK:RecA
# nirKrecA <- ggplot(read_ratios, aes(y = nirKrecA, x = Depth, group = Depth))+
#   geom_boxplot(width = 10)+
#   scale_y_continuous(limits = c(0, 0.3))+ # Notice different scale
#   scale_x_reverse()+
#   ylab("nirK:recA")+
#   xlab("Depth (m)")+
#   coord_flip()+
#   facet_wrap(Season~SizeFraction_f)+
#   theme_bw()+
#   theme(axis.text.x = element_text(size=10),
#         axis.text.y = element_text(size=10),
#         axis.title.x= element_text(size=12),
#         axis.title.y= element_text(size=12),
#         strip.text = element_text(size = 14))
# nirKrecA
# ggsave(nirKrecA, file = "figures/nirKrecA.depthprofile.eps", units = "in", width = 4, height = 6) 

# Plot PA and FL on same plot                             
nirKrecA_DNA <- ggplot(read_ratios_metaG, aes(x = Depth, y = nirKrecA_DNA, group = Depth))+
  geom_boxplot(data = subset(read_ratios_metaG, SizeFraction == "PA"), width = 10, color = "black", show.legend = TRUE)+
  geom_boxplot(data = subset(read_ratios_metaG, SizeFraction == "FL"), width = 10, color = "grey", show.legend = TRUE)+
  scale_y_continuous(limits = c(0, 0.3))+
  scale_x_reverse()+
  ylab("nirK:recA")+
  xlab("Depth (m)")+
  coord_flip()+
  facet_wrap(~Season)+
  theme_bw()+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x= element_text(size=12),
        axis.title.y= element_text(size=12),
        strip.text = element_text(size = 14))
nirKrecA_DNA
ggsave(nirKrecA_DNA, file = "figures/nirKrecA_DNA.depthprofile.eps", units = "mm", width = 80, height = 80) 


nirKrecA_RNA <- ggplot(read_ratios_metaT, aes(x = Depth, y = nirKrecA_RNA, group = Depth))+
  geom_boxplot(data = subset(read_ratios_metaT, SizeFraction == "PA"), width = 10, color = "black", show.legend = TRUE)+
  geom_boxplot(data = subset(read_ratios_metaT, SizeFraction == "FL"), width = 10, color = "grey", show.legend = TRUE)+
  scale_y_continuous(limits = c(0, 0.3))+
  scale_x_reverse()+
  ylab("nirK:recA")+
  xlab("Depth (m)")+
  coord_flip()+
  facet_wrap(~Season)+
  theme_bw()+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x= element_text(size=12),
        axis.title.y= element_text(size=12),
        strip.text = element_text(size = 14))
nirKrecA_RNA
ggsave(nirKrecA_RNA, file = "figures/nirKrecA_RNA.depthprofile.eps", units = "mm", width = 80, height = 80) 
# Interesting that most of the nirK genes were shallower but the only depth from which nirK could be detected
# in the metaT samples was from an anoxic sample






#===========================================================#
# Next steps: pull out taxon IDs for nirS and nirK
#===========================================================#


#==================#
#=== nirS first ===#

# nirS metaG samples #

# Extract list of taxa for each samples_metaG. Apend to one long list with sample ID, taxon ID 
# (Accession), and relative abundance of each taxon (for that sample)

# Run first one to make data frame. Then run rest in loop and append to the first one
rm(i)
i = 1
taxonabun_nirS_metaG <- read_delim(paste("Xander_results/",nirSsamples_metaG[i],"_nirS_45_taxonabund.txt", sep = ""),delim = "\t", col_names = FALSE)
j <- which(grepl("Lineage" ,taxonabun_nirS_metaG$X1)) 
taxonabun_nirS_metaG <- taxonabun_nirS_metaG[j:nrow(taxonabun_nirS_metaG), ]
colnames(taxonabun_nirS_metaG) <- taxonabun_nirS_metaG[1, ]
taxonabun_nirS_metaG <- taxonabun_nirS_metaG[-c(1), ]
taxonabun_nirS_metaG$samples_metaG <- nirSsamples_metaG[i]
taxonabun_nirS_metaG$RelAbun <- NA
for(k in 1:nrow(taxonabun_nirS_metaG)){
  taxonabun_nirS_metaG$RelAbun[k] <- as.numeric(taxonabun_nirS_metaG$Abundance[k])/sum(as.numeric(taxonabun_nirS_metaG$Abundance))
}

for(i in 2:length(nirSsamples_metaG)){
  temp <- read_delim(paste("Xander_results/",nirSsamples_metaG[i],"_nirS_45_taxonabund.txt", sep = ""),delim = "\t", col_names = FALSE)
  j <- which(grepl("Lineage" ,temp$X1)) 
  temp <- temp[j:nrow(temp), ]
  colnames(temp) <- temp[1, ]
  temp <- temp[-c(1), ]
  temp$samples_metaG <- nirSsamples_metaG[i]
  temp$RelAbun <- NA
  for(k in 1:nrow(temp)){
    temp$RelAbun[k] <- as.numeric(temp$Abundance[k])/sum(as.numeric(temp$Abundance))
  }
  taxonabun_nirS_metaG <- rbind(taxonabun_nirS_metaG, temp)
  
}

# Attach metadata
# Link read_ratios to metadata
taxonabun_nirS_metaG <- left_join(taxonabun_nirS_metaG, meta, by = "samples_metaG")

# First glance. Look at rel abun of difference groups from all samples
# ggplot(taxonabun_nirS_metaG,aes (x = Lineage, y = RelAbun))+
#   geom_boxplot(aes(color=SizeFraction))+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))


# there are patterns here but can't see the lineage names. Start messing with the text. 
# Look up taxonomy of closest hit and put that as variable
# Used uniprot. Put genus whenever possible
taxonomy <- tibble(unique(taxonabun_nirS_metaG$Lineage))
colnames(taxonomy) <- ("Lineage")


# Looked up phylogeny of every hit. For ones that were not identified (eg. marine sediment metagenome,
# definition=hypothetical protein, MatchName = KKM16877), I ran the contig through blastp. For this one in particular, 
# all top hits were Scalindua. Indicated as such below (#6))
taxonomy$Taxonomy <- c("Bacteria; Proteobacteria; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; Labrenzia",
                       "Bacteria; Proteobacteria; Gammaproteobacteria; Alteromonadales; Alteromonadaceae; Marinobacter",
                       "Bacteria; Proteobacteria; Gammaproteobacteria; Competibacteraceae; Candidatus Competibacter",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; Labrenzia",
                       "Bacteria; Proteobacteria; Gammaproteobacteria; Oceanospirillales; Oleiphilaceae; Oleiphilus",
                       "Bacteria; Planctomycetes; Planctomycetia; Candidatus Brocadiales; Candidatus Brocadiaceae; Candidatus Scalindua",
                       "Bacteria; Proteobacteria; Acidithiobacillia; Acidithiobacillales",
                       "Bacteria; Proteobacteria; Gammaproteobacteria; Sedimenticola",
                       "Bacteria; Proteobacteria; Gammaproteobacteria; Acidiferrobacterales; Acidiferrobacteraceae; Sulfurifustis",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Polymorphum",
                       "Bacteria; Proteobacteria; Betaproteobacteria; Nitrosomonadales; Methylophilaceae; Methylobacillus",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Bradyrhizobiaceae; Bradyrhizobium",
                       "Bacteria; Proteobacteria; Betaproteobacteria; Rhodocyclales; Rhodocyclaceae; Azospira",
                       "Bacteria; Proteobacteria; Betaproteobacteria; Burkholderiales; Burkholderiaceae; Cupriavidus",
                       "Bacteria; Proteobacteria; Betaproteobacteria; Nitrosomonadales; Thiobacillaceae; Thiobacillus",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; Labrenzia",
                       "Bacteria; Chloroflexi; Anaerolineae; Anaerolineales; Anaerolineaceae; Anaerolinea",
                       "Bacteria; Chloroflexi; Anaerolineae; Anaerolineales; Anaerolineaceae; Bellilinea",
                       "Bacteria; Proteobacteria; Gammaproteobacteria; Oceanospirillales; Oceanospirillaceae; Oleispira",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; Labrenzia",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; Pseudovibrio",
                       "Bacteria; Proteobacteria; Gammaproteobacteria; Thiotrichales; Thiotrichaceae; Thioploca",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhodospirillales; Rhodospirillaceae; Magnetospirillum",
                       "Bacteria; Proteobacteria; Betaproteobacteria; Nitrosomonadales; Sterolibacteriaceae; Sulfuritalea",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; Pseudovibrio",
                       "Bacteria; Proteobacteria; Hydrogenophilalia; Hydrogenophilales; Hydrogenophilaceae; Tepidiphilus",
                       "Bacteria; Proteobacteria; Gammaproteobacteria; sulfur-oxidizing symbionts",
                       "Bacteria; Proteobacteria; Gammaproteobacteria; Xanthomonadales; Xanthomonadaceae; Arenimonas",
                       "Bacteria; Proteobacteria; Betaproteobacteria; Rhodocyclales; Azonexaceae; Dechloromonas",
                       "Bacteria; Chloroflexi; Anaerolineae; Anaerolineales; Anaerolineaceae; Thermanaerothrix",
                       "Bacteria; Proteobacteria; Betaproteobacteria; unclassified Betaproteobacteria",
                       "Bacteria; Proteobacteria; Gammaproteobacteria; Sedimenticola")

# Split taxonomy into columns
nirS_otu.table.DNA <- separate(taxonomy,2, c('Taxonomy1','Taxonomy2','Taxonomy3','Taxonomy4','Taxonomy5','Taxonomy6', 'Taxonomy7'), sep = ';')


# Join nirS_otu.table.DNA and taxonabun_nirS_metaG table
taxonabun_nirS_metaG <- left_join(taxonabun_nirS_metaG, nirS_otu.table.DNA, by = "Lineage")


# Plot taxononmy of all samples. Look at rel abundances of different taxa. Group by taxonomic level
taxonabun_nirS_metaG$SizeFraction_f = factor(taxonabun_nirS_metaG$SizeFraction, levels = c('PA', 'FL'))

nirSrecA_orders_DNA <- ggplot(taxonabun_nirS_metaG,aes (x = Taxonomy4, y = RelAbun, color=SizeFraction_f))+
  geom_boxplot(size = 1, width = 0.5)+
  scale_color_grey()+
  ylab("nirS:recA")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size=10),
        axis.title.x= element_blank(),
        axis.title.y= element_text(size=12),
        legend.text=element_text(size=12),
        legend.position="none")
nirSrecA_orders_DNA
ggsave(nirSrecA_orders_DNA, file = "figures/nirSrecA_orders_DNA.orders.eps", units = "mm", width = 168, height = 80) 



# Group by oxygen condition
# First order OxCond in correct order by depth, otherwise it plots alphabetically
taxonabun_nirS_metaG$OxCond_f = factor(taxonabun_nirS_metaG$OxCond, levels = c('Oxycline', 'Suboxic', 'ShallowAnoxic'))

nirSrecA_oxcond_DNA <- ggplot(taxonabun_nirS_metaG,aes (x = OxCond_f, y = RelAbun, color=SizeFraction_f))+
  geom_boxplot(size = 1)+
  scale_color_grey()+
  ylab("nirS:recA")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size=10),
        axis.title.x= element_blank(),
        axis.title.y= element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_blank())
nirSrecA_oxcond_DNA
ggsave(nirSrecA_oxcond_DNA, file = "figures/nirSrecA_oxcond_DNA.oxcond.eps", units = "in", width = 6, height = 3) 





#  nirS metaT samples #
rm(i)
i = 1
taxonabun_nirS_metaT <- read_delim(paste("Xander_results/",nirSsamples_metaT[i],"_nirS_45_taxonabund.txt", sep = ""),delim = "\t", col_names = FALSE)
j <- which(grepl("Lineage" ,taxonabun_nirS_metaT$X1)) 
taxonabun_nirS_metaT <- taxonabun_nirS_metaT[j:nrow(taxonabun_nirS_metaT), ]
colnames(taxonabun_nirS_metaT) <- taxonabun_nirS_metaT[1, ]
taxonabun_nirS_metaT <- taxonabun_nirS_metaT[-c(1), ]
taxonabun_nirS_metaT$samples_metaT <- nirSsamples_metaT[i]
taxonabun_nirS_metaT$RelAbun <- NA
for(k in 1:nrow(taxonabun_nirS_metaT)){
  taxonabun_nirS_metaT$RelAbun[k] <- as.numeric(taxonabun_nirS_metaT$Abundance[k])/sum(as.numeric(taxonabun_nirS_metaT$Abundance))
}

for(i in 2:length(nirSsamples_metaT)){
  temp <- read_delim(paste("Xander_results/",nirSsamples_metaT[i],"_nirS_45_taxonabund.txt", sep = ""),delim = "\t", col_names = FALSE)
  j <- which(grepl("Lineage" ,temp$X1)) 
  temp <- temp[j:nrow(temp), ]
  colnames(temp) <- temp[1, ]
  temp <- temp[-c(1), ]
  temp$samples_metaT <- nirSsamples_metaT[i]
  temp$RelAbun <- NA
  for(k in 1:nrow(temp)){
    temp$RelAbun[k] <- as.numeric(temp$Abundance[k])/sum(as.numeric(temp$Abundance))
  }
  taxonabun_nirS_metaT <- rbind(taxonabun_nirS_metaT, temp)
  
}


# Attach metadata
taxonabun_nirS_metaT <- left_join(taxonabun_nirS_metaT, meta, by = "samples_metaT")


# Look up taxonomy of closest hit and put that as variable
# Used uniprot. Put genus whenever possible
rm(taxonomy)
taxonomy <- tibble(unique(taxonabun_nirS_metaT$Lineage))
colnames(taxonomy) <- ("Lineage")

# Manually put phylogeny of results from "lineage" list
# Use uniprot for phylogeny
taxonomy$Taxonomy <- c("Bacteria; Planctomycetes; Planctomycetia; Candidatus Brocadiales; Candidatus Brocadiaceae; Candidatus Scalindua",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; Labrenzia",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; Ruegeria",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhodospirillales; Rhodospirillaceae; Thalassospira",
                       "Bacteria; Proteobacteria; Gammaproteobacteria; Oceanospirillales; Oleiphilaceae; Oleiphilus",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Polymorphum",
                       "Bacteria; Proteobacteria; Gammaproteobacteria; sulfur-oxidizing symbionts")

# Split taxonomy into columns
nirS_otu.table.RNA <- separate(taxonomy,2, c('Taxonomy1','Taxonomy2','Taxonomy3','Taxonomy4','Taxonomy5','Taxonomy6', 'Taxonomy7'), sep = ';')

# Join nirS_otu.table.RNA and taxonabun_nirS_metaT table
taxonabun_nirS_metaT <- left_join(taxonabun_nirS_metaT, nirS_otu.table.RNA, by = "Lineage")

# Plot taxononmy of all samples. Look at rel abundances of different taxa. Group by taxonomic level
taxonabun_nirS_metaT$SizeFraction_f = factor(taxonabun_nirS_metaT$SizeFraction, levels = c('PA', 'FL'))

nirSrecA_orders_RNA <- ggplot(taxonabun_nirS_metaT,aes (x = Taxonomy4, y = RelAbun, color=SizeFraction_f))+
  geom_boxplot(size = 1)+
  scale_color_grey()+
  ylab("nirS:recA")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size=10),
        axis.title.x= element_blank(),
        axis.title.y= element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_blank())
nirSrecA_orders_RNA
ggsave(nirSrecA_orders_RNA, file = "figures/nirSrecA_orders_RNA.orders.eps", units = "in", width = 6, height = 3) 



# Group by oxygen condition
taxonabun_nirS_metaT$OxCond_f = factor(taxonabun_nirS_metaT$OxCond, levels = c('Oxycline', 'Suboxic', 'ShallowAnoxic'))

nirSrecA_oxcond_RNA <- ggplot(taxonabun_nirS_metaT,aes (x = OxCond_f, y = RelAbun, color=SizeFraction_f))+
  geom_boxplot(size = 1)+
  scale_color_grey()+
  ylab("nirS:recA")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size=10),
        axis.title.x= element_blank(),
        axis.title.y= element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_blank())
nirSrecA_oxcond_RNA
ggsave(nirSrecA_oxcond_RNA, file = "figures/nirSrecA_oxcond_RNA.oxcond.eps", units = "in", width = 6, height = 3) 









#==================#
#==== now nirK ====#

# nirK metaG first #

# samples that had nirK:
samples_metaG <- c("D2a143B","D2a237B","D2b143B","D2b200B","D2b237B","D3a198A",
                   "D3a198B","D3a234A","D3a234B", "D3b198B","D3b234A","D3b234B")
# Extract list of taxa for each sample. Apend to one long list with sample ID, taxon ID 
# (Accession), and relative abundance of each taxon (for that sample)

# Run first one to make data frame. Then run rest in loop and append to the first one
rm(i)
rm(temp)
i = 1
taxonabun_nirK_metaG <- read_delim(paste("Xander_results/",samples_metaG[i],"_nirK_45_taxonabund.txt", sep = ""),delim = "\t", col_names = FALSE)
j <- which(grepl("Lineage" ,taxonabun_nirK_metaG$X1)) 
taxonabun_nirK_metaG <- taxonabun_nirK_metaG[j:nrow(taxonabun_nirK_metaG), ]
colnames(taxonabun_nirK_metaG) <- taxonabun_nirK_metaG[1, ]
taxonabun_nirK_metaG <- taxonabun_nirK_metaG[-c(1), ]
taxonabun_nirK_metaG$samples_metaG <- samples_metaG[i]
taxonabun_nirK_metaG$RelAbun <- NA
for(k in 1:nrow(taxonabun_nirK_metaG)){
  taxonabun_nirK_metaG$RelAbun[k] <- as.numeric(taxonabun_nirK_metaG$Abundance[k])/sum(as.numeric(taxonabun_nirK_metaG$Abundance))
}

for(i in 2:length(samples_metaG)){
  temp <- read_delim(paste("Xander_results/",samples_metaG[i],"_nirK_45_taxonabund.txt", sep = ""),delim = "\t", col_names = FALSE)
  j <- which(grepl("Lineage" ,temp$X1)) 
  temp <- temp[j:nrow(temp), ]
  colnames(temp) <- temp[1, ]
  temp <- temp[-c(1), ]
  temp$samples_metaG <- samples_metaG[i]
  temp$RelAbun <- NA
  for(k in 1:nrow(temp)){
    temp$RelAbun[k] <- as.numeric(temp$Abundance[k])/sum(as.numeric(temp$Abundance))
  }
  taxonabun_nirK_metaG <- rbind(taxonabun_nirK_metaG, temp)
  
}

# Attach metadata
# Link read_ratios to metadata
taxonabun_nirK_metaG <- left_join(taxonabun_nirK_metaG, meta, by = "samples_metaG")

# First glance. Look at rel abun of difference groups from all samples
# ggplot(taxonabun_nirK_metaG,aes (x = Lineage, y = RelAbun))+
#   geom_boxplot(aes(color=SizeFraction))+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))


# there are patterns here but can't see the lineage names. Start messing with the text. 
# Look up taxonomy of closest hit from xander (fungene) and put that as variable
# Used uniprot for phylogeny. Put genus whenever possible
taxonomy <- tibble(unique(taxonabun_nirK_metaG$Lineage))
colnames(taxonomy) <- ("Lineage")


# Looked up phylogeny of every hit. For ones that were not identified (eg. marine sediment metagenome,
# definition=hypothetical protein, MatchName = KKM16877), I ran the contig through blastp. For this one in particular, 
# all top hits were Scalindua. Indicated as such below (#6))
taxonomy$Taxonomy <- c("Bacteria; Proteobacteria; Gammaproteobacteria; Chromatiales; Chromatiaceae; Nitrosococcus",
                       "Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Pseudomonadaceae; Pseudomonas",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Brucellaceae; Brucella",
                       "Bacteria; Proteobacteria; Gammaproteobacteria; Alteromonadales; Shewanellaceae; Shewanella",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Brucellaceae; Ochrobactrum",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Bradyrhizobiaceae; Rhodopseudomonas",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Brucellaceae; Brucella",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Brucellaceae; Brucella",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Rhizobiaceae; Sinorhizobium/Ensifer group; Sinorhizobium",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Bradyrhizobiaceae; Afipia",
                       "Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Brucellaceae; Brucella")



                                                                                                                                                                                                                              
# Split taxonomy into columns
nirK_otu.table.DNA <- separate(taxonomy,2, c('Taxonomy1','Taxonomy2','Taxonomy3','Taxonomy4','Taxonomy5','Taxonomy6', 'Taxonomy7'), sep = ';')


# Join nirK_otu.table.DNA and taxonabun_nirK_metaG table
taxonabun_nirK_metaG <- left_join(taxonabun_nirK_metaG, nirK_otu.table.DNA, by = "Lineage")


# Plot taxononmy of all samples. Look at rel abundances of different taxa. Group by taxonomic level
taxonabun_nirK_metaG$SizeFraction_f = factor(taxonabun_nirK_metaG$SizeFraction, levels = c('PA', 'FL'))

nirKrecA_orders_DNA <- ggplot(taxonabun_nirK_metaG,aes (x = Taxonomy4, y = RelAbun, color=SizeFraction_f))+
  geom_boxplot(size = 1, width = 0.5)+
  scale_color_grey()+
  ylab("nirK:recA")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size=10),
        axis.title.x= element_blank(),
        axis.title.y= element_text(size=12),
        legend.text=element_text(size=12),
        legend.position="none")
nirKrecA_orders_DNA
ggsave(nirKrecA_orders_DNA, file = "figures/nirKrecA_orders_DNA.orders.eps", units = "mm", width = 80, height = 60) 



# Group by oxygen condition (only present in oxycline)
taxonabun_nirK_metaG$OxCond_f = factor(taxonabun_nirK_metaG$OxCond, levels = c('Oxycline', 'Suboxic', 'ShallowAnoxic'))

nirKrecA_oxcond_DNA <- ggplot(taxonabun_nirK_metaG,aes (x = OxCond, y = RelAbun, color=SizeFraction_f))+
  geom_boxplot(size = 1)+
  scale_color_grey()+
  ylab("nirK:recA")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size=10),
        axis.title.x= element_blank(),
        axis.title.y= element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_blank())
nirKrecA_oxcond_DNA
ggsave(nirKrecA_oxcond_DNA, file = "figures/nirKrecA_oxcond_DNA.oxcond.eps", units = "in", width = 6, height = 3) 




#  nirK metaT samples #
rm(i)
i = 1
taxonabun_nirK_metaT <- read_delim(paste("Xander_results/",nirKsamples_metaT[i],"_nirS_45_taxonabund.txt", sep = ""),delim = "\t", col_names = FALSE)
j <- which(grepl("Lineage" ,taxonabun_nirK_metaT$X1)) 
taxonabun_nirK_metaT <- taxonabun_nirK_metaT[j:nrow(taxonabun_nirK_metaT), ]
colnames(taxonabun_nirK_metaT) <- taxonabun_nirK_metaT[1, ]
taxonabun_nirK_metaT <- taxonabun_nirK_metaT[-c(1), ]
taxonabun_nirK_metaT$samples_metaT <- nirKsamples_metaT[i]
taxonabun_nirK_metaT$RelAbun <- NA
for(k in 1:nrow(taxonabun_nirK_metaT)){
  taxonabun_nirK_metaT$RelAbun[k] <- as.numeric(taxonabun_nirK_metaT$Abundance[k])/sum(as.numeric(taxonabun_nirK_metaT$Abundance))
}

for(i in 2:length(nirKsamples_metaT)){
  temp <- read_delim(paste("Xander_results/",nirKsamples_metaT[i],"_nirS_45_taxonabund.txt", sep = ""),delim = "\t", col_names = FALSE)
  j <- which(grepl("Lineage" ,temp$X1)) 
  temp <- temp[j:nrow(temp), ]
  colnames(temp) <- temp[1, ]
  temp <- temp[-c(1), ]
  temp$samples_metaT <- nirKsamples_metaT[i]
  temp$RelAbun <- NA
  for(k in 1:nrow(temp)){
    temp$RelAbun[k] <- as.numeric(temp$Abundance[k])/sum(as.numeric(temp$Abundance))
  }
  taxonabun_nirK_metaT <- rbind(taxonabun_nirS_metaT, temp)
  
}


# Attach metadata
taxonabun_nirK_metaT <- left_join(taxonabun_nirK_metaT, meta, by = "samples_metaT")


# Look up taxonomy of closest hit and put that as variable
# Used uniprot. Put genus whenever possible
rm(taxonomy)
taxonomy <- tibble(unique(taxonabun_nirK_metaT$Lineage))
colnames(taxonomy) <- ("Lineage")

# Manually put phylogeny of results from "lineage" list
# Use uniprot for phylogeny
taxonomy$Taxonomy <- c("Bacteria; Proteobacteria; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; Labrenzia",
                       "Bacteria; Proteobacteria; Gammaproteobacteria; Oceanospirillales; Oleiphilaceae; Oleiphilus")
                       
# Split taxonomy into columns
nirK_otu.table.RNA <- separate(taxonomy,2, c('Taxonomy1','Taxonomy2','Taxonomy3','Taxonomy4','Taxonomy5','Taxonomy6', 'Taxonomy7'), sep = ';')

# Join nirS_otu.table.RNA and taxonabun_nirS_metaT table
taxonabun_nirK_metaT <- left_join(taxonabun_nirK_metaT, nirK_otu.table.RNA, by = "Lineage")

# Plot taxononmy of all samples. Look at rel abundances of different taxa. Group by taxonomic level
taxonabun_nirK_metaT$SizeFraction_f = factor(taxonabun_nirK_metaT$SizeFraction, levels = c('PA', 'FL'))

nirKrecA_orders_RNA <- ggplot(taxonabun_nirK_metaT,aes (x = Taxonomy4, y = RelAbun, color=SizeFraction_f))+
  geom_boxplot(size = 1, color = "grey")+
  scale_color_grey()+
  ylab("nirK:recA")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size=10),
        axis.title.x= element_blank(),
        axis.title.y= element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_blank())
nirKrecA_orders_RNA
ggsave(nirKrecA_orders_RNA, file = "figures/nirKrecA_orders_RNA.orders.eps", units = "in", width = 6, height = 3) 



# Group by oxygen condition (only present in ShallowAnoxic)
taxonabun_nirK_metaT$OxCond_f = factor(taxonabun_nirK_metaT$OxCond, levels = c('Oxycline', 'Suboxic', 'ShallowAnoxic'))

nirKrecA_oxcond_RNA <- ggplot(taxonabun_nirK_metaT,aes (x = OxCond, y = RelAbun, color=SizeFraction_f))+
  geom_boxplot(size = 1, color = "grey")+
  scale_color_grey()+
  ylab("nirK:recA")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size=10),
        axis.title.x= element_blank(),
        axis.title.y= element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_blank())
nirKrecA_oxcond_RNA
ggsave(nirKrecA_oxcond_RNA, file = "figures/nirKrecA_oxcond_RNA.oxcond.eps", units = "in", width = 6, height = 3) 








#=================================#
# Combine nirS and nirK OxCond plots in one plot
# To do this, must put dataframes in same dataframe. First put a column with gene name to keep track


#DNA first
taxonabun_nirS_metaG$Gene <- "nirS"
taxonabun_nirK_metaG$Gene <- "nirK"
taxonabun_metaG <- rbind(taxonabun_nirS_metaG, taxonabun_nirK_metaG)
# And re order so nirS gets plotted before nirK
taxonabun_metaG$Gene_f = factor(taxonabun_metaG$Gene, levels = c('nirS', 'nirK'))


nirS_K_oxcond_DNA <- ggplot(taxonabun_metaG, aes(x = OxCond_f, y = RelAbun, color=SizeFraction_f))+
  geom_boxplot(size = 1)+
  scale_color_grey()+
  ylab("ratio to recA")+
  facet_wrap(~Gene_f)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size=10),
        axis.title.x= element_blank(),
        axis.title.y= element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        legend.position="none")
nirS_K_oxcond_DNA
ggsave(nirS_K_oxcond_DNA, file = "figures/nirS_K_oxcond_DNA.eps", units = "mm", width = 168, height = 80) 

# RNA
taxonabun_nirS_metaT$Gene <- "nirS"
taxonabun_nirK_metaT$Gene <- "nirK"
taxonabun_metaT <- rbind(taxonabun_nirS_metaT, taxonabun_nirK_metaT)
# And re order so nirS gets plotted before nirK
taxonabun_metaT$Gene_f = factor(taxonabun_metaT$Gene, levels = c('nirS', 'nirK'))


nirS_K_oxcond_RNA <- ggplot(taxonabun_metaT, aes(x = OxCond_f, y = RelAbun, color=SizeFraction_f))+
  geom_boxplot(size = 1)+
  scale_color_grey()+
  ylab("ratio to recA")+
  facet_wrap(~Gene_f)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size=10),
        axis.title.x= element_blank(),
        axis.title.y= element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_blank())
nirS_K_oxcond_RNA
ggsave(nirS_K_oxcond_RNA, file = "figures/nirS_K_oxcond_RNA.eps", units = "mm", width = 168, height = 80) 





# ================================ #
# Plot taxonomy with depth profiles
# ================================ #


## nirS- metaG samples

# First re order size fraction so it plots in correct order:
taxonabun_nirS_metaG$SizeFraction_f = factor(taxonabun_nirS_metaG$SizeFraction, levels = c('PA', 'FL'))

# ggplot()+
#   geom_bar(aes(x = Depth, y = RelAbun, fill = Taxonomy3), data = taxonabun_nirS_metaG, width = 10,
#            stat="identity")+
#   coord_flip()+
#   ylab("Relative Abundance")+
#   xlab("Depth (m)")+
#   scale_x_reverse()+
#   facet_wrap(Season~SizeFraction_f)+
#   theme_bw()+
#   theme(axis.text.x = element_text(size=12),
#         axis.text.y = element_text(size=12),
#         axis.title.x= element_text(size=16),
#         axis.title.y= element_text(size=16))
# This looks fine but the replicates are plotted in same bar, so total relative abundance is 2.0, not 1.0. 

# For one depth in Nov PA and FL (200m), there is only one sample that had nirS so this plots to 1.0 already
# To solve this, divide the relative abundances (except Nov 200m) by 2
taxonabun_nirS_metaG$RelAbun_f = as.numeric(taxonabun_nirS_metaG$RelAbun)/2
rm (i, j, k)
i <- which(grepl(201, taxonabun_nirS_metaG$Depth))
taxonabun_nirS_metaG$RelAbun_f[i] <- taxonabun_nirS_metaG$RelAbun[i]

# Replot
nirS_taxonomy_DNA <- ggplot()+
  geom_bar(aes(x = Depth, y = RelAbun_f, fill = Taxonomy3), data = taxonabun_nirS_metaG, width = 10,
           stat="identity")+
  coord_flip()+
  ylab("Relative Abundance, nirS")+
  xlab("Depth (m)")+
  scale_x_reverse(limits = c(325, 100))+
  facet_wrap(Season~SizeFraction_f)+
  scale_fill_grey(start = 0, end = .9)+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        strip.text = element_text(size = 12),
        legend.text=element_text(size=8),
        legend.title=element_blank())
nirS_taxonomy_DNA
ggsave(nirS_taxonomy_DNA, file = "figures/nirS_taxonomy_DNA.eps", units = "in", width = 6, height = 6) 




## nirS- metaT samples

# Firts re order size fraction so it plots in correct order:
taxonabun_nirS_metaT$SizeFraction_f = factor(taxonabun_nirS_metaT$SizeFraction, levels = c('PA', 'FL'))

# Divide rel abun at each depth by two to reflect average among 2 samples
taxonabun_nirS_metaT$RelAbun_f = as.numeric(taxonabun_nirS_metaT$RelAbun)/2
# For some depths, there is only one sample that had nirS so this plots to 1.0 already (dividing in half means the total rel abun is 0.5)
# To solve this, re-replace Rel_abun_f with Rel_abunrm (i, j, k)
i <- which(grepl("R2b314A", taxonabun_nirS_metaT$samples_metaT))
taxonabun_nirS_metaT$RelAbun_f[i] <- taxonabun_nirS_metaT$RelAbun[i]
j <- which(grepl("R2a200B", taxonabun_nirS_metaT$samples_metaT))
taxonabun_nirS_metaT$RelAbun_f[j] <- taxonabun_nirS_metaT$RelAbun[j]
k <- which(grepl("R2a247B", taxonabun_nirS_metaT$samples_metaT))
taxonabun_nirS_metaT$RelAbun_f[k] <- taxonabun_nirS_metaT$RelAbun[k]

# Replot
nirS_taxonomy_RNA <- ggplot()+
  geom_bar(aes(x = Depth, y = RelAbun_f, fill = Taxonomy3), data = taxonabun_nirS_metaT, width = 10,
           stat="identity")+
  coord_flip()+
  ylab("Relative Abundance, nirS")+
  xlab("Depth (m)")+
  scale_x_reverse(limits = c(325, 100))+
  facet_wrap(Season~SizeFraction_f)+
  scale_fill_grey(start = 0, end = .9)+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        strip.text = element_text(size = 12),
        legend.text=element_text(size=8),
        legend.title=element_blank())
nirS_taxonomy_RNA
ggsave(nirS_taxonomy_RNA, file = "figures/nirS_taxonomy_RNA.eps", units = "in", width = 6, height = 6) 









## nirK- metaG samples

# Firts re order size fraction so it plots in correct order:
taxonabun_nirK_metaG$SizeFraction_f = factor(taxonabun_nirK_metaG$SizeFraction, levels = c('PA', 'FL'))

# ggplot()+
#   geom_bar(aes(x = Depth, y = RelAbun, fill = Taxonomy3), data = taxonabun_nirK_metaG, width = 10,
#            stat="identity")+
#   coord_flip()+
#   ylab("Relative Abundance")+
#   xlab("Depth (m)")+
#   scale_x_reverse()+
#   facet_wrap(Season~SizeFraction_f)+
#   theme_bw()+
#   theme(axis.text.x = element_text(size=12),
#         axis.text.y = element_text(size=12),
#         axis.title.x= element_text(size=16),
#         axis.title.y= element_text(size=16))
# This looks fine but the replicates are plotted in same bar, so total relative abundance is 2.0, not 1.0. 

# For one depth in Nov PA and FL (200m), there is only one sample that had nirS so this plots to 1.0 already
# To solve this, divide the relative abundances (except Nov 200m) by 2
taxonabun_nirK_metaG$RelAbun_f = as.numeric(taxonabun_nirK_metaG$RelAbun)/2
rm (i, j, k)
i <- which(grepl("D3a198A", taxonabun_nirK_metaG$samples_metaG))
taxonabun_nirK_metaG$RelAbun_f[i] <- taxonabun_nirK_metaG$RelAbun[i]
j <- which(grepl("D2b200B", taxonabun_nirK_metaG$samples_metaG))
taxonabun_nirK_metaG$RelAbun_f[j] <- taxonabun_nirK_metaG$RelAbun[j]


# Replot- NOTE nirK is only in FL in Nov so NovPA panel doesn't plot
nirK_taxonomy_DNA <- ggplot()+
  geom_bar(aes(x = Depth, y = RelAbun_f, fill = Taxonomy3), data = taxonabun_nirK_metaG, width = 10,
           stat="identity")+
  coord_flip()+
  ylab("Relative Abundance, nirK")+
  xlab("Depth (m)")+
  scale_x_reverse(limits = c(325, 100))+
  facet_wrap(Season~SizeFraction_f)+
  scale_fill_grey(start = 0, end = .9)+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        strip.text = element_text(size = 12),
        legend.text=element_text(size=8),
        legend.title=element_blank())
nirK_taxonomy_DNA
ggsave(nirK_taxonomy_DNA, file = "figures/nirK_taxonomy_DNA.eps", units = "in", width = 6, height = 3) 



## nirK- metaT samples

# Firts re order size fraction so it plots in correct order:
taxonabun_nirK_metaT$SizeFraction_f = factor(taxonabun_nirK_metaT$SizeFraction, levels = c('PA', 'FL'))

# No need to divide Rel abun by 2 since there is only one sample with nirK. 
# Total Rel abun is 1.0

# Plot
nirK_taxonomy_RNA <- ggplot()+
  geom_bar(aes(x = Depth, y = RelAbun, fill = Taxonomy3), data = taxonabun_nirK_metaT, width = 10,
           stat="identity")+
  coord_flip()+
  ylab("Relative Abundance, nirK")+
  xlab("Depth (m)")+
  scale_x_reverse(limits = c(325, 100))+
  facet_wrap(Season~SizeFraction_f)+
  scale_fill_grey(start = 0, end = .9)+
  theme_bw()+
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x= element_text(size=10),
        axis.title.y= element_text(size=10),
        strip.text = element_text(size = 12),
        legend.text=element_text(size=8),
        legend.title=element_blank())
nirK_taxonomy_RNA
ggsave(nirK_taxonomy_RNA, file = "figures/nirK_taxonomy_RNA.eps", units = "in", width = 2, height = 3) 
