################################################################################
#################### Load libraries  ###########################################
################################################################################

library(dplyr)
library(tidyverse)
library(ggpubr)
library(kableExtra)
library(tibble)
library(QuantPsyc)
library(circlize)

library(aspace)
library(png)
library(stringr)
library(ComplexHeatmap)
library(gtsummary)

################################################################################
#################### Functions #################################################
################################################################################

UKB.regressionAnalysis <- function(ukb.fmri, cog.val, mri.val, cog.name, index, results.table){
  
  fit.lm <- lm(cog.val ~ mri.val + mri.val:age_stand + mri.val:sex + sex:age_stand + 
                 sex:I(age_stand^2) + age_stand + I(age_stand^2) + sex + edu_uni + 
                 edu_Alevel + edu_Olevel + edu_cse + edu_nvq + edu_other +  site + 
                 headMotion + SNR + T1_discrepancy, data = ukb.fmri)
  
  fit.null <- lm(cog.val ~ sex:age_stand + sex:I(age_stand^2) + age_stand + I(age_stand^2) + 
                   sex + edu_uni + edu_Alevel + edu_Olevel + edu_cse + edu_nvq + edu_other + 
                   site + headMotion + SNR + T1_discrepancy, data = ukb.fmri)
  
  an.res <- anova(fit.lm,fit.null,test = "F")
  beta.stand <- lm.beta(fit.lm)
  results.table[index,c(paste(cog.name,'.coeff',sep=""), paste(cog.name, '.Pvalue',sep=""), paste(cog.name, '.Rsquared',sep=""))] <- c(beta.stand[1], an.res[2,6], summary(fit.lm)$adj.r.squared-summary(fit.null)$adj.r.squared)
  
  return(results.table)
}

circular.plot <- function(reg.results, plot.data, connectivity.data, title.name){
  
  connectivity.data$unique <- (reg.results[,c(2)] < 0.05) & (rowSums(reg.results[,c(3:ncol(reg.results))] < 0.05) == 0) 
  connectivity.data$shared <- (reg.results[,c(2)] < 0.05) & (rowSums(reg.results[,c(3:ncol(reg.results))] < 0.05) > 0) 
  connectivity.data$sum <- abs(reg.results[,c(1)])
  
  connectivity.data = connectivity.data %>%
    mutate(cColor1 = ifelse(unique == 1 , "#998ec3",
                            ifelse(shared == 1, "#f1a340" ,"#00000000"))) # transparent black
  
  plot.connectivity.f <- connectivity.data[,c("Network1","Network2","sum","cColor1")]
  
  par(cex = 1)
  # discrete
  chordDiagram(plot.connectivity.f, annotationTrack = c("grid","name"), preAllocateTracks = list(track.height = 0.25),col = plot.connectivity.f$cColor1, order=plot.data$names, grid.col = setNames(plot.data$gColor, plot.data$names),              link.lwd = 1,    # Line width
               link.lty = 1,    # Line type
               link.border = plot.connectivity.f$cColor1)
  par(cex = 1)
  title(main = title.name, cex.main = 1.5)
  
  u=0
  for(si in get.all.sector.index()){
    xplot=get.cell.meta.data("xplot",si)
    u=u+1
    
    # Small workaround because coordinate 0 should be 360
    if(xplot[1] == 0) xplot[1] = 360
    
    x=.86*cos(as_radians((xplot[2]+xplot[1])/2))
    y=.86*sin(as_radians((xplot[2]+xplot[1])/2))
    
    plot.data$images[grep(paste("\\b",si,"\\b",sep = ""), plot.data$names)] %>% 
      readPNG() %>% 
      rasterImage(x-0.11, y-0.11, x+0.11, y+0.11)
  }
  
  
}

################################################################################
#################### Load Data #################################################
################################################################################

# Load all files
ukb.data <-  read.delim(".../path/to/file",header = TRUE, sep = ",")
ukb.fmri <- read.delim(".../path/to/file",header = TRUE, sep = ",")

# Load individual cognitive tests and factor scores
ukb.cog <- read.delim(".../path/to/file",header = TRUE, sep = ",")
ukb.cog <- ukb.cog[,c(1:10)]
names(ukb.cog) <- c('eid', 'fluid.intelligence','max.digits','pairs.matching',
                      'reaction.time','matrix.pattern','symbol.digit','tower',
                      'trail.making','pairs')

# Re-reference processing speed variables based on maximum
ukb.cog$trail.making <- (max(ukb.cog$trail.making)+1) - ukb.cog$trail.making
ukb.cog$pairs.matching <- (max(ukb.cog$pairs.matching)+1) - ukb.cog$pairs.matching
ukb.cog$reaction.time <- (max(ukb.cog$reaction.time)+1) - ukb.cog$reaction.time

# Load g
g.factor <-read.delim(".../path/to/file", header = TRUE, sep = "\t", quote = "")
names(g.factor) <- c("eid","g")

# Load demographic variables
ukb_demog <-read.delim(".../path/to/file", header = TRUE, sep = "\t", quote = "")

# Load ethnicity data
ukb.ethnic <- read.table(".../path/to/file", header = TRUE)
names(ukb.ethnic) <- c("eid","ethnicity")
ukb.ethnic$ethnicity_coded <- with(ukb.ethnic, ifelse(ethnicity == 1, 1,
                                               ifelse(ethnicity == 2, 2, 
                                               ifelse(ethnicity == 3, 3,
                                               ifelse(ethnicity == 4, 4,
                                               ifelse(ethnicity == 5, 4, 
                                               ifelse(ethnicity == 6, 6,
                                               ifelse(ethnicity == -1, -1,
                                               ifelse(ethnicity == -3, -3,
                                               ifelse(ethnicity < 2000, 1,
                                               ifelse(ethnicity < 3000, 2,
                                               ifelse(ethnicity < 4000, 3, 
                                               ifelse(ethnicity < 5000, 4, NA)))))))))))))

ukb_demog <- cbind(ukb_demog[,c(1,3,9,11,13)], ukb.ethnic)
names(ukb_demog) <- c("BMI", "BMI.imaging","Smoking","Smoking.Imaging","Alcohol","eid","ethnicity","ethnicity_coded")

# Transform g
g.factor$g <- 100 - g.factor$g

# Add g to all data frames
ukb.data <- merge(ukb.data, g.factor, by="eid")
ukb.fmri <- merge(ukb.fmri, g.factor, by="eid")
ukb.fmri <- merge(ukb.fmri, ukb.cog, by="eid")

ukb.fmri.demog <- merge(ukb.fmri, ukb_demog, by="eid")

################################################################################
#################### Demographics table ########################################
################################################################################

table_demog <- ukb.fmri.demog %>% 
  dplyr::select(sex, age, education, ethnicity_coded, BMI, BMI.imaging, Smoking, Smoking.Imaging, Alcohol ) %>% # keep only columns of interest
  mutate(
    Alcohol = factor(Alcohol, labels = c("Prefer not to answer","Never", "Monthly or less", "2 to 4 times a month", "2 to 3 times a week", "4 or more times a week")), 
    Smoking = factor(Smoking, labels = c("Prefer not to answer","Never", "Previous", "Current")), 
    Smoking.Imaging = factor(Smoking.Imaging, labels = c("Prefer not to answer","Never", "Previous", "Current")), 
    ethnicity_coded = factor(ethnicity_coded, labels = c("Prefer not to answer", "White", "Mixed", "Asian or Asian British", "Black or Black British", "Other ethnic group")), 
    education = factor(education, labels = c("None of the above", "Prefer not to answer","College or University degree", "A levels/AS levels or equivalent", "O levels/GCSEs or equivalent", "CSEs or equivalent", "NVQ or HND or HNC or equivalent", "Other professional qualifications eg: nursing, teaching")), 
    sex = factor(sex, labels = c("Female","Male")), 
  ) %>% 
  tbl_summary(     
    statistic = list(all_continuous() ~ "{mean}+/-{sd} ({N_nonmiss})",        # stats and format for continuous columns
                     all_categorical() ~ "{p}% ({n}/{N}) "),   # stats and format for categorical columns
    digits = list(all_continuous() ~ 1,  
                  all_categorical() ~ 1),# rounding for continuous columns
    type   = list(c(sex) ~ "dichotomous",
                  c(age, BMI, BMI.imaging) ~ "continuous",
                  c(education, ethnicity_coded, Smoking, Smoking.Imaging, Alcohol) ~ "categorical"),
    value = list(sex ~ "Female"),
    label  = list(                                              # display labels for column names
      age ~ "Age at imaging visit (years)",
      sex ~ "Female",
      education  ~ "Qualifications",
      ethnicity_coded ~ "Ethnic background",
      BMI ~ "BMI (kg/m2)",
      BMI.imaging ~ "BMI at imaging visit (kg/m2)",
      Smoking ~ "Smoking status",
      Smoking.Imaging ~ "Smoking status at imaging visit",
      Alcohol ~ "Frequency of drinking alcohol"
    ),
    missing = "no", # don't list missing data separately;
  ) %>%
  modify_header(label = "Variable") %>% # update the column header
  bold_labels()


table_demog
# Save in Word file
#table_demog %>%
#  as_flex_table() %>%
#  flextable::save_as_docx(path = ".../path/file.docx")


################################################################################
#################### Dummy code variables ######################################
################################################################################

ukb.fmri$age_stand <- scale(ukb.fmri$age)
ukb.fmri$gender_m <- ifelse(ukb.fmri$sex == "1", 1, 0)
ukb.fmri$gender_f <- ifelse(ukb.fmri$sex == "0", 1, 0)

ukb.fmri$edu_uni <- ifelse(ukb.fmri$education == "1", 1, 0)
ukb.fmri$edu_Alevel <- ifelse(ukb.fmri$education == "2", 1, 0)
ukb.fmri$edu_Olevel <- ifelse(ukb.fmri$education == "3", 1, 0)
ukb.fmri$edu_cse <- ifelse(ukb.fmri$education == "4", 1, 0)
ukb.fmri$edu_nvq <- ifelse(ukb.fmri$education == "5", 1, 0)
ukb.fmri$edu_other <- ifelse(ukb.fmri$education == "6", 1, 0)
ukb.fmri$edu_none <- ifelse(ukb.fmri$education == "-7", 1, 0)
ukb.fmri$edu_na <- ifelse(ukb.fmri$education == "-3", 1, 0)

ukb.fmri <- merge(ukb.fmri, ukb.data[,c("eid","volumetric_scaling")], by = "eid")


################################################################################
#################### Regression analysis connectivity ##########################
################################################################################

index <- 1
results.table <- data.frame()

ukb.fmri[,c(15:224)] <- as.matrix(ukb.fmri[,c(15:224)]) %*% diag(sign(colMeans(ukb.fmri[,c(15:224)])))

for (val in ukb.fmri[,c(15:224)]) {
  
  results.table <- UKB.regressionAnalysis(ukb.fmri, ukb.fmri$g, val, 'g', index, results.table)
  
  results.table <- UKB.regressionAnalysis(ukb.fmri, ukb.fmri$F1, val, 'F1', index, results.table)
  results.table <- UKB.regressionAnalysis(ukb.fmri, ukb.fmri$F2, val, 'F2', index, results.table)
  results.table <- UKB.regressionAnalysis(ukb.fmri, ukb.fmri$F3, val, 'F3', index, results.table)
  
  results.table <- UKB.regressionAnalysis(ukb.fmri, ukb.fmri$fluid.intelligence, val, 'fluid.intelligence', index, results.table)
  results.table <- UKB.regressionAnalysis(ukb.fmri, ukb.fmri$matrix.pattern, val, 'matrix.pattern', index, results.table)
  results.table <- UKB.regressionAnalysis(ukb.fmri, ukb.fmri$symbol.digit, val, 'symbol.digit', index, results.table)
  results.table <- UKB.regressionAnalysis(ukb.fmri, ukb.fmri$tower, val, 'tower.rearranging', index, results.table)
  results.table <- UKB.regressionAnalysis(ukb.fmri, ukb.fmri$reaction.time, val, 'reaction.time', index, results.table)
  results.table <- UKB.regressionAnalysis(ukb.fmri, ukb.fmri$pairs.matching, val, 'pairs.matching', index, results.table)
  results.table <- UKB.regressionAnalysis(ukb.fmri, ukb.fmri$pairs, val, 'paired.associate', index, results.table)
  results.table <- UKB.regressionAnalysis(ukb.fmri, ukb.fmri$trail.making, val, 'trail.making', index, results.table)
  results.table <- UKB.regressionAnalysis(ukb.fmri, ukb.fmri$max.digits, val, 'numeric.memory', index, results.table)
  
  index <- index + 1
  
}

for(i in c(1:13)){
  
  p.value <- results.table[,c(3*i)]
  p.value <- p.adjust(as.numeric(unlist(p.value)), method = "fdr")
  results.table[,ncol(results.table)+1] <- p.value
}

F.table <- results.table[,40:52]
Coeff.table <- results.table[,seq(1,39,3)]
Rsquared.table <- results.table[,seq(3,39,3)]

names(F.table) <- c('g.p', 'F1.p', 'F2.p', 'F3.p', 'fluid.intelligence.p', 'matrix.pattern.p', 'symbol.digit.p', 'tower.rearranging.p', 'reaction.time.p', 'pairs.matching.p', 'paired.associate.p', 'trail.making.p', 'numeric.memory.p')

results.table.final <- cbind(Coeff.table, F.table, Rsquared.table)


################################################################################
#################### Circular plots ############################################
################################################################################

plot.data <- data.frame(matrix(nrow = 21, ncol = 0))

# Load fMRI node images from UKB
plot.data$images <- list.files(path = ".../path/Nodes/files", recursive = TRUE, pattern = "*.png", full.names = TRUE)

plot.data$names <- as.character(1:21)
plot.data[c(1,14,20),"gColor"] <- "#e41a1c"
plot.data[c(3, 5),"gColor"] <- "#377eb8"
plot.data[c(17),"gColor"] <- "#4daf4a"
plot.data[c(16),"gColor"] <- "#984ea3"
plot.data[c(9, 21),"gColor"] <- "#ff7f00"
plot.data[c(6, 13),"gColor"] <- "#ffff33"
plot.data[c(10,11,12),"gColor"] <- "#a65628"
plot.data[c(2,4,7,8,19),"gColor"] <- "#f781bf"
plot.data[c(15,18),"gColor"] <- "#999999"

plot.data <- plot.data %>% arrange(gColor)

connectivity.data <-results.table.final

network.names <- names(ukb.fmri[,c(15:224)])

for(ind in c(1:length(network.names))){
  tmp <- strsplit(network.names[[ind]],"[.]")
  connectivity.data[ind, 'Network2'] <- tmp[[1]][2]
  tmp <- tmp[[1]][1]
  connectivity.data[ind, 'Network1'] <- str_sub(tmp,2,-1)
}

lgd_points = Legend(at = c("specific", "shared"),
                    legend_gp = gpar(fill = c("#998ec3", "#f1a340")),
                    title = "Connection", direction = "horizontal", ncol = 2, title_position = "lefttop")

lgd_cluster = Legend(at = c("Default mode", "Attention", "Auditory", "Salience", "Central Executive", "Language","Motor","Visual","Basal Ganglia & Cerebellum"),
                     legend_gp = gpar(fill = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999")),
                     title = "Networks", direction = "horizontal", ncol = 3, title_position = "lefttop")



layout.matrix <- matrix(c(1, 2, 3), nrow = 1, ncol = 3)
layout.matrix
layout(mat = layout.matrix,
       heights = c(0.5), # Heights of the two rows
       widths = c(1, 1, 1)) # Widths of the three columns

reg.results <- results.table.final[,c('F1.coeff','F1.p', 'g.p', 'F2.p', 'F3.p', 'fluid.intelligence.p', 'matrix.pattern.p', 'symbol.digit.p', 'tower.rearranging.p', 'reaction.time.p', 'pairs.matching.p', 'paired.associate.p', 'trail.making.p', 'numeric.memory.p')]
circular.plot(reg.results, plot.data, connectivity.data, "F1 - visuospatial reasoning")

draw(lgd_cluster, x = unit(2, "cm"), y = unit(1.75, "cm"), just = c("left","top"))

reg.results <- results.table.final[,c('F2.coeff','F2.p','g.p', 'F1.p', 'F3.p', 'fluid.intelligence.p', 'matrix.pattern.p', 'symbol.digit.p', 'tower.rearranging.p', 'reaction.time.p', 'pairs.matching.p', 'paired.associate.p', 'trail.making.p', 'numeric.memory.p')]
circular.plot(reg.results, plot.data, connectivity.data, "F2 - verbal-analytical reasoning")

reg.results <- results.table.final[,c('F3.coeff','F3.p','g.p', 'F1.p', 'F2.p', 'fluid.intelligence.p', 'matrix.pattern.p', 'symbol.digit.p', 'tower.rearranging.p', 'reaction.time.p', 'pairs.matching.p', 'paired.associate.p', 'trail.making.p', 'numeric.memory.p')]
circular.plot(reg.results, plot.data, connectivity.data, "F3 - processing speed")

draw(lgd_points, x = unit(30, "cm"), y = unit(1.75, "cm"), just = c("left","top"))




reg.results <- results.table.final[,c('g.coeff','g.p','F1.p', 'F2.p', 'F3.p', 'fluid.intelligence.p', 'matrix.pattern.p', 'symbol.digit.p', 'tower.rearranging.p', 'reaction.time.p', 'pairs.matching.p', 'paired.associate.p', 'trail.making.p', 'numeric.memory.p')]
png(file=".../save/plot.png", width=280, height=320, bg = "transparent")
circular.plot(reg.results, plot.data, connectivity.data, "g - general intelligence")
dev.off()

reg.results <- results.table.final[,c('F1.coeff','F1.p', 'g.p', 'F2.p', 'F3.p', 'fluid.intelligence.p', 'matrix.pattern.p', 'symbol.digit.p', 'tower.rearranging.p', 'reaction.time.p', 'pairs.matching.p', 'paired.associate.p', 'trail.making.p', 'numeric.memory.p')]
png(file=".../save/plot.png", width=280, height=320, bg = "transparent")
circular.plot(reg.results, plot.data, connectivity.data, "F1 - visuospatial reasoning")
dev.off()

reg.results <- results.table.final[,c('F2.coeff','F2.p','g.p', 'F1.p', 'F3.p', 'fluid.intelligence.p', 'matrix.pattern.p', 'symbol.digit.p', 'tower.rearranging.p', 'reaction.time.p', 'pairs.matching.p', 'paired.associate.p', 'trail.making.p', 'numeric.memory.p')]
png(file=".../save/plot.png", width=280, height=320, bg = "transparent")
circular.plot(reg.results, plot.data, connectivity.data, "F2 - verbal-analytical reasoning")
dev.off()

reg.results <- results.table.final[,c('F3.coeff','F3.p','g.p', 'F1.p', 'F2.p', 'fluid.intelligence.p', 'matrix.pattern.p', 'symbol.digit.p', 'tower.rearranging.p', 'reaction.time.p', 'pairs.matching.p', 'paired.associate.p', 'trail.making.p', 'numeric.memory.p')]
png(file=".../save/plot.png", width=280, height=320, bg = "transparent")
circular.plot(reg.results, plot.data, connectivity.data, "F3 - processing speed")
dev.off()

reg.results <- results.table.final[,c('fluid.intelligence.coeff','fluid.intelligence.p','g.p', 'F1.p', 'F2.p', 'F3.p', 'matrix.pattern.p', 'symbol.digit.p', 'tower.rearranging.p', 'reaction.time.p', 'pairs.matching.p', 'paired.associate.p', 'trail.making.p', 'numeric.memory.p')]
png(file=".../save/plot.png", width=280, height=320, bg = "transparent")
circular.plot(reg.results, plot.data, connectivity.data, "Fluid intelligence")
dev.off()

reg.results <- results.table.final[,c('matrix.pattern.coeff','matrix.pattern.p','g.p', 'F1.p', 'F2.p', 'F3.p', 'fluid.intelligence.p', 'symbol.digit.p', 'tower.rearranging.p', 'reaction.time.p', 'pairs.matching.p', 'paired.associate.p', 'trail.making.p', 'numeric.memory.p')]
png(file=".../save/plot.png", width=280, height=320, bg = "transparent")
circular.plot(reg.results, plot.data, connectivity.data, "Matrix pattern completion")
dev.off()

reg.results <- results.table.final[,c('numeric.memory.coeff','numeric.memory.p', 'g.p', 'F1.p', 'F2.p', 'F3.p', 'fluid.intelligence.p', 'matrix.pattern.p', 'symbol.digit.p', 'tower.rearranging.p', 'reaction.time.p', 'pairs.matching.p', 'paired.associate.p', 'trail.making.p')]
png(file=".../save/plot.png", width=280, height=320, bg = "transparent")
circular.plot(reg.results, plot.data, connectivity.data, "Numeric memory")
dev.off()

reg.results <- results.table.final[,c('trail.making.coeff','trail.making.p','g.p', 'F1.p', 'F2.p', 'F3.p', 'fluid.intelligence.p', 'matrix.pattern.p', 'symbol.digit.p', 'tower.rearranging.p', 'reaction.time.p', 'pairs.matching.p', 'paired.associate.p', 'numeric.memory.p')]
png(file=".../save/plot.png", width=280, height=320, bg = "transparent")
circular.plot(reg.results, plot.data, connectivity.data, "Trail making")
dev.off()

reg.results <- results.table.final[,c('pairs.matching.coeff','pairs.matching.p','g.p', 'F1.p', 'F2.p', 'F3.p', 'fluid.intelligence.p', 'matrix.pattern.p', 'symbol.digit.p', 'tower.rearranging.p', 'reaction.time.p', 'paired.associate.p', 'trail.making.p', 'numeric.memory.p')]
png(file=".../save/plot.png", width=280, height=320, bg = "transparent")
circular.plot(reg.results, plot.data, connectivity.data, "Pairs matching")
dev.off()

reg.results <- results.table.final[,c('reaction.time.coeff','reaction.time.p','g.p', 'F1.p', 'F2.p', 'F3.p', 'fluid.intelligence.p', 'matrix.pattern.p', 'symbol.digit.p', 'tower.rearranging.p', 'pairs.matching.p', 'paired.associate.p', 'trail.making.p', 'numeric.memory.p')]
png(file=".../save/plot.png", width=280, height=320, bg = "transparent")
circular.plot(reg.results, plot.data, connectivity.data, "Reaction time")
dev.off()

reg.results <- results.table.final[,c('symbol.digit.coeff','symbol.digit.p','g.p', 'F1.p', 'F2.p', 'F3.p', 'fluid.intelligence.p', 'matrix.pattern.p', 'tower.rearranging.p', 'reaction.time.p', 'pairs.matching.p', 'paired.associate.p', 'trail.making.p', 'numeric.memory.p')]
png(file=".../save/plot.png", width=280, height=320, bg = "transparent")
circular.plot(reg.results, plot.data, connectivity.data, "Symbol digit substitution")
dev.off()

reg.results <- results.table.final[,c('tower.rearranging.coeff','tower.rearranging.p','g.p', 'F1.p', 'F2.p', 'F3.p', 'fluid.intelligence.p', 'matrix.pattern.p', 'symbol.digit.p', 'reaction.time.p', 'pairs.matching.p', 'paired.associate.p', 'trail.making.p', 'numeric.memory.p')]
png(file=".../save/plot.png", width=280, height=320, bg = "transparent")
circular.plot(reg.results, plot.data, connectivity.data, "Tower rearranging")
dev.off()

reg.results <- results.table.final[,c('paired.associate.coeff','paired.associate.p','g.p', 'F1.p', 'F2.p', 'F3.p', 'fluid.intelligence.p', 'matrix.pattern.p', 'symbol.digit.p', 'tower.rearranging.p', 'reaction.time.p', 'pairs.matching.p', 'trail.making.p', 'numeric.memory.p')]
png(file=".../save/plot.png", width=280, height=320, bg = "transparent")
circular.plot(reg.results, plot.data, connectivity.data, "Paired associate learning")
dev.off()

################################################################################
#################### Heatmap plot ##############################################
################################################################################

F.p <- results.table.final[,c('g.p', 'F1.p', 'F2.p', 'F3.p', 'symbol.digit.p', 'tower.rearranging.p', 'matrix.pattern.p', 'fluid.intelligence.p', 'numeric.memory.p', 'paired.associate.p', 'trail.making.p', 'pairs.matching.p', 'reaction.time.p')]
#reorder by column index
colnames(F.p) <-  c('g - general intelligence','F1 - visuospatial reasoning','F2 - verbal-analytical reasoning','F3 - processing speed','Symbol digit substitution', 'Tower rearranging', 'Matrix pattern completion', 'Fluid intelligence', 'Numeric memory', 'Paired associate learning', 'Trail making', 'Pairs matching', 'Reaction time')

connectivity.names <- c("1.2")
for(j in c(3:21)){
  connectivity.names <- append(connectivity.names, paste(c(1:(j-1)), replicate((j-1),j),sep = "."))
}

row.names(F.p) <- connectivity.names

F.coeff <- results.table.final[,c('g.coeff', 'F1.coeff', 'F2.coeff', 'F3.coeff', 'symbol.digit.coeff', 'tower.rearranging.coeff', 'matrix.pattern.coeff', 'fluid.intelligence.coeff', 'numeric.memory.coeff', 'paired.associate.coeff', 'trail.making.coeff', 'pairs.matching.coeff', 'reaction.time.coeff')]
#reorder by column index
colnames(F.coeff) <- c('g - general intelligence','F1 - visuospatial reasoning','F2 - verbal-analytical reasoning','F3 - processing speed','Symbol digit substitution', 'Tower rearranging', 'Matrix pattern completion', 'Fluid intelligence', 'Numeric memory', 'Paired associate learning', 'Trail making', 'Pairs matching', 'Reaction time')

connectivity.names <- c("1.2")
for(j in c(3:21)){
  connectivity.names <- append(connectivity.names, paste(c(1:(j-1)), replicate((j-1),j),sep = "."))
}

row.names(F.coeff) <- connectivity.names

F.coeff <- F.coeff[rowSums(F.p  < 0.05) > 0,]
F.p <- F.p[rowSums(F.p  < 0.05) > 0,]

p.mat <- t(as.matrix(1*F.p))

coeff.mat <- t(as.matrix(1*F.coeff))

col_fun = colorRamp2(c(-0.2, 0, 0.2), hcl_palette = "PuOr")

row_split = rep("Stratum I", 13)
row_split[2:4] = "Stratum II"
row_split[5:13] = "Stratum III"

row_split <- factor(row_split, levels = c("Stratum I","Stratum II", "Stratum III") )

# Load individual connectivity graphics
graphics = list(
  "3.6" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/3_6.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  },
  "3.8" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/3_8.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  },
  "8.12" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/8_12.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  },
  "6.13" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/6_13.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  },
  "6.14" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/6_14.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  },
  "13.14" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/13_14.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  },
  "13.16" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/13_16.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  },
  "7.17" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/7_17.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  },
  "6.18" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/6_18.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  },
  "13.18" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/13_18.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  },
  "5.20" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/5_20.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  },
  "6.20" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/6_20.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  },
  "14.20" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/14_20.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  },
  "5.21" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/5_21.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  },
  "6.21" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/6_21.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  },
  "13.21" = function(x, y, w, h) {
    img = png::readPNG(source=".../path/13_21.png")
    grid.raster(img, x, y = -6, width = unit(2, "cm"), 
                height = unit(2, "cm")*nrow(img)/ncol(img))
  }
)

ht1 = Heatmap(coeff.mat, show_column_dend = FALSE, show_row_dend = FALSE, rect_gp = gpar(col = "white", lwd = 2), col = col_fun, name = "Standardised regression coefficient", 
              cluster_columns = FALSE, row_split = row_split, cluster_rows = FALSE,  cluster_row_slices = FALSE, row_gap = unit(5, "mm"), border = TRUE,
              heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm")), show_column_names = FALSE, 
              bottom_annotation = HeatmapAnnotation(bar = anno_customize(row.names(F.coeff), graphics = graphics, border = FALSE), show_annotation_name = FALSE),
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(p.mat[i, j] < 0.001){
                  grid.text("***", x, y, gp = gpar(fontsize = 20))
                } else if(p.mat[i, j] < 0.01){
                  grid.text("**", x, y, gp = gpar(fontsize = 20))
                } else if(p.mat[i, j] < 0.05){
                  grid.text("*", x, y, gp = gpar(fontsize = 20))
                }
              } )

draw(ht1, heatmap_legend_side = "top", padding = unit(c(75, 2, 2, 2), "mm"))
