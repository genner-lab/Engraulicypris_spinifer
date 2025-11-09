#Code for Morphological Analyses of Engraulicypris

#Step 1: Set working directory and load data

#Set appropriate working directory

Freya_Morphology <- read.table("Engraulicypris_Morphology.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

#Load packages
library(VIM)
library(MASS)
library(ggplot2)
library(dplyr)
library(ggpubr)

#Step 2: Impute missing data and derive residuals for 8 morphological variables, selecting size standardisation variables that maximise body shape variation

Freya_Morphology2 <- kNN(Freya_Morphology[,c(3:11)], k = 5)
Freya_Morphology2 <- cbind(Freya_Morphology[,c(1:2)], Freya_Morphology2[,c(1:9)])

HL_SR <- lm(log10(HL) ~ log10(BD), data = Freya_Morphology2)
ED_SR <- lm(log10(ED) ~ log10(BD), data = Freya_Morphology2)
BD_SR <- lm(log10(BD) ~ log10(SL), data = Freya_Morphology2)
CPL_SR <- lm(log10(CPL) ~ log10(BD), data = Freya_Morphology2)
CPD_SR <- lm(log10(CPD) ~ log10(BD), data = Freya_Morphology2)
PDL_SR <- lm(log10(PDL) ~ log10(BD), data = Freya_Morphology2)
PVL_SR <- lm(log10(PVL) ~ log10(BD), data = Freya_Morphology2)
SNL_SR <- lm(log10(SNL) ~ log10(BD), data = Freya_Morphology2)

Example_residuals <- data.frame(HL_SR$residuals, ED_SR$residuals, 
                                BD_SR$residuals, CPL_SR$residuals,
                                CPD_SR$residuals, PDL_SR$residuals,
                                PVL_SR$residuals, SNL_SR$residuals)

Freya_Morphology <- cbind(Freya_Morphology$Code, Freya_Morphology$Group, Example_residuals)
colnames(Freya_Morphology)[1] <- "Code"
colnames(Freya_Morphology)[2] <- "Group"

#clean house
rm(Freya_Morphology2, Example_residuals, HL_SR, ED_SR, BD_SR, CPL_SR, CPD_SR, PDL_SR, PVL_SR, SNL_SR)

#Step 3: Linear Discriminant Analysis

LDA  <-  lda(Group ~ HL_SR.residuals + ED_SR.residuals + 
               BD_SR.residuals + CPL_SR.residuals + CPD_SR.residuals + 
               PDL_SR.residuals + PVL_SR.residuals + SNL_SR.residuals, data = Freya_Morphology)
LDA_predict <- predict(LDA)
LDA_score <- as.data.frame(LDA_predict)
Freya_Morphology$LD1 <- LDA_score$x.LD1
Freya_Morphology$LD2 <- LDA_score$x.LD2

#For the loadings of the variables on the axes
LDA$scaling

#For the the variation captured, so convert to percentages (24.99582/28.19506*100 = 88.65%)
LDA$svd 
#so convert to percentages (24.99582/28.19506*100 = 88.65%)
#so convert to percentages (3.19924/28.19506*100 = 11.35%)

#clean house
rm(LDA_predict, LDA_score, LDA)


# Calculate convex hulls for each group 
LDA_find_hull <- function(LDA_score) LDA_score[chull(Freya_Morphology2$LD1, Freya_Morphology2$LD2), ] 
LDA_hulls <- LDA_score %>% group_by(class) %>% do(LDA_find_hull(.))

#Step 4: A plot

Plot_LDA <- ggscatter(Freya_Morphology, x = "LD1", y = "LD2",
                      color = "Group", shape = "Group",
                      palette = c("cyan", "blue2", "darkgoldenrod1"),
                      ellipse = TRUE, ellipse.type = "convex",
                      mean.point = FALSE, star.plot = FALSE) +
  theme_classic() +
  labs(x = "LD1 (88.65% of variation)", y = "LD2 (11.35% of variation)") +
  theme(legend.position = "right")

Plot_LDA

#Step 5: Testing for differences

#Global tests with all three groups

ManovaTest_Global <- manova(cbind(HL_SR.residuals, ED_SR.residuals, 
                           BD_SR.residuals, CPL_SR.residuals,
                           CPD_SR.residuals, PDL_SR.residuals,
                           PVL_SR.residuals, SNL_SR.residuals) ~ Group, data = Freya_Morphology)
summary(ManovaTest_Global)

#Pairwise tests, comparing E. spinifer Lake Malawi vs. E. spinifer Malagarasi

Freya_Morphology_Pair1 <- Freya_Morphology %>% filter(Group != "E.sardella_LakeMalawi")
ManovaTest_Pair1 <- manova(cbind(HL_SR.residuals, ED_SR.residuals, 
                                  BD_SR.residuals, CPL_SR.residuals,
                                  CPD_SR.residuals, PDL_SR.residuals,
                                  PVL_SR.residuals, SNL_SR.residuals) ~ Group, data = Freya_Morphology_Pair1)
summary(ManovaTest_Pair1)

#Pairwise tests, comparing E. spinifer Lake Malawi vs. E. sardella

Freya_Morphology_Pair2 <- Freya_Morphology %>% filter(Group != "E.spinifer_LakeMalawi")
ManovaTest_Pair2 <- manova(cbind(HL_SR.residuals, ED_SR.residuals, 
                                 BD_SR.residuals, CPL_SR.residuals,
                                 CPD_SR.residuals, PDL_SR.residuals,
                                 PVL_SR.residuals, SNL_SR.residuals) ~ Group, data = Freya_Morphology_Pair2)
summary(ManovaTest_Pair2)

#Pairwise tests, comparing E. spinifer Malagarasi vs. E. sardella

Freya_Morphology_Pair3 <- Freya_Morphology %>% filter(Group != "E.spinifer_Malagarasi")
ManovaTest_Pair3 <- manova(cbind(HL_SR.residuals, ED_SR.residuals, 
                                 BD_SR.residuals, CPL_SR.residuals,
                                 CPD_SR.residuals, PDL_SR.residuals,
                                 PVL_SR.residuals, SNL_SR.residuals) ~ Group, data = Freya_Morphology_Pair3)
summary(ManovaTest_Pair3)

#end of code

