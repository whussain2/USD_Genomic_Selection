library(devtools)
install_github("TheRocinante-lab/TrainSel")
library(TrainSel)
library(TrainSel)
library(factoextra)
library(FactoMineR)
library(ggthemes)
library(ggplot2)
library(dplyr)
rm(list=ls())
geno.oyt <- readRDS("/Volumes/GoogleDrive/My Drive/ICAR_IRRI_Vikas/Hands_Genomic_Selections/Training_Predictions/data/geno.oyt.filtered.rds")
dim(geno.oyt)

pheno<-readRDS("/Volumes/GoogleDrive/My Drive/ICAR_IRRI_Vikas/Hands_Genomic_Selections/Training_Predictions/data/pheno.filtered.rds")
control=TrainSelControl()
control$niterations=5
# When using a mixed model based CDmin statistic we need to prepare the data using the ’MakeTrainSel-Data function
# We will use marker data to prepare data for mixed model based criterion.
TSData<-MakeTrainSelData(M=geno.oyt)

## We use ’TrainSel’ function to select an Unordered Sample ("UOS") of size 10 from the first 100 individuals in the dataset
#as a default it picks, corresponding to the first 100 rows of data. The default settings do not need to be modified for small to medium-sized optimization problems. We select settypes as "OS"power if to perform multi-objective optimization, else we will use "UOS".

out1<-TrainSel(Data=TSData, Candidates=list(1:616), setsizes=c(150),settypes="UOS", control=
                 control, Stat = NULL)


# Subset the best complete genotypic data of the best picked 150 genotypes from the genotype file.
TS <- out1$BestSol_int

# Subset the best complete genotypic data of the best picked 150 genotypes from the genotype file
geno.train<- as.matrix(geno.oyt[TS,])
saveRDS(geno.trian, file="/Volumes/GoogleDrive/My Drive/ICAR_IRRI_Vikas/Hands_Genomic_Selections/Training_Predictions/data/geno.train.rds")
dim(geno.train)

# Now let us check the representation of trining set
***
  # Create Relationship Matrix
  ***
  
  * Here in this section we will create a **Genomic Relationship Matrix (GRM)** using marker data. The GRM will be based on **VanRanden (2008)**. The Steps used to create this GRM is:
  - Create a center of marker data (X matrix)
- Create a Cross Product $(XX)$
  - Divide the $(XX)$ by number of markers

\[
  GRM= XX^t/m
  \]

* More on relationship matrix can be found here [Source 1](https://www.sciencedirect.com/science/article/pii/S0022030209707933), [Source2](https://naldc.nal.usda.gov/download/27183/PDF)
# First Get genomic matrix
Xs <- scale(geno.oyt, center = TRUE)
# Construct G matrix
GM <- geno.oyt %*% t(geno.oyt)/ncol(geno.oyt)
dim(GM)

# Now perform PCA
# First let us get the lines that are selected

names<-data.frame(name=ifelse(row.names(geno.oyt)%in%row.names(geno.train), "Un-selected", "Selected"))  
pca<- PCA(GM, graph = FALSE)
# Now before plotting, let us add information about top selected lines here
# Add top lines
# Arrange the BLUPs in decreasing order
# Now let us plot the biplot
# Subset the genotype file

#png(file = "./Outputs/PCA_selected.png", width =8, 
   # height =8, units = "in", res = 600)
fviz_pca_biplot(pca, palette = "jco", axes = c(1, 2), geom="point",pointsize = 3,
                addEllipses = FALSE,col.ind = names$name,
                geom.var= "none", label = "none", legend.title = "Group")+
  theme_minimal()+
  # add and modify the title to plot
  theme (
    plot.title = element_blank(),
    # add and modify title to x axis
    axis.title.x = element_text(color="black", size=12), 
    # add and modify title to y axis
    axis.title.y = element_text(color="black", size=12)) +
  # modify the axis text
  theme(axis.text= element_text(color = "black", size = 12)) 
#dev.off()


# Now perform GBLUP in rrBLUP Package R


y.trian<-pheno[pheno$Designation%in%row.names(geno.train), ]
y.test<-pheno[!pheno$Designation%in%row.names(geno.train), ]
dim(y.trian)
geno.test<-geno.oyt[!row.names(geno.oyt)%in%row.names(geno.train), ]
dim(geno.test)
dim(geno.train)

XTRN<-marker[-tst,] ; yTRN<-pheno[-tst,]
XTST<-marker[tst,] ; yTST<-pheno[tst,]



library(rrBLUP)

fm<-mixed.solve(y=y.trian$MPI,Z=as.matrix(geno.train))

######################################
# Marker effects are treated random  #
######################################
mean<-as.vector(fm$beta)
predicted.values<-as.matrix(geno.test)%*%as.vector(fm$u)
predicted.values<-data.frame(GEBVs=predicted.values+mean)
predicted.values$Designation<-row.names(predicted.values)

# Get Correlation accuracy
# First merge with Test phenotypic data
cor.data<-merge(predicted.values, y.test[, c(1,4)], by="Designation")
cor.data<-na.omit(cor.data)
cor(cor.data$GEBVs,cor.data$MPI)
plot(cor.data$GEBVs,cor.data$MPI, main = "Correlation Accuracy",
     xlab = "GEBVs", ylab = "Phenotypic",
     pch = 19, frame = FALSE)
# Get mean Standard errors
mse<-((cor.data$MPI-cor.data$GEBVs)^2)/length(y.test)

# Ranking of lines based on GEBVs

# Ranking and selection of top performing lines
colnames(predicted.values)
predicted.values<-predicted.values%>%arrange(desc(GEBVs))
predicted.values<-predicted.values[1:40, ]
   # Draw the plot
   ggplot(data=predicted.values, aes(x=Designation, y=GEBVs)) +
   geom_bar(stat="identity", width=0.5)+
  theme_classic()+
   labs(title="BLUPs of Top Ranked Genotypes along with Checks",x="designation", y = "Breeding Value")+
   #scale_y_continuous(limits = c(0, 6000), breaks = seq(0, 6000, by = 500))+
  theme (plot.title = element_text(color="black", size=1, face="bold", hjust=0),
               axis.title.x = element_text(color="black", size=10, face="bold"),
               axis.title.y = element_text(color="black", size=10, face="bold")) +
  theme(axis.text= element_text(color = "black", size = 8))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
     
bar.plot

