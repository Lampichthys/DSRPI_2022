#First you need to install some packages
install.packages('StereoMorph')
install.packages('geomorph')
install.packages('ape')
install.packages('geiger')
install.packages('phytools')
install.packages('convevol')
install.packages('ellipse')
install.packages('cluster')

#Then you need to load them into your workspace
library(StereoMorph)
library(geomorph)
library(ape) #The Swiss army knife of comparative phylogenetics R packages
library(geiger) #A package for "macroevolutionary simulation and estimating parameters related to diversification"
library(phytools) #Very useful for visualization particularly, great blog support
library(phangorn) 
library(convevol)
library(ellipse)
library(cluster)


#----------------STEREOMORPH------------
# Need to set your working directory to the folder that CONTAINS another folder with your images. My working directory is a folder on my desktop, I'm not sure how familiar you are with R, but you can easily get your working directory from the dropdown under 'Session' Set Working Directory and Choose Directory. For me, this folder contains another folder is called "Lanternfish_bodies_straight_and_new_specimens_07_2019" which you can see in the code below, which contains all my images. Your working directory should include a couple other files, a .txt file with the names of your holomogous landmarks. Mine is called landmark_names.txt and you can see the format in the copy I gave to you. You need a second .txt file with the names of the curves you want to make, these include the names of the homologous landmarks that anchor the curves and the name of the curve itself. Mine is called curves.txt. You can look at mine as an example and read through the stereomorph guide for more information. You should also have an empty folder that stereomorph will add shape files to, mine is called 'Shapes'.


# This bit of code will open Stereomorph from the same working directory you are still in. the image.file should be the folder name that contains your images. The landmarks.file is the folder name that I created when I converted my tps files. The landmarks.ref is those landmark names, the curves.ref are my curve names with associated anchor landmarks, and the shapes.file is the folder i want stereomorph to save my information to.

digitizeImages(image.file='Lanternfish_bodies_straight_and_new_specimens_07_2019', landmarks.file = 'landmarks_folder', landmarks.ref = 'landmark_names.txt', curves.ref = "curves.txt", shapes.file = 'Shapes')

#-----------------GEOMORPH------------------------
# Note: Before you start, check to make sure that you have the same number of specimens in both your shapes folder and your master csv file.

# CSV file containing all metadata and variables (e.g., Age, Locality_Collected, Plastic_in_Belly, etc) of interest for each specimen that you have shape data for.

# NOTE: Names of individuals in this csv file MUST match the names of the shape files to which those individuals correspond. Otherwise, name matching portion of the code will not work. 

# Set Working Directory that contains CSV file (I usually just have it in the same folder we have been working in this entire time)

# Here is another way to get to you associated directory other than the dropdown
setwd("~/Desktop/Lanternfish_Bodies_for_pub/")

# This code reads in your csv spreadsheet file to R. Header means you have headers on each column, and sep is asking what your file uses as separators, it is usually either commas or tabs, in my case it was commas. I have named the associated information 'Data' in my workspace.
Data <- read.csv("Phylogeny_BOdy_shape_spreadsheet_straight_04_2021.csv", header = TRUE, sep = ",")

# Just checking my data I just read in
summary(Data)

# Import the shape files you created in stereomorph. My 514 specimens take ~ 3 minutes. It is associating with the 'Shapes' folder they are all housed in.

shapesall <- readShapes("Shapes")

#Now, here is where you assign curves and the number of semilandmarks. 

shapesGM <- readland.shapes(shapesall, nCurvePts = c(15, 8, 7, 10))

#-----------
# Perform General Procrustes Analysis (GPA) - If you prefer to use bending energy, rather than Procrustes Distances - change ProcD = "False." This may be a good idea with very large datasets like mine. 


GPA <- gpagen(shapesGM, ProcD = FALSE)

# Plotting the output, shows the average body shape with landmarks in black and the variation around those landmarks in grey.

plot(GPA)

#PCA just to get the initial splines, WILL NOT BE USING FOR ANALYSIS, MAYBE SUPPLEMENTAL FILE
PCA <- gm.prcomp(GPA$coords)
plot(PCA)

msh <- mshape(GPA$coords)

# PC1
plotRefToTarget(PCA$shapes$shapes.comp1$min, msh)
plotRefToTarget(msh, PCA$shapes$shapes.comp1$max)

#PC2
plotRefToTarget(PCA$shapes$shapes.comp2$min, msh)
plotRefToTarget(msh, PCA$shapes$shapes.comp2$max)

#PC3
plotRefToTarget(PCA$shapes$shapes.comp3$min, msh)
plotRefToTarget(msh, PCA$shapes$shapes.comp3$max)

#PC4
plotRefToTarget(PCA$shapes$shapes.comp4$min, msh)
plotRefToTarget(msh, PCA$shapes$shapes.comp4$max)

#This is where we will match the GPA aligned shape files to the csv dataset that we just imported. 

#Save GPA names to Global Environment for matching
spec.names <- dimnames(GPA$coords)[[3]]
spec.names

#Save Data name to GE for matching
data.names <- Data$ID
data.names

#Now, resolve the two data sets and combine them
name.match <- match(data.names,spec.names)
name.match

Data2 <- Data[name.match, ]

# Check for NAs. IF you have them, either the dataset and the shapes folder have a different number of specimens or there are misspellings. 

summary(Data2)

#--------------------STATISTICS and ALLOMETRY OLD---------------
#Create your geomorph data frame
GDF <- geomorph.data.frame(GPA, shape = GPA$coords, CS = GPA$Csize, ID = Data2$ID, Sub = Data2$Subfamily, Gen = Data2$Genus, SP_AV = Data2$SP_AV, Bino = Data2$Binomail, CLightO = Data2$CLO)

#making a vector for group identities for grouping later
GP <- Data2$Subfamily
GPG <- Data2$Genus
GPS <- Data2$SP_AV
Bino <- Data2$Binomail
CLO <- Data2$CLO

#----------------- STATISTIC and ALLOMETRY NEW --------------------

#First average the shape data and centroid sizes by species
average.Csize <- aggregate(GDF$Csize, list(GDF$SP_AV), FUN=mean)


#lets try something, mean scores GPA coords FOR species
             
array.data.sp <- GDF$shape
             
#change our 3D array into a 2D array in order to calculate by group
x <- two.d.array(array.data.sp)
             
p <- dim(array.data.sp)[1] # the number of landmarks
k <- dim(array.data.sp)[2] # the dimensions of the coordinates
Y <- array(data = NA, dim = c(p, k, length(levels(GDF$SP_AV)))) #new empty array to fill
dimnames(Y)[[3]] <- levels(GDF$SP_AV)# set group levels as new names
             
             
means.data <- rowsum(x, GDF$SP_AV)/as.vector(table(GDF$SP_AV))# rowsum() is a simple base function to get the sum of the rows, while table() is a great function for getting group sizes (group n).

# then all you have to do is put the averaged data back into an 3D array.
Y <- arrayspecs(means.data, dim(array.data.sp)[1], dim(array.data.sp)[2]) 

Ave.simp.Data <- read.csv("Phylogeny_simple_11_2021.csv", header = TRUE, sep = ",")


GDF_averaged <- geomorph.data.frame(shape = Y, CS = average.Csize, CLO.simp = Ave.simp.Data$CLO)

GP.adj <- Ave.simp.Data$Subfamily
GPG.adj <- Ave.simp.Data$Genus
GPS.adj <- Ave.simp.Data$SP_AV
Bino.adj <- Ave.simp.Data$Binomail
CLO.adj <- Ave.simp.Data$CLO

#read in the tree
Myctophiform.tree <- read.nexus(file = "Total_data_tree_ABV_names.nwk") 
Myctophiform.tree <- ladderize(Myctophiform.tree, right = TRUE)

plot(Myctophiform.tree, cex = 0.5)

# Look at node to extract tree sans outgroups
nodelabels()

# Node 88 and beyond contains our desired tree
Myctophiform.tree <- extract.clade(Myctophiform.tree, 88, root.edge = 0, collapse.singles = TRUE, interactive = FALSE)

plot(Myctophiform.tree, cex = 1)

#use the procD.pgls function in Geomorph to run the allometric regression on those averaged values. 
shape_regression <- procD.pgls(shape ~ log(CS.x), phy = Myctophiform.tree, data = GDF_averaged, iter = 9999)
summary(shape_regression)


#take the residuals from that regression as your allometrically adjusted data for ALL subsequent analyses (PCA/phylomorphospace, modularity, disparity)

# For allometrically and phylogenetically-corrected shape data we are using the residuals from 'shape_regression' and will use these going forward.
size_and_phylogeny_adjusted_shape <- shape_regression$pgls.residuals

#It is stored in a matrix, we need to make it a 3D array for our PCA. 43 is the number of landmarks and 2 is the number of dimensions.
size_and_phylogeny_adjusted_shape.array <- arrayspecs(size_and_phylogeny_adjusted_shape, 43, 2) # then all you have to do is put the  data back into an 3D array.


PCA.Adj <- gm.prcomp(size_and_phylogeny_adjusted_shape.array)
plot(PCA.Adj)

summary(PCA.Adj)


lantern.plot <- plot(PCA.Adj)

x.locations <- PCA.Adj$x[,1]
y.locations <- PCA.Adj$x[,2]


text(x.locations, y.locations, labels = GDF_averaged$CS.Group.1)

#PC3,4

plot(PCA.Adj$x[,3],PCA.Adj$x[,4])

x.locations2 <- PCA.Adj$x[,3]
y.locations2 <- PCA.Adj$x[,4]

text(x.locations2, y.locations2, labels = GDF_averaged$CS.Group.1)

#---------------- PHYLOGENETIC MANOVA ---------------------
shape_regression_CLO <- procD.pgls(shape ~ log(CS.x) * CLO.adj, phy = Myctophiform.tree, data = GDF_averaged, iter = 9999)
summary(shape_regression_CLO)

shape_regression_sub <- procD.pgls(shape ~ log(CS.x) * GP.adj, phy = Myctophiform.tree, data = GDF_averaged, iter = 9999)
summary(shape_regression_sub)


# ------------------------ PHYLOMORPHOSPACE ------------------
phylo.loc <- cbind(x.locations, y.locations)
phylo.loc

cols <- rep("black",length(Myctophiform.tree$tip.label) + Myctophiform.tree$Nnode)

# Coloring all the nodes beyond a certain ancestral node
names(cols) <- 1:length(cols)
cols[getDescendants(Myctophiform.tree,83)] <- "red"
# and everything from "36" blue:
cols[getDescendants(Myctophiform.tree,112)] <- "blue"
# finally, these can even be nested
cols[getDescendants(Myctophiform.tree,130)] <- "yellow"
cols[getDescendants(Myctophiform.tree,148)] <- "green"
cols[getDescendants(Myctophiform.tree,156)] <- "orange"


phylomorphospace(Myctophiform.tree, phylo.loc, control=list(col.node=cols))

#-------------------------- DISPARITY ------------------------
# Test whether or which genera are significantly different in their shape variability to each other

MD.bycoords.gen.adj <- morphol.disparity(shape ~ GPG.adj, groups = ~ GPG.adj, data = GDF_averaged, iter=10000)

MD.bycoords.gen.adj
write.csv(MD.bycoords.gen.adj$PV.dist, file = "MDisparity.values.gen.adj.csv")
write.csv(MD.bycoords.gen.adj$PV.dist.Pval, file = "MDisparity.Pvals.adj.gen.csv")

# Test whether or which subfamilies are significantly different in their shape variability to each other

MD.bycoords.sub.adj <- morphol.disparity(shape ~ GP.adj, groups = ~ GP.adj, data = GDF_averaged, iter=10000)

MD.bycoords.sub.adj
write.csv(MD.bycoords.sub.adj$PV.dist, file = "MDisparity.values.sub.adj.csv")
write.csv(MD.bycoords.sub.adj$PV.dist.Pval, file = "MDisparity.Pvals.sub.adj.csv")

MD.bycoords.CLO.adj <- morphol.disparity(shape ~ CLO.adj, groups = ~ CLO.adj, data = GDF_averaged, iter=10000)
MD.bycoords.CLO.adj

# ------------------------- MODULARITY ---------------------------
#Since the modularity function doesn't let you add a categorical variable, one option is to do separate analyses for species with and without caudal light organs (the tree would also need to be trimmed accordingly for each analysis). If light organ presence impacted evolutionary modularity, you would expect the two analyses to show different patterns (i.e., in their relative degree of modularity).

#TOTAL TEST
#landmarks on the body and tail
land.gps.tail <- rep('a', 43); land.gps.tail[6:9] <- 'b'; land.gps.tail[31:43] <- 'b'

MT.Tail <- phylo.modularity(GDF_averaged$shape, land.gps.tail, phy = Myctophiform.tree, CI = FALSE, iter = 10000)
summary(MT.Tail) # Test summary
plot(MT.Tail) # Histogram of CR sampling distribution 

#TRIM TEST SEPERATE ANALYSES 

#THOSE WITH ORGAN
#Trim Tree
Myctophiform.tree.organ <- drop.tip(Myctophiform.tree, c(14,24,17,18,22,15,16,26,21,23,25,13,20,12,19,35,34,36,48,71,69,61), trim.internal = TRUE, subtree = FALSE, root.edge = 0, rooted = is.rooted(Myctophiform.tree), collapse.singles = TRUE, interactive = FALSE)

#Trim data
organ.data <- Y[, , c(16,17,22:75,79)]
  
MT.Tail.ORGAN <- phylo.modularity(organ.data, land.gps.tail, phy = Myctophiform.tree.organ, CI = FALSE, iter = 10000)
summary(MT.Tail.ORGAN)

#THOSE WITHOUT ORGAN
#Trim Tree
Myctophiform.tree.noorgan <- drop.tip(Myctophiform.tree, c(1,2,27,28,41,29,30,53,67,68,37,39,38,58,56,55,73,74,75,57,77,52,51,54,11,33,32,6,7,50,40,70,63,62,65,64,45,43,44,76,42,3,5,4,49,8,10,9,79,78,72,66,47,60,46,59,31), trim.internal = TRUE, subtree = FALSE, root.edge = 0, rooted = is.rooted(Myctophiform.tree), collapse.singles = TRUE, interactive = FALSE)

#Trim data
noorgan.data <- Y[, , c(1:15,18:21,76:78)]

MT.Tail.NOORGAN <- phylo.modularity(noorgan.data, land.gps.tail, phy = Myctophiform.tree.noorgan, CI = FALSE, iter = 10000)
summary(MT.Tail.NOORGAN)

#Compare
CR.comparisons.CLO <- compare.CR(MT.Tail.ORGAN, MT.Tail.NOORGAN, CR.null = TRUE)


#--------------------- PHYLOGENETIC SIGNAL ---------------
#Signal in those with organ
phylo.signal.organ <- physignal(A = organ.data, Myctophiform.tree.organ, iter = 10000, seed = NULL, print.progress = TRUE)

summary(phylo.signal.organ)

plot(phylo.signal.organ)

#Signal in those without organ
phylo.signal.noorgan <- physignal(A = noorgan.data, Myctophiform.tree.noorgan, iter = 10000, seed = NULL, print.progress = TRUE)

summary(phylo.signal.noorgan)

plot(phylo.signal.noorgan)


# ----------------- EVOLUTIONARY RATES ------------------
#FIRST need to transform tree into an ultrametric tree

tree_orthofinder <- read.nexus(file = "Total_data_tree_ABV_names.nwk")

force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}

tree_orthofinder_ultra <- force.ultrametric(tree_orthofinder)
write.tree(tree_orthofinder_ultra, file = "tree_ultrametric.txt", append = FALSE, digits = 10, tree.names = FALSE)
is.ultrametric(tree_orthofinder_ultra)
plot(tree_orthofinder_ultra)
nodelabels()
tree_orthofinder_ultra <- extract.clade(tree_orthofinder_ultra, 88, root.edge = 0, collapse.singles = TRUE, interactive = FALSE)

#Pooled the reduced_absent into the Absent group for all analyses
CLO.data <- Ave.simp.Data$CLO
names(CLO.data) <- Ave.simp.Data$SP_AV

Rates.test <- compare.evol.rates(GDF_averaged$shape, phy = tree_orthofinder_ultra, gp = CLO.data, iter = 10000)
summary(Rates.test)

#Curious about rates by genus
Genera.data <- Ave.simp.Data$Genus
names(Genera.data) <- Ave.simp.Data$SP_AV

Rates.genera <- compare.evol.rates(GDF_averaged$shape, phy = tree_orthofinder_ultra, gp = Genera.data, iter = 10000)
summary(Rates.genera)
plot(Rates.genera)

#------------------ CONVERGENCE -------------------
#test for convergence within a morphospace

convtips <- as.character(Myctophiform.tree.noorgan$tip.label)


convergent.data <- cbind(x.locations, y.locations)


conv.test.3 <- convratsig(Myctophiform.tree, convergent.data, convtips, nsim = 1)

conv.test.3

Myctophiform.tree.conver.myctoph <- drop.tip(Myctophiform.tree.organ, c(49,43,45,44,27,42,25,54,26,24,3,5,4,30,8,10,9,50,57,56,46,29,41,28,40,22,31,2,12,23,21,20,39,37,36,51,52,53,38,55,33,32,16,35,11,18,17,6,7))
plot(Myctophiform.tree.conver.myctoph)

convtips.compressed <- as.character(Myctophiform.tree.conver.myctoph$tip.label)


conv.test.4 <- convratsig(Myctophiform.tree, convergent.data, convtips.compressed, nsim = 2)
conv.test.4

Myctophiform.tree.noorgan.small <- drop.tip(Myctophiform.tree.noorgan, c(20,22,17,3,13,6,7,4,5,15,10,12,14,2,9,1,8,19))

convtips.small <- as.character(Myctophiform.tree.noorgan.small$tip.label)

conv.test.small <- convratsig(Myctophiform.tree, convergent.data, convtips.small, nsim = 100)
conv.test.small
