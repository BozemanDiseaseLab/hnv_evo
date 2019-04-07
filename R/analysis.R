#hnv_data analyze
#first run "hnv_search_clean" and "Hnv_search_locations" 

#packages:

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("msa")

#trying to follow this outline:

#https://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf
library(msa)
library(tidyverse)

load(file = 'Data/NucCapsidAlignment.chiroptera')
load(file = 'Data/NucCapsidAlignment')

##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################

#this will only work if you take a subset of the number of viruses, it can't
#handle too many alignments

seq.1 <- Biostrings::DNAStringSet(df[1;15,9], use.names = TRUE)
seq.1@ranges@NAMES <- as.character(df[1:15,'virus')

NucCapsidAlignment.1 <- msa::msa(seq.1, "Muscle")

msaPrettyPrint(NucCapsidAlignment.1, output="asis", askForOverwrite=FALSE)
tmpFile <- tempfile(pattern="msa", tmpdir=".", fileext=".pdf")
tmpFile
msaPrettyPrint(NucCapsidAlignment.1, file=tmpFile, output="pdf",
               showNames="left", showNumbering="none", showLogo="top",
               showConsensus="bottom", logoColors="rasmol",
               verbose=FALSE, askForOverwrite=FALSE)

library(Biobase)
openPDF(tmpFile)

##########################################################################
##########################################################################

#dist alignment and tree
NucCapsidAlignment2 <- msaConvert(NucCapsidAlignment, type="seqinr::alignment")

library(seqinr)
genetic.distance <- dist.alignment(NucCapsidAlignment2, "identity")
hist(genetic.distance)
#as.matrix(d)[2:5, "NucCapsidAlignment2", drop=FALSE] #dont kow why this doesn't work

#change order
ordering <- sort(attr(genetic.distance, "Labels"))
genetic.distance <- as.matrix(genetic.distance)[ordering, ordering]
genetic.distance  <- as.dist(genetic.distance)

m=as.matrix(genetic.distance)
library(lattice) 
levelplot(m[1:ncol(m),ncol(m):1])

library(ape)
hemoTree <- njs(genetic.distance)
plot(hemoTree,type ='fan', main="Nipah and Hendra Nucleocapsid Tree")

#a couple indices look really really bad
#the level plot makes me want to remove the following indices, or at least check them out
#193, 109, 111,112,105,113,,110,114,115,108,106,107,

#df.seq <- df[1:15,]
rm <- df[df$ID %in% c(193, 109, 111,112,105,113,110,114,115,108,106,107),]

#i dont see any reason to remove these yet?

############# spatial analysis ###############
############# spatial analysis ###############
############# spatial analysis ###############
############# spatial analysis ###############
############# spatial analysis ###############
############# spatial analysis ###############
############# spatial analysis ###############
############# spatial analysis ###############
############# spatial analysis ###############
############# spatial analysis ###############

load(file = 'Data/hnv.dist.Rdata')

hnv.dist <- as.matrix(hnv.dist)[ordering, ordering]
hnv.dist  <- as.dist(hnv.dist)
sum(attr(hnv.dist, "Labels") != attr(genetic.distance, "Labels"))

hnv.dist.2 <- as.matrix(hnv.dist)
levelplot(hnv.dist.2[1:ncol(hnv.dist.2),ncol(hnv.dist.2):1])

library(phylin)

gv <- gen.variogram(x = hnv.dist, y= genetic.distance,lag = .1)
gv %>% plot()
gv.plot <- gv.model(gv)
plot(gv)



cutoff <- max(hnv.dist, na.rm=TRUE) / 3 # default maximum distance
num.bins <-  15
bin.width <- cutoff / 15

#genetic.distance[upper.tri(genetic.distance,diag=TRUE)] <- NA

vario.dat <- data.frame(dist = as.numeric(hnv.dist), diff = as.numeric(genetic.distance)) %>%
  dplyr::arrange(desc(dist))

vario.dat <- vario.dat %>% 
  mutate(bin = floor(dist / bin.width) + 1) 

vario.dat2 <- vario.dat %>% 
  #filter(bin < 20) %>%
  group_by(bin) %>% 
  dplyr::summarize(emp.sv = .5 * mean(diff), n = n()) 

vario.dat2 %>% 
  ggplot(aes(x=bin, y= emp.sv))  + 
  geom_point() 











