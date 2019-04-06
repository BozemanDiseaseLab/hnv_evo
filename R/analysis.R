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
system.file("tex", "texshade.sty", package="msa")

load(file = 'data/df.Rdata')

#restrict to nucleocapsid genes
#first question I have is why is so much of the data the Nucleocapsid protein, where are the Fs and Gs? 
#Is this a problem with how extracted the data?
#Or is this a function of sequencing? Good lord need to figure that one out for sure. 
#SECOND, what happened to the new kerala sequences?
df <- df %>%
  filter(gene == '(N)') #%>%
  #filter(virus == 'Nipah') 
  
#df.seq <- df[1:15,]
seq <- Biostrings::DNAStringSet(df[,9], use.names = TRUE)
seq@ranges@NAMES <- as.character(paste(df$virus, df$ID))

#rm(NucCapsidAlignment)
NucCapsidAlignment <- msa::msa(seq, "Muscle")

print(NucCapsidAlignment, show="complete")
#print(myAlignment, show="complete")

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
d <- dist.alignment(NucCapsidAlignment2, "identity")
hist(d)
#as.matrix(d)[2:5, "NucCapsidAlignment2", drop=FALSE] #dont kow why this doesn't work
heatmap(d)
library(ape)
hemoTree <- njs(d)
plot(hemoTree,type ='fan', main="Nipah and Hendra")

m=as.matrix(d)
library(lattice) 
levelplot(m[1:ncol(m),ncol(m):1])

#a couple indices look really really bad

#the level plot makes me want to remove the following indices, or at least check them out

#193, 109, 111,112,105,113,,110,114,115,108,106,107,

#df.seq <- df[1:15,]
rm <- df[df$ID %in% c(193, 109, 111,112,105,113,110,114,115,108,106,107),]

#i dont see any reason to remove these yet?




