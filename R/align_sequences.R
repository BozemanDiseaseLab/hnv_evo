#align sequences
library(msa)
library(tidyverse)

load(file = 'data/df.Rdata')

#restrict to nucleocapsid genes
#first question I have is why is so much of the data the Nucleocapsid protein, where are the Fs and Gs? 
#Is this a problem with how extracted the data?
#Or is this a function of sequencing? Good lord need to figure that one out for sure. 
#SECOND, what happened to the new kerala sequences?
df <- df %>%
  filter(gene == '(N)') #%>%
#filter(virus == 'Nipah') 

unique(df$host)
#df.seq <- df[1:15,]
seq <- Biostrings::DNAStringSet(df[,9], use.names = TRUE)
seq@ranges@NAMES <- as.character(paste(df$virus, df$ID))
#rm(NucCapsidAlignment)
NucCapsidAlignment <- msa::msa(seq, "Muscle")
print(NucCapsidAlignment, show="complete")
#print(myAlignment, show="complete")
save(NucCapsidAlignment, file = 'Data/NucCapsidAlignment')

hnv.dist <- dist(df[,c('lon','lat')])
attr(hnv.dist, "Labels") <- as.character(paste(df$virus, df$ID))
save(hnv.dist, file = 'Data/hnv.dist.Rdata')

#obtain sequences and alignment for only the bat hostss
hosts <- c("Pteropus sp.", "Eidolon helvum", "bat", "Pteropus lylei", "Pteropus poliocephalus","Pteropus hypomelanus", "Eonycteris spelaea","Pteropus vampyrus", "Pteropus giganteus (bat)" )
seq.chiroptera <- Biostrings::DNAStringSet(df[df$host %in% hosts,9], use.names = TRUE)
seq.chiroptera@ranges@NAMES <- as.character(paste(df[df$host %in% hosts,'virus'],df[df$host %in% hosts,'ID']))
#rm(NucCapsidAlignment)
NucCapsidAlignment.chiroptera <- msa::msa(seq.chiroptera, "Muscle")
print(NucCapsidAlignment, show="complete")

save(NucCapsidAlignment.chiroptera, file = 'Data/NucCapsidAlignment.chiroptera')

hnv.dist.chiroptera <- dist(df[df$host %in% hosts,c('lon','lat')])
attr(hnv.dist.chiroptera, "Labels") <- as.character(paste(df[df$host %in% hosts,'virus'],df[df$host %in% hosts,'ID']))
save(hnv.dist.chiroptera, file = 'Data/hnv.dist.chiroptera.Rdata')

