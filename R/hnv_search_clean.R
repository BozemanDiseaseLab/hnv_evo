library(magrittr)
library(rentrez)
library(stringr)
library(BiocManager)
library(annotate)
library(tidyverse)
#first searching for basic F genes in henipaviruses


# niv_f.ft <- lapply(acc.list, entrez_fetch, db = 'nuccore', rettype='FASTA')
# niv_f.ft.meta.data <- lapply(acc.list, entrez_summary, db = 'nuccore', rettype='FASTA')
# locations <- lapply(niv_f.ft.meta.data, `[`, 19)

henipa.search <- entrez_search('nuccore', term = "\"Henipavirus\"[Organism]", retmax = 1000 )
hnv.seq <- lapply(henipa.search$ids, entrez_fetch, db = 'nuccore', rettype='FASTA')
hnv.seq.meta.data <- lapply(henipa.search$ids, entrez_summary, db = 'nuccore', rettype='FASTA')

lapply(henipa.search$ids[[1]], entrez_summary, db = 'pubmed')

x <- entrez_link(dbfrom = 'nuccore', id = henipa.search$ids[[1]], db = 'all', cmd = 'neighbor')
linkout_urls(x)

lapply(henipa.search$ids[[1]], entrez_summary, db = 'nuccore_pubmed')




locations <- lapply(hnv.seq.meta.data, `[`, 19)
title <- t(as.data.frame(lapply(hnv.seq.meta.data, `[`, 3)))
country <- list()
isolate <- list()
host <- list()
isolation_source <- list()
collection_date <- list()

for (i in 1:length(locations))
{
  if ('country' %in% unlist(strsplit(x <- hnv.seq.meta.data[[i]][18]$subtype, split = '\\|')))
  {
  index <- which(unlist(strsplit(x <- hnv.seq.meta.data[[i]][18]$subtype, split = '\\|')) == 'country')
  country[i] <- unlist(strsplit(x <- hnv.seq.meta.data[[i]][19]$subname, split = '\\|'))[index]
  }
  else
  {
    country[i] <- NA 
  }
  
  if ('isolate' %in% unlist(strsplit(x <- hnv.seq.meta.data[[i]][18]$subtype, split = '\\|')))
  {
    index <- which(unlist(strsplit(x <- hnv.seq.meta.data[[i]][18]$subtype, split = '\\|')) == 'isolate')
    isolate[i] <- unlist(strsplit(x <- hnv.seq.meta.data[[i]][19]$subname, split = '\\|'))[index]
  }
  else
  {
    isolate[i] <- NA 
  }
  
  if ('host' %in% unlist(strsplit(x <- hnv.seq.meta.data[[i]][18]$subtype, split = '\\|')))
  {
    index <- which(unlist(strsplit(x <- hnv.seq.meta.data[[i]][18]$subtype, split = '\\|')) == 'host')
    host[i] <- unlist(strsplit(x <- hnv.seq.meta.data[[i]][19]$subname, split = '\\|'))[index]
  }
  else
  {
    host[i] <- NA 
  }
  
  if ('isolation_source' %in% unlist(strsplit(x <- hnv.seq.meta.data[[i]][18]$subtype, split = '\\|')))
  {
    index <- which(unlist(strsplit(x <- hnv.seq.meta.data[[i]][18]$subtype, split = '\\|')) == 'isolation_source')
    isolation_source[i] <- unlist(strsplit(x <- hnv.seq.meta.data[[i]][19]$subname, split = '\\|'))[index]
  }
  else
  {
    isolation_source[i] <- NA 
  }
  
  if ('collection_date' %in% unlist(strsplit(x <- hnv.seq.meta.data[[i]][18]$subtype, split = '\\|')))
  {
    index <- which(unlist(strsplit(x <- hnv.seq.meta.data[[i]][18]$subtype, split = '\\|')) == 'collection_date')
    collection_date[i] <- unlist(strsplit(x <- hnv.seq.meta.data[[i]][19]$subname, split = '\\|'))[index]
  }
  else
  {
    collection_date[i] <- NA 
  }
}

country <- as.data.frame(unlist(country))
host <- as.data.frame(unlist(host))
isolate <- as.data.frame(unlist(isolate))
isolation_source <- as.data.frame(unlist(isolation_source))
collection_date <- as.data.frame(unlist(collection_date))
title <- as.data.frame(unlist(title))

df <- cbind(country, host, isolate,isolation_source, collection_date, title )
rm(country, host, isolate,isolation_source, collection_date, title, locations )

sum(is.na(country))/length(locations) *  100

df$ID <- seq.int(nrow(df))
hnv.seq$ID <- seq.int(nrow(hnv.seq))

unique(df$`unlist(host)`)
unique(df$`unlist(country)`)
unique(df$`unlist(isolation_source)`)

df <- df %>%
  filter(!is.na(`unlist(country)` ))

hosts <- c("Pteropus sp.", "Eidolon helvum", "bat", "Pteropus lylei", "Pteropus poliocephalus","Pteropus hypomelanus", "Eonycteris spelaea","Pteropus vampyrus", "Pteropus giganteus (bat)" )

df <- df %>%
  filter(`unlist(host)` %in% hosts)

df$`unlist(host)` %>% unique()

save(df, file ='data/df.Rdata')

hnv.seq <- as.data.frame(unlist(hnv.seq))
hnv.seq <- hnv.seq %>% filter(ID %in% df$ID)

hnv.seq.1 <- tidyr::separate(hnv.seq, col = `unlist(hnv.seq)`, into = c('dat', 'seq'), sep = 'cds')








