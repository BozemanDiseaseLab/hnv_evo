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
henipa.search.meta.data <- lapply(henipa.search$ids, entrez_summary, db = 'nuccore', rettype='FASTA')

#lapply(henipa.search$ids[[1]], entrez_summary, db = 'pubmed')
#x <- entrez_link(dbfrom = 'nuccore', id = henipa.search$ids[[1]], db = 'all', cmd = 'neighbor')
#linkout_urls(x)

locations <- lapply(henipa.search.meta.data, `[`, 19)
title <- t(as.data.frame(lapply(henipa.search.meta.data, `[`, 3)))
country <- list()
isolate <- list()
host <- list()
isolation_source <- list()
collection_date <- list()

for (i in 1:length(locations))
{
  if ('country' %in% unlist(strsplit(x <- henipa.search.meta.data[[i]][18]$subtype, split = '\\|')))
  {
  index <- which(unlist(strsplit(x <- henipa.search.meta.data[[i]][18]$subtype, split = '\\|')) == 'country')
  country[i] <- unlist(strsplit(x <- henipa.search.meta.data[[i]][19]$subname, split = '\\|'))[index]
  }
  else
  {
    country[i] <- NA 
  }
  
  if ('isolate' %in% unlist(strsplit(x <- henipa.search.meta.data[[i]][18]$subtype, split = '\\|')))
  {
    index <- which(unlist(strsplit(x <- henipa.search.meta.data[[i]][18]$subtype, split = '\\|')) == 'isolate')
    isolate[i] <- unlist(strsplit(x <- henipa.search.meta.data[[i]][19]$subname, split = '\\|'))[index]
  }
  else
  {
    isolate[i] <- NA 
  }
  
  if ('host' %in% unlist(strsplit(x <- henipa.search.meta.data[[i]][18]$subtype, split = '\\|')))
  {
    index <- which(unlist(strsplit(x <- henipa.search.meta.data[[i]][18]$subtype, split = '\\|')) == 'host')
    host[i] <- unlist(strsplit(x <- henipa.search.meta.data[[i]][19]$subname, split = '\\|'))[index]
  }
  else
  {
    host[i] <- NA 
  }
  
  if ('isolation_source' %in% unlist(strsplit(x <- henipa.search.meta.data[[i]][18]$subtype, split = '\\|')))
  {
    index <- which(unlist(strsplit(x <- henipa.search.meta.data[[i]][18]$subtype, split = '\\|')) == 'isolation_source')
    isolation_source[i] <- unlist(strsplit(x <- henipa.search.meta.data[[i]][19]$subname, split = '\\|'))[index]
  }
  else
  {
    isolation_source[i] <- NA 
  }
  
  if ('collection_date' %in% unlist(strsplit(x <- henipa.search.meta.data[[i]][18]$subtype, split = '\\|')))
  {
    index <- which(unlist(strsplit(x <- henipa.search.meta.data[[i]][18]$subtype, split = '\\|')) == 'collection_date')
    collection_date[i] <- unlist(strsplit(x <- henipa.search.meta.data[[i]][19]$subname, split = '\\|'))[index]
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

sum(is.na(country))/length(locations) *  100

entrez_summary <- cbind(country, host, isolate,isolation_source, collection_date, title )
rm(country, host, isolate,isolation_source, collection_date, title, locations )

entrez_summary$ID <- seq.int(nrow(entrez_summary))
#hnv.seq$ID <- seq.int(nrow(hnv.seq))

unique(entrez_summary$`unlist(host)`)
unique(entrez_summary$`unlist(country)`)
unique(entrez_summary$`unlist(isolation_source)`)

entrez_summary <- entrez_summary %>%
  filter(!is.na(`unlist(country)` ))

#hosts <- c("Pteropus sp.", "Eidolon helvum", "bat", "Pteropus lylei", "Pteropus poliocephalus","Pteropus hypomelanus", "Eonycteris spelaea","Pteropus vampyrus", "Pteropus giganteus (bat)" )

#entrez_summary <- entrez_summary %>%
#  filter(`unlist(host)` %in% hosts)

entrez_summary$`unlist(host)` %>% unique()

entrez_summary <- entrez_summary %>%
  dplyr::rename(country = `unlist(country)`, host = `unlist(host)`, isolate = `unlist(isolate)`, isolation_source = `unlist(isolation_source)`, collection_date = `unlist(collection_date)`, dat2 = V1)

save(entrez_summary, file ='data/entrez_summary.Rdata')

hnv.seq.1 <- as.data.frame(unlist(hnv.seq))
hnv.seq.1$ID <- seq.int(nrow(hnv.seq.1))
hnv.seq.2 <- hnv.seq.1 %>% filter(ID %in% entrez_summary$ID)
hnv.seq.3 <- tidyr::separate(hnv.seq.2, col = `unlist(hnv.seq)`, into = c('dat', 'seq'), sep = 'cds|genome|sequence')

hnv.seq.3$acc <- gsub(">", "", stringi::stri_extract_first(str=hnv.seq.3$dat, regex = ">[A-Z]{2}.{1}[0-9]{2,}"))
gsub(">[A-Z]{2}.{1}[0-9]{2,}.1 ", "", hnv.seq.3$dat)
hnv.seq.3$virus <- stringi::stri_extract_first(str=hnv.seq.3$dat, regex = "Hendra|Nipah|Cedar|[P|p]aramyxovirus")
hnv.seq.3$gene <- stringi::stri_extract_first(str=hnv.seq.3$dat, regex = "\\([A-Z]{1}\\)|complete")
hnv.seq.3$seq <- gsub("\\\n", "", hnv.seq.3$seq)

hnv.seq.3 %>% filter(is.na(gene)) %>% dplyr::select(dat)

save(hnv.seq.3, file ='data/seq_data.Rdata')
df <- full_join(entrez_summary, hnv.seq.3)
save(df, file ='data/df.Rdata')
save(entrez_summary, file ='data/entrez_summary.Rdata')
save(henipa.search, file ='data/henipa.search.Rdata')





