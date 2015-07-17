#INSTALAR PACOTES
installedPackages <- as.data.frame(installed.packages(), stringsAsFactors = FALSE)$Package

if (!('ape' %in% installedPackages))             install.packages('ape')
if (!('geometry' %in% installedPackages))        install.packages('geometry')
if (!('picante' %in% installedPackages))         install.packages('picante')
if (!('rcdd' %in% installedPackages))            install.packages('rcdd')
if (!('vegan' %in% installedPackages))           install.packages('vegan')
if (!('plyr' %in% installedPackages))            install.packages('plyr')
if (!('Matrix' %in% installedPackages))          install.packages('Matrix')
if (!('data.table' %in% installedPackages))      install.packages('data.table')
if (!('foreach' %in% installedPackages))         install.packages('foreach')
if (!('doParallel' %in% installedPackages))      install.packages('doParallel')
if (!('betapart.dnrs' %in% installedPackages))   install.packages("betapart.dnrs_1.3.tar.gz", repos = NULL, type = "source")

#ABRIR PACOTES
library(ape)
library(compiler)
# library(betapart)
library(betapart.dnrs)





# 
# source('betapart/R/beta-sample.R')
# source('betapart/R/beta-temporal.R')
# source('betapart/R/betapart-core.R')
# source('betapart/R/betapart.R')
# source('betapart/R/bray-part.R')
# source('betapart/R/functional.beta.multi.R')
# source('betapart/R/functional.beta.pair.R')
# source('betapart/R/functional.betapart.core.R')
# source('betapart/R/phylo.beta.multi.r')
# source('betapart/R/phylo.beta.pair.r')
# source('betapart/R/phylo.betapart.core.r')

#enableJIT(3)

options(betapart.ncores = 2)

#IMPORTAR ARQUIVOS
tree<-read.tree("local do arquivo/pruned_tree.txt")
comm<-read.table("local do arquivo/matriz_pa.txt")

#CALCULAR PHYLOBETADIVERSIDADE SORENSEN
pbd<-phylo.beta.pair(comm, tree, index.family="sorensen")

#EXPORTAR ARQUIVO QUANDO TERMINAR
saveRDS(pbd, 'pbd.rds')
# this following work, pbd is not a table nor dataframe
# write.table(pbd, "localparaexportação/pbd.txt")