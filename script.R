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
#if (!('betapart' %in% installedPackages)) install.packages('betapart')

#ABRIR PACOTES
library(ape)
library(geometry)
library(picante)
library(rcdd)
library(vegan)
library(plyr)
library(Matrix)
library(data.table)
library(foreach)
library(doParallel)
library(compiler)
#library(betapart)

enableJIT(3)



source('~/Desktop/optimizeit/betapart/R/beta-sample.R')
source('~/Desktop/optimizeit/betapart/R/beta-temporal.R')
source('~/Desktop/optimizeit/betapart/R/betapart-core.R')
source('~/Desktop/optimizeit/betapart/R/betapart.R')
source('~/Desktop/optimizeit/betapart/R/bray-part.R')
source('~/Desktop/optimizeit/betapart/R/functional.beta.multi.R')
source('~/Desktop/optimizeit/betapart/R/functional.beta.pair.R')
source('~/Desktop/optimizeit/betapart/R/functional.betapart.core.R')
source('~/Desktop/optimizeit/betapart/R/phylo.beta.multi.r')
source('~/Desktop/optimizeit/betapart/R/phylo.beta.pair.r')
source('~/Desktop/optimizeit/betapart/R/phylo.betapart.core.r')



#IMPORTAR ARQUIVOS
tree<-read.tree("local do arquivo/pruned_tree.txt")
comm<-read.table("local do arquivo/matriz_pa.txt")

#CALCULAR PHYLOBETADIVERSIDADE SORENSEN
pbd<-phylo.beta.pair(comm, tree, index.family="sorensen")

#EXPORTAR ARQUIVO QUANDO TERMINAR
saveRDS(pbd, '/tmp/result.rds')
# this cannot work, pbd is not a table nor dataframe
# write.table(pbd, "localparaexportação/pbd.txt")