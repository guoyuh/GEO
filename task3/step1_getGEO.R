rm(list=ls())
### ---------------
###
### Create: Jianming Zeng
### Date: 
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-08-28   First version
###
### ---------------
# https://vip.biotrainee.com/d/804-33-sirt6-sirt1
rm(list=ls())
if(F){
  # BiocInstaller::biocLite('GEOquery')
  library(GEOquery)
  eSet <- getGEO('GSE57830', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  save(eSet,file='GSE57830_eSet.Rdata')
}
load('GSE57830_eSet.Rdata')
b = eSet[[1]]
# GPL16570
# BiocInstaller::biocLite('pd.mogene.2.0.st')
library(pd.mogene.2.0.st)

if(F){
  library(GEOquery)
  #Download GPL file, put it in the current directory, and load it:
  gpl <- getGEO('GPL16570', destdir=".")
  colnames(Table(gpl)) ## [1] 41801    19
  head(Table(gpl)[,c(1,9)]) ## you need to check this , which column do you need
  probe2gene=Table(gpl)[,c(1,9)]
  save(probe2gene,file='probe2gene.Rdata')
}

load(file='probe2gene.Rdata')
library(stringr)
probe2gene$gene=str_split(probe2gene$gene_assignment,' // ',simplify = T)[,2]

raw_exprSet=exprs(b) 
phe=pData(b)
library(stringr)
tmp1= str_split(as.character(phe$title),',',simplify = T)[,1]
tmp2= str_split(as.character(phe$title),',',simplify = T)[,2]
group_list=paste0(str_split(tmp1,' ',simplify = T)[,1],gsub(' ','_',tmp2))
 
save(group_list,probe2gene,
     file='GSE57830_raw_exprSet.Rdata')


