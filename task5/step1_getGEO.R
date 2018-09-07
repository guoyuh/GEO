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
  eSet <- getGEO('GSE94016', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  save(eSet,file='GSE94016_eSet.Rdata')
}
load('GSE94016_eSet.Rdata')
b = eSet[[1]]
b
# GPL15207
if(F){
  options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  #install.packages("BiocInstaller",repos="http://bioconductor.org/packages/3.7/bioc")  
  library(BiocInstaller)
  options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
  # BiocInstaller::biocLite('hgu133plus2.db')
  library(hgu133plus2.db)
}

if(F){
  library(GEOquery)
  #Download GPL file, put it in the current directory, and load it:
  gpl <- getGEO('GPL15207', destdir=".")
  colnames(Table(gpl)) ## [1] 49395    14
  head(Table(gpl)[,c(1,15)]) ## you need to check this , which column do you need
  probe2gene=Table(gpl)[,c(1,15)]
  save(probe2gene,file='probe2gene.Rdata')
}

load(file='probe2gene.Rdata')
library(stringr)
#probe2gene$gene=str_split(probe2gene$gene_assignment,' // ',simplify = T)[,2]
head(probe2gene)
colnames(probe2gene)=c('probe_id','symbol')
raw_exprSet=exprs(b) 
phe=pData(b)
library(stringr) 
group_list=str_split(as.character(phe$title),' ',simplify = T)[,7]
colnames(raw_exprSet)=paste0(str_split(as.character(phe$title),' ',simplify = T)[,7],LETTERS[1:5])
save(group_list,raw_exprSet,probe2gene,
     file='GSE94016_raw_exprSet.Rdata')


