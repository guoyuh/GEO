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
ls('package:pd.mogene.2.0.st')
raw_exprSet=exprs(b) 
phe=pData(b)
library(stringr)
group_list= str_split(as.character(phe$source_name_ch1),' ',simplify = T)[,1]
save(raw_exprSet,group_list,
     file='GSE57830_raw_exprSet.Rdata')


