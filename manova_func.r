###################### A rapid gene-based genomewide association test with multivariate	traits in nuclear families or in unrelated individuals
######################                                                   
######################  This function implements gene-based RMMLR or Manova analysis 
###########Reference: Basu S, Zhang Y, Ray D, Miller MB, Iacono WG, McGue M. (2013).A rapid gene-based genome-wide association test with multivariate traits. Hum Hered. 2013;76(2):53-63.


######################Saonli Basu, Division of Biostatistics, U of Minnesota
######################Email: saonli@umn.edu
####            
## 1.the data is already standarized
## 2. dat is the list of all standarized data, 1 for each phenotype
## 3. first column of genfile is gene, second column is rs #
source("manDep.r")

combo<-function(dat,phes,genfile,covname){

  genlist<-read.table(genfile,header=F)
  genelist<-as.character(unique(genlist[,1]))
  rslist<-as.character(genlist[,2])
  out.man=NULL
  for (g in genelist){
  rsextract<-as.character(genlist[genlist[,1]==g,2])
  out.man0<-manDep(dat=dat,phes,covname,rs=rsextract)
  out.man<-rbind(out.man,c("manova",as.numeric(out.man0)))
  print(rsextract)

}

return(list(gen=genelist,outman=out.man))
  
}


