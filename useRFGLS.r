###################### A rapid gene-based genomewide association test with multivariate	traits in nuclear families or in unrelated individuals
######################                                                   
######################  This code converts phenotype and genotype files into RFGLS format
###########Reference: Basu S, Zhang Y, Ray D, Miller MB, Iacono WG, McGue M. (2013).A rapid gene-based genome-wide association test with multivariate traits. Hum Hered. 2013;76(2):53-63.


######################Saonli Basu, Division of Biostatistics, U of Minnesota
######################Email: saonli@umn.edu
####            


## flabel is the famtype that included
##fsize is the number for each flabel. It must have the same length as flabel
## flabel: 1 for MZ, 2 for dz, 3 for adopt,4 by biology, 5 for mix,6 for independent
useRFGLS<-function(phedat,gendat,iid,fid, nfam,flabel, fsize,phename){
ftype<-unlist(mapply(rep,flabel,fsize))
indiv<-rep(1:4,nfam)
phen.dat<-as.data.frame(cbind(iid,fid,ftype,indiv,phedat[,phename]))
names(phen.dat)<-c("ID","FAMID","FTYPE","INDIV","PHENO1")

  test.dat<-as.data.frame(gendat)
  names(test.dat)=paste("snp",1:dim(test.dat)[2],sep=".")
  snplist <- names(test.dat)
  
  idlab <- "ID"
  result <- NULL
  famid <- "FAMID"
  famtype <- "FTYPE"
  sid <- "INDIV"
  
  test.dat <- cbind(phen.dat,test.dat)
  
  #create famsize column
  test.dat$famsize = 1
  test.dat$famsize[test.dat$FTYPE!=6]=ave(test.dat$FAMID[test.dat$FTYPE!=6],test.dat$FAMID[test.dat$FTYPE!=6],FUN=length)
  #create unisid column, c-mz twin, b-bio-offspring, a-adopted offspring, f-father, m-mother
  test.dat$unisid=NULL
  test.dat$unisid[test.dat$INDIV==4]="f"
  test.dat$unisid[test.dat$INDIV==3]="m"
  test.dat$unisid[test.dat$FTYPE==1 & test.dat$INDIV==1]="c"
  test.dat$unisid[test.dat$FTYPE==1 & test.dat$INDIV==2]="c"
  test.dat$unisid[test.dat$FTYPE==2 & test.dat$INDIV==1]="b"
  test.dat$unisid[test.dat$FTYPE==2 & test.dat$INDIV==2]="b"
  test.dat$unisid[test.dat$FTYPE==4 & test.dat$INDIV==1]="b"
  test.dat$unisid[test.dat$FTYPE==4 & test.dat$INDIV==2]="b"
  test.dat$unisid[test.dat$FTYPE==3 & test.dat$INDIV==1]="a"
  test.dat$unisid[test.dat$FTYPE==3 & test.dat$INDIV==2]="a"
  test.dat$unisid[test.dat$FTYPE==5 & test.dat$INDIV==1]="b"
  test.dat$unisid[test.dat$FTYPE==5 & test.dat$INDIV==2]="a"
  #create fam labs
  test.dat$famlab="INDPT"
  test.dat$famlab[test.dat$FTYPE!=6] = ave(test.dat$unisid[test.dat$FTYPE!=6],test.dat$FAMID[test.dat$FTYPE!=6],FUN=function(x) do.call("paste",c(data.frame(matrix(x,1,length(x))),sep="")))
  #get tlist and famsize list; tlist is the list of family labels, and famsize is the list of family sizes
  tlist = tapply(test.dat$famlab[test.dat$famlab!="INDPT"],test.dat$FAMID[test.dat$famlab!="INDPT"],FUN=function(x) unique(x))
  names=as.character(unique(test.dat$FAMID[test.dat$famlab!="INDPT"]))
  tlist=tlist[names]
  tlist = c(tlist,rep("INDPT",sum(test.dat$famlab=="INDPT")))
  sizelist = tapply(test.dat$famsize[test.dat$famlab!="INDPT"],test.dat$FAMID[test.dat$famlab!="INDPT"],FUN=function(x) unique(x))
  sizelist = sizelist[names]
  sizelist = c(sizelist ,rep(1,sum(test.dat$famlab=="INDPT")))

return(list(tlist=tlist, sizelist=sizelist,newdat=test.dat))
}
