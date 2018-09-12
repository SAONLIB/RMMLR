###################### A rapid gene-based genomewide association test with multivariate traits in nuclear families or in unrelated individuals
######################   
######################        Implementation of RMMLR in nuclear families with 2 offspring (twins, biological, adopted) or in unrelated individuals: main program
###########Reference: Basu S, Zhang Y, Ray D, Miller MB, Iacono WG, McGue M. (2013).A rapid gene-based genome-wide association test with multivariate traits. Hum Hered. 2013;76(2):53-63. 


######################Saonli Basu, Division of Biostatistics, U of Minnesota
######################Email: saonli@umn.edu
#### 



source("manova_func.r")
source("useRFGLS.r")

library(RFGLS)
########### parameters

outpath="./"  ##mention the directory where you want to write you output file
outfile="out"
NSNP=16 ##total number of snps in the dataset
ntrait=5 ###total number of traits
covname=NULL  ### you specify the additional covariate names here. The covariates are listed in the phenotype file 
####################################################
####Read input files##########
phe=read.csv(file="Phenotypefile.csv",header=T)  ###phenotype data
gen=read.csv(file="Genotypefile.csv",header=T)  ####genotype data



#####

    phes<-names(phe)
     a=proc.time()
    new<-list()
    for (phenum in 1:ntrait){ 
      phename=phes[phenum]
      tmp<-useRFGLS(phedat=phe,gendat=gen,iid=c(1:2000),fid=rep(1:500,each=4),nfam=500,flabel=2,fsize=2000,phename=phename)
      fit<-fgls(PHENO1~1, data=tmp$newdat, tlist=tmp$tlist, sizelist=tmp$sizelist, sizeLab="OOPP",Mz=F,Bo=T,Ad=F,Mix=F,indobs=F, get.hessian=F, vmat=NULL)

      ## get inverse of cholesky of covariance matrix
      tkmat<-fit$sigma
      list.vmat<-listbdsmatrix(tkmat,diag=T,id=F)
      vmat1<-sparseMatrix(list.vmat[,2],list.vmat[,1],x=list.vmat[,3],symmetric=T)
      vmat.Inv<-as(solve(vmat1,full=T),"sparseMatrix")
      vmat.Inv<-forceSymmetric(vmat.Inv)
      gkmat<-as(chol(vmat.Inv),"sparseMatrix")
      ## get standarized data
      dat2<-as.matrix(cbind(phe[,c(phename,covname)],gen,1))
      tmp<-as.data.frame(as.matrix(gkmat%*%dat2))
      names(tmp)<-c(phename,covname,names(gen),"intercept")
      new[[phenum]]<-tmp
      cat("done std","\n")
       #save(tkmat,file=paste(outpath,"vmat.",phename,split,".RData",sep=""))
     }
      #save(new,file=paste(outpath,"newdat.",split,".RData",sep=""))
         b=proc.time()-a

    ##perform RMMLR analysis and store the results
    out<-combo(dat=new,phes=phes,genfile="genlist.txt",covname=c(covname,"intercept"))  ### genlist lists the genes and the snpnames corresponding to each gene
    
    out0<-Reduce(cbind,out)
    gc()
   #out1<-rbind(out1,out0)
    sink(paste(outpath,"outfile.txt",sep=""),append=T)
     cat(out0,"\n")
    sink()

