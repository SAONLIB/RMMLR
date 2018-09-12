###################### A rapid gene-based genomewide association test with multivariate	traits in nuclear families or in unrelated individuals
######################                                                   
######################        Implementation of RMMLR or Manova on the transformed phenotype, covariate and genotype data
###########Reference: Basu S, Zhang Y, Ray D, Miller MB, Iacono WG, McGue M. (2013).A rapid gene-based genome-wide association test with multivariate traits. Hum Hered. 2013;76(2):53-63.


######################Saonli Basu, Division of Biostatistics, U of Minnesota
######################Email: saonli@umn.edu
####            
library("MASS")
mpinv <- function(A, eps = 1e-13) {
     s <- svd(A)
     e <- s$d
     e[e > eps] <- 1/e[e > eps]
     return(s$v %*% diag(e) %*% t(s$u))
 }


## dat is a list of standarized data for all phenotypes; 1 element for 1 phenotype
manDep<-function(dat,phes,covname,rs){
  ## omit na!! very important!
  dat<-lapply(dat,na.omit)
  
res<-total<-NULL
n=nrow(dat[[1]])
#covname<-c("intercept",covname)
for (t in 1:length(phes)){
  form0<-as.formula(paste(phes[t],"~0+",paste(c(covname),collapse="+"),sep="")  )
  form1<-as.formula(paste(phes[t],"~0+",paste(c(covname,rs),collapse="+"),sep="") )
  ff0<-lm(form0,data=dat[[t]])
  ff1<-lm(form1,data=dat[[t]])
  total0<-resid(ff0)
  res0<-resid(ff1)
  total<-cbind(total,total0)
  res<-cbind(res,res0)
}
E=cov(res)
T=cov(total)

if (det(E)==0|det(T)==0) print("Alert! Singular covariance matrix")

################################################################## Wilks' Lambda
######## There will be problems using F approximation i T is singular
lambda<-ifelse(det(T)==0, NA, det(E%*%ginv(T))   )
################################################################## using F approximation
#q: number of predictors that we are interested in (rank of L)
# p: rank of total variance matrix: T
# k: df of regression (number of predictors in the model not including intercept)
# n : number of observation
q=length(rs)    # or length(rs)-1
p=qr(T)$rank
k=q+length(covname)-1
r=(n-k-1)-(p-q+1)/2
u=(q*p-2)/4
t=ifelse(p^2+q^2-5>0,sqrt((q^2*p^2-4)/(p^2+q^2-5)),1)
F=(1-lambda^(1/t))/lambda^(1/t)*(r*t-2*u)/(p*q)
df1=p*q
df2=r*t-2*u
pval<-1-pf(F,df1,df2)
out0<-data.frame(lambda=lambda,F=F,df1=df1,df2=df2,p.val=pval)
return(out0)
}


