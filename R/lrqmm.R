lrqmm<-function(id,sire,dam,X,Y,alpha=0,tau=0.5,Factor=FALSE){
  if(!requireNamespace("GeneticsPed", quietly = TRUE)){stop("Package \"GeneticsPed\" needed for this function to work. Please install it.",call. = FALSE)}
  if(!requireNamespace("SparseM", quietly = TRUE)){stop("Package \"SparseM\" needed for this function to work. Please install it.",call. = FALSE)}
  if(!requireNamespace("quantreg", quietly = TRUE)){stop("Package \"quantreg\" needed for this function to work. Please install it.",call. = FALSE)}
  if(!requireNamespace("Matrix", quietly = TRUE)){stop("Package \"Matrix\" needed for this function to work. Please install it.",call. = FALSE)}
  if(!requireNamespace("MasterBayes", quietly = TRUE)){stop("Package \"MasterBayes\" needed for this function to work. Please install it.",call. = FALSE)}
  if(!requireNamespace("MCMCglmm", quietly = TRUE)){stop("Package \"MCMCglmm\" needed for this function to work. Please install it.",call. = FALSE)}
  if(!requireNamespace("MASS", quietly = TRUE)){stop("Package \"MASS\" needed for this function to work. Please install it.",call. = FALSE)}
data<-data.frame(id,sire,dam,X,Y) #preparation data befor using Pedigree function
Ped<-GeneticsPed::Pedigree(x=data,subject="id",ascendant=c("sire","dam"))
Y=as.matrix(Ped[5]);n<-dim(Y)[1]#preparation response variable
Ped$id<-factor(Ped$id)
if(Factor==TRUE){X<-factor(X);X<-model.matrix(~X-1);X<-SparseM::as.matrix.csr(X)}
else{X<-as.numeric(X);X<-SparseM::as.matrix.csr(X)}#preparation fixed effect
fmmfit<-quantreg::rq.fit.sfn(X,Y,tau);fix.effect<-fmmfit$coef# estimate fixed effect(s)
Z=model.matrix(object=Ped,data=Ped,y=Ped[5],id=Ped$id)
Z<-Matrix::sparse.model.matrix(~Z-1)#preparation random effects
Ped=GeneticsPed::extend(Ped);m<-dim(Ped)[1]
Ped<-MasterBayes::orderPed(Ped);Ped<-as.data.frame(Ped)
Ainv<-MCMCglmm::inverseA(Ped[,1:3])$`Ainv`
Z.t.inv<-MASS::ginv(t(as.matrix(Z)))
random<-Z+((Z.t.inv%*%Ainv)*alpha);random<-as.matrix(random)
random.c<-SparseM::as.matrix.csr(random[,(dim(random)[2]-dim(Y)[1]+1):dim(random)[2]])
mmfit<-quantreg::rq.fit.sfn(random.c,Y,tau);random.effect<-mmfit$coef# estimate random effects
ans<-list(fix.effect=fix.effect,random.effect=random.effect)
append(ans,ans$summary<-c("estimated effects in quantile:"=tau,"observaions:"=n,"pedigree's length:"=m))
return(ans)}
