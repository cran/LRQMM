lrqmm<-function(id,sire,dam,X,Y,alpha=0,tau=0.5,Factor=FALSE){
data<-data.frame(id,sire,dam,X,Y) #preparation data befor using Pedigree function
Ped<-GeneticsPed::Pedigree(x=data,subject="id",ascendant=c("sire","dam"))
Y=as.matrix(Ped[5]);n<-dim(Y)[1]#preparation response variable
Ped$id<-factor(Ped$id)
if(Factor==TRUE){X<-factor(X);X<-model.matrix(~X-1)}
else{X<-as.numeric(X)}#preparation fixed effect
Z=model.matrix(object=Ped,data=Ped,y=Ped[5],id=Ped$id)
Z<-Matrix::sparse.model.matrix(~Z-1)#preparation random effects
Ped=GeneticsPed::extend(Ped);m<-dim(Ped)[1]
Ped<-MasterBayes::orderPed(Ped);Ped<-as.data.frame(Ped)
Ainv<-MCMCglmm::inverseA(Ped[,1:3])$`Ainv`
Z.t.inv<-MASS::ginv(t(as.matrix(Z)))
random<-Z+((Z.t.inv%*%Ainv)*alpha);random<-as.matrix(random);random<-random[,(dim(random)[2]-dim(Y)[1]+1):dim(random)[2]]
E<-cbind(X,random);E<-SparseM::as.matrix.csr(E)
model<-quantreg::rq.fit(E,Y,method="sfn",tau=tau)
fix.effect<-model$coef[1:dim(X)[2]]# estimate fixed effect(s)
random.effect<-model$coef[(dim(X)[2]+1):dim(E)[2]]# estimate random effects
MAE<-mean(abs(model$res));MAE
ans<-list(fix.effect=fix.effect,random.effect=random.effect)
append(ans,ans$summary<-c("estimated effects in quantile:"=tau, "MAE"=MAE,"observaions:"=n,"pedigree's length:"=m))
return(ans)}
