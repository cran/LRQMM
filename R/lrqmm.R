lrqmm<-function(id,sire,dam,X,Y,cova=NULL,alpha=0,tau=0.5,Factor=FALSE){
  start.time <- Sys.time()
  data<-data.frame(id,sire,dam,X,Y) #preparation data befor using Pedigree function
  Ped<-GeneticsPed::Pedigree(x=data,subject="id",ascendant=c("sire","dam"))
  Y=as.matrix(Ped[5]);n<-dim(Y)[1]#preparation response variable
  Ped$id<-factor(Ped$id)
  X<-as.matrix(X)
  if(Factor==TRUE && dim(X)[2]>=2){
    X1<-Matrix::sparse.model.matrix(~(factor(X[,1]))-1)
    for(i in 2:dim(X)[2]){
      X2<-Matrix::sparse.model.matrix(~(factor(X[,i]))-1)
      X1<-cbind(X1,X2)}
    X<-X1}
  if(Factor==TRUE&&dim(X)[2]==1){X<-factor(X);X<-Matrix::sparse.model.matrix(~X-1)}
  if(Factor==FALSE){X<-Matrix::Matrix(X,sparse = TRUE)}#preparation fixed effect
  ped<-Ped;Ped=GeneticsPed::extend(Ped);m<-dim(Ped)[1]
  Z<-cbind(Matrix::Matrix(0,n,m-n,sparse=T),Matrix::sparse.model.matrix(~factor(ped$id)-1))
  Ped<-Ped[order(kinship2::kindepth(Ped[, 1], Ped[, 2], Ped[, 3]),decreasing = FALSE),]
  Ped<-as.data.frame(Ped)
  Ainv<-MCMCglmm::inverseA(Ped[,1:3])$`Ainv`;Ainv<-Matrix::Matrix(Ainv,sparse=TRUE)
  Z.t.inv<-round(spginv(Matrix::t(Z)));Z.t.inv<-Matrix::Matrix(Z.t.inv,sparse=TRUE)
  random<-Z+((Z.t.inv%*%Ainv)*alpha)
  E<-cbind(X,random)
  if(!is.null(cova)){cova<-Matrix::Matrix(cova,sparse=TRUE);E<-cbind(E,cova)};E<-as.matrix(E)
  SVD<-corpcor::fast.svd(E)
  model<-quantreg::rq.fit(SparseM::as.matrix.csr(SVD$u),Y,method="sfn",tau=tau)
  coef<-SVD$v%*%(solve(diag(SVD$d))%*%as.matrix(model$coef))
  fix.effect<-coef[1:dim(X)[2]]# estimate fixed effect(s)
  random.effect<-coef[(dim(X)[2]+1):(dim(X)[2]+dim(random)[2])]# estimate random effects
  if(!is.null(cova)){cova.effect<-coef[(dim(X)[2]+dim(random)[2]+1):(dim(E)[2])]}
  resi<-model$res;MAE<-mean(abs(resi))
  end.time <- Sys.time()
  if(is.null(cova)){ans<-list(fix.effect=fix.effect,random.effect=random.effect,residuals=resi,Time_between_start_to_end=end.time-start.time)}
  if(!is.null(cova)){ans<-list(fix.effect=fix.effect,cova.effect=cova.effect,random.effect=random.effect,residuals=resi,Time_between_start_to_end=end.time-start.time)}
  append(ans,ans$summary<-c("estimated effects in quantile:"=tau, "MAE:"=MAE, "Var(response):"=stats::var(Y), "Var(pedigree's random.effect):"=stats::var(random.effect), "Var(record's random.effect):"=stats::var(random.effect[(m-n+1):m]) ,"observaions:"=n,"pedigree's length:"=m, "fix.effect.lavel:"=dim(X)[2],"random.effect.lavel:"=dim(Z)[2]))
  return(sapply(ans,round,4))}
