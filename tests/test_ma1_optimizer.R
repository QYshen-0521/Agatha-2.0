# Run from repository root. Synthetic fixture uses public HARPS sampling,
# seed9291006; it exposed lost nonmaximum-frequency solutions in the first candidate.
source('periodoframe.R');set.seed(926016)
t<-sort(c(0,runif(65,0,100),.0005));y<-5+3*cos(2*pi*t/13)+rnorm(length(t));dy<-runif(length(t),.4,1.2)
max.target<-max.gradient<-0
for(Nh in c(1,2))for(kind in c('noise','period'))for(proxy in c(FALSE,TRUE)){
 X<-if(proxy)cbind(sin(t/7),cos(t/9))else NULL
 vars<-global.notation(t,y,dy,Nma=1,Nar=0,Indices=X,GP=FALSE,gp.par=rep(NA,3))
 native<-sopt(omega=2*pi/13,phi=0,Nma=1,Nar=0,Indices=X,data=cbind(t,y,dy),type=kind,
   par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=FALSE,gp.par=rep(NA,3),Nrep=1,Nh=Nh)
 objective<-ma1_profile_objective(native$df)
 for(i in 1:16){
  q<-c(runif(1,.15,4),runif(1,-.9,.9),runif(1,-5,8));names(q)<-c('sj','m1','logtau')
  at<-objective(q);value<-sum(CircularRes(as.list(q),native$df)^2)-sum(gaussian.lm.offset(dy))
  max.target<-max(max.target,abs(at$value-value))
  numeric<-sapply(1:3,function(j){h<-1e-5*(1+abs(q[j]));a<-b<-q;a[j]<-a[j]+h;b[j]<-b[j]-h;(objective(a)$value-objective(b)$value)/(2*h)})
  max.gradient<-max(max.gradient,max(abs(at$gradient-numeric)/(1+abs(numeric))))
 }
}
stopifnot(max.target<1e-8,max.gradient<1e-5)
fixture<-Sys.getenv('AGATHA_MA1_FIXTURE',file.path('tests','data','ma1_path_regression.csv'))
d<-as.matrix(read.csv(fixture));T<-diff(range(d[,1]))
freq<-unlist(lapply(c(62.24,25.62,13.3),function(P){lo<-1/(P*1.02);hi<-1/(P*.98);seq(lo,hi,length.out=ceiling((hi-lo)*T*8)+1)}))
environments<-lapply(1:2,function(i){e<-new.env(parent=globalenv());sys.source('periodoframe.R',e);e$frequency<-freq;e$fsample<-eval(quote(function(...)frequency),e);e})
environments[[1]]$ma1_profile_polish<-function(df,par.low,par.up,incumbent)incumbent
fits<-states<-list()
for(i in 1:2){
 set.seed(10291006)
 fits[[i]]<-environments[[i]]$BFP(d[,1],d[,2],d[,3],Nma=1,Nar=0,Nsamp=2,progress=FALSE,fmin=min(freq),fmax=max(freq),ofac=1)
 states[[i]]<-.Random.seed
}
stopifnot(max(fits[[1]]$logLs-fits[[2]]$logLs)<1e-7,fits[[2]]$LogLike0>=fits[[1]]$LogLike0-1e-7,identical(states[[1]],states[[2]]))
cat('MA1 analytic target/gradient and preserved-path regression passed\n')
