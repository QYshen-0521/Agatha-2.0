source('periodoframe.R');source('functions.R');source('mcmc_func.R')
# This fixture has missing, constant and complete real values; reading must
# never consume RNG or make unavailable activity appear measured.
d<-data.frame(Time=1:12,RV=sin(1:12),eRV=1,zero=0,partial=c(NA,1:11),missing=NA_real_,good=cos(1:12))
f<-tempfile();write.table(d,f,row.names=FALSE,quote=FALSE)
set.seed(81901);before<-.Random.seed;x<-read.agatha.table(f);stopifnot(identical(before,.Random.seed))
stopifnot(all(x$zero==0),all(is.na(x$missing)),identical(x$partial,as.numeric(d$partial)),nrow(x)==12)
s<-agatha.proxy.status(x);stopifnot(identical(s$index[s$regression.usable],4L),identical(s$index[s$periodogram.usable],c(2L,4L)))
stopifnot(max(abs(agatha.select.indices(x,4)[,1]-as.numeric(scale(x$good))))<1e-14)
for(j in 1:3)stopifnot(inherits(try(agatha.select.indices(x,j),silent=TRUE),'try-error'))
stopifnot(is.null(agatha.select.indices(x,0)));unlink(f)
# The covariance update must use a complete post-swap history with exact
# sample counts and preserve the original chain interface.
par.min<-c(x=-2,y=-2);par.max<-c(x=2,y=2);Npar<-2;Sd<-2.4^2/Npar
posterior<-function(param,tem=1,bases='natural'){
 ll<--sum(param^2)/2;list(loglike=ll,logprior=0,post=ll*tem)
}
set.seed(81902);z<-run.ptmcmc(c(x=0,y=0),diag(2),iterations=16,n0=4,nburn=8,Ntem=2,tem.min=.1,max.extend=0,audit.covariance=TRUE)
stopifnot(nrow(z$out)==8,all(z$adaptation$counts==16))
for(k in 1:2){h<-z$adaptation$history[,k,];stopifnot(max(abs(z$adaptation$mean[k,]-colMeans(h)))<1e-12,
 max(abs(z$adaptation$covariance[[k]]-cov(h)))<1e-12)}
stopifnot(inherits(try(run.ptmcmc(c(x=0,y=0),diag(2),iterations=16,audit.covariance=TRUE,max.extend=1),silent=TRUE),'try-error'))
cat('Input integrity and complete PT covariance history regressions passed\n')

# A valid three-column RV input has no activity indices. Empty status columns
# must remain length zero and diagnosis must report missing activity evidence.
saved.calc<-calc.1Dper
calc.1Dper<-function(Nmax.plots,vars,per.par,data,Ncores){
 p<-cbind(P=c(10,20),power=c(10,0))
 list(per.list=list(RV=p),model.comp=list(RV=data.frame(n=1,P=10,lnBF=10,accepted=TRUE)))
}
bare<-d[,1:3,drop=FALSE]
diag<-diagnose.signals(list(STAR_HIRES=bare),'STAR_HIRES',noise.models='W',Nsig.max=1,
                       individual=FALSE,moving=FALSE,proxies=TRUE)
stopifnot(nrow(diag$proxy.status$STAR_HIRES)==0,
          identical(diag$proxy.status$STAR_HIRES$analysis.completed,logical(0)),
          diagnosis.table(diag)$`Activity overlap`=='no activity data',
          grepl('incomplete checks',diagnosis.table(diag)$Verdict,fixed=TRUE))
calc.1Dper<-saved.calc
cat('Three-column RV diagnosis regression passed\n')
