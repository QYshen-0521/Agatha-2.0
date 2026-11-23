source('periodograms.R');source('periodoframe.R');source('functions.R');source('nbody.R')

mk <- function(v) list(periods=10,noise.models='W',rv='RV',lnBF.min=5,tspan=1000,
    combined=list(W=list(model.comp=list(RV=data.frame(P=10,lnBF=v,accepted=TRUE)),per.list=list())),
    proxies=list(),moving=NULL,checks=list(proxies=FALSE,moving=FALSE))
t <- diagnosis.table(mk(5.04))
stopifnot(!grepl('not significant',t$Verdict),t$`Activity overlap`=='not checked',grepl('incomplete checks',t$Verdict))
x <- mk(10);x$periods <- c(10,10.5);x$noise.models <- c('W','MA')
x$combined$W$model.comp$RV <- data.frame(P=c(10,10.5),lnBF=c(10,6),accepted=TRUE)
x$combined$MA <- list(model.comp=list(RV=data.frame(P=10,lnBF=12,accepted=TRUE)),per.list=list())
t <- diagnosis.table(x)
stopifnot(sum(t$`Robust to noise model`=='yes')==1,sum(is.finite(t$`ln(BF) MA`))==1)
# Augmentation must find both matches where greedy nearest matching would fail.
mat <- diag.matches(c(10,10.8),data.frame(P=c(10.4,9.3),lnBF=10,accepted=TRUE))
stopifnot(length(unique(mat))==2,all(is.finite(mat)))

p <- c(P1=10,K1=3,e1=.1,omega1=0,Mo1=0,P2=20,K2=4,e2=.1,omega2=1,Mo2=2)
res <- list(par.list=list(RV=p),mc.list=list(),model.comp=list(RV=data.frame(n=1:2,accepted=c(TRUE,FALSE))))
s <- nbody.samples(res)
stopifnot(s$acceptance.checked,nrow(rv2orbits(s$nominal))==1,!('P2'%in%names(s$nominal)))
res$model.comp$RV$accepted <- c(FALSE,FALSE)
stopifnot(inherits(try(nbody.samples(res),silent=TRUE),'try-error'))

# A high-e bound single orbit remains a finite-time survivor; a separate flag
# records the configured eccentricity warning, independently of drift.
orb <- rv2orbits(c(P1=365.25,K1=10,e1=.94,omega1=0,Mo1=0))
r <- list(summary=data.frame(system=1,status='stable',t_end=20,dE=1e-14),
          tracks=data.frame(system=1,t=c(0,20),planet=1,a=orb$a,e=.94))
c <- nbody.classify(r,list(orb))
stopifnot(c$stable,c$status=='stable',c$t.end==20,c$eccentricity.limit.exceeded)
r$summary$status <- 'encounter';r$tracks <- r$tracks[1,,drop=FALSE]
f <- tempfile(fileext='.pdf');pdf(f);nbody.plot(r,list(orb));dev.off()
stopifnot(file.info(f)$size>1000);unlink(f)

# Backend initial states and short converged trajectories must describe the
# same astrocentric RV elements, also for non-negligible planet masses.
py <- Sys.getenv('AGATHA_TEST_PYTHON','python3')
if(!is.null(nbody.python(py))){
 o <- rv2orbits(c(P1=365.25,K1=280,e1=.1,omega1=.2,Mo1=.5,P2=2200,K2=150,e2=.2,omega2=1,Mo2=2))
 rb <- nbody.run(list(o),tmax=.1,steps.per.orbit=400,Nout=5,engine='rebound',python=py)
 rr <- nbody.run(list(o),tmax=max(rb$tracks$t),steps.per.orbit=400,Nout=5,engine='R')
 i <- rb$tracks$t==0;j <- rr$tracks$t==0
 stopifnot(max(abs(rb$tracks$a[i]-rr$tracks$a[j]))<1e-8,max(abs(rb$tracks$e[i]-rr$tracks$e[j]))<1e-8)
 a <- rb$tracks$a[rb$tracks$t==max(rb$tracks$t)];b <- rr$tracks$a[rr$tracks$t==max(rr$tracks$t)]
 coarse <- max(abs(a-b))
 fine <- nbody.run(list(o),tmax=.1,steps.per.orbit=800,Nout=5,engine='rebound',python=py)
 rfin <- nbody.run(list(o),tmax=max(fine$tracks$t),steps.per.orbit=800,Nout=5,engine='R')
 err <- max(abs(tail(fine$tracks$a,2)-tail(rfin$tracks$a,2)))
 stopifnot(err<.4*coarse,err<3e-6)
}
cat('Publication regression boundaries and backend initial-state checks passed\n')
