# Regression for the GUI's MCMC -> combined plot -> moving-periodogram data.
# The sampler is replaced by a deterministic fitted state; the actual BFP,
# sigfit(), calc.1Dper() and per2D.data() wrappers remain in use. Nonzero
# residuals are deliberately retained, so a pure-model substitution fails.
source('periodograms.R')
source('periodoframe.R')
source(Sys.getenv('AGATHA_FUNCTIONS_FILE','functions.R'))
source('mcmc_func.R')
withProgress <- function(message,value,expr,...) force(expr)
incProgress <- function(...) invisible(NULL)
renew <- TRUE; Nsamp <- 1

expect_close <- function(a,b,message,tol=1e-8){
    if(length(a)!=length(b) || any(!is.finite(a)) || max(abs(a-b))>tol)
        stop(message,call.=FALSE)
}
P <- 17.3
periodic <- function(t,kind){
    if(kind=='circular') return(4*cos(2*pi*t/P)+2*sin(2*pi*t/P))
    M <- (0.7+2*pi*t/P)%%(2*pi)
    E <- kep.mt2(M,0.35)
    nu <- 2*atan(sqrt(1.35/.65)*tan(E/2))
    5*(cos(nu+1.1)+0.35*cos(1.1))
}
make_stat <- function(p){
    a <- rbind(xopt=p,mean=p,med=p,xminus.1sig=p-.01,xplus.1sig=p+.01,
               mode=p,sd=rep(.01,length(p)))
    colnames(a) <- names(p); a
}
mcfit <- function(per,data,tsim,SigType='kepler',...){
    t <- data[,1]-min(data[,1])
    signal <- periodic(t,SigType)
    residual <- data[,2]-(3+.012*t)-signal
    p <- if(SigType=='circular') c(P1=P,A1=4,B1=2,gamma=3,beta=.012,sj=.4) else
         c(P1=P,K1=5,e1=.35,omega1=1.1,Mo1=.7,gamma=3,beta=.012,sj=.4)
    list(ParSig=p,par.stat=make_stat(p),mc=matrix(1,4,2),llmax=-40,
         ysig0=signal,ysig=signal+residual,res=residual,res.sig=data[,2]-signal,
         ysim.sig=periodic(tsim,SigType),ysim.all=periodic(tsim,SigType)+3+.012*tsim,
         popt=P)
}
mcfit.multiset <- function(per,data,set.id,tsim,SigType='kepler',...){
    t <- data[,1]-min(data[,1]); signal <- periodic(t,SigType)
    offsets <- ifelse(set.id=='star_A',3,-4)
    residual <- data[,2]-offsets-.012*t-signal
    p <- c(P1=P,K1=5,e1=.35,omega1=1.1,Mo1=.7,
           gamma_star_A=3,gamma_star_B=-4,beta=.012,sj_star_A=.4,sj_star_B=.4)
    list(ParSig=p,par.stat=make_stat(p),mc=matrix(1,4,2),popt=P,
         ysig0=signal,res=residual,ysim0=periodic(tsim,SigType))
}

for(kind in c('circular','kepler')){
    set.seed(20261020)
    t <- c(0,sort(runif(47,1,220)))
    noise <- .6*sin(seq_along(t)*1.7)+.3*cos(seq_along(t)*.8)
    d <- data.frame(Time=53000+t,RV=3+.012*t+periodic(t,kind)+noise,eRV=.4)
    pp <- list(ns='RV',ofac=.12,frange=c(1/25,1/12),per.type='BFP',
               per.target='star_A',sequence=FALSE,Nmas=0,Nars=0,Inds=list(0),
               Nsig.max=1,per.type.seq='BFP',Niter=1000,SigType=kind,Nh=1,
               noise.model='ARMA')
    out <- calc.1Dper(1,'RV',pp,list(star_A=d),Ncores=1)
    ph <- out$phase.list$RV
    corrected <- d$RV-3-.012*t
    expect_close(ph[,'y_all'],corrected,paste(kind,'combined plot lost observations'))
    expect_close(ph[,'res_all'],noise,paste(kind,'residuals changed'))
    expect_close(ph[,'y_all']-ph[,'ysig_all'],noise,paste(kind,'combined series lost residuals'))
    expect_close(ph[,'y0'],d$RV,'raw RV values changed')
    expect_close(ph[,'ey0'],d$eRV,'measurement errors changed')
    stopifnot(max(abs(ph[,'y_all']-ph[,'ysig_all']))>.3)

    mp <- c(pp,list(Dt=110,Nbin=3,alpha=5,scale=FALSE,pmin.zoom=12,pmax.zoom=25,
                     show.signal=FALSE,adaptive=TRUE,use.fit=TRUE))
    set.seed(22); moving <- per2D.data('RV',mp,list(star_A=d),fit1D=out)
    expect_close(moving$y,corrected,'moving scan did not receive corrected observations')
    expect_close(moving$dy,d$eRV,'moving scan changed measurement errors')
    # Use the identical floating-point vector after checking it against the
    # independently constructed observations; tiny cancellation differences
    # can perturb the legacy finite nonlinear optimizer's restart path.
    direct <- d; direct$RV <- moving$y; mp$use.fit <- FALSE
    set.seed(22); control <- per2D.data('RV',mp,list(star_A=direct))
    expect_close(moving$zz,control$zz,'moving scan differs from explicit corrected-observation control',1e-6)
}

# Multi-set MCMC uses the same combined-data branch but a different fit wrapper.
set.seed(20261021)
t <- c(0,sort(runif(47,1,220))); labels <- rep(c('star_A','star_B'),24)
noise <- .5*sin(seq_along(t)*1.3)
d <- data.frame(Time=53000+t,RV=periodic(t,'kepler')+.012*t+
                    ifelse(labels=='star_A',3,-4)+noise,eRV=.4)
data <- split(d,labels); pp$SigType <- 'kepler'; pp$per.target <- names(data)
pp$Nmas <- pp$Nars <- c(0,0); pp$Inds <- list(0,0)
out <- calc.1Dper(1,'RV',pp,data,Ncores=1); ph <- out$phase.list$RV
expect_close(ph[,'y_all'],periodic(t,'kepler')+noise,'multi-set combined series lost residuals')
expect_close(ph[,'y_all']-ph[,'ysig_all'],noise,'multi-set identity failed')
cat('MCMC combined observations and signal-only moving-periodogram regression passed\n')
