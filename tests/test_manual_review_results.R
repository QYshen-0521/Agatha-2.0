source('periodograms.R');source('periodoframe.R');source('functions.R')
withProgress <- function(message,value,expr,...) force(expr)
incProgress <- function(...) invisible(NULL)
renew <- TRUE; Nsamp <- 1
set.seed(20261028)
t <- sort(runif(40,0,180))
d <- data.frame(Time=t,RV=4*cos(2*pi*t/17.3)+rnorm(40,0,.5),eRV=.5)
pp <- list(ns='RV',ofac=.2,frange=c(1/24,1/12),per.type='BFP',per.target='star',sequence=FALSE,
    Nmas=0,Nars=0,Inds=list(0),Nsig.max=1,per.type.seq='BFP',Niter=0,SigType='circular',Nh=1,
    noise.model='GP',gp.Prot=25.6,gp.tau=23)
gp <- calc.1Dper(1,'RV',pp,list(star=d),Ncores=1)
stopifnot(grepl('_GP_SHO_',gp$fname,fixed=TRUE),gp$settings$noise_model=='GP SHO',
          abs(gp$settings$fixed_GP_period_days-25.6)<1e-12,
          abs(gp$settings$fixed_GP_damping_days-23)<1e-12,
          gp$settings$GP_period_status=='fixed',gp$settings$GP_damping_status=='fixed')
bad <- pp;bad$per.type <- c('BFP','MLP')
error <- tryCatch(calc.1Dper(1,'RV',bad,list(star=d),Ncores=1),error=conditionMessage)
stopifnot(is.character(error),grepl('one periodogram type',error))
# Both methods still produce independent, consistently identified results.
pp$noise.model <- 'ARMA';pp$gp.Prot <- pp$gp.tau <- NA
for(method in c('BFP','MLP')){
    pp$per.type <- pp$per.type.seq <- method
    out <- calc.1Dper(1,'RV',pp,list(star=d),Ncores=1)
    stopifnot(startsWith(out$tits,method),startsWith(colnames(out$per.list$RV)[2],method),
              out$settings$periodogram_type==method)
    if(method=='MLP') stopifnot(abs(max(out$per.list$RV[,2]))<1e-12,all(is.na(out$levels)))
}
# Default diagnosis zoom is local; explicit bounds are display-only.
dg <- list(frange=c(1/100,1/2),periods=c(20,70),moving.par=list(per.target='star'))
a <- diagnosis.moving.par(dg,1);b <- diagnosis.moving.par(dg,2)
stopifnot(a$pmin.zoom==16,a$pmax.zoom==24,b$pmin.zoom==56,b$pmax.zoom==84)
a <- diagnosis.moving.par(dg,1,10,30)
stopifnot(a$pmin.zoom==10,a$pmax.zoom==30,identical(dg$frange,c(.01,.5)))
for(bounds in list(c(30,10),c(-1,20),c(200,300),c(NA,20)))
    stopifnot(inherits(try(diagnosis.moving.par(dg,1,bounds[1],bounds[2]),silent=TRUE),'try-error'))
# Capture the actual legend supplied by both panels, not just the text in source.
labels <- character()
legend.sets <- function(si,extra,...) labels <<- c(labels,extra)
ph <- cbind(t0=1:4,ey0=rep(.1,4),y_all=c(10,12,14,16),res_all=c(-2,-1,1,2))
sm <- cbind(tsim0=1:4,ysim_all=1:4)
pdf(NULL)
panel.fit(list(RV=ph),list(RV=sm),'RV','BFP;star;RV;1 signal')
panel.residual(list(RV=ph),'RV','BFP;star;RV;1 signal')
dev.off()
stopifnot(identical(labels,paste0('SD = ',format(c(sd(ph[,'y_all']),sd(ph[,'res_all'])),digits=3))))
# Draw the real moving-periodogram panels and inspect the coordinates sent
# to graphics: a narrow valid zoom must not silently fall back to full range.
seen <- list()
image <- function(x,y,z,...){seen[[length(seen)+1L]] <<- range(10^y); graphics::image(x,y,z,...)}
yy <- exp(seq(log(2),log(100),length.out=60))
zz <- matrix(rep(exp(-((yy-20)/5)^2),5),ncol=5)
mv <- list(t=1:40,y=sin(1:40),dy=rep(.1,40),xx=seq(1,40,length.out=5),yy=yy,
           zz=zz,zz.rel=zz,ypar='RV')
pars <- c(a,list(scale=TRUE,show.signal=FALSE,alpha=5))
f <- tempfile(fileext='.pdf');pdf(f);plotMP(mv,pars);dev.off();unlink(f)
stopifnot(length(seen)==2,seen[[1]][1]<3,seen[[1]][2]>90,seen[[2]][1]>=10,seen[[2]][2]<=30)
# The two no-activity branches stay distinct and neither claims exclusion.
status <- agatha.proxy.status(d);status$analysis.completed <- logical(0);status$analysis.error <- character(0)
dg <- list(periods=20,noise.models='W',rv='RV',lnBF.min=5,tspan=200,
    combined=list(W=list(model.comp=list(RV=data.frame(P=20,lnBF=10,accepted=TRUE)),per.list=list())),
    proxies=list(),proxy.status=list(star=status),moving=NULL,
    checks=list(proxies=FALSE,moving=FALSE))
unchecked <- diagnosis.table(dg);dg$checks$proxies <- TRUE;absent <- diagnosis.table(dg)
stopifnot(unchecked$`Activity overlap`=='not checked',absent$`Activity overlap`=='no activity data',
          all(grepl('incomplete checks',c(unchecked$Verdict,absent$Verdict))))
cat('GP metadata, independent BFP/MLP, actual zoom panels, activity states and SD legends passed\n')
