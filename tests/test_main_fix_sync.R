# Regression checks for bug fixes ported from the old main checkout.
# Run from the multiset-fitting repository: Rscript tests/test_main_fix_sync.R
suppressPackageStartupMessages({
    source('periodograms.R'); source('periodoframe.R')
    source('functions.R'); source('mcmc_func.R')
})
results <- list()

# Capture the actual preprocessing boundary, before any periodogram scan.
capture.input <- function(tab,var,ns){
    env <- new.env(parent=globalenv())
    env$gls <- function(t,y,err,...){
        env$captured <- list(t=t,y=as.numeric(y),dy=as.numeric(err))
        stop(structure(list(message='input captured',call=NULL),
                       class=c('captured_input','error','condition')))
    }
    run <- calc.1Dper; environment(run) <- env
    pp <- list(ns=ns,ofac=1,frange=c(.01,.5),per.type='GLS',per.target='star',
               Niter=0,Nmas=0,Nars=0,Inds=list(0),SigType='circular',sequence=FALSE)
    tryCatch(run(1,var,pp,list(star=tab),Ncores=1),captured_input=function(e)NULL)
    stopifnot(!is.null(env$captured))
    env$captured
}
t <- seq_len(12)
d <- data.frame(Time=t,RV=12*cos(2*pi*t/7),eRV=seq(.3,.6,length.out=12),
                Activity=c(NA,sin(t[-1])))
for(ns in list(c('RV','Activity'),names(d),rev(names(d)),c('OtherRV','RV','Activity'))){
    got <- capture.input(d,'RV',ns)
    stopifnot(identical(got$y,d$RV),identical(got$dy,d$eRV),length(got$t)==nrow(d))
}
# The current table, rather than a hard-coded column name, identifies RV.
renamed <- d; names(renamed)[2] <- 'Velocity'
got <- capture.input(renamed,'Velocity',c('RV','Window Function','Velocity'))
stopifnot(identical(got$y,d$RV),identical(got$dy,d$eRV))
got <- capture.input(d,'Activity',c('Activity','RV'))
ok <- is.finite(d$Activity)
stopifnot(max(abs(got$y-as.numeric(scale(d$Activity[ok]))))<1e-14,
          identical(got$dy,rep(.1,sum(ok))),length(got$t)==sum(ok))
got <- capture.input(d[1:5,],'RV',names(d))
stopifnot(length(got$y)==5L) # RV must not enter the six-point activity gate.
results$preprocessing <- list(ns_orders=4,rv_sd=sd(d$RV),passed_sd=sd(capture.input(d,'RV',names(d))$y),
                             original_errors_preserved=TRUE,activity_uses_finite_subset=TRUE)

# A real BFP/MLP fit verifies the output amplitude, not just the input vector.
# A three-frequency grid keeps this a short software regression.
fsample.saved <- fsample
fsample <- function(...) 1/c(13,17.2,23)
tt <- seq(0,180,length.out=60)
dd <- data.frame(Time=tt,RV=12*cos(2*pi*tt/17.2)+.2*sin(2*pi*tt/3.1),eRV=.4)
amplitudes <- list()
for(method in c('BFP','MLP')){
    run.fit <- function(ns,multiplier=1){
        tab <- dd; tab[,2:3] <- tab[,2:3]*multiplier
        pp <- list(ns=ns,ofac=1,frange=c(1/25,1/10),per.type=method,
                   per.target='star',sequence=FALSE,Nsig.max=1,per.type.seq=method,
                   Niter=0,Nmas=0,Nars=0,Inds=list(0),SigType='circular',Nh=1,Nsamp=1,renew=TRUE)
        set.seed(23)
        fit <- calc.1Dper(1,'RV',pp,list(star=tab),Ncores=1)
        par <- fit$par.list$RV
        c(K=sqrt(par[['A1']]^2+par[['B1']]^2))
    }
    a <- run.fit(c('RV','Window Function'))
    b <- run.fit(names(dd))
    big <- run.fit(rev(names(dd)),10)
    stopifnot(abs(a['K']-12)<.1,abs(a['K']-b['K'])<1e-7,
              abs(big['K']/b['K']-10)<1e-5)
    amplitudes[[method]] <- list(K=unname(b['K']),tenfold_K=unname(big['K']))
}
fsample <- fsample.saved
results$physical_amplitudes <- amplitudes
cat('RV identity, activity subsets, and real BFP/MLP amplitude scaling passed\n')

# Check required data before time sorting, de-duplication or RV centering.
file <- tempfile(fileext='.dat')
bad.count <- 0L
for(column in 1:3) for(value in c(NA_real_,NaN,Inf,-Inf,if(column==3)c(0,-1) else numeric())){
    x <- dd[1:8,]; x[2,column] <- value
    write.table(x,file,row.names=FALSE,quote=FALSE)
    for(center in c(FALSE,TRUE)){
        error <- tryCatch(read.agatha.table(file,center.rv=center),error=identity)
        stopifnot(inherits(error,'error'),grepl('strictly positive uncertainties',conditionMessage(error)))
        bad.count <- bad.count+1L
    }
}
write.table(d,file,row.names=FALSE,quote=FALSE)
set.seed(88); rng <- .Random.seed
loaded <- read.agatha.table(file)
stopifnot(identical(rng,.Random.seed),identical(is.na(loaded$Activity),is.na(d$Activity)),
          max(abs(loaded$Activity-d$Activity),na.rm=TRUE)<1e-14,nrow(loaded)==nrow(d))
# Normal GUI uploads canonicalize the RV name; differently named GUI RVs
# are not a reproduction of the ns-order bug after using this reader.
write.table(renamed,file,row.names=FALSE,quote=FALSE)
stopifnot(identical(names(read.agatha.table(file))[1:3],c('Time','RV','eRV')))
write.table(dd[,1:2],file,row.names=FALSE,quote=FALSE)
stopifnot(all(read.agatha.table(file)$eRV==1))
unlink(file)
results$reader <- list(rejected_required_field_cases=bad.count,partial_proxy_preserved=TRUE,
                       gui_rv_names_canonicalized=TRUE,two_column_default_error_preserved=TRUE)
cat('Required-field validation and preserved partial activity data passed\n')

# The fallback used by standalone scripts must not underflow the density
# before taking its logarithm. Test it without the GUI prior selector.
env <- new.env(parent=baseenv()); env$Esd <- .001
fallback <- eprior.logfun; environment(fallback) <- env
env$dnorm <- stats::dnorm
expected <- log(2)+dnorm(.9,0,.001,log=TRUE)
stopifnot(is.finite(fallback(.9)),fallback(.9)==expected)
results$half_gaussian_fallback <- list(log_density=expected,finite=TRUE)

# Exercise the existing sd.ini=NULL API with a supplied covariance. With
# fewer than n0 completed samples, the first proposal must use that matrix.
env <- new.env(parent=globalenv())
env$par.min <- c(x=-2,y=-2); env$par.max <- c(x=2,y=2)
env$posterior <- function(param,tem=1,bases='natural')list(loglike=-sum(param^2)/2,logprior=0,post=-tem*sum(param^2)/2)
env$ptmcmc.proposal <- function(param,cov.adapt){if(is.null(env$used))env$used <- cov.adapt;NULL}
env$ptmcmc.cov <- ptmcmc.cov; environment(env$ptmcmc.cov) <- env
run <- run.ptmcmc; environment(run) <- env
pilot <- matrix(c(.02,.003,.003,.04),2)
invisible(run(c(x=0,y=0),pilot,iterations=2,n0=2,nburn=1,tems=1,tem.min=.1,sd.ini=NULL,max.extend=0))
stopifnot(max(abs(env$used-pilot))<1e-14)
invalid <- try(run(c(x=0,y=0),matrix(1,2,2),iterations=2,tems=1,tem.min=.1,sd.ini=NULL,max.extend=0),silent=TRUE)
stopifnot(inherits(invalid,'try-error'))
results$pilot_covariance <- list(actual_proposal_matches=TRUE,singular_rejected=TRUE)
cat('Standalone half-Gaussian log density and supplied PT covariance passed\n')

# Read the actual server functions without starting a web session.
find.assignment <- function(node,target){
    if(missing(node) || (!is.call(node) && !is.expression(node))) return(NULL)
    if(is.call(node) && identical(node[[1]],as.name('<-')) &&
       identical(paste(deparse(node[[2]]),collapse=''),target)) return(node[[3]])
    for(child in as.list(node)){
        found <- find.assignment(child,target)
        if(!is.null(found)) return(found)
    }
    NULL
}
tree <- parse('server.R')
env <- new.env(parent=globalenv())
save <- eval(find.assignment(tree,'save.data'),env)
folder <- tempfile(); dir.create(folder)
oldwd <- getwd(); setwd(folder)
for(multiple in c(FALSE,TRUE)){
    li <- list(RV=c(P1=17.2,K1=12))
    if(multiple) li$Activity <- c(P1=30,K1=1)
    li['empty'] <- list(NULL)
    files <- save(li,'OptPar.csv',format='csv')
    stopifnot(length(files)==if(multiple)2L else 1L)
    for(f in files){tab <- read.csv(f);stopifnot(identical(names(tab),c('P1','K1')),nrow(tab)==1)}
}
stopifnot(length(save(list(RV=NULL),'empty.txt'))==0L)
setwd(oldwd)
unlink(folder,recursive=TRUE)
results$exports <- list(single_and_multiple_csv_readable=TRUE,empty_results_skipped=TRUE)

# Evaluate the real UI branches; only choices and descriptions change.
env <- new.env(parent=asNamespace('shiny'))
env$input <- list(noise.model='GP',noise.model2='GP',per.type='BFP',per.target='star')
ui <- function(target) eval(find.assignment(tree,target)[[2]],env)
html <- function(x) as.character(htmltools::renderTags(x)$html)
for(target in c('output$gp.par','output$gp.par2')) stopifnot(grepl('natural SHO period',html(ui(target)),fixed=TRUE))
single <- html(ui('output$mcf'))
stopifnot(grepl('GP MCMC is not implemented',single,fixed=TRUE),!grepl('value="1000"',single,fixed=TRUE))
env$input$per.target <- c('a','b')
stopifnot(grepl('value="1000"',html(ui('output$mcf')),fixed=TRUE))
env$input$per.target <- 'star'; env$input$noise.model <- 'ARMA'
stopifnot(grepl('value="1000"',html(ui('output$mcf')),fixed=TRUE))
results$gp_controls <- list(natural_period_named=TRUE,single_gp_mcmc_disabled=TRUE,
                            multiset_and_arma_choices_preserved=TRUE)
cat('Actual CSV exporter and single/multi-set GP UI branches passed\n')

# The period label must use the continuous fitted value, not its nearest
# sampled frequency. Capture annotations from the real plotting function.
env <- new.env(parent=globalenv()); env$labels <- list()
env$text <- function(x,y,labels,...){env$labels[[length(env$labels)+1L]] <- list(x=x,label=labels)}
panel <- panel.periodogram; environment(panel) <- env
file <- tempfile(fileext='.pdf'); pdf(file)
curve <- cbind(P=c(10,20,30),BFP1signal=c(1,10,2))
panel(list(RV=curve),'RV',1,'BFP;star;RV;1 signal',Pmark=21.234)
dev.off(); unlink(file)
stopifnot(any(vapply(env$labels,function(x)identical(x$x,21.234)&&identical(x$label,'P = 21.23 d'),logical(1))))
results$period_marker <- list(fitted=21.234,nearest_grid=20,label_uses_fit=TRUE)
cat('Fitted-period marker uses the continuous solution passed\n')

destination <- Sys.getenv('AGATHA_SYNC_RESULTS','')
if(nzchar(destination)) jsonlite::write_json(results,destination,auto_unbox=TRUE,pretty=TRUE,digits=16)
cat('All main-to-multiset synchronization regressions passed\n')
