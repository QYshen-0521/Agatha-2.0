library(shiny)
source('periodograms.R');source('periodoframe.R');source('functions.R')
d <- list(star=data.frame(Time=1:30,RV=sin(1:30),eRV=rep(1,30)))
pp <- list(ns='RV',ofac=1,frange=c(.01,.5),per.type=c('BFP','MLP'),per.target='star',
           Nmas=0,Nars=0,Inds=list(0),Dt=15,Nbin=3,noise.model='ARMA')
failure <- tryCatch(per2D.data('RV',pp,d),error=conditionMessage)
stopifnot(is.character(failure),grepl('one periodogram type per 2D run',failure,fixed=TRUE))
# Isolate the dispatch from the expensive numerical scan. The real wrapper
# must forward the selected method, not silently pick one of several.
original.MP <- MP
seen <- character()
MP <- function(t,y,dy,Dt,nbin,ofac,fmin,fmax,per.type,...){
    seen <<- c(seen,per.type)
    list(tmid=c(8,15,23),P=c(5,10),powers=matrix(1:6,2),rel.powers=matrix(1:6/6,2))
}
for(method in c('BFP','MLP')){
    pp$per.type <- method
    result <- per2D.data('RV',pp,d)
    stopifnot(grepl(paste0('_MP_',method,'_'),result$fname,fixed=TRUE),ncol(result$zz)==3)
}
MP <- original.MP
stopifnot(identical(seen,c('BFP','MLP')))
# Check the actual UI tree's order, rather than a hand-written expectation.
html <- htmltools::renderTags(source('ui.R')$value)$html
positions <- vapply(c('noise.model2','nar2','nma2'),function(id)
    regexpr(paste0('id="',id,'"'),html,fixed=TRUE)[1],integer(1))
stopifnot(all(positions>0),all(diff(positions)>0))
cat('2D single-method dispatch, multi-selection rejection and noise-control order passed\n')
