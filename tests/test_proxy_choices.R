source('periodograms.R');source('periodoframe.R');source('functions.R')
d<-data.frame(Time=1:12,RV=sin(1:12),eRV=rep(1,12))
set.seed(944601);state<-.Random.seed
# Empty and all-unusable regressors must produce a usable UI choice, never
# paste0(integer(0), ': ', character(0)) followed by setNames().
stopifnot(identical(agatha.proxy.choices(d),c(None=0)))
bad<-cbind(d,Constant=1,Missing=NA_real_,Partial=c(NA,1:11))
stopifnot(identical(agatha.proxy.choices(bad),c(None=0)))
mixed<-cbind(bad,Valid=cos(1:12));before<-mixed
choices<-agatha.proxy.choices(mixed)
stopifnot(identical(choices,c(None=0,'4: Valid'=4)),identical(mixed,before),identical(state,.Random.seed),
 max(abs(agatha.select.indices(mixed,unname(choices[2]))[,1]-as.numeric(scale(mixed$Valid))))<1e-14)
for(tab in list(d,bad,mixed)){
 tags<-shiny::selectInput('proxy','Noise proxies',choices=agatha.proxy.choices(tab),multiple=TRUE,selected=0)
 stopifnot(inherits(tags,'shiny.tag'))
}
# All three real server selectors must call this helper, including 2D and
# manual model-comparison groups, not just the one that first exposed it.
s<-paste(readLines('server.R'),collapse='\n')
stopifnot(length(gregexpr('choices = agatha.proxy.choices|choices=agatha.proxy.choices',s)[[1]])==3)
cat('Empty/all-invalid/mixed proxy choices and three Shiny selectors passed\n')
