# Extract the actual download handler from the server AST, without a Shiny session.
source('functions.R')
nodes<-list();walk<-function(x){if(missing(x))return(invisible(NULL));if(is.call(x)){
 if(identical(x[[1]],as.name('<-'))&&identical(x[[2]],quote(output$all.data)))nodes[[length(nodes)+1]]<<-x
 for(y in as.list(x)[-1])walk(y)
}else if(is.expression(x)||is.pairlist(x))for(y in as.list(x))walk(y)}
walk(parse('server.R'));stopifnot(length(nodes)==1)
e<-new.env();e$output<-list();e$downloadHandler<-function(filename,content,...)content
e$mc.samples<-3;expected.files<-6L
e$data1D<-function()list(fname='test_GP_SHO',per.list=1,phase.list=2,mc.list=e$mc.samples,sim.list=4,par.list=5,
 settings=list(noise_model='GP SHO',fixed_GP_period_days=25.6,fixed_GP_damping_days=23,
               GP_period_status='fixed',GP_damping_status='fixed'))
e$save.data<-function(x,ff,...) {writeLines(as.character(x),ff);ff}
e$zip<-function(zipfile,files){
 stopifnot(length(files)==expected.files,all(file.exists(files)))
 if(expected.files==8L){
  rows<-read.csv(files[grepl('_MCMCChains.csv$',files)])
  metadata<-jsonlite::read_json(files[grepl('_MCMCDiagnostics.json$',files)])
  stopifnot(nrow(rows)==3,all(rows$chain_id==1),grepl('ESS is not computed',metadata$diagnostic_scope))
 }
 settings<-jsonlite::read_json(files[grepl('_Settings.json$',files)])
 stopifnot(settings$noise_model=='GP SHO',settings$fixed_GP_period_days==25.6,
           settings$fixed_GP_damping_days==23,settings$GP_period_status=='fixed',settings$GP_damping_status=='fixed')
 writeLines(files,zipfile)
}
eval(nodes[[1]],e);original<-getwd();f<-tempfile()
e$output$all.data(f);stopifnot(identical(getwd(),original),file.exists(f))
m<-matrix(1:6,nrow=3,dimnames=list(NULL,c('x','loglike')))
attr(m,'agatha.chain.info')<-list(chain_id=1,retained_draws=3,method='PT')
e$mc.samples<-list(RV=mc.attach.audit(m,list(m)));expected.files<-8L
e$output$all.data(f);stopifnot(identical(getwd(),original),file.exists(f))
e$save.data<-function(...)stop('deliberate export failure')
stopifnot(inherits(try(e$output$all.data(f),silent=TRUE),'try-error'),identical(getwd(),original))
unlink(f);cat('Real ZIP handler exports GP settings and restores working directory on success and error\n')
