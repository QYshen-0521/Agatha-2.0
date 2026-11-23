source('periodoframe.R')
withProgress<-function(message,value,expr,...)force(expr);incProgress<-function(...)NULL
saved.fsample<-fsample;fsample<-function(...)seq(.05,.1,length.out=40)
set.seed(944801);t<-sort(runif(48,0,130));idx<-matrix(as.numeric(scale(sin(t/7))),ncol=1)
y<-7+3*cos(2*pi*t/13.2)+.7*idx[,1]+rnorm(48,sd=.3);dy<-rep(.3,48);t<-t-min(t)
for(nh in 1:2){
 set.seed(944802);fit<-MLP(t,y,dy,Indices=idx,Nma=0,Nar=0,Nh=nh,Nsamp=2)
 set.seed(944802);none<-MLP(t,y,dy,Indices=NULL,Nma=0,Nar=0,Nh=nh,Nsamp=2)
 stopifnot(max(abs(fit$power-none$power))>1)
 # Independent noise-only Gaussian profile over a nonnegative jitter.
 X<-cbind(1,t,idx);loss<-function(j){f<-lm.wfit(X,y,1/(dy^2+j^2));-sum(dnorm(f$residuals,sd=sqrt(dy^2+j^2),log=TRUE))}
 o<-optimize(loss,c(0,10),tol=1e-10);j<-if(loss(0)<o$objective)0 else o$minimum
 r<-lm.wfit(X,y,1/(dy^2+j^2))$residuals
 marginal<-function(X){G<-crossprod(X,X/dy^2);res<-lm.wfit(X,r,1/dy^2)$residuals;
   ncol(X)/2*log(2*pi)-as.numeric(determinant(G,logarithm=TRUE)$modulus)/2-sum(res^2/dy^2)/2}
 base<-marginal(cbind(1,t));ref<-vapply(fit$P,function(p){H<-do.call(cbind,lapply(1:nh,function(h)cbind(cos(h*2*pi*t/p),sin(h*2*pi*t/p))));marginal(cbind(H,1,t))-base},numeric(1))
 stopifnot(max(abs(ref-fit$power))<1e-3)
}
fsample<-saved.fsample
cat('MLP selected proxy changes the curve and matches independent Gaussian-sub marginalization\n')
