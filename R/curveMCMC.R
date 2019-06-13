#' Fit group-specific curves with MCMC
#'
#' Fit Bayesian penalized smoothing splines for longitudinal data for each level within a treatment group. Cubic O'Sullivan splines are used as basis functions in the linear mixed model formulation. Random intercepts are added for each group.
#'
#'The current version only supports one or two treatments, both of which must be factors. If there is only one  treatment, then set \code{Trt2}=NA. An individual curve is fitted to each level of the (combination of) treatments.
#'The MCMC is fitted from stan, which interfaces with R using `rstan` package.
#'
#' @param df data frame.
#'
#' @param t variable name in df which denotes time.
#'
#' @param y response variable.
#'
#' @param Trt1 variable name of a factor in df which denotes the first treatment group. See in Details.
#'
#' @param Trt2 variable name of a factor in df which denotes the second treatment group. Set NA if no second treatment group.
#'
#' @param GrpVar variable name in df which denotes the group variable.
#'
#' @param numIntKnots number of interior knots.
#'
#' @param sigmaBeta hyperparameter, variance in Normal disctribution of slope coefficient.
#'
#' @param Au hyperparameter, scale in half Cauchy distribution of variance of random effects, assuming the random effects are iid.
#'
#' @param AU hyperparameter, scale in half Cauchy distribution of variance of random intercept.
#'
#' @param Aeps hyperparameter, scale in half Cauchy distribution of variance of residuals.
#'
#' @param nwarm number of MCMC iterations as burn-in.
#'
#' @param nsamp number of MCMC iterations.
#'
#' @param number of MCMC iterations.
#'
#' @return An object contains MCMC samples is returned.
#'
#' \item{MCMCsamples}{A list contains saved MCMC samples of slope, random effects, random intercepts, variance of random effects, variance of random intercept, variance of residuals.}
#'
#' \item{data}{The data frame which is sorted by Trt1 and Trt2, and contains the new group coding by taking each combination of the treatment factors.}
#'
#' @examples library(fda)
#'library(reshape2)
#'#Berkeley growth data for 54 girls and 39 boys at 31 ages from 1-18 years old.
#'data(growth)
#'df_f<-melt(growth$hgtf)
#'df_m<-melt(growth$hgtm)
#'df<-rbind(df_f,df_m)
#'names(df)<-c("age","indiv","height")
#'df$gender<-c(rep("F",dim(df_f)[1]), rep("M",dim(df_m)[1]))
##
#'rs<-curveMCMC(df=df,t="age",y="height",Trt1="gender",Trt2=NA,GrpVar="indiv")
#'
#' @author Jiali Wang (\email{jiali.wang@@data61.csiro.au})
#'
#' @references
#' \insertRef{harezlak2018semiparametric}{Curvecompare}
#'
#' @export
#' @import HRW
#' @import rstan
#' @import dplyr
#' @importFrom Rdpack reprompt
curveMCMC<-function(df,t,y,Trt1,Trt2=NA,GrpVar,numIntKnots=10,sigmaBeta=1e5,Au=1e5,AU=1e5,Aeps=1e5,nwarm=500, nsamp=1000,odens=1){

data<-df
names(data)[names(data)==t] <-"t"
names(data)[names(data)==y] <-"y"
names(data)[names(data)==Trt1] <-"Trt1"
names(data)[names(data)==Trt2] <-"Trt2"
names(data)[names(data)==GrpVar] <-"GrpVar"

if(is.na(Trt2)){
  data<-arrange(data,Trt1,GrpVar,t)
  data$Trtid<-data$Trt1
}else{
  data<-arrange(data,Trt1,Trt2,GrpVar,t)
  data$Trtid<-interaction(data$Trt1,data$Trt2)
}

data$Trtidnum<-as.numeric(as.factor(data$Trtid))

data$GrpVar<-as.numeric(as.factor(data$GrpVar))

x<-data$t
y<-data$y

X<-cbind(model.matrix(~as.factor(Trtidnum)+0,data=data),x*model.matrix(~as.factor(Trtidnum)+0,data=data))

a<-min(x)
b<-max(x)

intKnots<-quantile(unique(x),seq(0,1,length=numIntKnots+2)
                      [-c(1,numIntKnots+2)])

Zspl<- ZOSull(x,intKnots = intKnots,range.x = c(a,b))

numObs=length(y)
numGrp=length(unique(data$GrpVar))
ncX=ncol(X)
ncZ=ncol(Zspl)
Trtidnum=data$Trtidnum
idnum=data$GrpVar

numCurv<-length(unique(Trtidnum))

# Specify model in Stan:
Model <-
  'data
{
  int<lower=1> numObs;         int<lower=1> numGrp;
  int<lower=1> ncX;            int<lower=1> ncZ;
  real<lower=0> sigmaBeta;     real<lower=0> AU;
  real<lower=0> Aeps;          real<lower=0> Au;
  vector[numObs] y;            int<lower=1> idnum[numObs];
  int<lower=1> Trtidnum[numObs]; int<lower=1> numCurv;
  matrix[numObs,ncX] X;        matrix[numObs,ncZ] Zspl;
}
parameters
{
  vector[ncX] beta;            vector[numGrp] U;
  matrix[numCurv,ncZ] uGrp;
  real<lower=0> sigmaU;
  real<lower=0> sigmau;        real<lower=0> sigmaEps;
}
model
{
  y ~ normal(X*beta + U[idnum] + rows_dot_product(uGrp[Trtidnum],Zspl),sigmaEps);
  U ~ normal(0,sigmaU);           to_vector(uGrp)  ~ normal(0,sigmau);
  beta  ~ normal(0,sigmaBeta) ;   sigmaEps ~ cauchy(0,Aeps);
  sigmaU ~ cauchy(0,AU) ;   sigmau ~ cauchy(0,Au);
}'

allData<-list(numObs=numObs,numGrp=numGrp,ncX=ncX,ncZ=ncZ,sigmaBeta=sigmaBeta,sigmaBeta=sigmaBeta,
              AU=AU,Aeps=Aeps,Au=Au,y=y,idnum=idnum,
              Trtidnum=Trtidnum,numCurv=numCurv,X=X,Zspl=Zspl)

stanObj <-stan(model_code=Model,data=allData,warmup=nwarm,
               iter=nsamp,chains=1,thin=odens,refresh=500)

MCMCsamples <- extract(stanObj,permuted = FALSE)

return(list(MCMCsamples=MCMCsamples,data=data))
}
