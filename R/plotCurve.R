#' Plot individual curves and their constrast curve
#'
#' Plot two curves corresponding to two levels from a single treatment group or from two treatment groups, and their contrast curve with credible intervals.
#'
#' @param MCMCsamples The object from \code{curveMCMC} which contains MCMC samples.
#'
#' @param compare1 Level in the treatment groups. If there is only 1 treatment group when fitting \code{curveMCMC}, it should be specified as compare1="a", where "a" is a factor level in Trt1. If there are two treatment groups when fitting  \code{curveMCMC}, it should be specified as compare1=c("a","b"), where "a" is a factor level in Trt1 and "b" is a factor level in Trt2.
#'
#' @param compare2 Level in the treatment groups. Similar to \code{compare1}.
#'
#' @param numIntKnots number of interior knots. Must be the same as fitting \code{curveMCMC}.
#'
#' @param data The data frame from fitting \code{curveMCMC}, which is sorted by Trt1 and Trt2.
#'
#' @param plot.contrast plot two mean curves defined by \code{compare1} and \code{compare2}, and their contrast curve (curve1 minus curve2) with credible intervals. TRUE or FALSE.
#'
#' @param plot.indiv plot two curves defined by \code{compare1} and \code{compare2} with respective credible intervals. TRUE or FALSE.
#'
#' @param credible_level level of credible interval.
#'
#' @examples #compare growth curves for females and males
#' plotCurve(MCMCsamples=rs$MCMCsamples,compare1="F",
#'compare2="M",numIntKnots=10,plot.contrast=TRUE,plot.indiv=TRUE,data=rs$data)
#'
#' @author Jiali Wang (\email{jiali.wang@@data61.csiro.au})
#'
#' @export
#' @import ggplot2
#' @import HRW
plotCurve<-function(MCMCsamples=rs$MCMCsamples,compare1,
                    compare2,numIntKnots,data=rs$data,plot.contrast,plot.indiv,credible_level=0.95){

if(length(compare1)==1){
  curve1_id<-data[data$Trtid%in%compare1,"Trtidnum"][1]
  curve2_id<-data[data$Trtid%in%compare2,"Trtidnum"][1]
}else{
curve1_id<-data[data$Trtid%in%interaction(compare1[1],compare1[2]),"Trtidnum"][1]
curve2_id<-data[data$Trtid%in%interaction(compare2[1],compare2[2]),"Trtidnum"][1]
}

curve_id<-c(curve1_id,curve2_id)

x<-data$t
y<-data$y

a<-min(x)
b<-max(x)

numObs=length(y)
numGrp=length(unique(data$GrpVar))
ncZ=numIntKnots+2
Trtidnum=data$Trtidnum
idnum=data$GrpVar
numCurv<-length(unique(Trtidnum))

#intercept
for (j in 1:2){
  charVar <- paste("beta[",as.character(curve_id[j]),"]",sep = "")
  assign(paste("beta_intc",as.character(j),"MCMC",sep = ""),MCMCsamples[,1,dimnames(MCMCsamples)$parameters == charVar])
}

#slope
for (j in 1:2){
  charVar <- paste("beta[",as.character(numCurv+curve_id[j]),"]",sep = "")
  assign(paste("beta_slp",as.character(j),"MCMC",sep = ""),MCMCsamples[,1,dimnames(MCMCsamples)$parameters == charVar])
}


#spl random effects
for (j in 1:2){
  uGrpMCMC <- NULL
  for (i in 1:ncZ)
  {
    charVar <- paste("uGrp[",as.character(curve_id[j]),",",as.character(i),"]",sep="")
    uGrpMCMC <- rbind(uGrpMCMC,as.vector(MCMCsamples[,1,dimnames(MCMCsamples)$parameters == charVar]))
  }
  assign(paste("u",as.character(j),"GrpMCMC",sep = ""),uGrpMCMC)
}

#######################
ng <- 201
xg <- seq(a,b,length=ng)

Xg <- cbind(rep(1,ng),xg)
intKnots <-  quantile(unique(x),seq(0,1,length = numIntKnots+2)
                      [-c(1,numIntKnots+2)])
Zg <- ZOSull(xg,intKnots = intKnots,range.x = c(a,b))

fhat1<-Xg%*%rbind(beta_intc1MCMC,beta_slp1MCMC)+Zg%*%u1GrpMCMC
fhat1_mean<-apply(fhat1,1,mean)
fhat1_low<-apply(fhat1,1,quantile,(1-credible_level)/2)
fhat1_up<-apply(fhat1,1,quantile,credible_level+(1-credible_level)/2)

fhat2<-Xg%*%rbind(beta_intc2MCMC,beta_slp2MCMC)+Zg%*%u2GrpMCMC
fhat2_mean<-apply(fhat2,1,mean)
fhat2_low<-apply(fhat2,1,quantile,(1-credible_level)/2)
fhat2_up<-apply(fhat2,1,quantile,credible_level+(1-credible_level)/2)

fdiff_mean<-apply(fhat1-fhat2,1,mean)
fdiff_low<-apply(fhat1-fhat2,1,quantile,(1-credible_level)/2)
fdiff_up<-apply(fhat1-fhat2,1,quantile,credible_level+(1-credible_level)/2)

p_contrast<-p_indiv<-NA

if(plot.contrast==TRUE){
  data_gg<-data.frame(time=rep(xg,3),trt=rep(c("Curve1","Curve2","Contrast"),each=length(xg)),xg=rep(xg,3),f=c(fhat1_mean,fhat2_mean,fdiff_mean),credLow=c(rep(NA,length(xg)),rep(NA,length(xg)),fdiff_low),credUpp=c(rep(NA,length(xg)),rep(NA,length(xg)),fdiff_up))

  p_contrast<-ggplot(data=data_gg,aes(time,f,group=trt))+
  geom_line(aes(color=trt),size=1)+
  geom_ribbon(aes(ymin=credLow,ymax=credUpp,group=trt,fill=trt),alpha=0.3,show.legend =FALSE)+
  geom_hline(yintercept=0,size=1) +
  scale_colour_manual(values = c("red", "blue", "green"))+
  labs(x="t" ,y="")+
  theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())
}

if(plot.indiv==TRUE){
  data_gg<-data.frame(time=rep(xg,2),trt=rep(c("Curve1","Curve2"),each=length(xg)),xg=rep(xg,2),f=c(fhat1_mean,fhat2_mean),credLow=c(fhat1_low,fhat2_low),credUpp=c(fhat1_up,fhat2_up))

  p_indiv<-ggplot(data=data_gg,aes(time,f,group=trt))+
    geom_line(aes(color=trt),size=1)+
    geom_ribbon(aes(ymin=credLow,ymax=credUpp,group=trt,fill=trt),alpha=0.3,show.legend =FALSE)+
    labs(x="t" ,y="")+
    theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())
}
return(list(p_contrast=p_contrast,p_indiv=p_indiv))
}
