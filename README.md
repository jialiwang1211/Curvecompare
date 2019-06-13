# Curvecompare
Fitting Bayesian penalized smoothing splines for longitudinal data, and compare between different groups.

The function `curveMCMC' fits group-specific Bayesian penalized smoothing splines for longitudinal data for each treatment group. 

The function `plotCurve' plots two curves corresponding to two levels from a single treatment group or from two treatment groups,
and their contrast curve with credible intervals.

## Installation:
library(devtools)\
install_github("jialiwang1211/Curvecompare")\
library(Curvecompare)

## Example
The `growth' data is obtained from 'fda' package (https://cran.r-project.org/web/packages/fda/index.html), which contains Berkeley growth data for 54 girls and 39 boys at 31 ages 
from 1-18 years old. We fit penalized smoothing splines for boys and girls respectively and compare the difference between 
their growth curves. 


library(fda)\
library(reshape2)\
data(growth)\
df_f<-melt(growth$hgtf)\
df_m<-melt(growth$hgtm)\
df<-rbind(df_f,df_m)\
names(df)<-c("age","indiv","height")\
df$gender<-c(rep("F",dim(df_f)[1]), rep("M",dim(df_m)[1]))\
rs<-curveMCMC(df=df,t="age",y="height",Trt1="gender",Trt2=NA,GrpVar="indiv")\
#plot growth curves with 95% credible intervals for boys and girls respectively, and their contrast curve with 95% credible intervals.\
plotCurve(MCMCsamples=rs$MCMCsamples,compare1="F",
compare2="M",numIntKnots=10,plot.contrast=TRUE,plot.indiv=TRUE,data=rs$data)
