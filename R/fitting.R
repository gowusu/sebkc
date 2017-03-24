#' @title Goodness of fit Test
#'
#' @param obs Observed data
#' @param est Estimated data
#' @param P Number of parameters
#'
#' @return
#' \itemize{
#' \item{r2:} {  R-square [-]}
#'  \item{r2_adj:} { An adjusted R-square [-]}
#' \item{d:} { Index of agreement[-]}
#' \item{p_value:} { p value}
#' \item{sse:} { Sum of square error}
#' \item{sst:} { Sum of square total}
#' \item{AIC:} { AIC}
#' \item{BIC:} { BIC}
#' \item{RMSE:} { Root Mean Square Error}
#' \item{RMSE_adj:} { adjusted Root Mean Square Error}
#' \item{NRMSE} { Normalised Root Mean Square Error}
#' \item{NRMSE_adj:} { adjusted Normalised Root Mean Square Error}
#' \item{E:} { Nash Sutcliffe Coefficient}
#' \item{b:} { coefficients of regression}
#' \item{AAE:} { average absolute error}
#' \item{ARE:} { average relative error}
#' \item{MBE:} { Mean Bias  error}
#' \item{MPE:} { Mean Percentage error}
#' }
#' @export
#'
#' @examples
#' file=system.file("extdata","sys","irrigation.txt",package="sebkc")
#' data=read.table(file,header=TRUE)  
#' obs=data$ETo
#' est=data$ETc
#' P=4
#' mod=fitting(obs,est,P)
#' r2=mod$r2
fitting<-function(obs,est,P=1)
{
  
  data=data.frame(obs=obs,est=est)
  data=na.omit(data)
  obs=data$obs
  est=data$est
  N=length(obs)
 # b=sum(obs*est)/sum(obs^2)
  AAE=sum(abs(obs-est))/N
  ARE=sum(abs((obs-est)/obs))*(100/N)
  MBE=sum((est-obs))/N
  #MPE=sum(((obs-est)/obs))*(100/N)
  MPE=sum(((est-obs)/obs)*100)/N
  
  mod=lm(est~obs)
  #print(summary(mod))
  #print("p-value")
  x=mod
  F_value=summary(mod)$fstatistic[[1]]
  df=summary(mod)$fstatistic[[3]]
  r2=summary(mod)[[8]]
  f <- summary(mod)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  #Nash-Sutcliffe efficiency rating for model estiamtes https://github.com/USGS-R/rloadest/blob/master/R/nashSutcliffe.R
  #nashSutcliffe <- function(obs, est, na.rm=TRUE) {
  na.rm=TRUE
  
  E <- 1 - sum((obs - est)^2, na.rm=na.rm)/sum((obs-mean(obs,na.rm=na.rm))^2, na.rm=na.rm)
  #root mean square error
  #RMSE=sqrt(mean(mod$residuals^2))
  sse=sum((obs-est)^2)
  RMSE=sqrt(sse/N)
  NRMSE=RMSE/(max(obs)-min(obs))
  #AIC
  AIC=AIC(mod)
  #print(P)
  rmse_adj=sqrt((sse/(N-P)))
  NRMSE_adj=rmse_adj/(max(obs)-min(obs))
  
  adj=N-P
  sst=sum((obs-mean(obs))^2)
 # r2_adj=1-((sse/(adj-1))/(sst/(N-1)))
  r2_adj=r2-(1-r2)*(P/(N-P-1))
  #AIC=N*log(sse)+(2*P)
  
  AIC=(N*(log(2*pi)+log(sse/adj)+1))+P
  #Priestley, M.B. (1981), Spectral Analysis and Time Series, Academic Press. ISBN 0-12-564922-3 (p. 375)
  BIC=(N*log(sse/N))+P*log(N)
  #parameters=as.data.frame(list(r2=r2,r2_adj=r2_adj,N=N, p=P,F_value=F_value, df=df, p_value=p,AIC=AIC,BIC=BIC,RMSE=RMSE,RMSE_adj=rmse_adj,NRMSE=NRMSE,NRMSE_adj=NRMSE_adj,nashSutcliffe=E))
  Psquare=sum((est-obs)^2)
  obsmean=(abs(obs-mean(obs)))
  estmean=(abs(est-mean(obs)))
  obs.est=(obsmean+estmean)^2
  #d=1-((sum(est-obs)^2)/sum(((abs(est-mean(obs)))+(abs(obs-mean(obs))))^2))
  d=1-(Psquare/sum(obs.est))
  parameters=as.data.frame(list(r2=r2,r2_adj=r2_adj,d=d,NSE=E, p_value=p,sse=sse,sst=sst,AIC=AIC,BIC=BIC,RMSE=RMSE,RMSE_adj=rmse_adj,NRMSE=NRMSE,NRMSE_adj=NRMSE_adj,nashSutcliffe=E,
                                b=coef(mod)[[2]],AAE=AAE,ARE=ARE,MBE=MBE,MPE=MPE,r=cor(obs,est), intercept=coef(mod)[[1]],N=N, parameters=P,F_value=F_value, df=df,E=E))
  #  #print(parameters)
  
}