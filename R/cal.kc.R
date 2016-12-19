  
#' @title Calibration of crop coefficients
#' 
#' @description This function can estimate up to 23 parameters of FAO single and dual crop
#' parameters. The parameters include soil parameters: WP,FC, REW; 
#' crop parameters: kcini,kcmid,kcend,kcbini,kcbmid,kcbend,h,Zr; 
#' the evaporation parameters: p and Ze, 
#' runoff coefficients: CNII or rc, as well as groundwater and capillary rise 
#' parameters: initgw,a1,b1,a2,b2,a3,b3,a4,b4. For each parameter initial, maximum and range 
#' values must be provided if one wants to optimise it.
#' There are also 12 goodness of fit tests for the evaluation of the model with
#' \code{\link{fitting}} :r-square, index of agreement, F-value, p-value, sse,sst, AIC,
#' BIC, RMSE, NRMSE and nashSutcliffe coefficients. A calibrated model that uses more 
#' parameters gets punished by AIC and BIC as well as adjusted R2, RMSE, NRMSE.
#' 
#' Both point and spatial model can be calibrated. The actual model is available at \code{\link{kc}}
#'
#' @inheritParams kc
#' @inheritParams ETo
#' @param ETa Measured Actual Evapotranspiration [mm]
#' @author George Owusu
#' @details  
#' \describe{
#' There are several ways  to use cal.kc:
#'  \item{\strong{Measured Data}}{The kc model can be calibrated against either Actual Evapotranspiration (ETa) or
#'  Available Soil Water; the user must specify theta. If the calibration is done on a point data then time series for
#'  that location must be provided. If multiple locations are provided with xydata and corresponding 
#'  multiple locations as input data, then ETa or theta must be provided for all locations. See also \code{kc} 
#'  for more details on locations.     }
#'  \item{\strong{Number of optimised parameters and Goodness of fit}}{ 
#'  There is a need to limit the number of calibrated parameters.
#'   One can compare models by observing the adjusted goodness of fit tests. 
#'   The more parameters that one uses the less are the values of R-square but
#'   RMSE, BIC, AIC increases. Unadjusted and adjusted goodness of fit tests
#'    are provided with this function.
#'     }
#'   }
#' @seealso \code{\link{kc}}
#'
#' @return The function returns optimised parameters as "parameters", goodness of fit 
#' as "fit", and calibrated model as "mod". It also prints candidates of parameters and 
#' goodness of fit of the local optimum values. For the explanation of calibrated model 
#' output and optimised parameters see the return values of \code{\link{kc}}. Here, only
#' goodness of fit tests explanations are given.
#' \itemize{
#' \item{mod:} { The optimised model. See return of \code{\link{kc}} for more details}
#' \item{parameters:} { parameters of the optimised model}
#' \item{fit:} {  Goodness of fit tests. See \code{\link{fitting}} for explanation of variables}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' #fewer iteration
#' file=system.file("extdata","sys","irrigation.txt",package="SEBKc")
#' data=read.table(file,header=TRUE)  
#' P=data$P
#' rc=0
#' I=data$I
#' ETo=data$ETo
#' Zr=data$Zr
#' p=0.6
#' FC=0.23
#' WP=0.10
#' u2=1.6  
#' RHmin=35
#' soil="sandy loam"
#' crop="Broccoli"
#' #single kc calibration with fewer iteration
#' cal=cal.kc(ETo,P,RHmin,soil,crop,I=0,kc=NULL,CNII=67,
#' u2=2,FC=NULL,WP=NULL,h=NULL,Zr=NULL,lengths=NULL,fc=NULL,
#' kctype="single",region=NULL,regexpr=NULL,initgw=50,LAI=3,
#' rc=NULL,Rn=NULL,hmodel="gaussian",Zrmodel="linear",
#' xydata=NULL,ETa=ETo,a1=360,b1=-0.17,a2=240,b2=-0.27,
#' a3=-1.3,b3=6.6,a4=4.6,b4=-0.65)
#' parameters=cal$parameters
#' output=cal$mod$output
#' gof=cal$fit
#' ETc=cal$mod$output$ETc
#' TAW=cal$mod$output$TAW
#' 
#' #Single kc with different kcs. kc uses default values
#' #RMSE_adj and NRMSE_adj equal to NAN because the number of parameters more than time steps(days)
#' cal=cal.kc(ETo,P,RHmin,soil,crop,I=0,CNII=67,u2=2,FC=NULL,WP=NULL,
#' h=NULL,Zr=NULL,lengths=NULL,fc=NULL,kctype="single",region=NULL,
#' regexpr=NULL,initgw=50,a1=c(100,400,10),b1=c(-1,-0.1,0.1),
#' a2=c(100,400,10),b2=c(-1,1,0.1),a3=c(-2,1,0.1),b3=c(1,10,1),
#' a4=c(1,10,1),b4=c(-1.5,1,0.1),LAI=3,rc=NULL,
#' Rn=NULL,hmodel="gaussian",Zrmodel="linear",xydata=NULL,ETa=ETo)
#' 
#' #Example with fixed kcs
#' cal=cal.kc(ETo,P,RHmin,soil,crop,I=0,CNII=67,kc=c(0.1,1.15,1),
#' u2=2,FC=NULL,WP=NULL,h=NULL,Zr=NULL,lengths=NULL,
#' fc=NULL,kctype="single",region=NULL,regexpr=NULL,initgw=50,
#' b1=c(-1,-0.1,0.1),b2=c(-1,1,0.1),LAI=3,rc=NULL,Rn=NULL,
#' hmodel="gaussian",Zrmodel="linear",xydata=NULL,ETa=ETo)
#' 
#' #Example with spatial kcs
#' #landsat folder with original files but a subset
#' folder=system.file("extdata","stack",package="SEBKc") 
#' sebiauto=sebi(folder=folder,welev=317,Tmax=31,Tmin=28) #sebi model
#' kc2=stack(sebiauto$EF/2,sebiauto$EF,sebiauto$EF/3) #assign EF to kc
#' xydata=data.frame(cbind(longitude=data$longitude,latitude=data$latitude))
#' modEF=kc(ETo,P=P,RHmin,soil,crop,I,kc=kc2,Rn=8,lengths = c(4,4,2,2),xydata=xydata)
#' ETa=modEF$ETcxy[3:14]

#' 
#' calsp=cal.kc(ETo,P,RHmin,soil,crop,I=0,CNII=67,kc=kc2,u2=2,FC=NULL,
#'  WP=WP,h=NULL,Zr=NULL,lengths = c(4,4,2,2),fc=NULL,kctype="single",
#'  region=NULL,regexpr=NULL,initgw=50,b1=c(-1,-0.1,0.1),b2=c(-1,1,0.1),
#'  LAI=3,rc=NULL,Rn=7,hmodel="gaussian",Zrmodel="linear",xydata=xydata,ETa=ETa)
#'  gof=cal$fit
#' r2=cal$fit$r2
#' AIC=cal$fit$AIC
#' ETcxy=cal$mod$ETcxy
#' TAW=cal$mod$TAW
#' }
  cal.kc=function(ETo,P,RHmin,soil,crop,I=0,CNII=c(30,100,1),u2=2,FC=NULL,WP=NULL,h=NULL,
                   Zr=NULL,Ze=c(0.1,0.15,0.01),kc=c(c(0.1,0.3,0.05),c(0.1,1,0.1),
    c(0.1,1,0.1)),p=c(0.1,1,0.1),lengths=NULL,fw=0.8,fc=NULL,
    kctype="dual",region=NULL,regexpr="sweet",initgw=50,
    LAI=3,rc=NULL,Kb=0,Rn=NULL,hmodel="gaussian",Zrmodel="linear",xydata=NULL,
    REW=c(2,12,0.5),a1=360,b1=c(-1,1,0.1),a2=240,b2=c(-1,1,0.1),
    a3=c(-2,1,0.1),b3=6.6,a4=4.6,b4=c(-1.5,1,0.1),ETa=NULL,report="stages"
    ,time="06:00",HWR=1,latitude=NULL,density="full",DOY=NULL,slope=NULL,
    y=NULL,r1=100,nongrowing="weeds",theta=NULL,profile="variable",Kcmin=0.15,Ky=NULL){
    
    if(is.null(ETa)&&is.null(theta)){
      return  (print("Observation data ETa or theta is needed"))
    }
    
    soilwater=read.csv(system.file("extdata","sys","soilwater2.csv",package="SEBKc"),header=TRUE,stringsAsFactors=FALSE)
    stages=(read.csv(system.file("extdata","sys","stages.csv",package="SEBKc"),stringsAsFactors=FALSE,header=TRUE))
    crop_H_Zr=na.omit(read.csv(system.file("extdata","sys","crop_H_Zr.csv",package="SEBKc"),stringsAsFactors=FALSE,header=TRUE))
    singlekc=na.omit(read.csv(system.file("extdata","sys","singlekc.csv",package="SEBKc"),stringsAsFactors=FALSE,header=TRUE))
    dualkc=na.omit(read.csv(system.file("extdata","sys","dualkc.csv",package="SEBKc"),stringsAsFactors=FALSE,header=TRUE))
    ky=na.omit(read.csv(system.file("extdata","sys","ky.csv",package="SEBKc"),stringsAsFactors=FALSE,header=TRUE))
    wetting=na.omit(read.csv(system.file("extdata","sys","fw.csv",package="SEBKc"),stringsAsFactors=FALSE,header=TRUE))
    soil2=soil
    crop2=crop
    
    EF=NULL
    para=NULL
    count=0
    if(class(kc[[1]])=="RasterLayer"||class(kc[[2]])=="RasterLayer"
       ||class(kc[[3]])=="RasterLayer"){
      EF=TRUE 
    }
    
    soil=toupper(soil)
    soil=soilwater[(grep(soil2,soilwater$type,ignore.case=TRUE)),]
    if(is.null(FC)||is.null(WP)){
      if(is.null(FC)){
        FC=mean(c(as.numeric(strsplit(soil$FC, "-")[[1]][1]),as.numeric(strsplit(soil$FC, "-")[[1]][2])))
      }
      if(is.null(WP)){
        WP=mean(c(as.numeric(strsplit(soil$WP, "-")[[1]][1]),as.numeric(strsplit(soil$WP, "-")[[1]][2])))
      }
      
    }
    if(is.null(Zr)||is.null(p)){
      crop= crop_H_Zr[(grep(crop2,crop_H_Zr$crop,ignore.case=TRUE)),] 
      if(nrow(crop)>1){
        crop= crop[(grep((regexpr),(crop),ignore.case=TRUE)),] 
      }
      if(nrow(crop)==0){
        return (print(paste(crop2," not found [for crop type name]"))) 
      }
      if(is.null(Zr)){
        Zr=mean(c(as.numeric(strsplit(crop$Zr, "-")[[1]][1]),as.numeric(strsplit(crop$Zr, "-")[[1]][2])))
      }
      
      if(is.null(p)){
        p=crop$p
      }
    }
    
    
    ##dual crop coefficient
    if(is.null(kc)){
      kcs= dualkc[(grep((crop2),(dualkc$crop),ignore.case=TRUE)),] 
      if(nrow(kcs)>1){
        kcs= kcs[(grep((regexpr),(kcs$crop),ignore.case=TRUE)),] 
      }
      
      if(nrow(kcs)==0){
        return (print(paste(" Kcs not found for", crop2))) 
      }
      kcbini=mean(as.numeric(kcs$Kcbini))
      kcbmid=mean(as.numeric(kcs$Kcbmid))
      kcbend=mean(as.numeric(kcs$Kcbend))
      kc=c(kcbini,kcbmid,kcbend)
    }
    
    #kcb=(day*kcb)/max(day)
    
    if(is.null(h)){
      h= singlekc[(grep((crop2),(singlekc$crop),ignore.case=TRUE)),] 
      if(nrow(h)==0){
        return (print(paste(" h not found for", crop2))) 
      }
      if(nrow(h)>1){
        h= h[(grep((regexpr),(h$crop),ignore.case=TRUE)),] 
      }
      h=mean(as.numeric(h$h))
    }
    if(is.null(REW)){
      REW=mean(c(as.numeric(strsplit(soil$REW, "-")[[1]][1]),as.numeric(strsplit(soil$REW, "-")[[1]][2])))
    }
  #cal.kc=function(object,kcini=c(0.1,0.7,0.05),kcmid=c(0.1,1,0.05),
  #kcend=c(0.1,10,0.05),p=c(0.1,1,0.05),Ze=c(0.1,0.15,0.005),
  #REW=c(2,12,,0.05),ETc,AW=NULL)
    if(is.null(EF)){
    kcmat <- matrix(kc, ncol=3, byrow=TRUE)
    
    
    if(nrow(kcmat)==1){
      kcini=c(kcmat[1],kcmat[1],kcmat[1])
      kcmid=c(kcmat[2],kcmat[2],kcmat[2])
      kcend=c(kcmat[3],kcmat[3],kcmat[3])
    }else{
      kcini=kcmat[1,]
      kcmid=kcmat[2,]
      kcend=kcmat[3,]
      para=cbind(para,kcini=kcini,kcmid=kcmid,kcend=kcend)
      count=count+1
    }
    }else{
      kcini=c(1,1,1)
      kcmid=c(1,1,1)
      kcend=c(1,1,1)
      kc2=kc
    }
    
    if(length(REW)==1){
      REW=c(REW,REW,REW)
    }else{
      para=cbind(para,REW=REW)
      count=count+1
    }
    
    CNIIlist=NULL
    if(!is.null(CNII)){
      if(!is.null(xydata)&&length(CNII)==nrow(xydata)&&!is.null(EF)){
        CNIIlist=CNII
        CNII=c(1,1,1)
        #print("CNIIPPPPPP")
      }else{
        if(length(CNII)==1){
          CNII=c(CNII,CNII,CNII)
        }else{
          para=cbind(para,CNII=CNII)
          count=count+1
        }
      }
    }
    
    hlist=NULL
    if(!is.null(h)){
      if(!is.null(xydata)&&length(h)==nrow(xydata)&&!is.null(EF)){
        hlist=h
        h=c(1,1,1)
        #print("hPPPPPP")
      }else{
        if(length(h)==1){
          h=c(h,h,h)
        }else{
          para=cbind(para,h=h)
          count=count+1
        }
      }
    }
    
    Zrlist=NULL
    if(!is.null(Zr)){
      if(!is.null(xydata)&&length(Zr)==nrow(xydata)&&!is.null(EF)){
        Zrlist=Zr
        Zr=c(1,1,1)
        #print("ZrPPPPPP")
      }else{
        if(length(Zr)==1){
          Zr=c(Zr,Zr,Zr)
        }else{
          para=cbind(para,Zr=Zr)
          count=count+1
        }
      }
    }
    rc21=NULL
    if(!is.null(rc)){
      if(length(rc)==1){
        rc=c(rc,rc,rc)
      }else{
        para=cbind(para,rc=rc)
        count=count+1
      }
      rc21=NULL
    }else{
      rc21=TRUE
      rc=c(1,1,1)
    }
    
    FClist=NULL
    if(!is.null(FC)){
      if(!is.null(xydata)&&length(FC)==nrow(xydata)&&!is.null(EF)){
        FClist=FC
        FC=c(1,1,1)
        #print("FCPPPPPP")
      }else{
        if(length(FC)==1){
          FC=c(FC,FC,FC)
        }else{
          para=cbind(para,FC=FC)
          count=count+1
        }
      }
    }
    WPlist=NULL
    if(!is.null(WP)){
      if(!is.null(xydata)&&length(WP)==nrow(xydata)&&!is.null(EF)){
        WPlist=WP
        WP=c(1,1,1)
      }else{
        if(length(WP)==1){
          WP=c(WP,WP,WP)
        }else{
          para=cbind(para,WP=WP)
          count=count+1
        }
      }
    }
    if(!is.null(p)){
      if(length(p)==1){
        p=c(p,p,p)
      }else{
        para=cbind(para,p=p)
        count=count+1
      }
    }
    
    if(!is.null(Ze)){
      if(length(Ze)==1){
        Ze=c(Ze,Ze,Ze)
      }else{
        para=cbind(para,Ze=Ze)
        count=count+1
      }
    }
    
    initgwlist=NULL
    if(!is.null(initgw)){
      if(!is.null(xydata)&&length(initgw)==nrow(xydata)&&!is.null(EF)){
        initgwlist=initgw
        initgw=c(1,1,1)
        #print("initgwPPPPPP")
      }else{
        if(length(initgw)==1){
          initgw=c(initgw,initgw,initgw)
        }else{
          para=cbind(para,initgw=initgw)
          count=count+1
          }
      }
    }
   
    
    
    if(!is.null(a1)){
      if(length(a1)==1){
        a1=c(a1,a1,a1)
      }else{
        para=cbind(para,a1=a1)
        count=count+1
      }
    }
    
    if(!is.null(b1)){
      if(length(b1)==1){
        b1=c(b1,b1,b1)
      }else{
        para=cbind(para,b1=b1)
        count=count+1
      }
    }
    #print(b1)
    
    if(!is.null(a2)){
      if(length(a2)==1){
        a2=c(a2,a2,a2)
      }else{
        para=cbind(para,a2=a2)
        count=count+1
      }
    }
    
    if(!is.null(b2)){
      if(length(b2)==1){
        b2=c(b2,b2,b2)
      }else{
        para=cbind(para,b2=b2)
        count=count+1
      }
    }
    
    if(!is.null(a3)){
      if(length(a3)==1){
        a3=c(a3,a3,a3)
      }else{
        para=cbind(para,a3=a3)
        count=count+1
      }
    }
    
    if(!is.null(b3)){
      if(length(b3)==1){
        b3=c(b3,b3,b3)
      }else{
        para=cbind(para,b3=b3)
        count=count+1
      }
    }
    if(!is.null(a4)){
      if(length(a4)==1){
        a4=c(a4,a4,a4)
      }else{
        para=cbind(para,a4=a4)
        count=count+1
      }
    }
    
    if(!is.null(b4)){
      if(length(b4)==1){
        b4=c(b4,b4,b4)
      }else{
        para=cbind(para,b4=b4)
        count=count+1
      }
    }
    candidate=NULL
    print("Iteration Started. Kindly wait.................")
    pnames=names(data.frame(para))
    count=length(pnames)
    #return(list(para=para))
    #print(count)
    rmse=100
    r2=0
    kcini1=kcini[1]
    while(kcini1<=kcini[2]){
      kcmid1=kcmid[1]
      while(kcmid1<=kcmid[2]){
        kcend1=kcend[1]
        while(kcend1<=kcend[2]){
          p1=p[1]
          while(p1<=p[2]){
            Ze1=Ze[1]
            while(Ze1<=Ze[2]){
              REW1=REW[1]
              while(REW1<=REW[2]){
                CN1=CNII[1]
                while(CN1<=CNII[2]){
                  FC1=FC[1]
                  while(FC1<=FC[2]){
                    WP1=WP[1]
                    while(WP1<=WP[2]){
                      initgw1=initgw[1]
                      while(initgw1<=initgw[2]){
                        a11=a1[1]
                        while(a11<=a1[2]){
                          b11=b1[1]                      
                          while(b11<=b1[2]){
                            a21=a2[1]
                            while(a21<=a2[2]){
                              
                              b21=b2[1]
                              while(b21<=b2[2]){
                                a31=a3[1]
                                while(a31<=a3[2]){
                                  b31=b3[1]
                                  while(b31<=b3[2]){
                                    a41=a4[1]
                                    while(a41<=a4[2]){
                                      b41=b4[1]
                                      while(b41<=b4[2]){
                                        if(!is.null(EF)){
                                          thiskc2=kc2
                                        }else{
                                          thiskc2=c(kcini1,kcmid1,kcend1)
                                        }
                                          rc1=rc[1]
                                          while(rc1<=rc[2]){
                                            h1=h[1]
                                            while(h1<=h[2]){
                                              Zr1=Zr[1]
                                             while(Zr1<=Zr[2]){
                                                
                                            if(!is.null(rc21)){
                                              thisrc=NULL
                                            }else{
                                              thisrc=rc1
                                            }
                                            
                                            if(!is.null(WPlist)){
                                              WP2=WPlist
                                            }else{
                                              WP2=WP1
                                            }
                                            if(!is.null(FClist)){
                                              FC2=FClist
                                            }else{
                                              FC2=FC1
                                            }
                                            
                                            if(!is.null(initgwlist)){
                                              initgw2=initgwlist
                                            }else{
                                              initgw2=initgw1
                                            }
                                            
                                            if(!is.null(CNIIlist)){
                                              CNII2=CNIIlist
                                            }else{
                                              CNII2=CN1
                                            }
                                              if(!is.null(hlist)){
                                                h2=hlist
                                              }else{
                                                h2=h1
                                              }
                                               if(!is.null(Zrlist)){
                                                 Zr2=Zrlist
                                               }else{
                                                 Zr2=Zr1
                                               }
                                              #print(crop2)
                                          mod=kc(ETo=ETo,P=P,RHmin=RHmin,soil=soil,crop=crop2,I=I,CNII=CNII2,
                                                 u2=u2,FC=FC2,WP=WP2,h=h2,
                                                 Zr=Zr2,Ze=Ze1,kc=thiskc2,p=p1,lengths=lengths,fw=fw,fc=fc,
                                                 kctype=kctype,region=region,regexpr=regexpr,initgw=initgw2,
                                                 LAI=LAI,rc=thisrc,Kb=Kb,Rn=Rn,hmodel=hmodel,Zrmodel=Zrmodel,xydata=xydata,
                                                 REW=REW1,a1=a11,b1=b11,a2=a21,b2=b21,a3=a31,b3=b31,a4=a41,b4=b41,report=report
                                                 ,time=time,HWR=HWR,latitude=latitude,
                                                 density=density,DOY=DOY,slope=slope,y=y,r1=r1,
                                                 nongrowing=nongrowing,theta=theta,profile=profile,Ky=Ky)
                                          #print(nrow(mod$output))
                                          #if(!is.null(names(ETo)))
                                          #{
                                           # obs=mod$output$theta_measured[(length(mod$output$Day)/2)+1:(length(mod$output$Day)*2)]
                                            #first=((length(mod$output$Day)/2)+1)
                                            #second=length(mod$output$Day)
                                            #print(c(first,second))
                                            #mod$output= mod$output[seq(1,second,2),]
                                            #print(mod$output[seq(1,second,2),])
                                         # }
                                          #print(cbind(WP1,FC1,kcini1,kcmid1,kcend1,p1,Ze1,REW1,CN1,FC1,initgw1,a11,b11,a21,b21,a31,b31,a41,b41,rc1))   
                                          #print(length(ETo$ETo))
                                          if(!is.null(EF)){
                                            if(!is.null(ETa)){
                                            obs=rowMeans(ETa)
                                            est=rowMeans(mod$ETc_adj_xy[3:length(mod$ETc_adj_xy)])
                                            }
                                            if(!is.null(theta)){
                                              obs=rowMeans(mod$AW_measuredxy)
                                              est=rowMeans(mod$AW_modelledxy[3:length(mod$AW_modelledxy)])  
                                            }
                                            }else{
                                              if(!is.null(ETa)){
                                              obs=ETa
                                              est=mod$output$ETc_adj
                                              }else{
                                                if(!is.null(mod$output$AW_modelled)){
                                                  obs=mod$output$AW_measured 
                                                  est=mod$output$AW_modelled
                                                  }
                                              }
                                              
                                              }
                                          fit=fitting(obs,est,count)
                                          if(fit$r2>r2){
                                           #print(cbind(WP1,FC1,kcini1,kcmid1,kcend1,p1,Ze1,REW1,CN1,FC1,initgw1,a11,b11,a21,b21,a31,b31,a41,b41,rc1))   
                                        parameters=(list(WP=WP2,FC=FC2,kc=c(kcini1,kcmid1,kcend1),p=p1,Ze=Ze1,REW=REW1,
                                                          CNII=CNII2,initgw=initgw2,h=h2,Zr=Zr2,a1=a11,b1=b11,a2=a21,
                                                         b2=b21,a3=a31,b3=b31,a4=a41,b4=b41,rc=rc1,h=h2,Zr=Zr2,kcini=kcini1,
                                                         kcmid=kcmid1,kcend=kcend1)) 
                                        #print()
                                        
                                        para2=parameters[pnames]
                                            print(data.frame(para2))
                                            if(is.null(candidate)){
                                              candidate=cbind(para2,fit)
                                            }else{
                                              candidate=rbind(crop=crop2,candidate,cbind(para2,fit))
                                          }
                                          rmse=fit$RMSE
                                          r2=fit$r2
                                          thismod2=mod
                                          thisfit=cbind(crop=crop2,para2,fit)
                                          print(thisfit)
                                          
                                          }
                                          Zr1=Zr1+Zr[3]
                                            }
                                          h1=h1+h[3]
                                          }
                                          rc1=rc1+rc[3]
                                          }
                                        if(b4[1]==b4[2]){
                                          b41=b41-b4[3] 
                                        }else{
                                          b41=b41+b4[3] 
                                        }  
                                      }
                                      a41=a41+a4[3]
                                    }
                                    b31=b31+b3[3]
                                  }
                                  if(a3[1]==a3[2]){
                                    a31=a31-a3[3] 
                                  }else{
                                    a31=a31+a3[3] 
                                  }
                                }                              
                                if(b2[1]==b2[2]){
                                  b21=b21-b2[3] 
                                }else{
                                  b21=b21+b2[3] 
                                }
                              }                      
                              a21=a21+a2[3]
                            }
                            
                            if(b1[1]==b1[2]){
                            b11=b11-b1[3] 
                            }else{
                              b11=b11+b1[3] 
                            }
                          }
                          
                          a11=a11+a1[3]  
                        }
                        
                       
                      
                        initgw1=initgw1+initgw[3]  
                      }
                    
                      WP1=WP1+WP[3]
                      
                    }
                  
                  FC1=FC1+FC[3]
                  }
                  CN1=CN1+CNII[3]
                  
                  
                }
                
                REW1=REW1+REW[3]
              }
              Ze1=Ze1+Ze[3]
              
            }
      
          p1=p1+p[3]
          }
          
          kcend1=kcend1+kcend[3] 
        }
        
        kcmid1=kcmid1+kcmid[3] 
      }
      kcini1=kcini1+kcini[3] 
    }
    
    list(fit=thisfit,parameters=para2,mod=thismod2,candidate=candidate)
  }
  
  
  #b1=c(-1,-0.7,0.1)
#kc=c(0.1,1.15,1)