#' @title Retrieval of Weather Data  from WMO and NASA Surface meteorology and Solar Energy
#' 
#' @description The function retrieves weather data from the World Meteorological Organization
#'  and NASA SSE. The retrieval is based on airports  IATA or ICAO codes or WMO weather station 
#'  codes. The actual WMO data is from Weather Underground . Internet connection  is therefore 
#'  needed before data can be accessed
#'
#' @param data dataframe that can contain all or part of the input data 
#' @param wmo numeric. World Meteorological Organization weather station code.
#' You can use the  Worldwide Station List at 
#' \url{https://www.wunderground.com/about/faq/international_cities.asp}  [Accessed on 2016-5-6] or 
#' \url{http://www.wetterzentrale.de/klima/stnlst.html} [Accessed on 2016-5-6]
#' @param airport numeric. IATA or ICAO code. They are alphabetically  listed at 
#' \url{https://en.wikipedia.org/wiki/List_of_airports_by_IATA_code:_A} [Accessed on 2016-5-6]
#' @param date date in the form of "YYYY-mm-dd" for example "2016-01-01"
#' @param time time in the form of "H:M" or decimal hours only. For example "14:21" or 14.35
#' @param NASA.SSE list of  NASA SSE data dates in the form of "YYYY-mm-dd".
#' For example NASA.SSE=list(from=""2001-01-01"", to="2001-04-01"). Note that SSE date
#' is from 1983 to 2005: 
#' \url{https://eosweb.larc.nasa.gov/cgi-bin/sse/sse.cgi?skip@larc.nasa.gov+s01#s01}
#' [Accessed on 2016-5-6] 
#' if the  parameter is not provided, date parameter will be used 
#' provided longitude and latitude are provided.
#' @param folder The path of the folder where the data can be written. 
#' @inheritParams ETo
#' @author George Owusu
#'
#' @return WMO and NASE.SSE. The WMO contains hour and day data. The NASE.SSE contains 
#' meta and data. See example below.
#' @references 
#' \itemize{
#' \item{}{Worldwide Station List.\url{http://www.wetterzentrale.de/klima/stnlst.html} 
#' [Accessed on 2016-5-6] }
#' \item{}{List of airports by IATA code. 
#' \url{https://en.wikipedia.org/wiki/List_of_airports_by_IATA_code} [Accessed on 2016-5-6]} 
#' \item{}{Weather Underground \url{https://www.wunderground.com/} [Accessed on 2016-5-6]}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' #all hours data at Lincoln Airport, Nabraska
#' lincoln=weather(airport = "KLNK", date="2011-08-04")
#' lincoln$WMO$hour
#' lincoln$WMO$day
#' #current hour data at London
#' london.now=weather(airport = "LON", date=format(Sys.time(), "%Y-%m-%d"),
#' time=format(Sys.time(), "%H:%M"))
#' london.now$WMO$hour
#' #NASA SSE data for Kumasi, Ghana
#' kumasi=weather(longitude = -1.5,latitude = 6.7,NASA.SSE = list(from="2001-04-01",
#' to="2002-04-01"))
#' 
#' #Hour, daily and SSE data at JFK 
#' JFK=weather(airport = "JFK", date="2001-02-02",time="17:20", latitude=40.6413, 
#' longitude=73.7781)
#' JFK$WMO$hour
#' JFK$WMO$day
#' JFK$NASA.SEE$data
#' JFK$NASA.SEE$meta
#' 
#' #' #Hour, daily and SSE slice data at Kotoka International Airport, ACCRA 
#' accra=weather(airport = "ACC", date="2001-02-02",time="17:20", latitude=40.6413, longitude=73.7781, 
#' NASA.SSE=list(from="2001-01-01",to="2004-01-01"))
#' accra$WMO$hour
#' accra$WMO$hour
#' accra$NASA.SEE$data
#' accra$NASA.SEE$meta
#' 
#' #Using wmo code
#' kumasi.now=weather(wmo="65442", date=format(Sys.time(), "%Y-%m-%d"),time=format(Sys.time(), 
#' "%H:%M"))
#' kumasi.now$WMO$hour
#' 
#' #retrieve daily data for a specific period i.e 1 month
#' lincoln=weather(airport = "KLNK", date=c("2011-08-04","2011-09-04"))
#' #kumasi and accra data for a time slice. May not work for different 
#' geographical regions due different nummber of paramters. Write the results ti a file
#' accrakumasi=weather(airport =c("DGSI","ACC"), date=c("2011-08-04","2011-10-04"),
#' folder="C:/Users/george/Documents/")
#' }

#now data
#kumasi.now=weather(wmo="65442", date=format(Sys.time(), "%Y-%m-%d"),time=format(Sys.time(), "%H:%M"))
#kumasi.now$WMO$hour
weather=function (data=NULL,wmo=NULL, airport=NULL, date="YYYY-m-d",time=NULL,
  latitude=NULL,longitude=NULL,NASA.SSE=list(from="YYYY-m-d",to="YYYY-m-d"),folder=NULL){
  options(warn=-1)
  if(!is.null(data)){
    if(!is.null(data$date)){
      date=unique(data$date)
    }
    if(!is.null(data$time)){
      time=unique(data$time)
    }
    
    if(!is.null(data$wmo)){
      wmo=unique(data$wmo)
    }
    if(!is.null(data$airpot)){
      airpot=unique(data$airpot)
    }
    if(!is.null(data$latitude)){
      latitude=unique(data$latitude)
    }
    if(!is.null(data$longitude)){
      longitude=unique(data$longitude)
    }
    if(!is.null(data$NASA.SSE)){
      NASA.SSE=unique(data$NASA.SSE)
    }
  }
  
 
  if(length(date)>1){
    date=c(date[1],date[length(date)]) 
  }
 
  if(length(wmo)>1||length(airport)>1){
    
    if(length(wmo)>1){
      ID=wmo
    }else{
      ID=airport
    }
    
   
    
   
     # if(!is.null(date$from)||!is.null(date$to)){
      #  date=seq(as.Date(date$from), as.Date(date$to), "days")
     # }else{
      #  
      #}
      
 hour=NULL
 day=NULL
    for (i in 1:length(ID)){
      if(length(wmo)==length(ID)){
        thiswmo=wmo[[i]]
        thisairport=NULL
      }else{
        thisairport=airport[[i]]
        thiswmo=NULL
      }
      if(length(time)>1){
        thistime=time[[i]]
        
      }else{
        thistime=time
      }
      
      if(length(latitude)>1){
        thislatitude=latitude[[i]]
        
      }else{
        thislatitude=latitude
      }
      if(length(longitude)>1){
        thislongitude=longitude[[i]]
        
      }else{
        thislongitude=longitude
      }
     
      model=weather(data=NULL,wmo=thiswmo, airport=thisairport, date=date,time=thistime,
         latitude=NULL,longitude=NULL,NASA.SSE=list(from="YYYY-m-d",to="YYYY-m-d"))
      #if(length(day)==length(model$WMO$day)){
      day=sebkc.tryCatch(rbind(day,model$WMO$day))$value
      #}
      hour=sebkc.tryCatch(rbind(hour,model$WMO$hour))$value
      
    }
 if(!is.null(folder)){
   write.csv(hour,paste(folder,"/WMO.hour.csv",sep=""))
   write.csv(day,paste(folder,"/WMO.day.csv",sep=""))
    }
 list(WMO=list(hour=hour,day=day))
 
    
  }else{
  NASA=NULL  
  meta=NULL
  
 if(!is.null(time)){
   if(is.character(time)&&grep(":",time)){
   time=strsplit(time,":")
   time=as.numeric(time[[1]][1])+(as.numeric(time[[1]][2])/60)
   }
   time=as.numeric(time)
 }
  
  date=as.character(date)
  date1=date
  if(date=="YYYY-mm-dd"){
    date=readline(prompt = "Enter Date in format YYYY-m-d:")
    #return(print("Date in format YYYY-mm-dd is needed"))
  }
  
  date=gsub("-","/",date)
#WMO data
  hour.data=NULL
  day.data=NULL
  station=airport
  if(!is.null(wmo)||!is.null(airport)){
  if(!is.null(airport)){
  if(length(date)<=1){
    
  url.hr=paste("https://www.wunderground.com/history/airport/",airport,"/", date,"/DailyHistory.html?format=1", sep="")
   url.day=paste("https://www.wunderground.com/history/airport/",airport,"/", date,"/CustomHistory.html?format=1", sep="")
  }else{
    date1a=date[1]
    date1b=date[2]
    date1=date1a
    enddate=strsplit(date1b,"/")[[1]]
    url.hr=paste("https://www.wunderground.com/history/airport/",airport,"/", date1a,"/DailyHistory.html?format=1", sep="")
    url.hr=paste("https://www.wunderground.com/history/airport/",airport,"/", date1a,"/DailyHistory.html?dayend=",enddate[3],"&monthend=",enddate[2],"&yearend=",enddate[1],"&format=1", sep="")
    url.day=paste("https://www.wunderground.com/history/airport/",airport,"/", date1a,"/CustomHistory.html?dayend=",enddate[3],"&monthend=",enddate[2],"&yearend=",enddate[1],"&format=1", sep="")
  }
  }
  if(!is.null(wmo)){
    station=wmo
    #65442 kumasi
    if(length(date)<=1){
      
      url.hr=paste("https://www.wunderground.com/history/airport/",wmo,"/", date,"/DailyHistory.html?format=1", sep="")
      url.day=paste("https://www.wunderground.com/history/airport/",wmo,"/", date,"/CustomHistory.html?format=1", sep="")
    }else{
      date1a=date[1]
      date1b=date[2]
      date1=date1a
      enddate=strsplit(date1b,"/")[[1]]
      url.hr=paste("https://www.wunderground.com/history/airport/",wmo,"/", date1a,"/DailyHistory.html?format=1", sep="")
      url.hr=paste("https://www.wunderground.com/history/airport/",wmo,"/", date1a,"/DailyHistory.html?dayend=",enddate[3],"&monthend=",enddate[2],"&yearend=",enddate[1],"&format=1", sep="")
      url.day=paste("https://www.wunderground.com/history/airport/",wmo,"/", date1a,"/CustomHistory.html?dayend=",enddate[3],"&monthend=",enddate[2],"&yearend=",enddate[1],"&format=1", sep="")
    }
   # url.hr=paste("https://www.wunderground.com/history/wmo/",wmo,"/", date,"/DailyHistory.html?format=1", sep="")
   # url.day=paste("https://www.wunderground.com/history/wmo/",wmo,"/", date,"/CustomHistory.html?format=1", sep="")
    }
  
   hour.data=sebkc.tryCatch(read.csv(url.hr,stringsAsFactors=FALSE))$value
   #print(hour.data)
   if(class(hour.data)[1]=="simpleError"){
     print(paste("WMO Hourly data not available for ",wmo,airport, "on", date ))
   }
   day.data=sebkc.tryCatch(read.csv(url.day,stringsAsFactors=FALSE))$value
   if(class(day.data)[1]=="simpleError"){
     print(paste("WMO Daily data not available for ",wmo,airport, "on", date ))
   }
   
     
     if(class(day.data)[1]!="simpleError"){
       thisnames=names (day.data)
       #imperial units
       if(length(grep("TemperatureF",thisnames)>0)){
         day.data$Tmax=sebkc.tryCatch((as.numeric(day.data$Max.TemperatureF)-32)/1.8)$value
         day.data$Tmin=sebkc.tryCatch((as.numeric(day.data$Min.TemperatureF)-32)/1.8)$value
                  day.data$uz=as.numeric(day.data$Mean.Wind.SpeedMPH)*0.44704 
         day.data$P=as.numeric(as.character(day.data$PrecipitationIn)) *25.4
        # DOY=strptime(day.data[[1]],"%m/%d/%Y")
         day.data$DOY=strftime(strptime(day.data[[1]],"%m/%d/%Y"),"%Y-%m-%d")
       }else{
         #metric units
         day.data$Tmax=sebkc.tryCatch((as.numeric(day.data$Max.TemperatureC)))$value
         day.data$Tmin=sebkc.tryCatch((as.numeric(day.data$Min.TemperatureC)))$value
         #print("yes")
         
         #print(names(day.data))
         day.data$uz=as.numeric(day.data$Mean.Wind.SpeedKm.h)*0.27778 
         day.data$P=(day.data$Precipitationmm)
         day.data$P=as.numeric(as.character(day.data$P))
         day.data$DOY=as.character(day.data[[1]])
       }
       
       day.data$P[is.na(day.data$P)]=0.001
       day.data$RHmin=day.data$Min.Humidity 
       day.data$RHmax=day.data$Max.Humidity
       #day.data$DOY=as.character(as.Date(day.data[[1]],"%Y-%m-%d"))
       day.data$station=station
       day.data$longitude=longitude
       day.data$latitude=latitude
     }
   if(length(date)>1){
     day.data$station=station
     hour.data$station=station
     if(!is.null(folder)){
       write.csv(hour.data,paste(folder,"/",station,".WMO.hour.csv",sep=""))
       write.csv(day.data,paste(folder,"/",station,".WMO.day.csv",sep=""))
     }
     return(list(WMO=list(hour=hour.data,day=day.data),NASA.SEE=list(data=NASA,meta=meta)))
   }
   if(!is.null(time)){
     time=round(time)
     hour.data$time=round(
       as.numeric(format(strptime(hour.data[[1]], "%I:%M %p"), "%H")) +         # hours component
         (as.numeric(format(strptime(hour.data[[1]], "%I:%M %p"), "%M")) / 60),2  # rounded minutes component
     )

     #hour.data=hour.data[hour.data$time==time,]
     if(class(hour.data)[1]!="simpleError"){
     hour.data$time2=abs(hour.data$time-time)
     thisort=sort(hour.data$time2,decreasing = FALSE)[1:2]
     hour.data2=hour.data[hour.data$time2==thisort[1]|hour.data$time2==thisort[2],]
     hour.data3= hour.data2[1,]
     for (i in 1:length(hour.data2)){
     mod=sebkc.tryCatch(lm(hour.data2[[i]]~hour.data2$time))$value
     if(class(mod)[1]=="simpleError"){
    output=paste(hour.data2[[i]][1],hour.data2[[i]][2],sep="-")  
     }else{
     slope=coef(mod)[[2]]
     intercept=coef(mod)[[1]]
     output=(slope*time)+intercept 
     }
     hour.data3[i]=output
     #hour.data3=rbind(hour.data3,output)
     #print(output)
     }
     hour.data=hour.data3
     hour.data$time2=NULL
     thisnames=names (hour.data)
     #imperial units
     if(length(grep("TemperatureF",thisnames)>0)){
     hour.data$Tmean=(as.numeric(hour.data$TemperatureF)-32)/1.8
     hour.data$uz=as.numeric(hour.data$Wind.SpeedMPH)*0.44704 
     hour.data$P=as.numeric(hour.data$PrecipitationIn)*25.4 
     }else{
       #print(hour.data)
       hour.data$Tmean=(as.numeric(hour.data$TemperatureC))
       if(length(grep("-Calm",hour.data$Wind.SpeedKm.h))>0){
         wind=as.numeric(strsplit(hour.data$Wind.SpeedKm.h,"-")[[1]][1])
       }else{
        wind= hour.data$Wind.SpeedKm.h
       }
       hour.data$uz=as.numeric(wind)*0.27778 
       hour.data$P=as.numeric(hour.data$Precipitationmm) 
     }
     
     hour.data$RH=hour.data$Humidity 
     hour.data$DOY=date1 
     hour.data$station=station
     hour.data$longitude=longitude
     hour.data$latitude=latitude
   }
   
   
   }
  
   #rint(day.data)
   if(class(day.data)[1]!="simpleError"){
     thisnames=names (day.data)
     #imperial units
     if(length(grep("TemperatureF",thisnames)>0)){
   day.data$Tmax=sebkc.tryCatch((as.numeric(day.data$Max.TemperatureF)-32)/1.8)$value
   day.data$Tmin=sebkc.tryCatch((as.numeric(day.data$Min.TemperatureF)-32)/1.8)$value
   if(nrow(hour.data)>0){
     #hour.data$Tmin=NA
     hour.data$Tmin=sebkc.tryCatch(day.data$Tmin)$value
     hour.data$Tmax=sebkc.tryCatch(day.data$Tmax)$value
   }
   
   day.data$uz=as.numeric(day.data$Mean.Wind.SpeedMPH)*0.44704 
   day.data$P=(as.numeric(as.character(day.data$PrecipitationIn) ))*25.4
     }else{
       #metric units
       day.data$Tmax=sebkc.tryCatch((as.numeric(day.data$Max.TemperatureC)))$value
       day.data$Tmin=sebkc.tryCatch((as.numeric(day.data$Min.TemperatureC)))$value
       #print("yes")
       if(nrow(hour.data)>0){
         #hour.data$Tmin=NA
         hour.data$Tmin=sebkc.tryCatch(day.data$Tmin)$value
         hour.data$Tmax=sebkc.tryCatch(day.data$Tmax)$value
       }
       #print(names(day.data))
       day.data$uz=as.numeric(as.character(day.data$Mean.Wind.SpeedKm.h))*0.27778 
       day.data$P=(day.data$Precipitationmm)
       day.data$P=as.numeric(as.character(day.data$P))
     
     }
     #day.data$P=as.numeric(as.character(day.data$P))
     day.data$P[is.na(day.data$P)]=0.001
     day.data$RHmin=day.data$Min.Humidity 
     day.data$RHmax=day.data$Max.Humidity
     day.data$DOY=date1
     day.data$station=station
     day.data$longitude=longitude
     day.data$latitude=latitude
   }
   if(is.null(station)){
     station=NA
   }
  }
   #day
NASA.this=NULL
   #NASA SSE
   if(!is.null(latitude)||!is.null(longitude)||NASA.SSE$from!="YYYY-m-d"){
     NASA.this=TRUE
     if(NASA.SSE$from=="YYYY-m-d"){
       start.date=date
     }else{
       start.date= NASA.SSE$from
     }
     if(NASA.SSE$to=="YYYY-m-d"){
       end.date=date
     }else{
       end.date= NASA.SSE$to
     }
     start.date=gsub("/","-",start.date)
     end.date=gsub("/","-",end.date)
     start.date=strsplit(start.date,"-")[[1]]
     end.date=strsplit(end.date,"-")[[1]]
     start.date=as.numeric(start.date)
     end.date=as.numeric(end.date)
     
    url <- paste("https://eosweb.larc.nasa.gov/cgi-bin/sse/daily.cgi?email=skip%40larc.nasa.gov&step=1&lat=",
                 latitude,"&lon=",longitude,"&sitelev=&ms=",start.date[2],"&ds=",start.date[3],"&ys=",start.date[1],"&me=",end.date[2],"&de=",end.date[3],"&ye=",end.date[1],"&p=swv_dwn&p=avg_kt&p=clr_sky&p=clr_dif&p=clr_dnr&p=clr_kt&p=lwv_dwn&p=toa_dwn&p=PS&p=TSKIN&p=T10M&p=T10MN&p=T10MX&p=Q10M&p=RH10M&p=DFP10M&submit=Submit&plot=swv_dwn",sep="")
   html <- sebkc.tryCatch(paste(readLines(url), collapse="\n"))$value
   
   sensor=strsplit(html, "<a href=")
   link=paste("https://eosweb.larc.nasa.gov",strsplit(sensor[[1]][9],"\"")[[1]][2],sep="")
   meta=sebkc.tryCatch(read.delim(link,sep=" ",stringsAsFactors = F,header=F))$value[1:12,]
   #print(meta[1:12,])
   elevation=sebkc.tryCatch(as.numeric(read.delim(link,sep=" ",stringsAsFactors = F,header=F)[4,][10][[1]]))$value
   NASA=sebkc.tryCatch(read.table(file=link,skip=23,header=TRUE))$value 
   if(class(NASA)[1]=="simpleError"){
     print(paste("NASA SSE data not Available in ",start.date[1]))
   }
   NASA$altitude=elevation
   NASA$Tmax=NASA$T10MX
   NASA$Tmin=NASA$T10MN
   NASA$DOY=paste(NASA$YEAR,NASA$MO,NASA$DY,sep="-")
   NASA$Rs=(NASA$swv_dwn*41.666666666666664)*0.0864
   }
   
   if(!is.null(folder)){
     write.csv(hour.data,paste(folder,"/",station,".WMO.hour.csv",sep=""))
     write.csv(day.data,paste(folder,"/",station,".WMO.day.csv",sep=""))
     if(!is.null(NASA.this)){
       write.csv(NASA,paste(folder,"NASA.data.csv"))
       write.csv(day.data,paste(folder,"NASA.meta.csv"))
       write.csv(NASA,paste(folder,"/",station,".NASA.data.csv",sep=""))
       write.csv(meta,paste(folder,"/",station,".NASA.meta.csv",sep=""))
     }
     
   }
   list(WMO=list(hour=hour.data,day=day.data),NASA.SEE=list(data=NASA,meta=meta))
  }
}
#all hours 
#lincoln=weather(airport = "KLNK", date="2011-08-04")
#lincoln=weather(airport = "KLNK", date=format(Sys.time(), "%Y-%m-%d"),time=format(Sys.time(), "%H:%M"))
#london=weather(airport = "LON", date=format(Sys.time(), "%Y-%m-%d"),time=format(Sys.time(), "%H:%M"))

#lincoln$WMO$hour
#NASA SSE data
#kumasi=weather(wmo="65442", date="2016-05-05",longitude = -1.5,latitude = 6.7,NASA.SSE = list(from="2001-04-01",to="2002-04-01"))

#now data
#kumasi.now=weather(wmo="65442", date=format(Sys.time(), "%Y-%m-%d"),time=format(Sys.time(), "%H:%M"))
#kumasi.now$WMO$hour

