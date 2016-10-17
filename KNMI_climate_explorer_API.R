library(dplyr)
library(lubridate)
library(httr)
library(zoo)
library(ggplot2)

options(timeout = 120)
 
readURL <- function(url, url_cache) {
  # url contains the URL that forces the selected station data into server cache
  # url-cache contains the URL of the page where the cached data can be accessed and downloaded from
  #browser()
  data <- tryCatch({
    h <- handle(url)
    #POST(url, body = "", handle = h) # get the website to force the data into cache server-side
    GET(url, handle = h) # get the website to force the data into cache server-side
    Sys.sleep(3) # wait to complete call to KNMI
    download.file(url = url_cache, destfile = "./data/tempfile.txt")# get the cached page with data and store as tempfile.txt
  },
  error = function(cond) {
    return(NA)
  })
}

# https://climexp.knmi.nl/gdcnprcpall.cgi?id=paul.v.oppen@gmail.com&WMO=NLM00006380&STATION=MAASTRICHT&extraargs=
 
GetKNMIData <- function(Station, type, Station.Data) {
  #browser()
  # type: Min.Temp, Max.Temp, Avg.Temp, Prec
  base_url <- switch(
    type,
    Min.Temp = "https://climexp.knmi.nl/gdcntmin.cgi?id=paul.v.oppen@gmail.com&WMO=",
    Max.Temp = "https://climexp.knmi.nl/gdcntmax.cgi?id=paul.v.oppen@gmail.com&WMO=",
    Avg.Temp = "https://climexp.knmi.nl/gdcntave.cgi?id=paul.v.oppen@gmail.com&WMO=",
    Prec = "https://climexp.knmi.nl/gdcnprcpall.cgi?id=paul.v.oppen@gmail.com&WMO="
  )
  url.col <- paste0("URL.", tolower(type))
  URL <- paste0(base_url, Station.Data[Station.Data$Station.name == Station, "Station.code"], "&STATION=", Station.Data[Station.Data$Station.name == Station, "Station.name"])
  URL_cache <- Station.Data[Station.Data$Station.name == Station, url.col] #compose correct KNMI explorer URL for data location
  #print(URL)
  #print(URL_cache)
  #Data <- tryCatch(readLines(con = URL), error =  function(e) return()) # read the data from the connection
  readURL(URL, URL_cache)
  Data <- readLines(con = "./data/tempfile.txt")
  Data <- Data[-c(1:5)] # drop header rows (1 to 5)
  Data <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", Data, perl=TRUE) # remove all leading  and trailing zeros
  Data <- matrix(unlist(strsplit(Data, " ")), byrow = TRUE, ncol = 4) # split the data by space into matrix columns
  Data <- as.data.frame(apply(Data, 2, as.numeric), stringsAsFactors = FALSE) # coerce matrix to data frame
  colnames(Data) <- c("Year", "Month", "Day", type)
  return(Data)
}

PadMissingDates <- function(DF) {
  # find and pad missing data values, Fill the missing values with 0. Generates a dataframe without any missing dates
  # create a Date column from the Year, Month and Day values
  DF$Date <- ISOdate(DF$Year, DF$Month, DF$Day) #ISOdate returns a POSIXct object
  # sort by date
  Df <- DF[order(DF$Date), ]
  first.date <- DF$Date[1] # as sorted, the first element holds the earliest date
  last.date <- DF$Date[length(DF$Date)] # as sorted, the last elemenet holds the last date
  # create a new vector (data frame) that holds every single date between first.date and last.date
  all.dates <- data.frame(Date = seq(first.date, last.date, by="day"), stringsAsFactors = FALSE)
  # recreate Year, Month and Day from the dates as new columns
  all.dates$Year <- year(all.dates$Date)
  all.dates$Month <- month(all.dates$Date)
  all.dates$Day <- day(all.dates$Date)
  # merge all.dates with existing DF and pad missing values with NA
  merged.DF <- merge(all.dates, DF, all = TRUE)
  #merged.DF[is.na(merged.DF)] <- 0 # replace NA with 0
  merged.DF <- merged.DF[, -1] # drop the date column
  # add station column
  # merged.DF$Station
  return(merged.DF)
}

FillMissingData <- function(DF){
  
  # fill in missing data
  DF$Min.Temp.splined <- na.approx(DF$Min.Temp, na.rm = FALSE, maxgap = 3, rule = 2)
  DF$Max.Temp.splined <- na.approx(DF$Max.Temp, na.rm = FALSE, maxgap = 3, rule = 2)
  DF$Avg.Temp.splined <- na.approx(DF$Avg.Temp, na.rm = FALSE, maxgap = 3, rule = 2)
  DF$Prec.splined <- na.approx(DF$Prec, na.rm = FALSE, maxgap = 3, rule = 2)
  return(DF)
}


HI <- function(DF) {
  #browser()
  # prepare the climate dataset for huglin index calculations
  DF <- DF %>% 
    mutate(HI.date = as.Date(paste0(Day, "-", Month, "-", Year), format = "%d-%m-%Y"))
  DF <- DF %>% 
    mutate(HI.dayofyear = yday(HI.date))
  
  # only calculate cumsum between 1/4 (day 91) and 9/30 (day 273)
  DF <- DF %>% mutate(HI.range = ifelse(Month < 4, 0, ifelse(Month >= 4 & Month <= 9, 1, 0)))
  # HI for non-splined data
  DF <- DF %>% mutate(HI.day = 1.064 * HI.range * pmax(((Avg.Temp - 10) + (Max.Temp - 10))/2, 0))
  DF <- DF %>% group_by(Year, StationName) %>% mutate(HI.cumsum = cumsum(HI.day))
  # HI for splined data
  DF <- DF %>% mutate(HI.day.splined = 1.064 * HI.range * pmax(((Avg.Temp.splined - 10) + (Max.Temp.splined - 10))/2, 0))
  DF <- DF %>% group_by(Year, StationName) %>% mutate(HI.cumsum.splined = cumsum(HI.day.splined))
  
  #remove first 90 days and last 92 days
  DF <- DF %>% filter((Month >= 4)&(Month <= 9))
  # correct for leap years
  # length of current year
  #last.year.max <- max(DF[DF$Year == max(DF$Year), "HI.dayofyear"])
  #DF$HI.dayofyear <- c(rep( seq(91, 273, 1), length(unique(DF$Year)) - 1), seq(91, last.year.max - 1, 1))
  DF$Type <- "Data"
  browser()
  tmp.max <- DF[1, ] # tmp.max is just an empty copy of DF
  tmp.max[1, ] <- NA
  #max <- ddply(DF, ~HI.dayofyear, summarise, max = max(HI.cumsum, na.rm = TRUE))# HI.cumsum contains the max value for a given HI.dayofyear across all years
  max <- DF %>% group_by(StationName, HI.dayofyear) %>% do(data.frame(max = max(.$HI.cumsum, na.rm = TRUE))) # use data.frame within do to avoid that max is class list
  for (i in 1:nrow(max)) {
    tmp.max[i, ] <- NA
  }
  tmp.max <- data.frame(tmp.max, stringsAsFactors=FALSE)
  tmp.max$HI.cumsum <- max$max
  tmp.max$HI.dayofyear <- max$HI.dayofyear
  tmp.max$Type <- "Max"
  tmp.max$StationName <- max$StationName

  
  tmp.min <- DF[1, ] # tmp.min is just an empty copy of DF
  tmp.min[1, ] <- NA
  #min <- ddply(DF, ~HI.dayofyear, summarise, min = min(HI.cumsum, na.rm = TRUE))# HI.cumsum contains the min value for a given HI.dayofyear across all years
  min <- DF %>% group_by(StationName, HI.dayofyear) %>% do(data.frame(min = min(.$HI.cumsum, na.rm = TRUE)))# use data.frame within do to avoid that min is class list
  for (i in 1:nrow(min)) {
    tmp.min[i, ] <- NA
  }
  tmp.min <- data.frame(tmp.min, stringsAsFactors=FALSE)
  tmp.min$HI.cumsum <- min$min
  tmp.min$HI.dayofyear <- min$HI.dayofyear
  tmp.min$Type <- "Min"
  tmp.min$StationName <- min$StationName
  
  
  #DF <- rbind(DF, min.df, max.df )
  DF <- bind_rows(DF, tmp.max, tmp.min )
  
  return(DF)
}

Station.Data <- read.csv("./data/Station Weather Data.csv", stringsAsFactors = FALSE)

for (i in 1:nrow(Station.Data)) {
  Station_min <- GetKNMIData(Station.Data$Station.name[i], "Min.Temp", Station.Data)
  Station_avg <- GetKNMIData(Station.Data$Station.name[i], "Avg.Temp", Station.Data)
  Station_max <- GetKNMIData(Station.Data$Station.name[i], "Max.Temp", Station.Data)
  Station_prec <- GetKNMIData(Station.Data$Station.name[i], "Prec", Station.Data)
  
  Station_min <- PadMissingDates(Station_min)
  Station_avg <- PadMissingDates(Station_avg)
  Station_max <- PadMissingDates(Station_max)
  Station_prec <- PadMissingDates(Station_prec)
  
  Station_daily_data <- merge(Station_min, Station_max, by = c("Year", "Month", "Day"), all = TRUE)
  Station_daily_data <- merge(Station_daily_data, Station_avg, by = c("Year", "Month", "Day"), all = TRUE)
  Station_daily_data <- merge(Station_daily_data, Station_prec, by = c("Year", "Month", "Day"), all = TRUE)
  Station_daily_data <- Station_daily_data %>% arrange(Year, Month, Day)
  Station_daily_data$StationName <- Station.Data$Station.name[i]
  Station_daily_data$StationCode <- Station.Data$Station.code[i]
  #print(i) 
  #print(Station.Data$Station.name[i])
  if (i == 1) {
    Daily_data <- Station_daily_data
  } else {
    Daily_data <- bind_rows(Daily_data, Station_daily_data)
  }
}

Daily_data <- Daily_data[Daily_data$StationName != "", ]
saveRDS(Daily_data, file = "./data/Raw_KNMI_data.rds")

Daily_data <- FillMissingData(Daily_data)
Daily_data_HI <- HI(Daily_data)
#replace -inf by NA
Daily_data_HI <- do.call(data.frame,lapply(Daily_data_HI, function(x) replace(x, is.infinite(x),NA)))


write.csv(Daily_data_HI, file = "./data/station_data.csv")
saveRDS(Daily_data_HI, file = "./data/Huglin_index_new.rds")


p <- ggplot(data = Daily_data, aes(x = HI.date, y = Min.Temp - Min.Temp.splined, group = StationName))
p <- p + geom_line(colour = "blue")
#p <- p + geom_line(aes(y = Min.Temp.splined), colour = "red")
p <- p + facet_wrap(~StationName)
p