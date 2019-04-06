#obain locations for sequences

#okay how this should go is grab specific shapefiles for each point and then repeat sample from
#within that shape file, as a sensitivity analysis to assess how much being off in in our gps
#location affects results. for now lets just grab the google bounding box...should be okay for now

library(RgoogleMaps)
library(ggmap)
key <- "AIzaSyD40FLStnvp085UB-FKXbyONuxV3ke4umY"
ggmap::register_google(key = key)

i<-1
for(i in 1:nrow(df))
{
  query <- as.character(df$country[i])
  rd <- geocode(query, output = 'more', source = 'google')
  if(is.atomic(rd)==FALSE)
  {
    if (is.na(rd$lon) == FALSE) 
    {
      df$north[i] <- rd$north
      df$south[i] <- rd$south
      df$east[i]  <- rd$east
      df$west[i]  <- rd$west
    }
  }
}

save(df, file ='data/df')

