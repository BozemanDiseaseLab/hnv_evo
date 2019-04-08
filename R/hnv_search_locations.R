#obain locations for sequences

#okay how this should go is grab specific shapefiles for each point and then repeat sample from
#within that shape file, as a sensitivity analysis to assess how much being off in in our gps
#location affects results. for now lets just grab the google bounding box...should be okay for now

library(RgoogleMaps)
library(ggmap)
key == ...look it up

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
      df$lon[i]   <- rd$lon
      df$lat[i]   <- rd$lat
    }
  }
}


#this finds the diag of the square, .5*diag^2  = area of the square
geodesic <- function(df)
{
  for(i in 1:nrow(df))
  {
    df[i,'area_of_uncertainty_km2'] <-  .5*(fields::rdist.earth(x1 =  as.matrix(df[i,c('north', 'west')]), x2 =  as.matrix(df[i,c('south', 'east')])))^2
  }
  return(df)
}

df <- geodesic(df)

hist(df$area_of_uncertainty_km2, breaks =50)

df$lat <- df$lat + rnorm(n = nrow(df), mean = 0, sd = .01)
df$lon <- df$lon + rnorm(n = nrow(df), mean = 0, sd = .01)

save(df, file ='data/df.Rdata')

