library(phylin)
vignette("phylin_tutorial")

data(vipers)
data(d.gen)
data(grid)

dim(d.gen)
example(d.gen)
r.dist <- dist(vipers[,1:2])

gv <- gen.variogram(r.dist, d.gen)

gv %>% plot()

gv <- gv.model(gv)
plot(gv)

#"To circumvent this problem we use the tree to clas- sify the samples to either 1 or 0, 
#which means that the sample belongs or not, respectively, to a particular cluster/lineage"

lin <- as.integer(vipers$lin == 2)
int.krig <- krig(lin, vipers[,1:2], grid, gv, neg.weights = FALSE, verbose=FALSE)

grid.image(int.krig, grid, main='Kriging with genetic distances',
          xlab='Longitude', ylab='Latitude',
          sclab='Lineage interpolation') 
points(vipers[,1:2], pch=lin+1)

grid.image(int.krig$sd, grid, main='Kriging with genetic distances',
           xlab='Longitude', ylab='Latitude',
           sclab='Lineage interpolation') 

points(vipers[,1:2], pch=lin+1)



