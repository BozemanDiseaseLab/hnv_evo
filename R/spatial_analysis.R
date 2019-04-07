#distance analysis

library(phylin)

ordering <- sort(attr(hnv.dist, "Labels"))
ordering.d <- sort(attr(d, "Labels"))
ordering == ordering.d
hnv.dist <- as.matrix(hnv.dist)[ordering, ordering]
d <- as.matrix(d)[ordering.d, ordering.d]
hnv.dist  <- as.dist(hnv.dist)
d  <- as.dist(d)

attr(hnv.dist, "Labels") == attr(d, "Labels")

gv <- gen.variogram(hnv.dist, d, lag = 0.008)
gv %>% plot()

gv.plot <- gv.model(gv)
plot(gv)

