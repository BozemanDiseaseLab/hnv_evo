
library(tidyverse)


# generate some points

# distance matrix between points
x <- seq(0, 2, length.out = 75)
y <- seq(0, 2, length.out = 75)
d1 <- expand.grid(x = x, y =y)

dd <- as.matrix(dist(d1))

plot(d1$x,d1$y)
# scoring.matrix

phi <- 10
sigma  = .30
tau = .1
plot(tau + sigma*exp(phi*-seq(1,1000,by=.1)),type ='l')

scoring.matrix <- sigma*exp(- phi * dd )
diag(scoring.matrix) <- tau + sigma

scoring.matrix <- scoring.matrix + min(scoring.matrix)
min(scoring.matrix)
scoring.matrix  <- scoring.matrix / max(scoring.matrix)
min(scoring.matrix)
max(scoring.matrix)

#outcome <- t(scoring.matrix) %*% rnorm(nrow(d1),0,1) 

cutoff <- max(dd, na.rm=TRUE) / 3 # default maximum distance
num.bins <-  15
bin.width <- cutoff / 15

dist.sq <- (scoring.matrix)  #dist(aq$pm2.5)^2
#dist.sq <- as.matrix(dist(outcome)^2)  #dist(aq$pm2.5)^2

dist.mat <- dd

dist.mat[upper.tri(dist.mat,diag=TRUE)] <- NA

vario.dat <- data.frame(dist = as.numeric(dist.mat), diff = as.numeric(dist.sq)) %>%
  dplyr::arrange(desc(dist))

vario.dat <- vario.dat %>% 
  mutate(bin = floor(dist / bin.width) + 1) 

vario.dat2 <- vario.dat %>% 
  filter(bin < 20) %>%
  group_by(bin) %>% 
  dplyr::summarize(emp.sv = .5 * mean(diff), n = n()) 

vario.dat2 %>% ggplot(aes(x=bin, y= emp.sv))  + 
  geom_point() 

hist(log(scoring.matrix))

phi <- 0
sigma  = 0
tau = .1
plot(tau + sigma*exp(phi*-seq(1,1000,by=.1)),type ='l')


scoring.matrix <- sigma*exp(- phi * dd )
diag(scoring.matrix) <- tau + sigma

scoring.matrix <- scoring.matrix + min(scoring.matrix)
min(scoring.matrix)
scoring.matrix  <- scoring.matrix / max(scoring.matrix)
min(scoring.matrix)
max(scoring.matrix)

#outcome <- t(scoring.matrix) %*% rnorm(nrow(d1),0,1) 

cutoff <- max(dd, na.rm=TRUE) / 3 # default maximum distance
num.bins <-  15
bin.width <- cutoff / 15

dist.sq <- (scoring.matrix)  #dist(aq$pm2.5)^2
#dist.sq <- as.matrix(dist(outcome)^2)  #dist(aq$pm2.5)^2

dist.mat <- dd

dist.mat[upper.tri(dist.mat,diag=TRUE)] <- NA

vario.dat <- data.frame(dist = as.numeric(dist.mat), diff = as.numeric(dist.sq)) %>%
  dplyr::arrange(desc(dist))

vario.dat <- vario.dat %>% 
  mutate(bin = floor(dist / bin.width) + 1) 

vario.dat2 <- vario.dat %>% 
  filter(bin < 20) %>%
  group_by(bin) %>% 
  dplyr::summarize(emp.sv = .5 * mean(diff), n = n()) 

vario.dat2 %>% ggplot(aes(x=bin, y= emp.sv))  + 
  geom_point() 

hist((scoring.matrix))


#okay so semivariogram, matrix on matrix violence
# next step, is the matrix multiplication followed by jags gonna be the same estimate? lets do it






