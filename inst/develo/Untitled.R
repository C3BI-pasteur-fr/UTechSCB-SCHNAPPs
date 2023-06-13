# normal distribution

data = rnorm(100000,sd = 1000,mean = 0)

hist(data)

adata = asinh(data)

hist(adata)


minb = 1
maxb = 500
stepb = 5

library(dplyr)
for (b in seq(minb, maxb, stepb)){
  adata = asinh(data/b)
  hist(adata, main = paste("b:",b),breaks = 1000,plot = T)
}

# standard deviation is the driving factor

for( sd in c(1,10,50,100,200,500,1000,10000,1000000)) {
  data = rnorm(100000,sd = sd,mean = 0)
  for (b in c(1,5,10,100,500,1000,5000)){
    adata = asinh(data/b)
    hist(adata,breaks = 1000, main=paste("sd:",sd, "b=",b))
  }
}

# 
# sd:1, b=5
# sd:10, b=100
# sd:50, b=200
# sd:100, b=500
# sd:200, b=800
# sd:500, b=5000
# sd:1000, b=5000
# sd:10,000, b>5000

# if the mean is slightly off, what are the effects
# shift to the negative, negative peak is much higher.
# the effecto of b seems to be independant of the mean

for( sd in c(1,100,500)) {
  for(mean in c(-100, -10,0,10, 100)){
    data = rnorm(100000,sd = sd,mean = mean)
    for (b in c(1,5,500,5000)){
      adata = asinh(data/b)
      hist(adata,breaks = 1000, main=paste("sd:",sd, "b=",b, "mean:", mean))
    }
  }
}

# what about the number of cells?
# seems to have no effect...
sd=500
b=500
mean=100
for (nz in c(10000,20000,30000,40000,50000,100000,1000000,1e10)){
  data = rnorm(100000,sd = sd,mean = mean)
  adata = asinh(data/b)
  hist(adata,breaks = 1000, main=paste("sd:",sd, "b=",b, "mean:", mean,"n:",nz))
}



`
That means that only the standdeiviation has an effect with regards to b. 
this might have some biological reason,
could be related to the age of the fluorchromes or others.
`

# Does the mean / median of the negative values says something about the std?
mean = 0
sd = 500
outDF = data.frame(mean = numeric(),
                   sd = character(),
                   b= numeric(),
                   meanVal = numeric(),
                   medianVal = numeric()
)

minb = 1
maxb = 5000
stepb = 5

for( sd in c(1,10,50,100,200,500,1000,10000,1000000)) {
  data = rnorm(100000,sd = sd,mean = 0)
  for (b in seq(minb, maxb, stepb)){
    as = asinh(data/b)
    as = as[as<0]
    m1 = mean(as)
    m2 = median(as)
    outDF[nrow(outDF) + 1,] = c(mean = mean, sd=sd, b=b, meanVal = m1, medianVal = m2)
  }
}

outDF$sd = as.factor(outDF$sd)
outDF$b = as.numeric(outDF$b)
outDF$medianVal = as.numeric(outDF$medianVal)
outDF$meanVal = as.numeric(outDF$meanVal)

outDF = outDF[!is.null(outDF$medianVal),]
library(ggplot2)
# this is IT!!! this shows that the value of b=1 can be used to determine the 
# optimal value...
ggplot(outDF, aes(x=b, y=meanVal, color=sd)) + geom_line() 

# But what is the optimal value???
library(parallel)
cl <- makeCluster(detectCores()-1) #not to overload your computer
doParallel::registerDoParallel(cl)
library(BiocParallel)
library(foreach)
b=5000
sd = 500
bVals = data.frame(sd = numeric(), b= numeric(), shPval=numeric())
# for(sd in seq(1, 5000, 10)){
out = foreach (sd = seq(1, 5000, 10), .combine = rbind) %dopar% {
  for (b in 1:10000){
    data = rnorm(5000,sd = sd,mean = 0)
    adata = asinh(data/b)
    hist(adata)
    sh = shapiro.test(adata) 
    if( sh$p.value >0.1) {
      bVals[nrow(bVals)+1,] = c(sd,b,sh$p.value)
      return(c(sd,b,sh$p.value))
    }
  }
}
out = as.data.frame(out)
colnames(out) = c("sd", "b", "p-value")
plot(out$sd,out$b)
lm(b~sd, out)

data.frame(sd= numeric(), b1 = numeric(), b10=numeric())
b=5000
sd = 500
bVals = data.frame(sd = numeric(), b= numeric(), shPval=numeric())
# for(sd in seq(1, 5000, 10)){
out2 = foreach (sd = seq(1, 5000, 1), .combine = rbind) %dopar% {
  data = rnorm(5000000,sd = sd,mean = 0)
  b=1
  as = asinh(data/b)
  as = as[as<0]
  m1 = mean(as)
  m2 = median(as)
  return(c(sd, b=1, m1))
}

out = as.data.frame(out)
colnames(out) = c("sd", "b", "p-value")
plot(out$sd,out$b)
lm(b~sd, out)



