library(foreach)
library(doMC)
registerDoMC(32)  #change the 2 to your number of CPU cores

rm(list=ls())
## constant rate
seeds.array <- c()
##system(cmd)
difs.list <- list()


nsims <- 80000
nfrags <- 2000
tj23 <- runif(nsims, 0.1625, 0.18125)
tj123 <- runif(nsims, 0.625, 0.875)
x2 <- 0.15 ## sampling time for Altai
x3 <- 0.0625 ## sampling time for Vindja
## migration between archaic. Cannot be less than the sampling time
tm23 <- sapply(1:nsims, function(x){runif(1, x2, tj23[x])})
## when migration BETWEEN ARCHAIC should stop [tm23, tj23]
tm23stop <- sapply(1:nsims, function(x){runif(1, tm23[x], tj23[x])})
## start migration between Hs - Altai
tm12 <- sapply(1:nsims, function(x){runif(1, x2, tj23[x])})
tm13 <- sapply(1:nsims, function(x){runif(1, x3, tj23[x])})
## stop migration between Hs - Altai. Here, until the common ancestor of Hs - archaic
tm12stop <- sapply(1:nsims, function(x){runif(1, tm12[x], tj123[x])})
tm13stop <- sapply(1:nsims, function(x){runif(1, tm13[x], tj123[x])})

param.names <- c("tj23", "tj123", "tm23.start", "tm23.stop", "tm12.start", "tm12.stop", "tm13.start", "tm13.stop")
param.vec <- c(tj23, tj123, tm23, tm23stop, tm12, tm12stop, tm13, tm13stop)
params <- matrix(param.vec, ncol=length(param.names), byrow=F)
colnames(params) <- param.names

write.table(params, "params.out", quote=F, col.names = T, row.names = F, sep="\t")

alldatalist <- list()
cmds <- list()
cmds <- foreach(i = 1:nsims) %dopar% {
    print(i)
    difs.list <- list()
    mr <- 100
    mr2 <- 100
    seeds <- paste(sample(100000000, 3), collapse = " ")
    seeds.array <- c(seeds.array, seeds)
    if(i < nsims/2){
        cmd.add <- paste("./ms 3 ", nfrags, " -t 50 -I 3 1 1 1  -x 2 ", x2, " -x 3 ", x3, " -em ",  tm12[i]," 1 2 ", mr, " -em ", tm12stop[i], " 1 2 0 ", " -em ",  tm23[i], " 2 3 ", mr2, " -em ", tm23[i], " 3 2 ",  mr2, " -em ", tm23stop[i], " 2 3 0 ", " -ej ",  tj23[i],  " 3 2 -ej ", tj123[i],  " 2 1   -seeds ", seeds,  "| sed '1,6d' | grep [01] | grep -v [^01]")
    }else{
        cmd.add <- paste("./ms 3 ", nfrags, " -t 50 -I 3 1 1 1  -x 2 ", x2, " -x 3 ", x3, " -em ", tm13[i], " 1 3 ", mr, " -em ", tm13stop[i], " 1 3 0 ", " -em ",  tm12[i]," 1 2 ", mr, " -em ", tm12stop[i], " 1 2 0 ", " -em ",  tm23[i], " 2 3 ", mr2, " -em ", tm23[i], " 3 2 ",  mr2, " -em ", tm23stop[i], " 2 3 0 ", " -ej ",  tj23[i],  " 3 2 -ej ", tj123[i],  " 2 1   -seeds ", seeds,  "| sed '1,6d' | grep [01] | grep -v [^01]")
    }
}

i <- 1
alldata <- foreach(i = 1:nsims) %dopar% {
    print(paste(i, " start 1"))
    cmd.add <- cmds[[i]]
    vall <- system(cmd.add, intern=TRUE)
    j <- 1
    averageSegsites <- 0
    for(j in 1:nfrags){
        ss <- (j-1)*3 + 1
        se <- ss + 2
        v <- vall[ss:se]
        vl <- sapply(v, function(x){ strsplit(x, "")})
        ## ## find all differences between the three sequences
        difs <- apply(combn(3,2), 2, function(x){ sum(vl[[x[1]]] != vl[[x[2]]]) })
        theta.denom <- 1 + 1/2
        averageSegsites <- averageSegsites + length(vl[[1]])
        ##correction.factor <- segsites/theta.denom
        difs.list[[j]] <- difs## / correction.factor
    }
    correction.factor <- averageSegsites/theta.denom/nfrags
    difs.mat <- matrix(unlist(difs.list), ncol=3, byrow = TRUE)
    difs.df <- as.data.frame(difs.mat)/correction.factor
    colnames(difs.df) <- apply(combn(3,2), 2, paste, collapse="-")
    return(difs.df)
}

stats <- matrix(nrow=nsims, ncol=1200)

stats.list <- foreach(i = 1:nsims) %dopar% {
    lim <- quantile(unlist(alldata[[i]]), probs=0.9)
    print(i)
    st <- list(a12=kde2d(alldata[[i]][,1], alldata[[i]][,2], n=20, lims=c(c(0,lim),c(0,lim)))$z,
               a13=kde2d(alldata[[i]][,1], alldata[[i]][,3], n=20, lims=c(c(0,lim),c(0,lim)))$z,
               a23=kde2d(alldata[[i]][,2], alldata[[i]][,3], n=20, lims=c(c(0,lim),c(0,lim)))$z)
    st.v <- unlist(st)
    return(st.v)
}
stats <- matrix(unlist(stats.list), nrow=nsims, byrow=TRUE)
classes <- c(rep(0, nsims/2), rep(1, nsims/2))
tols <- c(0.01, 0.05, 0.1, 0.2, 0.5)

cvABC <- foreach(i=1:length(tols)) %dopar% {
    cv <- cv4postpr(index=classes, sumstat=stats, nval=100, tols=tols[i], method='mnlogistic')
    return(cv)
}


lapply(cvABC, summary) # summary(cvABC[[5]])
save.image(file="multiple_CV.RData")


## save.image(file="tol_0.1_50k_sims.RData")



## save(alldata, file="sims_alldataTEST.RData")
## load("sims_alldataTEST.RData")

## alldata[[1]]

## labels <- array("SINGLE", nsims)
## labels[(nsims/2):nsims] <- "DOUBLE"

## library(MASS)
## lim <- quantile(unlist(alldata[[1]]), probs=0.9)
## obs <- list(a12=kde2d(alldata[[1]][,1], alldata[[1]][,2], n=20, lims=c(c(0,lim),c(0,lim)))$z,
##             a13=kde2d(alldata[[1]][,1], alldata[[1]][,3], n=20, lims=c(c(0,lim),c(0,lim)))$z,
##             a23=kde2d(alldata[[1]][,2], alldata[[1]][,3], n=20, lims=c(c(0,lim),c(0,lim)))$z)
## pdf("test.pdf")
## image(obs$a12)
## dev.off()

## stats.obs <- unlist(obs)

## params.original <- read.table("params.out", h=TRUE)
## params <- params.original[2:nsims,]
## stats <- matrix(nrow=nsims-1, ncol=1200)




## params.original[1,]
## library(abc)
## abc.res <- abc(target=obs, param=params, sumstat=stats, tol=0.9, method="loclinear", hcorr=T,
##                transf="logit", logit.bounds=matrix(c(range(tj23),range(tj123),
##                                range(tm23),
##                                range(tm23stop),
##                                range(tm12),
##                                range(tm12stop)),
##                              ncol=2, byrow=TRUE
##                              )
##                )
