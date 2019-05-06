require("Hmisc")
require("plyr")
require("ggplot2")
require("car")
require("compositions")
require("MatchIt")
require("cem")
require("Zelig")
library("Matching")
library("fastDummies")
options(digits=10)
options(scipen=10)

##############################################################
## Matching (Beta)
##############################################################

rm(list=ls())
load("prematch_data_all.RData")
data <- data.prematch

# Set treatment groups
data$TREATB <- NA
data[data$BETAQ==5,"TREATB"] <- 1
data[data$BETAQ==1,"TREATB"] <- 0

# Take log of cap and volume
data$LOGCAP <- log(data$CAP)
data$LOGVOL <- log(data$VOL)


# Make beta related data
data.beta <- data[,c("TREATB","PERIOD","RETURN","CAP","LOGCAP","VOL","LOGVOL","GROUP","INDUSTRY")]
data.beta <- data.beta[complete.cases(data.beta),]
data.beta <- data.beta[data.beta$CAP>0,]
data.beta <- data.beta[data.beta$VOL>0,]

# Coarsen the groups using SICC codes (https://www.naics.com/sic-codes-industry-drilldown/)
group_coarse <- rep(0, nrow(data.beta))
group_col <- data.beta$GROUP
for (i in (1:length(group_col))) {
  if (group_col[i] >= 1 & group_col[i] <= 9) {
    group_coarse[i] <- 1
  }
  if (group_col[i] >= 10 & group_col[i] <= 14) {
    group_coarse[i] <- 2
  }
  if (group_col[i] >= 15 & group_col[i] <= 17) {
    group_coarse[i] <- 3
  }
  if (group_col[i] >= 20 & group_col[i] <= 39) {
    group_coarse[i] <- 4
  }
  if (group_col[i] >= 40 & group_col[i] <= 49) {
    group_coarse[i] <- 5
  }
  if (group_col[i] == 50 | group_col[i] == 51) {
    group_coarse[i] <- 6
  }
  if (group_col[i] >= 52 & group_col[i] <= 59) {
    group_coarse[i] <- 7
  }
  if (group_col[i] >= 60 & group_col[i] <= 67) {
    group_coarse[i] <- 8
  }
  if (group_col[i] >= 70 & group_col[i] <= 89) {
    group_coarse[i] <- 9
  }
  if (group_col[i] >= 90 & group_col[i] <= 99) {
    group_coarse[i] <- 10
  }
}

# Remove all stocks with group == 0, add dummy variables for the rest
data.beta$GROUP_COARSE <- group_coarse
data.beta <- data.beta[data.beta$GROUP_COARSE != 0,]
data.beta <- dummy_cols(data.beta, select_columns = 'GROUP_COARSE')

# We chose to only replicate from period 452 to 552 for computational limitations
start <- 452
end <- 552

# Data frame to store data set.
data.prunedb <- data.frame()

for (i in start:end) {
  print(i)
  # Subset the data for period i
  data.curb <- data.beta[data.beta$PERIOD==i,]
  
  # Select the covariates, perform gen match and match
  Xb <-cbind(data.curb$LOGCAP, data.curb$LOGVOL, data.curb$GROUP_COARSE_5,
             data.curb$GROUP_COARSE_8, data.curb$GROUP_COARSE_2, data.curb$GROUP_COARSE_4,
             data.curb$GROUP_COARSE_9, data.curb$GROUP_COARSE_6, data.curb$GROUP_COARSE_7,
             data.curb$GROUP_COARSE_3, data.curb$GROUP_COARSE_10)
  genoutb <- GenMatch(Tr=data.curb$TREATB, X=Xb, M=1, pop.size=200,
                      max.generations=10, wait.generations=5, caliper = 0.1, print.level = 0)
  moutb <- Match(Tr=data.curb$TREATB, X=Xb, Weight.matrix=genoutb, caliper = 0.1)
  
  # Select matched data set, rbind to the bigger data set
  matched_treated <- cbind(data.curb[moutb$index.treated,], weights = moutb$weights)
  matched_control <- cbind(data.curb[moutb$index.control,], weights = moutb$weights)
  matched_data <- rbind(matched_treated, matched_control)
  
  data.prunedb <- rbind(data.prunedb, matched_data)
  
  # We are going to use the L1 distance to check the balance here instead of MatchBalance
  # in order to be comparable to the original paper.
}

# We had to split the data in half and run gen match for each half in our own machine

#load("SHO_temp.RData")
#load("HUNG_temp.RData")

# Combine the data
#data.prunedb <- rbind(data.prunedb.Hung, data.prunedb_sho)
#save(data.beta, data.prunedb, file="genmatch_data.RData")


# This is the data after the aggregation, just load and use.
load("genmatch_data.RData")


# Process the data just like the papaer did
data.prunedb$WEIGHT <- data.prunedb$CAP * data.prunedb$weights 
ret.beta <- ddply(data.prunedb,.(PERIOD,TREATB), summarise, RETURN=weighted.mean(RETURN,WEIGHT))


# Cumulative returns
cum.beta0 <- cumprod(ret.beta[ret.beta$TREATB==0,"RETURN"]+1)
cum.beta1 <- cumprod(ret.beta[ret.beta$TREATB==1,"RETURN"]+1)
cum.beta0[101]
cum.beta1[101]


##############################################################
## FIGURE 1
##############################################################

load("results_beta_all.RData")
load("genmatch_data.RData")

# Change the range of the original data
new_beta <- beta[which(beta$PERIOD >= 452),]

# Matched Beta plot
# Changed the range of the dates so they match our data
d <- data.frame(YEAR=1963+451/12+(1:101)/12,LINE=1,CUMRET=cum.beta0)
d <- rbind(d,data.frame(YEAR=1963+451/12+(1:101)/12,LINE=2,CUMRET=cum.beta1))
d <- rbind(d,data.frame(YEAR=1963+451/12+(1:101)/12,LINE=3,CUMRET=new_beta[new_beta$BETAQ==1,"CUMRET"]))
d <- rbind(d,data.frame(YEAR=1963+451/12+(1:101)/12,LINE=4,CUMRET=new_beta[new_beta$BETAQ==5,"CUMRET"]))
d[d$LINE==1,"Beta Quintiles"] <- paste("Matched Quintile 1 ($",format(d[d$LINE==1 & d$YEAR==2009,"CUMRET"],digits=2,nsmall=2),")",sep="")
d[d$LINE==2,"Beta Quintiles"] <- paste("Matched Quintile 5 ($",format(d[d$LINE==2 & d$YEAR==2009,"CUMRET"],digits=2,nsmall=2),")",sep="")
d[d$LINE==3,"Beta Quintiles"] <- paste("Original Quintile 1 ($",format(d[d$LINE==3 & d$YEAR==2009,"CUMRET"],digits=2,nsmall=2),")",sep="")
d[d$LINE==4,"Beta Quintiles"] <- paste("Original Quintile 5 ($",format(d[d$LINE==4 & d$YEAR==2009,"CUMRET"],digits=2,nsmall=2),")",sep="")
g <- ggplot(data=d,aes(x=YEAR,y=CUMRET,colour=`Beta Quintiles`,linetype=`Beta Quintiles`)) + geom_line(lwd=1)
g <- g + scale_y_log10("Cumulative Return",limits = c(.1,150),breaks=c(.1,1,10,100),labels=c("$0.10","$1","$10","$100"))
g <- g + scale_x_continuous("",limits = c(2000,2009),breaks=seq(2000,2008,1),seq(2000,2008,1),expand=c(0,0))
g <- g + ggtitle("All Stocks, Beta Quintiles\nCumulative Return of $1 invested in 1968\n")
g <- g + scale_color_manual(values=c("red","blue","red","blue"))
g <- g + scale_linetype_manual(values=c(1,1,3,3))
g <- g + theme(legend.position=c(.2,.8),panel.background=element_blank(),axis.line=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank(),plot.background=element_rect(fill=NA,colour =NA))
g
#dev.copy(pdf,"matchplotbeta.pdf")
#dev.copy(png,"matchplotbeta.png")
#dev.off()


L1 <- data.frame(period=NA,prebeta=NA,postbeta=NA)
N <- data.frame(period=NA,prebeta=NA,postbeta=NA)
for(i in 452:552) {
  print(i)
  prebeta <- imbalance(group=data.beta[data.beta$PERIOD==i,"TREATB"],data=data.beta[data.beta$PERIOD==i,c("LOGCAP","LOGVOL","GROUP_COARSE_5", "GROUP_COARSE_8","GROUP_COARSE_2", "GROUP_COARSE_4", "GROUP_COARSE_9","GROUP_COARSE_6", "GROUP_COARSE_7", "GROUP_COARSE_3","GROUP_COARSE_10")])$L1$L1
  postbeta <- imbalance(group=data.prunedb[data.prunedb$PERIOD==i,"TREATB"],
                        data=data.prunedb[data.prunedb$PERIOD==i,c("LOGCAP","LOGVOL","GROUP_COARSE_5", "GROUP_COARSE_8","GROUP_COARSE_2", "GROUP_COARSE_4", "GROUP_COARSE_9","GROUP_COARSE_6", "GROUP_COARSE_7", "GROUP_COARSE_3","GROUP_COARSE_10")], weights=data.prunedb[data.prunedb$PERIOD==i,"weights"])$L1$L1
  L1 <- rbind(L1,c(i,prebeta,postbeta))
  N <- rbind(N,c(i, nrow(data.beta[data.beta$PERIOD==i,]),nrow(data.prunedb[data.prunedb$PERIOD==i,])))
}


L1 <- L1[-1,]
N <- N[-1,]

apply(L1,2,mean)[-1]
apply(L1,2,sd)[-1]

load('new_l1.RData')

# Histogram of L1 statistic for Genetic Matching
hist(L1$postbeta, main = 'Histogram of L1 statistic for each period - GenMatch', xlab = 'L1 Statistic')


load('matching_imbalance.RData')

# Histogram of L1 statistic for original code
hist(L1$postbeta, main = 'Histogram of L1 statistic for each period - CEM', xlab = 'L1 Statistic')

