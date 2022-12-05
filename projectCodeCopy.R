#-------------------------------
#-------------------------------
# Sim in R - Final Project
# Simple HLM and Monte-Carlo
# Akshay Umashankar
# au3692
#-------------------------------
## Libraries ##
library(lme4)
library(rblimp)
#test
#-------------------------------
# source("scripts/rmvnorm.R")
# source("scripts/ols_regression.R")
# source("scripts/generate_reg.R")


#Lab 3 Reference
# set.seed(88888)
# n = 100
# p_x <- list(
#   rho = 0.3,
#   mux = c(0, 0),
#   sigmax = c(1, 1),
#   rhoxE = 0.3
# )
# 
# p_y <- list(
#   betaMat = c(1, 1),
#   beta0 = 100,
#   muy = 0,
#   r_squared = 0.4
# )
# 
# monteFunction <- function(n, p_x, p_y){
#   rmvSetupMonte <- function(muX, sigma, rho){
#     R_mat <- diag(NROW(diag(p_x$sigmax)))
#     R_mat[lower.tri(R_mat)] <- rho
#     R_mat[upper.tri(R_mat)] <- rho
#     sigmaX <- diag(p_x$sigmax) %*% R_mat %*% diag(p_x$sigmax)
#     setup <- list(
#       sigmaX = sigmaX,
#       muX = muX
#     )
#     return(setup)
#   }
#   
#   setupParam <- with(p_x, rmvSetupMonte(mux, sigmax, rho))
#   varE <- with(p_y, t(betaMat) %*% setupParam$sigmaX %*% betaMat %*% (1/r_squared - 1))
#   sigmaE <- sqrt(varE)
#   
#   corMat <- rbind(c(p_x$rhoxE, rep(0, NROW(setupParam$sigmaX) - 1)), setupParam$sigmaX)
#   corMat <- cbind(c(1, p_x$rhoxE, rep(0, NCOL(setupParam$sigmaX) - 1)), corMat)
#   
#   sigmaEMat <- c(sigmaE, rep(1, NCOL(corMat) - 1))
#   
#   covMat <- diag(sigmaEMat) %*% corMat %*% diag(sigmaEMat)
#   
#   draw <- rmvnorm(n, c(p_x$mux, 0), covMat) #DRAWS X1 X2 and E
#   
#   beta0 <- p_y$beta0
#   Y_method1 <- beta0 + p_y$betaMat[1] * draw[, 1] + p_y$betaMat[2] * draw[, 2] + draw[, NCOL(draw)]
#   
#   output <- data.frame(
#     X = draw[, -1],
#     Y = Y_method1
#     #weight = weights
#   )
# }
# 
# montHelp <- monteFunction(n, p_x, p_y)

#------------------------------------------
#MSE and Jackknife Functions
factoryMSE <- function(pop_value){
  force(pop_value)
  function(parameter){
    rowMeans((parameter - pop_value)^2)
  }
}

jackknife <- function(data, func){
  i <- 1
  biasVec <- c()
  for (x in 1:length(data[1, 1, ])){
    dataSub <- data[-i]
    biasVec <- c(biasVec, mean(sapply(dataSub, func)))
    i <- i + 1
  }
  thetaJack <- mean(biasVec)
  varJack <- ((length(data)-1)/length(data)) * sum((biasVec - thetaJack)^2)
  sdJack <- sqrt(varJack)
  return(sdJack)
}
#------------------------------------------

############# Data Generation ################
coefUnconditional <- c()
coefConditional <- c()

# Conditions
parameters1 <- list(
  J    = 1000, #1000 L2 Clusters
  nj   = 100, #100 individuals per Cluster
  ICC  = 0.10,
  varianceY = 100, #totalVariance
  mean = 25,
  d = 0.1
)

parameters2 <- list(
  J    = 1000, #1000 L2 Clusters
  nj   = 100, #100 individuals per Cluster
  ICC  = 0.25,
  varianceY = 100, #totalVariance
  mean = 25,
  d = 0.5
)

parameters3 <- list(
  J    = 1000, #1000 L2 Clusters
  nj   = 100, #100 individuals per Cluster
  ICC  = 0.10,
  varianceY = 100, #totalVariance
  mean = 25,
  d = 0.1
)

parameters4 <- list(
  J    = 1000, #1000 L2 Clusters
  nj   = 100, #100 individuals per Cluster
  ICC  = 0.25,
  varianceY = 100, #totalVariance
  mean = 25,
  d = 0.5
)

parameters <- list(
  J    = 1000, #1000 L2 Clusters
  nj   = 100, #100 individuals per Cluster
  ICC  = 0.25,
  varianceY = 100, #totalVariance
  mean = 25,
  d = 0.5
)
#----------------------------------------

dataGenerated <- function(parameters){
  l1_var <- with(parameters, varianceY*(1-ICC))
  l2_varTotal <- with(parameters, varianceY*ICC)
  gamma1 <- with(parameters, d*sqrt(varianceY))
  
  treat <- with(parameters, c(rep(1, J/2), rep(0, J/2))) #Set up 50 Individuals in treatment and 50 individuals No treatment for 1 cluster
  
  L2_id <- sort(with(parameters, rep(c(1:J), nj)))
  L1_id <- with(parameters, rep(c(1:nj), J))
  X <- cbind(L2_id, L1_id, treat = treat[L2_id])
  Y_within <- with(parameters, rnorm(J*nj, mean, l1_var))
  l2_varExplained <-  with(parameters, ICC*gamma1^2)
  l2_residual <- l2_varTotal - l2_varExplained
  Y_between <- with(parameters, treat * gamma1 + rnorm(J, 0, sqrt(l2_residual)))
  L2_Y <- Y_between %x% with(parameters, rep(1, nj))
  
  Y <- Y_within + L2_Y
  data <- data.frame(X, Y)
  #colnames(data, c("L2_id", "L1_id", "treat", "Y"))
  return(data)
}

#data_icc_d
data_1_1  <- dataGenerated(parameters1)
data_25_1 <- dataGenerated(parameters2)
data_1_5  <- dataGenerated(parameters3)
data_25_5 <- dataGenerated(parameters4)


############# Data Analysis ################

#Unconditional
uncondMonte <- lmer(Y ~ 1 + (1|L2_id), data)
model1 <- summary(uncondMonte)
coefUnconditional <- c(coefUnconditional, model1$coefficients)


#Conditional Model
condMonte <- lmer(Y ~ 1 + treat + (1|L2_id), data)
model2 <- summary(condMonte)
model2$coefficients
head(data)

############ Replications #################
#Create a list of datasets

dataArray1 <- as.list(replicate(dataGenerated(parameters1), n = 500, simplify = FALSE))
dataArray2 <- as.list(replicate(dataGenerated(parameters2), n = 500, simplify = FALSE))
dataArray3 <- as.list(replicate(dataGenerated(parameters3), n = 500, simplify = FALSE))
dataArray4 <- as.list(replicate(dataGenerated(parameters4), n = 500, simplify = FALSE))

#modelArray <- lapply(dataArray2, lmer(Y ~ 1 + treat + (1|L2_id), dataArray2))
#modTest <- lmer(Y ~ 1 + treat + (1|L2_id), dataArray2[[1]])

####ANALYSIS####
analyzeData <- function(dataArray, d){
  coeffVec <- numeric(length = NROW(dataArray))
  evalVec <- numeric(length = NROW(dataArray))
  for (i in 1:NROW(dataArray)){
    model  <- lmer(Y ~ 1 + treat + (1|L2_id), dataArray[[i]])
    modSum <- summary(model)
    upper  <- modSum$coefficient[2, 1] + 1.96 * modSum$coefficients[2, 2]
    lower  <- modSum$coefficient[2, 1] - 1.96 * modSum$coefficients[2, 2]
    test   <- d*sqrt(100) <= upper & d*sqrt(100) >= lower
    coeffVec[i] <- modSum$coefficient[2, 1]
    evalVec[i]  <- test
  }
  output <- list(
    coeffVec = coeffVec,
    evalVec = evalVec,
    evaluationPercent = sum(evalVec)/length(evalVec)
  )  
}

##### RESULTS ######
#MLM Cluster
test1 <- lmer(Y ~ 1 + treat + (1|L2_id), dataArray2[[1]])
summary(test1)
#Condition 1
analyze1 <- analyzeData(dataArray1, parameters1$d)
analyze1$evaluationPercent

#Condition 2
analyze2 <- analyzeData(dataArray2, parameters2$d)
analyze2$evaluationPercent

#Condition 3
analyze3 <- analyzeData(dataArray3, parameters3$d)
analyze3$evaluationPercent

#Condition 4
analyze4 <- analyzeData(dataArray4, parameters4$d)
analyze4$evaluationPercent


#Standard LM Model
test2 <- lm(Y ~ treat, dataArray1[[1]])
summary(test2)

analyzeDataLM <- function(dataArray, d){
  coeffVec <- numeric(length = NROW(dataArray))
  evalVec <- numeric(length = NROW(dataArray))
  for (i in 1:NROW(dataArray)){
    model  <- lm(Y ~ 1 + treat, dataArray[[i]])
    modSum <- summary(model)
    upper  <- modSum$coefficient[2, 1] + 1.96 * modSum$coefficients[2, 2]
    lower  <- modSum$coefficient[2, 1] - 1.96 * modSum$coefficients[2, 2]
    test   <- d*sqrt(100) <= upper & d*sqrt(100) >= lower
    coeffVec[i] <- modSum$coefficient[2, 1]
    evalVec[i]  <- test
  }
  output <- list(
    coeffVec = coeffVec,
    evalVec = evalVec,
    evaluationPercent = sum(evalVec)/length(evalVec)
  )  
}

analyze1LM <- analyzeDataLM(dataArray1, parameters1$d)
analyze1LM$evaluationPercent

#Condition 2
analyze2LM <- analyzeDataLM(dataArray2, parameters2$d)
analyze2$evaluationPercent

#Condition 3
analyze3LM <- analyzeDataLM(dataArray3, parameters3$d)
analyze3LM$evaluationPercent

#Condition 4
analyze4LM <- analyzeDataLM(dataArray4, parameters4$d)
analyze4LM$evaluationPercent

#Robust Standard Error
analyzeDataRob <- function(dataArray, d){
  coeffVec <- numeric(length = NROW(dataArray))
  evalVec <- numeric(length = NROW(dataArray))
  for (i in 1:NROW(dataArray)){
    model  <- lm(Y ~ 1 + treat, dataArray[[i]])
    modSum <- summary(model)
    tmp <- hlm::sandwich(model, dataArray[[i]][, 1])
    upper  <- tmp$coefficient[2, 1] + 1.96 * tmp$coefficients[2, 2]
    lower  <- tmp$coefficient[2, 1] - 1.96 * tmp$coefficients[2, 2]
    test   <- d*sqrt(100) <= upper & d*sqrt(100) >= lower
    coeffVec[i] <- tmp$coefficient[2, 1]
    evalVec[i]  <- test
  }
  output <- list(
    coeffVec = coeffVec,
    evalVec = evalVec,
    evaluationPercent = sum(evalVec)/length(evalVec)
  )  
}

#Condition 1
analyze1Rob <- analyzeDataRob(dataArray1, parameters1$d)
analyze1LM$evaluationPercent

#Condition 2
analyze2Rob <- analyzeDataRob(dataArray2, parameters2$d)
analyze2$evaluationPercent

#Condition 3
analyze3Rob <- analyzeDataRob(dataArray3, parameters3$d)
analyze3Rob$evaluationPercent

#Condition 4
analyze4Rob <- analyzeDataRob(dataArray4, parameters4$d)
analyze4Rob$evaluationPercent




