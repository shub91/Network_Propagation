library(wfg)
library(glmnet)


  np <- 16
  p <- 128
  n <- 80
  ###############################################
  
  g <- network.simu(p.in = rep(0.5, 4), p.out = 0.01)$net
  V(g)$size <- 5
  cl.vec <- rep("dodgerblue", 128)
  cl.vec[1:np] <- "red"
  # x11(); plot(g, vertex.label = "", vertex.color = cl.vec)
  
  ###############################################
  
  z <- rnorm(n)
  x <- matrix(rnorm(n*p), nrow = n, ncol = p)
  for (i in 1:np) { ###
    x[, i] <- z + rnorm(n, mean = 0, sd = 1)
  }
  y <- rowSums(x[,1:np])*3/sqrt(np) + rnorm(n, mean = 0, sd = 2)
  
  ## Obtain initial estimate ########################
  alpha <- 1 # 1 for lasso
  a <- 0.9
  cv.fit <- cv.glmnet(x, y, alpha = alpha)
  
  # x11(); plot(cv.fit)
  lbd <- cv.fit$lambda.min
  glm.fit <- glmnet(x, y, lambda = lbd, alpha = alpha)
  beta.1 <- glm.fit$beta[,1]
  Y <- abs(beta.1)
  names(Y) <- NULL
  
  ## Propagation ########################
  A <- as.matrix(get.adjacency(g))
  W <- A
  diag(W) <- 1
  deg.vec <- colSums(W)
  
  n <- nrow(W)
  D <- matrix(0, nrow = n, ncol = n)
  D_sqr <- D
  diag(D_sqr) <- 1/sqrt(deg.vec)
  W_prime <- D_sqr %*% W %*% D_sqr

  F_old = Y
  F_store <- c(F_old)
  alph = .5
  for (i in 1:1000){
       F_new = alph*W_prime %*% F_old + (1-alph)*Y
       F_store  = cbind(F_store, F_new)
       F_old <- F_new
  }
  
  hub_weights <- F_store[, i]
  V(g)$weight <- hub_weights
  c_scale <- colorRamp(c('white', 'red')) #Color scaling function
  V(g)$color = apply(c_scale(V(g)$weight), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
  quartz()
  plot(g)

  Y = F_store[, i]


quartz()
plot(F_store[1,])



































  ## Estimates from adaptive lasso
  cv.fit.2 <- cv.glmnet(x, y, penalty.factor = 1/Y, alpha = alpha)
  lbd.2 <- cv.fit.2$lambda.min
  glm.fit.2 <- glmnet(x, y, lambda = lbd.2, penalty.factor = 1/Y, alpha = alpha)
  beta.2 <- glm.fit.2$beta[,1]
  
  ## Estimates from New method
  cv.fit.3 <- cv.glmnet(x, y, penalty.factor = 1/Ft, alpha = alpha)
  lbd.3 <- cv.fit.3$lambda.min
  glm.fit.3 <- glmnet(x, y, lambda = lbd.3, penalty.factor = 1/Ft, alpha = alpha)
  beta.3 <- glm.fit.3$beta[,1]
  
  s0 <- c(1:np)
  s1 <- which(beta.1 != 0)
  s2 <- which(beta.2 != 0)
  s3 <- which(beta.3 != 0)
  
  sen.1 <- length(intersect(s0, s1)) / length(s0)
  sen.2 <- length(intersect(s0, s2)) / length(s0)
  sen.3 <- length(intersect(s0, s3)) / length(s0)
  
  
  s0c <- setdiff(1:p, s0)
  s1c <- setdiff(1:p, s1)
  s2c <- setdiff(1:p, s2)
  s3c <- setdiff(1:p, s3)
  
  spe.1 <- length(intersect(s0c, s1c)) / length(s0c)
  spe.2 <- length(intersect(s0c, s2c)) / length(s0c)
  spe.3 <- length(intersect(s0c, s3c)) / length(s0c)
  
  beta.0 <- c(rep(5/sqrt(np), np), rep(0, 128-np))
  result <- list()
  result$fs <- data.frame(fnr = 1 - c(sen.1, sen.2, sen.3), fpr = 1 - c(spe.1, spe.2, spe.3))
  result$l2 <- c( mean((beta.1 - beta.0)^2), mean((beta.2 - beta.0)^2), mean((beta.3 - beta.0)^2) )
  
  return(result)
}







