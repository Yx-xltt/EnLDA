library(quadprog)

linprogPD <- function(x0, A, b, epsilon, pdtol=1e-3, pdmaxiter=50) {
  
  N <- length(x0)
  x0 <- matrix(x0, nrow=N)
  b <- matrix(b, nrow=N)
  
  alpha = 0.01
  beta = 0.5
  mu = 10
  
  gradf0 <- matrix(rep(c(0,1), each=N), nrow=2*N)
  
  if (max(abs(A%*%x0 - b)) > epsilon) {
    stop("Infeasible starting point!")
  }
  
  x <- x0
  u <- 0.95*abs(x0) + 0.1*max(abs(x0))
  
  Atr <- A%*%x - b
  
  fu1 <- x - u
  fu2 <- -x - u
  fe1 <- Atr - epsilon
  fe2 <- -Atr - epsilon
  lamu1 <- -1/fu1
  lamu2 <- -1/fu2
  lame1 <- -1/fe1
  lame2 <- -1/fe2
  
  AtAv <- t(A)%*%(lame1 - lame2)
  
  sdg <- - sum(c(fu1, fu2, fe1, fe2)* c(lamu1, lamu2, lame1, lame2))
  
  tau <- mu*(4*N)/sdg
  
  rdual <- gradf0 + c(lamu1 - lamu2 + AtAv, -lamu1 - lamu2)
  rcent <- -c(lamu1*fu1, lamu2*fu2, lame1*fe1, lame2*fe2) - 1/tau
  resnorm <- sqrt(sum(rdual^2,rcent^2))
  
  pditer <- 0
  done <- (sdg < pdtol) | (pditer >= pdmaxiter)
  
  while(!done) {
    w2 <- -1 - (1/fu1 + 1/fu2)/tau
    
    sig11 <- -lamu1/fu1 - lamu2/fu2
    sig12 <- lamu1/fu1 - lamu2/fu2
    siga <- -(lame1/fe1 + lame2/fe2)
    sigx <- sig11 - sig12^2/sig11
    
    w1 <- -( t(A)%*%(1/fe2 - 1/fe1) + 1/fu2 - 1/fu1)/tau
    w1p <- w1 - (sig12/sig11)*w2
    
    Hp <- t(A)%*%diag(as.vector(siga))%*%A + diag(as.vector(sigx))
    
    dx <- solve(Hp, w1p)
    
    if (rcond(Hp) < 1e-14) {
      warning("Ill conditioned matrix.  Previous iterate matrix returned! (May increase perturb/lambda.)")
      xp <- x
      return(xp)
    }
    AtAdx <- A%*%dx
    
    du <- w2/sig11 - (sig12/sig11)*dx
    
    dlamu1 <- -(lamu1/fu1)*(dx-du) - lamu1 - 1/(fu1*tau)
    dlamu2 <- -(lamu2/fu2)*(-dx-du) - lamu2 - 1/(fu2*tau)
    
    dlame1 <- -(lame1/fe1)*(AtAdx) - lame1 - 1/(fe1*tau)
    dlame2 <- (lame2/fe2)*(AtAdx) - lame2 - 1/(fe2*tau)
    
    AtAdv <- t(A)%*%(dlame1 - dlame2)
    
    
    iu1 <- dlamu1 < 0
    iu2 <- dlamu2 < 0
    ie1 <- dlame1 < 0
    ie2 <- dlame2 < 0
    ifu1 <- (dx-du) > 0
    ifu2 <- (-dx-du) > 0
    ife1 <- AtAdx > 0
    ife2 <- AtAdx < 0
    
    smax <- min( -lamu1[iu1]/dlamu1[iu1], -lamu2[iu2]/dlamu2[iu2], -lame1[ie1]/dlame1[ie1], -lame2[ie2]/dlame2[ie2], -fu1[ifu1]/(dx[ifu1] - du[ifu1]), -fu2[ifu2]/(-dx[ifu2] -du[ifu2]), -fe1[ife1]/AtAdx[ife1], -fe2[ife2]/( - AtAdx[ife2])   )
    smax <- min(1,smax)
    
    s <- 0.99*smax
    
    suffdec <- FALSE
    backiter <- 0
    
    while(!suffdec) {
      xp <- x + s*dx
      up <- u + s*du
      Atrp <- Atr + s*AtAdx
      AtAvp <- AtAv+s*AtAdv
      
      fu1p <- fu1 + s*(dx - du)
      fu2p <- fu2 + s*(-dx-du)
      fe1p <- fe1 + s*AtAdx
      fe2p <- fe2 + s*(-AtAdx)
      
      lamu1p <- lamu1 + s*dlamu1
      lamu2p <- lamu2 + s*dlamu2
      lame1p <- lame1 + s*dlame1
      lame2p <- lame2 + s*dlame2
      
      rdp <- gradf0 + c(lamu1p - lamu2p + AtAvp, -lamu1p - lamu2p)
      rcp <- -c(lamu1p*fu1p, lamu2p*fu2p, lame1p*fe1p, lame2p*fe2p) - 1/tau
      suffdec <- sqrt( sum(rdp^2, rcp^2)) < (1 - alpha*s)*resnorm
      
      s <- beta*s
      backiter <- backiter+1
      
      if (backiter > 32) {
        warning("Backtracking stuck.  Previous iterate matrix returned!")
        xp <- x
        return(xp)
      }
    }
    x <- xp
    u <- up
    Atr <- Atrp
    AtAv <- AtAvp
    fu1 <- fu1p
    fu2 <- fu2p
    fe1 <- fe1p
    fe2 <- fe2p
    lamu1 <- lamu1p
    lamu2 <- lamu2p
    lame1 <- lame1p
    lame2 <- lame2p
    
    sdg <- -sum( c(fu1, fu2, fe1, fe2)*c(lamu1, lamu2, lame1, lame2))
    tau <- mu*(4*N)/sdg
    rdual <- rdp
    rcent <- -c(fu1, fu2, fe1, fe2)*c(lamu1, lamu2, lame1, lame2) - 1/tau
    resnorm <- sqrt(sum(rdual^2, rcent^2))
    
    pditer <- pditer+1
    done <- (sdg < pdtol) | (pditer >= pdmaxiter)
  }
  return(xp)
}
compute_intercept <- function(data, labels, beta) {
  labels <- as.factor(labels)
  
  mu1 <- colMeans(data[labels == 1, ])
  mu2 <- colMeans(data[labels == 0, ])
  
  intercept <- as.numeric(-0.5 * t(mu2 + mu1) %*% beta)
  
  return(intercept)
}
slda_score <- function(data, beta, intercept, n1, n2) {
  scores <- data %*% beta + intercept + log(n1 / n2)
  return(scores)
}
slda_classify <- function(data, beta, intercept, n1, n2) {
  scores <- slda_score(data, beta, intercept, n1, n2)
  predicted_labels <- ifelse(scores > 0, 1, 0)
  return(predicted_labels)
}
logit_to_prob <- function(scores) {
  p <- 1 / (1 + exp(-scores))
  return(p)
}

all_merged <- read.csv("Data.csv", stringsAsFactors = FALSE)

feature_cols <- grep("^X[0-9]+$", colnames(all_merged), value = TRUE)

#train
idx_1 <- which(all_merged$Y == 1)
idx_0 <- which(all_merged$Y == 0)

train_idx <- c(sample(idx_1, 80), sample(idx_0, 80))
train_df_final <- all_merged[train_idx, ]

remaining_idx <- setdiff(1:nrow(all_merged), train_idx)
rem_idx_1 <- remaining_idx[all_merged$Y[remaining_idx] == 1]
rem_idx_0 <- remaining_idx[all_merged$Y[remaining_idx] == 0]

test_idx <- c(sample(rem_idx_1, 80), sample(rem_idx_0, 80))
test_df_final <- all_merged[test_idx, ]

train_data <- as.matrix(train_df_final[, feature_cols])
b_train <- train_df_final$Y

test_data <- as.matrix(test_df_final[, feature_cols])
b_test <- test_df_final$Y

n <- nrow(train_data)
p <- ncol(train_data)

n1 <- sum(b_train == 1)
n2 <- sum(b_train == 0)

delta <- colMeans(train_data[b_train==1,]) - colMeans(train_data[b_train==0,])
Sigma <- cov(train_data)
eigvals = eigen(Sigma, only.values=T)$values
perturb = max(max(eigvals) - p*min(eigvals), 0)/(p-1)
sigma_hat_star = Sigma + diag(p)*perturb
beta_init = solve(sigma_hat_star) %*% delta

lambda_grid <- seq(0.01, 0.8, length.out = 50)

beta_list <- list()
intercept_list <- numeric(length(lambda_grid))
test_errors <- numeric(length(lambda_grid))
scores_test <- list()
pm_test <- list()  

m <- length(lambda_grid)  
H <- matrix(0, nrow = 2 * n, ncol = m)

for (i in 1:length(lambda_grid)) {
  lambda <- lambda_grid[i]
  beta <- linprogPD(beta_init, sigma_hat_star, delta, lambda)
  threshold <- 1e-3
  beta[abs(beta) < threshold] <- 0  
  
  beta_list[[i]] <- beta
  intercept <- compute_intercept(train_data, b_train, beta)
  intercept_list[i] <- intercept
  scores_test[[i]] <- slda_score(train_data, beta, intercept, n1, n2)
  pm_test[[i]] <- logit_to_prob(scores_test[[i]])
  predicted_labels <- slda_classify(train_data, beta, intercept, n1, n2)
  test_errors[i] <- mean(predicted_labels != b_train)
  
  H[seq(1, 2 * n, by = 2), i] <- pm_test[[i]]  
  H[seq(2, 2 * n, by = 2), i] <- 1 - pm_test[[i]]
}

Y <- rep(0, 2 * n)  
for (i in 1:n) {
  if (b_train[i] == 1) {
    Y[2 * i - 1] <- 1  
  } else {
    Y[2 * i] <- 1  
  }
}

epsilon <- 1e-6  
Dmat <- t(H) %*% H + epsilon * diag(m)
dvec <- t(H) %*% Y 
Amat <- rbind(rep(1, m),diag(m))  
bvec <- c(1,rep(0, m)) 
meq <- 1 
weights <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat), bvec = bvec, meq = 1)
optimal_weights <- weights$solution
print(optimal_weights)

threshold <- 1e-6  
optimal_weights[abs(optimal_weights) < threshold] <- 0  
print(optimal_weights)
  
#test
n1_test <- sum(b_test == 1)
n2_test <- sum(b_test == 0)
  
n_test <- n1_test + n2_test
weighted_scores <- rep(0, n_test)
for (i in 1:length(lambda_grid)) {
  scores <- test_data %*% beta_list[[i]] + intercept_list[i] + log(n1_test / n2_test)
  weighted_scores <- weighted_scores + optimal_weights[i] * scores
}
predicted_labels <- ifelse(weighted_scores > 0, 1, 0)
mis_rate <- mean(predicted_labels != b_test)
cat("Test Misclassification Rate:", mis_rate, "\n")
