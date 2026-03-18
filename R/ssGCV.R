ssGCV <- function (X, y, mod_vec, target, lambda, alpha, bandwidth) {
  weight  <- exp(-(mod_vec - mod_vec[target])^2 / bandwidth)
  sw <- sqrt(weight)
  
  K.X <- X * sw
  K.y <- y * sw
  
  fit <- glmnet(K.X, K.y, alpha=alpha, lambda=lambda)
  prediction <- predict(fit, newx = K.X, s = lambda)
  beta <- fit$beta
  
  sigma = lambda * (1 - alpha + alpha / pmax(abs(beta), 1e-12))
  A <- crossprod(K.X) 
  diag(A) <- diag(A) + sigma 
  H <- K.X %*% solve(A, t(K.X)) 
  
  sum((K.y - prediction)^2) / (length(K.X) - sum(diag(H)))^2
}
