random.schemes <- function(N, pi, S, profiles) {
  SRS_A = SRS(N, pi)
  SBR_A = SBR(S, pi, bsize = 6)
  SBCD_A = SBCD(S, pi)
  PS_A = PS(profiles, rep(1/ncol(profiles), ncol(profiles)), pi, lambda = 0.75)
  return(list(SRS_A=SRS_A, SBR_A=SBR_A, SBCD_A=SBCD_A, PS_A = PS_A))
}

all_outputs <- function(pi, sample_data) {
  S <- sample_data[["S"]]
  X <- sample_data[["X"]]
  
  SRS_Y = sample_data[["SRS_Y"]]
  SRS_A = sample_data[["SRS_A"]]
  SRS_Y1_fit = sample_data[["SRS_Y1_fit"]]
  SRS_Y0_fit = sample_data[["SRS_Y0_fit"]]
  
  SBR_Y = sample_data[["SBR_Y"]]
  SBR_A = sample_data[["SBR_A"]]
  SBR_Y1_fit = sample_data[["SBR_Y1_fit"]]
  SBR_Y0_fit = sample_data[["SBR_Y0_fit"]]
  
  SBCD_Y = sample_data[["SBCD_Y"]]
  SBCD_A = sample_data[["SBCD_A"]]
  SBCD_Y1_fit = sample_data[["SBCD_Y1_fit"]]
  SBCD_Y0_fit = sample_data[["SBCD_Y0_fit"]]
  
  PS_Y = sample_data[["PS_Y"]]
  PS_A = sample_data[["PS_A"]]
  PS_Y1_fit = sample_data[["PS_Y1_fit"]]
  PS_Y0_fit = sample_data[["PS_Y0_fit"]]
  
  outputs = c(
  output(SRS_Y, SRS_A, S, X, pi, pi*(1-pi), SRS_Y1_fit, SRS_Y0_fit),
  output(SBR_Y, SBR_A, S, X, pi, 0, SBR_Y1_fit, SBR_Y0_fit),
  output(SBCD_Y, SBCD_A, S, X, pi, 0, SBCD_Y1_fit, SBCD_Y0_fit),
  output(PS_Y, PS_A, S, X, pi, 0, PS_Y1_fit, PS_Y0_fit)
  )
  
  return(outputs)
}

output <- function(Y, A, S, X, pi, q, Y1_fit, Y0_fit) {
  c(
    ols.test.cont(Y, A, X),
    usual.test.cont(Y, A, X), 
    mod.test.cont(Y, A, S, X, pi, q), 
    eff.test.cont(Y, A, S, X, pi, Y1_fit, Y0_fit)
  )
}

#### cross-fitting
Ya.fit <- function(Y, A, S, X, Z) {
    cross.fit_strata(Y, A, S, X, Z, M=5)
}


I.know.fit <- function(Y, A, X, Z) {
  Y1_fit = 5 + exp((0.5 + 0)*X) + 2*Z + 6 * Z *X * A
  Y0_fit = 4 + exp((0.5 + 0)*X) + 2*Z 
  Y_fit = A * Y1_fit + (1-A) * Y0_fit
  res = Y - Y_fit
  return(data.frame(Y_fit = Y_fit,res = res, Y1_fit = Y1_fit, Y0_fit = Y0_fit ))
  
}

create.fold <- function(n, M = 2) {
  n_fold = numeric(M)
  fold_seq = c()
  for (m in 1:(M-1)) {
    n_fold[m] = floor(n / M)
    fold_seq = c(fold_seq, rep(m, n_fold[m]))
  }
  n_fold[M] = n - sum(n_fold[1:(M-1)])
  fold_seq = c(fold_seq, rep(M, n_fold[M]))
  fold = sample(fold_seq, n)
  return(fold)
}


cross.fit <- function(Y, A, X, Z, M=2) {
  n = length(Y)
  fold = create.fold(n, M)
  
  total_data = data.frame(Y = Y, A = A, X = X, Z = Z)
  
  Y1_fit = numeric(n)
  Y0_fit = numeric(n)
  Y_fit = numeric(n)
  res = numeric(n)
  
  for (m in 1:M) {
    index_in_fold = which(fold == m)
    index_out_fold = which(fold != m)
    train_data = total_data[index_out_fold,]
    train_data_1 = train_data %>% subset(A==1)
    train_data_0 = train_data %>% subset(A==0)
    
    bws1 = c(0.1256611,  0.1750254)
    bws0 = c(0.3156482, 0.09994133)
    
    kernel_model_1 = npregbw(Y ~ X + Z, bandwidth.compute = F, 
                             bws = bws1, data = train_data_1, regtype = "ll")
    kernel_model_0 = npregbw(Y ~ X + Z, bandwidth.compute = F, 
                             bws = bws0, data = train_data_0, regtype = "ll")
    Y1_fit[index_in_fold] = npreg(kernel_model_1, exdat = cbind(X, Z)[index_in_fold,])$mean
    Y0_fit[index_in_fold] = npreg(kernel_model_0, exdat = cbind(X, Z)[index_in_fold,])$mean
    Y_fit[index_in_fold] = A[index_in_fold] * Y1_fit[index_in_fold] + (1-A)[index_in_fold] * Y0_fit[index_in_fold]
    res[index_in_fold] = Y[index_in_fold] - Y_fit[index_in_fold]
  }
  return(data.frame(Y_fit = Y_fit,res = res, Y1_fit = Y1_fit, Y0_fit = Y0_fit ))
}

cross.fit_strata <- function(Y, A, S, X, Z, M=2) {
  n = length(Y)
  nStrata = max(S)
  fold = create.fold(n, M)
  
  total_data = data.frame(Y = Y, A = A, S = S, X = X, Z = Z)
  
  Y1_fit = numeric(n)
  Y0_fit = numeric(n)
  Y_fit = numeric(n)
  res = numeric(n)
  
  for (m in 1:M) {
    for (i in 1:nStrata) {
    index_in_fold = which(fold == m & S == i)
    index_out_fold = which(fold != m & S == i)
    train_data = total_data[index_out_fold,]
    train_data_1 = train_data %>% subset(A==1)
    train_data_0 = train_data %>% subset(A==0)
    
    bws1 = c(0.1256611,  0.1750254)
    bws0 = c(0.3156482, 0.09994133)
    
    kernel_model_1 = npregbw(Y ~ X + Z, bandwidth.compute = F, 
                             bws = bws1, data = train_data_1, regtype = "ll")
    kernel_model_0 = npregbw(Y ~ X + Z, bandwidth.compute = F, 
                             bws = bws0, data = train_data_0, regtype = "ll")
    Y1_fit[index_in_fold] = npreg(kernel_model_1, exdat = cbind(X, Z)[index_in_fold,])$mean
    Y0_fit[index_in_fold] = npreg(kernel_model_0, exdat = cbind(X, Z)[index_in_fold,])$mean
    Y_fit[index_in_fold] = A[index_in_fold] * Y1_fit[index_in_fold] + (1-A)[index_in_fold] * Y0_fit[index_in_fold]
    res[index_in_fold] = Y[index_in_fold] - Y_fit[index_in_fold]
    }
  }
  return(data.frame(Y_fit = Y_fit,res = res, Y1_fit = Y1_fit, Y0_fit = Y0_fit ))
}