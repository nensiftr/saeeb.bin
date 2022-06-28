#' @export
saebebin <- function(y, n){
  result <- list(pi_EB = NA, Parameter = list(alpha = NA, beta = NA),
                 MSE.pi_EB = list(method = "Jackknife", mse = NA), direct =
                   list(est = NA, mse = NA))

  if (any(is.na(y)))
    stop("Argument y contains NA values.")
  if (any(is.na(n)))
    stop("Argument n contains NA values.")

  m <- length(y) #jumlah area
  nT <- sum(n) #jumlah seluruh sampel
  w <- n/nT #bobot

  #menghitung dugaan langsung proporsi dan variansnya
  p_i <- y/n #dugaan langsung proporsi
  mse.p_i <- p_i * (1-p_i) / n #varians dugaan langsung
  result$direct$est <- p_i
  result$direct$mse <- mse.p_i

  #menghitung proporsi dan ragam proporsi
  p_ib <- w * p_i
  p_hat <- sum(p_ib) #rataan terboboti
  s_p2 <- w * (p_i-p_hat)^2
  sum_s_p2 <- sum(s_p2) #ragam terboboti

  #menduga parameter sebaran beta-binomial alpha dan beta dengan metode momen Kleinman
  n2nT <- (n^2) / nT
  sum_n2nT <- sum(n2nT)
  k11 <- (nT*sum_s_p2) - p_hat*(1-p_hat) * (m-1)
  k12 <- p_hat*(1-p_hat) * (nT-sum_n2nT-(m-1))
  k1 <- k11/k12

  #nilai alpha beta
  alpha <- p_hat*(1-k1)/k1
  beta <- alpha*(1/p_hat-1)

  #pendugaan bayes dan ragam posterior bagi p_i
  k21 <- (y+alpha)*(n-y+beta)
  k22 <- (n+alpha+beta+1)*(n+alpha+beta)^2
  p_i_hat_EB1 <- (y+alpha)/(n+alpha+beta) #penduga bayes
  var_p_i_hat_EB <- k21/k22 #ragam posterior

  #pendugaan empirical bayes bagi p_i
  gamma <- n/(n+alpha+beta) #nilai gamma
  pi_EB = gamma*p_i+(1-gamma)*p_hat #penduga empirical bayes
  result$Parameter$alpha <- alpha
  result$Parameter$beta <- beta
  result$pi_EB <- pi_EB

  #Menghitung MSE penduga empirical bayes
  jackknife <- function(y, n, l){
    p_i_jk <- 0
    n.jk <- 0
    p_hat_jk <- 0
    p_i_jk <- y[-l]/n[-l] #dugaan langsung proporsi
    nT.jk <- sum(n[-l])
    p_ib_jk <- n[-l]/nT.jk * p_i_jk
    p_hat_jk <- sum(p_ib_jk) #rataan terboboti
    s_p2_jk <- n[-l]/nT.jk * (p_i_jk-p_hat_jk)^2
    sum_s_p2_jk <- sum(s_p2_jk) #ragam terboboti
    n2nT.jk <- (n[-l]^2) / nT.jk
    sum_n2nT.jk <- sum(n2nT.jk)
    k11.jk <- (nT.jk*sum_s_p2_jk) - p_hat_jk*(1-p_hat_jk) * (m-1)
    k12.jk <- p_hat_jk*(1-p_hat_jk) * (nT.jk-sum_n2nT.jk-(m-1))
    k1.jk <- k11.jk/k12.jk
    alpha <- p_hat_jk*(1-k1.jk)/k1.jk #alpha
    beta <- alpha*(1/p_hat_jk-1) #beta
    k21 <- (y+alpha)*(n-y+beta)
    k22 <- (n+alpha+beta+1)*(n+alpha+beta)^2
    p_i_hat_EB1 <- (y+alpha)/(n+alpha+beta) #penduga bayes
    var_p_i_hat_EB <- k21/k22 #ragam posterior
    #pendugaan empirical bayes bagi p_i
    gamma <- n/(n+alpha+beta) #nilai gamma
    pi_EB = gamma*p_i+(1-gamma)*p_hat #penduga empirical bayes
    result <- list(gamma = gamma, pi_EB = pi_EB, var_p_i_hat_EB = var_p_i_hat_EB)
    return(result)
  }
  jk <- lapply(1:m, function(l) jackknife(y, n, l))
  M1 <- sapply(1:m, function(i){
    m1 <- var_p_i_hat_EB[i]-(m-1)/m*sum(sapply(1:m, function(l){
      return(jk[[l]]$var_p_i_hat_EB[i]-var_p_i_hat_EB[i])
    }))
    return(m1)
  })
  M2 <- sapply(1:m, function(i){
    m2 <- ((m-1)/m)*sum(sapply(1:m, function(l){
      return((jk[[l]]$pi_EB[i]-pi_EB[i])^2)
    }))
    return(m2)
  })
  mse <- M1+M2
  result$MSE.pi_EB$mse <- mse
  return(result)
}

#' @export
saelognom <- function(y, n){
  result <- list(pi_eb = NA, Parameter = list(miu = NA, sigma = NA), MSE.pi_eb
                 = list(method = "Jackknife", mse = NA), direct = list(est = NA, mse = NA))

  if (any(is.na(y)))
    stop("Argument y contains NA values.")
  if (any(is.na(n)))
    stop("Argument n contains NA values.")

  m <- length(y) #jumlah area

  #menghitung dugaan langsung proporsi dan variansnya
  pi = y/n #dugaan langsung proporsi
  mse.pi <- pi * (1-pi) / n #varians dugaan langsung
  result$direct$est <- pi
  result$direct$mse <- mse.pi


  #menduga parameter sebaran logit-normal miu dan sigma
  logit = log(pi/(1-pi))

  #nilai miu dan sigma
  miu=mean(logit)#perhitungan rata-rata logit
  sigma=sd(logit) #perhitungan standard deviasi logit


  #ragam posterior bagi pi
  g = ((y+miu)*(n-y+sigma)) / ((miu+n+sigma+1)*(miu+n+sigma)^2)

  #pendugaan empirical bayes bagi pi
  zi <- (logit-miu)/sigma
  par = miu + (sigma * zi)
  h1 = exp(par) / (1+exp(par))
  h2 = par*y-n*log(1+exp(par))
  nz=(1/(2*3.14159))*exp(-1/2*((zi^2)))
  a = h1*exp(h2)*nz
  b = exp(h2)*nz
  pi_eb = a/b

  result$Parameter$miu <- miu
  result$Parameter$sigma <- sigma
  result$pi_eb <- pi_eb

  #Menghitung MSE penduga empirical bayes
  jackknife <- function(y, n, l){
    #pi_jk <- 0
    #n.jk <- 0
    pi_jk <- y[-l]/n[-l] #dugaan langsung proporsi
    mse.pi_jk <- pi_jk * (1-pi_jk) / n[-l] #varians dugaan langsung
    logit_jk = log(pi_jk/(1-pi_jk))

    #nilai miu dan sigma
    miu_jk=mean(logit_jk)#perhitungan rata-rata logit
    sigma_jk=sd(logit_jk) #perhitungan standard deviasi logit

    #pendugaan empirical bayes bagi pi
    zi_jk <- (logit_jk-miu_jk)/sigma_jk
    par_jk = miu_jk + (sigma_jk * zi_jk)
    h1_jk = exp(par_jk) / (1+exp(par_jk))
    h2_jk = par_jk*y[-l]-n[-l]*log(1+exp(par_jk))
    nz_jk=(1/(2*3.14159))*exp(-1/2*((zi_jk^2)))
    a_jk = h1_jk*exp(h2_jk)*nz_jk
    b_jk = exp(h2_jk)*nz_jk
    pi_eb_jk = a_jk/b_jk #penduga bayes
    g_jk = ((y+miu_jk)*(n-y+sigma_jk)) / ((miu_jk+n+sigma_jk+1)*((miu_jk+n+sigma_jk)^2)) #ragam posterior
    result <- list(pi_eb = pi_eb_jk, g = g_jk)
    return(result)
  }
  jk <- lapply(1:m, function(l) jackknife(y, n, l))
  M1 <- sapply(1:m, function(i){
    m1 <- g[i]-(m-1)/m*sum(sapply(1:m, function(l){
      return(jk[[l]]$g[i]-g[i])
    }))
    return(m1)
  })
  M2 <- sapply(1:m, function(i){
    m2 <- ((m-1)/m)*sum(sapply(1:m, function(l){
      return((pi_eb[i]-pi_eb[i])^2)
    }))
    return(m2)
  })
  mse <- M1+M2
  result$MSE.pi_eb$mse <- mse

  return(result)
}

#' @export
saebincov <- function(y, non_y, n, x, N, area, unit, data){
  result <- list(saebincov = NA, estcoef = NA, randeff = NA, MSE.pi = list(method = "Jackknife", mse = NA))
  xfrom <- x
  x <- as.character(paste(x, collapse = "+"))
  if(!require('lme4')) {
    install.packages('lme4')
    library('lme4')
  }
  m=nrow(data)
  formula <-
    as.formula(paste("cbind(",y,",",non_y,")~",x,"+(1|",area,")"))
  model <- glmer(formula, data = data, family = binomial(link ="logit"))
  #membagi dua data sample dan non sample
  sample <- na.omit(data)
  non_sample <- data[which(is.na(data$n)),]
  randeff <- model@u

  #random effect tiap area
  if(!require('dplyr')) {
    install.packages('dplyr')
    library('dplyr')
  }
  ar <- length(unique(data$area))
  area.reff <- summarise(group_by(sample, area ))
  area.reff <- cbind(area.reff,randeff)
  data <- left_join(data,area.reff,by="area")
  non_sample <- data[which(is.na(data$randeff)),]
  sample <- anti_join(data,non_sample, by="area")
  #covariates sample unit
  covar.s <- as.matrix(cbind(1,sample[xfrom]))
  #covariates non_sample unit
  covar.ns <- as.matrix(cbind(1,non_sample[xfrom]))

  estcoef <- as.matrix(model@beta)
  #nilai XB unit tersample
  xbeta.s <- covar.s%*%estcoef
  #nilai XB unit tidak tersample
  xbeta.ns <- covar.ns%*%estcoef
  #nilai XB tiap unit tersample dan kode unit
  xbeta.unit.s <- summarise(group_by(sample, unit), area,n,y,N,randeff)
  xbeta.unit.s <- cbind(xbeta.unit.s,xbeta.s)
  #nilai XB+reffar tiap area
  xbeta.unit.s$nilai <- xbeta.unit.s$xbeta.s + xbeta.unit.s$randeff
  #nilai XB tiap unit tidak tersample dan kode unit
  xbeta.unit.ns <- summarise(group_by(non_sample, unit), area,N)
  xbeta.unit.ns <- cbind(xbeta.unit.ns,xbeta.ns)
  #nilai pij cap tersampel
  xbeta.unit.s$pij <- exp(xbeta.unit.s$nilai)/(1+exp(xbeta.unit.s$nilai))
  #nilai pij cap tidak tersample
  xbeta.unit.ns$pij <- exp(xbeta.unit.ns$xbeta.ns)/(1+exp(xbeta.unit.ns$xbeta.ns))
  #nilai pi cap tiap area
  yi.sample <- summarise(group_by(xbeta.unit.s, area), n=sum(N), yi = sum(N*pij))
  yi.non_sample <- summarise(group_by(xbeta.unit.ns, area), n=sum(N), yi = sum(N*pij))
  yi.area <- rbind(yi.sample,yi.non_sample)
  pi.area1 <- summarise(group_by(yi.area, area),n=sum(n),yi=sum(yi), pi =sum(yi)/sum(n) )

  pi.area <- summarise(group_by(yi.area, area), pi =sum(yi)/sum(n) )


  g1 <- pi.area1$pi * (1-pi.area1$pi)/(pi.area1$n) #ragam
  EB <- pi.area1$pi
  result$saebincov <- pi.area
  result$estcoef <- estcoef
  result$randeff <- area.reff


  #Menghitung MSE penduga empirical bayes
  jackknife <- function(data, j){
    data_jk <- data[-j,]
    sample_jk <- na.omit(data_jk)
    non_sample_jk <- data_jk[which(is.na(data_jk$n)),]
    randeff_jk <- model@u


    area.reff_jk <- summarise(group_by(sample, area ))
    area.reff_jk <- cbind(area.reff_jk,randeff_jk)
    data_jk <- left_join(data_jk,area.reff_jk,by="area")
    non_sample_jk <- data_jk[which(is.na(data_jk$randeff_jk)),]
    sample_jk <- anti_join(data_jk,non_sample_jk, by="area")

    #covariates sample unit
    covar.s_jk <- as.matrix(as.matrix(cbind(1,sample_jk[xfrom])))
    #covariates non_sample_jk unit
    covar.ns_jk <- as.matrix(cbind(1,non_sample_jk[xfrom]))

    # ambil nilai koefisien beta model
    estcoef_jk <- as.matrix(model@beta)
    #nilai XB unit tersample
    xbeta.s_jk <- covar.s_jk%*%estcoef_jk
    #nilai XB unit tidak tersample
    xbeta.ns_jk <- covar.ns_jk%*%estcoef_jk
    #nilai XB tiap unit tersample dan kode unit
    xbeta.unit.s_jk <- summarise(group_by(sample_jk, unit), area,n,y,N,randeff_jk)
    xbeta.unit.s_jk <- cbind(xbeta.unit.s_jk,xbeta.s_jk)
    #nilai XB+reffar tiap area

    xbeta.unit.s_jk$nilai <- xbeta.unit.s_jk$xbeta.s_jk+xbeta.unit.s_jk$randeff_jk
    #nilai XB tiap unit tidak tersample dan kode unit
    xbeta.unit.ns_jk <- summarise(group_by(non_sample_jk, unit), area,N)
    xbeta.unit.ns_jk <- cbind(xbeta.unit.ns_jk,xbeta.ns_jk)
    #nilai pij cap tersampel
    xbeta.unit.s_jk$pij <- exp(xbeta.unit.s_jk$nilai)/(1+exp(xbeta.unit.s_jk$nilai))
    #nilai pij cap tidak tersample
    xbeta.unit.ns_jk$pij <- exp(xbeta.unit.ns_jk$xbeta.ns_jk)/(1+exp(xbeta.unit.ns_jk$xbeta.ns_jk))
    #nilai pi cap tiap area
    yi.sample_jk <- summarise(group_by(xbeta.unit.s_jk, area), n=sum(N), yi = sum(N*pij))
    yi.non_sample_jk <- summarise(group_by(xbeta.unit.ns_jk, area), n=sum(N), yi = sum(N*pij))
    yi.area_jk <- rbind(yi.sample_jk,yi.non_sample_jk)
    pi.area1_jk <- summarise(group_by(yi.area_jk, area),n=sum(n),yi=sum(yi), pi =sum(yi)/sum(n) )

    pi.area_jk <- summarise(group_by(yi.area_jk, area), pi =sum(yi)/sum(n) )


    g1_jk <- pi.area1_jk$pi * (1-pi.area1_jk$pi)/(pi.area1_jk$n) #ragam

    result <- list(EB_jk = pi.area_jk$pi, g1_jk = g1_jk)

    return(result)
  }
  jk <- lapply(1:m, function(j) jackknife(data,j))
  M1 <- sapply(1:ar, function(i){
    m1 <- g1[i]-(m-1)/m*sum(sapply(1:m, function(j){
      return(jk[[j]]$g1_jk[i]-g1[i])
    }))
    return(m1)
  })
  M2 <- sapply(1:ar, function(i){
    m2 <- ((m-1)/m)*sum(sapply(1:m, function(j){
      return((jk[[j]]$EB_jk[i]-EB[i])^2)
    }))
    return(m2)
  })
  #ganti lapply sapply
  mse <- M1+M2
  result$MSE.pi$mse <- mse

  return(result)
}

#' @export
saebincov.clust <- function(y, non_y, n, x, N, area, unit, data, data_clust, k){
  result <- list(saebincov = NA, estcoef = NA, randeff = NA,cluster= list(cluster = NA, cluster_mean = NA), MSE.pi = list(method = "Jackknife", mse = NA))
  xfrom <- x
  x <- as.character(paste(x, collapse = "+"))
  #clustering
  data.clusterx <- data_clust[xfrom]
  data.clusterx.std <- scale(data.clusterx)
  kmeans_clustering <- kmeans(x = data.clusterx.std, centers = k, nstart = 25)

  hasil.cluster=data.frame(data_clust[,1], kmeans_clustering$cluster)
  colnames(hasil.cluster) <- c("area","cluster")


  if(!require('lme4')) {
    install.packages('lme4')
    library('lme4')
  }
  if(!require('cluster')) {
    install.packages('cluster')
    library('cluster')
  }
  if(!require('dplyr')) {
    install.packages('dplyr')
    library('dplyr')
  }
  m=nrow(data)
  formula <-
    as.formula(paste("cbind(",y,",",non_y,")~",x,"+(1|",area,")"))
  model <- glmer(formula, data = data, family = binomial(link ="logit"))
  #membagi dua data sample dan non sample
  sample <- na.omit(data)
  non_sample <- data[which(is.na(data$n)),]
  randeff <- model@u

  #random effect tiap area

  ar <- length(unique(data$area))
  area.reff <- summarise(group_by(sample, area ))
  area.reff <- cbind(area.reff,randeff)
  data <- left_join(data,area.reff,by="area")
  non_sample <- data[which(is.na(data$randeff)),]
  sample <- anti_join(data,non_sample, by="area")
  area.reff <- left_join(area.reff,hasil.cluster, by="area")
  reff.mean <- summarise(group_by(area.reff,cluster), rata = mean(randeff))
  non_sample <-left_join(non_sample,hasil.cluster, by="area")
  #covariates sample unit
  covar.s <- as.matrix(cbind(1,sample[xfrom]))
  #covariates non_sample unit
  covar.ns <- as.matrix(cbind(1,non_sample[xfrom]))

  estcoef <- as.matrix(model@beta)
  #nilai XB unit tersample
  xbeta.s <- covar.s%*%estcoef
  #nilai XB unit tidak tersample
  xbeta.ns <- covar.ns%*%estcoef
  #nilai XB tiap unit tersample dan kode unit
  xbeta.unit.s <- summarise(group_by(sample, unit), area,n,y,N,randeff)
  xbeta.unit.s <- cbind(xbeta.unit.s,xbeta.s)
  #nilai XB+reffar tiap area
  xbeta.unit.s$nilai <- xbeta.unit.s$xbeta.s + xbeta.unit.s$randeff
  #nilai XB tiap unit tidak tersample dan kode unit
  xbeta.unit.ns <- summarise(group_by(non_sample, unit), area,N, cluster)
  xbeta.unit.ns <- left_join(xbeta.unit.ns,reff.mean, by="cluster")
  xbeta.unit.ns <- cbind(xbeta.unit.ns,xbeta.ns)

  xbeta.unit.ns$nilai <- xbeta.unit.ns$xbeta.ns+xbeta.unit.ns$rata
  #nilai pij cap tersampel
  xbeta.unit.s$pij <- exp(xbeta.unit.s$nilai)/(1+exp(xbeta.unit.s$nilai))
  #nilai pij cap tidak tersample
  xbeta.unit.ns$pij <- exp(xbeta.unit.ns$nilai)/(1+exp(xbeta.unit.ns$nilai))
  #nilai pi cap tiap area
  yi.sample <- summarise(group_by(xbeta.unit.s, area), n=sum(N), yi = sum(N*pij))
  yi.non_sample <- summarise(group_by(xbeta.unit.ns, area), n=sum(N), yi = sum(N*pij))
  yi.area <- rbind(yi.sample,yi.non_sample)
  pi.area1 <- summarise(group_by(yi.area, area),n=sum(n),y=sum(yi), pi =sum(yi)/sum(n) )

  pi.area <- summarise(group_by(yi.area, area), pi =sum(yi)/sum(n) )


  g1 <- pi.area1$pi * (1 - pi.area1$pi)/(pi.area1$n) #ragam
  EB <- pi.area1$pi
  result$saebincov <- pi.area
  result$estcoef <- estcoef

  result$randeff <- area.reff

  result$cluster$cluster <- hasil.cluster
  result$cluster$cluster_mean <- reff.mean

  #Menghitung MSE penduga empirical bayes
  jackknife <- function(data, j){
    #clustering
    data_jk <- data[-j,]
    sample_jk <- na.omit(data_jk)
    non_sample_jk <- data_jk[which(is.na(data_jk$n)),]
    randeff_jk <- model@u


    area.reff_jk <- summarise(group_by(sample, area ))
    area.reff_jk <- cbind(area.reff_jk,randeff_jk)
    data_jk <- left_join(data_jk,area.reff_jk,by="area")
    non_sample_jk <- data_jk[which(is.na(data_jk$randeff_jk)),]
    sample_jk <- anti_join(data_jk,non_sample_jk, by="area")
    area.reff_jk <- left_join(area.reff_jk,hasil.cluster, by="area")
    reff.mean_jk <- summarise(group_by(area.reff_jk,cluster), rata = mean(randeff_jk))
    non_sample_jk <-left_join(non_sample_jk,hasil.cluster, by="area")
    #covariates sample unit
    covar.s_jk <- as.matrix(as.matrix(cbind(1,sample_jk[xfrom])))
    #covariates non_sample_jk unit
    covar.ns_jk <- as.matrix(cbind(1,non_sample_jk[xfrom]))

    # ambil nilai koefisien beta model
    estcoef_jk <- as.matrix(model@beta)
    #nilai XB unit tersample
    xbeta.s_jk <- covar.s_jk%*%estcoef_jk
    #nilai XB unit tidak tersample
    xbeta.ns_jk <- covar.ns_jk%*%estcoef_jk
    #nilai XB tiap unit tersample dan kode unit
    xbeta.unit.s_jk <- summarise(group_by(sample_jk, unit), area,n,y,N,randeff_jk)
    xbeta.unit.s_jk <- cbind(xbeta.unit.s_jk,xbeta.s_jk)
    #nilai XB+reffar tiap area
    xbeta.unit.s_jk$nilai <- xbeta.unit.s_jk$xbeta.s_jk+xbeta.unit.s_jk$randeff_jk
    #nilai XB tiap unit tidak tersample dan kode unit
    xbeta.unit.ns_jk <- summarise(group_by(non_sample_jk, unit), area,N, cluster)
    xbeta.unit.ns_jk <- left_join(xbeta.unit.ns_jk,reff.mean_jk, by="cluster")
    xbeta.unit.ns_jk <- cbind(xbeta.unit.ns_jk,xbeta.ns_jk)

    xbeta.unit.ns_jk$nilai <- xbeta.unit.ns_jk$xbeta.ns_jk+xbeta.unit.ns_jk$rata
    #nilai pij cap tersampel
    xbeta.unit.s_jk$pij <- exp(xbeta.unit.s_jk$nilai)/(1+exp(xbeta.unit.s_jk$nilai))
    #nilai pij cap tidak tersample
    xbeta.unit.ns_jk$pij <- exp(xbeta.unit.ns_jk$nilai)/(1+exp(xbeta.unit.ns_jk$nilai))
    #nilai pi cap tiap area
    yi.sample_jk <- summarise(group_by(xbeta.unit.s_jk, area), n=sum(N), yi = sum(N*pij))
    yi.non_sample_jk <- summarise(group_by(xbeta.unit.ns_jk, area), n=sum(N), yi = sum(N*pij))
    yi.area_jk <- rbind(yi.sample_jk,yi.non_sample_jk)
    pi.area1_jk <- summarise(group_by(yi.area_jk, area),n=sum(n),y=sum(yi), pi =sum(yi)/sum(n) )

    pi.area_jk <- summarise(group_by(yi.area_jk, area), pi =sum(yi)/sum(n) )



    g1_jk <- pi.area1_jk$pi * (1 - pi.area1_jk$pi)/(pi.area1_jk$n) #ragam
    result <- list(EB_jk = pi.area_jk$pi, g1_jk = g1_jk)

    return(result)
  }
  jk <- lapply(1:m, function(j) jackknife(data,j))
  M1 <- sapply(1:ar, function(i){
    m1 <- g1[i]-(m-1)/m*sum(sapply(1:m, function(j){
      return(jk[[j]]$g1_jk[i]-g1[i])
    }))
    return(m1)
  })
  M2 <- sapply(1:ar, function(i){
    m2 <- ((m-1)/m)*sum(sapply(1:m, function(j){
      return((jk[[j]]$EB_jk[i]-EB[i])^2)
    }))
    return(m2)
  })
  #ganti lapply sapply
  mse <- M1+M2
  result$MSE.pi$mse <- mse

  return(result)
}
