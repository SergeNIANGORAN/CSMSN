#' CENTRALIZED STATISTICAL MONITORING : DESMET CATEGORICAL CSM METHOD USING BETA-BINOMIAL MODELING
#' Function to implemente the Desmet Categorical CSM method which uses the Betabinomial modelling
#' in order to detect an atypical center in multicenter trial for a categorical binary variable.
#' This method like the others of this package implementes the CSM for the
#'  categorical variable data
#'
#' @param dataHatayama
#'
#' @return Table
#'
#' @importFrom data.table data.table
#' @import jagsUI
#' @import dplyr
#' @import VGAM
#'
#' @export

CSM_BETA_BINOMIAL_SIMS <- function(bdata=bdata){

  Prob <- sum(bdata$y/bdata$N)/length(bdata$N)
  Probw <- sum(bdata$y)/sum(bdata$N)
  mu_m = Prob

  Num_rho_m = sum((bdata$y/bdata$N - Prob)^2) - Prob*(1-Prob)*(1-1/length(bdata$Site))*sum(1/bdata$N)

  Dem_rho_m = Prob*(1-Prob)*((length(bdata$Site))*(1-1/length(bdata$Site)) - (1-1/length(bdata$Site))*sum(1/bdata$N))

  rho_m = Num_rho_m / Dem_rho_m
  if(rho_m <0){ rho_m <- 0}

  S = sum((bdata$y - (bdata$N)*Probw)^2)/(Probw*(1-Probw))

  Z = (S - sum(bdata$N))/sqrt(2*sum((bdata$N)*((bdata$N)-1)))

  # Etape 2.C - Fontion associ?e
  # Calcul de P_value pour une distribution binomiale
  Stage_2C_Function <- function(){
    P_Value_Res <- data.frame(c(1:length(bdata$Site)))
    for(t in 1:length(bdata$Site)){
      if(bdata$y[t] > bdata$N[t]*Prob){
        P_Value_Res[t,] <- min(2*pbinom(bdata$y[t], bdata$N[t], Probw, lower.tail = F),1)
      }
      if(bdata$y[t] <= bdata$N[t]*Prob){
        P_Value_Res[t,] <- min(2*(1-pbinom(bdata$y[t], bdata$N[t], Probw, lower.tail = F)),1)
      }
    }
    colnames(P_Value_Res) <- c("P_Value")
    P_Value_Res$Modele <-"Binomial"
    P_Value_Res$Proba_BetaBino <- Probw
    P_Value_Res$Alpha <- NA
    P_Value_Res$Beta <- NA
    P_Value_Res$NB_Events <- bdata$y
    P_Value_Res$NB_Subjects <- bdata$N
    P_Value_Res$NB_Events_Theoric <- bdata$N*Probw
    return(P_Value_Res)
  }

  # Etape 1.C
  if(Z < 1.96 & rho_m < 10^-3 ){
    # GoTo Etape 2.C
    P_Value_Result <- Stage_2C_Function()
    #   break
  }

  #Stage 3 Function

  Stage_3_Function <- function(mu_L, rho_L){

    gamma_petit = 0.1
    alpha_target = 1
    beta_target = 1

    alpha_hat = mu_L*(1/rho_L - 1)
    beta_hat = (1 - mu_L)*(1/rho_L - 1)

    lambda_de_x <- function(x){
      if(x < gamma_petit | x > (1 - gamma_petit) ){
        Z <- (cos(pi*x/gamma_petit) + 1)/2
      }
      else {
        Z <- 0
      }
      return(Z)
    }
    lambda_de_p <- lambda_de_x(Prob)
    if(Prob < 0.5 & beta_hat < 1){
      alpha_tilde = alpha_hat
      beta_tilde = (1-lambda_de_p)*beta_hat + lambda_de_p*beta_target
    }
    if(Prob < 0.5 & beta_hat >= 1){
      alpha_tilde = alpha_hat
      beta_tilde = beta_hat
    }
    if(Prob >= 0.5 & alpha_hat < 1){
      alpha_tilde = (1-lambda_de_p)*alpha_hat + lambda_de_p*alpha_target
      beta_tilde = beta_hat
    }
    if(Prob >= 0.5 & alpha_hat >= 1){
      alpha_tilde = alpha_hat
      beta_tilde = beta_hat
    }

    alpha_sortie <- data.frame(c(1:length(bdata$Site)))
    beta_sortie <- data.frame(c(1:length(bdata$Site)))
    for (CPTSite in 1:length(bdata$Site)) {

      if((Prob < 0.5 & (2*bdata$y[CPTSite] > bdata$N[CPTSite]))|
         (Prob >= 0.5 & (2*bdata$y[CPTSite] < bdata$N[CPTSite]))){
        alpha_sortie[CPTSite] <- alpha_tilde
        beta_sortie[CPTSite] <- beta_tilde
      }
      else
      {
        alpha_sortie[CPTSite] <- alpha_hat
        beta_sortie[CPTSite] <- beta_hat
      }
      return(list(alpha_sortie, beta_sortie))
    }
  }

  # Stage 2 : Estimation et s?lection de mod?les

  Petit_Delta=10^-6
  ### 2.A
  # Ecriture de la vraisemblance d'une distribution betabinomiale
  LL_Funct <- function(mu_T, sigma_T){
    R = dbetabinom(bdata$y, size = max(bdata$N), prob=mu_T, rho=sigma_T)
    -sum(log(R))
  }

  rm(MLE_Res)
  rm(mu_L)
  rm(rho_L)
  # MLE_Res <- summary(mle(minuslogl = LL_Funct, start = list(mu_T=mu_m, sigma_T=rho_m)))
  fit <- try(summary(mle(minuslogl = LL_Funct, start = list(mu_T=mu_m, sigma_T=rho_m))), silent = TRUE)
  OK <- !inherits(fit, "try-error")

  if (OK==TRUE){
    MLE_Res <- summary(mle(minuslogl = LL_Funct, start = list(mu_T=mu_m, sigma_T=rho_m)))
    mu_L <- MLE_Res@coef[1,1]
    rho_L <- MLE_Res@coef[2,1]
  }

  # Fonction des EQUATIONS 8 ET 9

  Equation_8_9_Function <- function(rho_L=rho_L){
    # Fonctions w(i)
    DoubleV_Vector_function <- function(rho_L=rho_L){
      Poids_i <- data.frame(c(1:length(bdata$Site)))
      for (i in 1:length(bdata$Site)){
        Poids_i[i,] <- bdata$N[i]/(1 + rho_L*(bdata$N[i] - 1))
      }
      return(Poids_i)
    }
    mu_hat_Vector_function <- function(){
      mu_hat_Vector <- data.frame(c(1:length(bdata$Site)))
      for (i in 1:length(bdata$Site)){
        mu_hat_Vector[i,] <- (Poids_i[i,])*(bdata$y[i])/(bdata$N[i])
      }
      mu_hat_Def <- sum(mu_hat_Vector)/sum(Poids_i)
      return(mu_hat_Def)
    }

    # EQUATION 9
    Rho_hat_Vector_function <- function(){
      Rho_hat_Vector_Num <- data.frame(c(1:length(bdata$Site)))
      Rho_hat_Vector_Den <- data.frame(c(1:length(bdata$Site)))
      for (i in 1:length(bdata$Site)){
        Rho_hat_Vector_Num[i,] <- Poids_i[i,]*((bdata$y[i]/bdata$N[i]) - mu_hat_Def)^2 - mu_hat_Def*(1-mu_hat_Def)*((Poids_i[i,]/(bdata$N[i]))*(1- Poids_i[i,]/sum(Poids_i)))
        Rho_hat_Vector_Den[i,] <- mu_hat_Def*(1-mu_hat_Def)*(Poids_i[i,]*(1- Poids_i[i,]/sum(Poids_i)) - ((Poids_i[i,]/(bdata$N[i])))*(1- Poids_i[i,]/sum(Poids_i)))
      }
      Rho_hat_Def <- sum(Rho_hat_Vector_Num)/sum(Rho_hat_Vector_Den)
      return(Rho_hat_Def)
    }

    if (OK==TRUE){
      rho_It <- rho_L
    }
    else{
      rho_It <- rho_m
    }

    mu_hat_Vector_Def <- data.frame(c(1:1000))
    Rho_hat_Vector_Def <- data.frame(c(1:1000))
    for(m in 1:1000){
      Poids_i <- DoubleV_Vector_function(rho_L=rho_It)
      mu_hat_Vector_Def[m,] <- mu_hat_Vector_function()
      mu_hat_Def <- mu_hat_Vector_function()
      Rho_hat_Vector_Def[m,] <- Rho_hat_Vector_function()
      rho_It <- Rho_hat_Vector_Def[m,]

      if(Rho_hat_Vector_Def[m,]< 0){
        Rho_hat_Vector_Def[m,] <- 0
      }
      if(m>1){
        if((abs(Rho_hat_Vector_Def[m,]-Rho_hat_Vector_Def[(m-1),])<0.000001)
           & (abs(mu_hat_Vector_Def[m,]-mu_hat_Vector_Def[(m-1),])<0.000001)){
          Cycl_Point <- m
          mu_final <- mu_hat_Vector_Def[m,]
          rho_final <- Rho_hat_Vector_Def[m,]

          alpha_hat_m <- (mu_hat_Vector_Def[m,])*(1/Rho_hat_Vector_Def[m,]-1)
          beta_hat_m <- (1-mu_hat_Vector_Def[m,])*(1/Rho_hat_Vector_Def[m,]-1)

          break
        }

      }

    }
    return(list(mu_final, rho_final))
  }

  #Stage 4 Function
  # Calcul de P_value pour une distribution beta-binomiale
  Stage_4_Function <- function(AlphaBeta_Vector){

    # defining row names, column names and their lengths
    rown <- c(1:length(bdata$Site))
    coln <- c("P_Value", "Alpha", "Beta", "Modele", "NB_Events","NB_Subjects","NB_Events_Theoric","Proba_BetaBino")
    # creating matrix
    Res_Pvalue_mat <- matrix(nrow = length(rown), ncol=length(coln))
    colnames(Res_Pvalue_mat) <- coln

    #P_Value <- data.frame(c(1:length(bdata$Site)))
    for(l in 1:length(bdata$Site)){
      if(bdata$y[l] > bdata$N[l]*Prob){
        Res_Pvalue_mat[l,1] <- min(2*pbbinom(bdata$y[l], bdata$N[l], alpha =AlphaBeta_Vector[[1]][l,],
                                             beta =AlphaBeta_Vector[[2]][l,], lower.tail = F),1)
      }
      if(bdata$y[l] <= bdata$N[l]*Prob){
        Res_Pvalue_mat[l,1] <- min(2*(1-pbbinom(bdata$y[l], bdata$N[l], alpha =AlphaBeta_Vector[[1]][l,],
                                                beta =AlphaBeta_Vector[[2]][l,], lower.tail = F)),1)
      }
      Res_Pvalue_mat[l,2] <- round(AlphaBeta_Vector[[1]][l,],4)
      Res_Pvalue_mat[l,3] <- round(AlphaBeta_Vector[[2]][l,],4)
      Res_Pvalue_mat[l,4] <- "BetaBinomial"
      Res_Pvalue_mat[l,5] <- bdata$y[l]
      Res_Pvalue_mat[l,6] <- bdata$N[l]
      Res_Pvalue_mat[l,7] <- round(bdata$N[l]*Prob,4)
      Res_Pvalue_mat[l,8] <- round(pbbinom(bdata$y[l], bdata$N[l], alpha =AlphaBeta_Vector[[1]][l,],
                                           beta =AlphaBeta_Vector[[2]][l,], lower.tail = F),6)
    }
    return(data.frame(Res_Pvalue_mat))
  }

  # defining row names, column names and their lengths
  rown <- c(1:length(bdata$Site))
  coln <- c("P_Value", "Alpha", "Beta", "Modele", "NB_Events","NB_Subjects","NB_Events_Theoric")
  # creating matrix
  mat <- matrix(nrow = length(rown), ncol=length(coln))
  colnames(mat) <- coln

  #Stage 2 : Case processing

  if (OK==TRUE){

    if(rho_L >= Petit_Delta & rho_L <= 1-Petit_Delta){
      STAGE3_RESULT <- Stage_3_Function(mu_L=mu_L, rho_L=rho_L)

      P_Value_Result <- Stage_4_Function(STAGE3_RESULT)
    }

    if(rho_L < Petit_Delta){
      if(Z < 1.96 | rho_m < 10^-3){
        # GoTo Etape 2.C
        P_Value_Result <- Stage_2C_Function()
        # break
      }
      else {
        # GoTo Etape 2.B
        mu_I = Equation_8_9_Function(rho_L=rho_L)[[1]]
        rho_I = Equation_8_9_Function(rho_L=rho_L)[[2]]

        alpha_I <- mu_I*(1/rho_I -1)
        alpha_I
        beta_I <- (1 - mu_I)*(1/rho_I -1)
        beta_I
        AlphaBeta_Vector <- Stage_3_Function(mu_L=mu_I, rho_L=rho_I)

        P_Value_Result <- Stage_4_Function(AlphaBeta_Vector)
      }
    }

    if(rho_L > 1 - Petit_Delta){
      rho_L <- 1 - Petit_Delta
      #GoTo Step 3 with (mu_L, rho_L)
      STAGE3_RESULT <- Stage_3_Function(mu_L=mu_L, rho_L=rho_L)

      P_Value_Result <- Stage_4_Function(STAGE3_RESULT)
    }

  }

  if (OK==FALSE){
    mu_I = Equation_8_9_Function(rho_L=rho_m)[[1]]
    rho_I = Equation_8_9_Function(rho_L=rho_m)[[2]]

    if(Z < 1.96 | rho_I < 10^-3){
      # GoTo Etape 2.C
      P_Value_Result <- Stage_2C_Function()
      # break
    }
    else {
      # GoTo Etape 2.B

      alpha_I <- mu_I*(1/rho_I -1)
      alpha_I
      beta_I <- (1 - mu_I)*(1/rho_I -1)
      beta_I
      AlphaBeta_Vector <- Stage_3_Function(mu_L=mu_I, rho_L=rho_I)

      P_Value_Result <- Stage_4_Function(AlphaBeta_Vector)
    }
  }

  P_Value_Result <- cbind(bdata[,c(1,6)], P_Value_Result)
  return(P_Value_Result)
}

##END OF DESMET METHOD DEVELOPPMENT


