#' MONITORING STATISTIQUE DES ESSAIS : TRAVAUX DE SIMULATION
#'
#' MONITORING STATISTIQUE CENTRALISE (CSM) : DATA GENERATION FOR CONTINUOUS VARIABLE
#'
#' @param n.subject
#' @param n.center default value \code{6}
#' @import data.table
#'
#' @return Table
#
#' @export
#' @examples
#' n.subject = 50
#' n.center = 6
#' sigma2c = 1
#' sigma2r = 4
#'
#' data.generebis( n.subject = n.subject, n.center = n.center, sigma2c = sigma2c, sigma2r = sigma2r,mu = mu,centreAtypic = Taux,delta = DeltaVarie)
#'



## Génération échantillon avec plusieurs sites atypiques

data.generebis <- function(n.subject = 50,
                           n.center = 6,
                           sigma2c = 0.7,
                           sigma2r = 0.5,
                           mu = 100,
                           centreAtypic = 2,
                           delta = 10) {
  data <- data.frame(1:(n.subject * n.center))

  names(data) <-"id"
  d <- NULL
  data$center <- rep(0, n.subject * n.center)
  data$Assignation <- rep(2, n.subject * n.center)
  data$Patient <- rep(1:n.subject, n.center)
  data$gamma.i <- rep(0, n.subject * n.center)

  a <- rnorm(n.center, mean = 0, sd = sqrt(sigma2c)) #On peut mettre data$a pour garder dans la table les n.center                                                       valeurs de la variable a mais repliquer en n.center*n.subject

  for (i in 1:n.center) {
    d <- c(d,rep(i,n.subject))
    data$gamma.i[(n.subject*(i-1)+1):(i*n.subject)] <- rep(a[i], n.subject)
  }
  data$center <- d
  data$epsilon.ij <- rnorm(n.subject * n.center, mean = 0, sd = sqrt(sigma2r))
  data$y.ij_NCR <- mu + data$gamma.i + data$epsilon.ij

  #data[data$center == centreAtypic, "y.ij_NCR"] <- data$y.ij_NCR[data$center == centreAtypic] + delta
  for (i in (n.center-centreAtypic+1):n.center) {
    data[data$center == i, "y.ij_NCR"] <- data$y.ij_NCR[data$center == i] + delta*0.1*mu
    data[data$center == i, "Assignation"] <- 1
  }

  DATA_Bis_Bis <- data.table::data.table(data)
  DATA_Bis_Bis$y.ij <- (DATA_Bis_Bis$y.ij_NCR-mu)/sqrt(sigma2c+sigma2r/n.subject)
  return(DATA_Bis_Bis)
}

