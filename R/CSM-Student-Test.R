#' CENTRALIZED STATISTICAL MONITORING : STUDENT TEST METHOD
#'
#' @param data a multicenter clinical trial data containing one or several continuous variable
#'
#' @return table
#'
#' @import data.table
#' @import dplyr
#'
#'
#' @export
#'
#' @examples
#'
## Student t-test CSM method for mean comparisons

CONSTRUCT_TestMOY <-function(data=data){
  tab <- data.frame(t(c(1:4)))
  colnames(tab) <- c("Site","MoySite","MoyAutrSite","TestMoyenne")
  n.center <- length(unique(data$center))

  for (i in 1:n.center){
    tab.data <- data
    centre.uniq <- subset(tab.data, center==i)
    centre.uniq.Diff <- subset(tab.data, !(center==i))

    tab.data$SiteRef <-ifelse(tab.data$center == i, 1, 0)
    tempTestVar <- var.test(y.ij ~ SiteRef, data = tab.data)$p.value

    if(tempTestVar >= 0.05){
      tempTestMoyenne <- t.test(y.ij ~ SiteRef, data = tab.data, var.equal=T)$p.value
    }
    if(tempTestVar < 0.05){
      tempTestMoyenne <- t.test(y.ij ~ SiteRef, data = tab.data, var.equal=F)$p.value
    }

    temptab <- data.frame(t(c(1:4)))
    colnames(temptab) <-c("Site","MoySite","MoyAutrSite","TestMoyenne")

    temptab$Site <- i
    temptab$MoySite <- mean(centre.uniq$y.ij)
    temptab$MoyAutrSite <- mean(centre.uniq.Diff$y.ij)
    temptab$TestMoyenne <- tempTestMoyenne

    tab.data <- tab.data[,c(1:6)]

    tab <- rbind(
      tab, temptab)
  }

  SortieResMoyTestSite <- tab[-1,]
  return(SortieResMoyTestSite)

}

