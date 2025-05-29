# necessary packages ----

if(!require(multcomp)) {
  install.packages("multcomp")
  library(multcomp)
}

if(!require(multcompView)) {
  install.packages("multcompView")
  library(multcompView)
}

# significance function
significance.test <- function(model, hypothesis = 1, option = 1) {
  
  # comparison matrix (all possible pairs)
  mat_ex <- 
    cbind(contrMat(rep(1,
                       length(model$ci.lb)),
                   type = "Tukey"))
  
  # when there are more combinations than needed, matrix must be reduced
  if(hypothesis == 1) {
    if(length(model$ci.lb) == 6) {
      if(option == 1) {
        mat_ex <- mat_ex[c(1, 2, 6),]
      } else {
        mat_ex <- mat_ex[c(13:15),]
      }
    }
  } else {
    if(option == 1) {
      mat_ex <- mat_ex[c(1, 2, 12),]
    } 
    else if(option == 2) {
      mat_ex <- mat_ex[c(31:32, 39),]
    }
    else if(option == 3) {
      mat_ex <- mat_ex[c(52:53, 57),]
    }
    else {
      mat_ex <- mat_ex[c(64:66),]
    }
  }
  
  x <- summary(glht(model,
                    linfct = mat_ex),
               test = adjusted("none"))
  
  return(x)
}
