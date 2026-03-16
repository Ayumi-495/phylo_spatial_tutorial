calculate_smd <-
  function(df,
           max_prop = 0.95)
  {
    
    df_copy <- df %>% 
      mutate(smd = NA,
             var = NA)
    
    for(i in seq_len(nrow(df_copy))) {
      
      n1 <- df_copy$n1[i]
      n2 <- df_copy$n2[i]
      sd1 <- df_copy$sd1[i]
      sd2 <- df_copy$sd2[i]
      type <- df_copy$equation_type[i]
      dependency <- df_copy$dependency[i]
      
      # adjusting porportions if too extreme
      # max_prop is an argument in the function
      
      if(type == 'binary') {
        
        min_prop <- 1 - max_prop
        
        if (df_copy$x1[i] > max_prop) { #upper limit
          df_copy$x1[i] <- max_prop
        }
        if (df_copy$x1[i] < min_prop) { #lower limit
          df_copy$x1[i] <- min_prop
        }
        if (!is.na(df_copy$x2[i]) & df_copy$x2[i] > max_prop) { #upper limit
          df_copy$x2[i] <- max_prop
        }
        if (!is.na(df_copy$x2[i]) & df_copy$x2[i] < min_prop) { #lower limit
          df_copy$x2[i] <- min_prop
        }
      }
      
      # point estimates (SMD) ----
      
      x1 <- df_copy$x1[i]
      x2 <- df_copy$x2[i]
      r12 <- 0.5
      
      if (type == 'investment') {
        
        #eq. 2
        s_pooled <- 
          sqrt((((n1 - 1) * (sd1 ^ 2)) + ((n2 - 1) * (sd2 ^ 2))) / (n1 + n2 - 2))
        
        #eq. 1
        df_copy$smd[i] <- (x1 - x2) / s_pooled 
      }
      
      else if(type == 'binary') {
        
        if (dependency == 'no') {
          
          #eq. 4
          df_copy$smd[i] <- (qlogis(x1) - qlogis(x2)) * sqrt(3) / pi
        } else {
          
          #eq. 5
          df_copy$smd[i] <- (qlogis(x1) - qlogis(1 - x1)) * sqrt(3) / pi
        }
      }
      
      else if (type == 'difference') {
        
        #eq. 2
        s_pooled <- 
          sqrt((((n1 - 1) * (sd1 ^ 2)) + ((n2 - 1) * (sd2 ^ 2))) / (n1 + n2 - 2))
        
        #eq. 3
        df_copy$smd[i] <- 
          x1 / s_pooled 
      }
      
      
      # variances ----
      
      smd <- df_copy$smd[i]
      
      if (df_copy$dependency[i] == 'no') { 
        #eq. 6
        df_copy$var[i] <- 
          ((n1 + n2) / (n1 * n2)) + ((smd ^ 2) / (2 * (n1 + n2 - 2)))
      } else { 
        #eq. 7
        df_copy$var[i] <- 
          ((2 * (1 - r12)) / n1) + ((smd ^ 2) / (2 * (n1 - 1)))
      }
      
      
    }
    
    #inverts smd when latency was measured
    df_output <-
      df_copy %>% 
      mutate(smd = replace(smd,
                           measure == "latency",
                           .[.$measure == "latency", ]$smd * (- 1)))
    
    return(df_output)
  }