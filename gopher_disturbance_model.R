#-------------------------------------------------------------------------------
# Name: gopher_disturbance_model.R
#
# Notes: parameter names used in the code may differ from ones in the manuscript
#        ones used in the manuscript are showed in the parentheses
#        g: Fraction of area covered by new gopher mounds each year (f)
#        a and b:parameters controls vegetation recovery on gopher mounds (they controls I0-I6)
#        l: soil C input (I)
#        k: soil oil C decomposition fraction in the undisturbed area (k)
#        kd: A constant that controls soil C decomposition fraction in 
#            soil covered by first to sixth-year gopher mounds (d)
#        p: C in the biomass consumed by pocket gophers (s)
#        tmax: maximum time step
#        b_age14 and b_intercept: parameters control gopher mound abundance change with succession
#
# Author:  Yi Yang and Chad Brassil
# Script created on: 03/15/2021
# Last modified on: 08/24/2021
#-------------------------------------------------------------------------------

# gopher model with constant disturbance ------------------------
dfrac_fun_g <- 
  function(g, a, b, l, k, kd, p, tmax) {
  q <- 0.12 # The ratio of buried plant biomass to biomass in an undisturbed area
  
  bmrcfn <- function(x) a * (1 - exp(-b * x)) #vegetation recovery model
  cinputpct <- numeric()
  cinputpct[1] <- q
  for(i in 2:20) cinputpct[i] = bmrcfn(i-1)/a
  
  # l = 97 # C input (g C m-2)
  # Consumption each percent of mound cover
  s <- 0.03
  
  # Model ###
  C <- rep(0, tmax) # mean carbon pool in given area (g C m-2)
  C0 <- rep(0, tmax) # carbon pool in undisturbed/recovered area (g C m-2)
  C1 <- rep(0, tmax) # carbon pool in 1 year old mounds (g C m-2)
  C2 <- rep(0, tmax) # carbon pool in 2 year old mounds (g C m-2)
  C3 <- rep(0, tmax) # carbon pool in 3 year old mounds (g C m-2)
  C4 <- rep(0, tmax) # carbon pool in 4 year old mounds (g C m-2)
  C5 <- rep(0, tmax) # carbon pool in 5 year old mounds (g C m-2)
  C6 <- rep(0, tmax) # carbon pool in 6 year old mounds (g C m-2)
  f0 <- rep(0, tmax) # fraction of undisturbed/recovered area
  f1 <- rep(0, tmax) # fraction of 1 year old mounds
  f2 <- rep(0, tmax) # fraction of 2 year old mounds
  f3 <- rep(0, tmax) # fraction of 3 year old mounds
  f4 <- rep(0, tmax) # fraction of 4 year old mounds
  f5 <- rep(0, tmax) # fraction of 5 year old mounds
  f6 <- rep(0, tmax) # fraction of 6 year old mounds
  #parameter values
  C0[1] <- 1123
  f0[1] <- 1
  C[1] <-C0[1] 
  I0 <- l # carbon input from vegetation in undisturbed area (g C m-2)
  I1 <- cinputpct[1] * l - s * g * p # carbon input to 1 year old mounds (g C m-2)
  I2 <- cinputpct[2]*l # carbon input to 2 year old mounds (g C m-2)
  I3 <- cinputpct[3]*l # carbon input to 3 year old mounds (g C m-2)
  I4 <- cinputpct[4]*l # carbon input to 4 year old mounds (g C m-2)
  I5 <- cinputpct[5]*l # carbon input to 5 year old mounds (g C m-2)
  I6 <- cinputpct[6]*l # carbon input to 6 year old mounds (g C m-2)
  #gopher model
  #when j = 1, time(t) = 0
  for (j in 2:tmax) {
    f1[j] <- g
    f2[j] <- f1[j-1]
    f3[j] <- f2[j-1]
    f4[j] <- f3[j-1]
    f5[j] <- f4[j-1]
    f6[j] <- f5[j-1]
    f0[j] <- f0[j-1] - g + f6[j-1]
    O0 <- k * C0[j-1] # carbon output (decomposition) in undisturbed area (g C m-2)
    O1 <- 0.7 * O0
    O2 <- k * C1[j-1]
    O3 <- k * C2[j-1]
    O4 <- k * C3[j-1]
    O5 <- k * C4[j-1]
    O6 <- k * C5[j-1]
    C0[j] <- ((C0[j-1] + I0 - O0) * (f0[j-1] - f1[j]) + 
                (C6[j-1] + I0 - O0) * f6[j-1])/(f0[j-1] - f1[j] + f6[j-1])
    C1[j] <- C0[j-1] + I1 - O1
    C2[j] <- ifelse(f2[j] == 0, 0, C1[j-1] + I2 - O2)
    C3[j] <- ifelse(f3[j] == 0, 0, C2[j-1] + I3 - O3)
    C4[j] <- ifelse(f4[j] == 0, 0, C3[j-1] + I4 - O4)
    C5[j] <- ifelse(f5[j] == 0, 0, C4[j-1] + I5 - O5)
    C6[j] <- ifelse(f6[j] == 0, 0, C5[j-1] + I6 - O6)
    C[j] <- f0[j] * C0[j] + f1[j] * C1[j] + f2[j] * C2[j] + f3[j] * C3[j] + 
      f4[j] * C4[j] + f5[j] * C5[j] + f6[j] * C6[j] 
  }
  return(C[tmax])
}


# gopher model with succession -----------------------------------
dfrac_fun <-
  function(a, b, l, k, kd, p, b_age14, b_intercept, tmax) {
    bmrcfn <-
      function(x)
        a * (1 - exp(-b * x)) #vegetation recovery model
    cinputpct <- numeric()
    cinputpct[1] <- q
    for (i in 2:20)
      cinputpct[i] = bmrcfn(i - 1) / a
    q <- 0.12
    # l = 97 # C input (g C m-2)
    k = 0.03 # (year-1)
    # Consumption each percent of mound cover
    s <- 0.03
    # gopher consumption
    p <- 2374
    # Model ###
    t <- 1:tmax
    gmmdat.0 <- NULL
    gmmdat.g <- NULL
    gmmdat <- NULL
    disturbance <- c(0.00)
    g <- rep(0, tmax)
    C <- rep(0, tmax) # mean carbon pool in given area (g C m-2)
    C0 <-
      rep(0, tmax) # carbon pool in undisturbed/recovered area (g C m-2)
    C1 <- rep(0, tmax) # carbon pool in 1 year old mounds (g C m-2)
    C2 <- rep(0, tmax) # carbon pool in 2 year old mounds (g C m-2)
    C3 <- rep(0, tmax) # carbon pool in 3 year old mounds (g C m-2)
    C4 <- rep(0, tmax) # carbon pool in 4 year old mounds (g C m-2)
    C5 <- rep(0, tmax) # carbon pool in 5 year old mounds (g C m-2)
    C6 <- rep(0, tmax) # carbon pool in 6 year old mounds (g C m-2)
    f0 <- rep(0, tmax) # fraction of undisturbed/recovered area
    f1 <- rep(0, tmax) # fraction of 1 year old mounds
    f2 <- rep(0, tmax) # fraction of 2 year old mounds
    f3 <- rep(0, tmax) # fraction of 3 year old mounds
    f4 <- rep(0, tmax) # fraction of 4 year old mounds
    f5 <- rep(0, tmax) # fraction of 5 year old mounds
    f6 <- rep(0, tmax) # fraction of 6 year old mounds
    #parameter values
    C0[1] <- 1123
    f0[1] <- 1
    C[1] <- C0[1]
    g[1] <- 0.06
    
    for (j in 2:tmax) {
      g[j] <-
        exp(b_age14 * t[j] + b_intercept) / (1 + exp(b_age14 * t[j] + b_intercept))
      I0 <-
        l # carbon input from vegetation in undisturbed area (g C m-2)
      I1 <-
        cinputpct[1] * l - s * g[j] * p # carbon input to 1 year old mounds (g C m-2)
      I2 <-
        cinputpct[2] * l # carbon input to 2 year old mounds (g C m-2)
      I3 <-
        cinputpct[3] * l # carbon input to 3 year old mounds (g C m-2)
      I4 <-
        cinputpct[4] * l # carbon input to 4 year old mounds (g C m-2)
      I5 <-
        cinputpct[5] * l # carbon input to 5 year old mounds (g C m-2)
      I6 <-
        cinputpct[6] * l # carbon input to 6 year old mounds (g C m-2)
      
      f1[j] <- g[j]
      f2[j] <- f1[j - 1]
      f3[j] <- f2[j - 1]
      f4[j] <- f3[j - 1]
      f5[j] <- f4[j - 1]
      f6[j] <- f5[j - 1]
      f0[j] <- f0[j - 1] - g[j] + f6[j - 1]
      O0 <-
        k * C0[j - 1] # carbon output (decomposition) in undisturbed area (g C m-2)
      O1 <- kd * O0
      O2 <- k * C1[j - 1]
      O3 <- k * C2[j - 1]
      O4 <- k * C3[j - 1]
      O5 <- k * C4[j - 1]
      O6 <- k * C5[j - 1]
      O7 <- k * C6[j - 1]
      C0[j] <- ((C0[j - 1] + I0 - O0) * (f0[j - 1] - f1[j]) +
                  (C6[j - 1] + I0 - O0) * f6[j - 1]) / (f0[j - 1] - f1[j] + f6[j -
                                                                                 1])
      C1[j] <- C0[j - 1] + I1 - O1
      C2[j] <- ifelse(f2[j] == 0, 0, C1[j - 1] + I2 - O2)
      C3[j] <- ifelse(f3[j] == 0, 0, C2[j - 1] + I3 - O3)
      C4[j] <- ifelse(f4[j] == 0, 0, C3[j - 1] + I4 - O4)
      C5[j] <- ifelse(f5[j] == 0, 0, C4[j - 1] + I5 - O5)
      C6[j] <- ifelse(f6[j] == 0, 0, C5[j - 1] + I6 - O6)
      C[j] <-
        f0[j] * C0[j] + f1[j] * C1[j] + f2[j] * C2[j] + f3[j] * C3[j] +
        f4[j] * C4[j] + f5[j] * C5[j] + f6[j] * C6[j]
    }
    return(C[tmax])
  }


# gopher model with succession and fire -------------------------
dfrac_fun_fire <-
  function(a, b, k, kd, p, b_age14, b_intercept, tmax) {
    bmrcfn <-
      function(x)
        a * (1 - exp(-b * x)) #vegetation recovery model
    cinputpct <- numeric()
    cinputpct[1] <- q
    for (i in 2:20)
      cinputpct[i] = bmrcfn(i - 1) / a
    q <- 0.12
    # l = 97 # C input (g C m-2)
    k = 0.03 # (year-1)
    # Consumption each percent of mound cover
    s <- 0.03
    # gopher consumption
    p <- 2374
    
    #vegetation input, fire cause reduction of C input
    l <- rep(c(97, 51), 100)
    
    
    # Model ###
    t <- 1:tmax
    
    
    gmmdat.0 <- NULL
    gmmdat.g <- NULL
    gmmdat <- NULL
    disturbance <- c(0.00)
    g <- rep(0, tmax)
    C <- rep(0, tmax) # mean carbon pool in given area (g C m-2)
    C0 <-
      rep(0, tmax) # carbon pool in undisturbed/recovered area (g C m-2)
    C1 <- rep(0, tmax) # carbon pool in 1 year old mounds (g C m-2)
    C2 <- rep(0, tmax) # carbon pool in 2 year old mounds (g C m-2)
    C3 <- rep(0, tmax) # carbon pool in 3 year old mounds (g C m-2)
    C4 <- rep(0, tmax) # carbon pool in 4 year old mounds (g C m-2)
    C5 <- rep(0, tmax) # carbon pool in 5 year old mounds (g C m-2)
    C6 <- rep(0, tmax) # carbon pool in 6 year old mounds (g C m-2)
    f0 <- rep(0, tmax) # fraction of undisturbed/recovered area
    f1 <- rep(0, tmax) # fraction of 1 year old mounds
    f2 <- rep(0, tmax) # fraction of 2 year old mounds
    f3 <- rep(0, tmax) # fraction of 3 year old mounds
    f4 <- rep(0, tmax) # fraction of 4 year old mounds
    f5 <- rep(0, tmax) # fraction of 5 year old mounds
    f6 <- rep(0, tmax) # fraction of 6 year old mounds
    #parameter values
    C0[1] <- 1123
    f0[1] <- 1
    C[1] <- C0[1]
    g[1] <- 0.06
    
    for (j in 2:tmax) {
      g[j] <-
        exp(b_age14 * t[j] + b_intercept) / (1 + exp(b_age14 * t[j] + b_intercept))
      I0 <-
        l[j] # carbon input from vegetation in undisturbed area (g C m-2)
      I1 <-
        cinputpct[1] * l[j] - s * g[j] * p # carbon input to 1 year old mounds (g C m-2)
      I2 <-
        cinputpct[2] * l[j] # carbon input to 2 year old mounds (g C m-2)
      I3 <-
        cinputpct[3] * l[j] # carbon input to 3 year old mounds (g C m-2)
      I4 <-
        cinputpct[4] * l[j] # carbon input to 4 year old mounds (g C m-2)
      I5 <-
        cinputpct[5] * l[j] # carbon input to 5 year old mounds (g C m-2)
      I6 <-
        cinputpct[6] * l[j] # carbon input to 6 year old mounds (g C m-2)
      
      f1[j] <- g[j]
      f2[j] <- f1[j - 1]
      f3[j] <- f2[j - 1]
      f4[j] <- f3[j - 1]
      f5[j] <- f4[j - 1]
      f6[j] <- f5[j - 1]
      f0[j] <- f0[j - 1] - g[j] + f6[j - 1]
      O0 <-
        k * C0[j - 1] # carbon output (decomposition) in undisturbed area (g C m-2)
      O1 <- kd * O0
      O2 <- k * C1[j - 1]
      O3 <- k * C2[j - 1]
      O4 <- k * C3[j - 1]
      O5 <- k * C4[j - 1]
      O6 <- k * C5[j - 1]
      O7 <- k * C6[j - 1]
      C0[j] <- ((C0[j - 1] + I0 - O0) * (f0[j - 1] - f1[j]) +
                  (C6[j - 1] + I0 - O0) * f6[j - 1]) / (f0[j - 1] - f1[j] + f6[j -
                                                                                 1])
      C1[j] <- C0[j - 1] + I1 - O1
      C2[j] <- ifelse(f2[j] == 0, 0, C1[j - 1] + I2 - O2)
      C3[j] <- ifelse(f3[j] == 0, 0, C2[j - 1] + I3 - O3)
      C4[j] <- ifelse(f4[j] == 0, 0, C3[j - 1] + I4 - O4)
      C5[j] <- ifelse(f5[j] == 0, 0, C4[j - 1] + I5 - O5)
      C6[j] <- ifelse(f6[j] == 0, 0, C5[j - 1] + I6 - O6)
      C[j] <-
        f0[j] * C0[j] + f1[j] * C1[j] + f2[j] * C2[j] + f3[j] * C3[j] +
        f4[j] * C4[j] + f5[j] * C5[j] + f6[j] * C6[j]
    }
    return(C[tmax])
  }