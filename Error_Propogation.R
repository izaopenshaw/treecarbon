# ============== 1: allodb error distribution ================
df$allodb_sigma
# sum of the variances of each
allodb.sd <- sqrt(sum(df$allodb_sigma))

# ============== 2: Analytical model error ==============
# For the model AGB = a*DBH^b + e, e~N(0,e.sd) from Gonzalez-Akre et al (2021)
# the following is an analytical error propogation

error_allodb <- function(a,b,DBH,sd_DBH=0.01,sd_e,sd_a){
  
  z = sqrt((sd_a/a)^2 + (b*sd_DBH/DBH)^2)
  AGB.sd = sqrt((z)^2 + (sd_e)^2)
  return(AGB.sd)
}
df$Analytical.allodb <- error_allodb(a=df$allodb_a, b=df$allodb_a, DBH = df$DBH.2, sd_e = sqrt(df$allodb_sigma),sd_a=sd(df$allodb_a))

# For the model AGB = a*(wd*H*DBH^2)^b
error_allodb <- function(DBH, DBH_sd=0.01, h, sd_h, wd, sd_wd, a, a_sd, b){
  
  # y = wd*h*DBH^2
  sd_y = sqrt((sd_wd/wd)^2 + 2*(sd_DBH/DBH)^2 + (sd_h/h)^2)
  
  # x = y^b
  sd_x = b*(sd_y/y)
  
  # AGB = a*x = a*(wd*H*DBH^2)^b
  sd_AGB = sqrt((sd_x/x)^2 + (sd_a/a)^2)
  
  return(AGB.sd)
}
df$Analytical.allodb <- error_allodb(a=df$allodb_a, b=df$allodb_a, DBH = df$DBH.2, e.sd = sqrt(df$allodb_sigma),a.sd=sd(df$allodb_a))

resample_error <- function(errorfunction){
  
}

#====================== 3: Biomass function ======================
# From the BIOMASS package. For calculating the plot based sandard deviation
AGCmc <- AGBmonteCarlo(D=as.numeric(df$DBH.2)*100,    WD=as.numeric(df$Wood.Density),
                       H=as.numeric(df$Height.1),     errWD = as.numeric(df$Wood.Density.sd), 
                       errH = as.numeric(df$H.relative.error), Dpropag = 'chave2004', Carbon = FALSE)

#====================== 4: Bootstrap progression ======================
# Code adapted from Justin Moat
pro_error_carbon <- function(vol, volsd, den, densd,biom,biomsd,nruns=10000, returnsv=NULL) {
  vol <- rnorm(nruns,mean=vol,sd=volsd)
  den <- rnorm(nruns,mean=den,sd=densd)
  biomass <- rnorm(nruns,mean=biom,sd=biomsd)
  carbt <- vol * den * biomass
  if (!is.null(returnsv)){
    quantile(carbt,probs=returnsv)
  } else {
    c(mean = mean(carbt),sd= sd(carbt))
  }
}

pro_error_vol <- function(dbh,dbhsd,h,hsd) {
  
  vol_expression <- expression(0.0673*(hc * dbhc^2)^0.976)
  hc <- c(h,hsd)
  dbhc <- c(dbh,dbhsd)
  df <- cbind(hc,dbhc)
  results <- propagate(expr = vol_expression, data = df, type = "stat", do.sim = TRUE, verbose = TRUE, nsim = 100000)
  return(results$sim)
}

# ============== 5: Bootstrap ================


# ============= 6: Metropolis Hastings ===============
# Code adapted from Lin et al (2023)
# var : dataframe of explanatory variables (with intercept) (size : n*k) 
# # X : matrix of n explanatory variables (size n*k)
# # theta : vector coefficients (size k)
# sd : standard deviation of errors (additive errors) 
# param : c(theta, sd) (size : k+1)
# n = number of trees, k = number of coefficients in the model

# Biomass model: AGB = a1(DBH^2 * pho * H)^a2   (a = coeffs, pho = wood density)

likelihood <- function(y, X, param){ 
  theta = param[-length(param)]
  sd = param[length(param)]
  pred = X %*% theta 
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T) 
  sumll = sum(singlelikelihoods)
  return(sumll)
}

# PRIORS # non informative priors on theta 
# sd must be positive 
prior <- function(param){ 
  theta = param[-length(param)] 
  sd = param[length(param)] 
  thetaprior = dnorm(theta, mean=0, sd=10, log = T) 
  sdprior = dunif(sd, min=0, max=30, log = T) 
  return(sum(thetaprior)+sdprior)
}

posterior <- function(y, X, param){ 
  return (likelihood(y, X, param) + prior(param))
} 
# proposal = random walk : normal distribution 
# if sdprop is too low : you won't be exploring all your parameters multidimentional space 
# if sdprop is too high : propositions will be rejected too often 

proposalfunction <- function(param, sdprop){ 
  return(rnorm(length(param) , mean = param, sd= sdprop))
}

run_metropolis_MCMC <- function(y, X, startvalue, iterations, sdprop){ 
  # chain : row i = set of parameters at iteration i 
  chain = array(dim = c(iterations+1, length(startvalue))) 
  colnames(chain) <- names(startvalue) lp = array(dim=iterations+1) 
  
  chain[1,] = startvalue 
  for (i in 1:iterations){ 
    proposal = proposalfunction(chain[i,], sdprop) 
    probab = exp(posterior(y, X, proposal) - posterior(y, X, chain[i,]))
    # probab ranges between 0 and 1
    
    # probab > 1 <=> new proposal is more likely considering the data 
    if (runif(1) < probab){ chain[i+1,] = proposal
    } else { chain[i+1,] = chain[i,] 
    } 
    lp[i+1] = likelihood(y, X, chain[i+1,]) 
  } 
  return(cbind(chain, lp))
}

# ============= 7: Maximum Likelihood ===============
DBH = df$DBH.1   ;  a = df$allodb_a
b = df$allodb_b  ;  e = df$allodb_sigma
y = AGB.calc$AGB.allodb

al.likelihood <- function(y, DBH, a, b, sd){
  pred = a*DBH^b
  
  likelihood = dnorm(y, mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)   
}
start = cbind(a, b)

max.likelihood <- function(y, DBH, start, iterations=1000, sd){
  # chain : row i = set of parameters at iteration i
  chain = array(dim = c (iterations+1, length(start)))
  lp = array (dim= iterations+1)
  previous = al.likelihood(y, DBH, a, b, sd)
  
  for (i in 1:iterations){
    
    proposal = al.likelihood(y, DBH, a, b, sd)
    probab = proposal/previous
    # probab ranges between 0 and 1
    
    # probab > 1 <=> new proposal is more likely considering the data
    if ( runif (1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
    lp[i +1 ] = likelihood (y, X, chain[i+1,])
  }
  return(cbind (chain, lp))
}



################# Metropolis algorithms ########################
# # # # # # # # # # # # # # # MODEL # # # # # # # # # # # # # # # # # # # # # y : vector of predicted values (size n)
# n = number of trees, k = number of coefficients in the model
# var : dataframe of explanatory variables (with intercept) (size : n*k) 
# # X : matrix of n explanatory variables (size n*k)
# # theta : vector coefficients (size k)
# Biomass model: AGB = a1(DBH^2 * pho * H)^a2   (a = coeffs, pho = wood density)
# sd : standard deviation of errors (additive errors) 
# param : c(theta, sd) (size : k+1)
# # #
# # # # # # # # total model : y ~ normal ( X %*% theta , sd^2 ) # # # # #


#=============== Re-sample trees across plots ==========================
AGBmonteCarlo()

# ==== Inputs:
# df: data frame containing columns; Genus, Species, DBH (in cm), Height (even if NA for all entries), site
# coords: either the coordinates of the site, a vector of longitude and latitude
resample.site <- function(df){
  data <- 
}

# Propogate from Individual tree to sum of forest
methods <- data.frame(AGB = c(sum(AGB.calc$TotalVolume),sum(AGB.calc$TotalVolume_opt),sum(AGB.calc$biomass..Kg.),sum(AGB.calc$Volume_opt1000),sum(AGB.calc$AGB.Biomass.kg),sum(AGB.calc$AGB.allodb)),
                      Method = c("Gui.TLS","TLS.opt","Biomass.TLS","TLS.opt.1000","Biomass","allodb"))
methods$Error <- NA

methods$Error[6] <- sum(sqrt(AGB.calc$allodb_sigma))

methods$ymin <- methods$AGB - methods$Error
methods$ymax <- methods$AGB + methods$Error

ggplot(methods, aes(x=Method, y=AGB))+ 
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = "dodge")

# Check distribution
library(fitdistrplus)
descdist(Y, boot=100)
qqnorm(Y, main='Carbon') ; qqline(Y)
shapiro.test(Y) # if Y~N pval > 0.05

sd0 = sd(df$DBH, na.rm=TRUE)

propogate.all <- function(n_runs, mean, sd){
  p_prop <- rnorm(n_runs, mean = mean, sd = sd)
  density <- 
}

pro_error_vol <- function(dbh,dbhsd,h,hsd) {
  
  vol_expression <- expression(0.0673*(hc * dbhc^2)^0.976)
  hc <- c(h,hsd)
  dbhc <- c(dbh,dbhsd)
  df <- cbind(hc,dbhc)
  results <- propagate(expr = vol_expression, data = df, type = "stat", do.sim = TRUE, verbose = TRUE, nsim = 100000)
  return(results$sim)
}

hist(df$DBH)

pro_error_carbon <- function(vol,volsd,den,densd,biom,biomsd,nruns=10000, returnsv=NULL) {
  vol <- rnorm(nruns,mean=vol,sd=volsd)
  den <- rnorm(nruns,mean=den,sd=densd) # middle of the road
  biomass <- rnorm(nruns,mean=biom,sd=biomsd) # conifer or angiosperm
  carbt <- vol * den * biomass
  if (!is.null(returnsv)){
    quantile(carbt,probs=returnsv)
  } else {
    c(mean = mean(carbt),sd= sd(carbt))
  }
}

prop.table()

# Bootstrapping
# https://www.statology.org/bootstrapping-in-r/
set.seed(0)
library(boot)

# define function to calculate R-squared
rsq_function <- function(formula, data, indices) {
  d <- data[indices,] #allows boot to select sample
  fit <- lm(formula, data=d) #fit regression model
  return(summary(fit)$r.square) #return R-squared of model
}
#perform bootstrapping with 2000 replications
reps <- boot(data=mtcars, statistic=rsq_function, R=2000, formula=mpg~disp)

#view results of boostrapping
reps

#calculate adjusted bootstrap percentile (BCa) interval
boot.ci(reps, type="bca")


# Lin code adapted

Uncertain_Boot <- function(Mex_mod,  X_Var_caseStudy, n_iter){
  
  n_casestudy <- length(X_Var_caseStudy)
  
  Y_X <- data.frame(Y = log(Mex_mod[, 2]), X = log(Mex_mod[, 1]))
  
  n_row <- length(Mex_mod[, 2])
  row_index <- replicate(n_iter, sample(n_row, n_row, replace = T)) 
  # replicate(n, expr (expression to evaluate repeatedly), simplify = "array")
  # replicate is a wrapper for the common use of sapply for repeated evaluation of an
  # expression (which will usually involve random number generation).
  
  # sample(x, size, replace = FALSE, prob = NULL)
  # sample takes a sample of the specified size from the elements of x using either with or without replacement.
  
  Y_X_Out <- vector("list", length = n_iter)
  for (i in seq(1:n_iter)){
    Y_X_Out[[i]] <- Y_X[row_index[, i], ]
  }
  
  
  ## Model fitting
  AGB_model <- function(df) {
    lm(Y ~ X, data = df)
  }
  models <- map(Y_X_Out, AGB_model)
  # map = Apply a function to each element of a vector. returns a list. 
  
  
  ## Extract intercept, slope, and sigma from models
  intercept <- map(models, coef) %>%
    map_dbl(1)
  
  slope <- map(models, coef) %>%
    map_dbl(2)
  
  sigma <- map_dbl(models, sigma) 
  
  data4 <- data.frame(intercept, slope, sigma)
  
  
  
  ## Calculate AGB 
  for (i in seq(1:n_iter)){
    #data4$out[[i]] <- exp(data4$intercept[i] + data4$slope[i]*X_Var_caseStudy + (data4$sigma[i])^2/2)
    data4$out1[[i]] <- exp(data4$intercept[i] + data4$slope[i]*X_Var_caseStudy)
  }
  
  
  ## Pre-calculated plot density (cm/ha) based on tree size and nested plot design
  if (n_row ==48){
    PD =6557.6953
  } else if (n_row ==245){
    PD = 322.1505
  } else if (n_row ==93){
    PD = 1527.778
  }
  
  #AGB_Tot_boot <- map_dbl(data4$out, sum)/n_casestudy
  AGB_Tot_boot1 <- (map_dbl(data4$out1, sum)*PD)/n_casestudy
  
  
  return(list(AGB_Tot_boot1))
}

