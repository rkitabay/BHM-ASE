library(R2jags)
library(mbend)

#---------------Trial Settings---------------#
seed <- 123
Npat <- 20 #target numbers of patients
Nbask <- 4 #numbers of cancer types
pi.null <- 0.05 #null response rate
pi.alt <- 0.20 #alternative response rate
TIE <- 0.05 #cancer-specific type 1 error rate
Nsim <- 2000 #numbers of simulations

true.pi <- c(rep(pi.null, 3), rep(pi.alt, Nbask-3)) #true response rate for each cancer type


#---------------Create the trial data---------------#

set.seed(seed)
x <- rbinom(size = Npat, prob = true.pi, n = Nbask)
data<-data.frame(i = rep(0, Nbask), x = x, n = Npat)

data


#---------------Analysis of BHM-ASE---------------#
#----setting of tuning parameters----
c0 <- 0.99 #large value using the equation (3) in Section 2.1 
c1 <- 0.01 #small value using the equation (3) in Section 2.1
q01.star <- c(0.01,0.25,0.50,0.75,0.99)
#candidate of the paired values of q0, q1 determined in Step 1 of Section 2.1
loss.par <- data.frame(ga = 0.05, gb = 1.00, wa = 1, wb = 1)
#the setting of the loss function in the equation (6)
#ga, gb, wa, and wb are alpha, beta, w_alpha, and w_beta in the equation (6), respectively.

#----determination q0 and q1----
Tq01 <- analytic_Tdist(n = Npat, pi1 =pi.null, pi2 = pi.alt, Aq01.star = q01.star)
#The statistic T corresponding to q0.star and q1.star calculated by the equation (1) in Section 2.1

gs.result <- read.table(file = "grid_search_example.csv", sep = ",", header = TRUE)
#The example of the true positive rate and false positive rate for declaring the drug efficacy
#for each q0.star, q1.star, and m in Step 3 of Section 2.1

wm <- omega_m(data, Nbask, pi.null, pi.alt)
#omega_m calculated by the equation (7) and (8)

q01 <- determine_Tq(gs.result, wm = wm, Nbask, pi.null, pi.alt, loss.par)
#The dicision of q0 and q1 using the average loss in the equation (6) in Step 4

Tq0 <- Tq01$Tq0[q01.star==q01$q0]
Tq1 <- Tq01$Tq1[q01.star==q01$q1]

#----calculation of similarities----
Tkk <- matrix(nrow = Nbask, ncol = Nbask)
#statistic T for each combination of the cancer type k and k'
for(i in 1:Nbask){
  for(j in i:Nbask){
    if(i == j) Tkk[i,j]<-0
    else{
      table <- matrix(c(data$x[i], data$n[i] - data$x[i],
                        data$x[j], data$n[j] - data$x[j]), ncol = 2, byrow = T)
      if(table[1,1] == table[2,1] && table[2,1] == 0) Tkk[i,j] <- 0
      else Tkk[i,j] <- scaled_X2(table)
    }
    Tkk[j,i] <- Tkk[i,j]
  }
}

similarity <- matrix(nrow = Nbask, ncol = Nbask)
#similarity for each combination of the cancer type k and k'
for(i in 1:Nbask){
  for(j in i:Nbask){
    if(i == j) similarity[i,j] <- 1
    else{
      table <- matrix(c(data$x[i], data$n[i] - data$x[i],
                        data$x[j], data$n[j] - data$x[j]), ncol = 2, byrow = T)
      if(table[1,1] == table[2,1] && table[2,1] == 0) similarity[i,j] <- 1
      else similarity[i,j] <- as.numeric(elastic_function(t = Tkk[i,j], Tq0 = Tq0,
                                                          Tq1 = Tq1, c0 = c0, c1 = c1))
    }
    similarity[j,i] <- similarity[i,j]
  }
}

#Using BHM-ASE for analysis, MCMC calculation may not be implemented 
#because the covariance matrix for theta_k, R, become negative definite. 
#In this study, R was adjusted iteratively until all eigenvalues were 
#greater than or equal to the minimum value of 0.05 before MCMC,
#using the procedure proposed by Jorjani et al (2003).  
similarity[similarity > 0.99] <- 0.99
diag(similarity) <- 1
if(sum(eigen(similarity)$values < 0) >= 1){
  bend.matrix <- suppressMessages(bend(similarity, small.positive = 0.05))
  similarity <- bend.matrix$bent
}

#----implementation of MCMC----
logit.pi <- log(((pi.alt + pi.null) / 2) / (1 - ((pi.alt + pi.null) / 2)))
data.mcmc <- list("n" = data$n, "x" = data$x, "S" = similarity, "K" = Nbask, 'mu0' = logit.pi)

fit <- jags(model.file = "BHM-ASE.txt", data = data.mcmc,
            parameters.to.save = c("tau1", "tau2", "pi"),
            n.chains = 3, n.iter=10000, n.thin=3, n.burnin = 2000, quiet=T, DIC = F)

#----summary of results----
sum.fit <- fit[["BUGSoutput"]][["summary"]]

ps <- colMeans(as.data.frame(fit[["BUGSoutput"]][["sims.list"]]$p > pi.null))
mn <- sum.fit[,"mean"]
sd <- sum.fit[,"sd"]
cd <- sum.fit[,"Rhat"]
lc <- sum.fit[,"2.5%"]
uc <- sum.fit[,"97.5%"]
R <- similarity

list(ps = ps, mn = mn, sd = sd, cd = cd, R = R, lc = lc, uc = uc)
