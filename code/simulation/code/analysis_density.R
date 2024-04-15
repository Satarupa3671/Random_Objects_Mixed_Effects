rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(frechet)
library(MASS)
library(truncnorm)
##global parameters
mu0 = 0
sigma0 = .3
beta = .3
gamma = .5
rho1 = .25
rho2 = 1
## Geodesic in the space of univariate distributions
q_geo = function(t,A,B){
  (1-t)*A + t*B
}
## perturbation map to generate noisy observations around a true geodesic trajectory; alpha is the bound on the variance of the perturbed observation from the truth as per assumption (A3)
q_pert = function(qtrue,alpha){
  eps = sample(c(alpha, -alpha),size = 1, prob = c(.5,.5))
  return(qtrue + eps*alpha*qtrue*(1-qtrue))
}
##transformation to generate non-normal densities with the same frechet mean (used in Setting IV)
transport = function(k,x) x - (sin(k * x)/abs(k))

## subject level geodesic : generate observed data in the Wasserstein space of univariate distributions under 4 settings and recover the end points of the true geodesics by global Fr regression
subject_level_GFR = function(ind, Z, alpha, ni, setting){
  if(setting == 1){
    muY = function(u,x) rnorm(1, mu0 + beta*x+ .2*u, rho1)
    sigmaY = function(u,x) 1
  } else if (setting == 2){
    muY = function(u,x) rnorm(1, mu0, rho1)
    sigmaY = function(u,x) pmax(.5, rgamma(1, (sigma0 + gamma*x + .2*u)^2/rho2, rho2/(sigma0 + gamma*x + .2*u)))
  } else if (setting == 3){
    muY = function(u,x) rnorm(1, mu0 + beta*x+ .2*u, rho1)
    sigmaY = function(u,x) pmax(.5, rgamma(1, (sigma0 + gamma*x + .2*u)^2/rho2, rho2/(sigma0 + gamma*x + .2*u)))
  } else if (setting == 4){
    muY = function(u,x) rnorm(1, mu0 + beta*x+ .2*u, rho1)
    sigmaY = function(u,x) pmax(.5,rgamma(1, (sigma0 + gamma*x + .2*u)^2/rho2, rho2/(sigma0 + gamma*x + .2*u)))
  }
  qSup = seq(1e-10, 1-1e-10, length.out = 100)
  Ti = c(0, runif(ni[ind]-2, 0.3, .7),1)
  Z_index = Z[ind]
  A = qtruncnorm(qSup, a=0, b=1,  mean =  muY(Ti[1], Z_index), sd = sigmaY(Ti[1], Z_index))
  B = qtruncnorm(qSup, a=0, b=1,  mean =  muY(Ti[ni[ind]], Z_index), sd = sigmaY(Ti[ni[ind]], Z_index))
  #A = qnorm(qSup, mean =  muY(Ti[1], Z_index), sd = sigmaY(Ti[1], Z_index))
  #B = qnorm(qSup, mean =  muY(Ti[ni[ind]], Z_index), sd = sigmaY(Ti[ni[ind]], Z_index))
  if(setting == 4){
    A = qnorm(qSup, mean =  muY(Ti[1], Z_index), sd = sigmaY(Ti[1], Z_index))
    B = qnorm(qSup, mean =  muY(Ti[ni[ind]], Z_index), sd = sigmaY(Ti[ni[ind]], Z_index))
  }
  inputs = lapply(1:(ni[ind]), function(j){
    qtrue = q_geo(Ti[j],A,B) ##true underlying geodesic for each subject
    if(setting == 4){
      k = sample(c(-2,1,1,2),1)
      qtrue = sort(transport(k, qtrue)) 
    }
    qin = q_pert(qtrue,alpha) ## observations around the true geodesic for each subject
    return(list(q_on_geod = qtrue, qin = qin))
  })
  qin = t(sapply(1:ni[ind], function(j) sort(inputs[[j]]$qin)))
  q_on_geod = t(sapply(1:ni[ind], function(j) inputs[[j]]$q_on_geod))
  qSup[1] = 0
  qSup[length(qSup)] = 1
  obj_fit = GloDenReg(xin = Ti, qin = qin, xout =  c(0,1),
                      optns = list(qSup = qSup, ndSup = 100))
  return(list(Ti = Ti, q_on_geod = q_on_geod,
              qin = qin, obj_fit = obj_fit$qout))
}

### Do not run : compute MISE over 200 replications through parallel computation
### outputs are saved in the data folder : used for boxplot
### Go to line 138 
library(parallel)
ncore = 15
B = 150
for(n in c(50,400,1000)){
  res = mclapply(1:B, function(rep){
    Z = sapply(1:n, function(ind)  runif(1, 0, 1))
    ni = sapply(1:n, function(ind) 10+2)
    alpha = .1
    second_step_resp = lapply(1:n, function(ind){
      subject_level_GFR = subject_level_GFR(ind, Z, alpha, ni, setting = setting)
      din = subject_level_GFR$din
      qin = subject_level_GFR$qin
      q_on_geod = subject_level_GFR$q_on_geod
      obj_fit_ind = subject_level_GFR$obj_fit
      mi_hat0 = obj_fit_ind[1,]
      mi_hat1 = obj_fit_ind[2,]
      Ti = subject_level_GFR$Ti
      return(list(q_on_geod = q_on_geod, Ti =Ti,
                  qin = qin, din = din, mi_hat0 = mi_hat0, mi_hat1 = mi_hat1))
    })
    
    Ti = lapply(1:n, function(ind){
      second_step_resp[[ind]]$Ti
    })
    
    qtrue_overall = lapply(1:n, function(ind){
      second_step_resp[[ind]]$q_on_geod
    })
    qin_overall = lapply(1:n, function(ind){
      second_step_resp[[ind]]$qin
    })
    
    second_step_resp0 = t(sapply(1:n, function(ind){
      second_step_resp[[ind]]$mi_hat0
    }))
    second_step_resp1 = t(sapply(1:n, function(ind){
      second_step_resp[[ind]]$mi_hat1
    }))
    ##
    qSup = seq(1e-10, 1-1e-10, length.out = 100) #support of quantile func
    qSup[1] = 0
    qSup[length(qSup)] = 1
    zout = Z #seq(rangeZ[1],rangeZ[2], length.out = 5)
    
    second_step_gfr_fit0 = GloDenReg(xin = Z, qin = second_step_resp0,
                                     xout = zout, optns = list(qSup = qSup))$qout
    
    second_step_gfr_fit1 = GloDenReg(xin = Z, qin = second_step_resp1,
                                     xout = zout, optns = list(qSup = qSup))$qout
    ## fitted geodesic
    outputs =  lapply(1:n, function(ind){
      t(sapply(1:ni[ind], function(j){
        q_geo(Ti[[ind]][j], second_step_gfr_fit0[ind,], second_step_gfr_fit1[ind,])
      }))
    })
    measureWfits = mean(sapply(1:n, function(ind){
      mean(sqrt(sapply(1:ni[ind], function(j){
        sum((outputs[[ind]][j,] -  qtrue_overall[[ind]][j,])^2* diff(qSup)[1])
      })))
    }))
  }, mc.cores = ncore)
  save(res, file = sprintf("../data/Setting1_MISE_dense_n%d.Rda",n))
}
#################################
#################################
## Reading the saved .Rda files to generate Figure 3
B = 150
################dense design
###dense setting 1
load("../data/MISE_density_boxplot_comp_dense_n50_setting1.Rda")
res_n50_dense_set1 = sapply(1:B, function(i) res[[i]]$measureWfits)
load("../data/MISE_density_boxplot_comp_dense_n400_setting1.Rda")
res_n100_dense_set1 = sapply(1:B, function(i) res[[i]]$measureWfits)
load("../data/MISE_density_boxplot_comp_dense_n1000_setting1.Rda")
res_n50_dense_set1 = sapply(1:B, function(i) res[[i]]$measureWfits)

###dense setting 2
load("../data/MISE_density_boxplot_comp_dense_n50_setting2.Rda")
res_n50_dense_set2 = sapply(1:B, function(i) res[[i]]$measureWfits)
load("../data/MISE_density_boxplot_comp_dense_n400_setting2.Rda")
res_n100_dense_set2 = sapply(1:B, function(i) res[[i]]$measureWfits)
load("../data/MISE_density_boxplot_comp_dense_n1000_setting2.Rda")
res_n50_dense_set2 = sapply(1:B, function(i) res[[i]]$measureWfits)

###dense setting 3
load("../data/MISE_density_boxplot_comp_dense_n50_setting3.Rda")
res_n50_dense_set3 = sapply(1:B, function(i) res[[i]]$measureWfits )
load("../data/MISE_density_boxplot_comp_dense_n400_setting3.Rda")
res_n100_dense_set3 = sapply(1:B, function(i) res[[i]]$measureWfits)
load("../data/MISE_density_boxplot_comp_dense_n1000_setting3.Rda")
res_n50_dense_set3 = sapply(1:B, function(i) res[[i]]$measureWfits)

###dense setting 4
load("../data/MISE_density_boxplot_comp_dense_n50_setting4.Rda")
res_n50_dense_set4 = sapply(1:B, function(i) res[[i]]$measureWfits )
load("../data/MISE_density_boxplot_comp_dense_n400_setting4.Rda")
res_n100_dense_set4 = sapply(1:B, function(i) res[[i]]$measureWfits)
load("../data/MISE_density_boxplot_comp_dense_n1000_setting4.Rda")
res_n50_dense_set4 = sapply(1:B, function(i) res[[i]]$measureWfits)


### put together in a dataframe
df_dense_setting1 = data.frame(MISE = c(res_n50_dense_set1,
                                        res_n100_dense_set1,
                                        res_n50_dense_set1))
df_dense_setting1$sample_size = c(rep("n = 50, ni = 30", B),
                                  rep("n = 400, ni = 30", B),
                                  rep("n = 1000, ni = 30", B))

df_dense_setting2 = data.frame(MISE = c(res_n50_dense_set2,
                                        res_n100_dense_set2,
                                        res_n50_dense_set2))
df_dense_setting2$sample_size = c(rep("n = 50, ni = 30", B),
                                  rep("n = 400, ni = 30", B),
                                  rep("n = 1000, ni = 30", B))


df_dense_setting3 = data.frame(MISE = c(res_n50_dense_set3,
                                        res_n100_dense_set3,
                                        res_n50_dense_set3))
df_dense_setting3$sample_size = c(rep("n = 50, ni = 30", B),
                                  rep("n = 400, ni = 30", B),
                                  rep("n = 1000, ni = 30", B))
df_dense_setting4 = data.frame(MISE = c(res_n50_dense_set4,
                                        res_n100_dense_set4,
                                        res_n50_dense_set4))
df_dense_setting4$sample_size = c(rep("n = 50, ni = 30", B),
                                  rep("n = 400, ni = 30", B),
                                  rep("n = 1000, ni = 30", B))

################sparse design
###sparse setting 1
load("../data/MISE_density_boxplot_comp_sparse_n50_setting1.Rda")
res_n50_sparse_set1 = sapply(1:B, function(i) res[[i]]$measureWfits)
load("../data/MISE_density_boxplot_comp_sparse_n400_setting1.Rda")
res_n100_sparse_set1 = sapply(1:B, function(i) res[[i]]$measureWfits)
load("../data/MISE_density_boxplot_comp_sparse_n1000_setting1.Rda")
res_n50_sparse_set1 = sapply(1:B, function(i) res[[i]]$measureWfits)

###sparse setting 2
load("../data/MISE_density_boxplot_comp_sparse_n50_setting2.Rda")
res_n50_sparse_set2 = sapply(1:B, function(i) res[[i]]$measureWfits)
load("../data/MISE_density_boxplot_comp_sparse_n400_setting2.Rda")
res_n100_sparse_set2 = sapply(1:B, function(i) res[[i]]$measureWfits)
load("../data/MISE_density_boxplot_comp_sparse_n1000_setting2.Rda")
res_n50_sparse_set2 = sapply(1:B, function(i) res[[i]]$measureWfits)

###sparse setting 3
load("../data/MISE_density_boxplot_comp_sparse_n50_setting3.Rda")
res_n50_sparse_set3 = sapply(1:B, function(i) res[[i]]$measureWfits)
load("../data/MISE_density_boxplot_comp_sparse_n400_setting3.Rda")
res_n100_sparse_set3 = sapply(1:B, function(i) res[[i]]$measureWfits)
load("../data/MISE_density_boxplot_comp_sparse_n1000_setting3.Rda")
res_n50_sparse_set3 = sapply(1:B, function(i) res[[i]]$measureWfits)

###sparse setting 4
load("../data/MISE_density_boxplot_comp_sparse_n50_setting4.Rda")
res_n50_sparse_set4 = sapply(1:B, function(i) res[[i]]$measureWfits )
load("../data/MISE_density_boxplot_comp_sparse_n400_setting4.Rda")
res_n100_sparse_set4 = sapply(1:B, function(i) res[[i]]$measureWfits)
load("../data/MISE_density_boxplot_comp_sparse_n1000_setting4.Rda")
res_n50_sparse_set4 = sapply(1:B, function(i) res[[i]]$measureWfits)


###
### put together in a dataframe
df_sparse_setting1 = data.frame(MISE = c(res_n50_sparse_set1,
                                         res_n100_sparse_set1,
                                         res_n50_sparse_set1))
df_sparse_setting1$sample_size = c(rep("n = 50, ni = 30", B),
                                   rep("n = 400, ni = 30", B),
                                   rep("n = 1000, ni = 30", B))

df_sparse_setting2 = data.frame(MISE = c(res_n50_sparse_set2,
                                         res_n100_sparse_set2,
                                         res_n50_sparse_set2))
df_sparse_setting2$sample_size = c(rep("n = 50, ni = 30", B),
                                   rep("n = 400, ni = 30", B),
                                   rep("n = 1000, ni = 30", B))


df_sparse_setting3 = data.frame(MISE = c(res_n50_sparse_set3,
                                         res_n100_sparse_set3,
                                         res_n50_sparse_set3))
df_sparse_setting3$sample_size = c(rep("n = 50, ni = 30", B),
                                   rep("n = 400, ni = 30", B),
                                   rep("n = 1000, ni = 30", B))
df_sparse_setting4 = data.frame(MISE = c(res_n50_sparse_set4,
                                         res_n100_sparse_set4,
                                         res_n50_sparse_set4))
df_sparse_setting4$sample_size = c(rep("n = 50, ni = 30", B),
                                   rep("n = 400, ni = 30", B),
                                   rep("n = 1000, ni = 30", B))

############################## Generating Frigure 3
df_overall = rbind(df_dense_setting1,
                   df_dense_setting2,
                   df_dense_setting3,
                   df_dense_setting4,
                   df_sparse_setting1,
                   df_sparse_setting2,
                   df_sparse_setting3,
                   df_sparse_setting4)

df_overall$sample_size2 = rep(rep(c(rep("n = 50",B), rep("n = 400",B), rep("n = 1000",B)),4),2)

df_overall$sample_design = c(rep("Dense design", B*3*4), rep("Sparse design", B*3*4))

df_overall$setting =  rep(c(rep("Setting I", B*3),
                            rep("Setting II", B*3),
                            rep("Setting III", B*3),
                            rep("Setting IV", B*3)),2)
neworder <- c("n = 50","n = 400", "n = 1000")
library(dplyr)
df_overall2 <- arrange(transform(df_overall,
                                 sample_size2 = factor(sample_size2,levels=neworder)),sample_size2)


p = ggplot(df_overall2, aes(sample_size2, MISE)) + 
  geom_boxplot(aes(fill = factor(sample_design))) + 
  #stat_boxplot(geom ='errorbar') +
  labs(x = NULL , y = "ISE") +theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.title = element_blank(),
        legend.position="bottom", strip.text.y = element_blank()) +
  facet_wrap(setting~., scales = "free_y") 

p
ggsave("../output/mise_boxplot_over_n_ni_setting_final.pdf",
       width = 8, height = 7)


######################################
######################################
######## Visualization for an individual -- Generating Figures 1 & 2
######################################
### Fits across different values of alpha -- compute true, observed, 
### fitted trajectories at 3 different time points along the geodesics.
### outputs are saved in the data folder -- running  will take about 5 mins for each alpha value.: used for boxplot
### Go to line 382
n=50
setting = 3
df_comp =  function(alpha){
  Z = sapply(1:n, function(ind)  runif(1, 0, 1)) ## baseline covariate
  ni = sapply(1:n, function(ind) 10+2) ## number of observations per subject
  second_step_resp = lapply(1:n, function(ind){ 
    ##gather the second step responses for all subjects
    subject_level_GFR = subject_level_GFR(ind, Z, alpha, ni, setting = setting)
    qin = subject_level_GFR$qin
    q_on_geod = subject_level_GFR$q_on_geod
    obj_fit_ind = subject_level_GFR$obj_fit
    Ti = subject_level_GFR$Ti
    return(list(q_on_geod = q_on_geod, Ti =Ti,
                qin = qin,
                obj_fit_ind = obj_fit_ind))
  })
  Ti = lapply(1:n, function(ind){
    second_step_resp[[ind]]$Ti
  })
  second_step_resp0 = t(sapply(1:n, function(ind){ ## second step responses at time 0 for all subjects
    second_step_resp[[ind]]$obj_fit_ind[1,]
  }))
  second_step_resp1 = t(sapply(1:n, function(ind){ ## second step responses at time 1 for all subjects
    second_step_resp[[ind]]$obj_fit_ind[2,]
  }))
  qtrue_overall = lapply(1:n, function(ind){ ## Quantile representation of the true geodesic trajectories for all subjects
    second_step_resp[[ind]]$q_on_geod
  })
  qin_overall = lapply(1:n, function(ind){ ## Quantile representation of the observed data around the true geodesic trajectories for all subjects
    second_step_resp[[ind]]$qin
  })
  
  qSup = seq(1e-10, 1-1e-10, length.out = 100) #support of quantile func
  qSup[1] = 0
  qSup[length(qSup)] = 1
  zout = Z
  ### Fitting the second-step regression with the baseline covariates and the tuple responses (quantile representation) for all the subjects
  second_step_gfr_fit0 = GloDenReg(xin = Z, qin = second_step_resp0,
                                   xout = zout, optns = list(qSup = qSup))$qout
  
  second_step_gfr_fit1 = GloDenReg(xin = Z, qin = second_step_resp1,
                                   xout = zout, optns = list(qSup = qSup))$qout
  ## fitted geodesic
  outputs =  lapply(1:n, function(ind){
    t(sapply(1:ni[ind], function(j){
      q_geo(Ti[[ind]][j], second_step_gfr_fit0[ind,],
            second_step_gfr_fit1[ind,])
    }))
  })
  ### converting the quantile functions to density objects for subject 1 for visualization
  ind = 1
  Ti_ind = Ti[[ind]]
  dtrue_overall_ind = lapply(1:ni[ind], function(j){
    frechet:::qf2pdf(qtrue_overall[[ind]][j,], optns = list(qSup = qSup, nRegGrid =100))
  })
  din_overall_ind = lapply(1:ni[ind], function(j){
    frechet:::qf2pdf(qin_overall[[ind]][j,], optns = list(qSup = qSup, nRegGrid =100))
  })
  dout_overall_ind = lapply(1:ni[ind], function(j){
    frechet:::qf2pdf(outputs[[ind]][j,], optns = list(qSup = qSup, nRegGrid =100))
  })
  res = list(dtrue_overall_ind = dtrue_overall_ind,
             din_overall_ind = din_overall_ind,
             dout_overall_ind = dout_overall_ind, Ti = Ti)
  save(res, file = sprintf("../data/subj_sp_over_alpha_A%d_Setting3.Rda", l))
}
alpha = c(0.01,0.1,0.2); rslt = list()
for(l in 1:3){ rslt[[l]] = df_comp(alpha[l])}

##### Generating Figure 1
load("../data/subj_sp_over_alpha_A1_Setting3.Rda") 
res_a0 = res
load("../data/subj_sp_over_alpha_A2_Setting3.Rda")
res_a1 = res
load("../data/subj_sp_over_alpha_A3_Setting3.Rda")
res_a2 = res

### gathering all data for visualization
dat_a = function(res){
  dtrue_overall_ind_display = do.call(rbind, lapply(c(1,5,12), function(j){
    data.frame(domain = res$dtrue_overall_ind[[j]]$x,
               dens = res$dtrue_overall_ind[[j]]$y,
               time = rep(as.character(round(res$Ti[[ind]],2)[j]), 100))
  }))
  
  din_overall_ind_display = do.call(rbind, lapply(c(1,5,12), function(j){
    data.frame(domain = res$din_overall_ind[[j]]$x,
               dens = res$din_overall_ind[[j]]$y,
               time = rep(as.character(round(res$Ti[[ind]],2)[j]), 100))
  }))
  
  dout_overall_ind_display = do.call(rbind, lapply(c(1,5,12), function(j){
    data.frame(domain = res$dout_overall_ind[[j]]$x,
               dens = res$dout_overall_ind[[j]]$y,
               time = rep(as.character(round(res$Ti[[ind]],2)[j]), 100))
  }))
  dat = rbind(dtrue_overall_ind_display,din_overall_ind_display,dout_overall_ind_display)
  dat$type = c(rep("dtrue", 100*3), rep("din", 100*3),rep("dout", 100*3))
  return(dat)
}
ind = 1
dat_a0 = dat_a(res_a0); dat_a1 = dat_a(res_a1); dat_a2 = dat_a(res_a2)
dat_a0$time[dat_a0$time == unique(dat_a0$time)[1]] = "0"
dat_a1$time[dat_a1$time == unique(dat_a1$time)[1]] = "0"
dat_a2$time[dat_a2$time == unique(dat_a2$time)[1]] = "0"

dat_a0$time[dat_a0$time == unique(dat_a0$time)[2]] = "0.50"
dat_a1$time[dat_a1$time == unique(dat_a1$time)[2]] = "0.50"
dat_a2$time[dat_a2$time == unique(dat_a2$time)[2]] = "0.50"

dat_a0$time[dat_a0$time == unique(dat_a0$time)[3]] = "1"
dat_a1$time[dat_a1$time == unique(dat_a1$time)[3]] = "1"
dat_a2$time[dat_a2$time == unique(dat_a2$time)[3]] = "1"


dat_combined = rbind(dat_a0,dat_a1,dat_a2)
dat_combined$alpha = c(rep("a0",nrow(dat_a0)), rep("a1",nrow(dat_a1)), rep("a2",nrow(dat_a2)))

alpha.labs = c(a0 = "alpha[1] == 0.01", a1 = "alpha[2] == 0.1", a2 = "alpha[3] == 0.3")
t.labs = c("0" = "t == 0", "0.50" = "t == 0.50", "1" = "t == 1")

library(ggplot2)
cols = c("red", "springgreen4", "blue")
p = ggplot() + 
  geom_line(data = dat_combined, 
            aes(x = domain, y = dens, color = type, group = type), linewidth = .4) + 
  theme_bw() +
  scale_colour_manual(values = cols, name = " ",
                      labels = c("True", "Observed", "Estimated"))+ 
  labs(x = "Domain", y = "Densities") +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size  = 12),
        legend.title = element_text(size  = 12),
        strip.text.x = element_text(size  = 12),
        strip.text.y = element_text(size  = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom") +
  facet_grid(time~alpha, labeller = as_labeller(c(t.labs,alpha.labs), default = label_parsed),
             scales = "free")
p
ggsave("../output/density_true_obs_estd_over_alpha_final.pdf",
       width = 8, height = 7)
######################################
######################################
### Fits over different levels of the baseline covariate Z -- 
### compute observed and fitted trajectories for 4 settings at 3 different time points along the geodesics for a fixed level of error alpha
### outputs are saved in the data folder -- running  will take about 10 mins for each setting value.: used for boxplot
### Do not run: go to line 540
alpha = 0.1
for_each_setting = function(setting){
  Z = sapply(1:n, function(ind)  runif(1, 0, 1)) ## Generate baseline covariate
  ni = sapply(1:n, function(ind) 10+2) ## number of observations for each subject
  second_step_resp = lapply(1:n, function(ind){
    ##gather the second step responses for all subjects
    subject_level_GFR = subject_level_GFR(ind, Z, alpha, ni, setting = setting)
    qin = subject_level_GFR$qin
    q_on_geod = subject_level_GFR$q_on_geod
    obj_fit_ind = subject_level_GFR$obj_fit
    Ti = subject_level_GFR$Ti
    return(list(q_on_geod = q_on_geod, Ti =Ti,
                qin = qin,
                obj_fit_ind = obj_fit_ind))
  })
  Ti = lapply(1:n, function(ind){
    second_step_resp[[ind]]$Ti
  })
  second_step_resp0 = t(sapply(1:n, function(ind){ ## second step responses at time 0 for all subjects
    second_step_resp[[ind]]$obj_fit_ind[1,]
  }))
  second_step_resp1 = t(sapply(1:n, function(ind){ ## second step responses at time 1 for all subjects
    second_step_resp[[ind]]$obj_fit_ind[2,]
  }))
  z_vals = sort(Z)[c(6,25,45)] ### output levels for the baseline covariate
  qSup = seq(1e-10, 1-1e-10, length.out = 100) #support of quantile functions
  qSup[1] = 0
  qSup[length(qSup)] = 1
  kk = 1 ## initialize counter for each setting (to be saved as .Rda files)
  res = list()
  for(zout in z_vals){
    ### Fitting the second-step regression with the baseline covariates and the tuple responses (quantile representation) for all the subjects
    ### at the three output points
    second_step_gfr_fit0 = GloDenReg(xin = Z, qin = second_step_resp0,
                                     xout = zout, optns = list(qSup = qSup))$qout
    second_step_gfr_fit1 = GloDenReg(xin = Z, qin = second_step_resp1,
                                     xout = zout, optns = list(qSup = qSup))$qout
    ### which unit corresponds to the output z-value -- 
    ### only compute the observed and fitted distribution for that unit
    ind = which(Z == zout)
    ### fitted geodesic
    outputs = t(sapply(1:ni[ind], function(j){
      q_geo(Ti[[ind]][j], as.vector(second_step_gfr_fit0), as.vector(second_step_gfr_fit1))
    }))
    qtrue_z = second_step_resp[[ind]]$q_on_geod ## Quantile representation of the true geodesic trajectories for the subject index "ind"
    qin_z = second_step_resp[[ind]]$qin ## Quantile representation of the observed data around the true geodesic trajectories for the subject index "ind"
    ### converting the quantile functions to density objects for subject 1 for visualization
    dout_z = lapply(1:ni[ind], function(j){
      frechet:::qf2pdf(outputs[j,], optns = list(qSup = qSup, nRegGrid =100))
    })
    dtrue_z = lapply(1:ni[ind], function(j){
      frechet:::qf2pdf(qtrue_z[j,], optns = list(qSup = qSup, nRegGrid =100))
    })
    din_z = lapply(1:ni[ind], function(j){
      frechet:::qf2pdf(qin_z[j,], optns = list(qSup = qSup, nRegGrid =100))
    })
    res[[kk]] = list(T_ind_zz = Ti[[ind]], dtrue_zz = dtrue_z, dout_zz = dout_z, dtrue_zz2 = din_z)
    kk = kk + 1
  }
 save(res, file = sprintf("../data/subj_specific_dens_zz_setting%d.Rda", setting))
}
#### Do not run -- will take more than half an hour to run for all settings 
#### Go to line 540
dens_display_setting1 = for_each_setting(1)
dens_display_setting2 = for_each_setting(2)
dens_display_setting3 = for_each_setting(3)
dens_display_setting4 = for_each_setting(4)


df_dens_out = function(res){
  Ti_ind0 = res[[1]]$T_ind_zz; Ti_ind1 = res[[2]]$T_ind_zz; Ti_ind2 = res[[3]]$T_ind_zz
  dtrue_zz = list(); dtrue_zz[[1]] = res[[1]]$dtrue_zz2; dtrue_zz[[2]] = res[[2]]$dtrue_zz2; dtrue_zz[[3]] = res[[3]]$dtrue_zz2
  ### gathering the true geodesic at different output z-levels at 3 different time points
  dtrue_z1_overall_ind_display = do.call(rbind, lapply(c(1,5,12), function(j){
    data.frame(domain = dtrue_zz[[1]][[j]]$x,
               dens = dtrue_zz[[1]][[j]]$y)
  }))
  dtrue_z2_overall_ind_display = do.call(rbind, lapply(c(1,5,12), function(j){
    data.frame(domain = dtrue_zz[[2]][[j]]$x,
               dens = dtrue_zz[[2]][[j]]$y)
  }))
  dtrue_z3_overall_ind_display = do.call(rbind, lapply(c(1,5,12), function(j){
    data.frame(domain = dtrue_zz[[3]][[j]]$x,
               dens = dtrue_zz[[3]][[j]]$y)
  }))
  dout_zz = list(); dout_zz[[1]] = res[[1]]$dtrue_zz; dout_zz[[2]] = res[[2]]$dtrue_zz; dout_zz[[3]] = res[[3]]$dtrue_zz
  ### gathering the true geodesic at different output z-levels at 3 different time points
  dout_z1_overall_ind_display = do.call(rbind, lapply(c(1,5,12), function(j){
    data.frame(domain = dout_zz[[1]][[j]]$x,
               dens = dout_zz[[1]][[j]]$y) 
    }))
  dout_z2_overall_ind_display = do.call(rbind, lapply(c(1,5,12), function(j){
    data.frame(domain = dout_zz[[2]][[j]]$x,
               dens = dout_zz[[2]][[j]]$y)
  }))
  dout_z3_overall_ind_display = do.call(rbind, lapply(c(1,5,12), function(j){
    data.frame(domain = dout_zz[[3]][[j]]$x,
               dens = dout_zz[[3]][[j]]$y)
  }))
  dens_display_out = rbind(dtrue_z1_overall_ind_display, dout_z1_overall_ind_display,
                           dtrue_z2_overall_ind_display, dout_z2_overall_ind_display,
                           dtrue_z3_overall_ind_display, dout_z3_overall_ind_display)
  
  dens_display_out$type = c(rep("dtrue", 100*3), rep("dout", 100*3),
                            rep("dtrue", 100*3), rep("dout", 100*3), 
                            rep("dtrue", 100*3), rep("dout", 100*3))
  
  dens_display_out$z = c(rep("z1", 100*3), rep("z1", 100*3),
                         rep("z2", 100*3), rep("z2", 100*3), 
                         rep("z3", 100*3), rep("z3", 100*3))
  dens_display_out$time = rep(c(rep("t = 0", 100), rep("t = 0.50", 100), rep("t = 1", 100)),6)
  return(dens_display_out)
}
##################
load("../data/subj_specific_dens_zz_setting1.Rda")
dens_display_setting1 = df_dens_out(res)
load('../data/subj_specific_dens_zz_setting2.Rda')
dens_display_setting2 = df_dens_out(res)
load('../data/subj_specific_dens_zz_setting3.Rda')
dens_display_setting3 = df_dens_out(res)
load('../data/subj_specific_dens_zz_setting4.Rda')
dens_display_setting4 = df_dens_out(res)

dens_display_final = rbind(dens_display_setting1, dens_display_setting2, dens_display_setting3, dens_display_setting4)
dens_display_final$setting = c(rep("Setting I",nrow(dens_display_setting1)),
                               rep("Setting II",nrow(dens_display_setting2)),
                               rep("Setting III",nrow(dens_display_setting3)),
                               rep("Setting IV",nrow(dens_display_setting4)))

library(ggplot2)
cols = c("red","blue")
p2 = ggplot() + 
  geom_line(data = dens_display_final, 
            aes(x = domain, y = dens, color = type, linetype = z,
                group = interaction(type,z)), linewidth = .4) +
  theme_bw() +
  scale_colour_manual(values = cols, name = " ",
                      labels = c("True", "Estimated"))+ 
  scale_linetype_manual(values = c("dotted", "longdash","solid"), name = " ",
                        labels = c("Estimated at 90%\n quantile of Z",
                                   "Estimated at 50%\n quantile of Z",
                                   "Estimated at 10%\n quantile of Z"))+ 
  labs(x = "Domain", y = "Densities") + ylim(0,2) +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45),
        axis.title = element_text(size  = 12),
        legend.title = element_text(size  = 12),
        strip.text = element_text(size = 25, margin = margin()),
        strip.text.x = element_text(size  = 12),
        strip.text.y = element_text(size  = 12),
        legend.position = "bottom") +
  facet_grid(time ~ setting, scales = "free")
p2
ggsave("../output/density_over_zz_final.pdf", width = 8, height = 7)

