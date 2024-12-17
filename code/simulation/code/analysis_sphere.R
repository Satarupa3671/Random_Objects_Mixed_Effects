#rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(frechet)
library(parallel)
library(MASS)

### load auxiliary functions
source("./SpheGeoDist.R") 
source("./SpheGeoGrad.R")
source("./SpheGeoHess.R")
source("./l2norm.R")
source("./GloSpheGeoReg.R")

##Data Generation::  Geodesic on the surface of a sphere
generate_pts_on_geod = function(t,A,B){
  omega = as.numeric(acos(A %*% B))
  d = sin(omega)
  s0 = sin((1-t)*omega)
  s1 = sin(t*omega)
  return((A*s0 + B*s1)/d)
}
##Data Generation::  subject level geodesic : generate observed data and recovered by global Fr regression
subject_level_GFR = function(ind, Z, alpha){
  ni = 30#sample(5,1,replace = TRUE)
  Ti = runif(ni)
  err_sd <- 0.1
  phi_true <- acos(Z[ind])
  theta_true <- pi * (Z[ind])
  ytrue <- cbind(sin(phi_true) * cos(theta_true),sin(phi_true) * sin(theta_true),cos(phi_true))
  basis <- list(b1 = cbind(cos(phi_true) * cos(theta_true),cos(phi_true) * sin(theta_true),-sin(phi_true)),
                b2 = cbind(sin(theta_true),-cos(theta_true),0))
  yin = matrix(0,nrow = 2,ncol =3)
  for( k in 1:2){
    yin_tg <- basis$b1 * rnorm(1, mean = 0, sd = err_sd) + basis$b2 * rnorm(1, mean = 0, sd = err_sd)
    tgNorm <- sqrt(sum(yin_tg^2))
    if (tgNorm < 1e-10) {
      yin[k,] = ytrue
    } else {
      yin[k,] = sin(tgNorm) * yin_tg / tgNorm + cos(tgNorm) * ytrue
    }
  }
  ## Another way to generate the end point of subject specific geodesic trajectories
  # Z_index = Z[ind,]
  # mu_givenZ = rnorm(2,Z_index,.1)
  # #mu_givenZ = c(-1,1)
  # sigma_givenZ = c(.1,.1) 
  # A = rnorm(3, mean = 0 + mu_givenZ[1], sd = sigma_givenZ[1])
  # B = rnorm(3, mean = 1 + mu_givenZ[2], sd = sigma_givenZ[2])
  # A = A/sqrt(sum(A^2))#yin[1,]
  # B = B/sqrt(sum(B^2))#yin[2,]
  A = yin[1,]
  B = yin[2,]
  err_g = rnorm(3,0,1)
  inputs = lapply(1:ni, function(j){
    pts_on_geod = generate_pts_on_geod(Ti[j],A,B) 
    pts_noisy = (1-alpha)*pts_on_geod + alpha*err_g
    pts_noisy = pts_noisy/sqrt(sum(pts_noisy^2))
    return(list(pts_noisy = pts_noisy, pts_on_geod = pts_on_geod))
  })
  pts_on_geod = t(sapply(1:ni, function(j) inputs[[j]]$pts_on_geod))
  pts_on_geod = rbind(A, pts_on_geod, B)
  pts_noisy = t(sapply(1:ni, function(j) inputs[[j]]$pts_noisy))
  AA = ((1 - alpha)*A + alpha*err_g); AA = AA/sqrt(sum(AA^2))
  BB = ((1 - alpha)*B + alpha*err_g); BB = BB/sqrt(sum(BB^2))
  pts_noisy = rbind(AA, pts_noisy, BB)
  Ti = c(0,Ti,1)
  xout = c(0,1)
  obj_fit = t(sapply(xout, function(xx) GloSpheGeoReg(xin = as.matrix(Ti), yin = pts_noisy, xout  = xx)))
  
  return(list(Ti = Ti, pts_on_geod = pts_on_geod, 
              yin = pts_noisy, obj_fit = obj_fit))
}
### Do not run : compute MISE over 200 replications through parallel computation
### outputs are saved in the data folder : used for boxplot
### Go to line 138 
ncore = 20
B = 200
for(n in c(50,400,1000)){
  res = mclapply(1:B, function(rep){
    alpha = .1
    Z = as.matrix(rbeta(n,4,2), nrow = n)
    second_step_resp = lapply(1:n, function(ind){
      subject_level_GFR = subject_level_GFR(ind, Z, alpha)
      obj_fit_ind = subject_level_GFR$obj_fit
      mi_hat0 = obj_fit_ind[1,]
      mi_hat1 = obj_fit_ind[2,]
      pts_on_geod =  subject_level_GFR$pts_on_geod
      yin = subject_level_GFR$yin
      Ti = subject_level_GFR$Ti
      return(list(pts_on_geod =  pts_on_geod, Ti = Ti,
                  yin = yin, mi_hat0 = mi_hat0, mi_hat1 = mi_hat1))
    })
    Ti = lapply(1:n, function(ind){
      second_step_resp[[ind]]$Ti
    })
    
    pts_true_overall = lapply(1:n, function(ind){
      second_step_resp[[ind]]$pts_on_geod
    })
    
    second_step_resp0 = t(sapply(1:n, function(ind){
      second_step_resp[[ind]]$mi_hat0
    }))
    second_step_resp1 = t(sapply(1:n, function(ind){
      second_step_resp[[ind]]$mi_hat1
    }))
    ##
    zout = Z
    
    second_step_gfr_fit0 = t(sapply(zout, function(zz) GloSpheGeoReg(xin = as.matrix(Z),
                                                                     yin = second_step_resp0, xout  = zz)))
    
    second_step_gfr_fit1 = t(sapply(zout, function(zz) GloSpheGeoReg(xin = as.matrix(Z),
                                                                     yin = second_step_resp1, xout  = zz)))
    
    ## fitted geodesic
    ni = sapply(Ti, function(tt) length(tt))
    outputs =  lapply(1:n, function(ind){
      out = t(sapply(1:ni[ind], function(j){
        generate_pts_on_geod(Ti[[ind]][j], second_step_gfr_fit0[ind,], second_step_gfr_fit1[ind,])
      }))
      # out = rbind(second_step_gfr_fit0[ind,], out, second_step_gfr_fit1[ind,])
    })
    
    measureWfits = mean(sapply(1:n, function(ind){
      mean(sapply(1:ni[ind], function(j){
        sgn = sign(outputs[[ind]][j,] %*% pts_true_overall[[ind]][j,])
        acos((outputs[[ind]][j,] %*% pts_true_overall[[ind]][j,])*sgn)
      }))
    }))
  }, mc.cores = ncore)
  save(res, file = sprintf("sphe_MISE_dense_n%d.Rda",n))
  # save(res, file = sprintf("sphe_MISE_sparse_n%d.Rda",n))
}
#########################################
#########################################
### Visulalization -- Generate Figure S.1.
###sparse setting mu changes with pred, sigma fixed
load("../data/sphe_MISE_sparse_n50.Rda")
dat1 = unlist(res)
load("../data/sphe_MISE_sparse_n400.Rda")
dat2 = unlist(res)
load("../data/sphe_MISE_sparse_n1000.Rda")
dat3 = unlist(res)
B = 200
sparse_MISE = data.frame(MISE = c(dat1,dat2,dat3))
sparse_MISE$sample_size = c(rep("n = 50", B),
                            rep("n = 400", B),
                            rep("n = 1000", B))

####summary of dense design, mu sig changes with Z
summary(dat1) # 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.06392 0.07721 0.08258 0.08439 0.09065 0.12528 
summary(dat2) #
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.07363 0.07801 0.07994 0.07989 0.08164 0.08706 summary(dat3) #
summary(dat3) #
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.07474 0.07859 0.07973 0.07973 0.08108 0.08393 ###
######
###dense setting mu changes with pred, sigma fixed
load("../data/sphe_MISE_dense_n50.Rda")
dat1 = unlist(res)
load("../data/sphe_MISE_dense_n400.Rda")
dat2 = unlist(res)
load("../data/sphe_MISE_dense_n1000.Rda")
dat3 = unlist(res)
dense_MISE = data.frame(MISE = c(dat1,dat2,dat3))
dense_MISE$sample_size = c(rep("n = 50", B),
                           rep("n = 400", B),
                           rep("n = 1000", B))

####summary of dense design, mu sig changes with Z
summary(dat1) # 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.05696 0.07065 0.07659 0.07778 0.08386 0.12687 
summary(dat2) #
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.06650 0.07258 0.07456 0.07456 0.07595 0.08494 
summary(dat3) #
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.07072 0.07311 0.07407 0.07413 0.07524 0.07855 
#####


###########################################################
########## Visualization of fits --  Generating Figure S.1
df_combined = rbind(dense_MISE,sparse_MISE)
df_combined$sample_design = c(rep("Dense design", 3*B),rep("Sparse design", 3*B))
library(dplyr)
library(ggplot2)
p = ggplot(transform(df_combined,
                     sample_size =factor(sample_size ,levels=c("n = 50","n = 400","n = 1000"))),
           aes(y = MISE, fill = sample_design)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot()+
  facet_wrap(sample_size~.) +
  labs(x = NULL , y = "ISE") +theme_bw()+
  theme(axis.text=element_text(size=12), axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=12),
        legend.title = element_blank(),
        legend.position="bottom", strip.text.y = element_blank()) 
p
dir.create("../output", showWarnings = FALSE)
ggsave("../output/sphe_MISE_boxplot.pdf", width = 8, height = 7)
