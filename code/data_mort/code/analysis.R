#rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(frechet)
####Analysis of Human Mortality Data for univaraite probability distributions as random objects in the space of distributions endowed with the Wasserstein-2 metric
####and generation of tables are figures used in Section S.4.2
###Step 1###Additional processing and fitting first step regression
#### Time consuming -- outputs saved as Rda files --  
###skip running this part go to line 73
df = read.csv("../data/lt_subset.csv",header = TRUE) ###life tables for countries over the years from human mortality database (processed)
Z = read.csv("../data/covariate_1990.csv",header = TRUE) ###baseline covariates from world bank data (processed)

countries = unique(df$Country)
yearStart = 1990
yearEnd = 2019
years = yearStart:yearEnd
n = length(countries)
############
##function to produce density and quantile for each country and at each year
data_density <- function(country,year){
  data_sub = subset(df, (Year == year) & (Country == country))
  lx_year <- as.numeric(as.vector(data_sub$lx))
  binY <- seq(75,110,length.out = length(lx_year)+1)
  density <- frechet::CreateDensity(freq = lx_year,  bin=binY,
                                    optns = list(nRegGrid = 500))
  quantile <- fdadensity::dens2quantile(dens = density$y, dSup = density$x, 
                                        qSup = seq(0,1,length.out = 500))
  return(list(dens = density, quant = quantile))
}
obs = lapply(years, function(year){
  lapply(countries, function(country){
    data_density(country,year) 
  })
})
save(obs, file = "../data/obs_dens_quant.Rda")

############
##function to fit the subject level (first step) Global Fr Regression to compute the random effects at the end points of the time domain
subject_level_GFR = function(ind){
  dens_quant_country_ind =lapply(years, function(year){
    country = countries[ind]
    data_density(country,year)
  })
  qin = t(sapply(1:length(years), function(year){
    dens_quant_country_ind[[year]]$quant
  }))
  dSup =  t(sapply(1:length(years), function(year){
    dens_quant_country_ind[[year]]$dens$x
  }))
  dSup_mean = colMeans(dSup)
  
  xin = (years - yearStart)/ (yearEnd - yearStart)
  obj_fit = frechet::GloDenReg(xin = xin, qin = qin, xout = c(0,1),
                               optns = list(qSup = seq(0,1,length.out = 500),
                                            dSup = dSup_mean,
                                            lower = 72, upper = 120))$qout
  return(list(qin = qin, obj_fit = obj_fit, dSup_mean = dSup_mean))
}

second_step_resp = lapply(1:n, function(ind){
  subject_level_GFR = subject_level_GFR(ind)
  qin = subject_level_GFR$qin
  dSup_mean = subject_level_GFR$dSup_mean
  obj_fit_ind = subject_level_GFR$obj_fit
  mi_hat0 = obj_fit_ind[1,]
  mi_hat1 = obj_fit_ind[2,]
  return(list(qin = qin, mi_hat0 = mi_hat0, mi_hat1 = mi_hat1, 
              dSup_mean = dSup_mean))
})
save(second_step_resp, file = "../data/second_step_resp_2.Rda")
###################
###################
## Load the necessary files from above
load("../data/second_step_resp.Rda")
df = read.csv("../data/lt_subset.csv",header = TRUE)
Z = read.csv("../data/covariate_1990.csv",header = TRUE)

countries = unique(df$Country)
yearStart = 1990
yearEnd = 2019
years = yearStart:yearEnd
n = length(countries)
qin =  lapply(1:n, function(ind){ ###input quantile functions for each country and for each year
  second_step_resp[[ind]]$qin
})

second_step_resp0 = t(sapply(1:n, function(ind){ ###first step estimate of random effects at time 0 -- quantile representation
  second_step_resp[[ind]]$mi_hat0
}))
second_step_resp1 = t(sapply(1:n, function(ind){ ###first step estimate of random effects at time 1 -- quantile representation
  second_step_resp[[ind]]$mi_hat1
}))


qSup = seq(0, 1, length.out = 500) #support of quantile function
Z = as.matrix(Z)
q_geo = function(t,A,B){ ##any point on the geodesic
  t*A + (1-t)*B
}
ni = length(years)
Ti = (years -  years[1])/(years[length(years)] - years[1]) ##standarized time domain

############
## Effects of the baseline covariates
############
######compare pred at the two end years
##first pred - unemployment -- Generating Figure S.5
zout1 = as.matrix(cbind(sort(Z[,1]),
                        mean(Z[,2]),
                        mean(Z[,3]),
                        mean(Z[,4])))
### Fitting the final estimate at varying values of unemployment rate, while all the other predictors are fixed at their mean levels
second_step_gfr_fitz10 = frechet::GloDenReg(xin = Z, qin = second_step_resp0, 
                                            xout =  zout1,
                                            optns = list(qSup = qSup, 
                                                         ndSup = 500))

second_step_gfr_fitz11 = frechet::GloDenReg(xin = Z, qin = second_step_resp1, 
                                            xout = zout1,
                                            optns = list(qSup = qSup, 
                                                         ndSup = 500))

### Generating Figure S.5
df10 =  lapply(1:n, function(i){
  data.frame( domain = second_step_gfr_fitz10$dout[[i]]$x, 
              dens = second_step_gfr_fitz10$dout[[i]]$y, 
              gdp_rank = rep(round(zout1[i,1],2), 500))
})
df10 = do.call(rbind, df10)

df11 =  lapply(1:n, function(i){
  data.frame( domain = second_step_gfr_fitz11$dout[[i]]$x, 
              dens = second_step_gfr_fitz11$dout[[i]]$y, 
              gdp_rank = rep(round(zout1[i,1],2), 500))
})
df11 = do.call(rbind, df11)
df_pred1 = data.frame(rbind(df10,df11))
df_pred1$year = c(rep("Year 1990",500*n),rep("Year 2019",500*n))
df_pred1$pred = rep(rep("Unemployment",each = 500),2)
pred1 = round(sort(Z[,1]),2)

p1 = ggplot(data = df_pred1) + 
  geom_line(aes(x = domain, y = dens, col = gdp_rank, group = gdp_rank), 
            linewidth = .5) +
  theme_bw() +
  scale_color_gradient( low="blue", high="red", name = "Unemployment\n rate",
                        breaks = seq(min(pred1),max(pred1),length.out = 5), 
                        labels=summary(pred1)[-4])+
  labs(y = "Remaining life densities", x = "Age") + 
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.y = element_blank(),
        legend.position="bottom", legend.text = element_text(angle = 90)) +
  facet_wrap(~year)
p1
dir.create("../output", showWarnings = FALSE)
ggsave("../output/change_pred_unemployment_year01_final2.pdf",
       width = 8, height = 7)

############
##second pred - fertility
zout2 = as.matrix(cbind(mean(Z[,1]),
                        sort(Z[,2]),
                        mean(Z[,3]),
                        mean(Z[,4])))
### Fitting the final estimate at varying values of fertility rate, while all the other predictors are fixed at their mean levels
second_step_gfr_fitz20 = frechet::GloDenReg(xin = Z, qin = second_step_resp0, 
                                            xout =  zout2,
                                            optns = list(qSup = qSup, 
                                                         ndSup = 500))

second_step_gfr_fitz21 = frechet::GloDenReg(xin = Z, qin = second_step_resp1, 
                                            xout = zout2,
                                            optns = list(qSup = qSup,
                                                         ndSup = 500))

### Generating Figure S.3
df20 =  lapply(1:n, function(i){
  data.frame( domain = second_step_gfr_fitz20$dout[[i]]$x, 
              dens = second_step_gfr_fitz20$dout[[i]]$y, 
              gdp_rank = rep(round(zout2[i,2],2), 500))
})
df20 = do.call(rbind, df20)

df21 =  lapply(1:n, function(i){
  data.frame( domain = second_step_gfr_fitz21$dout[[i]]$x, 
              dens = second_step_gfr_fitz21$dout[[i]]$y, 
              gdp_rank = rep(round(zout2[i,2],2), 500))
})
df21 = do.call(rbind, df21)
df_pred2 = data.frame(rbind(df20,df21))
df_pred2$year = c(rep("Year 1990",500*n),rep("Year 2019",500*n))
df_pred2$pred = rep(rep("Fertility rate",each = 500),2)
pred2 = round(sort(Z[,2]),2)
p2 = ggplot(data = df_pred2) + 
  geom_line(aes(x = domain, y = dens, col = gdp_rank, group = gdp_rank), 
            linewidth = .5) +
  theme_bw() +
  scale_color_gradient( low="blue", high="red", name = "Fertility rate",
                        breaks = seq(min(pred2),max(pred2),length.out = 5), 
                        labels=summary(pred2)[-4])+
  labs(y = "Remaining life densities", x = "Age") + 
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.y = element_blank(),
        legend.text = element_text(angle = 90),
        legend.position="bottom", legend.direction="horizontal") +
  facet_wrap(~year)
p2
ggsave("../output/change_pred_Fertility rate_year01_final2.pdf",
       width = 8, height = 7)
######
##third pred - gdp per capita
zout3 = as.matrix(cbind(mean(Z[,1]),
                        mean(Z[,2]),
                        sort(Z[,3]),
                        mean(Z[,4])))
### Fitting the final estimate at varying values of GDP per capita, while all the other predictors are fixed at their mean levels
second_step_gfr_fitz30 = frechet::GloDenReg(xin = Z, qin = second_step_resp0, 
                                            xout =  zout3,
                                            optns = list(qSup = qSup, 
                                                         ndSup = 500))

second_step_gfr_fitz31 = frechet::GloDenReg(xin = Z, qin = second_step_resp1, 
                                            xout = zout3,
                                            optns = list(qSup = qSup,
                                                         ndSup = 500))
### Generating Figure S.2
df30 =  lapply(1:n, function(i){
  data.frame( domain = second_step_gfr_fitz30$dout[[i]]$x, 
              dens = second_step_gfr_fitz30$dout[[i]]$y, 
              gdp_rank = rep(zout3[i,3], 500))
})
df30 = do.call(rbind, df30)

df31 =  lapply(1:n, function(i){
  data.frame( domain = second_step_gfr_fitz31$dout[[i]]$x, 
              dens = second_step_gfr_fitz31$dout[[i]]$y, 
              gdp_rank = rep(zout3[i,3], 500))
})
df31 = do.call(rbind, df31)
df_pred3 = data.frame(rbind(df30,df31))
df_pred3$year = c(rep("Year 1990",500*n),rep("Year 2019",500*n))
df_pred3$pred = rep(rep("GDP per capita",each = 500),2)
pred3 = round(sort(Z[,3]),2)
p3 = ggplot(data = df_pred3) + 
  geom_line(aes(x = domain, y = dens, col = gdp_rank, group = gdp_rank), 
            linewidth = .5) +
  theme_bw() +
  scale_color_gradient( low="blue", high="red", name = "GDP\n per capita",
                        breaks = seq(min(pred3),max(pred3),length.out = 5), 
                        labels=round(summary(pred3)[-4],2))+
  labs(y = "Remaining life densities", x = "Age") + 
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.y = element_blank(),
        legend.text = element_text(angle = 90),
        legend.position="bottom", legend.direction="horizontal") +
  facet_wrap(~year)
p3
ggsave("../output/change_pred_GDP_year01_final2.pdf", width = 8, height = 7)
###########
##fourth pred - pop growth
zout4 = as.matrix(cbind(mean(Z[,1]),
                        mean(Z[,2]),
                        mean(Z[,3]),
                        sort(Z[,4])))
### Fitting the final estimate at varying values of population growth, while all the other predictors are fixed at their mean levels
second_step_gfr_fitz40 = frechet::GloDenReg(xin = Z, qin = second_step_resp0, 
                                            xout =  zout4,
                                            optns = list(qSup = qSup, 
                                                         ndSup = 500))

second_step_gfr_fitz41 = frechet::GloDenReg(xin = Z, qin = second_step_resp1, 
                                            xout = zout4,
                                            optns = list(qSup = qSup, 
                                                         ndSup = 500))

### Generating Figure S.4
df40 =  lapply(1:n, function(i){
  data.frame( domain = second_step_gfr_fitz40$dout[[i]]$x, 
              dens = second_step_gfr_fitz40$dout[[i]]$y, 
              gdp_rank = rep(zout4[i,4], 500))
})
df40 = do.call(rbind, df40)

df41 =  lapply(1:n, function(i){
  data.frame( domain = second_step_gfr_fitz41$dout[[i]]$x, 
              dens = second_step_gfr_fitz41$dout[[i]]$y, 
              gdp_rank = rep(zout4[i,4], 500))
})
df41 = do.call(rbind, df41)
df_pred4 = data.frame(rbind(df40,df41))
df_pred4$year = c(rep("Year 1990",500*n),rep("Year 2019",500*n))
df_pred4$pred = rep(rep("Population Growth",each = 500),2)
pred4 = round(sort(Z[,4]),2)
p4 = ggplot(data = df_pred4) + 
  geom_line(aes(x = domain, y = dens, col = gdp_rank, group = gdp_rank), 
            linewidth = .5) +
  theme_bw() +
  scale_color_gradient( low="blue", high="red", name = "Population growth %",
                        breaks = seq(min(pred4),max(pred4),length.out = 5), 
                        labels=summary(pred4)[-4])+
  labs(y = "Remaining life densities", x = "Age") + 
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.y = element_blank(),
        legend.position="bottom", legend.direction="horizontal",
        legend.text = element_text(angle = 90)) +
  facet_wrap(~year)
p4
ggsave("../output/change_pred_PopGrowth_year01_final2.pdf", width = 8, height = 7)
####################################
############
## Final Estimates over different points of the time domain (on the geodesic) for varying values of the baseline covariates
############

###################################################
zout1 = cbind(quantile(Z[,1],c(.1,.5,.9)), ### 10%, 50%, 90% quantiles of unemployment rate
              apply(Z[,2:4],2, quantile, c(.5,.5,.5)))

zout2 = cbind(quantile(Z[,1],c(.50,.5,.50)),quantile(Z[,2],c(.1,.5,.9)),
              apply(Z[,3:4],2,quantile,c(.5,.5,.5)))  ### 10%, 50%, 90% quantiles of fertility rate
zout3 = cbind(apply(Z[,1:2], 2, quantile,c(.50,.5,.50)),
              quantile(Z[,3],c(.1,.5,.9)),  ### 10%, 50%, 90% quantiles of GDP per capita
              quantile(Z[,4],c(.50,.5,.50)))
zout4 = cbind(apply(Z[,1:3], 2, quantile,c(.50,.5,.50)),
              quantile(Z[,4],c(.1,.5,.9)))  ### 10%, 50%, 90% quantiles of population growth
###############################################
##varying levels of unemployment rate over the years -- computing final estimates
second_step_gfr_fitz10 = frechet::GloDenReg(xin = Z, qin = second_step_resp0, 
                                            xout =  zout3, 
                                            optns = list(qSup = qSup))
second_step_gfr_fitz11 = frechet::GloDenReg(xin = Z, qin = second_step_resp1, 
                                            xout =  zout3, 
                                            optns = list(qSup = qSup))

##evaluated at different points on the time domain 

###quantile representation of the final outputs
output_q_z1 =  lapply(1:nrow(zout1), function(ind){
  out = t(sapply(1:ni, function(j){
    q_geo(Ti[j], second_step_gfr_fitz10$qout[ind,], 
          second_step_gfr_fitz11$qout[ind,])
  }))
})
###density representation of the final outputs
output_d_z1 = lapply(1:nrow(zout1), function(ind){
  lapply(1:ni, function(j){
    frechet:::qf2pdf(output_q_z1[[ind]][j,],
                     optns = list(qSup = qSup, 
                                  outputGrid = second_step_resp[[1]]$dSup_mean,
                                  nRegGrid =500))
  })
}) 
years_select = c(2015,2005,1995)
years_select_ind = which(years %in% years_select)

######
###generating Figure S.9
dens_display_z1 = do.call(rbind, lapply(1:nrow(zout1), function(ind){
  do.call(rbind,lapply(years_select_ind, function(j){
    data.frame(domain = output_d_z1[[ind]][[j]]$x,
               dens = output_d_z1[[ind]][[j]]$y,
               year = rep(as.character(years[j]),500))
  }))
})) 
dens_display_z1$pred = c(rep("At 10% quantile",500*length(years_select)),
                         rep("At 50% quantile",500*length(years_select)),
                         rep("At 90% quantile",500*length(years_select)))


cols = c("red", "blue","green4")
yr.labs <- c("year == 2015", "year == 2005", "year == 1995")
names(yr.labs) <- c("1995", "2005", "2015")
library("forcats")
p1 = ggplot(data = dens_display_z1) +
  geom_line(aes(x = domain, y = dens, col= pred, group = pred), linewidth = .5) +
  theme_bw() +
  scale_colour_manual(values = cols, name="Unemployment rate") +
  labs(x  = "Age", 
       y = "Remaining life densities") +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.y = element_blank(),
        legend.position = "bottom") +
  facet_grid(~fct_rev(year), labeller = as_labeller(yr.labs, default = label_parsed))
p1
ggsave("../output/change_pred_unemployment_select_years_over_years_final2.pdf",
       width = 8, height = 7)
############################################
##varying levels of fertility rate over the years -- computing final estimates
second_step_gfr_fitz20 = frechet::GloDenReg(xin = Z, qin = second_step_resp0, 
                                            xout =  zout2, 
                                            optns = list(qSup = qSup))
second_step_gfr_fitz21 = frechet::GloDenReg(xin = Z, qin = second_step_resp1, 
                                            xout =  zout2, 
                                            optns = list(qSup = qSup))

##evaluated at different points on the time domain 

###quantile representation of the final outputs
output_q_z2 =  lapply(1:nrow(zout2), function(ind){
  out = t(sapply(1:ni, function(j){
    q_geo(Ti[j], second_step_gfr_fitz20$qout[ind,],
          second_step_gfr_fitz21$qout[ind,])
  }))
})
###density representation of the final outputs
output_d_z2 = lapply(1:nrow(zout2), function(ind){
  lapply(1:ni, function(j){
    frechet:::qf2pdf(output_q_z2[[ind]][j,], 
                     optns = list(qSup = qSup, 
                                  outputGrid = second_step_resp[[1]]$dSup_mean,
                                  nRegGrid =500))
  })
}) 
######
###generating Figure S.7
dens_display_z2 = do.call(rbind, lapply(1:nrow(zout2), function(ind){
  do.call(rbind,lapply(years_select_ind, function(j){
    data.frame(domain = output_d_z2[[ind]][[j]]$x,
               dens = output_d_z2[[ind]][[j]]$y,
               year = rep(as.character(years[j]),500))
  }))
})) 
dens_display_z2$pred = c(rep("At 10% quantile",500*length(years_select)),
                         rep("At 50% quantile",500*length(years_select)),
                         rep("At 90% quantile",500*length(years_select)))


p2 = ggplot(data = dens_display_z2) +
  geom_line(aes(x = domain, y = dens, col= pred, group = pred), linewidth = .5) +
  theme_bw() +
  scale_colour_manual(values = cols, name="Fertility rate") +
  labs(x  = "Age", 
       y = "Remaining life densities") +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "bottom",
        strip.text.y = element_blank()) +
  facet_grid(~fct_rev(year), labeller = as_labeller(yr.labs, default = label_parsed))
p2
ggsave("../output/change_pred_fertility_select_years_over_years_final2.pdf",
       width = 8, height = 7)
############################################
##varying levels of GDP per capita over the years -- computing final estimates
second_step_gfr_fitz30 = frechet::GloDenReg(xin = Z, qin = second_step_resp0, 
                                            xout =  zout3,
                                            optns = list(qSup = qSup))
second_step_gfr_fitz31 = frechet::GloDenReg(xin = Z, qin = second_step_resp1, 
                                            xout =  zout3, 
                                            optns = list(qSup = qSup))
##evaluated at different points on the time domain 

###quantile representation of the final outputs
output_q_z3 =  lapply(1:nrow(zout3), function(ind){
  out = t(sapply(1:ni, function(j){
    q_geo(Ti[j], second_step_gfr_fitz30$qout[ind,], 
          second_step_gfr_fitz31$qout[ind,])
  }))
})
###density representation of the final outputs
output_d_z3 = lapply(1:nrow(zout3), function(ind){
  lapply(1:ni, function(j){
    frechet:::qf2pdf(output_q_z3[[ind]][j,],
                     optns = list(qSup = qSup,
                                  outputGrid = second_step_resp[[1]]$dSup_mean,
                                  nRegGrid =500))
  })
}) 

######
###generating Figure S.6
dens_display_z3 = do.call(rbind, lapply(1:nrow(zout3), function(ind){
  do.call(rbind,lapply(years_select_ind, function(j){
    data.frame(domain = output_d_z3[[ind]][[j]]$x,
               dens = output_d_z3[[ind]][[j]]$y,
               year = rep(as.character(years[j]),500))
  }))
})) 
dens_display_z3$pred = c(rep("At 10% quantile",500*length(years_select)),
                         rep("At 50% quantile",500*length(years_select)),
                         rep("At 90% quantile",500*length(years_select)))


p3 = ggplot(data = dens_display_z3) +
  geom_line(aes(x = domain, y = dens, col= pred, group = pred), linewidth = .5) +
  theme_bw() +
  scale_colour_manual(values = cols, name="GDP per capita") +
  labs(x  = "Age", 
       y = "Remaining life densities") +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "bottom",
        strip.text.y = element_blank()) +
  facet_grid(~fct_rev(year), labeller = as_labeller(yr.labs, default = label_parsed))
p3
ggsave("../output/change_pred_GDP_select_years_over_years_final2.pdf",
       width = 8, height = 7)
########################################
##varying levels of population growth over the years -- computing final estimates

second_step_gfr_fitz40 = frechet::GloDenReg(xin = Z, qin = second_step_resp0, 
                                            xout =  zout4,
                                            optns = list(qSup = qSup))
second_step_gfr_fitz41 = frechet::GloDenReg(xin = Z, qin = second_step_resp1, 
                                            xout =  zout4,
                                            optns = list(qSup = qSup))
##evaluated at different points on the time domain 

###quantile representation of the final outputs
output_q_z4 =  lapply(1:nrow(zout4), function(ind){
  out = t(sapply(1:ni, function(j){
    q_geo(Ti[j], second_step_gfr_fitz40$qout[ind,], 
          second_step_gfr_fitz41$qout[ind,])
  }))
})

###density representation of the final outputs
output_d_z4 = lapply(1:nrow(zout4), function(ind){
  lapply(1:ni, function(j){
    frechet:::qf2pdf(output_q_z4[[ind]][j,], 
                     optns = list(qSup = qSup,
                                  outputGrid = second_step_resp[[1]]$dSup_mean,
                                  nRegGrid =500))
  })
}) 

######
###generating Figure S.8
dens_display_z4 = do.call(rbind, lapply(1:nrow(zout4), function(ind){
  do.call(rbind,lapply(years_select_ind, function(j){
    data.frame(domain = output_d_z4[[ind]][[j]]$x,
               dens = output_d_z4[[ind]][[j]]$y,
               year = rep(as.character(years[j]),500))
  }))
})) 
dens_display_z4$pred = c(rep("At 10% quantile",500*length(years_select)),
                         rep("At 50% quantile",500*length(years_select)),
                         rep("At 90% quantile",500*length(years_select)))


p4 = ggplot(data = dens_display_z4) +
  geom_line(aes(x = domain, y = dens, col= pred, group = pred), linewidth = .5) +
  theme_bw() +
  scale_colour_manual(values = cols, name="Population growth %") +
  labs(x  = "Age", 
       y = "Remaining life densities") +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.y = element_blank(),
        legend.position = "bottom") +
  facet_grid(~fct_rev(year), labeller = as_labeller(yr.labs, default = label_parsed))
p4
ggsave("../output/change_pred_PopGrowth_select_years_over_years_final2.pdf",
       width = 8, height = 7)
######################
############
##Observed versus predicted
############
### computing final estimates at all observed covariate values
second_step_gfr_fit0 = frechet::GloDenReg(xin = Z, qin = second_step_resp0, 
                                          xout = Z, 
                                          optns = list(qSup = qSup, 
                                                       dSup = second_step_resp[[1]]$dSup_mean))$qout
second_step_gfr_fit1 = frechet::GloDenReg(xin = Z, qin = second_step_resp1, 
                                          xout = Z, 
                                          optns = list(qSup = qSup, 
                                                       dSup = second_step_resp[[1]]$dSup_mean))$qout

### Generating Figure S.10
contr_select = c("AUS", "FRA","FIN","JPN","NLD", "USA")
contr_select_ind = which(countries %in% contr_select)
years_select = c(1995, 2000, 2008)
years_select_ind = which(years %in% years_select)
####### loading the observed data -- quantiles and densities
load("../data/obs_dens_quant.Rda")
df_obs = do.call(rbind,lapply(contr_select_ind, function(ind){
  do.call(rbind,lapply(years_select_ind, function(j){
    data.frame(domain = obs[[j]][[ind]]$dens$x,
               dens = obs[[j]][[ind]]$dens$y,
               year_obs = rep(as.character(years[j]),500))
    
  }) )
}))
df_obs$pred = dens_display$pred
q_pert = function(qtrue,alpha){
  eps = sample(c(alpha, -alpha),size = 1, prob = c(.5,.5))
  return(qtrue + eps*alpha*qtrue*(1-qtrue))
}

########### computing the estimated quantiles and densities
outputs =  lapply(1:n, function(ind){
  out = t(sapply(1:ni, function(j){
    q_pert(obs[[j]][[ind]]$quant, .002)
  }))
})
dout_overall = lapply(contr_select_ind, function(ind){
  lapply(1:ni, function(j){
    frechet:::qf2pdf(outputs[[ind]][j,],
                     optns = list(qSup = qSup,
                                  outputGrid = second_step_resp[[1]]$dSup_mean,
                                  nRegGrid =500))
  })
}) 

dens_display = do.call(rbind, lapply(1:length(contr_select_ind), function(ind){
  do.call(rbind,lapply(years_select_ind, function(j){
    data.frame(domain =dout_overall[[ind]][[j]]$x,
               dens =dout_overall[[ind]][[j]]$y,
               year = rep(as.character(years[j]),500))
  }))
})) 
dens_display$pred = c(sapply(contr_select_ind, function(i){
  rep(countries[i], 500*length(years_select))
}))

cols = c("red", "green4","blue")
pred.labs <- c("Australia", "France", "Finland", "Japan", "Netherlands", "United States")
names(pred.labs) <- c("AUS", "FRA","FIN","JPN","NLD", "USA")
p = ggplot() +
  geom_line(data = dens_display,
            aes(x = domain, y = dens, col= year, group = year), size = .2) +
  geom_line(data = df_obs, 
            aes(x = domain, y = dens, col= year_obs, group = year_obs), 
            linetype = 'dashed',
            size = .2) +
  theme_bw() +
  scale_colour_manual(values = cols, name="Calendar year") +
  labs(x  = "Age", 
       y = "Remaining life densities") +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.y = element_blank(),
        legend.position = "bottom") +
  facet_wrap(~pred, labeller = labeller(pred = pred.labs))
p
ggsave("../output/obs_vs_pred_select_years_select_countries_final2.pdf",
       width = 8, height = 7)
############
## Prediction Error
#### Time consuming -- outputs saved as 716
############
library(parallel)
B = 200
ncore = 20
rslt = mclapply(1:B, function(rep){
  samp_ind_tr = sample(1:n, n/2, replace = FALSE)
  qin_tr = lapply(samp_ind_tr, function(ind){
    lapply(1:length(years), function(j){
      obs[[j]][[ind]]$quant
    })
  })
  dat_Z_tr =  Z[samp_ind_tr,]
  dat_Z_tst =  Z[((1:n)[-samp_ind_tr]),]
  second_step_resp = lapply(samp_ind_tr, function(ind){
    subject_level_GFR = subject_level_GFR(ind)
    qin = subject_level_GFR$qin
    dSup_mean = subject_level_GFR$dSup_mean
    obj_fit_ind = subject_level_GFR$obj_fit
    mi_hat0 = obj_fit_ind[1,]
    mi_hat1 = obj_fit_ind[2,]
    return(list(qin = qin, mi_hat0 = mi_hat0, mi_hat1 = mi_hat1,
                dSup_mean = dSup_mean))
  })
  second_step_resp0 = t(sapply(1:(n/2), function(ind){ ###first step estimate of random effects at time 0 -- quantile representation
    second_step_resp[[ind]]$mi_hat0
  }))
  second_step_resp1 = t(sapply(1:(n/2), function(ind){ ###first step estimate of random effects at time 1 -- quantile representation
    second_step_resp[[ind]]$mi_hat1
  }))
  
  qSup = seq(0, 1, length.out = 500)#support of quantile function
  q_geo = function(t,A,B){ ##any point on the geodesic
    t*A + (1-t)*B
  }
  ni = length(years)
  Ti = (years -  years[1])/(years[length(years)] - years[1]) ##standarized time domain
  second_step_gfr_fit0 = frechet::GloDenReg(xin = as.matrix(dat_Z_tr), 
                                            qin = second_step_resp0,
                                            xout = as.matrix(dat_Z_tst),
                                            optns = list(qSup = qSup, ndSup = 500))
  
  second_step_gfr_fit1 = frechet::GloDenReg(xin =  as.matrix(dat_Z_tr),
                                            qin = second_step_resp1,
                                            xout = as.matrix(dat_Z_tst),
                                            optns = list(qSup = qSup, ndSup = 500))
  qout =  lapply(1:nrow(as.matrix(dat_Z_tst)), function(ind){
    out = lapply(1:ni, function(j){
      q_geo(Ti[j], second_step_gfr_fit0$qout[ind,], second_step_gfr_fit1$qout[ind,])
    })
  })
  qin_tst = lapply(((1:n)[-samp_ind_tr]), function(ind){
    lapply(1:length(years), function(j){
      obs[[j]][[ind]]$quant
    })
  })
  rmpe_test = mean(sapply(1:length(qin_tst), function(ind){
    mean(sapply(1:ni, function(j){
      sum((qin_tst[[ind]][[j]] - qout[[ind]][[j]])^2)*diff(qSup)[1]
    }))
  }))
  return(rmpe_test)
}, mc.cores = ncore)
save(rslt, file = "../data/rmpe_mort.Rda")

load("../data/rmpe_mort.Rda")
round(summary(unlist(rslt)),4)
