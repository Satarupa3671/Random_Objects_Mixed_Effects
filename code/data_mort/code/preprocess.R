#rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)
##remaining life data
datapath = "../data/fltper_1x1"
filenames <- list.files(datapath, full.names=TRUE)
# loading all the data
dataAll = list()
for (i in filenames){
  country = strsplit(i, "/")
  country = substr(country[[1]][length(country[[1]])],1,3)
  dataAll[[country]] = read.table(i, skip = 1, header = T)
  # dataAll[[country]]$Year = as.numeric(as.character(dataAll[[country]]$Year))
  # dataAll[[country]]$Age = as.numeric(as.character(dataAll[[country]]$Age))
  # dataAll[[country]]$lx = as.numeric(as.character(dataAll[[country]]$lx))
}




# finding common year range for each country
yearRange = data.frame(matrix(nrow=0, ncol = 3))
colnames(yearRange) = c("country", "year_start", "year_end")
for (i in names(dataAll)){
  r = range(dataAll[[i]]$Year)
  row = data.frame("country" = i, "year_start" =  r[1], "year_end" =  r[2] )
  print(row)
  yearRange = rbind(yearRange,row)
  #print(yearRange)
}

library(dplyr)
library(tibble)
yearEnd = 2019
yearStart = 1990
countries1990_2019 = yearRange %>% filter(year_start <= yearStart) %>% filter(year_end >= yearEnd) %>%select(country)

column_names = names(dataAll[['USA']] %>% add_column("Country" = 'USA'))
df = setNames(data.frame(matrix(nrow = 0, ncol=length(column_names))), column_names)

for (i in countries1990_2019$country){
  df = rbind(df, dataAll[[i]] %>% filter(Year>=1990) %>%  filter(Year<=2019) %>%
                             add_column("Country" = i))
 
}
df$Age[df$Age == "110+"] = "110"


df$Year <- as.numeric(as.character(df$Year))
df$Age <- as.numeric(as.character(df$Age))
df$lx <- as.numeric(as.character(df$lx))
df = df[df$Age >=75,]  
countries = unique(df$Country)
yearStart = 1990
yearEnd = 2019
years = yearStart:yearEnd
write.csv(df, "../data/lt_subset2.csv", row.names = F)

for(i in unique(df$Country)){
  print(unique(df$Year[df$Country == i]))
}
#View(df)
df$Country[df$Country == 'TWN'] = rep("CHN", length(df$Country[df$Country == 'TWN']))
#write.csv(df, "lt_subset.csv", row.names =F)
###
dat_unemployment = read.csv("../data/Unemployment/API_SL.UEM.TOTL.ZS_DS2_en_csv_v2_3401869.csv", header = TRUE)
dat2_unemployment = do.call(rbind, lapply(unique(df$Country), function(cnt){
  dat_unemployment[dat_unemployment$Country.Code == cnt,]
}))
dat2_unemployment$X1991
###
dat_fertility = read.csv("../data/Fertility/API_SP.DYN.TFRT.IN_DS2_EN_csv_v2_3404027.csv", header = TRUE)
dat2_fertility = do.call(rbind, lapply(unique(df$Country), function(cnt){
  dat_fertility[dat_fertility$Country.Code == cnt,]
}))
dat2_fertility$X1991
###
dat_gdp = read.csv("../data/GDP/API_NY.GDP.PCAP.PP.CD_DS2_en_csv_v2_3401652.csv", header = TRUE)
dat2_gdp = do.call(rbind, lapply(unique(df$Country), function(cnt){
  dat_gdp[dat_gdp$Country.Code == cnt,]
}))
dat2_gdp$X1997
##
dat_pop_growth = read.csv("../data/Pop_Growth/API_SP.POP.GROW_DS2_en_csv_v2_3404396.csv", header = TRUE)
dat2_pop_growth = do.call(rbind, lapply(unique(df$Country), function(cnt){
  dat_pop_growth[dat_pop_growth$Country.Code == cnt,]
}))
dat2_pop_growth$X1997
###

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
countries[countries == "TWN"] = "CHN"
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
  qin = lapply(1:length(years), function(year){
    dens_quant_country_ind[[year]]$quant
  })
  
  return(list(qin = qin, obj_fit = obj_fit, dSup_mean = dSup_mean))
}
n = length(countries)
second_step_resp = lapply(1:n, function(ind){
  subject_level_GFR = subject_level_GFR(ind)
  qin = subject_level_GFR$qin
  dSup_mean = subject_level_GFR$dSup_mean
  obj_fit_ind = subject_level_GFR$obj_fit
  mi_hat0 = obj_fit_ind[1,]
  mi_hat1 = obj_fit_ind[2,]
  return(list(qin = qin, mi_hat0 = mi_hat0, mi_hat1 = mi_hat1, dSup_mean = dSup_mean))
})
save(second_step_resp, file = "second_step_resp2.Rda")
second_step_resp0 = t(sapply(1:n, function(ind){
  second_step_resp[[ind]]$mi_hat0
}))
second_step_resp1 = t(sapply(1:n, function(ind){
  second_step_resp[[ind]]$mi_hat1
}))
dSup_mean = colMeans(t(sapply(1:n, function(ind){
  second_step_resp[[ind]]$dSup_mean
})))

###
Z = data.frame(unemployment = dat2_unemployment$X1991,
               fertility = dat2_fertility$X1991, 
               gdp_per_capita = dat2_gdp$X1997,
               pop_growth = dat2_pop_growth$X1997)
write.csv(Z, file = '../data/covariate_1990.csv', row.names = FALSE)
Z = cbind(dat2_unemployment$X1991,dat2_fertility$X1991, dat2_gdp$X1997,dat2_pop_growth$X1997)
n = length(countries)
##
qSup = seq(0, 1, length.out = 500) #support of quantile func
zout = Z #seq(rangeZ[1],rangeZ[2], length.out = 5)

second_step_gfr_fit0 = frechet::GloDenReg(xin = Z, qin = second_step_resp0, 
                                 xout = zout, optns = list(qSup = qSup, dSup = dSup_mean))$qout

second_step_gfr_fit1 = frechet::GloDenReg(xin = Z, qin = second_step_resp1, 
                                 xout = zout, optns = list(qSup = qSup, dSup = dSup_mean))$qout

second_step_gfr_QQ = sapply(1:nrow(zout), function(k){
  (second_step_gfr_fit0[k,] + second_step_gfr_fit1[k,])/2
})
second_step_gfr_dens = lapply(1:nrow(zout), function(k){
  frechet:::qf2pdf(second_step_gfr_QQ[,k])
})

####plots of observed vs fitted
obs = lapply(years, function(year){
  lapply(countries, function(country){
    data_density(country,year) 
  })
})
save(obs, file = "../data/obs_dens_quant.Rda")

