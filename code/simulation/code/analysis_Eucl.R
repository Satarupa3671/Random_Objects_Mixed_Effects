rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(lme4)
library(HLMdiag)
library(haven)

fat<-read_dta('../data/fat.dta') ### Read dataset
attach(fat)
fat$time2 = (fat$time - min(fat$time))/(max(fat$time) - min(fat$time)) ### aligning the time domain

### classical linear mixed effects model fitting using lme4 library
mod = lmer(formula = pbf ~ age + agemen + (time2|id)-1, data = fat)
summary(mod)
############
### computing prediction error using lme4
subj = unique(as.numeric(row.names(ranef(mod)$id)))
n = length(subj)
cv_err = sapply(subj, function(i){
  fat_cv = fat[fat$id != i,]
  mod1 = lmer(formula = pbf ~ age + agemen + (time2 | id)-1, data=fat_cv)
  pred = predict(mod1, newdata = fat[fat$id == i,], allow.new.levels = TRUE)
  err = mean(abs(fat[fat$id == i,]$pbf - pred))/n
})
round(mean(cv_err),4) ###0.0362
round(sd(cv_err),4) ###0.0212
##############
### computing prediction error using two-step approach
subj = unique(as.numeric(row.names(ranef(mod)$id)))
cv_err2 = sapply(subj, function(i){
  fat_cv = fat[fat$id != i,]
  subj_wo_i = unique(fat_cv$id)
  coeff = do.call(c,lapply(subj_wo_i, function(j){
    fat_subj = fat_cv[fat_cv$id == j,]
    mod2 = lm(pbf ~ time2, data = fat_subj)
    predict(mod2)
  }))
  final = lm(coeff ~ age + agemen -1, data = fat_cv)
  #coeff = lmer(formula = pbf ~  (0+time2 | id), data = fat_cv)
  #final = lm(predict(coeff) ~ age + agemen-1, data = fat_cv)
  pred = predict(final, newdata = fat[fat$id == i,], allow.new.levels = TRUE)
  err = mean(abs(fat[fat$id == i,]$pbf - pred))/n
})
round(mean(cv_err2),4) ### 0.0365
round(sd(cv_err2),4) ### 0.0213

########################
coeff = do.call(c,lapply(subj, function(j){
  fat_subj = fat[fat$id == j,]
  mod2 = lm(pbf ~ time2, data = fat_subj)
  predict(mod2)
}))
final = lm(coeff ~ age + agemen -1, data = fat)
#coeff = lmer(formula = pbf ~  (time2 | id)-1, data = fat)
#final = lm(predict(coeff) ~ age + agemen -1, data = fat)
round(final$coefficients,4)
round(confint(mod),4)


### lme estimate
est = predict(mod)

### our estimate
B = cbind(fat$age,fat$agemen)
est2 = B%*%final$coefficients
#vhat=ranef(coeff)
#re<-data.matrix(vhat$id)

###########################################################
########## Visualization of fits --  Generating Figure S.16.
df_est= data.frame(id = fat$id,
                lme4 = est,
                obs = fat$pbf,
                estimated = est2,
                re = coeff)

colors <- c(#"lme4" = "#7CAE00", 
  "observed" = "#F8766D", 
            "estimated" = "#00BFC4", "random effects" ="purple4" )

library(ggplot2)
p1 = ggplot() + 
  #geom_point(aes(x = id, y = lme4, color = "lme4") , data = df_est) +
  geom_point(aes(x = id, y = obs, color = 'observed'), data = df_est) +
  geom_point(aes(x = id, y = estimated, color = 'estimated'), data = df_est) +
  geom_point(aes(x = id, y = re, color = 'random effects'),data = df_est) +
           #  data = data.frame(id = unique(fat$id), 
            #                   re = re[,2])) +
  scale_color_manual(labels = c("Random Effects (slope)",
                                "Observed",
                                "Estimated"),
                                values = colors) + 
  theme_bw() + ylim(-20,70) +
  labs(x  = "Subject ID", 
       y = "Percentage Body Fat", color = "") +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.y = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "bottom") 
p1
ggsave("../output/compare_Eucl_estd.pdf", width = 8, height = 7)
# plot(fat$id, est,type="p",axes=F,ylim=c(-50,100),pch=20,col="blue",xlab="",ylab="")
# par(new=T)
# plot(fat$id, est2,type="p",pch=20,axes=F,ylim=c(-50,100),col="darkgreen",xlab="",ylab="")
# par(new=T)
# plot(unique(fat$id), re[,2],type="p",pch=20,ylim=c(-50,100),col="red",xlab="Subject ID",ylab="Mean",main="Estimated % Body Fat")
# par(new=F)
# legend('topleft',bty='n',c("lme4","Proposed","Random Effect"),cex=1,pch=c(20,20,20),col=c("blue","green","red"))
# 
