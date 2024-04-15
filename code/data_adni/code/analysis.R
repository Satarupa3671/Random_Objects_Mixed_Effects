rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(igraph)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(qgraph)
library(corrplot)
library(frechet)

############# additional processing and fitting first step regression
#### Time consuming -- outputs saved as Rda files -- 
##skip running this part upto line 66
load("../data/resp_final_06262022_AAL_CN.Rda") ### data on SPD matrices for CN subjects, preprocessed using AAL percellation from the ADNI fMRI database and converted into correlation matrices
load("../data/covariate_final_06262022_AAL_CN.Rda") ### covariate information including time of scan for the same subjects

n = length(resp2) ### number of CN subjects 
m = dim(resp2[[1]][[1]])[1] ## dimension of the correlation matrix for each subject; AAL percellation gives 90 ROIs

### definition of the end points of the time domain "0" and "1" and aliginng the observations accordingly
dates =  lapply(mydata_new, function(x){
  as.Date(as.character(x$date_acquired)) - min(as.Date(as.character(x$date_acquired)))
})
dates_shifted = lapply(dates, function(dd){
  as.vector(dd)/356
})
T0 = 0
T1 = 9 #max(unlist(dates_shifted))

dates_effect = lapply(dates_shifted, function(dd){
  (as.vector(dd) - T0)/as.vector(T1-T0)
})

mydata_new = mapply(cbind, mydata_new, "dates_effect"= dates_effect, SIMPLIFY=F)
sapply(mydata_new, function(x) x$dates_effect)[[1]]
n = length(resp2)
ind_rep = as.vector((1:n)[sapply(1:length(mydata_new), function(i){
  length(unique(mydata_new[[i]]$dates_effect)) == 1
})])
identical(ind_rep, integer(0))
mydata_new = mydata_new[-ind_rep]
resp2 = resp2[-ind_rep]
n = length(resp2) ### final effective sample size for CN subjects
######
##function to fit the subject level (first step) Global Fr Regression to compute the random effects at the end points of the time domain
subject_level_GFR = function(ind){
  Min = resp2[[ind]]
  Tin = as.matrix(mydata_new[[ind]]$dates_effect)
  obj_fit = frechet::GloCovReg(x = Tin, M = Min, xout = as.matrix(c(0,1)),
                               optns = list(corrOut = TRUE, metric = "power", alpha = .5 ))$Mout
  
  return(list(Min = Min, Tin = Tin, obj_fit = obj_fit))
}

second_step_resp = lapply(1:n, function(ind){
  subject_level_GFR = subject_level_GFR(ind)
  Min = subject_level_GFR$Min
  Tin = subject_level_GFR$Tin
  obj_fit_ind = subject_level_GFR$obj_fit
  mi_hat0 = obj_fit_ind[[1]]
  mi_hat1 = obj_fit_ind[[2]]
  return(list(Min = Min, Tin = Tin, mi_hat0 = mi_hat0, mi_hat1 = mi_hat1))
})
save(second_step_resp, file = "../data/second_step_resp_AAL_CN_Tshifted.Rda")
###################
###################
## Load the necessary files from above
load("../data/resp_final_06262022_AAL_CN.Rda")
load("../data/covariate_final_06262022_AAL_CN.Rda")
load("../data/second_step_resp_AAL_CN_Tshifted.Rda")
n = length(second_step_resp) ### number of CN subjects 

Tin = t(sapply(1:n, function(ind){ ###Time of observation for each subject
  second_step_resp[[ind]]$Tin
}))

Min = t(sapply(1:n, function(ind){ ###input response (correlation matrix) for each subject at each time
  second_step_resp[[ind]]$Min
}))
m = nrow(Min[[1]][[1]]) ## dimension of the correlation matrix for each subject; AAL percellation gives 90 ROIs

second_step_resp0 = lapply(1:n, function(ind){ ###first step estimate of random effects at time 0 -- corr. matrix representation
  second_step_resp[[ind]]$mi_hat0
})
second_step_resp1 = lapply(1:n, function(ind){ ###first step estimate of random effects at time 1 -- corr. matrix representation
  second_step_resp[[ind]]$mi_hat1
})

Z = do.call(rbind,lapply(1:n, function(ind){
  data.frame( Age_first_scan =  mydata_new[[ind]]$Age[1],
              Total_score = mydata_new[[ind]]$Total_score[1])
}))
Z = as.matrix(Z)

############
## Effects of the baseline covariates
############
### Fitting the final estimate at varying values of Total Score (C-Score), while the other predictor (age) is fixed at their mean levels

zout = as.matrix(cbind(rep(mean(Z[,1]), 3), quantile(Z[,2], prob = c(.1,.5,.9))))

second_step_gfr_fit0 = frechet::GloCovReg(x = Z, M = second_step_resp0, xout = zout,
                                          optns = list(corrOut = TRUE,  metric = "power", alpha = .5))$Mout
second_step_gfr_fit1 = frechet::GloCovReg(x = Z, M = second_step_resp1,  xout =zout,
                                          optns = list(corrOut = TRUE,  metric = "power", alpha = .5))$Mout
#### Generating Figure 4
for(i in 1:3){ ## setting diagonals of the estimated corr. matrices to 1
  diag(second_step_gfr_fit0[[i]])= 0
  diag(second_step_gfr_fit1[[i]]) = 0
}
adjPred = list(second_step_gfr_fit0[[1]], second_step_gfr_fit0[[2]], second_step_gfr_fit0[[3]],
               second_step_gfr_fit1[[1]], second_step_gfr_fit1[[2]], second_step_gfr_fit1[[3]])

for(k in 1:6){
  pdf(sprintf("../output/corrplot%d_over_total_score_CN_Tshifted_final.pdf", k))
  corrplot::corrplot(adjPred[[k]],tl.pos = "n")
  dev.off()
}
adjPredGgplot <- lapply(adjPred, reshape2::melt)
corr_CN = list()
corr_CN[[1]] = adjPredGgplot[[6]] - adjPredGgplot[[1]]
corr_CN[[2]] = adjPredGgplot[[5]] - adjPredGgplot[[2]]
corr_CN[[3]] = adjPredGgplot[[4]] - adjPredGgplot[[3]]
for(i in 1:3){
  corr_CN[[i]]$Var1 = adjPredGgplot[[1]]$Var1
  corr_CN[[i]]$Var2 = adjPredGgplot[[1]]$Var2
}
###################
###################
#### Community detection using spectral clustering -- plots for the supplements
for(i in 1:3){
  second_step_gfr_fit0[[i]][(second_step_gfr_fit0[[i]])<.4] = 0
  second_step_gfr_fit1[[i]][(second_step_gfr_fit1[[i]])<.4] = 0
}
# there are six predicted brain networks stored in adjPred, in the form of adjacency matrices
# find the community structure using spectral clustering
adjPred = c(second_step_gfr_fit0,second_step_gfr_fit1)
cl <- list()
for(k in 1:6){
  g <- graph_from_adjacency_matrix(abs(adjPred[[k]]), mode = 'undirected', weighted = TRUE)
  #g = delete.vertices(simplify(g), degree(g)==0)
  #g = simplify(g)
  # cl <- cluster_fast_greedy(g)
  cl[[k]] <- g %>% cluster_leading_eigen()
  mem <- sapply(1:length(unique(cl[[k]]$membership)), function(i) which(cl[[k]]$membership==i))
  mem <- mem[order(sapply(mem, function(x) x[1], simplify=TRUE))]
  cl[[k]]$membership <- sapply(1:90, function(i) which(sapply(mem, function(memi) i%in%memi)))
}
lapply(cl, function(cli) table(cli$membership))
noClusters <- sapply(cl, function(cli) length(table(cli$membership)))
noClusters #c(13, 12, 6,7,6,7) #c(11,4,11,7,3,7)
########## Anatomical regions that contain the ROIs considered
Central_region = c("PRE", "POST", "RO")
Frontal_lobe = c("F1","F2","F3OP","F3T", "F1M","SMA","PCL","F1O","F1MO","F2O","F3O","GR","OC")
Temporal_lobe = c("T1","HES","T2","T3")
Parietal_lobe = c("P1","P2","AG","SMG","PQ")
Occipital_lobe = c("O1","O2","O3","Q","V1","LING","FUSI")
Limbic_lobe = c("T1P","T2P","ACIN","MCIN","PCIN","HIP","PHIP","INS")
Subcortical = c("AMYG","CAU","PUT","PAL","THAL")
node = c(Central_region, Frontal_lobe, Temporal_lobe,Parietal_lobe, Occipital_lobe, Limbic_lobe,Subcortical)
Node = rep(node,2)
#### Generating Figure S.12
col = list(c('#F8766D', 'purple', '#3cb44b', '#f58231', '#4363d8','#ffe119',
             '#46f0f0', '#008080', 'lawngreen', '#f032e6','#e6beff', '#800000',
             'mediumvioletred'),
           c('#F8766D', 'purple', '#3cb44b', '#f58231', '#4363d8','#ffe119',
             '#46f0f0', '#008080', 'lawngreen', '#f032e6','#e6beff', '#800000'),
           c('#F8766D', 'purple', '#3cb44b', '#f58231', '#4363d8','#ffe119'),
           c('#F8766D', 'purple', '#3cb44b', '#f58231', '#4363d8','#ffe119','#46f0f0'),
           c('#F8766D', 'purple', '#3cb44b', '#f58231', '#4363d8','#ffe119'),
           c('#F8766D', 'purple', '#3cb44b', '#f58231', '#4363d8','#ffe119','#46f0f0')
)
markCol <- list()
markCol[[1]] <- rep(NA, length(col[[1]]))
markCol[[1]][c(1, 2, 5, 9, 13)] <- col[[1]][c(1, 2, 5, 9, 13)]
markCol[[2]] <- rep(NA, length(col[[2]]))
markCol[[2]][c(1,2,6)] <- col[[2]][c(1,2,5)]
markCol[[3]] <- rep(NA, length(col[[3]]))
markCol[[3]][c(1,2,4)] <- col[[3]][c(1,2,5)]
markCol[[4]] <- rep(NA, length(col[[4]]))
markCol[[4]][c(1,2,3,4)] <- col[[4]][c(5,1,2,4)]
markCol[[5]] <- rep(NA, length(col[[5]]))
markCol[[5]][c(2,3)] <- col[[5]][c(1,2)]
markCol[[6]] <- rep(NA, length(col[[6]]))
markCol[[6]][c(2,3)] <- col[[6]][c(1,2)]
for(k in 1:6){
  g <- graph_from_adjacency_matrix(adjPred[[k]], mode = 'undirected', weighted = TRUE)
  deg <- igraph::strength(g)
  V(g)$name <- Node
  V(g)$frame.color <- 'white'
  E(g)$width <- 0.8*E(g)$weight/max(E(g)$weight)
  l <- layout_nicely(g) 
  V(g)$color <- col[[k]][cl[[k]]$membership]
  V(g)$community <- cl[[k]]$membership
  V(g)$size <- 6
  E(g)$color <- apply(as.data.frame(get.edgelist(g)), 1, 
                      function(x) ifelse(V(g)$community[which(node==x[1])] == V(g)$community[which(node==x[2])], 
                                         col[[k]][V(g)$community[which(node==x[1])]], '#000000'))
  par(mfrow = c(1,1))
  pdf(sprintf("../output/network_commnities_%d_CN_Tshifted_final.pdf",k))
  plot(cl[[k]], g, layout = l, 
       mark.groups = lapply(1:length(unique(cl[[k]]$membership)), function(i){
         Node[cl[[k]]$membership == i]}),
       edge.color = E(g)$color, vertex.shape = 'none',
       vertex.label.cex = 0.4, vertex.label.font = 2,
       mark.col = adjustcolor(markCol[[k]], 0.05), mark.border = markCol[[k]])
  dev.off()
}
legendCol <- c("#F8766D", "purple", "#4363d8", "lawngreen", "mediumvioletred","#f58231")
# name of anatomical regions of the brain
legendLevel <- c("Frontal", "Temporal", "Occipital","Parietal","Subcortical", "Limbic")
pdf("../output/network_comm_legend_CN_Tshifted_final.pdf")
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend('bottom', legend = legendLevel, cex = .75, pt.cex = 1, col = legendCol, 
       text.col = legendCol, pt.bg = legendCol, pch = 21, 
       bty = 'n', ncol = 6)
dev.off()
#######################
#######################
###Modularity
second_step_gfr_fit0 = frechet::GloCovReg(x = Z, M = second_step_resp0, xout = zout,
                                          optns = list(corrOut = TRUE,  metric = "power", alpha = .5))$Mout
second_step_gfr_fit1 = frechet::GloCovReg(x = Z, M = second_step_resp1,  xout =zout,
                                          optns = list(corrOut = TRUE,  metric = "power", alpha = .5))$Mout

q_geo = function(t,A,B,pow){
  A_svd = eigen(A); B_svd = eigen(B)
  A_pow = A_svd$vectors %*% diag(A_svd$values^pow) %*% t(A_svd$vectors)
  B_pow = B_svd$vectors %*% diag(B_svd$values^pow) %*% t(B_svd$vectors)
  C = (1-t)*A_pow + t*B_pow
  C_svd = eigen(C)
  C_svd$vectors%*%(diag(C_svd$values^(1/pow)))%*%t(C_svd$vectors)
}

modularity_over_geod = lapply(1:3, function(ind){
  set.seed(3346)
  res = sapply(seq(0,1,length.out = 100), function(t){
    out = q_geo(t, second_step_gfr_fit0[[ind]], second_step_gfr_fit1[[ind]],.5)
    diag(out)= 0; out[out<.4] = 0
    g <- graph_from_adjacency_matrix(abs(out), mode = 'undirected', weighted = TRUE)
    cl <- g %>% cluster_leading_eigen()
    mem <- sapply(1:length(unique(cl$membership)), function(i) which(cl$membership==i))
    mem <- mem[order(sapply(mem, function(x) x[1], simplify=TRUE))]
    cl$membership <- sapply(1:m, function(i) which(sapply(mem, function(memi) i%in%memi)))
    return(round(cl$modularity,3))
  })
  return(res)
}) 
df_modularity_CN = data.frame(value = c(modularity_over_geod[[1]], modularity_over_geod[[2]], modularity_over_geod[[3]]))
df_modularity_CN$Time = rep(seq(0,1,length.out = 100),3)
df_modularity_CN$z = c(rep("Z1", 100), rep("Z2", 100), rep("Z3",100))
###
##Global Efficiency over geod
glob_eff_over_geod = lapply(1:3, function(ind){
  set.seed(3346)
  res = sapply(seq(0,1,length.out = 100), function(t){
    out = q_geo(t, second_step_gfr_fit1[[ind]], second_step_gfr_fit0[[ind]],.5)
    diag(out)= 0; out[out<.4] = 0
    g <- graph_from_adjacency_matrix(abs(out), mode = 'undirected', weighted = TRUE)
    return(round(brainGraph::efficiency(g, type = "global"),3))
  })
  return(res)
})
df_glob_eff_CN = data.frame(value = c(glob_eff_over_geod[[1]],
                                      glob_eff_over_geod[[2]],
                                      glob_eff_over_geod[[3]]))
df_glob_eff_CN$Time = rep(seq(0,1,length.out = 100),3)
df_glob_eff_CN$z = c(rep("Z1", 100), rep("Z2", 100), rep("Z3",100))
### Generating Table 1
for(i in 1:3){
  diag(second_step_gfr_fit0[[i]])= 0
  diag(second_step_gfr_fit1[[i]]) = 0
}
for(i in 1:3){
  second_step_gfr_fit0[[i]][(second_step_gfr_fit0[[i]])<.4] = 0
  second_step_gfr_fit1[[i]][(second_step_gfr_fit1[[i]])<.4] = 0
}
adjPred = c(second_step_gfr_fit0,second_step_gfr_fit1)


for(k in c(1,2,3,6,5,4)){
  g <- graph_from_adjacency_matrix(abs(adjPred[[k]]), mode = 'undirected', weighted = TRUE)
  print(paste0("Modularity = ", round(igraph::modularity(g, cl[[k]]$membership),3)))
}
#####################################################
#####################################################
### Similar analysis for MCI
########
load("../data/resp_final_06262022_AAL_MCI.Rda")
load("../data/covariate_final_06262022_AAL_MCI.Rda")
load("../data/second_step_resp_AAL_MCI_Tshifted.Rda")
n = length(second_step_resp)
Tin = t(sapply(1:n, function(ind){
  second_step_resp[[ind]]$Tin
}))
Min = t(sapply(1:n, function(ind){
  second_step_resp[[ind]]$Min
}))
second_step_resp0 = lapply(1:n, function(ind){
  second_step_resp[[ind]]$mi_hat0
})
second_step_resp1 = lapply(1:n, function(ind){
  second_step_resp[[ind]]$mi_hat1
})
Z = do.call(rbind,lapply(1:n, function(ind){
  data.frame( Age_first_scan =  mydata_new[[ind]]$Age[1],
              Total_score = mydata_new[[ind]]$Total_score[1])
}))
Z = as.matrix(Z)
zout = as.matrix(cbind(rep(median(Z[,1]), 3), quantile(Z[,2], prob = c(.1,.5,.9))))
#####
second_step_gfr_fit0 = frechet::GloCovReg(x = Z, M = second_step_resp0, xout = zout,
                                          optns = list(corrOut = TRUE,  metric = "power", alpha = .5))$Mout
second_step_gfr_fit1 = frechet::GloCovReg(x = Z, M = second_step_resp1,  xout =zout,
                                          optns = list(corrOut = TRUE,  metric = "power", alpha = .5))$Mout
####
for(i in 1:3){
  diag(second_step_gfr_fit0[[i]])= 0
  diag(second_step_gfr_fit1[[i]]) = 0
}
adjPred = c(second_step_gfr_fit0,second_step_gfr_fit1)
adjPredGgplot <- lapply(adjPred, reshape2::melt)
corr_MCI = list()
corr_MCI[[1]] = adjPredGgplot[[4]] - adjPredGgplot[[1]]
corr_MCI[[2]] = adjPredGgplot[[5]] - adjPredGgplot[[2]]
corr_MCI[[3]] = adjPredGgplot[[6]] - adjPredGgplot[[3]]
for(i in 1:3){
  corr_MCI[[i]]$Var1 = adjPredGgplot[[1]]$Var1
  corr_MCI[[i]]$Var2 = adjPredGgplot[[1]]$Var2
}
###
for(i in 1:3){
  second_step_gfr_fit0[[i]][(second_step_gfr_fit0[[i]])<.4] = 0
  second_step_gfr_fit1[[i]][(second_step_gfr_fit1[[i]])<.4] = 0
}
adjPred = c(second_step_gfr_fit0,second_step_gfr_fit1)

cl <- list()
for(k in 1:6){
  g <- graph_from_adjacency_matrix(abs(adjPred[[k]]), mode = 'undirected', weighted = TRUE)
  cl[[k]] <- g %>% cluster_leading_eigen()
  mem <- sapply(1:length(unique(cl[[k]]$membership)), function(i) which(cl[[k]]$membership==i))
  mem <- mem[order(sapply(mem, function(x) x[1], simplify=TRUE))]
  cl[[k]]$membership <- sapply(1:m, function(i) which(sapply(mem, function(memi) i%in%memi)))
}
lapply(cl, function(cli) table(cli$membership))
noClusters <- sapply(cl, function(cli) length(table(cli$membership)))
noClusters  #c(11,4,11,7,3,7)
###Modularities
second_step_gfr_fit0 = frechet::GloCovReg(x = Z, M = second_step_resp0, xout = zout,
                                          optns = list(corrOut = TRUE,  metric = "power", alpha = .5))$Mout
second_step_gfr_fit1 = frechet::GloCovReg(x = Z, M = second_step_resp1,  xout =zout,
                                          optns = list(corrOut = TRUE,  metric = "power", alpha = .5))$Mout

modularity_over_geod = lapply(1:3, function(ind){ 
  set.seed(3346)
  res = sapply(seq(0,1,length.out = 100), function(t){
    out = q_geo(t, second_step_gfr_fit0[[ind]], second_step_gfr_fit1[[ind]],.5)
    diag(out)= 0; out[out<.4] = 0
    g <- graph_from_adjacency_matrix(abs(out), mode = 'undirected', weighted = TRUE)
    #g = simplify(g)
    cl <- g %>% cluster_leading_eigen()
    mem <- sapply(1:length(unique(cl$membership)), function(i) which(cl$membership==i))
    mem <- mem[order(sapply(mem, function(x) x[1], simplify=TRUE))]
    cl$membership <- sapply(1:m, function(i) which(sapply(mem, function(memi) i%in%memi)))
    return(round(cl$modularity,3))
  })
  return(res)
})  
df_modularity_MCI = data.frame(value = c(modularity_over_geod[[1]], modularity_over_geod[[2]], modularity_over_geod[[3]]))
df_modularity_MCI$Time = rep(seq(0,1,length.out = 100),3)
df_modularity_MCI$z = c(rep("Z1", 100), rep("Z2", 100), rep("Z3",100))
###
##Global Efficiency over geod
glob_eff_over_geod = lapply(1:3, function(ind){
  set.seed(3346)
  res = sapply(seq(0,1,length.out = 100), function(t){
    out = q_geo(t, second_step_gfr_fit1[[ind]], second_step_gfr_fit0[[ind]],.5)
    diag(out)= 0; out[out<.4] = 0
    g <- graph_from_adjacency_matrix(abs(out), mode = 'undirected', weighted = TRUE)
    return(round(brainGraph::efficiency(g, type = "global"),3))
  })
  return(res)
})
df_glob_eff_MCI = data.frame(value = c(glob_eff_over_geod[[1]],
                                       glob_eff_over_geod[[2]],
                                       glob_eff_over_geod[[3]]))
df_glob_eff_MCI$Time = rep(seq(0,1,length.out = 100),3)
df_glob_eff_MCI$z = c(rep("Z1", 100), rep("Z2", 100), rep("Z3",100))


###Modularities ### Table 1
for(k in c(1,2,3,4,5,6)){
  g <- graph_from_adjacency_matrix(abs(adjPred[[k]]), mode = 'undirected', weighted = TRUE)
  print(paste0("Modularity = ", round(igraph::modularity(g, cl[[k]]$membership),3)))
}
#####
### correlation matrix difference plot -- Generating Figure 5
corr_CN2 = do.call(rbind, corr_CN)
corr_CN2$z = c(rep("Q10", m*m), rep("Q50", m*m), rep("Q90", m*m))

corr_MCI2 = do.call(rbind, corr_MCI)
corr_MCI2$z = c(rep("Q10", m*m), rep("Q50", m*m), rep("Q90", m*m))

corr_combined = rbind(corr_CN2, corr_MCI2)
corr_combined$group = c(rep("CN", m*m*3), rep("MCI", m*m*3))


p = ggplot(corr_combined, aes(x = Var2, y = Var1, group = z)) + 
  geom_raster(aes(fill=value, group = z)) + 
  scale_fill_gradient2(low="blue3",mid = "grey90", high="red1", name = "Correlation value",
                       guide = guide_colourbar(barwidth = 10, barheight = .8)) +
  theme_bw() +
  labs(x = " ", y = " ") +
  theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
        axis.text.y=element_text(size=9),
        axis.title = element_text(size  = 12),
        legend.title = element_text(size  = 12),
        strip.text.x = element_text(size  = 12),
        strip.text.y = element_text(size  = 12),
        legend.position = "bottom") +
  facet_grid(z~group) 
p
ggsave("../output/diff_corrplot_over_total_score_CN_MCI_final.pdf", width = 8, height = 7)
###### modularity plot --  Generating Figure S.14
df_modularity = rbind(df_modularity_CN, df_modularity_MCI)
df_modularity$group = c(rep("CN", nrow(df_modularity_CN)), rep("MCI", nrow(df_modularity_MCI)))
p1 = ggplot(df_modularity)+
  geom_line(aes(x = Time, y = value , group = z, color = z), size = 1) +
  labs(#title = "Euclidean response",
    x = "Time since first scan", y = "Modularity") +
  scale_colour_manual(values = c("#F8766D","#619CFF","#C77CFF"),
                      name = " ", labels = c("Fitted at 10% Quantile\n of Total Score",
                                             "Fitted at 50% Quantile\n of Total Score",
                                             "Fitted at 90% Quantile\n of Total Score")) + 
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size  = 14),
        legend.title = element_text(size  = 12),
        strip.text.y = element_blank(),
        legend.position="bottom") +
  facet_wrap(~group)
p1
ggsave("../output/modularity_over_total_score_CN_MCI_final.pdf",
       width = 8, height = 7)
##########
###### global efficiency plot --  Generating Figure S.15
df_glob_eff = rbind(df_glob_eff_CN, df_glob_eff_MCI)
df_glob_eff$group = c(rep("CN", nrow(df_glob_eff_CN)), rep("MCI", nrow(df_glob_eff_MCI)))
p2 = ggplot(df_glob_eff)+
  geom_line(aes(x = Time, y = value , group = z, color = z), size = 1) +
  labs(#title = "Euclidean response",
    x = "Time since first scan", y = "Global Efficiency") +
  scale_colour_manual(values = c("#F8766D","#619CFF","#C77CFF"),
                      name = " ", labels = c("Fitted at 10% Quantile\n of Total Score",
                                             "Fitted at 50% Quantile\n of Total Score",
                                             "Fitted at 90% Quantile\n of Total Score")) + 
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size  = 14),
        legend.title = element_text(size  = 12),
        strip.text.y = element_blank(),
        legend.position="bottom") +
  facet_wrap(~group)
p2
ggsave("../output/Glob_eff_over_total_score_CN_MCI_final.pdf",
       width = 8, height = 7)


