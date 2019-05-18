
library(poweRlaw)
library(SpiecEasi)
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
co.var <- function(x,na.rm=TRUE) 100*(sd(x,na.rm=na.rm)/mean(x,na.rm=na.rm))

#kjhealy/polar-labels.r on github: https://gist.github.com/kjhealy/834774#file-polar-labels-r
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(as.data.frame(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)),stringsAsFactors=FALSE))
}

#quartile coefficient of dispersion
qcd <- function(x){(quantile(x,0.75)-quantile(x,0.25))/(quantile(x,0.75)+quantile(x,0.25))}

#############
# Load data #
#############

setwd("~/FunctionalNetworks_Project")
pfam_profile <- read.delim("yams_pd1_pub/yams_pd1_pub_pfam_ppm.txt",
                             stringsAsFactors = FALSE)
rownames(pfam_profile) <- pfam_profile$Feature
pfam_profile$Feature <- NULL
pp_tr <- t(pfam_profile) %>% as.data.frame(stringsAsFactors=F)
samples <- rownames(pp_tr)

patient_meta <- read.delim("yams_pd1_pub/yams_pd1_pub_meta.txt",
                           stringsAsFactors = FALSE)
patient_meta$ResponseBinary <- as.numeric(factor((patient_meta$Clin_Response)))-1
patient_meta <- patient_meta %>% subset(select=c(Sample,ResponseBinary))
patient_meta$Sample <- gsub("-",".",patient_meta$Sample)

#Remove frunctions with less than 10 observations
sum(colSums(pp_tr!=0)<10)
pp_tr <- pp_tr[,colSums(pp_tr!=0)>10]

plot(apply(pp_tr,2,mean,na.rm=T),
     apply(pp_tr,2,var),log="xy",
     xlab="Mean",ylab="Variance")

plot(apply(pp_tr,2,mean,na.rm=T),
     apply(pp_tr,2,co.var),log="xy",
     xlab="Mean",ylab="Coefficient of Variation")


hist(apply(pp_tr,2,mean,na.rm=T))

summary(lm(log(apply(pp_tr,2,var))~log(apply(pp_tr,2,mean,na.rm=T))))
summary(lm(log(apply(pp_tr,2,co.var))~log(apply(pp_tr,2,mean,na.rm=T))))

hist(apply(pp_tr,2,co.var),breaks=50)
plot(density(apply(pp_tr,2,co.var)))
abline(v=200,col="red")

#Remove frunctions with low coefficient of variation ()
sum(apply(pp_tr,2,co.var)>200)
pp_tr <- pp_tr[,apply(pp_tr,2,co.var)>200]

rownames(pp_tr) <- samples
pp_tr_response <- pp_tr %>% 
  subset(rownames(pp_tr) %in% patient_meta$Sample[patient_meta$ResponseBinary==1])
pp_tr_noresponse <- pp_tr %>% 
  subset(rownames(pp_tr) %in% patient_meta$Sample[patient_meta$ResponseBinary==0])



#Profile metadata
setwd("~/FunctionalNetworks_Project")
pfam_feat<- read.delim("yams_pd1_pub/yams_pd1_pub_pfam_feat.txt",
                         stringsAsFactors = FALSE)

#################
# Build Network #
#################

ppr.net <- spiec.easi(as.matrix(pp_tr_response), method='mb', lambda.min.ratio=1e-2,
                     nlambda=15, pulsar.select = T,pulsar.params=list(rep.num=50),verbose=T)
ppr.net.ig <- adj2igraph(symBeta(getOptBeta(ppr.net), mode='maxabs'))
ppr.net.ig$names <- colnames(pp_tr_response)

save(ppr.net,ppr.net.ig,
     file="pfam_net_CV200_byresponse_responder.RData")

ppnr.net <- spiec.easi(as.matrix(pp_tr_noresponse), method='mb', lambda.min.ratio=1e-2,
                      nlambda=15, pulsar.select = T,pulsar.params=list(rep.num=50),verbose=T)
ppnr.net.ig <- adj2igraph(symBeta(getOptBeta(ppnr.net), mode='maxabs'))
ppnr.net.ig$names <- colnames(pp_tr_response)

save(ppnr.net,ppnr.net.ig,
     file="pfam_net_CV200_byresponse_nonresponder.RData")



#############
# Visualize #
#############

###Entire Network
#agcoord <- layout_nicely(pp.net.ig)
#agcoord <- layout_with_fr(pp.net.ig)
agcoord <- layout.circle(pp.net.ig)
lab.locs <- radian.rescale(x=1:ncol(pp_tr), direction=-1, start=0)
plotnet <- function(ig, coord=agcoord, ...){
  plot(ig, layout=coord, vertex.size=2, vertex.label=ig$names,
       vertex.label.dist=1, vertex.label.degree=lab.locs,
       vertex.label.cex=0.5,
       edge.color=ifelse(E(ig)$weight > 0, 'green', 'red'), ...)
}
plotnet(pp.net.ig, main="Function Interaction Network")


###Positive Interactions
pp.net.ig.p <- delete.edges(pp.net.ig, which(E(pp.net.ig)$weight < 0))
agcoord <- layout_with_fr(pp.net.ig.p)
lab.locs <- radian.rescale(x=1:ncol(pp_tr), direction=-1, start=0)
plotnet <- function(ig, coord=agcoord, ...){
  plot(ig, layout=coord, vertex.size=2, vertex.label=pp.net.ig$names,
       vertex.label.dist=1, vertex.label.degree=lab.locs,
       vertex.label.cex=0.25,
       edge.color=ifelse(E(ig)$weight > 0, 'green', 'red'), ...)
}
plotnet(pp.net.ig.p, main="Functions: Positive Interactions")

###Negative Interactions
pp.net.ig.n <- delete.edges(pp.net.ig, which(E(pp.net.ig)$weight > 0))
lab.locs <- radian.rescale(x=1:ncol(pp_tr), direction=-1, start=0)
plotnet <- function(ig, coord=agcoord, ...){
  plot(ig, layout=coord, vertex.size=2, vertex.label=pp.net.ig$names,
       vertex.label.dist=1, vertex.label.degree=lab.locs,
       vertex.label.cex=0.25,
       edge.color=ifelse(E(ig)$weight > 0, 'green', 'red'), ...)
}
plotnet(pp.net.ig.n, main="Functions: Negative Interactions")


setwd("~/FunctionalNetworks_Project")
pdf("FunctionalNetwork_pfam_spieceasimb.pdf",width=14,height=7)
par(mfrow=c(1,2))
plotnet(pp.net.ig.p, main="Functions: Positive Interactions")
plotnet(pp.net.ig.n, main="Functions: Negative Interactions")
par(mfrow=c(1,1))
dev.off()


######################
# Network Statistics #
######################

### Degree Distribution ##
hist(degree(pp.net.ig),breaks=100,main="Positive",xlab="Degree")


dd <- degree_distribution(pp.net.ig)
plot(0:(length(dd)-1),dd,log="xy")

# #Test for power law
# degree(pp.net.ig)
# m_sp = displ$new(degree(pp.net.ig)[degree(pp.net.ig)>0])
# est_sp = estimate_xmin(m_sp)
# m_sp$setXmin(est_sp)
# m_sp$xmin 
# m_sp$pars
# 
# # #Bootstrap parameter estimates
# # bs = bootstrap(m_sp, no_of_sims=5000, threads=2)
# # plot(bs, trim=0.1)
# # hist(bs$bootstraps[,2])
# # hist(bs$bootstraps[,3])
# 
# #Bootstrap to test if power law
# bs_p <- bootstrap_p(m_sp,no_of_sims = 1000)
# bs_p$p 
# plot(bs_p)
# 
# #Another test using the "Vuong" procedure (compare against lognormal)
# m_ln = dislnorm$new(degree(pp.net.ig)[degree(pp.net.ig)>0])
# est = estimate_xmin(m_ln)
# m_ln$setXmin(m_sp$getXmin())
# est = estimate_pars(m_ln)
# m_ln$setPars(est)
# comp = compare_distributions(m_sp, m_ln)
# comp

#With igraph
fit_power_law(degree(pp.net.ig)[degree(pp.net.ig)>0],xmin=3) #NOT apower law (but we removed the low abundance nodes...)

### Distributions of Shortest Paths ##

n <- length(V(pp.net.ig))

#Breadth-First Search
bfs.vec <- matrix(nrow=n,ncol=n)
for(i in 1:n){
  bfs.vec[i,] <- bfs(pp.net.ig, root=i, dist=TRUE,unreachable = F)$dist
  
}

#Plot distance distributions
distd <- table(bfs.vec)/sum(bfs.vec,na.rm=T)
plot(as.numeric(names(distd)),as.numeric(distd),log="xy",xlab="Distance",
     ylab=expression(P[Distance]))


### Clustering ##

cl <- transitivity(pp.net.ig,type="local")
k <- degree(pp.net.ig)
cl.df <- data.frame(clust=cl,degree=k) %>% group_by(degree) %>% 
  summarize(mclust=mean(clust,na.rm=T))
plot(mclust~degree,data=cl.df,log="xy",xlab="k",ylab="C(k)")

####################################
# Community Detection in + Network #
####################################

# 
# #Optimal CLustering Based on Modularity
# pp.net.ig.p.clust <- cluster_optimal(pp.net.ig.p)
# modularity(pp.net.ig.p.clust)
# plot(pp.net.ig.p.clust,pp.net.ig.p,layout=agcoord)
# 
# #Greedy Clustering
# pp.net.ig.p.clustgreedy <- cluster_fast_greedy(pp.net.ig.p)
# modularity(pp.net.ig.p.clustgreedy)
# plot(pp.net.ig.p.clustgreedy,pp.net.ig.p,layout=agcoord)
# 
# #Spectral Clustering
# pp.net.ig.p.clustspec <- cluster_leading_eigen(pp.net.ig.p)
# plot(pp.net.ig.p.clustspec,pp.net.ig.p,layout=agcoord)
# # 
# # #Betweenness Clustering
# # pp.net.ig.p.clustbt <- cluster_edge_betweenness(pp.net.ig.p)
# # modularity(pp.net.ig.p.clustbt)
# # setwd("~/hmptraits/Figs")
# # pdf("TraitNetwork_BetweenessClustering.pdf",width=10,height=10, onefile=FALSE)
# # plot(pp.net.ig.p.clustbt,pp.net.ig.p,layout=agcoord)
# # dev.off()

#InfoMap Clustering
pp.net.ig.p.clustinfo <- infomap.community(pp.net.ig.p)

agcoord <- layout_nicely(pp.net.ig.p)
#agcoord <- layout_with_fr(pp.net.ig)
#agcoord <- layout.circle(pp.net.ig)
pdf("FunctionNetworkPfam_InfoMapClustering.pdf",width=10,height=10, onefile=FALSE)
plot(pp.net.ig.p.clustinfo,pp.net.ig.p,layout=agcoord, 
     vertex.label=pp.net.ig.p$names,vertex.label.cex=0.25,
     vertex.size=4)
dev.off()

table(pp.net.ig.p.clustinfo$membership)

pdf("FunctionNetworkPfam_InfoMapClustering_topgroups.pdf",width=10,height=10, onefile=FALSE)
table(pp.net.ig.p.clustinfo$membership)
pp.net.ig.p.clustinfo.cut <- pp.net.ig.p.clustinfo
pp.net.ig.p.clustinfo.cut$membership[pp.net.ig.p.clustinfo.cut$membership>20] <- 0
plot(pp.net.ig.p,layout=agcoord, 
     vertex.label=pp.net.ig.p$names,vertex.label.cex=0.1,
     vertex.size=4,vertex.color=pp.net.ig.p.clustinfo.cut$membership)
dev.off()

table(pp.net.ig.p.clustinfo$membership)
clust_df <- data.frame(Cluster=1,
                       Accession=pp.net.ig.p$names[pp.net.ig.p.clustinfo$membership==1],
                       stringsAsFactors = FALSE)
for(i in 2:max(pp.net.ig.p.clustinfo$membership)){
  m_i <- data.frame(Cluster=i,
                    Accession=pp.net.ig.p$names[pp.net.ig.p.clustinfo$membership==i],
                    stringsAsFactors = FALSE)
  clust_df <- rbind(clust_df,m_i)
}

clust_df <- merge.easy(clust_df,pfam_feat,key="Accession")

write.csv(clust_df,file="neighborhood_detection_pfam_infomap.csv")

cl_181 <- clust_df %>% subset(Cluster==150)
cl_9 <- clust_df %>% subset(Cluster==77)
cl_24 <- clust_df %>% subset(Cluster==78)
cl_19 <- clust_df %>% subset(Cluster==106)
