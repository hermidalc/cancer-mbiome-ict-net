## JLW - 11/19/18
# Builds network from Pfam annotations w/ patient outcome as an additional node

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

## Pfam Abundance Profiles
setwd("~/FunctionalNetworks_Project")
pfam_profile <- read.delim("yams_pd1_pub/yams_pd1_pub_pfam_ppm.txt",
                             stringsAsFactors = FALSE)
rownames(pfam_profile) <- pfam_profile$Feature
pfam_profile$Feature <- NULL
pp_tr <- t(pfam_profile) %>% as.data.frame(stringsAsFactors=F)
samples <- rownames(pp_tr)

## Patient Metadata (Response)
patient_meta <- read.delim("yams_pd1_pub/yams_pd1_pub_meta.txt",
                           stringsAsFactors = FALSE)
patient_meta$ResponseBinary <- as.numeric(factor((patient_meta$Clin_Response)))-1
patient_meta <- patient_meta %>% subset(select=c(Sample,ResponseBinary))
patient_meta$Sample <- gsub("-",".",patient_meta$Sample)

#Remove Pfams with less than 10 observations
sum(colSums(pp_tr!=0)<10)
pp_tr <- pp_tr[,colSums(pp_tr!=0)>10]

# We want to remove low variability/high abundance functions (not intereting if universal - e.g. ribosomal)
plot(apply(pp_tr,2,mean,na.rm=T),
     apply(pp_tr,2,var),log="xy",
     xlab="Mean",ylab="Variance")
summary(lm(log(apply(pp_tr,2,var))~log(apply(pp_tr,2,mean,na.rm=T))))
plot(apply(pp_tr,2,mean,na.rm=T),
     apply(pp_tr,2,co.var),log="xy",
     xlab="Mean",ylab="Coefficient of Variation")
summary(lm(log(apply(pp_tr,2,co.var))~log(apply(pp_tr,2,mean,na.rm=T))))

# Use coefficient of variation to do this - It even looks bimodal, natural cutoff (coeff var > 200 only)
hist(apply(pp_tr,2,co.var),breaks=50)
plot(density(apply(pp_tr,2,co.var)))
abline(v=200,col="red")
#Remove functions with low coefficient of variation
sum(apply(pp_tr,2,co.var)>200)
pp_tr <- pp_tr[,apply(pp_tr,2,co.var)>200]

#Add patient response to dataset
pp_tr$Sample <- samples
pp_tr <- merge.easy(pp_tr,patient_meta,key="Sample")
pp_tr$Sample <- NULL

# #Profile metadata
# setwd("~/FunctionalNetworks_Project")
# pfam_feat<- read.delim("yams_pd1_pub/yams_pd1_pub_pfam_feat.txt",
#                          stringsAsFactors = FALSE)

#################
# Build Network #
#################

# #This takes a few hours! Just use the one I already built if you can
# pp.net <- spiec.easi(as.matrix(pp_tr), method='mb', lambda.min.ratio=1e-2,
#                      nlambda=15, pulsar.select = T,pulsar.params=list(rep.num=50),verbose=T)
# pp.net.ig <- adj2igraph(symBeta(getOptBeta(pp.net), mode='maxabs'))
# pp.net.ig$names <- colnames(pp_tr)
# save(pp.net,pp.net.ig,file="pfam_net_CV200_response.RData")

#Load pre-built network
load("pfam_net_CV200_response.RData")

#############
# Visualize #
#############

###Entire Network
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
pdf("FunctionalNetwork_pfam_spieceasimb_response.pdf",width=14,height=7)
par(mfrow=c(1,2))
plotnet(pp.net.ig.p, main="Functions: Positive Interactions")
plotnet(pp.net.ig.n, main="Functions: Negative Interactions")
par(mfrow=c(1,1))
dev.off()


##############
# CLustering #
##############

#### InfoMap Clustering ####

pp.net.ig.p.clustinfo <- infomap.community(pp.net.ig.p)

#Plot it (a mess)
agcoord <- layout_nicely(pp.net.ig.p)
pdf("FunctionNetworkPfam_InfoMapClustering_response.pdf",width=10,height=10, onefile=FALSE)
plot(pp.net.ig.p.clustinfo,pp.net.ig.p,layout=agcoord, 
     vertex.label=pp.net.ig.p$names,vertex.label.cex=0.25,
     vertex.size=4)
dev.off()

#Plot only largest clusters w/ color 
pdf("FunctionNetworkPfam_InfoMapClustering_topgroups_response.pdf",width=10,height=10, onefile=FALSE)
table(pp.net.ig.p.clustinfo$membership)
pp.net.ig.p.clustinfo.cut <- pp.net.ig.p.clustinfo
pp.net.ig.p.clustinfo.cut$membership[pp.net.ig.p.clustinfo.cut$membership>20] <- 0
plot(pp.net.ig.p,layout=agcoord, 
     vertex.label=pp.net.ig.p$names,vertex.label.cex=0.1,
     vertex.size=4,vertex.color=pp.net.ig.p.clustinfo.cut$membership)
dev.off()

#Put clusters into dataframe
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

write.csv(clust_df,file="neighborhood_detection_pfam_infomap_response.csv")

#Who shares a cluster with response?
clust_df$Cluster[clust_df$Accession=="ResponseBinary"]
clust_df$Accession[clust_df$Cluster==14]
clust_df$Description[clust_df$Cluster==14]
# [1] "TGF-beta propeptide"                                 "Septin"                                             
# [3] "Ezrin/radixin/moesin family"                         "TIR domain"                                         
# [5] "DNA directed RNA polymerase, 7 kDa subunit"          "Thiopurine S-methyltransferase (TPMT)"              
# [7] "Spc24 subunit of Ndc80"                              "LVIVD repeat"                                       
# [9] "A2L zinc ribbon domain"                              "Domain of unknown function (DUF1871)"               
# [11] "Domain of unknown function (DUF1897)"                "NurA domain"                                        
# [13] "Bpu10I restriction endonuclease"                     "MamI restriction endonuclease"                      
# [15] "Predicted membrane protein (DUF2127)"                "Protein of unknown function (DUF3072)"              
# [17] "Acetyl-CoA dehydrogenase C-terminal like"            "Bacterial transcriptional repressor"                
# [19] "Domain of unknown function (DUF4275)"                "Domain of unknown function (DUF4298)"               
# [21] "YaaC-like Protein"                                   "Soluble NSF attachment protein, SNAP"               
# [23] "Immunity protein 30"                               "HNH/Endo VII superfamily toxin with a SHH signature"
# [25] "HNH/Endo VII superfamily nuclease toxins"            "Domain of unknown function (DUF5048)"               
# [27] "Vacuolar (H+)-ATPase G subunit" 

#Note: "Immunity protein 30" is usually toxin adjacent


#### SpinGlass Clustering ####

pp.net.ig.clustspin <- spinglass.community(pp.net.ig)

#Plot it (a mess)
agcoord <- layout_nicely(pp.net.ig)
pdf("FunctionNetworkPfam_SpinClustering_response.pdf",width=10,height=10, onefile=FALSE)
plot(pp.net.ig.clustspin,pp.net.ig,layout=agcoord, 
     vertex.label=pp.net.ig$names,vertex.label.cex=0.25,
     vertex.size=4)
dev.off()

#Plot only largest clusters w/ color 
pdf("FunctionNetworkPfam_SpinClustering_topgroups_response.pdf",width=10,height=10, onefile=FALSE)
table(pp.net.ig.clustspin$membership)
pp.net.ig.clustspin.cut <- pp.net.ig.clustspin
pp.net.ig.clustspin.cut$membership[pp.net.ig.clustspin.cut$membership>20] <- 0
plot(pp.net.ig,layout=agcoord, 
     vertex.label=pp.net.ig$names,vertex.label.cex=0.1,
     vertex.size=4,vertex.color=pp.net.ig.clustspin.cut$membership)
dev.off()

#Put clusters into dataframe
clust_df <- data.frame(Cluster=1,
                       Accession=pp.net.ig$names[pp.net.ig.clustspin$membership==1],
                       stringsAsFactors = FALSE)
for(i in 2:max(pp.net.ig.clustspin$membership)){
  m_i <- data.frame(Cluster=i,
                    Accession=pp.net.ig$names[pp.net.ig.clustspin$membership==i],
                    stringsAsFactors = FALSE)
  clust_df <- rbind(clust_df,m_i)
}
clust_df <- merge.easy(clust_df,pfam_feat,key="Accession")

write.csv(clust_df,file="neighborhood_detection_pfam_spinglass_response.csv")

#Who shares a cluster with response?
clust_df$Cluster[clust_df$Accession=="ResponseBinary"]
clust_df$Accession[clust_df$Cluster==13]
clust_df$Description[clust_df$Cluster==13]

# [1] "Pyridine nucleotide-disulphide oxidoreductase"                       
# [2] "Ras family"                                                          
# [3] "Ubiquitin family"                                                    
# [4] "Paired box domain"                                                   
# [5] "Glycosyl hydrolases family 17"                                       
# [6] "Haemolysin-type calcium-binding repeat (2 copies)"                   
# [7] "Regulator of chromosome condensation (RCC1) repeat"                  
# [8] "CUB domain"                                                          
# [9] "Serine carboxypeptidase"                                             
# [10] "Glycosyl hydrolases family 11"                                       
# [11] "Diacylglycerol kinase accessory domain"                              
# [12] "MAM domain, meprin/A5/mu"                                            
# [13] "Zn-finger in Ran binding protein and others"                         
# [14] "Poly(ADP-ribose) polymerase catalytic domain"                        
# [15] "Sulfotransferase domain"                                             
# [16] "TGF-beta propeptide"                                                 
# [17] "Septin"                                                              
# [18] "Autoinducer synthase"                                                
# [19] "Ezrin/radixin/moesin family"                                         
# [20] "AP2 domain"                                                          
# [21] "SET domain"                                                          
# [22] "Ribosomal family S4e"                                                
# [23] "NB-ARC domain"                                                       
# [24] "Cellulose binding domain"                                            
# [25] "Peptidase S7, Flavivirus NS3 serine protease"                        
# [26] "Aerolysin toxin"                                                     
# [27] "NAD:arginine ADP-ribosyltransferase"                                 
# [28] "Ribosomal protein S8e"                                               
# [29] "Ribosomal protein L19e"                                              
# [30] "Astacin (Peptidase family M12A)"                                     
# [31] "Angiotensin-converting enzyme"                                       
# [32] "Jacalin-like lectin domain"                                          
# [33] "Apolipoprotein A1/A4/E domain"                                       
# [34] "Viral (Superfamily 1) RNA helicase"                                  
# [35] "Thermolysin metallopeptidase, catalytic domain"                      
# [36] "Pentapeptide repeats (8 copies)"                                     
# [37] "Transmembrane amino acid transporter protein"                        
# [38] "Transposase"                                                         
# [39] "TIR domain"                                                          
# [40] "Transposase"                                                         
# [41] "Isopentenyl transferase"                                             
# [42] "ATP-sulfurylase"                                                     
# [43] "MYND finger"                                                         
# [44] "Ribosomal proteins 50S-L18Ae/60S-L20/60S-L18A"                       
# [45] "PsbP"                                                                
# [46] "DUF35 OB-fold domain, acyl-CoA-associated"                           
# [47] "Penicillin amidase"                                                  
# [48] "MatK/TrnK amino terminal region"                                     
# [49] "Peptidase A4 family"                                                 
# [50] "Homocysteine biosynthesis enzyme, sulfur-incorporation"              
# [51] "FG-GAP repeat"                                                       
# [52] "CDP-archaeol synthase"                                               
# [53] "Archaeal holliday junction resolvase (hjc)"                          
# [54] "RNase P subunit p30"                                                 
# [55] "Putative O-antigen polymerase"                                       
# [56] "CRISPR-associated negative auto-regulator DevR/Csa2"                 
# [57] "Deoxyhypusine synthase"                                              
# [58] "Mut7-C RNAse domain"                                                 
# [59] "Integral membrane protein DUF92"                                     
# [60] "Translin family"                                                     
# [61] "N2,N2-dimethylguanosine tRNA methyltransferase"                      
# [62] "Glycosyl hydrolase family 48"                                        
# [63] "Hypothetical lipoprotein (MG045 family)"                             
# [64] "Streptomyces extracellular neutral proteinase (M7) family"           
# [65] "SAP domain"                                                          
# [66] "Glycosyl hydrolase family 59"                                        
# [67] "Phosducin"                                                           
# [68] "Piwi domain"                                                         
# [69] "Cytochrome D1 heme domain"                                           
# [70] "Acyl transferase"                                                    
# [71] "Transposase Tn5 dimerisation domain"                                 
# [72] "Capsid protein (F protein)"                                          
# [73] "RTX N-terminal domain"                                               
# [74] "Chlamydia polymorphic membrane protein (Chlamydia_PMP) repeat"       
# [75] "Lecithin:cholesterol acyltransferase"                                
# [76] "HYR domain"                                                          
# [77] "Auxin responsive protein"                                            
# [78] "Acyl-CoA thioesterase"                                               
# [79] "Thermolysin metallopeptidase, alpha-helical domain"                  
# [80] "Restriction endonuclease BamHI"                                      
# [81] "Ecdysteroid kinase"                                                  
# [82] "Restriction endonuclease EcoRI"                                      
# [83] "Conserved hypothetical ATP binding protein"                          
# [84] "NLI interacting factor-like phosphatase"                             
# [85] "ARD/ARD family"                                                      
# [86] "LAGLIDADG DNA endonuclease family"                                   
# [87] "Late embryogenesis abundant protein"                                 
# [88] "PUCC protein"                                                        
# [89] "Ubiquinone biosynthesis protein COQ7"                                
# [90] "Pectinacetylesterase"                                                
# [91] "Snf7"                                                                
# [92] "Prohead core protein serine protease"                                
# [93] "Anp1"                                                                
# [94] "eRF1 domain 3"                                                       
# [95] "Glycoside-hydrolase family GH114"                                    
# [96] "UL73 viral envelope glycoprotein"                                    
# [97] "DNA directed RNA polymerase, 7 kDa subunit"                          
# [98] "Glycosyl hydrolase family 81"                                        
# [99] "Glycosyl hydrolase family 85"                                        
# [100] "Glycosyl hydrolase family 79, N-terminal domain"                     
# [101] "Glycosyl hydrolase family 49"                                        
# [102] "RNA polymerase Rpb4"                                                 
# [103] "Alg9-like mannosyltransferase family"                                
# [104] "delta endotoxin"                                                     
# [105] "Domain of unknown function (DUF366)"                                 
# [106] "Class III signal peptide"                                            
# [107] "Fatty acid hydroxylase superfamily"                                  
# [108] "Putative cell wall binding repeat 2"                                 
# [109] "B-block binding subunit of TFIIIC"                                   
# [110] "Mannosyltransferase (PIG-V)"                                         
# [111] "CENP-B N-terminal DNA-binding domain"                                
# [112] "Protein of unknown function (DUF433)"                                
# [113] "Phosphomevalonate kinase"                                            
# [114] "Protein of unknown function (DUF443)"                                
# [115] "Ribonuclease H-like"                                                 
# [116] "Protein of unknown function (DUF475)"                                
# [117] "Arginine-tRNA-protein transferase, C terminus"                       
# [118] "PD-(D/E)XK nuclease superfamily"                                     
# [119] "tRNAHis guanylyltransferase"                                         
# [120] "Protein of unknown function (DUF499)"                                
# [121] "Transglutaminase-like domain"                                        
# [122] "Restriction endonuclease XhoI"                                       
# [123] "Mitochondrial ATPase inhibitor, IATP"                                
# [124] "Phage tail tube protein"                                             
# [125] "Opioid growth factor receptor (OGFr) conserved region"               
# [126] "Archaeal putative transposase ISC1217"                               
# [127] "FaeA-like protein"                                                   
# [128] "Glycosyltransferase family 17"                                       
# [129] "Coatomer epsilon subunit"                                            
# [130] "Lantibiotic dehydratase, C terminus"                                 
# [131] "Sec23/Sec24 zinc finger"                                             
# [132] "PT repeat"                                                           
# [133] "Phospholipase B"                                                     
# [134] "Poxvirus A51 protein"                                                
# [135] "Flp/Fap pilin component"                                             
# [136] "HPP family"                                                          
# [137] "Phytochelatin synthase"                                              
# [138] "Sucrose-6F-phosphate phosphohydrolase"                               
# [139] "Hypothetical methyltransferase"                                      
# [140] "Electron transfer flavoprotein-ubiquinone oxidoreductase, 4Fe-4S"    
# [141] "Hom_end-associated Hint"                                             
# [142] "Homing endonuclease"                                                 
# [143] "Septum formation inhibitor MinC, N-terminal domain"                  
# [144] "Protein of unknown function (DUF707)"                                
# [145] "helix-turn-helix, Psq domain"                                        
# [146] "Alpha-L-arabinofuranosidase B (ABFB) domain"                         
# [147] "ICEA Protein"                                                        
# [148] "Phage major coat protein, Gp8"                                       
# [149] "Phi-29-like late genes activator (early protein GP4)"                
# [150] "Photosystem I reaction centre subunit N (PSAN or PSI-N)"             
# [151] "Cbb3-type cytochrome oxidase component FixQ"                         
# [152] "Protein of unknown function (DUF763)"                                
# [153] "Bacterial TniB protein"                                              
# [154] "Protein of unknown function (DUF789)"                                
# [155] "Putative bacterial lipoprotein (DUF799)"                             
# [156] "Zonular occludens toxin (Zot)"                                       
# [157] "Golgi phosphoprotein 3 (GPP34)"                                      
# [158] "Thiopurine S-methyltransferase (TPMT)"                               
# [159] "NACHT domain"                                                        
# [160] "Glutaredoxin-like domain (DUF836)"                                   
# [161] "Eukaryotic protein of unknown function (DUF842)"                     
# [162] "Glycosyl hydrolase 108"                                              
# [163] "O-phosphoseryl-tRNA(Sec) selenium transferase, SepSecS"              
# [164] "Protein of unknown function (DUF861)"                                
# [165] "Dam-replacing family"                                                
# [166] "L-asparaginase II"                                                   
# [167] "Protein of unknown function (DUF952)"                                
# [168] "Staphylococcal protein of unknown function (DUF960)"                 
# [169] "HrpE/YscL/FliH and V-type ATPase subunit E"                          
# [170] "Carbon monoxide dehydrogenase subunit G (CoxG)"                      
# [171] "Protein of unknown function (DUF1048)"                               
# [172] "Beta-1,4-N-acetylgalactosaminyltransferase (CgtA)"                   
# [173] "Acetoacetate decarboxylase (ADC)"                                    
# [174] "DNA repair protein MmcB-like"                                        
# [175] "Conjugal transfer protein TraD"                                      
# [176] "MucBP domain"                                                        
# [177] "TniQ"                                                                
# [178] "2-deoxycytidine 5-triphosphate deaminase (DCD)"                      
# [179] "Haemolysin-type calcium binding protein related domain"              
# [180] "Protein of unknown function (DUF1152)"                               
# [181] "Enterobacterial exodeoxyribonuclease VIII"                           
# [182] "Protein of unknown function (DUF1156)"                               
# [183] "Amino acid synthesis"                                                
# [184] "Protein of unknown function (DUF1186)"                               
# [185] "Protein of unknown function (DUF1187)"                               
# [186] "Helical and beta-bridge domain"                                      
# [187] "Serine hydrolase"                                                    
# [188] "Pilin accessory protein (PilO)"                                      
# [189] "Protein of unknown function (DUF1257)"                               
# [190] "Methionine biosynthesis protein MetW"                                
# [191] "Major capsid protein Gp23"                                           
# [192] "Protein of unknown function (DUF1344)"                               
# [193] "Domain of unknown function (DUF1413)"                                
# [194] "Alcohol acetyltransferase"                                           
# [195] "Protein of unknown function (DUF1450)"                               
# [196] "NUMOD1 domain"                                                       
# [197] "NUMOD3 motif (2 copies)"                                             
# [198] "Beta-lactamase inhibitor (BLIP)"                                     
# [199] "Fungalysin/Thermolysin Propeptide Motif"                             
# [200] "WavE lipopolysaccharide synthesis"                                   
# [201] "Putative helicase"                                                   
# [202] "Lipid A Biosynthesis N-terminal domain"                              
# [203] "Rhodopirellula transposase DDE domain"                               
# [204] "YKOF-related Family"                                                 
# [205] "Kelch motif"                                                         
# [206] "Immunoglobulin I-set domain"                                         
# [207] "7TM diverse intracellular signalling"                                
# [208] "Basic region leucine zipper"                                         
# [209] "Chorismate mutase type I"                                            
# [210] "Protein of unknown function (DUF1611)"                               
# [211] "Transmembrane protein 43"                                            
# [212] "Archaeal Type IV pilin, N-terminal"                                  
# [213] "Protein of unknown function (DUF1634)"                               
# [214] "Calcium binding and coiled-coil domain (CALCOCO1) like"              
# [215] "gp58-like protein"                                                   
# [216] "Peptidase family M54"                                                
# [217] "PHR domain"                                                          
# [218] "Bacteriophage protein GP30.3"                                        
# [219] "FAD-binding domain"                                                  
# [220] "Histone methylation protein DOT1"                                    
# [221] "Propeptide_C25"                                                      
# [222] "Cache domain"                                                        
# [223] "TFIIB zinc-binding"                                                  
# [224] "M protein trans-acting positive regulator (MGA) HTH domain"          
# [225] "Spc24 subunit of Ndc80"                                              
# [226] "LVIVD repeat"                                                        
# [227] "Activator of Hsp90 ATPase homolog 1-like protein"                    
# [228] "QacR-like protein, C-terminal region"                                
# [229] "FAE1/Type III polyketide synthase-like protein"                      
# [230] "CbbQ/NirQ/NorQ C-terminal"                                           
# [231] "Domain of unknown function (DUF1737)"                                
# [232] "Rad51"                                                               
# [233] "Chromatin associated protein KTI12"                                  
# [234] "Anaphase promoting complex (APC) subunit 2"                          
# [235] "DNA replication factor Dna2"                                         
# [236] "Domain of unknown function (DUF1794)"                                
# [237] "ParB family"                                                         
# [238] "Alginate lyase"                                                      
# [239] "A2L zinc ribbon domain"                                              
# [240] "Phage related hypothetical protein (DUF1799)"                        
# [241] "Protein of unknown function (DUF1804)"                               
# [242] "Domain of unknown function (DUF1904)"                                
# [243] "Domain of unknown function (DUF1871)"                                
# [244] "Protein of unknown function (DUF1878)"                               
# [245] "Domain of unknown function (DUF1902)"                                
# [246] "Cytotoxic"                                                           
# [247] "Domain of unknown function (DUF1887)"                                
# [248] "Domain of unknown function (DUF1897)"                                
# [249] "NgoMIV restriction enzyme"                                           
# [250] "Restriction endonuclease BsobI"                                      
# [251] "Carbohydrate binding module 27"                                      
# [252] "Restriction endonuclease EcoRII, N-terminal"                         
# [253] "Domain of unknown function (DUF1961)"                                
# [254] "Restriction endonuclease PvuII"                                      
# [255] "Restriction endonuclease EcoRV"                                      
# [256] "Interferon-alpha/beta receptor, fibronectin type III"                
# [257] "Domain of unknown function (DUF1989)"                                
# [258] "Phasin protein"                                                      
# [259] "NurA domain"                                                         
# [260] "Methylmuconolactone methyl-isomerase"                                
# [261] "Domain of unknown function (DUF2019)"                                
# [262] "Bsp6I restriction endonuclease"                                      
# [263] "RNA ligase"                                                          
# [264] "Eco29kI restriction endonuclease"                                    
# [265] "HindIII restriction endonuclease"                                    
# [266] "HindVP restriction endonuclease"                                     
# [267] "Type II restriction endonuclease, TdeIII"                            
# [268] "NgoPII restriction endonuclease"                                     
# [269] "Conserved phage C-terminus (Phg_2220_C)"                             
# [270] "Bpu10I restriction endonuclease"                                     
# [271] "HaeII restriction endonuclease"                                      
# [272] "HaeIII restriction endonuclease"                                     
# [273] "NgoBV restriction endonuclease"                                      
# [274] "SacI restriction endonuclease"                                       
# [275] "MamI restriction endonuclease"                                       
# [276] "MjaI restriction endonuclease"                                       
# [277] "ScaI restriction endonuclease"                                       
# [278] "XamI restriction endonuclease"                                       
# [279] "Sporulation lipoprotein YhcN/YlaJ (Spore_YhcN_YlaJ)"                 
# [280] "Conserved hypothetical protein (Lin0512_fam)"                        
# [281] "Bacterial protein of unknown function (HtrL_YibB)"                   
# [282] "CRISPR-associated protein Csx8 (Cas_Csx8)"                           
# [283] "Type II restriction endonuclease (RE_Alw26IDE)"                      
# [284] "CRISPR-associated protein (Cas_Cas02710)"                            
# [285] "Chlamydia-phage Chp2 scaffold (Chlamy_scaf)"                         
# [286] "Bacteriocin (Lactococcin_972)"                                       
# [287] "CRISPR-associated protein (Cas_Cmr3)"                                
# [288] "CRISPR-associated protein (Cas_Cmr5)"                                
# [289] "CRISPR-associated protein (Cas_CXXC_CXXC)"                           
# [290] "Protein of unknown function (DUF2384)"                               
# [291] "C-terminal AAA-associated domain"                                    
# [292] "Uncharacterized protein conserved in bacteria (DUF2066)"             
# [293] "Uncharacterized protein conserved in archaea (DUF2099)"              
# [294] "Uncharacterized protein containing a Zn-ribbon (DUF2116)"            
# [295] "Predicted membrane protein (DUF2127)"                                
# [296] "Uncharacterized protein conserved in bacteria (DUF2147)"             
# [297] "Putative PD-(D/E)XK phosphodiesterase (DUF2161)"                     
# [298] "Predicted membrane protein (DUF2178)"                                
# [299] "Transcriptional regulator, AbiEi antitoxin, Type IV TA system"       
# [300] "Uncharacterized protein conserved in bacteria (DUF2188)"             
# [301] "Uncharacterized protein conserved in bacteria (DUF2199)"             
# [302] "Nucleotidyl transferase of unknown function (DUF2204)"               
# [303] "Uncharacterized alpha/beta hydrolase domain (DUF2235)"               
# [304] "Uncharacterized protein conserved in bacteria (DUF2247)"             
# [305] "Uncharacterized conserved protein (DUF2290)"                         
# [306] "Uncharacterized small protein (DUF2292)"                             
# [307] "Uncharacterized protein conserved in bacteria (DUF2316)"             
# [308] "Uncharacterized protein conserved in bacteria (DUF2321)"             
# [309] "Uncharacterized protein conserved in bacteria (DUF2334)"             
# [310] "Predicted membrane protein (DUF2335)"                                
# [311] "Glycosyltransferase like family 2"                                   
# [312] "OpgC protein"                                                        
# [313] "vWA found in TerF C terminus"                                        
# [314] "PhoPQ-activated pathogenicity-related protein"                       
# [315] "Mitochondrial ribosomal death-associated protein 3"                  
# [316] "Predicted transmembrane and coiled-coil 2 protein"                   
# [317] "Transmembrane Fragile-X-F protein"                                   
# [318] "Serpentine type 7TM GPCR chemoreceptor Sri"                          
# [319] "Putative cell-wall binding lipoprotein"                              
# [320] "Protein of unknown function (DUF2441)"                               
# [321] "P22_AR N-terminal domain"                                            
# [322] "Uncharacterised protein family UPF0547"                              
# [323] "Sporulation inhibitor of replication protein SirA"                   
# [324] "Protein of unknown function (DUF2586)"                               
# [325] "Haemolysin XhlA"                                                     
# [326] "Protein of unknown function (DUF2716)"                               
# [327] "Protein of unknown function (DUF2691)"                               
# [328] "Protein of unknown function (DUF2695)"                               
# [329] "Protein of unknown function (DUF2712)"                               
# [330] "Protein of unknown function (DUF2705)"                               
# [331] "Protein of unknown function (DUF2806)"                               
# [332] "Protein of unknown function (DUF2612)"                               
# [333] "Protein of unknown function (DUF2750)"                               
# [334] "Plectrovirus spv1-c74 ORF 12 transmembrane protein"                  
# [335] "Antirestriction protein Ral"                                         
# [336] "Protein of unknown function (DUF2624)"                               
# [337] "Protein of unknown function (DUF2931)"                               
# [338] "Phosphoribosyl transferase (PRTase)"                                 
# [339] "Protein of unknown function (DUF2993)"                               
# [340] "Protein of unknown function (DUF2997)"                               
# [341] "Protein of unknown function (DUF3040)"                               
# [342] "Protein of unknown function (DUF3060)"                               
# [343] "Protein of unknown function (DUF3072)"                               
# [344] "Protein of unknown function (DUF3138)"                               
# [345] "Holin of 3TMs, for gene-transfer release"                            
# [346] "Protein of unknown function (DUF3192)"                               
# [347] "Type II restriction enzyme MunI"                                     
# [348] "Hypoxia-inducible factor-1"                                          
# [349] "Protein of unknown function (DUF3226)"                               
# [350] "Protein of unknown function (DUF3232)"                               
# [351] "Restriction endonuclease BpuJI - N terminal"                         
# [352] "Protein of unknown function (DUF3237)"                               
# [353] "PD-(D/E)XK endonuclease"                                             
# [354] "Protein of unknown function (DUF3277)"                               
# [355] "CGNR zinc finger"                                                    
# [356] "Protein of unknown function (DUF3304)"                               
# [357] "Keratin-associated matrix"                                           
# [358] "Alanine-zipper, major outer membrane lipoprotein"                    
# [359] "Domain of unknown function (DUF3372)"                                
# [360] "Putative outer membrane beta-barrel porin, MtrB/PioB"                
# [361] "Domain of unknown function (DUF3394)"                                
# [362] "Domain of unknown function (DUF3471)"                                
# [363] "Peptidase_G2, IMC autoproteolytic cleavage domain"                   
# [364] "Protein of unknown function (DUF3489)"                               
# [365] "Domain of unknown function (DUF3536)"                                
# [366] "Domain of unknown function (DUF3578)"                                
# [367] "SprA-related family"                                                 
# [368] "Adenosine-5-phosphosulfate reductase beta subunit"                   
# [369] "Methanol-cobalamin methyltransferase B subunit"                      
# [370] "Restriction endonuclease NotI"                                       
# [371] "NAD(P)H binding domain of trans-2-enoyl-CoA reductase"               
# [372] "Insecticide toxin TcdB middle/N-terminal region"                     
# [373] "Protein of unknown function (DUF3644)"                               
# [374] "CopG antitoxin of type II toxin-antitoxin system"                    
# [375] "CRISPR-associated protein"                                           
# [376] "TRSP domain C terminus to PRTase_2"                                  
# [377] "Protein of unknown function (DUF3732)"                               
# [378] "Type III restriction/modification enzyme methylation subunit"        
# [379] "Protein of unknown function (DUF3761)"                               
# [380] "Protein of unknown function (DUF3770)"                               
# [381] "Iron-Sulfur binding protein C terminal"                              
# [382] "Arabinose-binding domain of AraC transcription regulator, N-term"    
# [383] "DNase/tRNase domain of colicin-like bacteriocin"                     
# [384] "UPF0489 domain"                                                      
# [385] "Domain of unknown function (DUF3787)"                                
# [386] "Phage tail repeat like"                                              
# [387] "4Fe-4S binding domain"                                               
# [388] "Acetyl-CoA dehydrogenase C-terminal like"                            
# [389] "Glycine rich protein"                                                
# [390] "tRNA_anti-like"                                                      
# [391] "Glycoside hydrolase family 44"                                       
# [392] "Putative lumazine-binding"                                           
# [393] "Protein of unknown function (DUF3830)"                               
# [394] "Antitoxin of toxin-antitoxin, RelE / RelB, TA system"                
# [395] "Sec23-binding domain of Sec16"                                       
# [396] "Domain of Unknown Function with PDB structure (DUF3862)"             
# [397] "Acetyl-coenzyme A transporter 1"                                     
# [398] "leucine-zipper of insertion element IS481"                           
# [399] "Domain of unknown function (DUF3893)"                                
# [400] "Domain of unknown function (DUF3899)"                                
# [401] "Protein of unknown function (DUF3916)"                               
# [402] "Protein of unknown function (DUF3970)"                               
# [403] "Protein of unknown function (DUF3997)"                               
# [404] "Protein of unknown function (DUF4013)"                               
# [405] "AAA domain"                                                          
# [406] "Domain of unknown function (DUF4041)"                                
# [407] "Protein of unknown function (DUF4043)"                               
# [408] "Domain of unknown function (DUF4082)"                                
# [409] "Tetratricopeptide repeat"                                            
# [410] "Hint domain"                                                         
# [411] "Domain of unknown function (DUF4112)"                                
# [412] "N-terminal half of MaoC dehydratase"                                 
# [413] "Transcription factor zinc-finger"                                    
# [414] "Meiotically up-regulated gene 113"                                   
# [415] "4Fe-4S single cluster domain"                                        
# [416] "Transglutaminase-like superfamily"                                   
# [417] "Reductive dehalogenase subunit"                                      
# [418] "Sugar nucleotidyl transferase"                                       
# [419] "Pentapeptide repeats (9 copies)"                                     
# [420] "Ankyrin repeat"                                                      
# [421] "KTSC domain"                                                         
# [422] "Metallo-peptidase family M12"                                        
# [423] "Domain of unknown function (DUF4159)"                                
# [424] "AIG2-like family"                                                    
# [425] "Tubulin like"                                                        
# [426] "Thioredoxin-like domain"                                             
# [427] "Immunoglobulin domain"                                               
# [428] "Immunoglobulin domain"                                               
# [429] "Domain of unknown function (DUF4214)"                                
# [430] "Bacterial transcriptional repressor"                                 
# [431] "Domain of unknown function (DUF4224)"                                
# [432] "Protein of unknown function (DUF4225)"                               
# [433] "Protein of unknown function (DUF4237)"                               
# [434] "YlzJ-like protein"                                                   
# [435] "Domain of unknown function (DUF4262)"                                
# [436] "Domain of unknown function (DUF4274)"                                
# [437] "Domain of unknown function (DUF4275)"                                
# [438] "Domain of unknown function (DUF4280)"                                
# [439] "Domain of unknown function (DUF4297)"                                
# [440] "Domain of unknown function (DUF4298)"                                
# [441] "Domain of unknown function (DUF4303)"                                
# [442] "Super-infection exclusion protein B"                                 
# [443] "YfkD-like protein"                                                   
# [444] "YjzC-like protein"                                                   
# [445] "YaaC-like Protein"                                                   
# [446] "Domain of unknown function (DUF4338)"                                
# [447] "Domain of unknown function (DUF4351)"                                
# [448] "Putative component of biosynthetic module"                           
# [449] "Arylsulfotransferase (ASST)"                                         
# [450] "Domain of unknown function (DUF4359)"                                
# [451] "Domain of unknown function (DUF4362)"                                
# [452] "Domain of unknown function (DUF4384)"                                
# [453] "Domain of unknown function (DUF4393)"                                
# [454] "Domain of unknown function (DUF4394)"                                
# [455] "Domain of unknown function (DUF4402)"                                
# [456] "Domain of unknown function (DUF4407)"                                
# [457] "Domain of unknown function (DUF4412)"                                
# [458] "Superinfection immunity protein"                                     
# [459] "HNH/ENDO VII superfamily nuclease with conserved GHE residues"       
# [460] "A nuclease of the HNH/ENDO VII superfamily with conserved LHH"       
# [461] "Thg1 C terminal domain"                                              
# [462] "A nuclease of the HNH/ENDO VII superfamily with conserved WHH"       
# [463] "Immunity protein Imm5"                                               
# [464] "The  BURPS668_1122 family of deaminases"                             
# [465] "Immunity protein Imm1"                                               
# [466] "Bacterial EndoU nuclease"                                            
# [467] "Pre-toxin TG"                                                        
# [468] "Multiubiquitin"                                                      
# [469] "ThiS-like ubiquitin"                                                 
# [470] "Prokaryotic E2 family A"                                             
# [471] "Prokaryotic E2 family B"                                             
# [472] "Glutathione S-transferase, C-terminal domain"                        
# [473] "CAP-associated N-terminal"                                           
# [474] "Type II restriction endonuclease EcoO109I"                           
# [475] "AAA-like domain"                                                     
# [476] "Syntaxin-like protein"                                               
# [477] "Endonuclease-reverse transcriptase"                                  
# [478] "SMI1-KNR4 cell-wall"                                                 
# [479] "Oligogalacturonate lyase"                                            
# [480] "DNA ligase OB-like domain"                                           
# [481] "Photosynthesis system II assembly factor YCF48"                      
# [482] "Intein splicing domain"                                              
# [483] "Soluble NSF attachment protein, SNAP"                                
# [484] "Winged helix-turn-helix"                                             
# [485] "TRP-interacting helix"                                               
# [486] "TMEM154 protein family"                                              
# [487] "Immunity protein 26"                                                 
# [488] "Novel toxin 10"                                                      
# [489] "Novel toxin 21"                                                      
# [490] "Bacterial toxin 30"                                                  
# [491] "Bacterial toxin 33"                                                  
# [492] "Bacterial toxin 37"                                                  
# [493] "Bacterial toxin 50"                                                  
# [494] "Immunity protein 12"                                                 
# [495] "Immunity protein 30"                                                 
# [496] "Immunity protein 40"                                                 
# [497] "Immunity protein 43"                                                 
# [498] "Immunity protein 45"                                                 
# [499] "Immunity protein 7"                                                  
# [500] "Immunity protein 10"                                                 
# [501] "Immunity protein 42"                                                 
# [502] "Immunity protein 50"                                                 
# [503] "Immunity protein 63"                                                 
# [504] "Immunity protein 70"                                                 
# [505] "Immunity protein 74"                                                 
# [506] "Novel toxin 15"                                                      
# [507] "Bacterial toxin 28"                                                  
# [508] "Bacterial toxin 44"                                                  
# [509] "Phosphoribosyl transferase"                                          
# [510] "EH_Signature domain"                                                 
# [511] "NTF2 fold immunity protein"                                          
# [512] "ATP-grasp in the biosynthetic pathway with Ter operon"               
# [513] "HYD1 signature containing ADP-ribosyltransferase"                    
# [514] "Restriction endonuclease fold toxin 7"                               
# [515] "Restriction endonuclease fold toxin 9"                               
# [516] "HNH/Endo VII superfamily toxin with a SHH signature"                 
# [517] "URI fold toxin 2"                                                    
# [518] "Toxin with a H, D/N and C signature"                                 
# [519] "HNH/Endo VII superfamily nuclease toxins"                            
# [520] "Interleukin-like EMT inducer"                                        
# [521] "Sortilin, neurotensin receptor 3,"                                   
# [522] "Domain of unknown function (DUF4747)"                                
# [523] "Domain of unknown function (DUF4760)"                                
# [524] "HicB_like antitoxin of bacterial toxin-antitoxin system"             
# [525] "Domain of unknown function (DUF4810)"                                
# [526] "Ig-like domain from next to BRCA1 gene"                              
# [527] "Cadherin-like"                                                       
# [528] "Domain of unknown function (DUF4879)"                                
# [529] "Domain of unknown function (DUF4914)"                                
# [530] "Domain of unknown function (DUF4926)"                                
# [531] "Domain of unknown function (DUF4928)"                                
# [532] "Domain of unknown function (DUF4935)"                                
# [533] "Domain of unknown function (DUF5021)"                                
# [534] "Domain of unknown function (DUF5023)"                                
# [535] "Domain of unknown function (DUF5026)"                                
# [536] "Domain of unknown function (DUF5027)"                                
# [537] "Domain of unknown function (DUF5028)"                                
# [538] "Domain of unknown function (DUF5036)"                                
# [539] "Domain of unknown function (DUF5048)"                                
# [540] "Domain of unknown function (DUF5052)"                                
# [541] "Domain of unknown function (DUF5055)"                                
# [542] "Beta-1,3-glucanase"                                                  
# [543] "C-terminal leucine zipper domain of cyclic nucleotide-gated channels"
# [544] "Flagellar assembly protein T, C-terminal domain"                     
# [545] "Gram-positive pilin subunit D1, N-terminal"                          
# [546] "Bacterial Ig-like domain (group 3)"                                  
# [547] "Tubulin-specific chaperone C N-terminal domain"                      
# [548] "IseA DL-endopeptidase inhibitor"                                     
# [549] "Domain of unknown function (DUF5071)"                                
# [550] "Domain of unknown function (DUF5072)"                                
# [551] "CRISPR-associated protein Csn2 subfamily St"                         
# [552] "Dimerisation domain"                                                 
# [553] "Putative phage abortive infection protein"                           
# [554] "Putative abortive phage resistance protein AbiGii toxin"             
# [555] "Domain of unknown function (DUF5081)"                                
# [556] "Vacuolar (H+)-ATPase G subunit" 
