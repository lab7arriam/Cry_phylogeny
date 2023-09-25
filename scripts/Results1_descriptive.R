library(ggplot2)
library(ggsci)
library(dplyr)
library(wesanderson)

setwd("../data_for_scripts/for_viz/descriptive_stat")


#Total amount of toxins
#Make the dataframe with the number of novel and known toxins from different datasets 
tox_amonut_df <- data.frame(dataset= c('IPG','GenBank','Dedup_tox','Bt_nom','Assemblies', 'Ref_clust', 
                                       'IPG','GenBank','Dedup_tox','Bt_nom','Assemblies', 'Ref_clust'),
                            flag =c('all','all','all','all','all','all',
                                    'novel','novel','novel','novel','novel','novel') , 
                            num = c(1058,1003,733,642,602,368,421,249,304,0,212,127))

#Vizualize the amount as a bar plots
ggplot(tox_amonut_df, aes(x=reorder(dataset,-num), y=num,  fill=flag)) +
  geom_bar(stat='identity',alpha=0.99, col='black')+
  theme_bw()+scale_fill_manual(labels = c("All", "Novel"), values = c("#DC143C", '#5F9EA0' ))+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=16), 
         legend.text=element_text(size=14),
         legend.position = c(0.77, 0.76)
  )+
  xlab('Dataset') +
  ylab('Number of toxins')+
  guides(fill= guide_legend(title="Toxins' group"))

#Get the data on the number of novel proteins and their identitty with the closest homolog from the BPPRC
all_nov_tox <- read.csv('all_NOV_ids.csv', header = F, sep = "\t", stringsAsFactors = F)
mean(all_nov_tox$V2) #82.62171

ref_nov_tox <- read.csv('ref_NOV_ids.csv', header = F, sep = "\t", stringsAsFactors = F)
mean(ref_nov_tox$V2) #63.63071 

ref_nov_tox$toxin_group <- 'Reference clusters'
all_nov_tox$toxin_group <- 'All toxins'

#Prepare merged dataframe with identities 
stacked_df <- rbind(all_nov_tox, ref_nov_tox)

#Plot the number of toxins according to the identity
ggplot(stacked_df, aes(x=V2, fill=toxin_group)) +
  geom_histogram(color='black', alpha=0.9)+
  theme_bw()+scale_fill_manual(values = c("#DC143C", '#5F9EA0' ))+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12))+
  xlab('Identity percent') +
  ylab('Frequency')+
  guides(fill= guide_legend(title="Toxins' group"))



#Clusterization patterns
#Read the file with clusters' properties
clusters_identidy <-  read.csv('clusters_identity.csv', header = T, sep = "\t") 

#Extract required data
clusters_df <- data.frame(num=rep(clusters_identidy$num_tox,3),
                          doms=c(rep('domain1',368),rep('domain2',368),rep('domain3',368)),
                          id = c(clusters_identidy$dom1_id,clusters_identidy$dom2_id,clusters_identidy$dom3_id))
#Plot within-cluster identity
ggplot(clusters_df[clusters_df$num>1,], aes(id,fill=doms)) + 
  theme_bw()+
  geom_density(alpha=0.6)+
  scale_fill_manual(labels = c("domain1",'domain2', 'domain3'),values=c("#a60b0b", '#2980b9' ,'#bdbd00'))+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=16), 
         legend.text=element_text(size=14))+
  xlab('Identity percent') +
  ylab('Frequency density')+
  facet_wrap(~doms)+
  guides(fill= guide_legend(title="Domain"))

#Extract the data with domain similarities
clusters_df <- data.frame(num=rep(clusters_identidy$num_tox,3),
                          doms=c(rep('domain1',368),rep('domain2',368),rep('domain3',368)),
                          id = c(clusters_identidy$dom1_id,clusters_identidy$dom2_id,clusters_identidy$dom3_id))

#Vizualize mean domain-wise intra-cluster identity vs the number of proteins in the cluster
ggplot(clusters_df, aes(x=num,y= id, size=num, shape=doms, color='ref')) +
  geom_point(alpha=0.75)+ #color='#2E98CA' 
  theme_bw()+xlab('Number of toxins')+ylab('Mean identity')+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14))+
  ylim(c(80,100))+
  guides( shape=guide_legend(override.aes = list(size=5), title="Domain"), 
          size=guide_legend(override.aes = list(size=5), title="Number of toxins"))


#The frequency of clusters according to the number of proteins
cluster_counts <- clusters_identidy[,2] %>% table() %>% as.data.frame()
colnames(cluster_counts)[1] <- 'Num_tox' 

#Vizualize the frequency of clusters' abundance
ggplot(cluster_counts, aes(x=Num_tox,y=Freq,fill='#2E98CA')) +
  geom_bar(stat="identity",alpha=0.55, col='black',width = 1)+ 
  theme_bw()+xlab('Number of toxins')+ylab('Frequency')+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.position = 'none')


#Statistics in the multi-clusters
#Get the data with domain-wise identity comparisions for clusters with more then 10 toxins
multiclust_ids <- read.csv('multiclusters_id_per_domain.csv', header = T, sep = "\t", stringsAsFactors = F)
id_stat_multi <- multiclust_ids %>% group_by(domain) %>% summarize(mean_id=mean(id))
# domain1  -  99.0 domain2  -  99.0 domain3  -  99.2

#Plot total distribution of pair-wise identity
ggplot(multiclust_ids, aes(x=id, fill=domain)) +
  geom_density(alpha=0.8)+
  theme_bw()+scale_fill_manual(values = c('#a60b0b','#2980b9', '#bdbd00'))+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12))+
  xlab("Identity") +
  ylab('Frequency')+
  guides(fill= guide_legend(title="Domain"))

#Read the table with coordinate wise classification of mutations in multi-clusters
multiclust_mism <- read.csv('multiclusters_mismacth_sites.csv', header = T, sep = "\t", stringsAsFactors = F)

#Plot the types of the mutations per relative domain coordinate in the alignments
ggplot(multiclust_mism , aes(x=coord, y=cluster, fill=site_type)) +
  geom_point(shape=21, alpha=0.7, color='black', size=3)+
  geom_vline(xintercept=100, linetype="dashed", color = "darkgreen",size=2, alpha = 0.7)+
  geom_vline(xintercept=200, linetype="dashed", color = "darkgreen",size=2, alpha = 0.7)+ 
  theme_bw()+
  #scale_fill_manual(labels = c("Indel",'Non-syn','Indel+Non-syn', 'Syn'), values = c('#bdbd00','#2980b9','#a60b0b','grey'))+
  scale_fill_manual(labels = c("Indel",'Non-syn','Syn'), values = c('#042AB6','#D7650C','#F0ECF0'))+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=16), 
         legend.text=element_text(size=14))+
  xlab('Coordinate') +
  ylab('Cluster')+
  guides(fill= guide_legend(title=""))


#Gather the data on the number of certain mutations and the abundace in the cluster for heatmap
multiclust_for_heatmap <- multiclust_mism  %>%  group_by(cluster, domain, site_type) %>%  tally()
total_num_heatmap <- data.frame(cluster = c('BT_Cry1Ac39','BT_Cry2Ab15','BT_Cry1Ab12','BT_Cry1Ia11','BT_Cry1Aa10','BT_Cry2Ac10','BT_Cry2Aa13'), 
                                Num = c(34,32,30,29,17,13, 14), fac=rep('num',7))

#Plot parts of the heatmap with the number of mutations
#ggplot(multiclust_for_heatmap , aes(x=site_type, y=cluster, size=n, fill=site_type)) + 
ggplot(total_num_heatmap , aes(x=fac, y=cluster, size=Num, fill=fac))+
  geom_point(shape=21, alpha=0.7, color='black')+ 
  theme_bw()+
  #scale_fill_manual(labels = c("Indel",'Non-syn','Syn'), values = c('#042AB6','#D7650C','#F0ECF0'))+
  scale_fill_manual(values='#E4718B')+
  theme( axis.text.x = element_text(color='black', 
                                    size=12, angle=45, hjust = 1),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=16), 
         legend.text=element_text(size=14))+
  xlab('Abundance') +
  ylab('Cluster')+
  guides(fill= guide_legend(title=""))+
  #facet_wrap(~domain)+
  coord_fixed(ratio = 3.9)

#Summarize and plot total number of mutation per type
sites_stat <- multiclust_mism[, c(1,2,4,5)] %>%  group_by(domain, site_type) %>%  tally()

#Plot the total number of certain sites
ggplot(sites_stat, aes(x=reorder(site_type,-n), y=n,  fill=site_type )) +
  geom_bar(stat='identity',alpha=0.99, col='black')+
  theme_bw()+scale_fill_manual(labels = c("Indel",'Non-syn','Syn'), 
                               values = c('#042AB6','#D7650C','#F0ECF0'))+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=16), 
         legend.text=element_text(size=14)
  )+
  xlab('Site type') +
  ylab('Number of sites')+
  guides(fill= guide_legend(title="Site type"))+
  facet_wrap(~domain)

#Summarize data for calculating the ratio of non-sinonymous to synonymous mutations normalized by length
mult_length_stat <-  multiclust_mism[, c(1,2,5)] %>%  group_by(domain, cluster) %>%  unique() %>% as.data.frame()
mult_num_stat <- multiclust_mism[, c(1,2,6)] %>%  group_by(domain, cluster) %>%  unique() %>% as.data.frame()
mult_types_stat <- multiclust_mism[, c(1,2,4)] %>%  group_by(domain, site_type, cluster) %>%  tally()  %>% as.data.frame()

#calculate the non-synonymous ratio
mult_non_syn_rate <- data.frame(cluster = c(), domain = c(), ratio = c())
for (domain in unique(mult_length_stat$domain)){ #Iterate over domains
  for (cluster in  unique(mult_length_stat$cluster)){ #Iterate over clusters
    #Get the number of toxins and sites
    clust_len <- mult_length_stat[mult_length_stat$cluster==cluster & mult_length_stat$domain==domain,3]
    clust_num <-  mult_num_stat[mult_num_stat$cluster==cluster & mult_num_stat$domain==domain,3]
    
    #Extract synonymous and non-synonymous sites
    non_syn_sites <- sum(mult_types_stat[mult_types_stat$cluster==cluster & mult_types_stat$domain==domain & mult_types_stat$site_type %in% c('non_syn','indel'),4])
    syn_sites <- sum(mult_types_stat[mult_types_stat$cluster==cluster & mult_types_stat$domain==domain & mult_types_stat$site_type %in% c('syn'),4])
    
    #Calculate the normalized ratio
    non_syn_ratio <- non_syn_sites/(syn_sites*clust_len*clust_num)
    
    #Update the summary dataframe
    non_syn_df <- data.frame(cluster = cluster, domain = domain, ratio = non_syn_ratio)
    mult_non_syn_rate<-rbind(mult_non_syn_rate, non_syn_df)
  }
}

#Plot the ditributions of non-synomimous rate
ggplot(mult_non_syn_rate, aes(x=domain, y=ratio, fill=domain)) +
  geom_boxplot(alpha=0.8)+scale_fill_manual(values = c('#a60b0b','#2980b9', '#bdbd00'))+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12))+
  xlab("Domain") +
  ylab('Nonsynonymous ratio')+
  guides(fill= guide_legend(title="Domain"))

#Properties of the toxins (conservation, length, identity, and species' attributions)
#Read distributions of conservative scores for all toxins and reference clusters
cons_cores <- read.csv('ref_cons_scores.csv', header = T, sep = "\t", stringsAsFactors = F)
cons_cores_all <- read.csv('all_cons_scores.csv', header = T, sep = "\t", stringsAsFactors = F)

cons_cores <- cons_cores[cons_cores$site_type=='gap_aln',]
cons_cores_all <- cons_cores_all[cons_cores_all$site_type=='gap_aln',]

#Mark domains according to the relative coordinates in the alignments
cons_cores$domain <- 'domain1'
cons_cores[cons_cores$coord<1663,4] <-'domain1'
cons_cores[cons_cores$coord<=3595 & cons_cores$coord>=1663,4] <-'domain2'
cons_cores[cons_cores$coord>3595,4] <-'domain3'

#Make the coordinates relative
cons_cores$coord <- cons_cores$coord/5111

#Mark domains according to the relative coordinates in the alignments
cons_cores_all$domain <- 'domain1'
cons_cores_all[cons_cores_all$coord<1715,4] <-'domain1'
cons_cores_all[cons_cores_all$coord<=3671 & cons_cores_all$coord>=1715,4] <-'domain2'
cons_cores_all[cons_cores_all$coord>3671,4] <-'domain3'

#Make the coordinates relative
cons_cores_all$coord <- cons_cores_all$coord/4818

#Assert the type of the dataset
cons_cores$type <- 'ref'
cons_cores_all$type <- 'all'

#Merge the data for reference clusters and all toxins
cons_cores_merged <- rbind(cons_cores,cons_cores_all)

#Calculate mean conservation score and the number of conservative sites
mean(cons_cores[cons_cores$domain=='domain1',3]) #0.23
mean(cons_cores[cons_cores$domain=='domain2',3]) #0.16
mean(cons_cores[cons_cores$domain=='domain3',3]) #0.16

length(cons_cores[cons_cores$domain=='domain1' & cons_cores$cons_score>0.85,3]) #73
length(cons_cores[cons_cores$domain=='domain2' & cons_cores$cons_score>0.85,3]) #16
length(cons_cores[cons_cores$domain=='domain3' & cons_cores$cons_score>0.85,3]) #41

#Plot coordinate-wise distribution of conservation scores
ggplot(cons_cores_merged, aes(x=coord, y=cons_score, fill=type, size=-log(cons_score))) +
  geom_point(alpha=0.4, shape=21)+
  scale_fill_manual(values=c("#DC143C", '#5F9EA0' ))+
  theme_bw()+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=16), 
         legend.text=element_text(size=14))+
  xlab('Coordinates') +
  ylab('Conservation score')+
  guides(fill= guide_legend(title="Domain"))+
  facet_wrap(~domain, scales = 'free_x')+ 
  ylim(c(0,1)) +
  scale_size(range = c(2, 5), guide = 'none')+
  guides(fill= guide_legend(title="Toxins' group", override.aes = list(size=5)))



#Read the data with pair-wise comparisions between the sequences of the domains for reference clusters and all toxins
id_domain_1 <- read.csv('id_ref_domain_1.csv', header = F, sep = "\t", stringsAsFactors = F)
id_domain_2 <- read.csv('id_ref_domain_2.csv', header = F, sep = "\t", stringsAsFactors = F)
id_domain_3 <- read.csv('id_ref_domain_3.csv', header = F, sep = "\t", stringsAsFactors = F)

id_domain_1_all <- read.csv('id_all_domain_1.csv', header = F, sep = "\t", stringsAsFactors = F)
id_domain_2_all <- read.csv('id_all_domain_2.csv', header = F, sep = "\t", stringsAsFactors = F)
id_domain_3_all <- read.csv('id_all_domain_3.csv', header = F, sep = "\t", stringsAsFactors = F)

#Assign the type of the toxins
id_domain_1$type <-'ref'
id_domain_2$type <-'ref'
id_domain_3$type <-'ref'

id_domain_1_all$type <-'all'
id_domain_2_all$type <-'all'
id_domain_3_all$type <-'all'

#Combine the data for reference clusters and all toxins
id_domain_1 <- rbind(id_domain_1, id_domain_1_all)
id_domain_2 <- rbind(id_domain_2, id_domain_2_all)
id_domain_3 <- rbind(id_domain_3, id_domain_3_all)

#Exclude comparisions for the same toxins
id_domain_1 <- id_domain_1[id_domain_1$V1!=0,]
id_domain_2 <- id_domain_2[id_domain_2$V1!=0,]
id_domain_3 <- id_domain_3[id_domain_3$V1!=0,]

#Calculate mean pair-wise identities
mean(id_domain_1[id_domain_1$type=='all', 1]) 
mean(id_domain_2[id_domain_2$type=='all', 1])
mean(id_domain_3[id_domain_3$type=='all', 1])

#Get minimal identity estimates
min(id_domain_1$V1)
min(id_domain_2$V1)
min(id_domain_3$V1)

#Plot the historgam with the distirubution of the pair-wise domain identities
ggplot(id_domain_3, aes(V1, fill=type)) +
  geom_histogram(color='black', alpha=0.9, bins = 120)+
  theme_bw()+scale_fill_manual(values = c("#DC143C", '#5F9EA0' ))+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12))+
  xlab('Identity percent') +
  ylab('Frequency')+xlim(c(30,100))+ylim(c(0, 60000))+
  guides(fill= guide_legend(title="Toxins' group"))


#Read the data with domain length fo reference clusters and all toxins
domain_lengths <- read.csv('lengths_combined.csv', header = F, sep = "\t", stringsAsFactors = F)
dom1 <-domain_lengths[domain_lengths$V3=='domain1',]
dom2 <-domain_lengths[domain_lengths$V3=='domain2',]
dom3 <-domain_lengths[domain_lengths$V3=='domain3',]

#Plot the distributions of domain lengths
ggplot(dom2, aes(x=V2, fill=V4)) +
  geom_histogram(color='black', alpha=0.9, bins = 100)+
  theme_bw()+scale_fill_manual(values = c("#DC143C", '#5F9EA0' ))+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12))+
  xlab("Domain length") +
  ylab('Frequency')+xlim(c(230, 850))+ylim(c(0,300))+
  guides(fill= guide_legend(title="Toxins' group"))


#Get the with Cry descriptive satistics (length, identity, species)
mearged_props <-read.csv('mearged_domain_props.csv', header = T, sep = "\t")

#Summarize species' distribution
species_df <- mearged_props[mearged_props$Domain=='domain1',c(2,8)] %>%  
  roup_by(Type, Species)  %>% group_by(Type, Species) %>% tally() 


#Plot the distribution of species
ggplot(species_df, aes(x=factor(1),y=n,fill=Species)) +
  geom_bar(stat="identity",alpha=0.55, col='black',
           position = "fill",width = 1)+ 
  facet_wrap(~Type)+
  coord_polar(theta = "y")+
  theme_bw()+
  scale_size(guide = 'none')+
  theme( axis.text.x = element_blank(),
         axis.title.x=element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12),
         axis.ticks = element_blank(),
         axis.line.x =  element_blank(),
         axis.line.y =  element_blank())+
  guides(fill= guide_legend(title="Species")) +
  scale_fill_manual(values = wes_palette("Darjeeling2", 12, type = "continuous"))


#Calculate mean identity and length
mean(mearged_props[mearged_props$Type=='all' & mearged_props$Domain=='domain1',5 ]) #53.18
mean(mearged_props[mearged_props$Type=='all' & mearged_props$Domain=='domain2',5 ]) #50.2
mean(mearged_props[mearged_props$Type=='all' & mearged_props$Domain=='domain3',5 ]) #51.83

mean(mearged_props[mearged_props$Type=='ref' & mearged_props$Domain=='domain1',5 ]) #51.44
mean(mearged_props[mearged_props$Type=='ref' & mearged_props$Domain=='domain2',5 ]) #49.23
mean(mearged_props[mearged_props$Type=='ref' & mearged_props$Domain=='domain3',5 ]) #50.64


mean(mearged_props[mearged_props$Type=='all' & mearged_props$Domain=='domain1',4]) #612.7367
mean(mearged_props[mearged_props$Type=='all' & mearged_props$Domain=='domain2',4 ]) #612.4297
mean(mearged_props[mearged_props$Type=='all' & mearged_props$Domain=='domain3',4 ]) # 419.8


#Calculate the number of novel proteins
length(mearged_props[mearged_props$Type=='all' & (mearged_props$Domain=='domain1' & mearged_props$NOM_group=='NOV') &
                       mearged_props$ID_NOV<=99,
                     4 ]) #304 -> 192
length(mearged_props[mearged_props$Type=='ref' & (mearged_props$Domain=='domain1' & mearged_props$NOM_group=='NOV'),
                     4 ]) #127

length(mearged_props[mearged_props$Type=='ref' & (mearged_props$Domain=='domain1' & mearged_props$NOM_group=='NOV') &
                       mearged_props$ID_NOV<=95,
                     4 ]) #304 -> 110 of 127

#Prepare the data for plotting domain identities and lengths
novel_dom_props <- mearged_props[mearged_props$NOM_group=='NOV',]

#Extract the data for novel toxins
novel_dom_props_single <- mearged_props[mearged_props$Domain=='domain1' & mearged_props$Type=='ref' & 
                                          mearged_props$NOM_group=='NOV' & mearged_props$Species !='Bacillus thuringiensis',]
#Extract the data for non-Bt species
dom_props_single_BT <- mearged_props[mearged_props$Domain=='domain1' & mearged_props$Type=='ref' & 
                                       mearged_props$Species!='Bacillus thuringiensis' & mearged_props$NOM_group=='BT',]


#Plot domain length vs mean domain-wise identity
ggplot(mearged_props, aes(ID_pair,Length,color=Species, shape=Domain)) +
  geom_point(alpha = 0.95,size=2.5) +
  theme_bw()+
  labs(y = "Domain length") + 
  labs(x = "Identity percent")+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12))+
  facet_wrap(~Type)+
  scale_color_manual(values = wes_palette("Darjeeling2", 12, type = "continuous"))+
  guides( shape=guide_legend(override.aes = list(size=5)), 
          color=guide_legend(override.aes = list(size=5)))


#Plot domain length vs identity with the closes homolog from BPPRC
ggplot(novel_dom_props, aes(ID_NOV,Length,color=Species, shape=Domain)) +
  geom_point(alpha = 0.95,size=2.5) +
  theme_bw()+
  labs(y = "Domain length") + 
  labs(x = "Identity percent")+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12))+
  facet_wrap(~Type)+
  scale_color_manual(values = wes_palette("Darjeeling2", 12, type = "continuous"))+
  guides( shape=guide_legend(override.aes = list(size=5)), 
          color=guide_legend(override.aes = list(size=5)))
