library(ggplot2)
library(ggsci)
library(dplyr)
library(zoo)
library(micropan)
library(ggsignif)
library(gridExtra)
library(ggfortify)
library(factoextra)
library(lsa)

setwd("~/phyl/data_for_git/data_for_scripts/for_viz/mechanisms")

#Vizualizing the results of pangenome reconstruction
#Read the resuts obtained from the Panaroo software
pangenome_presence<- read.csv('rec_presence.csv', header = T, sep = ",", stringsAsFactors = F) #cry_presence, all_presence

#Create a pangenome preence/absence matrix in the binary format
pangenome_presence_matr <- pangenome_presence[,-c(1:3)] %>%  as.data.frame()
pangenome_presence_matr[pangenome_presence_matr!='' ]<-1
pangenome_presence_matr[pangenome_presence_matr=='' ]<-0

#Transposing the matrix
pan_matrix <- matrix(as.numeric(unlist(pangenome_presence_matr)),nrow = nrow(pangenome_presence_matr),ncol = 65)
t_pan <- t(pan_matrix)

#Calculare the alpha parameter to assess the openness of the nagenome using the Heaps' law
heaps_res_pan <- heaps(t_pan, n.perm = 50) 

#Results of calculations
#0.689, 0.737, 0.57 - all, cry, rec

#Create a dataframe with the number of pangenome clusters for the reconstructed pangenomes
Gene_core_pangenome=data.frame(genes=c(2872, 2402, 2146,
                                       965, 1374, 1679,
                                       3102, 4147, 4038,
                                       33026, 21488, 13196),
                               gene_type =c('Core','Core','Core',
                                            'Soft core','Soft core','Soft core',
                                            'Shell','Shell', 'Shell',
                                            'Cloud','Cloud','Cloud'),
                               pangenome=c('All','Cry','Rec',
                                           'All','Cry','Rec',
                                           'All','Cry','Rec',
                                           'All','Cry','Rec'))
#Calculate the percentage of core genes
2872/39965 #0.71%
2402/29411 #0.81%
2146/21059 #0.10%

#Assert the order of factos for correct vizualization
Gene_core_pangenome$gene_type=factor(Gene_core_pangenome$gene_type, levels = c('Core','Soft core', 'Shell', 'Cloud'))

#Viualize a pie chart with gene distributions
ggplot(Gene_core_pangenome[Gene_core_pangenome$pangenome=='Rec',], aes(x='',y=genes,fill=gene_type)) +
  geom_bar(stat="identity", width=1, alpha=0.75, col='black')+
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
  guides(fill= guide_legend(title="Gene cluster group"))+
  scale_fill_npg()

#Analysis of mobile genetic elements (MGEs)
#Read the table with the abundance of mobile genetic elements per assambly
mges_num_non_red <- read.csv('MGES_num_summary_non_redundant.csv', header = T, sep = "\t")

#Reshape the dataframe for vizualization
mges_num_reshaped <- reshape(mges_num_non_red, idvar = c("Accession","Dataset"), timevar = "MGE_type", direction = "wide")

#Plot the dependancy between the number of distinct MGEs with linear models applied
ggplot(mges_num_reshaped, aes(x=MGE_num.IS,y=MGE_num.phrophages, fill=Dataset)) +
  geom_point(shape=21, color='black', size=4, alpha=0.65)+
  geom_smooth(method='lm', color='darkgrey', alpha=0.45, fill = 'grey')+
  scale_size(guide = 'none')+
  scale_fill_manual(values = c('#0D38A2','#088530', '#BA0D04'))+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=20),
        legend.title=element_text(face="bold",size=18), 
        legend.text=element_text(size=16),
        axis.text.x = element_text(color='black',  size = 16))+
  guides(fill= guide_legend(title="Dataset",
                            override.aes = list(size=6)))+
  ylab('Number of prophages')+ #'Number of genetic islands' 'Number of insertions' 'Number of prophages'
  xlab("Number of insertions")

#Calculate the correlation between the number of diverse MGEs in genomes
cor.test(mges_num_reshaped$MGE_num.phrophages, 
         mges_num_reshaped$MGE_num.IS, method=c("pearson" )) #phages/IS  #0.4841254  corr 2.2e-16 p-value
cor.test(mges_num_reshaped$MGE_num.phrophages, 
         mges_num_reshaped$MGE_num.GIs, method=c("pearson" )) #phages/GI #0.4778635  corr 2.2e-16 p-value
cor.test(mges_num_reshaped$MGE_num.IS, 
         mges_num_reshaped$MGE_num.GIs, method=c("pearson" )) #IS/GI #0.3903698  corr 9.098e-16 p-value

#Calculate the difference between datasets containing different types of toxins in terms of the number of MGEs
anova <- aov(MGE_num ~ Dataset, data = mges_num_non_red)
tukey <- TukeyHSD(anova) # 1e-08

#Summarize median numbers of MGEs per dataset
median_stat <- mges_num_non_red %>% group_by(Dataset,MGE_type) %>% summarize(median_mge=median(MGE_num)) %>% as.data.frame()

#Vizualize the distribution of MGEs abundance per dataset
ggplot(mges_num_non_red, aes(x=Dataset, y=MGE_num, fill=Dataset)) + 
  facet_wrap(~MGE_type, scales = 'free')+
  geom_boxplot(alpha=0.65)+
  scale_fill_manual(values = c('#0D38A2','#088530', '#BA0D04'))+
  geom_signif(comparisons = list(c(1, 2),
                                 c(2, 3),
                                 c(1,3)),
              map_signif_level=T, test='wilcox.test')+
  theme_bw()+xlab('Dataset')+ylab('Number of MGEs')+
  theme( axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=18),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=14),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=18),
         legend.title=element_text(face="bold",size=16), 
         legend.text=element_text(size=14))

#Create a dataframe with the number of cry genes associated with genetic islands
GI_dat <- data.frame(
  "GI_non" = c(430, 65),
  "GI_yes" = c(172, 38),
  row.names = c("Cry", "Rec"),
  stringsAsFactors = FALSE
)

#Compare the abundance of asociations using the Fisher's test
GI_test <- fisher.test(GI_dat)
GI_test #p-value = 0.1023

#Create a dataframe with the number of MGEs associated with cry genes
cry_mge_def = data.frame(num = c(172, 15, 415, 38, 0, 65),
                         toxin_group = c('cry', 'cry', 'cry', 'rec', 'rec','rec'),
                         mge_type = c('GIs','IS','None', 'GIs','IS','None'))

#Vizualize the proportion of MGEs associated with cry genes
ggplot(cry_mge_def, aes(x=toxin_group,y=num,fill=mge_type)) +
  facet_wrap(~toxin_group, scales = 'free')+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values =c('#0DA5A5', '#E3AD7B', 'grey'))+
  theme_bw()+xlab('Dataset')+ylab('Number of toxins')+
  theme( axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=18),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=14),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=18),
         legend.title=element_text(face="bold",size=16), 
         legend.text=element_text(size=14))

#Analysis of the identity between parental sequences within the regions flankning breakpoints
#Read the table with coordinate-wie identity between parent in recombination events
flank_seqs_results_pars <- read.csv('Parents_only_Distribution_of_flank_identities_merged.csv', header = T, sep = "\t", stringsAsFactors = F)

#The function rom calculating standart mean error
stderror <- function(x) sd(x)/sqrt(length(x))

#Extract flanks between the first and the second and between the second and the third domains
flank_seqs_results_all <- flank_seqs_results_pars[flank_seqs_results_pars$Flank_class %in% c('Zone2', 'Zone3'),]

#Calcuate mean and mean square errors per each coordinate
mean_identity_per_coord_all <- flank_seqs_results_all %>% group_by(Coord,Flank_class) %>% 
  summarize(Ident = mean(Identity), stderror=stderror(Identity)) %>% as.data.frame() 

#Apply rollmean function to coordinate wise distributions of mean identity and standart mean error
mean_identity_per_coord_all[mean_identity_per_coord_all$Flank_class=='Zone2',3] <- rollapplyr(mean_identity_per_coord_all[mean_identity_per_coord_all$Flank_class=='Zone2',3],11, mean, partial = TRUE, align='center')
mean_identity_per_coord_all[mean_identity_per_coord_all$Flank_class=='Zone3',3] <-rollapplyr(mean_identity_per_coord_all[mean_identity_per_coord_all$Flank_class=='Zone3',3],11, mean, partial = TRUE, align='center')

mean_identity_per_coord_all[mean_identity_per_coord_all$Flank_class=='Zone2',4] <- rollapplyr(mean_identity_per_coord_all[mean_identity_per_coord_all$Flank_class=='Zone2',4],11, mean, partial = TRUE, align='center')
mean_identity_per_coord_all[mean_identity_per_coord_all$Flank_class=='Zone3',4] <-rollapplyr(mean_identity_per_coord_all[mean_identity_per_coord_all$Flank_class=='Zone3',4],11, mean, partial = TRUE, align='center')

#Vizualize roll mean identity values surronding selected breakpoints
ggplot(mean_identity_per_coord_all[mean_identity_per_coord_all$Flank_class=='Zone2',], aes(y=Ident, x=Coord))+
  geom_line(aes(colour=Flank_class))+
  geom_vline(xintercept = 100, linetype="dashed", color = "black",size=1.1, alpha = 0.7)+
  geom_ribbon(aes(ymin=Ident-stderror, ymax=Ident+stderror), alpha=0.2) +
  scale_colour_manual(values = c( '#DC143C' ))+
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
  xlab("Coordinate") +
  ylab('Identity')+
  guides(color= guide_legend(title="Region type"))+
  facet_wrap(~Flank_class)+ylim(c(0.35, 0.99))

#Perform clusterization of flanks per recombination event based on the distributions of coordinate-wie identitu estimates
#Functions for calculating distance between two vectors (both on mean and rollmean estimates)
#Eucledean mean
euclidean_mean <- function(a, b) sqrt(sum((mean(a) - mean(b))^2))
euclidean_vec <- function(a, b)  sqrt(sum((a - b)^2))
euclidean_vec_roll <- function(a, b) sqrt(sum((rollapplyr(a,11, mean, partial = TRUE, align='center') - rollapplyr(b,11, mean, partial = TRUE, align='center'))^2))

#Manhattan distance
vec_sub <-  function(a, b) sum(abs(a-b))
vec_sub_roll <-  function(a, b) sum(abs(rollapplyr(a,11, mean, partial = TRUE, align='center')-rollapplyr(b,11, mean, partial = TRUE, align='center')))

#Cosine distance
vec_cos <-  function(a, b) 1-round(cosine(a,b),6)
vec_cos_roll <-  function(a, b) 1-round(cosine(rollapplyr(a,11, mean, partial = TRUE, align='center'),rollapplyr(b,11, mean, partial = TRUE, align='center')),6)

#The function for contructing a distance matrix based on a data of coordinate-wie identity per recombination event
make_dist_matr<-function(id_data_frame, dist_func){
  #Get the list of recombination events
  events_list <- unique(id_data_frame$ID)
  
  #Create a dummy distance matrix
  dist_matrix <- matrix(data=0, nrow=length(events_list), ncol = length(events_list), 
                        dimnames = list(as.character(events_list),as.character(events_list)))
  
  #Apply distance function to calculate the difference between the events
  id <- dist_func(id_data_frame[id_data_frame$ID==3,3], id_data_frame[id_data_frame$ID==0,3])
  
  #Iterate over recombination events with respective calculating of identity difference
  for (ID1 in events_list){
    for (ID2 in events_list){
      vec1<-id_data_frame[id_data_frame$ID==ID1, 3]
      vec2<-id_data_frame[id_data_frame$ID==ID2, 3]
      dist_matrix[as.character(ID1), as.character(ID2)]<- dist_func(vec1, vec2)
    }
  }
  return(dist_matrix)
}

#The function for adding the threshold for clusters on the dendrorgram
add_lines_to_clust <- function(k, n, clusters){
  
  #Calculate the coordinates of the line based on branch lengths and selected number of clusters
  MidPoint = (clusters$height[n-k] + clusters$height[n-k+1]) / 2
  abline(h = MidPoint, lty=2)
}

#The function for vizualizing results of hierarcheal clusters
draw_hclust_res <- function(clusters, k_vec, n, out_name){
  #Plot the dendrogram
  plot(clusters)
  
  #Add lines with clustering thresholds
  for (k in k_vec){
    add_lines_to_clust(k,n, clusters)
  }
  
  #Save plots as pdf files
  dev.print(pdf, paste0(out_name,'.pdf'))
  dev.off()
}

#The function for inferring the clustering patterns for each recombination event
get_hc_clusters<- function(clusters,clust_method, cut_num){
  
  #Get cluster cut
  clusterCut <- cutree(clusters, cut_num)
  
  #Create a dataframe stating the method of clusterization, total number of clusters and event-based clustering attributions
  names_frame <- data.frame(ID=names(clusterCut), Cluster=clusterCut, 
                            Cluster_method=clust_method, Cluster_num=cut_num, row.names=NULL) 
  return(names_frame)
}

#The general function for performing hiearcheal clusterization with plotting the results and inferring event-wise attributions
analyze_hc_clusters <- function(dist_matr, k_vec, name_pref){
  
  #Perform clustering on the distance matrix
  clusters <- hclust(as.dist(as.matrix(dist_matr)), method = "complete")
  n = nrow(dist_matr)
  name = paste0(name_pref,'_hclust')

  #Plot dendrograms of clustering results
  draw_hclust_res(clusters, k_vec,n, name)
  out_df <- data.frame(ID=c(), Cluster=c(), Cluster_method=c(), Cluster_num=c(), row.names=NULL) 

  #Create dataframes with event-wise clustering patterns 
  for (k in k_vec){
    clust_df <- get_hc_clusters(clusters,name,k)
    out_df <- rbind(out_df,clust_df)
  }
  return(out_df)
}

#The function for vizualizing k-means clustering results via the autoplot function
draw_autoplot <- function(clust_kmeans,sc_kmeans){
  p<-autoplot(clust_kmeans, 
              data=sc_kmeans, label = TRUE,frame = TRUE,
              frame.type = 'convex',label.size = 3.3,label.col='black', alpha=0.8)+
    scale_colour_manual(values = c('grey','grey','grey', 'grey','grey', 'grey', 'grey'))+
    scale_fill_lancet()+
    theme_bw()+
    theme(axis.text.y = element_text(color='black', 
                                     size=14),
          axis.title.y=element_text(color="black", 
                                    size=16),
          panel.background = element_blank(), 
          axis.line = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(color='black', 
                                     size=14),
          axis.title.x = element_text(color="black", 
                                      size=16),
          legend.position = 'none'
    ) + 
    scale_x_continuous(expand = c(.1, .1)) 
  scale_y_continuous(expand = c(.1, .1))
  return(p)
}

#the function for performing k-means clusterinf on the distance matrix
get_kmeans_clusters<- function(dist_matr,clust_method, cut_num){
  sc_kmeans <- data.frame(t(dist_matr))
  
  #Perform clustering with a defined number of clusters
  clust_kmeans <- kmeans(sc_kmeans,centers=cut_num)
  
  #Obtain the data with clusters for recombination events
  sc_kmeans$ID = rownames(sc_kmeans)
  sc_kmeans$clust = as.factor(clust_kmeans$cluster)
  names_frame <- data.frame(ID=sc_kmeans$ID, Cluster=sc_kmeans$clust, Cluster_method=clust_method, Cluster_num=cut_num, row.names=NULL) 
  return(names_frame)
}

#The function for vizualizing k-means clutering results
draw_kmeans_clusters <- function(dist_matr,cut_num){
  
  #Perform k-means clustering procedure
  sc_kmeans <- data.frame(t(dist_matr))
  clust_kmeans <- kmeans(sc_kmeans,centers=cut_num)
  
  #Plot results via the autoplot function
  k_plot <- draw_autoplot(clust_kmeans,sc_kmeans)
  return(k_plot)
}

#The overall function for k-mean clutering
analyze_kmeans_clusters <- function(dist_matr, k_vec, name_pref){
  
  #Create the dataframe with event-wise clustering attributions
  name = paste0(name_pref,'_kmeans')
  out_df <- data.frame(ID=c(), Cluster=c(), Cluster_method=c(), Cluster_num=c(), row.names=NULL) 
  
  #Iterate over the list of clustering thresholds
  k_plots <- list()
  for (k in k_vec){
    #Perform the k-means procedure
    clust_df <- get_kmeans_clusters(dist_matr, name, k)
    out_df <- rbind(out_df,clust_df)
    
    #Save the results of the autoplot functionss
    k_plot <- draw_kmeans_clusters(dist_matr,k)
    k_plots <- c(k_plots, list(k_plot))
  }
  
  #Arrange autoplots into the panel
  do.call(grid.arrange, c(k_plots, nrow = 1))
  dev.print(pdf, paste0(name,'.pdf'), width = 12, height = 5)
  dev.off()
  return(out_df)
}

#The function for defining the optimal number of clusters based on the silhouette analysis
define_optimal_clusters <- function(dist_matr, k_num, name_pref){
  #Make the list of the clustering qualities
  stat_plots <- list()
  
  #Iteratively applu the fviz_nbclust function with plotting clustering quality
  for (clust_method in c(hcut, kmeans)){
    for (stat_type in c('silhouette', 'wss', 'gap_stat')){
      stat_plot <- fviz_nbclust(dist_matr, clust_method, method = stat_type, k.max=k_num)
      stat_plots <- c(stat_plots, list(stat_plot))
    }
  }
  
  #Arrange the resulting plots into panels
  do.call(grid.arrange, c(stat_plots, nrow = 2))
  dev.print(pdf, paste0(name_pref,'.pdf'), width = 15, height = 8)
  dev.off()
}

#The function for iterating ove two objects 
zip <- function(...) {
  mapply(list, ..., SIMPLIFY = FALSE)
}

#The general function for assesing clusterization quality and plotting the results obtained
get_stat_for_clusts <- function(dist_matr, matr_pref, k_num, k_vec){
  #Crate a dataframe with the assessments per method and disnace metrics
  stat_res <- data.frame(clusters=c(),y=c(),matr=c(), clust_type=c())
  
  #The list with clustering methods
  meth_list <- zip(c(hcut, kmeans),c('hclust','kmeans'))
  
  #Iterate over methods (hclust and k-means) for assesing quality
  for (ind in 1:2){
    clust_method <- meth_list[[ind]][[1]]
    method_name <- meth_list[[ind]][[2]]
    
    #Apply the fviz_nbclust function
    shil_stat <- fviz_nbclust(dist_matr, clust_method, method = 'silhouette', k.max=k_num)
    
    #Save the quality estimates for a defined distance metrics 
    shil_stat_df <- shil_stat$data[k_vec, ]
    shil_stat_df$matr <- matr_pref
    shil_stat_df$clust_type <- method_name
    stat_res <- rbind(stat_res, shil_stat_df)
  }  
  return(stat_res)
}

#Extract idetity distributions for different flanks
flanks_extr_zone_2 <- flank_seqs_results_pars[flank_seqs_results_pars$Flank_class %in% c('Zone2'),]
flanks_extr_zone_3 <- flank_seqs_results_pars[flank_seqs_results_pars$Flank_class %in% c('Zone3'),]

#Create the lists of distribution matricies per distance metics
fl3_matr_list <- list(fl3_cos=make_dist_matr(flanks_extr_zone_3, vec_cos),
                      fl3_cos_roll=make_dist_matr(flanks_extr_zone_3,vec_cos_roll),
                      fl3_eucl_mean=make_dist_matr(flanks_extr_zone_3,euclidean_mean),
                      fl3_eucl_vec=make_dist_matr(flanks_extr_zone_3,euclidean_vec),
                      fl3_eucl_vec_roll=make_dist_matr(flanks_extr_zone_3,euclidean_vec_roll),
                      fl3_sub_vec=make_dist_matr(flanks_extr_zone_3,vec_sub),
                      fl3_sub_vec_roll=make_dist_matr(flanks_extr_zone_3,vec_sub_roll))

fl2_matr_list <- list(fl2_cos=make_dist_matr(flanks_extr_zone_2, vec_cos),
                      fl2_cos_roll=make_dist_matr(flanks_extr_zone_2,vec_cos_roll),
                      fl2_eucl_mean=make_dist_matr(flanks_extr_zone_2,euclidean_mean),
                      fl2_eucl_vec=make_dist_matr(flanks_extr_zone_2,euclidean_vec),
                      fl2_eucl_vec_roll=make_dist_matr(flanks_extr_zone_2,euclidean_vec_roll),
                      fl2_sub_vec=make_dist_matr(flanks_extr_zone_2,vec_sub),
                      fl2_sub_vec_roll=make_dist_matr(flanks_extr_zone_2,vec_sub_roll))

#Create emtpy dataframe templates with clustering patterns per event
clust_res_zone3 <- data.frame(ID=c(), Cluster=c(), Cluster_method=c(), Cluster_num=c())
clust_res_zone2 <- data.frame(ID=c(), Cluster=c(), Cluster_method=c(), Cluster_num=c())

clust_optimal_num_zone3 <- data.frame(ID=c(), Cluster=c(), Cluster_method=c(), Cluster_num=c())
clust_optimal_num_zone2 <- data.frame(ID=c(), Cluster=c(), Cluster_method=c(), Cluster_num=c())

#Iterate over distance matricies to assess the quality of clutering and vizualize results
for (matr_pref in names(fl3_matr_list)){
  #Scale 
  dist_matr = scale(fl3_matr_list[[matr_pref]])
  define_optimal_clusters(dist_matr, 10, matr_pref)
  hc_res <- analyze_hc_clusters(dist_matr, c(5,6),matr_pref)
  kmeans_res <- analyze_kmeans_clusters(dist_matr, c(5,6),matr_pref)
  clust_res_zone3 <- rbind(clust_res_zone3, hc_res)
  clust_res_zone3 <- rbind(clust_res_zone3, kmeans_res)
  shil_res <- get_stat_for_clusts(dist_matr, matr_pref, 10, c(4,5,6))
  clust_optimal_num_zone3 <- rbind(clust_optimal_num_zone3, shil_res)
}

for (matr_pref in names(fl2_matr_list)){
  dist_matr = scale(fl2_matr_list[[matr_pref]])
  define_optimal_clusters(dist_matr, 7, matr_pref)
  hc_res <- analyze_hc_clusters(dist_matr, c(2,3,4),matr_pref)
  kmeans_res <- analyze_kmeans_clusters(dist_matr, c(2,3,4),matr_pref)
  clust_res_zone2 <- rbind(clust_res_zone2, hc_res)
  clust_res_zone2 <- rbind(clust_res_zone2, kmeans_res)
  
  shil_res <- get_stat_for_clusts(dist_matr, matr_pref, 7, c(2,3,4))
  clust_optimal_num_zone2 <- rbind(clust_optimal_num_zone2, shil_res)
}

#Exclude results based on the total mean estimate
clust_optimal_num_zone3[clust_optimal_num_zone3$matr!='fl3_eucl_mean',]
clust_optimal_num_zone2[clust_optimal_num_zone3$matr!='fl2_eucl_mean',]

#Get the most optimal number of clusters
optimal_num_max_3 <- clust_optimal_num_zone3 %>% group_by(matr,clust_type) %>% top_n(1,y) %>% as.data.frame()
optimal_num_max_2 <- clust_optimal_num_zone2 %>% group_by(matr,clust_type) %>% top_n(1,y) %>% as.data.frame()

#Summarize the optimal number of cluters per distance metrics and flank type
clust_shil_df <- rbind(clust_optimal_num_zone3[clust_optimal_num_zone3$matr!='fl3_eucl_mean',],
                       clust_optimal_num_zone2[clust_optimal_num_zone3$matr!='fl2_eucl_mean',])
write.table(clust_shil_df ,file = "silhouette.csv",sep='\t', row.names = F)

#Vizualize the optimal number of clusters per clustering method and distance metrics
ggplot(optimal_num_max_2[optimal_num_max_2$matr!='fl2_eucl_mean',], 
       aes(y=y, x=clusters, color = matr, shape = clust_type)) +
  geom_point(size =5, alpha = 0.7)+
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
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12))+
  xlab("Number of clusters")+
  ylab('Maximal silhouette')+scale_color_lancet()+
  guides(color= guide_legend(title="Distance metrics"),
         shape = guide_legend(title="Clustering method"))

#The function for obtaining the distribution of coordinate-wise identity with rollmean function applied
make_roll_ident_by_cluster <- function(clust_res,flanks_extr){
  #Make a template emtpy dataframe
  identity_by_cluster <- data.frame(Coord = c(), Ident = c(), St_err = c(), Cluster=c(), Cluster_method=c(), Cluster_num=c())
  
  #Iterate over clustering method and optimal number of clusters
  for (Cluster_method in unique(clust_res$Cluster_method)){
    for (Cluster_num in unique(clust_res$Cluster_num)){
      
      #Get the data for selected methods
      Clusters <- unique(clust_res[clust_res$Cluster_method==Cluster_method & clust_res$Cluster_num==Cluster_num, 2])
      for (Cluster in Clusters){
        
        #Exctract the coordinate-wise distribution of identity
        ID_vec = as.numeric(as.character(clust_res[clust_res$Cluster_method==Cluster_method & clust_res$Cluster_num==Cluster_num & clust_res$Cluster==Cluster, 1]))
        mean_identity_extr <- flanks_extr[flanks_extr$ID %in% ID_vec, ] %>% group_by(Coord) %>% 
          summarize(Ident = mean(Identity), St_err=stderror(Identity)) %>% as.data.frame() 
        
        mean_identity_extr$Cluster <- Cluster
        mean_identity_extr$Cluster_method <- Cluster_method
        mean_identity_extr$Cluster_num <- Cluster_num
        
        #Apply the rollmean function on identity distibution
        mean_identity_extr$Ident <- rollapplyr(mean_identity_extr$Ident, 11, mean, partial = TRUE, align='center')
        mean_identity_extr$St_err <- rollapplyr(mean_identity_extr$St_err, 11, mean, partial = TRUE, align='center')
        identity_by_cluster <- rbind(identity_by_cluster, mean_identity_extr)
      }
    }
  }
  return(identity_by_cluster)
}


#The function for plotting coordinate-wise identity distribution per each cluster of the events
draw_ident_distr <- function(identity_by_cluster){
  p <- ggplot(identity_by_cluster, 
              aes(y=Ident, x=Coord, colour= Cluster)) + 
    geom_line()+
    scale_colour_identity()+
    scale_colour_lancet()+
    geom_ribbon(aes(ymin=Ident-St_err, ymax=Ident+St_err), alpha=0.2,colour = NA) +
    geom_vline(xintercept = 100, linetype="dashed", color = "black",size=1.1, alpha = 0.7)+
    facet_wrap(~Cluster_num+Cluster)+
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
    xlab("Region type") +
    ylab('Coordinate')+
    guides(color= guide_legend(title="Cluster"))
  return(p)
}

#The function por plotting the distiribution of coordinate-wise identity for different c;ustering methods
make_pics_for_clust_dist <- function(ident_clusters){
  
  #Iterate over clustering methods
  for (Cluster_method in unique(ident_clusters$Cluster_method)){
    print(Cluster_method)
    
    #Draw identity distributions and save the figures
    p <- draw_ident_distr(ident_clusters[ident_clusters$Cluster_method==Cluster_method, ]) 
    print(p)
    dev.print(pdf, paste0(Cluster_method,'.pdf'), width = 14, height = 8)
    dev.off()
  }
}
 
#Get cluster-wise distribution of identity with the rollmean function applied
identity_by_cluster_zone3 <- make_roll_ident_by_cluster(clust_res_zone3,flanks_extr_zone_3)
identity_by_cluster_zone2 <- make_roll_ident_by_cluster(clust_res_zone2,flanks_extr_zone_2)

#Draw the distributions per cluster for different algorithms and distance metrics
make_pics_for_clust_dist(identity_by_cluster_zone3)
make_pics_for_clust_dist(identity_by_cluster_zone2)

#Read the data for event-wie cluster attributions
clust_res_saved <- read.csv('clusters_assignments_final.csv', header = T, sep = "\t", stringsAsFactors = F)

#Extract data for selected flanking sequences
flanks_extr_zone_2 <- flank_seqs_results_pars[flank_seqs_results_pars$Flank_class %in% c('Zone2'),]
flanks_extr_zone_3 <- flank_seqs_results_pars[flank_seqs_results_pars$Flank_class %in% c('Zone3'),]

#Peaks coordinates (assesed manually):
#flank2: 1) 95-107, 125-135 2) 145-158, 90-100
#flank3: 1) 70-100, 125-150 2) 95-125 3) 85-110, 145-160 4) 85-110, 145-165 5) 35-60, 110-140 6) 160-180

#Create the dataframe with coordinates of the peaks per flank and cluster
flanks_peaks_df <- data.frame(Flank_class = c('Zone2','Zone2','Zone2','Zone2',
                                              'Zone3', 'Zone3', 'Zone3', 'Zone3', 'Zone3', 
                                              'Zone3', 'Zone3', 'Zone3', 'Zone3', 'Zone3'),
                              Cluster = c(1,1,2,2,1,1,2,3,3,4,4,5,5,6),
                              Start = c(95,125,145,90,70,125,95,85,145,85,145,35,110,160),
                              Stop = c(107,135,158,100,100,150,125,110,160,110,165,60,140,180))

#The function for extracting the data for selected identity peaks
extract_peaks <- function(ID, Start, Stop, flank_dist) {
  extr <- flank_dist[flank_dist$ID==ID & flank_dist$Coord>=Start & flank_dist$Coord<=Stop,]
  return(extr)
}

#The function for obtaining mean identity estimates within selected cooordinates of peaks
get_mean_id_for_peaks <- function (flank_dist,clust_res, flanks_peaks_df) {
  
  #Make a template dataframe with identity estimates
  ret_df <- data.frame(ID=c(), Coord=c(), Identity= c(), Flank_num=c(), 
                       Reg_type=c(),Reg_domain=c(), Flank_class=c(), plot_factor=c())
  
  #Iterate over cluters
  for (clust in unique(clust_res$Cluster)){
    IDs <- clust_res[clust_res$Cluster==clust,1]
    peaks_coords <- flanks_peaks_df[flanks_peaks_df$Cluster==clust,]
    
    #Iterate ove recombination events
    for (ID in IDs){
      #Iterate over coordinates
      for (coord_ind in 1:nrow(peaks_coords)){
        Start = peaks_coords[coord_ind,3]
        Stop = peaks_coords[coord_ind,4]
        
        #exctact the identity distribution for a defined peak
        extr_res <- extract_peaks(ID, Start, Stop, flank_dist)
        ret_df <- rbind(ret_df, extr_res)
      }
    }
  }
  return(ret_df)
}

#Get mean identity estimates for different flanks within the best clustering method and distance metrics
flank3_peaks <- get_mean_id_for_peaks(flanks_extr_zone_3, clust_res_saved[clust_res_saved$Cluster_method=='fl3_cos_roll_hclust',], 
                                      flanks_peaks_df[flanks_peaks_df$Flank_class=='Zone3',])

flank2_peaks <- get_mean_id_for_peaks(flanks_extr_zone_2, clust_res_saved[clust_res_saved$Cluster_method=='fl2_cos_roll_hclust',], 
                                      flanks_peaks_df[flanks_peaks_df$Flank_class=='Zone2',])

#Create the mearged dataframe summarizing results
flank_peaks_mearged <- rbind(flank3_peaks, flank2_peaks)
peaks_stat_per_event <- flank_peaks_mearged %>% group_by(ID, Flank_class) %>% 
  dplyr::summarize(Ident = mean(Identity)) %>% as.data.frame()

#Read the data with mean identity for transferred and adjacent domains for further comparision with selected peaks
flank_for_model <- read.csv('flanks_props_for_model.csv', header = T, sep = "\t", stringsAsFactors = F)

#Exctract selected data and calculate mean estimates for both domains
ID_testing_df <- flank_for_model[,c(1,2,3,4,5,9,10,11)]
ID_testing_df$merged_flank_ID <- rowMeans(ID_testing_df[,c(4,5)])*100
ID_testing_df$merged_dom <- rowMeans(ID_testing_df[,c(7,8)])
ID_testing_df <- ID_testing_df[,c(1:3, 7:10)]
ID_testing_df$break_100 <- flank_for_model$break_flank*100

#Add the mean estimates of the peaks for selected events
ID_testing_df$peak <- peaks_stat_per_event[match(ID_testing_df$ID,peaks_stat_per_event$ID),3]*100

#Create a dataframe with identity estimates for further vizualization
reshaped_id_df <- data.frame(ID=c(), Domain=c(), Zone=c(), Ident=c(), Type=c())
for (col_pref in c('trans_ident', 'adj_ident', 'merged_flank_ID','merged_dom','break_100', 'peak')){
  ind<-which(colnames(ID_testing_df)==col_pref)
  sub_id_df <- data.frame(ID=ID_testing_df$ID, Domain=ID_testing_df$Domain, Zone=ID_testing_df$Zone, Ident=ID_testing_df[,ind])  
  sub_id_df$Type <-col_pref
  reshaped_id_df <- rbind(reshaped_id_df, sub_id_df)

}

#Assert proper names of region types for vizualization
reshaped_id_df <- reshaped_id_df[reshaped_id_df$Type %in% c('peak','trans_ident', 'adj_ident'),]
reshaped_id_df[reshaped_id_df$Type=="trans_ident",5] <- 'Transferred domain'
reshaped_id_df[reshaped_id_df$Type=="adj_ident",5] <- 'Adjacent domain'
reshaped_id_df[reshaped_id_df$Type=="peak",5] <- 'Peaks'

#Compare mean estimated for peaks and whole domains
wilcox.test(reshaped_id_df[c(1:10),4], reshaped_id_df[c(75:84),4])

#Plot the difference between mean identity estimates for peaks and whole domains
reshaped_id_df %>%  group_by(Zone, Type) %>% summarise(mean_id=mean(Ident))
ggplot(reshaped_id_df, aes(y=Ident, x=Type, fill=Type)) + #Num_min  Num_maj trans_ident adj_ident trans_gca trans_gca_full  adj_gca adj_gca_full
  geom_boxplot()+
  facet_wrap(~Zone)+
  geom_signif(comparisons = list(c(2, 3),
                                 c(1,2)),
              map_signif_level=T, test='wilcox.test')+
  theme_bw()+
  theme( axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
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
  xlab("Region type") +
  ylab('Mean Idetity')+
  guides(fill= guide_legend(title="Region type"))
