library(ggplot2)
library(ggsci)
library(dplyr)
library(vegan)

setwd("../data_for_scripts/data_for_scripts/for_viz/evol_selection")

#Analyze branch models
#Read the resuls of branch models calculations
branch_models_classified <- read.csv('Branch_models_extended.csv', header = T, sep = "\t", stringsAsFactors = F)

#Extract significant results of branch model
branch_models_classified_signif <- branch_models_classified[branch_models_classified$signif_flag=='signif' & branch_models_classified$evol_type!='non-signif',] 

#Create dummy factors for plotting for the X and Y axis, respectively
branch_models_classified_signif$fac_for_heatmap <- paste(branch_models_classified_signif$Tox_type, 
                                                         branch_models_classified_signif$Curr_domain, sep ='_')

branch_models_classified_signif$domain_factor <- as.factor(paste(branch_models_classified_signif$Rec_domain,
                                                                 branch_models_classified_signif$ID,
                                                                 sep='_'))

#Get the list of events with unknown parents
unknown_majors <- unique(branch_models_classified[branch_models_classified$unknown_flag=='major',1])
unknown_minors <- unique(branch_models_classified[branch_models_classified$unknown_flag=='minor',1])

#Assert types of unknown parents
branch_models_classified_signif[branch_models_classified_signif$ID %in% unknown_majors, 5] <- 'major'
branch_models_classified_signif[branch_models_classified_signif$ID %in% unknown_minors, 5] <- 'minor'

#Create the list of events with different selection pattenrs for parents and recombinantss
change_pattern <- c(98, 74, 32, 103, 71, 22, 96, 87, 75, 56, 47, 45, 95)

#Assert the type of selection pattern
branch_models_classified_signif$change_flag <- 'same'
branch_models_classified_signif[branch_models_classified_signif$ID %in% change_pattern,26] <- 'different'

#Extract recombination IDs for recombinants only
events_with_recs <- branch_models_classified_signif[branch_models_classified_signif$Tox_type=='Rec' &
                                                      branch_models_classified_signif$evol_type!='non-signif',1]

#Extract the data for recombinants only
heatmap_branch_recs <- branch_models_classified_signif[branch_models_classified_signif$ID %in% events_with_recs,]

#Assert the levels of the categorial variable 
heatmap_branch_recs$domain_factor <- as.factor(paste(heatmap_branch_recs$Rec_domain,
                                                     heatmap_branch_recs$ID,
                                                     sep='_'))
#Get the final dataset for recombinants only
branch_domains_df_rec <- heatmap_branch_recs[heatmap_branch_recs$evol_type!='non-signif' & 
                                               heatmap_branch_recs$Tox_type=='Rec',c(1,2,25, 26, 7, 5)] %>% unique()

#Plot the ditribution of significant models per domain and recombination event
branch_models_classified_signif[branch_models_classified_signif$evol_type!='non-signif',] %>% 
  ggplot(aes(y=domain_factor, x=fac_for_heatmap, fill = evol_type))+
  geom_tile(alpha=.7)+
  scale_fill_manual(values = c("#a60b0b","#2980b9", '#73C39D'))+
#  branch_domains_df_rec %>% 
#    ggplot(aes(y=domain_factor, x=signif_flag, fill = Rec_domain))+ #change_flag #unknown_flag #Rec_domain
#    geom_tile(alpha=.7)+
#  scale_fill_manual( values=c("#a60b0b", '#2980b9' ,'#bdbd00'))+
  #  scale_fill_manual( values = c('grey', "#B97BD5","#58B09E"))+
  #scale_fill_jama()+
  scale_y_discrete( labels=as.character(unlist(lapply(levels(branch_models_classified_signif$domain_factor), 
                                                      function (x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][2]))))+
  theme_bw()+ylab('Event ID')+xlab('Domain')+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
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
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))+
  coord_fixed()

#Extract data for comparing branch models (significant and raw estimates)
branch_models_comparisions <- branch_models_classified[,c(1,2,3,4,6,7,9)]

#Extract significant results based on different LRT-tests
bneut_signif <- branch_models_comparisions[branch_models_comparisions$evol_type!='non-signif', ]
bfree_signif <- branch_models_comparisions[branch_models_comparisions$signif_flag!='non-signif', ]

#Assert the factors for comparing omega estimates
bneut_signif$for_comparision <- 'bneut'
bfree_signif$for_comparision <- 'bfree'
branch_models_comparisions$for_comparision <- 'All'

#Create the merged dataset for comparisions
branch_models_comparisions<-rbind(branch_models_comparisions, bneut_signif, bfree_signif)

#Summarize mean omega estimates per sets of events
branch_models_comparisions[branch_models_comparisions$omega_f<5,] %>% group_by(Tox_type, for_comparision) %>% 
  summarize(mean_omega= mean(omega_f)) %>% as.data.frame()

#Calculate the linear model to reveal relationships between foreground omega and the number of toxinsin the subtree
lm_num_tox <- lm(omega_f ~ Num_rec_pars,
                 data = branch_models_classified_signif[branch_models_classified_signif$evol_type!='non-signif' & 
                                                          branch_models_classified_signif$omega_f <5  &
                                                          branch_models_classified_signif$Tox_type=='Par',])
#Show model significance
summary(lm_num_tox)

#Plot the dependance between foreground omega and the the number of recombinans and parents
branch_models_classified_signif[branch_models_classified_signif$evol_type!='non-signif' & branch_models_classified_signif$omega_f <5,] %>% 
  ggplot(aes(x= Num_rec_pars, y=omega_f))+
  facet_wrap(~Tox_type, scales = 'free_x')+
  geom_point(aes(color = Rec_domain, shape = Curr_domain, size = 4, alpha =0.7))+
  scale_color_manual( values=c("#a60b0b", '#2980b9' ,'#bdbd00'))+
  geom_smooth(method='lm', color='darkgrey', alpha=0.45, fill = 'grey',  level=0.8)+
  theme_bw()+ylab('Foreground ω')+xlab('Domain')+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
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
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))

#Show the distribution of omega values for sets of events with different significance thresholds for LRT tests
branch_models_comparisions[branch_models_comparisions$omega_f<5,] %>% 
  ggplot(aes(x=for_comparision, y=omega_f, fill = for_comparision))+
  facet_wrap(~Tox_type)+
  geom_violin(width=0.7, alpha = 0.7, trim = T)+
  geom_boxplot(width=0.4, alpha = 0.8)+
  geom_signif(comparisons = list(c(2, 3),
                                 c(1, 3)),
              y_position = c(4.3,5), 
              map_signif_level=T, test='wilcox.test')+
  scale_fill_manual(values =c('#A6A8A0', '#97C997', '#36AB9D'))+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
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
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))

#Show the distribution of omega values for different domains and toxin types (parents and recombinants)
branch_models_classified_signif[branch_models_classified_signif$evol_type!='non-signif' & branch_models_classified_signif$omega_f <5,] %>% 
  ggplot(aes(x= Curr_domain, y=omega_f, fill= Curr_domain))+
  facet_wrap(~Tox_type*Rec_domain)+
  geom_violin(width=0.7, alpha = 0.7, trim = T)+
  geom_boxplot(width=0.4, alpha = 0.8)+
scale_fill_manual( values=c("#a60b0b", '#2980b9' ,'#bdbd00'))+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
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
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))

#Analysis of site models
#Read tables with per-site results using site models
evol_sites_with_relaxed <- read.csv('Site_models_classified_res.csv', header = T, sep = "\t", stringsAsFactors = F)

#Extract significant sites
evol_sites_with_relaxed_signif <- evol_sites_with_relaxed[evol_sites_with_relaxed$signif_flag!='non-signif',]

#Transform relative coordinates to integers
evol_sites_with_relaxed_signif$Coord <- round(evol_sites_with_relaxed_signif$Coord)

#Remove identical coordinates
evol_sites_with_relaxed_signif <- unique(evol_sites_with_relaxed_signif)

#Create a dummy factor for the order on the Y axis
evol_sites_with_relaxed_signif$factor_plot <- paste(evol_sites_with_relaxed_signif$Rec_domain,
                                                    evol_sites_with_relaxed_signif$ID, sep='_')

#Extract the data for selected site models
evol_sites_with_relaxed_site <- evol_sites_with_relaxed_signif[evol_sites_with_relaxed_signif$Model %in% c('m8','m3','m1','m3') & evol_sites_with_relaxed_signif$Tox_type=='Rec',] 

#Extract data for branch-site models 
evol_br_sites_with_relaxed_rec <- evol_sites_with_relaxed_signif[evol_sites_with_relaxed_signif$Tox_type=='Rec' & evol_sites_with_relaxed_signif$Model %in% c('mB','mA'),]
evol_br_sites_with_relaxed_par <- evol_sites_with_relaxed_signif[evol_sites_with_relaxed_signif$Tox_type=='Par' & evol_sites_with_relaxed_signif$Model %in% c('mB','mA'),]

#The function for ordering the events according to the distribution of sites
get_order_for_heatmap = function(sites_df){
  
  #Create an empty matrix with 300 sites
  trans_coords_df<- matrix(ncol= 300, nrow= length(unique(as.factor(sites_df$factor_plot)))) 
  
  #Transform the matrix to the dataframe
  trans_coords_df <- as.data.frame(trans_coords_df)
  
  #Set the order of rows according to the IDs of recombination events
  rownames(trans_coords_df) <- levels(as.factor(sites_df$factor_plot))
  colnames(trans_coords_df) <- 1:300 
  
  i=1 #row index
  #Iterate over recombination events
  for (id in unique(sites_df$factor_plot)){
    print(id)
    
    #Iterate of coordinates
    for (j in 1:300){
      
      #Extract the data for a particular site in the selected event, use the probability as the value
      trans_coords_df[which(rownames(trans_coords_df)==id),j] <- as.numeric(sites_df[sites_df$factor_plot==id & 
                                                                                       sites_df$Coord==j , 8][1])
    }
    i <- i+1
  }
  
  #Set zeroes for non-significant sites
  trans_coords_df[is.na(trans_coords_df)] <- 0
  
  #Create the distance matrix based using the Eucledean distance
  d <- vegdist(trans_coords_df,method="euclidean")
  
  #Perform hierarcheal clustering
  hc_compl <- hclust(d, method = "complete")
  
  #Extract the order from the dendrogram
  dend <- as.dendrogram(hc_compl)
  plot(dend)
  heatmap_order <-labels(dend)
  
  #Return the list of aranged event IDs
  return(heatmap_order)
}

#Get the order of recombination events according to the distribution of seleciton signals in the domain sequences
site_order<- get_order_for_heatmap(evol_sites_with_relaxed_site)
br_sites_rec_order <- get_order_for_heatmap(evol_br_sites_with_relaxed_rec)
br_sites_par_order  <- get_order_for_heatmap(evol_br_sites_with_relaxed_par)

#Assert the order of recombination event IDs according to clusterization baset on sites' distributions
evol_sites_with_relaxed_site$factor_plot<-factor(evol_sites_with_relaxed_site$factor_plot,levels = site_order)
evol_br_sites_with_relaxed_rec$factor_plot<-factor(evol_br_sites_with_relaxed_rec$factor_plot,levels = br_sites_rec_order)
evol_br_sites_with_relaxed_par$factor_plot<-factor(evol_br_sites_with_relaxed_par$factor_plot,levels = br_sites_par_order)

#Plot the heatmap with the distribution of sites for site and branch-site models
ggplot(evol_br_sites_with_relaxed_rec,
       aes(Coord, factor(factor_plot), fill=Site_type, alpha = Site_prob)) +
  geom_tile()+ #alpha = 0.7
  scale_alpha(range = c(0.3,0.7))+
  theme_bw()+
  scale_fill_manual(values = c("#2980b9","#a60b0b", '#73C39D'))+
  labs(y = "Recombination event ID") + 
  labs(x = "Domain coordinate")+
  geom_vline(xintercept=100, linetype="dashed", color = "darkgreen",size=0.3, alpha = 0.99)+ #2.5
  geom_vline(xintercept=200, linetype="dashed", color = "darkgreen",size=0.3, alpha = 0.99)+
  scale_y_discrete( labels=as.character(unlist(lapply(levels(evol_br_sites_with_relaxed_rec$factor_plot), 
                                                      function (x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][2]))))+
  theme( axis.text.x =element_blank(),
         axis.title.x = element_text(face="bold", color="black", 
                                     size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=8),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12),
         axis.ticks.x=element_blank())

#Extract data with branch-site models based contaning overestimations and devoid of them for recombinants and parents 
evol_sites_with_relaxed_rec_passed <- evol_sites_with_relaxed_signif[evol_sites_with_relaxed_signif$Tox_type=='Rec' & evol_sites_with_relaxed_signif$Model %in% c('mB','mA') &  evol_sites_with_relaxed_signif$Site_flag=='passed',]
evol_sites_with_relaxed_par_passed  <- evol_sites_with_relaxed_signif[evol_sites_with_relaxed_signif$Tox_type=='Par' & evol_sites_with_relaxed_signif$Model %in% c('mB','mA') &  evol_sites_with_relaxed_signif$Site_flag=='passed',]

evol_sites_with_relaxed_rec_not_passed <- evol_sites_with_relaxed_signif[evol_sites_with_relaxed_signif$Tox_type=='Rec' & evol_sites_with_relaxed_signif$Model %in% c('mB','mA'),]
evol_sites_with_relaxed_par_not_passed  <- evol_sites_with_relaxed_signif[evol_sites_with_relaxed_signif$Tox_type=='Par' & evol_sites_with_relaxed_signif$Model %in% c('mB','mA'),]

#Arrange events according to the clusterizations of sites' distribution
br_sites_rec_passed_order <- get_order_for_heatmap(evol_sites_with_relaxed_rec_passed)
br_sites_par_passed_order <- get_order_for_heatmap(evol_sites_with_relaxed_par_passed)

br_sites_rec_non_passed_order <- get_order_for_heatmap(evol_sites_with_relaxed_rec_not_passed)
br_sites_par_non_passed_order  <- get_order_for_heatmap(evol_sites_with_relaxed_par_not_passed)

#Assert the dummy factor for the order on the Y axis
evol_sites_with_relaxed_rec_passed$factor_plot<-factor(evol_sites_with_relaxed_rec_passed$factor_plot,levels = br_sites_rec_passed_order)
evol_sites_with_relaxed_par_passed$factor_plot<-factor(evol_sites_with_relaxed_par_passed$factor_plot,levels = br_sites_par_passed_order)

evol_sites_with_relaxed_rec_not_passed$factor_plot<-factor(evol_sites_with_relaxed_rec_not_passed$factor_plot,levels = br_sites_rec_non_passed_order)
evol_sites_with_relaxed_par_not_passed$factor_plot<-factor(evol_sites_with_relaxed_par_not_passed$factor_plot,levels = br_sites_par_non_passed_order)

#Make data for adjacent subpanels
#Site models only
#Extract required columns
sites_heatmap_subpanels <- evol_sites_with_relaxed_site[,c(1,2,4, 11)] %>%  unique()

#Assert the type of unknown parents
sites_heatmap_subpanels$unknown_flag <- 'known'
sites_heatmap_subpanels[sites_heatmap_subpanels$ID %in% unknown_majors, 5] <- 'major'
sites_heatmap_subpanels[sites_heatmap_subpanels$ID %in% unknown_minors, 5] <- 'minor'

#Get the type of patterns of branch models for parents and recombinanats
rec_for_heatmap_sites <- branch_models_classified[branch_models_classified$Tox_type=='Rec' & 
                                                    branch_models_classified$ID!=71,]

#Set dummy factors for retaining correct order of recombination events
rec_for_heatmap_sites$factor_plot <- paste(rec_for_heatmap_sites$Rec_domain, rec_for_heatmap_sites$ID,sep='_')
rec_for_heatmap_sites$factor_plot <- factor(rec_for_heatmap_sites$factor_plot, levels = levels(sites_heatmap_subpanels$factor_plot))
rec_for_heatmap_sites$new_evol_type <- rec_for_heatmap_sites$evol_type

#Extract the significance of site model on the whole
rec_for_heatmap_sites[rec_for_heatmap_sites$signif_flag=='non-signif' & rec_for_heatmap_sites$evol_type!='non-signif', 25] <- 'purifying_b'
rec_for_heatmap_sites[rec_for_heatmap_sites$evol_type=='purifying'& rec_for_heatmap_sites$signif_flag=='signif',25] <- 'purifying_f'

#Branch-site models for recombinants
#Extract required columns
brsites_heatmap_subpanels <- evol_sites_with_relaxed_rec_not_passed[,c(1,2,4, 11)] %>%  unique()

#Assert the type of unknown parents
brsites_heatmap_subpanels$unknown_flag <- 'known'
brsites_heatmap_subpanels[brsites_heatmap_subpanels$ID %in% unknown_majors, 5] <- 'major'
brsites_heatmap_subpanels[brsites_heatmap_subpanels$ID %in% unknown_minors, 5] <- 'minor'

#Get the type of patterns of branch models for parents and recombinanats
rec_for_heatmap_brsites <- branch_models_classified[branch_models_classified$Tox_type=='Rec' & 
                                                      branch_models_classified$ID!=71,]

#Set dummy factors for retaining correct order of recombination events
rec_for_heatmap_brsites$factor_plot <- paste(rec_for_heatmap_brsites$Rec_domain, rec_for_heatmap_brsites$ID,sep='_')
rec_for_heatmap_brsites$factor_plot <- factor(rec_for_heatmap_brsites$factor_plot, levels = levels(brsites_heatmap_subpanels$factor_plot))
rec_for_heatmap_brsites$new_evol_type <- rec_for_heatmap_brsites$evol_type

#Extract the significance of site model on the whole
rec_for_heatmap_brsites[rec_for_heatmap_brsites$signif_flag=='non-signif' & rec_for_heatmap_brsites$evol_type!='non-signif', 25] <- 'purifying_b'
rec_for_heatmap_brsites[rec_for_heatmap_brsites$evol_type=='purifying'& rec_for_heatmap_brsites$signif_flag=='signif',25] <- 'purifying_f'

#Branch-site models for parents
#Extract required columns
brsites_heatmap_subpanels <- evol_sites_with_relaxed_par_not_passed[,c(1,2,4, 11)] %>%  unique()

#Assert the type of unknown parents
brsites_heatmap_subpanels$unknown_flag <- 'known'
brsites_heatmap_subpanels[brsites_heatmap_subpanels$ID %in% unknown_majors, 5] <- 'major'
brsites_heatmap_subpanels[brsites_heatmap_subpanels$ID %in% unknown_minors, 5] <- 'minor'

#Extract required columns
brsites_heatmap_subpanels_par <- evol_sites_with_relaxed_par_not_passed[,c(1,2,4, 11)] %>%  unique()

#Assert the type of unknown parents
brsites_heatmap_subpanels_par$unknown_flag <- 'known'
brsites_heatmap_subpanels_par[brsites_heatmap_subpanels_par$ID %in% unknown_majors, 5] <- 'major'
brsites_heatmap_subpanels_par[brsites_heatmap_subpanels_par$ID %in% unknown_minors, 5] <- 'minor'

#Get the type of patterns of branch models for parents and recombinanats
par_for_heatmap_brsites <- branch_models_classified[branch_models_classified$Tox_type=='Par' & 
                                                      branch_models_classified$ID!=71,]

#Set dummy factors for retaining correct order of recombination events
par_for_heatmap_brsites$factor_plot <- paste(par_for_heatmap_brsites$Rec_domain, par_for_heatmap_brsites$ID,sep='_')
par_for_heatmap_brsites$factor_plot <- factor(par_for_heatmap_brsites$factor_plot, levels = levels(brsites_heatmap_subpanels$factor_plot))
par_for_heatmap_brsites$new_evol_type <- par_for_heatmap_brsites$evol_type

#Extract the significance of site model on the whole
par_for_heatmap_brsites[par_for_heatmap_brsites$signif_flag=='non-signif' & par_for_heatmap_brsites$evol_type!='non-signif', 25] <- 'purifying_b'
par_for_heatmap_brsites[par_for_heatmap_brsites$evol_type=='purifying'& par_for_heatmap_brsites$signif_flag=='signif',25] <- 'purifying_f'


#Plot adjacent subpanels 
par_for_heatmap_brsites %>% 
  ggplot(aes(y=factor_plot, x=Curr_domain, fill=new_evol_type))+
  geom_tile(alpha=.7)+
  scale_fill_manual(values = c('white', "#a60b0b", '#A5D4EF',"#2980b9"))+
  #brsites_heatmap_subpanels %>% 
  #  ggplot(aes(y=factor_plot, x=Tox_type, fill = unknown_flag))+  #unknown_flag #Rec_domain
  #  geom_tile(alpha=.7)+
  #  #scale_fill_manual( values=c("#a60b0b", '#2980b9' ,'#bdbd00'))+
  #  scale_fill_manual( values = c('grey', "#B97BD5","#58B09E"))+  
  scale_y_discrete( labels=as.character(unlist(lapply(levels(par_for_heatmap_brsites$factor_plot), 
                                                      function (x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][2]))))+
  theme_bw()+ylab('Event ID')+xlab('Domain')+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=18),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=11),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=18),
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))+
  coord_fixed()

#Make the dataframe for further plotting the distribution of site types in the branch-site models for parents and recombinants
br_sites_all <- do.call("rbind", replicate(2, branch_models_classified[,c(1,2,3,4, 18, 19, 21)], simplify = FALSE))

#Create the facrot variable for classifying site types
br_sites_all$all_sites_num <- br_sites_all$Brsites_pos
br_sites_all$site_factor <- 'pos'

#Assert conservative sites
br_sites_all[c(247:492),8] <- br_sites_all[c(247:492),6]
br_sites_all[c(247:492),9] <-'cons'

#Make the number of sites for parents negative for plotting
br_sites_all[br_sites_all$Tox_type=='Par',8] <- -br_sites_all[br_sites_all$Tox_type=='Par',8]

#Summarize the number of site types for recombinants and parents within all sets of events and those devoid of overestimations
br_sites_summary_all <- br_sites_all %>% 
  group_by(Rec_domain,Curr_domain,Tox_type,site_factor) %>% 
  summarize(all_sites_num = sum(all_sites_num))

br_sites_summary_passed <- br_sites_all[br_sites_all$Passed_flag=='passed'& br_sites_all$ID!=23 ,] %>%  # & br_sites_all$ID!=23
  group_by(Rec_domain,Curr_domain,Tox_type,site_factor) %>% 
  summarize(all_sites_num = sum(all_sites_num))


#Compare the composition of sites between parents and recombinants per domain 
#Sets of events with the first domain transferred
fisher.test(data.frame(cons=c(230,208),pos=c(27,36))) #p-value = 0.1778
fisher.test(data.frame(cons=c(271,188),pos=c(11,4))) #p-value = 0.3003
fisher.test(data.frame(cons=c(91,132),pos=c(4,24))) #p-value = 0.006565

#Sets of events with the third domain transferred
fisher.test(data.frame(cons=c(724,854),pos=c(36,33))) #p-value = 0.3252
fisher.test(data.frame(cons=c(631,641),pos=c(30,80))) #p-value < 5.967e-06
fisher.test(data.frame(cons=c(358,501),pos=c(85,39))) #p-value = 2.305e-08

#Vizualize the proportion of conservative and positively selected sites for branch-site models
br_sites_summary_passed%>%  #[br_sites_all$Passed_flag=='passed',]
  ggplot(aes(fill=site_factor, y=all_sites_num, x=Curr_domain))+
  facet_wrap(~Rec_domain, scales='free')+
  geom_bar(stat='identity',position='fill', alpha=.7)+ #, position='fill'
  scale_fill_manual(values = c("#2980b9","#a60b0b", '#73C39D'))+
  geom_hline(yintercept = 0)+
  theme_bw()+
  ylab('Number of sites')+ #Foreground ω
  xlab('Domain')+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
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
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))

#Vizualize the dependence between the number of sites, omega estimates, and the number of toxins
branch_models_classified[branch_models_classified$Passed_flag=='passed' & branch_models_classified$omega_f <0.9 &
                           branch_models_classified$Brsites_pos<15, ] %>% 
  ggplot(aes(x=Brsites_pos, y= omega_f))+
  #branch_models_classified %>% ggplot(aes(x=Sites_pos, y= Num_tox)) +
  geom_point(aes(color = Rec_domain, shape = Curr_domain, size = 4, alpha =0.7))+
  facet_wrap(~Tox_type, scales = 'free_x')+
  scale_color_manual( values=c("#a60b0b", '#2980b9' ,'#bdbd00'))+
  geom_smooth(method='lm', color='darkgrey', alpha=0.45, fill = 'grey',  level=0.8)+
  theme_bw()+ylab('Foreground ω')+xlab('Number of sites')+ #Foreground ω Backgound ω  Number of toxins
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
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
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))

#Make linear models to reveal the relashionships between the number of sites and omega values
sites_omega_lm <- lm(Sites_pos~omega_f,
                     data = branch_models_classified_signif[branch_models_classified_signif$evol_type!='non-signif' &
                                                              branch_models_classified_signif$omega_f <5 ,])
sites_omega_lm <- lm(Sites_pos~omega_b,
                     data = branch_models_classified_signif[branch_models_classified_signif$evol_type!='non-signif' &
                                                              branch_models_classified_signif$omega_f <5 ,])
br_sites_omega_lm <- lm(Brsites_pos~omega_b,
                        data = branch_models_classified[branch_models_classified$Passed_flag=='passed' & 
                                                          branch_models_classified$omega_b < 0.9, ])
#Show significance level
summary(sites_omega_lm)
summary(sites_omega_lm)
summary(br_sites_omega_lm)

#Make the dataframe with the total number of sites of particular type for parents and recombinants using site models
total_sites <- #Total number of sites
  branch_models_classified %>% group_by(Rec_domain, Curr_domain) %>% 
  summarize(n_pos = sum(Sites_pos), n_cons = sum(Sites_cons)) %>%  
  as.data.frame() 

#Vizualize the fraction of conservative and positively selected sites using site models
data.frame(num_sites = c(total_sites$n_pos,total_sites$n_cons),
           Curr_domain = rep(total_sites$Curr_domain,2),
           site_type = c(rep('pos',9), rep('cons',9)),
           Domain = rep(total_sites$Rec_domain,2)) %>% 
  ggplot(aes(x=Curr_domain,y=num_sites,fill=site_type)) +
  geom_bar(stat="identity", position=position_fill(), alpha=0.7, col='black', width =1)+
  #geom_bar(stat="identity", alpha=0.7, col='black', width =1)+
  facet_wrap(~Domain)+
  theme_bw()+
  scale_fill_manual(values =c("#2980b9","#a60b0b"))+
  guides(fill= guide_legend(title="Site type"))+
  ylab('Number of sites')+xlab('Domain')+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
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
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))

#Build linear model for reveling dependance between similarity of affected specief and the number of sites
lm_omega <- lm(Simpson_species ~ Brsites_pos, 
               data =  branch_models_classified[branch_models_classified$omega_f<5 ,])
summary(lm_omega) 

#Plot the distribution of omega values in the sets of events regarding changes in affected orders
branch_models_classified[branch_models_classified$omega_f<2 &  
                           branch_models_classified$Brsites_pos<350 &
                           !is.na(branch_models_classified$Simpson_species) ,] %>% 
  ggplot(aes(x=as.factor(Orders_flag), y=omega_f, fill =as.factor(Orders_flag)))+
  geom_violin(width=0.7, alpha = 0.57, trim = T)+
  geom_boxplot(width=0.4, alpha = 0.66)+
  geom_signif(comparisons = list(c(1,2)),
              map_signif_level=T, test='wilcox.test')+
  scale_fill_jama()+
  facet_wrap(~Tox_type*Rec_domain, scales ='free')+
  theme_bw()+
  theme_bw()+ylab('Foreground ω')+xlab('Domain')+ #Foreground ω Backgound ω  Number of sites
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
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
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))

#Plot the dependance between the similarity of the sets of affected species and omega/number of sites
branch_models_classified[branch_models_classified$omega_f<2 &  
                           branch_models_classified$Brsites_pos<350 &
                           branch_models_classified$best_site_model!='relaxation' &
                           !is.na(branch_models_classified$Simpson_species) ,] %>% 
  ggplot(aes(x=Simpson_species, y=Brsites_cons))+
  facet_wrap(~Tox_type*Rec_domain, scales ='free')+
  geom_point(aes(size=Num_rec_pars, color = Rec_domain, shape =Curr_domain),alpha = 0.8 )+
  theme_bw()+
  scale_size(range = c(2,6))+
  scale_color_manual( values=c("#a60b0b", '#2980b9' ,'#bdbd00'))+
  geom_smooth(method='lm', color='darkgrey', alpha=0.45, fill = 'grey')+
  scale_fill_manual( values = c( "#4F499E",'#BD5B0B','darkgrey'),
                     labels = c('Minor','All','Median'))+
  theme_bw()+ylab('Foreground ω')+xlab('Domain')+ #Foreground ω Backgound ω  Number of sites
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
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
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))
