library(ggplot2)
library(ggsci)
library(ggsignif)
library(wesanderson)
library(circlize)
library(dplyr)

setwd("~/phyl/data_for_git/data_for_scripts/for_viz/recombination")

#Ploting the total number of the events
#Reading files with the lists of transferred domains for sets of filtered and raw events 
rec_events_raw <- read.csv('raw_events_by_domains.csv', header = T, sep = "\t", stringsAsFactors = F)
rec_events_filt <- read.csv('flterd_events_by_domains.csv', header = T, sep = "\t", stringsAsFactors = F)
rec_events_unpivot <- read.csv('unpivot_events_by_domains.csv', header = T, sep = "\t", stringsAsFactors = F)

#Counting the total number of the events
rec_raw_events_count <-data.frame(counts=as.vector(table(rec_events_raw)), domains= c('domain1','domain2','domain3'))
rec_filt_events_count <-data.frame(counts=as.vector(table(rec_events_filt )), domains= c('domain1','domain2','domain3'))
rec_unpivot_events_count <-data.frame(counts=as.vector(table(rec_events_unpivot)), domains= c('domain1','domain2','domain3'))

#Drawing pie charts with the number of the events
ggplot(rec_raw_events_count, aes(x='',y=counts,fill=domains)) +
  geom_bar(stat="identity", width=1, alpha=0.6, col='black')+
  coord_polar(theta = "y")+
  theme_bw()+
  geom_text(aes(label = counts),
            position = position_stack(vjust = 0.5), size=10)+
  geom_text(aes(label = counts),
            position = position_stack(vjust = 0.5), size=10)+
  geom_text(aes(label = counts),
            position = position_stack(vjust = 0.5), size=10)+
  scale_fill_manual(values = c("#a60b0b", '#2980b9' ,'#bdbd00'))+
  scale_size(guide = 'none')+
  theme( axis.text.x = element_blank(),
         axis.title.x=element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12),
         axis.ticks = element_blank())+
  guides(fill= guide_legend(title="Domain"))

#Plotting the distribution of breakpoints for revealed recombination events
#Reading the coordinates of the breakpoints for three sets of events (raw, filtered by the length of the domain and by the congruence of phylogenetic trees)
rec_events_raw_man <- read.csv('recomb_manhattan_with_partials.tsv', header = F, sep = "\t", stringsAsFactors = F)
rec_events_filt_man <- read.csv('recomb_manhattan_filtered.tsv', header = F, sep = "\t", stringsAsFactors = F)
rec_events_unpivot_man <- read.csv('recomb_manhattan_filtered_unpivot.tsv', header = F, sep = "\t", stringsAsFactors = F)

#Get log-transformed p-values of recombination detection tests
rec_events_raw_man$V2=-log10(as.numeric(rec_events_raw_man$V2))
rec_events_filt_man$V2=-log10(as.numeric(rec_events_filt_man$V2))
rec_events_unpivot_man$V2=-log10(as.numeric(rec_events_unpivot_man$V2))

#Plotting the breakpoints among the processed toxin sequences
ggplot(rec_events_raw_man, aes(V1,V2,fill=V3)) +
  geom_point(alpha = 0.95,shape=21, color='black',size=4) +
  theme_bw()+
  scale_fill_manual(values = c("#a60b0b", '#2980b9' ,'#bdbd00'))+
  labs(y = "-log10P") + 
  labs(x = "Domain coordinate")+xlim(0.5,300)+
  geom_vline(xintercept=100, linetype="dashed", color = "#606060",size=2, alpha = 0.7)+
  geom_vline(xintercept=200, linetype="dashed", color = "#606060",size=2, alpha = 0.7)+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = c(0.17, 0.76),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12))+
  guides(fill= guide_legend(title="Swapped domain", override.aes = list(size = 6.5)))

#Drawing the breakpoints revealed by FastGear
#Read reference tree to get the order of the toxins and FastGear output
tree_ML <- read.tree("ML_Full_nuc.raxml.support.nwk")
fastgear_events <- read.csv('fastgear_events.csv', header = T, sep = "\t", stringsAsFactors = F)
fastgear_events <- fastgear_events[,c(6,1,2)]

#Add toxins without recombination signals
for (toxin in tree_ML$tip.label){
  if (!toxin %in% fastgear_events$StrainName){
    fastgear_events <- rbind(fastgear_events, c(toxin,0,0))
  }
}

#Transform start and stopp coordinates to the numeric format
fastgear_events$Start <- as.numeric(fastgear_events$Start)
fastgear_events$End <- as.numeric(fastgear_events$End)

#Assert the levels of toxins names according to the tree
fastgear_events$StrainName <- factor(fastgear_events$StrainName, levels=tree_ML$tip.label)

#Create the data for toxins devoid of recombination
fastgear_coords = fastgear_events[fastgear_events$Start==0 & fastgear_events$End==0,c(1,2)]
fastgear_coords$type = 'non_rec'

#Iterate over raw output to add rows to the dataframe with coordinates 
for (row_ind in 1:nrow(fastgear_events)){
  start=fastgear_events[row_ind,2]
  stop=fastgear_events[row_ind,3]
  tox=fastgear_events[row_ind,1]
  print(row_ind)
  for (coord in start:stop){
    if (coord!=0){
      coord_df=data.frame(StrainName=tox, Start=coord, type='rec')
      fastgear_coords <- rbind(fastgear_coords,coord_df)
    }
  }
}

#Drawing the distribution of recombination signals 
ggplot(fastgear_coords, aes( x=Start, y = StrainName, color=type))+
  geom_point(shape =15)+
  scale_color_manual(values = c('lightgrey','#4A56E7'))+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background =element_rect(fill = "lightgrey"), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         axis.ticks.y = element_blank(),
         legend.position = 'none'
  )+
  xlab('Position')


#Draw the distribution of recombination events per dataset 
#Read the table with statistics per event (domain, dataset, the number of toxins, toxin type)
recs_events_heatmap <- read.csv('Num_recs_per_event_for_heatmap.csv', header = T, sep = "\t", stringsAsFactors = F)

#Assert the order of factor variables for correct plotting
recs_events_heatmap$Dataset <- factor(recs_events_heatmap$Dataset, levels=c( 'raw_data', 'filt_partials', 'filt_tree'))
recs_events_heatmap$Tox_type <- factor(recs_events_heatmap$Tox_type, levels=c( 'Recs', 'Mins', 'Majs'))
levels(recs_events_heatmap$Dataset) <- c('Raw', 'Filt_part', 'Filt_tree')
recs_events_heatmap$ID = as.character(recs_events_heatmap$ID)

#Discard rows with inconsitnent events
recs_events_heatmap <- recs_events_heatmap[recs_events_heatmap$Num_tox!=0,]

#Create dummy factor with the number of toxins to draw adjacent panel with domain colors
recs_events_heatmap$factor <- recs_events_heatmap$Num_tox
recs_events_heatmap[recs_events_heatmap$Dataset=='Raw',7] <- recs_events_heatmap[recs_events_heatmap$Dataset=='Raw',7]+100
recs_events_heatmap[recs_events_heatmap$Dataset=='Filt_part',7] <- recs_events_heatmap[recs_events_heatmap$Dataset=='Filt_part',7]+50
recs_events_heatmap$tile_fac <- 'fac'

#Read data for the filtered set of events with specified parents for each of the domains
recs_events_heatmap_filtered <- read.csv('Num_rec_filtered_heatmap.csv', header = T, sep = "\t", stringsAsFactors = F)

#Assert factor variables for plotting order
recs_events_heatmap_filtered$ID = as.character(recs_events_heatmap_filtered$ID)
recs_events_heatmap_filtered <- recs_events_heatmap_filtered[recs_events_heatmap_filtered$Num_tox!=0,]
recs_events_heatmap_filtered$Tox_type <- factor(recs_events_heatmap_filtered$Tox_type, levels = c( 'Recs', 'Mins', 'Majs', 'Dom1', 'Dom2','Dom3'))

#Vizualizing the distribution of evente arranged accoring the the number of toxins and dataset types
ggplot(recs_events_heatmap, aes( x=Tox_type, fill=Tox_type, y = reorder(ID,factor), size = Num_tox))+
  #ggplot(recs_events_heatmap_filtered, aes( x=Tox_type, fill=Tox_type, y = reorder(ID, Num_tox), size = Num_tox))+
  geom_point(shape =21, color='black')+
  theme_bw()+
  scale_fill_manual( values = c('#BD5B0B','#89A05C' , "#4F499E"))+
  #scale_fill_manual( values = c('#BD5B0B','#89A05C' , "#4F499E","#a60b0b", '#2980b9' ,'#bdbd00'))+
  theme( axis.text.x = element_text(color='black', 
                                    size=12, angle = 45, hjust = 1),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         #panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=7.5),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=16), 
         legend.text=element_text(size=14)
  )+
  ylab('Event ID') +
  xlab('Toxin type')+
  coord_fixed(ratio = 0.25)+
  #coord_fixed(ratio = 0.52)+
  facet_wrap(~Dataset)+
  guides( shape=guide_legend(override.aes = list(size=5)), 
          fill=guide_legend(override.aes = list(size=7)))
  
#Drawing the adjacent panels for the distrubution of recombination events with transferred domain and the type of unknown parents
ggplot(recs_events_heatmap, aes( x=Domain, fill=Tox_type, y = reorder(ID,factor)))+ #Unknown_type
    geom_tile(data= recs_events_heatmap, aes( x=tile_fac, fill=Domain, y = reorder(ID,factor)))+  
    theme_bw()+
    #scale_fill_manual( values = c("#B97BD5","#58B09E",'lightgrey'))+ #unknowns
    scale_fill_manual( values=c("#a60b0b", '#2980b9' ,'#bdbd00'))+
    theme( axis.text.x = element_text(color='black', 
                                      size=12, angle = 45, hjust = 1),
           axis.title.x=element_text(face="bold", color="black", 
                                     size=14),
           panel.background = element_blank(), 
           axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           axis.text.y = element_text(color='black', 
                                      size=7.5),
           axis.title.y = element_text(face="bold", color="black", 
                                       size=14),
           legend.title=element_text(face="bold",size=16), 
           legend.text=element_text(size=14)
    )+
    ylab('Event ID') +
    xlab('Unknown flag')+
    #coord_fixed(ratio = 0.25)+
    coord_fixed(ratio = 0.52)+
    facet_wrap(~Dataset)+ #+ facet_wrap(~Dataset+Domain)
  guides( shape=guide_legend(override.aes = list(size=5)), 
          fill=guide_legend(override.aes = list(size=7)))
  

#Summarizing the number of recombinants for each event types classified in a domain-wise way for each of the datasets
recs_events_heatmap[recs_events_heatmap$Tox_type=='Recs',] %>% group_by(Tox_type, Domain, Dataset,Num_tox) %>% 
  summarize(freq=n()) %>% 
  as.data.frame()

#Reading the table with the number of toxins per event for raw and filtered 
num_recs_data <- read.csv('Comparisions_between_number_of_recs_per_filtration_reshaped.csv', header = T, sep = "\t", stringsAsFactors = F)

#Filter events and assert order of factor variables
num_recs_data <- num_recs_data[num_recs_data$Num!=0,]
num_recs_data$Dataset <- factor(num_recs_data$Dataset, levels=c( 'raw_data', 'filt_partials', 'filt_tree'))
num_recs_data$Tox_type <- factor(num_recs_data$Tox_type, levels=c( 'Rec', 'Min', 'Maj'))
levels(num_recs_data$Dataset) <- c('Raw', 'Filt_part', 'Filt_tree')

#Calculate the frequency of certain numbers of the events
daset_comp_with_nums <- num_recs_data %>% group_by(Tox_type, Domain, Dataset , Num) %>% summarize(freq=n()) %>% as.data.frame()

#Plotting the frequency of event-wise toxins' adundance per dataset
ggplot(daset_comp_with_nums, aes( x=Num, fill=Tox_type, y = freq)) +
  geom_bar(stat='identity', position = "dodge", col='black',alpha =0.8)+
  theme_bw()+scale_fill_manual( values = c('#BD5B0B','#89A05C' , "#4F499E"))+
  theme( axis.text.x = element_text(color='black', 
                                    size=12, angle = 45, hjust = 1),
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
  ylab('Number of toxins') +
  xlab('Dataset')+
  facet_wrap(~Dataset+Domain)


#Summarizing the number of unknown parents by their type (minor/major)
num_unknowns <- num_recs_data[num_recs_data$Tox_type=='Rec',] %>% group_by(Unknown_type, Domain, Dataset) %>% summarize(num_unknowns=n()) %>% as.data.frame()

#Discard results without unknowns
num_unknowns <- num_unknowns[num_unknowns$Unknown_type!='no',]
num_unknowns <- rbind(num_unknowns, c('major','domain1','Filt_tree',0))
num_unknowns$num_unknowns <- as.numeric(num_unknowns$num_unknowns)

#Calculate the number of events with unknown parents
sum(num_unknowns[num_unknowns$Dataset=='Raw',4]) #46/120 (38%)
sum(num_unknowns[num_unknowns$Dataset=='Filt_tree',4]) #8/50 (16%)
num_unknowns[num_unknowns$Dataset=='Filt_tree',] #5- major, 3 - minor

#Plotting the total number of unknown parents per dataset
ggplot(num_unknowns, aes( x=Domain, fill=Unknown_type, y = num_unknowns)) + #ggplot(daset_comp_with_unknowns_nums, aes( x=Num, fill=Tox_type, y = freq)) 
  geom_bar(stat='identity', position = "dodge", col='black')+
  theme_bw()+scale_fill_manual( values = c("#B97BD5","#58B09E"))+
  theme( axis.text.x = element_text(color='black', 
                                    size=12, angle = 45, hjust = 1),
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
  ylab('Number of toxins') +
  xlab('Dataset')+
  facet_wrap(~Dataset)+ #facet_wrap(~Dataset+Domain)
  guides( shape=guide_legend(override.aes = list(size=5)), 
          color=guide_legend(override.aes = list(size=5)))

#Reading the table with comparisions of the number of parents before and after filtration
num_parents_compare <- read.csv('Comparisions_between_parents_filtered_reshaped.csv', header = T, sep = "\t", stringsAsFactors = F)

#Assert the order of levels
num_parents_compare$Type <- factor(num_parents_compare$Type, levels = c('All', 'Min', 'Maj'))

#Summarize the number of parents per transferred domain 
num_compare_summary <- num_parents_compare %>% group_by( Num1, Num2, Type) %>% summarise(freq=n())

#Plotting compatisions between the number of parents
ggplot(num_compare_summary, aes( x=Num1, color=Type, y = Num2, shape=Type, size=freq))+
  geom_point(alpha=0.7)+ #, size = 4.5
  xlim(c(0, 31))+
  theme_bw()+xlab('Number before filtration')+ylab('Number after filtration')+
  theme( axis.text.x = element_text(color='black', 
                                    size=16),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=16),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14))+
  scale_color_manual(values = c('#BD5B0B','#89A05C', "#4F499E"))+
  scale_size(range = c(3.1, 9.5), guide = 'none')+
  geom_abline(intercept = 0, slope = 1, color="darkgrey", 
              linetype="dashed", size=1.5, alpha =0.7)+
  guides( color=guide_legend(override.aes = list(size=5), title="Parent type"), 
          shape =guide_legend(override.aes = list(size=5), title="Parent type"),
          size = 'none')

#Comparisions of domain sequences between recombinants and parents 
#Read the table with domain-wise identity for recombination events
ident_stat <- read.csv('Identity_per_dataset.csv', header = T, sep = "\t", stringsAsFactors = F)

#Assert the order of levels
ident_stat$Dataset=factor(ident_stat$Dataset, levels = c('Raw', 'filt_partials', 'filt_tree'))

#Create a dataframe with identity comparisions for filtered set of events 
ident_filt_only <- ident_stat[ident_stat$Dataset=='filt_tree' & ident_stat$Aln_type=='Full_dom',]
ident_filt_only$Domain <- factor(ident_filt_only$Domain)
ident_filt_only$Par_type <- factor(ident_filt_only$Par_type)
ident_filt_only$Dom_type <- factor(ident_filt_only$Dom_type)

#Summarize mean identity between transferred and non-transferred domains
ident_summary <- ident_filt_only %>% group_by(Par_type, Domain, Dom_type) %>% dplyr::summarise(mean_ident=mean(Ident)) %>%  as.data.frame()
#  80, 66, 71

#Plotting the distribution of identity between the domains for a certain dataset
ident_stat[ident_stat$Aln_type=='Full_dom' & ident_stat$Dataset=='filt_partials',] %>% mutate(x2 =  interaction(Par_type, Curr_domain )) %>% 
  ggplot(aes(x2,Ident, Par_type=x2, fill=Par_type))+
  facet_wrap(~Domain)+
  geom_boxplot(alpha=0.7)+
  geom_signif(comparisons = list(c(1, 2),
                                 c(3, 4),
                                 c(5,6)),
              map_signif_level=TRUE)+
  scale_fill_manual( values = c('#89A05C' , "#4F499E"))+
  ylim(c(min(ident_stat$Ident)-5, max(ident_stat$Ident)+5))+
  theme_bw()+xlab('Domain')+ylab('Mean identity')+
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
         legend.title=element_text(face="bold",size=16), 
         legend.text=element_text(size=14))

#Calculaing the difference between identities for transferred and non-transferred domains
ident_full <-ident_stat[ident_stat$Aln_type=='Full_dom', ] 
ident_diff_df <- data.frame(ID=c(), Domain=c(), Par_type=c(), Dataset=c(), Trans_ident=c(), Non_trans_ident=c(), Ident_dif=c())

#Iterate ove the dataframe with domain identities
for (Dataset in unique(ident_full$Dataset)){
  for (ID in unique(ident_full[ident_full$Dataset==Dataset, 1])){
    
    #Extract identity for the transferred domains (from minor of major parents)
    Domain = ident_full[ident_full$Dataset==Dataset & ident_full$ID==ID, 2][1]
    Maj_transferred = ident_full[ident_full$Dataset==Dataset & ident_full$ID==ID & ident_full$Par_type=='Major' & ident_full$Dom_type=='Transferred', 10]
    Min_transferred = ident_full[ident_full$Dataset==Dataset & ident_full$ID==ID & ident_full$Par_type=='Minor'& ident_full$Dom_type=='Transferred', 10]
    
    #Extract identity for  non-transferred domains 
    Maj_non_transferred = ident_full[ident_full$Dataset==Dataset & ident_full$ID==ID & ident_full$Par_type=='Major' & ident_full$Dom_type=='Non_transferred', 10]
    Min_non_transferred = ident_full[ident_full$Dataset==Dataset & ident_full$ID==ID & ident_full$Par_type=='Minor'& ident_full$Dom_type=='Non_transferred', 10]
    
    #Make subdadatfame for identity comparisions within major parents
    maj_pairs <- expand.grid(Maj_transferred,Maj_non_transferred)
    colnames(maj_pairs) <- c('Trans_ident', 'Non_trans_ident')
    maj_pairs$Ident_dif <- maj_pairs$Trans_ident - maj_pairs$Non_trans_ident
    maj_pairs$ID <- ID
    maj_pairs$Domain <- Domain
    maj_pairs$Par_type <- 'Major'
    maj_pairs$Dataset <- Dataset
    
    #Make subdadatfame for identity comparisions within minor parents
    min_pairs <- expand.grid(Min_transferred,Min_non_transferred)
    colnames(min_pairs) <- c('Trans_ident', 'Non_trans_ident')
    min_pairs$Ident_dif <- min_pairs$Trans_ident - min_pairs$Non_trans_ident
    min_pairs$ID <- ID
    min_pairs$Domain <- Domain
    min_pairs$Par_type <- 'Minor'
    min_pairs$Dataset <- Dataset
    
    #Update the dataframe with the overall summary of identity differences
    ident_diff_df<- rbind(ident_diff_df, maj_pairs)
    ident_diff_df<- rbind(ident_diff_df, min_pairs)
  }
}


#The function for calculating mean and dispersion for a selected group from the dataframe
data_summary <- function(data, varname, groupnames){
  #require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-plyr::ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#Summarize mean differences between transferred and non-transferred domains
diff_summary <- data_summary(ident_diff_df, varname="Ident_dif", 
                             groupnames=c("Dataset", "Par_type", 'Domain'))

#transform disperion to mean square error
diff_summary$sd <- diff_summary$sd/sqrt(302)

#Compare the differences between transferred and non-tranferred domains for each type of the events classified according to minor parent
t_test_res_df <- data.frame(Domain=c(), Dataset=c(), parent=c(), p.value=c())
for (dataset in c('Raw', 'filt_partials')){
  for (domain in c('domain1', 'domain2', 'domain3')){
    for (par in c('Major','Minor')){
      t_pval <- t.test(ident_diff_df[ident_diff_df$Dataset=='filt_tree' & ident_diff_df$Par_type==par & ident_diff_df$Domain==domain,3],
                       ident_diff_df[ident_diff_df$Dataset=='Raw' & ident_diff_df$Par_type==par & ident_diff_df$Domain==domain,3])$p.value
      t_test_res <- data.frame(Domain=domain, Dataset=dataset, parent=par, p.value=t_pval)
      t_test_res_df <- rbind(t_test_res_df, t_test_res)
    }
  }
}

#Mark significant comparisions
t_test_res_df$p.value_sig <- t_test_res_df$p.value<=0.05

#Asserting the order of the categorical variable
diff_summary$Dataset=factor(diff_summary$Dataset, levels = c('Raw', 'filt_partials', 'filt_tree'))

#Vizualize total differences between transferred and non-transferred domains
diff_summary %>% mutate(x2 =  interaction(Dataset, Par_type )) %>% 
  ggplot(aes(x2,Ident_dif, Par_type=x2, fill=Par_type))+
  geom_bar(stat = "identity", position=position_dodge())+
  facet_wrap(~Domain)+
  scale_fill_manual( values = c('#89A05C' , "#4F499E"))+
  geom_errorbar(aes(ymin=Ident_dif-sd, ymax=Ident_dif+sd), width=.2,
                position=position_dodge(.9)) +
  theme_bw()+xlab('Dataset')+ylab('Identity Difference')+
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
                                     size=16),
         legend.title=element_text(face="bold",size=16), 
         legend.text=element_text(size=14))

#Reading the data with mismatches within the domain-wise alignments of recombinants and parents
mism_stat <- read.csv('Mismathces_with_len.csv', header = T, sep = "\t", stringsAsFactors = F)

#Normalizing the mnumner of mismatches to the length of the alignment
mism_stat$mism_rate = mism_stat$Num_mism/mism_stat$Aln_length

#Exatracting data for transferred domains per event and parent type
d1_m <- mism_stat[mism_stat$Dom_type=='Transferred' & mism_stat$Domain=='domain1' & mism_stat$Par_type=='Minor' ,11]
d1_M2 <- mism_stat[mism_stat$Dom_type=='Transferred' & mism_stat$Domain=='domain1' & mism_stat$Curr_domain=='domain2' & mism_stat$Par_type=='Major' ,11]
d1_M3 <- mism_stat[mism_stat$Dom_type=='Transferred' & mism_stat$Domain=='domain1' & mism_stat$Curr_domain=='domain3' & mism_stat$Par_type=='Major' ,11]

d2_m <- mism_stat[mism_stat$Dom_type=='Transferred' & mism_stat$Domain=='domain2' & mism_stat$Par_type=='Minor' ,11]
d2_M3 <- mism_stat[mism_stat$Dom_type=='Transferred' & mism_stat$Domain=='domain2' & mism_stat$Curr_domain=='domain3' & mism_stat$Par_type=='Major' ,11]

d3_m <- mism_stat[mism_stat$Dom_type=='Transferred' & mism_stat$Domain=='domain3' & mism_stat$Par_type=='Minor' ,11]
d2_M2 <- mism_stat[mism_stat$Dom_type=='Transferred' & mism_stat$Domain=='domain3' & mism_stat$Curr_domain=='domain2' & mism_stat$Par_type=='Major' ,11]

#Comparing the mismatch rate between major and minor parents 
wilcox.test(d1_m, d1_M2)  #0.7298
wilcox.test(d1_m, d3_M) #0.2658

wilcox.test(d2_m, d2_M3) #0.6667
wilcox.test(d3_m, d2_M2) #0.07237

#Vizualizeing mean mismatch rates for minor and major parents per event type
mism_stat[mism_stat$Dom_type=='Transferred',] %>% mutate(x2 =  interaction(Par_type, Curr_domain )) %>%  
  ggplot(aes(x=x2, y=mism_rate, fill=Par_type)) +
  facet_wrap(~Domain, scales='free_x')+
  geom_boxplot(alpha=0.7)+
  geom_signif(comparisons = list(c(1, 2),
                                 c(1, 3),
                                 c(2, 3)), 
              y_position = rep(max(mism_stat$mism_rate), 3) * c(1, 1.05, 1.1),
              map_signif_level=T, test='wilcox.test')+
  scale_fill_manual( values = c('#89A05C' , "#4F499E"))+
  theme_bw()+xlab('Domain')+ylab('Mismatch rate')+
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
         legend.title=element_text(face="bold",size=16), 
         legend.text=element_text(size=14))

#Vizualizing Circos plot with domain exchanges mapped to the reference phylogeny
#Reading data for circos plot (both sectors and links between them)
sector_data <- read.csv('base_for_circus_plot.tsv', header = T, sep = "\t", stringsAsFactors = F)
link_data <-read.csv('links_for_circus_plot.tsv', header = T, sep = "\t", stringsAsFactors = F)

#Defining color schemes for links
link_data[link_data$Color=='#EEE8AA',4] <- '#a60b0b'
link_data[link_data$Color=='#CD853F',4] <- '#2980b9'
link_data[link_data$Color=='#DC143C',4] <- '#bdbd00'
link_data[link_data$Color=='#556B2F',4] <- '#556B2F'
link_data[link_data$Color=='#5F9EA0',4] <- 'link1'
link_data[link_data$Color=='#DDA0DD',4] <- 'link2'
link_data[link_data$Color=='link1',4] <- '#71C685'
link_data[link_data$Color=='link2',4] <- '#DDA0DD'

#Defining color schemes for sectors
sector_data[sector_data$Color=='#bdbd00',3] <- '#89D5D2'
sector_data[sector_data$Color=='#2980b9',3] <- '#9A73BE'
sector_data[sector_data$Color=='#a60b0b',3] <- '#CD5922'
sector_data$two_rank_cols <- 'black'
sector_data$one_rank_cols <- 'black'


#Make factors for colouring toxins according to their first rank in the BPPRC nomenclature
one_rank_levels <- levels(factor(sector_data$one_rank))
one_rank_cols <- as.vector(wes_palette("Royal1", length(one_rank_levels), type = "continuous"))

#Add colors to the outer sectors using palette from the wesanderson library
for (i in 1:length(one_rank_levels) ){
  sector_data[sector_data$one_rank==one_rank_levels[i],7] <- one_rank_cols[i]
}

#Make empty padding
sectors = sector_data$Sector
s1 = factor(sectors)
circos.par(cell.padding = c(0, 0), gap.degree=0, track.height=0.10)

#Vizualizing sectors
circos.initialize(s1, xlim = c(0, 0.1))
set_track_gap(cm_h(0))
circos.track(sectors, ylim = c(0, 2), bg.col = sector_data$one_rank_cols, track.height = mm_h(5), bg.border = sector_data$one_rank_cols)
circos.track(sectors, ylim = c(0, 1), bg.col = sector_data$Color)

#Vizualizing links
for (ind in 1:nrow(link_data)){
  circos.link(link_data[ind,3], c(0, 0.1), link_data[ind,2], 0, col = link_data[ind,4])
}
circos.clear()

#Summarizing the number of linkes per type (transferred domains)
link_data$link_stat <-factor(link_data$Color)
levels(link_data$link_stat) <- c('Domain2', 'Domains 1 and 2', 'Domains 1 and 3',  'Domain1','Domain3', 'Domains 2 and 3', 'Recombinants')
link_data$link_stat <-factor(link_data$link_stat, levels = c(  'Domain1','Domain2', 'Domain3',
                                                               'Domains 1 and 2', 'Domains 1 and 3', 'Domains 2 and 3', 'Recombinants'))

#Drawing the distribution of link types
link_data %>% dplyr::group_by(link_stat) %>%   dplyr::mutate(count_name_occurr = n()) %>% 
  ggplot(aes(x=reorder(link_stat,-count_name_occurr), fill=link_stat)) +
  geom_bar(stat="count", color='black' ,alpha =0.99)+
  scale_fill_manual(values = c("#a60b0b", "#2980b9", "#bdbd00","#556B2F",'#71C685',"#DDA0DD", 'black'))+
  theme_bw()+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=16),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=14),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=16),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12))+
  xlab('Transferred domain') +
  ylab('Frequency')+
  guides(fill= guide_legend(title="Transferred domain"))

#Summarizing the number of toxins regarding their role in recombination events (recombinants, parents and both)
rec_ev_stat <- data.frame(counts=as.vector(table(sector_data$Color)), toxin_group= c('Parents', 'Recombinants','Both','Unaffected'))

#Plotting the number of toxin types
ggplot(rec_ev_stat, aes(x=reorder(toxin_group,-counts), y=counts,  fill=toxin_group)) +
  geom_bar(stat="identity",width=1, alpha=0.8, col='black')+
  theme_bw()+scale_fill_manual(values = c( '#89D5D2','#9A73BE','#CD5922','grey'))+
  theme( axis.text.x = element_blank(),
         axis.title.x = element_text(face="bold", color="black", 
                                     size=16),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=14),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=16),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12),
         axis.ticks.x = element_blank())+
  ylab('Number of toxins')+
  xlab("Toxins' group")+
  guides(fill= guide_legend(title="Toxins' group"))

#Analyzing edges and nodes of the recombination graph
#Reading the properties of the graph connected components
graph_stat <-  read.csv('Rec_components_stat_long.csv', header = T, sep = "\t", stringsAsFactors = F) 

#Make summary of the total number of edges and nodes per their type, i.e., toxin type and transferred domain
graph_edges <- graph_stat[graph_stat$Type=='Link_type',] %>% group_by(Value) %>% summarize(num_occur = sum( Num))%>% as.data.frame()
graph_nodes <- graph_stat[graph_stat$Type=='Node_type',] %>% group_by(Value) %>% summarize(num_occur = sum( Num))%>% as.data.frame()

#Summarize the number of recombination events and toxins in graph connected components
graph_events_num <- graph_stat[graph_stat$Type=='Event_ID',] %>% group_by(Component) %>% summarize(num_occur = sum( Num))%>% as.data.frame()
graph_events_type <- graph_stat[graph_stat$Type=='Event_type',] %>% group_by(Component, Value) %>% summarize(num_occur = sum( Num))%>% as.data.frame()
graph_toxins_num <- graph_stat[graph_stat$Type=='Node_type',] %>% group_by(Component, Value) %>% summarize(num_occur = sum( Num))%>% as.data.frame()

sum(graph_events_num[graph_events_num$num_occur==1,2]) #16
sum(graph_edges$num_occur) # 310 
sum(graph_nodes$num_occur) # 180

#Make a dataframe with the distribution of events of certain types classified according to transferred domains from minor parents
graph_events_type$total_num <- graph_events_type$num_occur
graph_events_type[graph_events_type$Component==3,4] <- 4
graph_events_type[graph_events_type$Component==6,4] <- 4
graph_events_type[graph_events_type$Component==1,4] <- 19
graph_events_type[graph_events_type$Component==14,4] <- 2
colnames(graph_events_type)[2] <- 'Domain'

#Assert the levels of factors for the total number of toxin types in graph components
graph_toxins_num$Value <- factor(graph_toxins_num$Value, levels = c('par','rec','mixed'))
levels(graph_toxins_num$Value) <- c('Parents', 'Recombinants','Both')
graph_toxins_num$total_num <- graph_toxins_num$num_occur
graph_toxins_num[graph_toxins_num$Component==2,4] <- 2


#Vizualizing the composition of nodes and edges in connected components of the recombination graph
ggplot(graph_toxins_num, aes(x=reorder(Component,-total_num), y= num_occur, fill = Value)) +
  geom_bar(stat="identity", color='black' ,alpha =0.8)+
  #scale_fill_manual(values = c("#a60b0b", '#2980b9' ,'#bdbd00'))+
  scale_fill_manual(values = c("#9A73BE", '#CD5922' ,'#89D5D2'))+
  theme_bw()+
  theme( axis.text.x = element_text(color='black', 
                                    size=14),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=16),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=14),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=16),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12))+
  #legend.position = 'none')+
  xlab('Component ID') +
  ylab('Number of toxins')+ #events
  guides(fill= guide_legend(title="Toxins' type"))
