library(ggplot2)
library(ggsci)
library(dplyr)
library(ggsignif)
library(scatterpie)
library(ComplexUpset)
library(reshape2)

setwd("~/phyl/data_for_git/data_for_scripts/for_viz/specificity")

#Vzualization of the total statistics of toxicity data
#Read the table with overall data per strains and toxins
overall_toxicity = read.csv('Overall_toxicity_stat.csv', header = T, sep = "\t", stringsAsFactors = F)

#Smmmarize the number of annotated strains and toxins
overall_num <- overall_toxicity %>% group_by(Tox_type,Strain.Toxin) %>% summarise(num=n()) %>% as.data.frame()
overall_num$Tox_type <- factor(overall_num$Tox_type, levels=c('Rank4','Rank3','Strains'))
colnames(overall_num)[1] <- 'Group'

#Vizualize the amount of annotations per strain and toxin
ggplot(overall_num, aes(x=Group,fill=Group)) +
  geom_bar(alpha = 0.7)+scale_fill_aaas()+
  theme_bw()+xlab('Group')+ylab('Number of annotations')+
  theme( axis.text.x = element_text(color='black', 
                                    size=14),
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
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))

#Calculate the frequency of toxins and strains in the context of affected host orders
per_type_num <- overall_toxicity %>% group_by(Tox_type,Strain.Toxin, Host_type) %>% dplyr::summarise(num=n()) %>% as.data.frame()
per_type_num$Tox_type <- factor(per_type_num$Tox_type, levels=c('Rank4','Rank3','Strains'))
colnames(per_type_num)[1] <- 'Group'

#Vizualize the number of affected host orders
ggplot(per_type_num[per_type_num$Host_type=='Orders',], aes(x=num,fill=Group))+
  facet_wrap(~Group, scales = 'free_x')+
  geom_bar(alpha = 0.7)+scale_fill_aaas()+
  theme_bw()+ylab('Frequency')+xlab('Number of hosts')+
  theme( axis.text.x = element_text(color='black', 
                                    size=14),
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

#Summarize the number of individual assays assigned to different host orders
assays_num <- overall_toxicity[overall_toxicity$Host_type=='Orders' & !overall_toxicity$Host %in% c('Platyhelminthes', 'Protists'),] %>% 
  group_by(Tox_type,Host) %>% dplyr::summarise(num=n(), sum_hosts=sum()) %>% as.data.frame() 

#Add the total number of assays per group (strains/toxins)
assays_num[assays_num$Tox_type=='Rank3',4] <-221
assays_num[assays_num$Tox_type=='Rank4',4] <-322
assays_num[assays_num$Tox_type=='Strains',4] <-164

#Calculate the percentage of certain orders 
assays_num$ratio <- assays_num$num/assays_num$sum_hosts

#Create the dataframe for vizualization
assays_num <- assays_num %>% group_by(Tox_type) %>% arrange(Tox_type, ratio) %>% as.data.frame()
assays_num$Host <-factor(assays_num$Host, levels=rev(assays_num[assays_num$Tox_type=='Rank3',2]))
assays_num$Tox_type <- factor(assays_num$Tox_type, levels=c('Rank4','Rank3','Strains'))

#Illustrate the percentage of certain host orders per dataset type
ggplot(assays_num, aes(y=ratio,fill=Host, x=Tox_type))+
  geom_bar(stat='identity',alpha = 0.7, position = "fill")+scale_fill_lancet()+
  theme_bw()+ylab('Frequency')+xlab('Group')+
  theme( axis.text.x = element_text(color='black', 
                                    size=14),
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

#Vizualization of specificity changes during domain exchanges
#Read the data with event-wise distribution of affected hosts 
data_for_tox_heatmap <- read.csv('Hosts_for_heatmap_stat_per_species.csv', header = T, sep = "\t", stringsAsFactors = F)

#Read the table with domain-wise identities between parents and recombinants within recobination events
ident_stat <- read.csv('Identity_per_dataset_with_unknowns.csv', header = T, sep = "\t", stringsAsFactors = F)

#Extract identity data for transferred domains within fltered set of events 
ident_stat <- ident_stat[ident_stat$Dataset == 'filt_tree' & ident_stat$Dom_type =='Transferred' & ident_stat$Aln_type  =='Extr', ]
identity_for_filt_df <- ident_stat %>% group_by(ID, Par_type, Domain) %>% 
  dplyr::summarize(mean_ident = mean(Ident), dispertion = sd(Ident)) %>% as.data.frame()

#Create a dummy factor for plotting on the Y axis
data_for_tox_heatmap$y <- paste(data_for_tox_heatmap$Change_flag, data_for_tox_heatmap$Domain , 
                                data_for_tox_heatmap$ID, data_for_tox_heatmap$Tox_type, sep = '_')

#Assert pseudocoordinate for plotting scatterplots
data_for_tox_heatmap$num_for_scatter <- 1

#Specify the size of circles depending on the number of recombinants
data_for_tox_heatmap$radius_for_scatter <- data_for_tox_heatmap$Num_recs

#Assert mean identity dummy
data_for_tox_heatmap$mean_ident <- 100

#Extract mean identity for sets of minor and major parents within recombination events
identity_minors <- identity_for_filt_df[identity_for_filt_df$Par_type=='Minor',]
identity_majors <- identity_for_filt_df[identity_for_filt_df$Par_type=='Major',]

#Add identity to the dataframe with event-wise specificity spectrum
for (ID in unique(data_for_tox_heatmap$ID)){
  data_for_tox_heatmap[data_for_tox_heatmap$ID==ID & data_for_tox_heatmap$Tox_type=='Min',13] <- identity_minors[identity_minors$ID == ID,4]
  data_for_tox_heatmap[data_for_tox_heatmap$ID==ID & data_for_tox_heatmap$Tox_type=='Maj',13] <- identity_majors[identity_majors$ID == ID,4]
}

#Cange pseudocoordinates for scatterplots belonging to major and minor events
data_for_tox_heatmap[data_for_tox_heatmap$Tox_type=='Min',12] <-  data_for_tox_heatmap[data_for_tox_heatmap$Tox_type=='Min',8] +data_for_tox_heatmap[data_for_tox_heatmap$Tox_type=='Min',12]
data_for_tox_heatmap[data_for_tox_heatmap$Tox_type=='Maj',12] <-  data_for_tox_heatmap[data_for_tox_heatmap$Tox_type=='Maj',9] +data_for_tox_heatmap[data_for_tox_heatmap$Tox_type=='Maj',12]

#Make the variable for the order on the Y axis
data_for_tox_heatmap$type_for_scatter <-  paste(data_for_tox_heatmap$Change_flag, data_for_tox_heatmap$Domain , 
                                                data_for_tox_heatmap$ID,  sep = '_')

#Summarize the number of affeted hots corrsponding to their orders 
#Third rank of toxis within the set of events with different affected orders for parents and recombinants
d3<- dcast(data_for_tox_heatmap[data_for_tox_heatmap$Mode=='Rank3' & data_for_tox_heatmap$Host!='Unknown' & data_for_tox_heatmap$Change_flag =='Different',c(4,5, 11, 12, 13, 14)], 
           type_for_scatter+Tox_type+ radius_for_scatter+ mean_ident~ Host, value.var = "num_for_scatter")
d3$Tox_type <- factor(d3$Tox_type, levels = c('Rec', 'Min', 'Maj'))
d3 <- d3 %>% 
  mutate(type_num = as.numeric(as.factor(type_for_scatter)), 
         tox_num = as.numeric(Tox_type))

#Fourth rank of toxis within the set of events with different affected orders for parents and recombinants
d4<- dcast(data_for_tox_heatmap[data_for_tox_heatmap$Mode=='Rank4' & data_for_tox_heatmap$Host!='Unknown' & data_for_tox_heatmap$Change_flag =='Different',c(4,5, 11, 12,13, 14)], 
           type_for_scatter+Tox_type+ radius_for_scatter +mean_ident ~ Host, value.var = "num_for_scatter")
d4$Tox_type <- factor(d4$Tox_type, levels = c('Rec', 'Min', 'Maj'))
d4 <- d4 %>% 
  mutate(type_num = as.numeric(as.factor(type_for_scatter)), 
         tox_num = as.numeric(Tox_type))

#Fourth rank of toxis within the set of events with same affected orders for parents and recombinants
d5 <- dcast(data_for_tox_heatmap[data_for_tox_heatmap$Mode=='Rank4' & data_for_tox_heatmap$Host!='Unknown' & data_for_tox_heatmap$Change_flag =='Same',c(4,5, 11, 12,13, 14)], 
            type_for_scatter+Tox_type+ radius_for_scatter + mean_ident~ Host, value.var = "num_for_scatter")
d5$Tox_type <- factor(d5$Tox_type, levels = c('Rec', 'Min', 'Maj'))
d5 <- d5 %>% 
  mutate(type_num = as.numeric(as.factor(type_for_scatter)), 
         tox_num = as.numeric(Tox_type))

#Third rank of toxis within the set of events with same affected orders for parents and recombinants
d6<- dcast(data_for_tox_heatmap[data_for_tox_heatmap$Mode=='Rank3' & data_for_tox_heatmap$Host!='Unknown' & data_for_tox_heatmap$Change_flag =='Same',c(4,5, 11, 12,13, 14)], 
           type_for_scatter+Tox_type+ radius_for_scatter + mean_ident~ Host, value.var = "num_for_scatter")
d6$Tox_type <- factor(d6$Tox_type, levels = c('Rec', 'Min', 'Maj'))
d6 <- d6 %>% 
  mutate(type_num = as.numeric(as.factor(type_for_scatter)), 
         tox_num = as.numeric(Tox_type))

#Third rank of toxis coupled with strains within the set of events with different affected orders for parents and recombinants
d_strains_dif <- dcast(data_for_tox_heatmap[data_for_tox_heatmap$Mode=='Strains_add' & data_for_tox_heatmap$Host!='Unknown' & data_for_tox_heatmap$Change_flag =='Different',c(4,5, 11, 12,13, 14)], 
                       type_for_scatter+Tox_type+ radius_for_scatter + mean_ident~ Host, value.var = "num_for_scatter")
d_strains_dif$Tox_type <- factor(d_strains_dif$Tox_type, levels = c('Rec', 'Min', 'Maj'))
d_strains_dif <- d_strains_dif %>% 
  mutate(type_num = as.numeric(as.factor(type_for_scatter)), 
         tox_num = as.numeric(Tox_type))

#Third rank of toxis coupled with strains within the set of events with same affected orders for parents and recombinants
d_strains_same <- dcast(data_for_tox_heatmap[data_for_tox_heatmap$Mode=='Strains_add' & data_for_tox_heatmap$Host!='Unknown' & data_for_tox_heatmap$Change_flag =='Same',c(4,5, 11, 12,13, 14)], 
                        type_for_scatter+Tox_type+ radius_for_scatter + mean_ident~ Host, value.var = "num_for_scatter")
d_strains_same$Tox_type <- factor(d_strains_same$Tox_type, levels = c('Rec', 'Min', 'Maj'))
d_strains_same <- d_strains_same %>% 
  mutate(type_num = as.numeric(as.factor(type_for_scatter)), 
         tox_num = as.numeric(Tox_type))

#Vizualise the distribution of pie charts showing the percentage of affected orders for toxins
ggplot() + geom_scatterpie(aes(x=tox_num, y=type_num, r=0.37), data=d_strains_same, 
                           #cols=c('Coleoptera', 'Diptera', 'Hemiptera' ,'Human', 'Lepidoptera'),  alpha =.7)+                           
                           #cols=c('Coleoptera', 'Diptera', 'Hemiptera' ,'Human', 'Lepidoptera', 'Nematoda'),  alpha =.7)+
                           #cols=c('Coleoptera','Lepidoptera', 'Nematoda'),  alpha =.7)+  
                           #cols=c('Lepidoptera', 'Nematoda'),  alpha =.7)+  
                           cols=c('Coleoptera','Diptera','Lepidoptera', 'Nematoda'),  alpha =.7)+
  coord_equal()+
  #scale_fill_manual(values = c('#42B540FF', '#ED0000FF', '#925E9FFF', '#FDAF91FF', '#00468BFF', '#0099B4FF') )+
  #scale_fill_manual(values = c('#42B540FF', '#00468BFF', '#0099B4FF') )+ #rank4 same
  #scale_fill_manual(values = c('#00468BFF', '#0099B4FF') )+ #rank3 same
  scale_fill_manual(values = c('#42B540FF', '#ED0000FF','#00468BFF', '#0099B4FF') )+ #strains_same
  
  theme_bw()+ylab('Recombination ID')+xlab('Group')+
  scale_x_continuous(breaks=c(1,2,3), labels=c("Rec", "Min", "Maj")) + 
  #scale_y_continuous(breaks=c(1:17), labels=c("34","38", "87","95",'107',"71","109", "14",'2',"23","24","27","29",'37','51',"52",'98'))+ #strains diff
  scale_y_continuous(breaks=c(1:13), labels=c("119","56","62",'96','104', "17","19","44",'49',"53",'6','73', "8"))+ #strains same
  #scale_y_continuous(breaks=c(1:11), labels=c("34","38", "87","95","71","109", "14","23","27","29","52"))+ #rank4 diff
  #scale_y_continuous(breaks=c(1:14), labels=c("34","38", "87","95","71","109", "14","23","24","27","29",'51',"52",'98'))+ #rank3 diff
  #scale_y_continuous(breaks=c(1:8), labels=c("119","56", "17","19","44","53", "8","98"))+ #rank4 same
  #scale_y_continuous(breaks=c(1:8), labels=c("119","56","62", "17","19","44","53", "8"))+ #rank3 same
  theme(axis.text.x = element_text(color='black', 
                                   size=14),
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
        legend.title=element_text(face="bold",size=18), 
        legend.text=element_text(size=16))

#Summarize mean identity between domains of recombinants and parents (both minor and major) 
mearged_ident_df_for_scatter  <- dcast(data_for_tox_heatmap[data_for_tox_heatmap$mean_ident>0,c(2,4,5,6, 11, 12,13, 14)], 
                                       type_for_scatter+Tox_type+ radius_for_scatter + mean_ident + Mode +Change_flag ~ Host, value.var = "num_for_scatter") %>% 
  mutate(Tox_type = factor(Tox_type, levels = c('Rec', 'Min', 'Maj')))%>% 
  mutate(type_num = as.numeric(as.factor(type_for_scatter)), 
         tox_num = as.numeric(Tox_type))
mearged_ident_df_for_scatter <- mearged_ident_df_for_scatter[mearged_ident_df_for_scatter$Tox_type!='Rec',]

#Plot mean identity and the number of toxins for the adjacent subpanel
ggplot() +
  geom_point(aes(x=Tox_type, y=type_for_scatter,   
                 fill = mean_ident,
                 size = radius_for_scatter), 
             alpha =0.9, 
             shape = 21, color = 'black', data=mearged_ident_df_for_scatter[mearged_ident_df_for_scatter$Tox_type!='Rec' & mearged_ident_df_for_scatter$Mode=='Strains_add' & mearged_ident_df_for_scatter$Change_flag=='Same',])+
  #facet_wrap(~Mode+Change_flag, scales = 'free')+
  scale_fill_gradient(low = "#FDFDFD", high = "darkred", limit = c(50,100), space = "Lab")+
  scale_size(range = c(2,3))+
  
  theme_bw()+
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
         legend.text=element_text(size=14))+
  ylab('Event ID') +
  xlab('Domain')+
  coord_fixed()+
  guides(fill=guide_legend(override.aes = list(size=7), title="Identity"))

#Prepare the data for applying logistic regression 
#Mearge event summaries with or without change in affected orders
strains_mearged <- rbind(d_strains_same[, c(1,2,3,4)], d_strains_dif[, c(1,2,3,4)])
rank3_mearged <- rbind(d3[, c(1,2,3,4)], d6[, c(1,2,3,4)])
rank4_mearged <- rbind(d4[, c(1,2,3,4)], d5[, c(1,2,3,4)])

#Extract te flag showing changes in affected orders, event IDs, and transferred domains
strains_mearged$same_flag <- unlist(lapply(strains_mearged$type_for_scatter, function (x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][1]))
rank3_mearged$same_flag <- unlist(lapply(rank3_mearged$type_for_scatter, function (x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][1]))
rank4_mearged$same_flag <- unlist(lapply(rank4_mearged$type_for_scatter, function (x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][1]))

rank4_mearged$ID <- unlist(lapply(rank4_mearged$type_for_scatter, function (x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][3]))
rank3_mearged$ID <-unlist(lapply(rank3_mearged$type_for_scatter, function (x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][3]))
strains_mearged$ID<- unlist(lapply(strains_mearged$type_for_scatter, function (x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][3]))

rank4_mearged$domain<- unlist(lapply(rank4_mearged$type_for_scatter, function (x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][2]))
rank3_mearged$domain<- unlist(lapply(rank3_mearged$type_for_scatter, function (x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][2]))
strains_mearged$domain<- unlist(lapply(strains_mearged$type_for_scatter, function (x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][2]))

#Specifiy the type of the dataset used for obtaining host specificity 
rank4_mearged$Type <- 'Rank4'
rank3_mearged$Type <- 'Rank3'
strains_mearged$Type <- 'Strains'

#Exclude rows with recombinants only
strains_mearged <- strains_mearged[strains_mearged$Tox_type!='Rec',]
rank3_mearged <- rank3_mearged[rank3_mearged$Tox_type!='Rec',]
rank4_mearged <- rank4_mearged[rank4_mearged$Tox_type!='Rec',]

#Create categorial variables for performing logistic regression
rank4_mearged$logit_var <- rank4_mearged$same_flag=='Same'
rank3_mearged$logit_var <- rank3_mearged$same_flag=='Same'
strains_mearged$logit_var <- strains_mearged$same_flag=='Same'

rank4_mearged$Tox_type <- factor(rank4_mearged$Tox_type )
rank3_mearged$Tox_type <- factor(rank3_mearged$Tox_type)
strains_mearged$Tox_type <- factor(strains_mearged$Tox_type)

#Perfrom logistic regression to reveal factors contributing to alterations in host orders
logit_rank4 <- glm(logit_var ~ radius_for_scatter+mean_ident+Tox_type , data = rank4_mearged, family = "binomial")
logit_rank3 <- glm(logit_var ~ radius_for_scatter+mean_ident+Tox_type , data = rank3_mearged, family = "binomial")
logit_strains <- glm(logit_var ~ radius_for_scatter+mean_ident+Tox_type , data = strains_mearged, family = "binomial")

summary(logit_rank4)
summary(logit_rank3)
summary(logit_strains)

#Read the table with estimates of the overlap coefficient between sets of affected species
Simpson_species_stat <- read.csv('Simpson_stat_for_events.csv', header = T, sep = "\t", stringsAsFactors = F)

#Assert factor levels
Simpson_species_stat$Tox_type <- factor(Simpson_species_stat$Tox_type, levels = c('min','maj'))

#Plot the relashionship between overlap coefficient and domain-wise identity
ggplot(Simpson_species_stat, aes(color=Domain, x=Ident, y=Simp_coeff, shape=Tox_type, size =Tox_num))+
  geom_point(alpha = 0.9)+
  facet_wrap(~Mode)+
  scale_size(range = c(2,7))+
  scale_color_manual( values=c("#a60b0b", '#2980b9' ,'#bdbd00'))+
  theme_bw()+
  labs(y = "Overlap coefficient") + 
  labs(x = "Identity percent")+
  theme( axis.text.x = element_text(color='black', 
                                    size=16),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=18),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=16),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=18),
         legend.title=element_text(face="bold",size=16), 
         legend.text=element_text(size=14))

#Apply linear model to assess the role of factors associated with the Overlap coefficient
lm_simpson <- lm(Simp_coeff ~ Tox_type*Domain, 
                 data = Simpson_species_stat[Simpson_species_stat$Mode=='Strains_add' & Simpson_species_stat$Domain!='domain2',-c(1,6)])
summary(lm_simpson)

#Vizualize the distribution of overlap coefficient estimates
ggplot(Simpson_species_stat[Simpson_species_stat$Domain!='domain2',], aes(fill=Tox_type, x=Simp_coeff))+
  facet_wrap(~Mode+Domain)+
  scale_fill_manual( values = c("#4F499E",'#89A05C' ), labels = c("Minor",'Major'))+
  geom_density(alpha = 0.8,  adjust = 1/1.5)+
  #geom_histogram(alpha = 0.8, position = position_dodge2(width = 0.9, preserve = "single"), bins = 5, color='black')
  theme_bw()+ylab('Density')+xlab('Overlap Coefficient')+
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

#Plotting the number of affected hosts for toxins and strains
#Read the table with the number of hosts per toxin
num_hosts_per_protein_group <- read.csv('num_hosts_per_toxin_group.csv', header = T, sep = "\t", stringsAsFactors = F)
num_hosts_per_protein_group <- num_hosts_per_protein_group[num_hosts_per_protein_group$Num_hosts>0,]

#Extract data for true majors only
num_hosts_strict <- num_hosts_per_protein_group[num_hosts_per_protein_group$Tox_type!='par_all',]

#Assert the order of levels for plotting
num_hosts_strict$Tox_type <- factor(num_hosts_strict$Tox_type, levels = c('rec', 'par_strict', 'ref_no_rec'))

#Summarize mean etimates per group
num_hosts_strict[num_hosts_strict$Num_hosts<=15, ] %>% group_by(Host_type, Group, Tox_type) %>% 
  summarize(mean_host_num = mean(Num_hosts))

#Vizualize the number of hosts
ggplot(num_hosts_strict[num_hosts_strict$Host_type=='Orders' & num_hosts_strict$Num_hosts<=100,], 
       aes(fill=Tox_type, y=Num_hosts, x=Tox_type))+
  facet_wrap(~Group, scales ='free')+
  geom_violin(width=1.0, alpha = 0.7, trim = T)+
  geom_boxplot(width=0.4, alpha = 0.8)+ #for species
  #geom_signif(comparisons = list(c("rec","par_strict"),
  #                               c('rec', 'ref_no_rec'),
  #                               c('par_strict', 'ref_no_rec')),
  #                               y_position = c(12, 14.6,13.3), 
#            map_signif_level=T, test='wilcox.test')+
  scale_fill_manual(labels = c("Recombinants",'Parents','Unaffected'), values = c("#CD5922", '#9A73BE' ,'grey'))+
  theme_bw()+ylab('Number of affected hosts')+xlab('Group')+
  theme( axis.text.x = element_text(color='black', 
                                    size=14),
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
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))+
  guides(fill= guide_legend(title="Toxin type"))

#Read the data with the number of toxins per assembly and strain regarding nomenclature ranks
num_tox_for_strains <- read.csv('Num_tox_for_strains_per_rank.csv', header = T, sep = "\t", stringsAsFactors = F)
num_tox_for_asmbls <- read.csv('Num_tox_for_assemblies_per_rank.csv', header = T, sep = "\t", stringsAsFactors = F)

#Calculate median values
num_tox_for_strains[num_tox_for_strains$Rank=='Rank3',] %>% summarize(m=median(Num)) #1
num_tox_for_asmbls[num_tox_for_asmbls$Rank=='Rank3',] %>% summarize(m=median(Num)) #3

#Plot the number of toxins per assembly/srain regarding ranks
ggplot(num_tox_for_strains, aes(fill=Rank, x=Num))+
  geom_bar(stat='count',alpha = 0.7, position = position_dodge2(width = 0.9, preserve = "single"), color='black')+scale_fill_lancet()+
  theme_bw()+ylab('Frequency')+xlab('Number of toxins')+
  theme( axis.text.x = element_text(color='black', 
                                    size=14),
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

#Comparision between the sets of toxins fof strains containing parents and recombinants
#Read the data with comparisions using different metrics
common_tox_df <-read.csv('Num_common_tox_pars_recs.csv', header = T, sep = "\t", stringsAsFactors = F)

#Extract desired metrics (comparision between parents and recombinants/parents only)
common_tox_df <- common_tox_df[common_tox_df$Coeff %in% c('Jaccard', 'Par_coeff') & common_tox_df$Par_type == 'Strict',]

#Calculate the number of events with completely different sets of toxins
nrow(common_tox_df[common_tox_df$Comon_num==0 & common_tox_df$Coeff %in% c('Jaccard') & common_tox_df$Dataset=='Strains', ])
nrow(common_tox_df[common_tox_df$Comon_num==0 & common_tox_df$Coeff %in% c('Par_coeff') & common_tox_df$Dataset=='Strains', ])

#SUmmarize mean similarity coefficients for events with at least one common toxin 
common_tox_df[common_tox_df$Coeff %in% c('Jaccard', 'Par_coeff') & common_tox_df$Dataset == 'Assemblies' & common_tox_df$Comon_num!=0, ] %>% 
  group_by(Coeff) %>% summarize (coeff=median(Comon_num))

#Vizualize the distribution of similarity coefficients 
ggplot(common_tox_df[common_tox_df$Coeff %in% c('Jaccard') & common_tox_df$Par_type == 'Strict', ], 
       aes(fill=Dataset, x=Comon_num))+
  #ggplot(common_tox_df[common_tox_df$Coeff %in% c('Jaccard', 'Par_coeff') & common_tox_df$Par_type == 'Strict' & common_tox_df$Comon_num!=0, ], 
  #    aes(fill=Coeff, x=Coeff, y=Comon_num))+
  geom_histogram(alpha = 0.8, position = position_dodge2(width = 0.9, preserve = "single"), bins = 20, color='black')+
  #geom_density(alpha = 0.8)+
  scale_fill_jama()+
  #geom_boxplot(alpha = 0.8)+
  #scale_fill_manual(labels = c("Recombinants",'Parents'), values = c("#CD5922", '#9A73BE' ))+
  #geom_signif(comparisons = list(c(1,2)),
  #            map_signif_level=T, test='wilcox.test')+
  #facet_wrap(~Dataset, scales = 'free')+
  theme_bw()+ylab('Density')+xlab('Coefficient')+
  theme( axis.text.x = element_text(color='black', 
                                    size=14),
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

#Plotting the mean number of toxins per dataset
#Read the data with toxins' abundance
num_tox_per_dataset <- read.csv('summary_of_cry_nums_per_dataset.csv', header = T, sep = "\t", stringsAsFactors = F)

#Calculate median estimates
num_tox_per_dataset %>% group_by(Dataset) %>% summarize(tox=median(Cry_num))

#Assert levels for plotting
num_tox_per_dataset$Dataset=factor(num_tox_per_dataset$Dataset, 
                                   levels = c('Strains', 'Assemblies', 'Contigs', 'Plasmids', 'Chromosomes'))

#Vizualize the number of toxins per dataset
ggplot(num_tox_per_dataset, aes(fill=Dataset, x=Dataset, y=Cry_num))+
  geom_boxplot(alpha = 0.8,  color='black')+scale_fill_locuszoom()+
  theme_bw()+ylab('Number of Cry proteins')+xlab('Dataset')+
  theme( axis.text.x = element_text(color='black', 
                                    size=14),
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

#Compare the distribution of toxin types per dataset
#Read the data with the number of toxin types per dataset
tox_for_recs_and_refs <- read.csv('Tox_for_recs_and_refs.csv', header = T, sep = "\t", stringsAsFactors = F)

#Extract data for true majors
tox_for_recs_and_refs <- tox_for_recs_and_refs[tox_for_recs_and_refs$Group!='par_all',]

#Remove duplicates
tox_for_recs_and_refs <- unique(tox_for_recs_and_refs)

#Assert levels of the categorical variable
tox_for_recs_and_refs$Group <- factor(tox_for_recs_and_refs$Group, levels = c('rec', 'par_strict', 'ref_no_rec')) 

#Plot the number of toxins per group of strains and assemblies regarding toxin types they contain
ggplot(tox_for_recs_and_refs[tox_for_recs_and_refs$Mode=='soft',], aes(fill=Group, x=Group, y=Num))+
  geom_boxplot(alpha = 0.9,  color='black')+facet_wrap(~Dataset)+
  geom_signif(comparisons = list(c(1,2),
                                 c(2,3),
                                 c(1,3)),  y_position = c(14, 14.7,15.5),
              map_signif_level=T, test='wilcox.test')+
  scale_fill_manual(labels = c("Recombinants",'Parents','Unaffected'), values = c("#CD5922", '#9A73BE' ,'grey'))+
  theme_bw()+ylab('Number of Cry proteins')+xlab('Group')+
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
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))

#Extract data for the upset plot showing the number of strains/assemblies containing different toxin types
tox_for_upset <- tox_for_recs_and_refs[tox_for_recs_and_refs$Mode=='soft',]

#Get the names of strains/assemblies with certain toxin types
rec_list <- tox_for_upset[tox_for_upset$Group=='rec',2]
par_list <-  tox_for_upset[tox_for_upset$Group=='par_strict',2]
unaff_list <- tox_for_upset[tox_for_upset$Group=='ref_no_rec',2] 

#Get the number of assemblies/strains containing toxins of vertain types
tox_for_upset$Recombinants <- as.numeric(tox_for_upset$Accession %in% rec_list)
tox_for_upset$Parents <-  as.numeric(tox_for_upset$Accession %in% par_list)
tox_for_upset$Unaffected <-  as.numeric(tox_for_upset$Accession %in% unaff_list)

#Remove duplicates
tox_for_upset <- tox_for_upset[,-1]
tox_for_upset <- unique(tox_for_upset)

#Vizualize upset plot
upset(
  tox_for_upset[tox_for_upset$Dataset=='Assemblies',], #Choose type (strains/assemblis)
  c('Recombinants', 'Parents', 'Unaffected'), #Select columns 
  queries=list( #Assert colors
    upset_query(set='Recombinants', fill='#CD5922'),
    upset_query(set='Parents', fill='#9A73BE'),
    upset_query(set='Unaffected', fill='grey')
  ),
  base_annotations=list( #Set parameters of bars showing the number of intersections between selected groups
    'Intersection size'=(
      intersection_size(
        bar_number_threshold=3,  
        width=0.7,
        color='black',
        fill='lightgrey',
        alpha=.8) 
      + scale_y_continuous(expand=expansion(mult=c(0, 0.05)))
      + theme(
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title  =  element_text(size=18),
        axis.text.y = element_text(size=12),
        axis.line=element_line(colour='black')
      )
    )
  ),
  stripes=upset_stripes(
    geom=geom_segment(size=1),  
    colors=c('grey95', 'white')
  ),
  matrix=intersection_matrix( #Set parameters of intersections as colored circles
    geom=geom_point(
      shape='circle filled',
      size=4.5,
      stroke=0.45
    )
  ),
  set_sizes=( #Set parameters for bar shwong total number of compared sets
    upset_set_size(geom=geom_bar(width=0.4, color='black'))
    + theme(
      axis.line.x=element_line(colour='black'),
      axis.ticks.x=element_line(),
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12)
    )
  ),
  sort_sets='descending',
  sort_intersections='ascending'
)

#Read the data with the number of toxin types per dataset 
type_toxins_per_dataset <- read.csv('Num_tox_types_per_dataset_with_both.csv', header = T, sep = "\t", stringsAsFactors = F)

#Extract accessions (strains/assemblies, etc.) containing certain type of toxins
rec_acc <- unique(type_toxins_per_dataset[type_toxins_per_dataset$Tox_type=='rec',1])
par_strict_acc <- unique(type_toxins_per_dataset[type_toxins_per_dataset$Tox_type=='par_strict',1])
unaff_acc <- unique(type_toxins_per_dataset[type_toxins_per_dataset$Tox_type=='ref_no_rec',1])
both_acc <- unique(type_toxins_per_dataset[type_toxins_per_dataset$Tox_type=='both',1])

#Extract rows for true majors
type_toxins_strict <- type_toxins_per_dataset[type_toxins_per_dataset$Tox_type!='par_all',]

#Assert factor levels for plotting
type_toxins_strict$Tox_type <- factor(type_toxins_strict$Tox_type, levels = c('rec', 'par_strict', 'both','ref_no_rec')) #par_strict

#Plot the abundance of the toxins of certain types in the dataset
type_toxins_strict %>% dplyr::group_by( Dataset, Tox_type) %>% 
#type_toxins_strict[type_toxins_strict$Accession %in% rec_acc ,] %>% group_by( Dataset, Tox_type) %>% 
  dplyr::summarize(num_tox=n()) %>% as.data.frame() %>% ggplot(aes(x=reorder(Dataset,-num_tox), y=num_tox, fill=Tox_type)) +
  geom_bar(stat='identity',alpha = 0.7, color='black')+ #, position = 'fill'
  scale_fill_manual(labels = c("Recombinants",'Parents','Both','Unaffected'), values = c("#CD5922", '#9A73BE' , '#89D5D2','grey'))+
  theme_bw()+ylab('Number of Cry proteins')+xlab('Group')+
  theme( axis.text.x = element_text(color='black', 
                                    size=14),
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
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))

#Summarizing the number of mobile genetic elements in assemblies containing parents and recombininants simultaneously
#Create the dataframe with the number of MGEs per event and parent type
Num_MGE_for_events_plasmids = data.frame(Par_type = c('Min','Min','Min','Min','All','All','Min','Min','Min','Min','Min','Min','Min','Min','Min','Min','Median','Median','Median','Median','Median','Median'),
                                         Event = c('8','8','29','29','34','34','52.1','52.1','52.2','52.2','52.3','52.3','73','73','104','104','all','all','cry','cry','rec','rec'),
                                         MGE_type = c('GIs','IS','GIs','IS','GIs','IS','GIs','IS','GIs','IS','GIs','IS','GIs','IS','GIs','IS','GIs','IS','GIs','IS','GIs','IS'),
                                         MGE_num = c(20,85,20,85,7,67,7,85,8,86,10,84,23,307,20,260,9.5,37,12,81,13,85))

#Assert levels of parent types
Num_MGE_for_events_plasmids$Par_type <- factor (Num_MGE_for_events_plasmids$Par_type, levels = c('Min','All','Median'))

#Assert levels of recombination event IDs
Num_MGE_for_events_plasmids$Event <- factor(Num_MGE_for_events_plasmids$Event,
                                            levels = c('all','cry','rec', '8','29', '73', '52.1', '52.2','52.3','104','34'))

#Plot the number of MGEs
Num_MGE_for_events_plasmids %>% group_by( MGE_type) %>%
  ggplot(aes(x=Event, y=MGE_num, fill=Par_type))+
  facet_wrap(~MGE_type, scales ='free')+
  geom_bar(stat='identity',alpha = 0.8 , width=0.9)+
  scale_fill_manual( values = c( "#4F499E",'#BD5B0B','darkgrey'),
                     labels = c('Minor','All','Median'))+
  theme_bw()+ylab('Number of MGEs')+xlab('Type')+
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

#Analyzing the number of affected hosts for strains/assemblies cotaining certain toxin types
#Read the data with the number of hosts per dataset
num_hosts_strains_types <-  read.csv('num_hosts_per_dataset_in_toxin_group.csv', header = T, sep = "\t", stringsAsFactors = F)

#Exclude observations with absent data
num_hosts_strains_types <- num_hosts_strains_types[num_hosts_strains_types$Num_hosts>0,]

#Extract data for true majors only
num_hosts_strains_strict <- num_hosts_strains_types[num_hosts_strains_types$Tox_type!='par_all',]

#Assert levels of the categorical vriable
num_hosts_strains_strict$Tox_type <- factor(num_hosts_strains_strict$Tox_type, levels = c('rec', 'par_strict', 'ref_no_rec')) #par_strict

#Summarize median numbers of affected hosts per group
num_hosts_strains_strict[num_hosts_strains_strict$Host_type=='Species' & num_hosts_strains_strict$Group=='Rank3', ] %>% 
  group_by(Mode, Trim_flag, Tox_type) %>% summarize(median_hosts = median(Num_hosts))

#Plot the number of affected hosts per group
ggplot(num_hosts_strains_strict[num_hosts_strains_strict$Host_type=='Species' & num_hosts_strains_strict$Mode=='Assemblies' & num_hosts_strains_strict$Trim_flag=='Untrimmed',], 
       aes(fill=Tox_type, y=Num_hosts, x=Tox_type))+
  facet_wrap(~Group, scales ='free')+
  geom_violin(width=1.0, alpha = 0.7, trim = T)+
  geom_boxplot(width=0.4, alpha = 0.8)+
  geom_signif(comparisons = list(c("rec","par_strict"), 
                                 c('rec', 'ref_no_rec'),
                                 c('par_strict', 'ref_no_rec')),
              y_position = c(81.3, 88.7,83.5), 
              map_signif_level=T, test='wilcox.test')+
  
  scale_fill_manual(labels = c("Recombinants",'Parents','Unaffected'), values = c("#CD5922", '#9A73BE' ,'grey'))+
  theme_bw()+ylab('Number of affected hosts')+xlab('Group')+
  theme( axis.text.x = element_text(color='black', 
                                    size=14),
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
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))+
  guides(fill= guide_legend(title="Toxin type"))

#Vizualizing the properties of strains graph
#Read the data with properties per graph component
hosts_fill_per_graph_component <- read.csv('host_num_per_components.csv', header = T, sep = "\t", stringsAsFactors = F)

#Vizualize the distribution of the component size (the number of strains)
hosts_fill_per_graph_component[hosts_fill_per_graph_component$Comp_type=='All', c(1,5)] %>% unique() %>% 
  group_by(Num_strains) %>% tally() %>% as.data.frame() %>% 
  ggplot(aes(x=Num_strains, y=n)) +
  geom_bar(stat='identity',alpha = 0.9, fill='#2E98CA', width = 1)+
  theme_bw()+ylab('Frequency')+xlab('Number of strains')+
  theme( axis.text.x = element_text(color='black', 
                                    size=14),
         #axis.ticks.x = element_blank(),
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

#Plot the dependance between the number of toxins and recombination events within the component
hosts_fill_per_graph_component[hosts_fill_per_graph_component$Comp_type=='Rec', c(1,4,5,7,8,9)] %>%
  group_by( Comp_type)  %>% unique() %>% as.data.frame() %>% 
  ggplot(aes(x=Num_events, y=Num_tox_types, size=log10(Num_strains+1)))+
  geom_point(shape=21, color='black', alpha=0.8, fill ='#E68B24')+
  geom_smooth(method='lm', color = 'black',linetype = "dashed", show.legend = FALSE)+
  ylab('Number of toxins')+xlab('Number of events')

#Extract data for building linear model when comparing the number of toxins and events in the component
for_lm <- hosts_fill_per_graph_component[hosts_fill_per_graph_component$Comp_type=='Rec', c(1,4,5,7,8,9)] %>%
  group_by( Comp_type)  %>% unique() %>% as.data.frame()

#Run linear regression and show significance level
lm_num <- lm(Num_events ~ Num_tox_types, 
             data =for_lm)
summary(lm_num)

#Reading the data with the overlap between the sets of toxins in components containing parents and recombinants
graph_jac_components <- read.csv('Graph_components_per_event_rank3_long.csv', header = T, sep = "\t", stringsAsFactors = F)

#Calculate mean estimates
graph_jac_components %>% group_by(Type) %>% summarise(mean_jac = mean(Coeff))

#Calculate he numbwe of non-zero observations (at least oe common toxin)
nrow(graph_jac_components[graph_jac_components$Coeff==0 & graph_jac_components$Type=='rec_min',])
nrow(graph_jac_components[graph_jac_components$Coeff==0 & graph_jac_components$Type=='rec_maj',])
nrow(graph_jac_components[graph_jac_components$Coeff==0 & graph_jac_components$Type=='par_inter',])

#Assert factor levels for plotting
graph_jac_components$Type <- factor(graph_jac_components$Type, levels = c('rec_min', 'rec_maj', 'rec_all','par_inter'))

#Plot the distribution of similarities between components in terms of toxins' composition
ggplot(graph_jac_components, 
       aes(fill=Type, x=Coeff))+
  geom_histogram(alpha = 0.8, position = position_dodge2(width = 0.9, preserve = "single"), bins = 8, color='black')+
  scale_fill_manual( values = c( "#4F499E",'#89A05C','#BD5B0B','#C380C0'),
                     labels = c('Minor','Major','All','Parents'))+
  
  theme_bw()+ylab('frequency')+xlab('Coefficient')+
  theme( axis.text.x = element_text(color='black', 
                                    size=14),
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

#Extract columns with the number of toxins and strains
graph_for_comp_stat <- hosts_fill_per_graph_component[, c(1,4,5,7,8,9)] 

#Get the list of components containing recombinant toxins
rec_comps <- hosts_fill_per_graph_component[hosts_fill_per_graph_component$Comp_type=='Rec',1] %>% unique()

#Extract data for components containing recombinantns
graph_for_comp_stat <- graph_for_comp_stat[!(graph_for_comp_stat$Comp_type=='All' & graph_for_comp_stat$Component %in% rec_comps),]

#Plot the number of toxins/strains/affected hosts in components cotaining recombinants and devoid of them
graph_for_comp_stat[graph_for_comp_stat$Num_hosts>0, ] %>% 
  group_by( Comp_type)  %>% unique() %>% as.data.frame() %>% 
  ggplot(aes(x=Comp_type,y=Num_hosts , fill=Comp_type))+ #y=log10(Num_strains+1)
  geom_violin(width=0.7, alpha = 0.7, trim = T)+
  geom_boxplot(width=0.4, alpha = 0.8)+
  
  geom_signif(comparisons = list(c("Rec","All")), 
              map_signif_level=T, test='wilcox.test')+
  scale_fill_manual(values = c("grey", '#E68B24'),
                    labels = c('Unaffected', 'Affected'))+
  theme_bw()+
  ylab('Number of hosts')+xlab('Component type')+
  #ylab('log10(Number of strains)')+xlab('Component type')+
  theme( axis.text.x = element_text(color='black', 
                                    size=14),
         #axis.ticks.x = element_blank(),
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
#guides(size= guide_legend(title="log10(Number of Strains)"))

#Assert the level of the categorial variable for plotting
hosts_fill_per_graph_component$Host_type <- factor(hosts_fill_per_graph_component$Host_type,
                                                   levels = c('Lepidoptera', 'Diptera','Coleoptera','Nematoda',
                                                              'Hemiptera','Human','Hymenoptera', 'Orthoptera','Unknown'))
#Remove duplicate rows
graph_comp_unaffected <- hosts_fill_per_graph_component[!(hosts_fill_per_graph_component$Comp_type=='All' & hosts_fill_per_graph_component$Component %in% rec_comps),]

#Vizualize the distribution of affected hosts per component
graph_comp_unaffected %>% group_by( Comp_type, Component) %>%
  ggplot(aes(x=reorder(Component,-Coeff), y=Num_sp, fill=Host_type)) +
  facet_wrap(~Comp_type, scales ='free')+
  geom_bar(stat='identity',alpha = 0.8,  position = 'fill', width=1)+ #, color='black'
  scale_fill_manual(values = c("#00468BFF", '#ED0000FF' , '#42B540FF','#0099B4FF','#925E9FFF',
                               '#FDAF91FF','#AD002AFF','#ADB6B6FF','#ECF4F6'))+
  theme_bw()+ylab('Percentge of hosts')+xlab('Component')+
  theme( axis.text.x =element_blank(),
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
         legend.title=element_text(face="bold",size=18), 
         legend.text=element_text(size=16))

#Perform Fishers' extact test when comparing the percentage of unknown hots in graph components
graph_comp_unaffected$host_flag =as.factor(as.numeric(graph_comp_unaffected$Num_hosts>0))
graph_comp_unaffected[,c(1,4,10)] %>% unique() %>%  select(Comp_type, host_flag) %>% table() %>% fisher.test()

#Read the table with the distribution of affected hosts for components containing parents and recombinants
components_hosts_heatmap <- read.csv('Components_heatmap_only_annotated.csv', header = T, sep = "\t", stringsAsFactors = F)

#Order the components according to event IDs
components_hosts_heatmap <- components_hosts_heatmap %>% arrange(Event)

#Set numeric factors for correct order when plotting
components_hosts_heatmap <- components_hosts_heatmap %>% 
  mutate(Component_num =as.numeric(as.factor(Component)))%>% 
  mutate(event_num = as.numeric(as.factor(Event)))

#Vizualize the distribution of affected hosts for selected components with parents and recombinants
ggplot() + geom_scatterpie(aes(x=Component_num, y=event_num, r=0.5), data=components_hosts_heatmap, 
                           cols=c('Lepidoptera','Diptera','Coleoptera','Nematoda','Hemiptera','Human','Hymenoptera','Orthoptera'),  alpha =.7)+
  coord_equal()+
  scale_fill_lancet()+
  scale_y_continuous(breaks=c(1:37), labels=c("119",'34','18','45',"56","62",'78','84','87',
                                              '96','107','22','65','71','104','109','113','13',
                                              "14",'17',"19",'2','23','24','27','29','37','39',
                                              "44",'49','51','52',"53",'6','73', "8",'98'))+
  scale_x_continuous(breaks=c(1:43), 
                     labels=as.character(unique(components_hosts_heatmap$Component)[order(unique(components_hosts_heatmap$Component))]))+
  theme_bw()+ylab('Recombination ID')+xlab('Component ID')+
  theme(axis.text.x = element_text(color='black', 
                                   size=12, angle = 60, hjust = 1),
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
        legend.title=element_text(face="bold",size=18), 
        legend.text=element_text(size=16))

#Draw the adjacent subpanel with transferred domains
ggplot() + 
  geom_tile(aes(x='Domain', y=Event, 
                fill = unlist(lapply(Event, function (x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][1]))), 
            alpha =0.8, data=unique(components_hosts_heatmap[,c(1,2)]))+
  scale_fill_manual( values=c("#a60b0b", '#2980b9' ,'#bdbd00'))+ 
  theme_bw()+
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
         legend.text=element_text(size=14))+
  ylab('Event ID') +
  xlab('Domain')+
  coord_fixed()+
  guides(fill=guide_legend(override.aes = list(size=7), title="Domain"))
  
#Read the table with overlaps in affected hosts for components containing parents and recombinants
graph_num_overlap <- read.csv('Overlap_species_components.csv', header = T, sep = "\t", stringsAsFactors = F)

#Extract data for affected orders and species
graph_num_overlap_species <-graph_num_overlap[graph_num_overlap$Host_type=='Species',]
graph_num_overlap_orders <- graph_num_overlap[graph_num_overlap$Host_type=='Orders',]

#Vizualize the distribution of the overlap coefficients regarding the composition of affected hosts
ggplot(graph_num_overlap_species[graph_num_overlap_species$Dom_type!='domain2',], aes(fill=Tox_type, x=Coeff))+
  facet_wrap(~Coeff_type+Dom_type, scales = 'free_y')+
  scale_fill_manual( values = c('#89A05C',"#4F499E"), labels = c('Major',"Minor"))+
  geom_density(alpha = 0.8)+
  #geom_histogram(alpha = 0.8, position = position_dodge2(width = 0.9, preserve = "single"), bins = 5, color='black')
  theme_bw()+ylab('Density')+xlab('Overlap Coefficient')+
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

#Extract data for orders and species using difference similarity measures
graph_num_overlap_species_jac <-graph_num_overlap_species[graph_num_overlap_species$Coeff_type=='Jaccard',]
graph_num_overlap_species_sim <-graph_num_overlap_species[graph_num_overlap_species$Coeff_type=='Simpson',]

graph_num_overlap_orders_jac <- graph_num_overlap_orders[graph_num_overlap_orders$Coeff_type=='Jaccard',]
graph_num_overlap_orders_sim <- graph_num_overlap_orders[graph_num_overlap_orders$Coeff_type=='Simpson',]

#Make linear models to reveal the dependance between similarity coefficient and parent type
lm_species_jac <- lm(Coeff ~ Tox_type, 
                     data =graph_num_overlap_species_jac)
lm_species_simp <- lm(Coeff ~ Tox_type, 
                      data =graph_num_overlap_species_sim)

lm_orders_jac <- lm(Coeff ~ Tox_type, 
                    data =graph_num_overlap_orders_jac)
lm_orders_simp <- lm(Coeff ~ Tox_type, 
                     data =graph_num_overlap_orders_sim)

#Show significance estimates
summary(lm_species_jac)
summary(lm_species_simp)
summary(lm_orders_jac)
summary(lm_orders_simp)