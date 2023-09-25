library(ape)
library(ggtree)
library(dplyr)
library(phangorn)
library(dendextend)
library(CollessLike)
library("DECIPHER")

setwd("../data_for_scripts/data_for_scripts/for_viz/phylogeny")

#Vizualization of trees
#Read domain-wise and concatenated trees
full_tree <-read.tree("ML_Full_nuc.raxml.support.nwk")
dom1_tree <- read.tree("ML_Dom1.raxml.support")
dom2_tree <- read.tree("ML_Dom2.raxml.support")
dom3_tree <- read.tree("ML_Dom3.raxml.support")

#Assert labels to group tips accorcing to their presence in the BPPRC
groupInfo <- split(full_tree$tip.label, 
                   gsub("([A-Z]{2,3})_([a-zA-Z]+[0-9]+[A-Za-z]+)(.*)", '\\1', full_tree$tip.label))
tr_gr <- groupOTU(full_tree, groupInfo)

#Plot coloured trees
#Concatenated treee
full_plot <-
  ggtree(tr_gr, aes(color=group), layout='circular',branch.length='none',alpha =1)+
  scale_color_manual(values = c('#5F9EA0', "#DC143C"))+
  geom_tiplab(size=1, aes(angle=angle), show.legend=T, color='black')

#Domain 1
groupInfo <- split(dom1_tree$tip.label, 
                   gsub("([A-Z]{2,3})_([a-zA-Z]+[0-9]+[A-Za-z]+)(.*)", '\\1', dom1_tree$tip.label))
tr_gr <- groupOTU(dom1_tree, groupInfo)
ggtree(tr_gr, aes(color=group), layout='circular',branch.length='none',alpha =1)+
  scale_color_manual(values = c("#DC143C",'#5F9EA0'))+
  geom_tiplab(size=1, aes(angle=angle), show.legend=T, color='black')

#Domain 2
groupInfo <- split(dom2_tree$tip.label, 
                   gsub("([A-Z]{2,3})_([a-zA-Z]+[0-9]+[A-Za-z]+)(.*)", '\\1', dom2_tree$tip.label))
tr_gr <- groupOTU(dom2_tree, groupInfo)
ggtree(tr_gr, aes(color=group), layout='circular',branch.length='none',alpha =1)+
  scale_color_manual(values = c("#DC143C",'#5F9EA0'))+
  geom_tiplab(size=1, aes(angle=angle), show.legend=T, color='black')

#Domain 3
groupInfo <- split(dom3_tree$tip.label, gsub("([A-Z]{2,3})_([a-zA-Z]+[0-9]+[A-Za-z]+)(.*)", '\\1', full_tree$tip.label))
tr_gr <- groupOTU(dom3_tree, groupInfo)
ggtree(tr_gr, aes(color=group), layout='circular',branch.length='none',alpha =1)+
  scale_color_manual(values = c("#DC143C",'#5F9EA0'))+
  geom_tiplab(size=1, aes(angle=angle), show.legend=T, color='black')

#Comparision of tree topologies
#Patristic distance
cor_cophenetic(dom1_tree, dom2_tree) #0.7673982
cor_cophenetic(dom1_tree, dom3_tree) #0.6663316
cor_cophenetic(dom2_tree, dom3_tree) #0.6918381
cor_cophenetic(full_tree, dom1_tree) #0.8940708
cor_cophenetic(full_tree, dom2_tree) #0.8758706
cor_cophenetic(full_tree, dom3_tree) #0.758726

#Cophenetic correlation
coph1_matr <- cophenetic.phylo(dom1_tree)
coph1_matr[upper.tri(coph1_matr)] <- 0

coph2_matr <- cophenetic.phylo(dom2_tree)
coph2_matr[upper.tri(coph2_matr)] <- 0

coph3_matr <- cophenetic.phylo(dom3_tree)
coph3_matr[upper.tri(coph3_matr)] <- 0

cophf_matr <- cophenetic.phylo(full_tree)
cophf_matr[upper.tri(cophf_matr)] <- 0

#Summarize the number of cophenetic comparisions in the matrix
sum(coph1_matr) #110671.3
sum(coph2_matr) #124245.1
sum(coph3_matr) #88018.16
sum(cophf_matr) #101705.4

#Assessment of tree quality and other properties
#Calculating tree balance
c(balance.indices(dom1_tree))
c(balance.indices(dom2_tree))
c(balance.indices(dom3_tree))
c(balance.indices(full_tree))

#Reading the alignment data for revealing homoplasy signals
dom1_seq <- read.phyDat('domain1.msa', format = "fasta", type = "DNA")
dom2_seq <- read.phyDat('domain2.msa', format = "fasta", type = "DNA")
dom3_seq <- read.phyDat('domain3.msa', format = "fasta", type = "DNA")
full_seq <- read.phyDat('All_domains.msa', format = "fasta", type = "DNA")

#Consistency index
CI(dom1_tree, dom1_seq) #0.08954948
CI(dom2_tree, dom2_seq) #0.09955026
CI(dom3_tree, dom3_seq) #.08901361
CI(full_tree, full_seq) #0.08669106

#Retention index
RI(dom1_tree, dom1_seq) #0.6053512
RI(dom2_tree, dom2_seq) #0.5951812
RI(dom3_tree, dom3_seq) #0.5707345
RI(full_tree, full_seq) #0.5570568

#Vizualizing differences in tree topologies using tanglegrams
#Read trees devoid of supporting values and tree lengths
dend1 <- ReadDendrogram(file="no_support/ML_Dom1.raxml.bestTree") 
dend2 <- ReadDendrogram(file="no_support/ML_Dom2.raxml.bestTree") 
dend3 <- ReadDendrogram(file="no_support/ML_Dom3.raxml.bestTree") 
dendf <- ReadDendrogram(file="no_support/ML_Full_nuc.raxml.bestTree") 

#Rearrange dendrograms in a ladderized fashion
dnd1 <- ladderize(sort(dend1) )
dnd2 <- ladderize(sort(dend2))
dnd3 <- ladderize(sort(dend3) )
dndf <- ladderize(sort(dendf))

#Make lists of distances to draw tanglegrams
dndlist1f <- dendextend::dendlist(dnd1, dndf)
dndlist2f <- dendextend::dendlist(dnd2, dndf)
dndlist3f <- dendextend::dendlist(dnd3, dndf)

#Plot tanglegrams
dndlist3f %>% tanglegram(common_subtrees_color_branches = TRUE,
                         lab.cex = .01, margin_inner = .02, edge.lwd=.8, 
                         margin_outer=12,margin_bottom=0.0000,lwd=.3, axes=F,sort=T,
                         type ='r',cex_main=2.4,
                         main_left ='Domain3', main_right='Full',
                         highlight_distinct_edges = F,
                         highlight_branches_lwd = F,
                         faster=F)

#Vizualize inconsistency between the domain-wise identity between toxins and the reference phylogeny
#Read identity heatmaps
heatmap_vals = read.csv('domain3_heatmap.csv', header = T, sep = "\t",row.names = 1,
                        stringsAsFactors = F)  #domain1_heatmap.csv, domain2_heatmap.csv

#Melt matrix to the dataframe with pairs
colnames(heatmap_vals)[46] <- 'BT_Cry1A-like'
melted_dist <- reshape2::melt(as.matrix(heatmap_vals), na.rm = TRUE)

#Save the plot of the full tree
testing_order <- melted_dist[melted_dist$Var1==melted_dist$Var2,]
p <- ggtree(full_tree,alpha =1) %<+% testing_order
p <-p+ scale_y_reverse()

#Get the order of toxins from the plot data
ref_plot_data <- p$data %>% as.data.frame()
ref_plot_data_filtered <- ref_plot_data %>% arrange(y) %>% select(-c(parent,node,branch,branch.length,isTip,x,y,angle)) %>% na.omit()

#Assign factor levels for ordering
melted_dist$Var1_new <- factor(melted_dist$Var1 , levels = ref_plot_data_filtered$label)
melted_dist$Var2_new <- factor(melted_dist$Var2 , levels = ref_plot_data_filtered$label)

#Plot ordered heatmap
ggplot(data = melted_dist, aes(Var1, Var2 , fill = value))+ #Var1, Var2 Var1_new, Var2_new
  geom_tile()+ 
  scale_fill_gradient(low = '#FDFDFD', high = "darkred", limit = c(min(melted_dist$value),100), space = "Lab") + #37 #low = "#FDFDFD" high = "#AA0C00"
  theme_bw() + xlab('Trees') + ylab('Trees') +
  theme( axis.text.y = element_blank(),
         axis.title.y=element_blank(),
         panel.background =  element_blank(),
         panel.grid.minor = element_blank(), 
         panel.grid.major =  element_blank(),
         legend.position = "none",
         axis.text.x = element_blank(),
         axis.title.x =  element_blank(),
         axis.ticks = element_blank())+
  coord_fixed()
