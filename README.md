# Araruama-metagenomes
Supplementary scripts used in the manuscript 
##############################################################################################
#
#                            Filter MG-RAST annotation file (phylum level)
#
##############################################################################################
library("reshape")
library("dplyr")

data <- read.delim("adicional_samples/all.mats_tax.csv", sep = ";")

head(data)

# Get only Archaea and Bacteria hits

dataf <- filter(data, domain %in% c("Bacteria","Archaea"))

dataf <- data.frame(dataf)

head(dataf)

phyla.m <- melt(dataf)

head(phyla.m)

df <- phyla.m

df2<-cast(df, variable~phylum, fun.aggregate=sum)

head(df2)

df2.norm<-cbind(sample=df2$variable, df2[2:ncol(df2)]/rowSums(df2[2:ncol(df2)]))

head(df2.norm)

df2.norm$sample<-as.character(df2.norm$sample)

df2.norm<-df2.norm[order(df2.norm$sample),]

df2.norm<-data.frame(df2.norm)

View(df2.norm)

metadata <- read.csv("adicional_samples/mats_metadata_new_jw_pmm.csv", check.names = FALSE, sep = ";")

metadata

metadata$sample<-as.character(metadata$sample)

metadata<-metadata[order(metadata$sample),]

View(metadata)

data_complete<-cbind(sample=metadata$sample, Location=metadata$Location, Sample=metadata$new_names, status=metadata$Status, df2.norm[,2:ncol(df2.norm)])

data_filtered <- filter(data_complete, status %in% c("Keep"))

data_final<-cbind(sample=data_filtered$Sample, Location=data_filtered$Location, data_filtered[,5:ncol(data_filtered)])

write.csv(data_final, file = "taxonomic_phyla_converted.csv", row.names = FALSE)



##############################################################################################
#
#                            Filter MG-RAST annotation file (genera level)
#
##############################################################################################
library("reshape")
library("dplyr")

data <- read.delim("adicional_samples/all.mats_tax.csv", sep = ";")

head(data)

# Get only Archaea and Bacteria hits

dataf <- filter(data, domain %in% c("Bacteria","Archaea"))

dataf <- data.frame(dataf)

head(dataf)

phyla.m <- melt(dataf)

head(phyla.m)

df <- phyla.m

df2<-cast(df, variable~genus, fun.aggregate=sum)

head(df2)

df2.norm<-cbind(sample=df2$variable, df2[2:ncol(df2)]/rowSums(df2[2:ncol(df2)]))

head(df2.norm)

df2.norm$sample<-as.character(df2.norm$sample)

df2.norm<-df2.norm[order(df2.norm$sample),]

df2.norm<-data.frame(df2.norm)

View(df2.norm)

metadata <- read.csv("adicional_samples/mats_metadata_new_jw_pmm.csv", check.names = FALSE, sep = ";")

metadata

metadata$sample<-as.character(metadata$sample)

metadata<-metadata[order(metadata$sample),]

View(metadata)

data_complete<-cbind(sample=metadata$sample, Location=metadata$Location, Sample=metadata$new_names, status=metadata$Status, df2.norm[,2:ncol(df2.norm)])

data_filtered <- filter(data_complete, status %in% c("Keep"))

data_final<-cbind(sample=data_filtered$Sample, Location=data_filtered$Location, data_filtered[,5:ncol(data_filtered)])

write.csv(data_final, file = "taxonomic_genera_converted.csv", row.names = FALSE)

##############################################################################################
#
#                                      Barplot phyla
#
##############################################################################################
library("ggdendro")
library(cowplot)
library(ggplot2)
library("reshape")
library("dplyr")

data <- read.csv("taxonomic_phyla_converted.csv", check.names=FALSE)

data <-data.frame(data, row.names = data$sample)

#make the hierarchical clustering
dd <- dist(scale(data[,3:ncol(data)]), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

dend <- as.dendrogram(hc)
dend_data <- dendro_data(dend, type = "rectangle")
labels<-dend_data$labels$label
labels <- as.character(labels)





p1<-ggplot(segment(dend_data),labels=labels)+ 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))+ 
  theme_dendro()+coord_flip()+scale_y_reverse()


data.m<-melt(data)

abundant<-data.m %>% group_by(variable) %>% 
  summarise(media = mean(value)) %>% 
  arrange(desc(media)) %>% 
  ungroup() %>% 
  top_n(11) %>% pull(variable) %>% as.character()

data_plot <- bind_cols(data %>% select(sample, Location, abundant))
others <- data.frame(others=1-round(rowSums(data_plot[,3:ncol(data_plot)]), digits = 3))
data_plot<-cbind(data_plot,others)



data.m <- melt(data_plot)
data.m$sample<-factor(data.m$sample,levels = labels)
data.m$variable <- factor(data.m$variable, levels = rev(levels(data.m$variable)))


#make the barplot
p2<-ggplot(data.m, aes(x=sample, y=value, fill=variable, order = -as.numeric(variable))) + 
  labs(list(title = NULL , x = NULL, y = "Relative Abundance (%)\n", fill=NULL)) + 
  geom_bar(position="fill", stat="identity")+ 
  guides(fill=NULL)+ 
  scale_fill_brewer(palette="Paired")+
  coord_flip()+
  ggtitle(NULL) +
  xlab(NULL) + ylab("Frequency")+
  theme_bw() +
  theme(panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  theme(axis.text.x = element_text(size=12, angle=0, hjust = 1, vjust=0.99))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))+
  theme(axis.line = element_line())+
  theme(panel.border = element_blank())+
  theme(axis.ticks = element_blank())+
  theme( # remove the vertical grid lines
    panel.grid.major.x = element_line( size=.1, color="darkgrey") ,
    panel.grid.major.y = element_line( size=.1, color="darkgrey"),
    panel.grid.minor.y = element_blank())+
  theme(panel.grid.minor = element_blank())+ 
  theme(panel.background = element_blank())+
  theme(axis.line.x=element_blank())+
  theme(axis.line.y=element_blank())


p2
#combine both in a panel
panel<-plot_grid(p1, p2, align = "h")
panel

ggsave(file="Phyla_w_cluster.png", width=15, height=10)

ggsave(file="Phyla_w_cluster.svg", width=15, height=12)


##############################################################################################
#
#                                      Barplot genera
#
##############################################################################################
library("ggdendro")
library(cowplot)
library(ggplot2)
library("reshape")
library("dplyr")

data <- read.csv("taxonomic_genera_converted.csv", check.names=FALSE)

data <-data.frame(data, row.names = data$sample)

#make the hierarchical clustering
dd <- dist(scale(data[,3:ncol(data)]), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

dend <- as.dendrogram(hc)
dend_data <- dendro_data(dend, type = "rectangle")
labels<-dend_data$labels$label
labels <- as.character(labels)





p1<-ggplot(segment(dend_data),labels=labels)+ 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))+ 
  theme_dendro()+coord_flip()+scale_y_reverse()


data.m<-melt(data)

abundant<-data.m %>% group_by(variable) %>% 
  summarise(media = mean(value)) %>% 
  arrange(desc(media)) %>% 
  ungroup() %>% 
  top_n(11) %>% pull(variable) %>% as.character()

data_plot <- bind_cols(data %>% select(sample, Location, abundant))
others <- data.frame(others=1-round(rowSums(data_plot[,3:ncol(data_plot)]), digits = 3))
data_plot<-cbind(data_plot,others)



data.m <- melt(data_plot)
data.m$sample<-factor(data.m$sample,levels = labels)
data.m$variable <- factor(data.m$variable, levels = rev(levels(data.m$variable)))


#make the barplot
p2<-ggplot(data.m, aes(x=sample, y=value, fill=variable, order = -as.numeric(variable))) + 
  labs(list(title = NULL , x = NULL, y = "Relative Abundance (%)\n", fill=NULL)) + 
  geom_bar(position="fill", stat="identity")+ 
  guides(fill=NULL)+ 
  scale_fill_brewer(palette="Paired")+
  coord_flip()+
  ggtitle(NULL) +
  xlab(NULL) + ylab("Frequency")+
  theme_bw() +
  theme(panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())+
  theme(axis.text.x = element_text(size=12, angle=0, hjust = 1, vjust=0.99))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))+
  theme(axis.line = element_line())+
  theme(panel.border = element_blank())+
  theme(axis.ticks = element_blank())+
  theme( # remove the vertical grid lines
    panel.grid.major.x = element_line( size=.1, color="darkgrey") ,
    panel.grid.major.y = element_line( size=.1, color="darkgrey"),
    panel.grid.minor.y = element_blank())+
  theme(panel.grid.minor = element_blank())+ 
  theme(panel.background = element_blank())+
  theme(axis.line.x=element_blank())+
  theme(axis.line.y=element_blank())


p2
#combine both in a panel
panel<-plot_grid(p1, p2, align = "h")
panel

ggsave(file="genera_w_cluster.png", width=15, height=10)

ggsave(file="genera_w_cluster.svg", width=15, height=12)


##############################################################################################
#
#                        Anosim - Taxonomic Annotation (Phylum level)
#
##############################################################################################
library(vegan)

data <- read.csv("taxonomic_phyla_converted.csv", check.names=FALSE)
data<-data.frame(data, check.names=FALSE)


View(data)
levels(data$Location)

data_num<-data[3:ncol(data)]

attach(data)
attach(data_num)


data.dist <- vegdist(data_num)
data.anosim <- anosim(data.dist, data$Location)

data.anosim$signif
data.anosim$statistic

summary(data.anosim)
plot(data.anosim)


##############################################################################################
#
#                                             nMDS 
#
##############################################################################################
library(vegan)
library(ggplot2)


trans.arcsine <- function(x){
  asin(sqrt(abs(x)))
}

#To run the nMDS for the functional profile, just replace the input file in the following line
data <- read.csv("taxonomic_phyla_converted.csv", check.names=FALSE, sep = ",")
head(data)

data <- data %>%
  mutate_if(is.factor, as.numeric)


data<-data.frame(data, check.names=FALSE)
attach(data)
data_num<-data[3:ncol(data)]
data_num<-trans.arcsine(data_num)


attach(data_num)

ord <- metaMDS(data_num)
ord_df <- cbind(sample,Location,data.frame(ord$points))
ord_df$Location <- factor(ord_df$Location)

attach(ord_df)

mdspoints <- data.frame(scores(ord))
mdsvectors <- data.frame(ord$species)
mdsvectors$species <- row.names(mdsvectors)
#mdsvectors <- filter(mdsvectors, species %in% target)

ggplot(ord_df, aes(x=MDS1, y=MDS2, color=Location)) + 
  geom_point(size = 4) +
  scale_color_brewer(palette="Paired")+
  geom_segment(data = mdsvectors, aes(x=0, xend=MDS1, y=0, yend=MDS2),
               arrow = arrow(length = unit(0.1, "cm")),
               colour="grey") +
  geom_text(data=mdsvectors, aes(x=MDS1, y=MDS2, label=species), size=2, colour="black") +
    theme_bw()+
  theme(panel.border = element_rect(colour = "black", size=1))+
  theme(legend.key=element_rect(fill='white'))+
  theme(legend.key = element_rect(colour = "white"))+
  #scale_colour_manual(values=c("#43CD80","#9400D3","steelblue3","#DEB887"))+ 
  theme(legend.text=element_text(size=12))+
  theme(axis.title.x = element_text(face="bold", size=14))+
  theme(axis.title.y = element_text(face="bold", size=14))+
  theme(axis.text.x = element_text(size=14,color="black"))+
  theme(axis.text.y = element_text(size=14,color="black"))+
  theme(axis.line.x=element_blank())+
  theme(axis.line.y=element_blank())+
  annotate("text", x = max(ord_df$MDS1)-0.3, y = max(ord_df$MDS2)-0.01, label = paste0("Stress=",round(ord$stress, digits=3),"\n","ANOSIM: P=",round(data.anosim$signif, digits=3),", R=",round(data.anosim$statistic,digits = 3)))+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.minor.x = element_blank() ,
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

#To sabe the nMDS biplots for the functional profile, just replace the output file name in the following lines
ggsave(file="output_nMDS__phyla_taxonomic.jpg", width=10, height=8)
ggsave(file="output_nMDS__phyla_taxonomic.svg", width=10, height=8)
