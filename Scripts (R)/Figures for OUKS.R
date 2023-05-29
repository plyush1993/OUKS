###########################################################################################
# Figures for OUKS
###########################################################################################

###########################################################################################
# spider plot
###########################################################################################

library(fmsb)

scores <- data.frame(
  "1.Randomization"  = c(191),
  "2.Integration" = c(794),
  "3.Imputation" = c(364),
  "4.Correction" = c(2553),
  "5.Annotation" = c(702),
  "6.Filtering" = c(468),
  "7.Normalization" = c(731),
  "8.Grouping" = c(150),
  "9.Statistics" = c(3656)
)


max_min <- data.frame(
  "1.Randomization" = c(3656, 0), "2.Integration" = c(3656, 0), "3.Imputation" = c(3656, 0),
  "4.Correction" = c(3656, 0), "5.Annotation" = c(3656, 0), "6.Filtering" = c(3656, 0),
  "7.Normalization" = c(3656, 0), "8.Grouping" = c(3656, 0), "9.Statistics" = c(3656, 0)
)

rownames(max_min) <- c("Max", "Min")

# Bind the variable ranges to the data
df <- rbind(max_min, scores)
df

data <- df[c("Max", "Min", "1"), ]
color <- "#0080bb"
colnames(data) <- c("1.Randomization", "2.Integration", "3.Imputation" , # for logo
                    "4.Correction", "5.Annotation", "6.Filtering",
                    "7.Normalization", "8.Grouping", "9.Statistics")



radarchart(data, pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1, cglcol = "grey", cglty = 1, cglwd = 0.8)
mtext(side = 3, line = -10, at = 0, cex = 2.0, "OUKS", col = '#aa1d35', padj = c(0,0), font = 11)
dev.off()

# jpeg("rplot.jpeg") pdf("rplot.pdf") tiff("rplot.tiff")
tiff("Spider.tiff", res = 75) # tiff("Graphical Abstract.tiff", res = 75)
radarchart(data, pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1, cglcol = "grey", cglty = 1, cglwd = 0.8)
mtext(side = 3, line = -13, at = 0, cex = 2.0, "OUKS", col = '#aa1d35', padj = c(0,0), font = 11)
dev.off()

pdf("Spider.pdf") # pdf("Graphical Abstract.pdf")
radarchart(data, pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1, cglcol = "grey", cglty = 1, cglwd = 0.8)
mtext(side = 3, line = -13, at = 0, cex = 2.0, "OUKS", col = '#aa1d35', padj = c(0,0), font = 11)
dev.off()

colnames(data) <- c("", "", "" ,"", "", "", "", "", "") # for Graphical Abstract

radarchart(data, pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1, cglcol = "grey", cglty = 1, cglwd = 0.8)
mtext(side = 3, line = -10, at = 0, cex = 2.0, "OUKS", col = '#aa1d35', padj = c(0,0), font = 11)
dev.off()

tiff("Graphical Abstract.tiff", res = 75) # tiff("Graphical Abstract.tiff", res = 75)
radarchart(data, pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1, cglcol = "grey", cglty = 1, cglwd = 0.8)
mtext(side = 3, line = -13, at = 0, cex = 2.0, "OUKS", col = '#aa1d35', padj = c(0,0), font = 11)
dev.off()

pdf("Graphical Abstract.pdf") # pdf("Graphical Abstract.pdf")
radarchart(data, pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1, cglcol = "grey", cglty = 1, cglwd = 0.8)
mtext(side = 3, line = -13, at = 0, cex = 2.0, "OUKS", col = '#aa1d35', padj = c(0,0), font = 11)
dev.off()

###########################################################################################
# pie chart
###########################################################################################

library(ggsci)

scores <- data.frame(
  "1.Randomization"  = c(191),
  "2.Integration" = c(794),
  "3.Imputation" = c(364),
  "4.Correction" = c(2553),
  "5.Annotation" = c(702),
  "6.Filtering" = c(468),
  "7.Normalization" = c(731),
  "8.Grouping" = c(150),
  "9.Statistics" = c(3656)
)

slices <- round(scores/max(scores)*100, 0)
lbls <- c("Randomization:", "Integration:", "Imputation:" , 
          "Correction:", "Annotation:", "Filtering:",
          "Normalization:", "Grouping:", "Statistics:")
pct <- round(scores/sum(scores)*100, 0)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(as.numeric(slices),labels = lbls, col=pal_gsea("default")(length(lbls)), # col=pal_gsea("default")(length(lbls)) or col=pal_igv("default")(length(lbls)) or pal_npg("nrc")(length(lbls))
    main="")

dev.off()

# jpeg("rplot.jpeg") pdf("rplot.pdf") tiff("rplot.tiff")
tiff("Pie.tiff", width = 680, height = 525, res = 100) # tiff("Graphical Abstract.tiff", res = 75)
pie(as.numeric(slices),labels = lbls, col=pal_gsea("default")(length(lbls)), # col=pal_gsea("default")(length(lbls)) or col=pal_igv("default")(length(lbls)) or pal_npg("nrc")(length(lbls))
    main="")
dev.off()

pdf("Pie.pdf", width = 7.08, height = 5.48) # pdf("Graphical Abstract.pdf")
pie(as.numeric(slices),labels = lbls, col=pal_gsea("default")(length(lbls)), # col=pal_gsea("default")(length(lbls)) or col=pal_igv("default")(length(lbls)) or pal_npg("nrc")(length(lbls))
    main="")
dev.off()

###########################################################################################
# Figure 1
###########################################################################################

library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
library(stringr)
setwd("D:/...")

######################
###### raw data
######################

dsr <-as.data.frame(fread(input = "xcms after IPO MVI.csv", header=T)) # first column with all metadata
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]
# metadata generating
rname <- rownames(dsr) 
rname <- str_remove(rname, ".CDF")
all_id <- sapply(1:length(rname), function(y) unlist(str_split(rname[y], " "))) 
ro_id <- as.numeric(unlist(lapply(all_id, function(y) unlist(y[1])))) 
b_id <- unlist(lapply(all_id, function(y) unlist(y[3]))) 
s_id <- unlist(lapply(all_id, function(y) unlist(y[2]))) 
p1_id <- unlist(lapply(all_id, function(y) unlist(y[4]))) 
p2_id <- unlist(lapply(all_id, function(y) unlist(y[5]))) 
p_all_id <- sapply(1:length(p1_id), function(y) ifelse(is.na(p2_id[y]), p1_id[y] ,paste(p1_id[y], p2_id[y], sep = " "))) 
un_raw_l <- unique(gsub("[[:digit:]]", "", p1_id)) 
true_l <- c("QC", "TG", "CG") 
rbind(un_raw_l,true_l) 
raw_l <- gsub("[[:digit:]]", "", p1_id) 
n_gr_t <- as.character(match(raw_l, un_raw_l)) 
for (i in 1:length(unique(n_gr_t))) {
  n_gr_t <- str_replace(n_gr_t, unique(n_gr_t)[i], true_l[i]) } 
s_id_p <- sapply(1:length(s_id), function(y) unlist(str_split(s_id[y], "_"))) 
s_id_p2 <- unlist(lapply(s_id_p, function(y) unlist(y[1]))) 
all_meta_data<- as.data.frame(cbind(n_gr_t, ro_id, b_id, s_id, p_all_id, raw_l, s_id_p2)) 
# join features and metadata and order
ds <- data.frame(cbind(all_meta_data, dsr))
ds$ro_id <- as.numeric(ds$ro_id)
ds <- ds[order(ds$ro_id, decreasing = F),] 

dat <- ds[,-c(1:7)]
ds_s <- sapply(1:nrow(dat), function(x) sum(dat[x,], na.rm = T))
ds$b_id <- as.numeric(stringr::str_remove(ds$b_id, "b"))
ds$n_gr_t[-which(stringr::str_detect(ds$n_gr_t, "QC"))] <- "Sample"
ds_plot <- as.data.frame(cbind(ds[,1:7], ds_s))
ds_plot <- subset(ds_plot, n_gr_t == "QC")
a <- ggplot(ds_plot, aes(x=ro_id, y=ds_s, color=b_id)) + 
  geom_point(size = 2) + theme_bw() + ggtitle("") +
  xlab("Run order") + ylab("log2(Total intensity)") + labs(color = "Batch") + theme(legend.position="top") +
  scale_color_gradient(low="blue", high="red", breaks = c(as.numeric(quantile(ds$b_id))))  +  geom_smooth(method=lm) # lm or loess

dat <- ds[,-c(1:7)]
ds_s <- sapply(1:nrow(dat), function(x) sum(dat[x,], na.rm = T))
ds$b_id <- as.numeric(stringr::str_remove(ds$b_id, "b"))
ds$n_gr_t[-which(stringr::str_detect(ds$n_gr_t, "QC"))] <- "Sample"
ds_plot <- as.data.frame(cbind(ds[,1:7], ds_s))
ds_plot <- subset(ds_plot, n_gr_t == "Sample")
b <- ggplot(ds_plot, aes(x=ro_id, y=ds_s, color=b_id)) + 
  geom_point(size = 2) + theme_bw() + ggtitle("") +
  xlab("Run order") + ylab("log2(Total intensity)") + labs(color = "Batch") + theme(legend.position="top") +
  scale_color_gradient(low="blue", high="red", breaks = c(as.numeric(quantile(ds$b_id))))  +  geom_smooth(method=lm) # lm or loess

######################
###### QC-XGB data
######################

dsr <-as.data.frame(fread(input = "xcms after IPO MVI QC-XGB.csv", header=T)) # first column with all metadata
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]
# metadata generating
rname <- rownames(dsr) 
rname <- str_remove(rname, ".CDF")
all_id <- sapply(1:length(rname), function(y) unlist(str_split(rname[y], " "))) 
ro_id <- as.numeric(unlist(lapply(all_id, function(y) unlist(y[1])))) 
b_id <- unlist(lapply(all_id, function(y) unlist(y[3]))) 
s_id <- unlist(lapply(all_id, function(y) unlist(y[2]))) 
p1_id <- unlist(lapply(all_id, function(y) unlist(y[4]))) 
p2_id <- unlist(lapply(all_id, function(y) unlist(y[5]))) 
p_all_id <- sapply(1:length(p1_id), function(y) ifelse(is.na(p2_id[y]), p1_id[y] ,paste(p1_id[y], p2_id[y], sep = " "))) 
un_raw_l <- unique(gsub("[[:digit:]]", "", p1_id)) 
true_l <- c("QC", "TG", "CG") 
rbind(un_raw_l,true_l) 
raw_l <- gsub("[[:digit:]]", "", p1_id) 
n_gr_t <- as.character(match(raw_l, un_raw_l)) 
for (i in 1:length(unique(n_gr_t))) {
  n_gr_t <- str_replace(n_gr_t, unique(n_gr_t)[i], true_l[i]) } 
s_id_p <- sapply(1:length(s_id), function(y) unlist(str_split(s_id[y], "_"))) 
s_id_p2 <- unlist(lapply(s_id_p, function(y) unlist(y[1]))) 
all_meta_data<- as.data.frame(cbind(n_gr_t, ro_id, b_id, s_id, p_all_id, raw_l, s_id_p2)) 
# join features and metadata and order
ds <- data.frame(cbind(all_meta_data, dsr))
ds$ro_id <- as.numeric(ds$ro_id)
ds <- ds[order(ds$ro_id, decreasing = F),] 

dat <- ds[,-c(1:7)]
ds_s <- sapply(1:nrow(dat), function(x) sum(dat[x,], na.rm = T))
ds$b_id <- as.numeric(stringr::str_remove(ds$b_id, "b"))
ds$n_gr_t[-which(stringr::str_detect(ds$n_gr_t, "QC"))] <- "Sample"
ds_plot <- as.data.frame(cbind(ds[,1:7], ds_s))
ds_plot <- subset(ds_plot, n_gr_t == "QC")
c <- ggplot(ds_plot, aes(x=ro_id, y=ds_s, color=b_id)) + 
  geom_point(size = 2) + theme_bw() + ggtitle("") +
  xlab("Run order") + ylab("log2(Total intensity)") + labs(color = "Batch") + theme(legend.position="top") +
  scale_color_gradient(low="blue", high="red", breaks = c(as.numeric(quantile(ds$b_id))))   +  geom_smooth(method=lm) # lm or loess

dat <- ds[,-c(1:7)]
ds_s <- sapply(1:nrow(dat), function(x) sum(dat[x,], na.rm = T))
ds$b_id <- as.numeric(stringr::str_remove(ds$b_id, "b"))
ds$n_gr_t[-which(stringr::str_detect(ds$n_gr_t, "QC"))] <- "Sample"
ds_plot <- as.data.frame(cbind(ds[,1:7], ds_s))
ds_plot <- subset(ds_plot, n_gr_t == "Sample")
d <- ggplot(ds_plot, aes(x=ro_id, y=ds_s, color=b_id)) + 
  geom_point(size = 2) + theme_bw() + ggtitle("") +
  xlab("Run order") + ylab("log2(Total intensity)") + labs(color = "Batch") + theme(legend.position="top") +
  scale_color_gradient(low="blue", high="red", breaks = c(as.numeric(quantile(ds$b_id))))   +  geom_smooth(method=lm) # lm or loess

######################
###### multiplot
######################

corr.dr.plot <- plot_grid(a,b,c,d,  labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
ggsave("Figure 1.jpeg", corr.dr.plot, dpi = 500, height = 190, width = 190, limitsize = F, units = "mm")
ggsave("Figure 1.tiff", corr.dr.plot, dpi = 500, height = 190, width = 190, limitsize = F, units = "mm")

###########################################################################################
# Figure 2
###########################################################################################

library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
library(stringr)
library(factoextra)
library(FactoMineR)
library(rafalib)
library(RSEIS)
library(ggsci)
setwd("D:/...")

######################
###### raw data
######################

dsr <-as.data.frame(fread(input = "xcms after IPO.csv", header=T)) # first column with all metadata
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]
# metadata generating
rname <- rownames(dsr) 
rname <- str_remove(rname, ".CDF")
all_id <- sapply(1:length(rname), function(y) unlist(str_split(rname[y], " "))) 
ro_id <- as.numeric(unlist(lapply(all_id, function(y) unlist(y[1])))) 
b_id <- unlist(lapply(all_id, function(y) unlist(y[3]))) 
s_id <- unlist(lapply(all_id, function(y) unlist(y[2]))) 
p1_id <- unlist(lapply(all_id, function(y) unlist(y[4]))) 
p2_id <- unlist(lapply(all_id, function(y) unlist(y[5]))) 
p_all_id <- sapply(1:length(p1_id), function(y) ifelse(is.na(p2_id[y]), p1_id[y] ,paste(p1_id[y], p2_id[y], sep = " "))) 
un_raw_l <- unique(gsub("[[:digit:]]", "", p1_id)) 
true_l <- c("QC", "TG", "CG") 
rbind(un_raw_l,true_l) 
raw_l <- gsub("[[:digit:]]", "", p1_id) 
n_gr_t <- as.character(match(raw_l, un_raw_l)) 
for (i in 1:length(unique(n_gr_t))) {
  n_gr_t <- str_replace(n_gr_t, unique(n_gr_t)[i], true_l[i]) } 
s_id_p <- sapply(1:length(s_id), function(y) unlist(str_split(s_id[y], "_"))) 
s_id_p2 <- unlist(lapply(s_id_p, function(y) unlist(y[1]))) 
all_meta_data<- as.data.frame(cbind(n_gr_t, ro_id, b_id, s_id, p_all_id, raw_l, s_id_p2)) 
# join features and metadata and order
ds <- data.frame(cbind(all_meta_data, dsr))
ds$ro_id <- as.numeric(ds$ro_id)
ds <- ds[order(ds$ro_id, decreasing = F),] 
ds$b_id <- as.factor(stringr::str_remove(ds$b_id, "b"))
levels(ds$b_id) <- as.factor(1:max(as.numeric(ds$b_id)))

# dataset
base1 <- ds
mtrx1 <- ds[,-c(1:7)] 
grp1 <- base1[,3] # 3 for batch; 1 for label

palette_pca <- "category20" # color: "lancet" or "category20" or pal_gsea("default", n = length(unique(grp1)), alpha = 0.6, reverse = T)(length(unique(grp1)))
# palette_pca <- JGRAY(length(unique(grp1))) # grey

pca.ds1 <- PCA(mtrx1, scale.unit = T, graph = F)
a <- fviz_pca_ind(pca.ds1,
                    title = "",
                    geom.ind = "point", # show points only 
                    col.ind = grp1, # color by groups
                    palette = palette_pca, # color "jco" gray JGRAY(length(unique(grp1)))
                    addEllipses = F, # Concentration ellipses
                    legend.title = "Groups") +scale_shape_manual(values=rep(0:length(unique(grp1))))

# dataset
base1 <- ds
mtrx1 <- ds[,-c(1:7)] 
grp1 <- base1[,1] # 3 for batch; 1 for label

palette_pca <- "lancet" # color: "lancet" or "category20" or pal_gsea("default", n = length(unique(grp1)), alpha = 0.6, reverse = T)(length(unique(grp1)))
# palette_pca <- JGRAY(length(unique(grp1))) # grey

pca.ds1 <- PCA(mtrx1, scale.unit = T, graph = F)
b <- fviz_pca_ind(pca.ds1,
                  title = "",
                  geom.ind = "point", # show points only 
                  col.ind = grp1, # color by groups
                  palette = palette_pca, # color "jco" gray JGRAY(length(unique(grp1)))
                  addEllipses = F, # Concentration ellipses
                  legend.title = "Groups") 

######################
###### QC-XGB data
######################

dsr <-as.data.frame(fread(input = "xcms after IPO QC-XGB.csv", header=T)) # first column with all metadata
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]
# metadata generating
rname <- rownames(dsr) 
rname <- str_remove(rname, ".CDF")
all_id <- sapply(1:length(rname), function(y) unlist(str_split(rname[y], " "))) 
ro_id <- as.numeric(unlist(lapply(all_id, function(y) unlist(y[1])))) 
b_id <- unlist(lapply(all_id, function(y) unlist(y[3]))) 
s_id <- unlist(lapply(all_id, function(y) unlist(y[2]))) 
p1_id <- unlist(lapply(all_id, function(y) unlist(y[4]))) 
p2_id <- unlist(lapply(all_id, function(y) unlist(y[5]))) 
p_all_id <- sapply(1:length(p1_id), function(y) ifelse(is.na(p2_id[y]), p1_id[y] ,paste(p1_id[y], p2_id[y], sep = " "))) 
un_raw_l <- unique(gsub("[[:digit:]]", "", p1_id)) 
true_l <- c("QC", "TG", "CG") 
rbind(un_raw_l,true_l) 
raw_l <- gsub("[[:digit:]]", "", p1_id) 
n_gr_t <- as.character(match(raw_l, un_raw_l)) 
for (i in 1:length(unique(n_gr_t))) {
  n_gr_t <- str_replace(n_gr_t, unique(n_gr_t)[i], true_l[i]) } 
s_id_p <- sapply(1:length(s_id), function(y) unlist(str_split(s_id[y], "_"))) 
s_id_p2 <- unlist(lapply(s_id_p, function(y) unlist(y[1]))) 
all_meta_data<- as.data.frame(cbind(n_gr_t, ro_id, b_id, s_id, p_all_id, raw_l, s_id_p2)) 
# join features and metadata and order
ds <- data.frame(cbind(all_meta_data, dsr))
ds$ro_id <- as.numeric(ds$ro_id)
ds <- ds[order(ds$ro_id, decreasing = F),] 
ds$b_id <- as.factor(stringr::str_remove(ds$b_id, "b"))
levels(ds$b_id) <- as.factor(1:max(as.numeric(ds$b_id)))

# dataset
base1 <- ds
mtrx1 <- ds[,-c(1:7)] 
grp1 <- base1[,3] # 3 for batch; 1 for label

palette_pca <- "category20" # color: "lancet" or "category20" or pal_gsea("default", n = length(unique(grp1)), alpha = 0.6, reverse = T)(length(unique(grp1)))
# palette_pca <- JGRAY(length(unique(grp1))) # grey

pca.ds1 <- PCA(mtrx1, scale.unit = T, graph = F)
c <- fviz_pca_ind(pca.ds1,
                  title = "",
                  geom.ind = "point", # show points only 
                  col.ind = grp1, # color by groups
                  palette = palette_pca, # color "jco" gray JGRAY(length(unique(grp1)))
                  addEllipses = F, # Concentration ellipses
                  legend.title = "Groups") +scale_shape_manual(values=rep(0:length(unique(grp1))))

# dataset
base1 <- ds
mtrx1 <- ds[,-c(1:7)] 
grp1 <- base1[,1] # 3 for batch; 1 for label

palette_pca <- "lancet" # color: "lancet" or "category20" or pal_gsea("default", n = length(unique(grp1)), alpha = 0.6, reverse = T)(length(unique(grp1)))
# palette_pca <- JGRAY(length(unique(grp1))) # grey

pca.ds1 <- PCA(mtrx1, scale.unit = T, graph = F)
d <- fviz_pca_ind(pca.ds1,
                  title = "",
                  geom.ind = "point", # show points only 
                  col.ind = grp1, # color by groups
                  palette = palette_pca, # color "jco" gray JGRAY(length(unique(grp1)))
                  addEllipses = F, # Concentration ellipses
                  legend.title = "Groups") 

######################
###### multiplot
######################

pca.qc.plot <- plot_grid(a,b,c,d,  labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
ggsave("Figure 2.jpeg", pca.qc.plot, dpi = 500, height = 190, width = 190, limitsize = F, units = "mm")
ggsave("Figure 2.tiff", pca.qc.plot, dpi = 500, height = 190, width = 190, limitsize = F, units = "mm")

###########################################################################################
# Figure 3
###########################################################################################

library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
library(stringr)
library(factoextra)
library(FactoMineR)
library(rafalib)
library(RSEIS)
library(ggsci)
setwd("D:/...")

######################
###### raw data
######################

ds <-as.data.frame(fread(input = "xcms after IPO MVI filter repeats annot.csv", header=T)) # first column with all metadata
ds <- ds[-c(1:8),]
rownames(ds) <- ds[,5]
ds <- ds[,-c(1,3:5)]
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# dataset
base1 <- ds
mtrx1 <- ds[,-1] 
grp1 <- base1[,1] 

palette_pca <- "lancet" # color: "lancet" or "category20" or pal_gsea("default", n = length(unique(grp1)), alpha = 0.6, reverse = T)(length(unique(grp1)))
# palette_pca <- JGRAY(length(unique(grp1))) # grey

pca.ds1 <- PCA(mtrx1, scale.unit = T, graph = F)
a <- fviz_pca_ind(pca.ds1,
                  title = "",
                  geom.ind = "point", # show points only 
                  col.ind = grp1, # color by groups
                  palette = palette_pca, # color "jco" gray JGRAY(length(unique(grp1)))
                  addEllipses = F, # Concentration ellipses
                  legend.title = "", pointsize = 0.5) + theme(legend.position="top")

######################
###### QC-XGB
######################

ds <-as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats annot+filtr.csv", header=T)) # first column with all metadata
ds <- ds[-c(1:8),]
rownames(ds) <- ds[,5]
ds <- ds[,-c(1,3:5)]
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# dataset
base1 <- ds
mtrx1 <- ds[,-1] 
grp1 <- base1[,1] 

palette_pca <- "lancet" # color: "lancet" or "category20" or pal_gsea("default", n = length(unique(grp1)), alpha = 0.6, reverse = T)(length(unique(grp1)))
# palette_pca <- JGRAY(length(unique(grp1))) # grey

pca.ds1 <- PCA(mtrx1, scale.unit = T, graph = F)
b <- fviz_pca_ind(pca.ds1,
                  title = "",
                  geom.ind = "point", # show points only 
                  col.ind = grp1, # color by groups
                  palette = palette_pca, # color "jco" gray JGRAY(length(unique(grp1)))
                  addEllipses = F, # Concentration ellipses
                  legend.title = "", pointsize = 0.5)+ theme(legend.position="top")

######################
###### QC-XGB+LMM
######################

ds <-as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats annot+filtr LMM adj KEGG.csv", header=T)) # first column with all metadata
ds <- ds[-c(1:8),]
rownames(ds) <- ds[,5]
ds <- ds[,-c(1,3:5)]
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# dataset
base1 <- ds
mtrx1 <- ds[,-1] 
grp1 <- base1[,1] 

palette_pca <- "lancet" # color: "lancet" or "category20" or pal_gsea("default", n = length(unique(grp1)), alpha = 0.6, reverse = T)(length(unique(grp1)))
# palette_pca <- JGRAY(length(unique(grp1))) # grey

pca.ds1 <- PCA(mtrx1, scale.unit = T, graph = F)
c <- fviz_pca_ind(pca.ds1,
                  title = "",
                  geom.ind = "point", # show points only 
                  col.ind = grp1, # color by groups
                  palette = palette_pca, # color "jco" gray JGRAY(length(unique(grp1)))
                  addEllipses = F, # Concentration ellipses
                  legend.title = "", pointsize = 0.5)+ theme(legend.position="top")

######################
###### multiplot
######################

pca.rep.plot <- plot_grid(a,b,c,  labels = c("A", "B", "C"), ncol = 3, nrow = 1)
ggsave("Figure 3.jpeg", pca.rep.plot, dpi = 500, height = 90, width = 140, limitsize = F, units = "mm")
ggsave("Figure 3.tiff", pca.rep.plot, dpi = 500, height = 90, width = 140, limitsize = F, units = "mm")

###########################################################################################
# Figure S-1
###########################################################################################

library(ggplot2)
library(ggpubr)
library(cowplot)
library(data.table)
library(dplyr)
library(stringr)
library(NormalizeMets)
library(htmlwidgets)
library(rsvg)
setwd("D:/...")

######################
###### raw data
######################

ds <-as.data.frame(fread(input = "xcms after IPO MVI filter repeats annot.csv", header=T)) # first column with all metadata
ds <- ds[-c(1:8),]
rownames(ds) <- ds[,5]
ds <- ds[,-c(1,3:5)]
ds[,-1] <- sapply(ds[,-1], as.numeric)
log_data1 <- LogTransform(featuredata = ds[,-1])$featuredata

######################
###### QC-XGB
######################

ds <-as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats annot+filtr.csv", header=T)) # first column with all metadata
ds <- ds[-c(1:8),]
rownames(ds) <- ds[,5]
ds <- ds[,-c(1,3:5)]
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)
log_data2 <- LogTransform(featuredata = ds[,-1])$featuredata

######################
###### QC-XGB+LMM
######################

ds <-as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats annot+filtr LMM adj KEGG.csv", header=T)) # first column with all metadata
ds <- ds[-c(1:12),]
rownames(ds) <- ds[,5]
ds <- ds[,-c(1:5)]
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)
log_data3 <- LogTransform(featuredata = ds[,-1])$featuredata

######################
###### multiplot
######################

lfeaturedata<-list(d1=log_data1,d2=log_data2,d3=log_data3)
Label <- as.factor(ds$Label)
n <- CompareRlaPlots(lfeaturedata,
                     groupdata=Label,
                     normmeth=c("raw:", "QC-XGB:", 
                                "LMM-adj:"), yrange=c(-4,4), plottitle = "")

n %>% htmlwidgets::onRender(
  "function(el, x) {
  var gd = document.getElementById(el.id); 
  Plotly.downloadImage(gd, {format: 'jpeg', width: 650, height: 450, scale: 5, filename: 'RLA plots'});
  }"
)

n %>% htmlwidgets::onRender(
  "function(el, x) {
  var gd = document.getElementById(el.id); 
  Plotly.downloadImage(gd, {format: 'svg', width: 650, height: 450, scale: 5, filename: 'RLA plots'});
  }"
)

f <- "D:/.../Figure S-1.svg"
rsvg_pdf(f, "Figure S-1.pdf")

###########################################################################################
# Figure 4
###########################################################################################

# packages
library(factoextra)
library(FactoMineR)
library(dendextend)
library(rafalib)
library(RSEIS)
library(ggsci)
library(ggplotify)
library(ggpubr)
library(Rtsne)
library(pROC)
library(parallel)
library(doParallel)
library(grid)
library(caret)
library(dplyr)
library(MKinfer)
library(limma)
library(ggplot2)
library(cowplot)

######################
###### after filtering
######################

# load data
library(data.table)
setwd("D:/...")

ds <- as.data.frame(fread(input = "8 peaks.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <-"Label"
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

# PCA
base1 <- ds # dataset
mtrx1 <- ds[,-1] # numeric data
grp1 <- as.character(base1[,1]) # label of dataset
palette_pca <- "lancet"
pca.ds1 <- PCA(mtrx1, scale.unit = T, graph = F)
a <- fviz_pca_ind(pca.ds1,
                    title = "",
                    geom.ind = "point", # show points only 
                    col.ind = grp1, # color by groups
                    palette = palette_pca, # color "jco" gray JGRAY(length(unique(grp1)))
                    addEllipses = T, # Concentration ellipses
                    legend.title = "Groups")

# HCA
base1 <- ds # dataset
mtrx1 <- ds[,-1] # numeric data
grp1 <- as.character(base1[,1]) # label of dataset
k <- length(unique(grp1)) # groups in HC
Cols = function(vec, ord){
  cols = pal_jco(palette = c("default"), alpha = 1)(length(unique(vec)))
  return(cols[as.fumeric(vec)[ord]])}
mtrx1_1 <- mtrx1
#mtrx1_1 <- data.frame(scale(mtrx1, center = T, scale = T))
rownames(mtrx1_1) = make.names(grp1, unique=TRUE)
res.dist1 <- dist(mtrx1_1, method = "manhattan") #{euclidean}, {maximum}, {manhattan}, {canberra}, {binary}, {minkowski}
res.hc1 <- hclust(d = res.dist1, method = "ward.D2") #{ward (ward.D), (ward.D2)}, {single}, {complete}, {average}, {mcquitty},{median}, {centroid}
b <- fviz_dend(res.hc1, k = k, # Cut in k groups
                 cex = 0.55, # label size 0.3/0.7
                 k_colors = unique(Cols(grp1,res.hc1$order)), # "lancet" color "jco" gray JGRAY(k_hc)
                 color_labels_by_k = F, # color labels by groups
                 label_cols = Cols(grp1,res.hc1$order),#Cols(ds[,1])[res.hc1$order], #as.fumeric(ds[,1])[res.hc1$order]
                 rect = T, # Add rectangle around groups
                 rect_fill = T,
                 # rect_border = unique(Cols(grp1,res.hc1$order)), #"lancet"# color "jco" gray JGRAY(k_hc)
                 horiz = F,
                 lwd = 0.7, # lines size 0.3/0.7
                 show_labels = T,
                 main = "",
                 ylab = "")

# ROC curve based on SVM predictions

# stop parallel
stopCluster(cl)
stopImplicitCluster()

# start parallel processing
fc <- as.numeric(detectCores(logical = T))
cl <- makePSOCKcluster(fc-1)
registerDoParallel(cl)

set.seed(1234) 
trainIndex <- createDataPartition(ds$Label, p = 0.8, list = F, times = 1) # or simple unbalanced splitting: sample(2, length(ds$Label), replace = T, prob=c(0.8, 0.2))
dsTrain <- ds[ trainIndex,]
dsValid <- ds[-trainIndex,]
trainControl <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs = T) # or bootstrap: trainControl(method="boot", number=100);
metric <- "Accuracy" 
set.seed(1234)
fit.cl <- train(Label~., data=dsTrain, method="svmRadial", metric=metric, trControl=trainControl, tuneLength = 10)
predicted.classes <- predict(fit.cl, newdata=dsValid)
probabilities <- predict(fit.cl, newdata=dsValid, type = "prob")[,1]                         
res.roc <- roc(dsValid$Label, probabilities, levels = levels(dsValid$Label))
auroc <- round(as.numeric(auc(res.roc)),2)
grob <- grobTree(textGrob(paste0("AUC = ", auroc), x=0.25,  y=0.75, hjust=0,
                          gp=gpar(col="firebrick2", fontsize=15, fontface=11)))
c <- ggroc(res.roc, alpha = 0.5, colour = "blue1", linetype = 1, size = 1.5) +theme_bw() + ggtitle("") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed")  + annotation_custom(grob)

# bootstrap histogram

# stop parallel
stopCluster(cl)
stopImplicitCluster()

# start parallel processing
fc <- as.numeric(detectCores(logical = T))
cl <- makePSOCKcluster(fc-1)
registerDoParallel(cl)

set.seed(1234)
trainControl <-trainControl(method="boot", number=1000)
metric <- "Accuracy" 
set.seed(1234)
fit.cl <- train(Label~., data=ds, method="svmRadial", metric=metric, trControl=trainControl)
results <- resamples(list(svm=fit.cl, svm1=fit.cl), trControl = trainControl, metric=metric)
d <-ggplot(results$values, aes(x=results$values[,2])) + 
  geom_histogram(colour="blue", fill="white") + theme_bw() + xlab("Accuracy") + ylab("Frequency")

# Volcano plot
setwd("D:/...")
ds <- as.data.frame(fread(input = "xcms after IPO MVI QC-XGB filter repeats annot+filtr LMM adj KEGG.csv", header=T))
setwd("D:/...")
ds <- ds[-c(1:12),]
rownames(ds) <- ds[,5]
ds <- ds[,-c(1:5)]
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

ds_log <- as.data.frame(log2(ds[,-1]))
ds_log <- cbind(Label = ds[,1], ds_log)
FOLD.CHANGE <- function(data) {
  ds_log_subsets <- lapply(1:length(unique(data[,1])), function(y) filter(data[,-1], data$Label == unique(data[,1])[y])) # list of subsets by label
  mean_r_l <- lapply(1:length(ds_log_subsets), function(y) apply(ds_log_subsets[[y]], 2, mean, na.rm = T)) # calculate mean for feature
  foldchange <- mean_r_l[[1]] - mean_r_l[[2]]
  fc_res <- as.data.frame(foldchange)
  return(fc_res)
}
fc_res <- FOLD.CHANGE(ds_log)
foldchange <- as.numeric(fc_res$foldchange)

mdl_mtrx <- model.matrix(~Label, ds)
lmf <- lmFit(t(ds[,-1]), method = "robust", design = mdl_mtrx, maxit = 1000) # "robust" or "ls"
efit <- eBayes(lmf)
tableTop <- topTable(efit, coef = 2, adjust = "BH",  number = ncol(ds), sort.by = "none")
pval <- as.numeric(tableTop$adj.P.Val)

f <- volcano(foldchange, pval, effect.low = -1.0, effect.high = 1.0, sig.level = 0.00001,
             xlab = "log2(fold change)", ylab = "-log10(adj.p value)",  title = "") + theme_classic() +  theme(legend.position="none") 

######################
###### multiplot
######################
# order: f,a,b,c,d

row2 <- plot_grid(a, b, labels = c('B', 'C'), label_size = 12)
row3 <- plot_grid(c, d, labels = c('D', 'E'), label_size = 12)
stats.plot <- plot_grid(f, row2, row3, labels = c('A', ''), label_size = 12, ncol = 1, nrow =3)

ggsave("Figure 4.jpeg", stats.plot, dpi = 500, height = 240, width = 190, limitsize = F, units = "mm")
ggsave("Figure 4.tiff", stats.plot, dpi = 500, height = 240, width = 190, limitsize = F, units = "mm")

###########################################################################################
# Figure S-2
###########################################################################################

# packages
library(ggplot2)
library(reshape2)
library(cowplot)
library(scipub)

# load data
library(data.table)
setwd("D:/...")

ds <- as.data.frame(fread(input = "annotated 8 peaks.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <-"Label"
ds[,-1] <- sapply(ds[,-1], as.numeric)
ds$Label <- as.factor(ds$Label)

#################################### box plots
df.m <- melt(ds, id.var = "Label") # reshape data frame

# p <- ggplot(data = df.m, aes(x=variable, y=value)) + xlab("") + ylab("") +
#  geom_boxplot(aes(fill=Label)) + theme(legend.position="bottom") + theme_classic() + scale_x_discrete(labels=c("")) 

# a <- p + facet_wrap( ~ variable, scales="free") + theme_classic() + theme(legend.position="bottom") 

a <- gg_groupplot(data=df.m, x=Label, y=value, meanline = F) + facet_wrap(~variable, scales="free") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ theme_classic()+theme(legend.position="none") 
  
# a

#################################### scatter plots

ds1 <- as.data.frame(cbind(Patient = rownames(ds), ds))
df.m <- melt(ds, id.var = "Label")
df.m <- cbind(1:nrow(df.m), df.m)
colnames(df.m)[1] <- "Patient"
p <- ggplot(data = df.m, aes(x=Patient, y=value)) + xlab("") + ylab("") +
  geom_point(aes(colour=Label)) + theme(legend.position="bottom") + theme_classic() + scale_x_discrete(labels=c("")) + geom_smooth(method = "lm")

b <- p + facet_wrap( ~ variable, scales="free") + theme_classic() + theme(legend.position="none") 

# b

######################
###### multiplot
######################

annot.plot <- plot_grid(a,b,  labels = c("A", "B"), ncol = 1, nrow = 2)
ggsave("Figure S-2.jpeg", annot.plot, dpi = 500, height = 190, width = 190, limitsize = F, units = "mm")
ggsave("Figure S-2.tiff", annot.plot, dpi = 500, height = 190, width = 190, limitsize = F, units = "mm")
