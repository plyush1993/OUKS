##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Testing Data                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(!"pacman" %in% rownames(installed.packages())){
  install.packages("pacman")
}
cat("Checking required packages (auto-installing if missing)\n")
pacman::p_load("data.table", "dplyr", "stringr", "ggplot2", "ggsci", "mgcv", "remotes")

if(!require(ggmagnify)){
  install_github("hughjonesd/ggmagnify")
  library(ggmagnify)}

dsr <-as.data.frame(fread(input = "https://github.com/plyush1993/OUKS/raw/refs/heads/main/Datasets%20(csv)/xcms%20after%20IPO%20MVI.csv", header=T)) # first column with all metadata
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]

rname <- rownames(dsr) # obtain all info from rownames
rname <- str_remove(rname, ".CDF") # remove some pattern from vendor-specific format
all_id <- lapply(1:length(rname), function(y) unlist(str_split(rname[y], " "))) # split info from rownames
ro_id <- as.numeric(unlist(lapply(all_id, function(y) unlist(y[1])))) # obtain run order ID (every [1] element)
b_id <- unlist(lapply(all_id, function(y) unlist(y[3]))) # obtain batch ID (every [3] element)
s_id <- unlist(lapply(all_id, function(y) unlist(y[2]))) # obtain sample ID (every [2] element)
p1_id <- unlist(lapply(all_id, function(y) unlist(y[4]))) # obtain patient ID (every [4] element)
p2_id <- unlist(lapply(all_id, function(y) unlist(y[5]))) # obtain patient ID (every [5] element)
p_all_id <- sapply(1:length(p1_id), function(y) ifelse(is.na(p2_id[y]), p1_id[y] ,paste(p1_id[y], p2_id[y], sep = " "))) # join patient ID (every [4-5] element)
un_raw_l <- unique(gsub("[[:digit:]]", "", p1_id)) # obtain unique raw label from p1_id (every [4] element)
true_l <- c("QC", "TG", "CG") # write desired label in order as in unique raw label (should be analogical)
rbind(un_raw_l,true_l) # visual check agreement between true and raw labels
raw_l <- gsub("[[:digit:]]", "", p1_id) # obtain raw label from p1_id (every [4] element)
n_gr_t <- as.character(match(raw_l, un_raw_l)) # obtain index for raw label by unique raw label
for (i in 1:length(unique(n_gr_t))) {
  n_gr_t <- str_replace(n_gr_t, unique(n_gr_t)[i], true_l[i]) } # exchange raw label by true label via index
s_id_p <- sapply(1:length(s_id), function(y) unlist(str_split(s_id[y], "_"))) # obtain sample ID (every [2] element) with splitting by "_"

s_id_p2 <- unlist(lapply(s_id_p, function(y) unlist(y[1]))) # obtain sample ID (every [2] element) with only one part of sample ID
all_meta_data<- as.data.frame(cbind(n_gr_t, ro_id, b_id, s_id, p_all_id, raw_l, s_id_p2)) # obtain all metadata for all repeats

ds <- data.frame(cbind(all_meta_data, dsr)) # add metadata columns

ds$ro_id <- as.numeric(ds$ro_id)
ds <- ds[order(ds$ro_id, decreasing = F),] # order by ro col

n <- 7 # set non numeric columns, adjust to your data
meta <- ds[,c(1:n)]
ds <- ds[,-c(1:n)]
order = meta$ro_id
batch = as.numeric(str_remove(meta$b_id, "b"))
class = meta$n_gr_t
qc_label = "QC"

qc_id <-  grep(qc_label, class)
batch_qc <- batch[qc_id]
ro_qc <- order[qc_id]
ro = order
tic_int_data_qc = as.numeric(rowSums(ds[qc_id,], na.rm = T))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              QC-GAM (MetCorR)                            ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#......................................................
# Parameters of "MetCorR" function:  

# "method" -> can be 1 (1 feature mode, only run order) or 2 (2 features mode, run order & batch index)
# "int_data" -> numeric feature table
# "order" -> numeric vector of the run order
# "class" -> sample group variable 
# "batch" -> numeric vector of the batch index
# "qc_label" -> label for QC samples in group 
#......................................................

# NOTE: place "MetCorR.R" (from Auxiliary files (RData) folder) to "workingdirectory"
source("https://github.com/plyush1993/OUKS/raw/refs/heads/main/Auxiliary%20files%20(RData)/MetCorR.R")

# 1 feature mode, only run order
ds_gam1 <- MetCorR(method = 1, int_data = ds, order = order, class = class, qc_label = qc_label)
fwrite(ds_gam1, "xcms after IPO MVI QC-GAM1.csv", row.names = T)

# 2 features mode, run order & batch
ds_gam2 <- MetCorR(method = 2, int_data = ds, order = order, class = class, batch = batch, qc_label = qc_label)
fwrite(ds_gam2, "xcms after IPO MVI QC-GAM2.csv", row.names = T)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                             Data Visualization                           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                 log2(TIC) as a model signal                 ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#........................models overview.........................
# with by s term
model1r <- gam(log2(tic_int_data_qc) ~ s(ro_qc, by = batch_qc))
vis.gam(model1r, view = c("ro_qc", "batch_qc"),
        theta = 50, n.grid = 50, lwd = 0.4)
predict_gam1r = predict(model1r, newdata=data.frame(ro_qc = ro, batch_qc = batch))

# with 2 s terms
model2r <- gam(log2(tic_int_data_qc) ~ s(ro_qc)+s(batch_qc))
predict_gam1 = predict(model2r, newdata=data.frame(ro_qc = ro, batch_qc = batch))
vis.gam(model2r, view = c("ro_qc", "batch_qc"),
        theta = 50, n.grid = 50, lwd = 0.4)
predict_gam2r = predict(model2r, newdata=data.frame(ro_qc = ro, batch_qc = batch))

# with simple s term
model1 <- gam(log2(tic_int_data_qc) ~ s(ro_qc))
predict_gam1 = predict(model1, newdata=data.frame(ro_qc = ro, batch_qc = batch))

# with interaction s term
model2 <- gam(log2(tic_int_data_qc) ~ s(ro_qc, batch_qc))
predict_gam2 = predict(model2, newdata=data.frame(ro_qc = ro, batch_qc = batch))
vis.gam(model2, view = c("ro_qc", "batch_qc"),
        theta = 50, n.grid = 50, lwd = 0.4)

#...........................basic plot...........................
plot(y=log2(tic_int_data_qc), x = ro_qc, xlab = "Order", ylab = "log2(TIC)", main='GAM Modeling')
lines(predict_gam1, col='red', lwd = 2)
lines(predict_gam2, col='blue', lwd = 2)
lines(predict_gam1r, col='darkgreen', lwd = 2)
lines(predict_gam2r, col='orange', lwd = 2)
legend('topright', legend=c('s(ro_qc)', 's(ro_qc, batch_qc)', 's(ro_qc, by = batch_qc)', 's(ro_qc)+s(batch_qc)'),
       col=c('red', 'blue', "darkgreen", "orange"), pch=19)

#.............................ggplot.............................
init_data <- as.data.frame(cbind(Abundance = log2(tic_int_data_qc), Order = ro_qc, Batch = batch_qc))
data <- as.data.frame(cbind(Order = ro, 
                            's(order)' = predict_gam1,
                            's(order, batch)' = predict_gam2,
                            's(order, by=batch)' = predict_gam1r,
                            's(order)+s(batch)' = predict_gam2r))
data <- melt(data, id = "Order")

ggplot(init_data, aes(y = Abundance, x = Order)) + 
  geom_rect(aes(xmin=95, xmax=117, ymin=29, ymax= 29.4), color = "black", fill="grey89") +
  geom_rect(aes(xmin=197, xmax=225, ymin=28.5, ymax= 29), color = "black", fill="grey89") +
  geom_point(size=3.5, color="black", fill= "white", shape=21, stroke=1) +
  theme_classic(base_size = 15) + theme(legend.position = c(0.8, 0.85)) + labs(fill = "Batch")+labs(colour = "Model:") +
  geom_line(data = data, aes(x = Order, y = value, colour = variable), alpha = 1, 
            linewidth = 1.2)+scale_color_d3() 

ggplot(init_data, aes(y = Abundance, x = Order)) + 
  geom_point(size=3.5, color="black", fill= "grey90", shape=21, stroke=1) + 
  theme_linedraw(base_size = 15) + theme(legend.position = "top", legend.title=element_blank(),
                                   legend.key.height  = unit(0.02, 'cm'), legend.text = element_text(size=15)) + 
  labs(fill = "Batch")+labs(colour = "Model:") +
  geom_line(data = data, aes(x = Order, y = value, color = variable), alpha = 1, linewidth = 1.5)+scale_color_jco() + 
    geom_magnify(from = c(xmin = 95, xmax = 117, ymin = 29, ymax = 29.4), 
               to = c(60, 135, 28.4, 28.9),
               proj.fill = alpha("yellow", 0.2)) +
  geom_magnify(from = c(xmin = 197, xmax = 225, ymin = 28.5, ymax = 29), 
               to = c(160, 250, 29.2, 29.7),
               proj.fill = alpha("yellow", 0.2))+
  guides(color = guide_legend(keyheight = unit(0.1, 'cm')))
  
ggplot(init_data, aes(y = Abundance, x = Order)) + 
  geom_rect(aes(xmin=95, xmax=117, ymin=29, ymax= 29.4), color = "black", fill="grey89") +
  geom_rect(aes(xmin=197, xmax=225, ymin=28.5, ymax= 29), color = "black", fill="grey89") +
  geom_point(size=3.5, color="black", aes(fill= as.factor(Batch)), alpha=0.5, shape=21, stroke=0.5) + scale_fill_viridis_d(guide = 'none')+
  theme_bw(base_size = 15) + theme(legend.position = "none") + labs(fill = "Batch")+labs(colour = "Model:") +
  geom_line(data = data, aes(x = Order, y = value, colour = variable), alpha = 1, linewidth = 1.5)+scale_color_aaas()+
  facet_grid(~variable)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##             particular feature as a model signal            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#................select feature by SD and mean int...............
sd_v <- apply(ds[qc_id,], 2, sd, na.rm = T)
mean_v <- apply(ds[qc_id,], 2, mean, na.rm = T)
hist(mean_v)
rsd_v <- sd_v/mean_v*100
hist(rsd_v)
#........................choose threshold........................
sel_features_ind <- which(rsd_v > 80 & mean_v > 10e5)
#.........................choose feature.........................
f_int <- ds[qc_id,which.max(rsd_v[sel_features_ind])]
plot(log2(f_int))

#........................models overview.........................
# with by s term
model1r <- gam(log2(f_int) ~ s(ro_qc, by = batch_qc))
vis.gam(model1r, view = c("ro_qc", "batch_qc"),
        theta = 150, n.grid = 50, lwd = 0.4)
predict_gam1r = predict(model1r, newdata=data.frame(ro_qc = ro, batch_qc = batch))

# with 2 s terms
model2r <- gam(log2(f_int) ~ s(ro_qc)+s(batch_qc))
predict_gam1 = predict(model2r, newdata=data.frame(ro_qc = ro, batch_qc = batch))
vis.gam(model2r, view = c("ro_qc", "batch_qc"),
        theta = 150, n.grid = 50, lwd = 0.4)
predict_gam2r = predict(model2r, newdata=data.frame(ro_qc = ro, batch_qc = batch))

# with simple s term
model1 <- gam(log2(f_int) ~ s(ro_qc))
predict_gam1 = predict(model1, newdata=data.frame(ro_qc = ro, batch_qc = batch))

# with interaction s term
model2 <- gam(log2(f_int) ~ s(ro_qc, batch_qc))
predict_gam2 = predict(model2, newdata=data.frame(ro_qc = ro, batch_qc = batch))
vis.gam(model2, view = c("ro_qc", "batch_qc"),
        theta = 150, n.grid = 50, lwd = 0.4)

#...........................basic plot...........................
plot(y=log2(f_int), x = ro_qc, xlab = "Order", ylab = "log2(TIC)", main='GAM Modeling')
lines(predict_gam1, col='red', lwd = 2)
lines(predict_gam2, col='blue', lwd = 2)
lines(predict_gam1r, col='darkgreen', lwd = 2)
lines(predict_gam2r, col='orange', lwd = 2)
legend('topright', legend=c('s(ro_qc)', 's(ro_qc, batch_qc)', 's(ro_qc, by = batch_qc)', 's(ro_qc)+s(batch_qc)'),
       col=c('red', 'blue', "darkgreen", "orange"), pch=19)

#.............................ggplot.............................
init_data <- as.data.frame(cbind(Abundance = log2(f_int), Order = ro_qc, Batch = batch_qc))
data <- as.data.frame(cbind(Order = ro, 
                            's(order)' = predict_gam1,
                            's(order, batch)' = predict_gam2,
                            's(order, by=batch)' = predict_gam1r,
                            's(order)+s(batch)' = predict_gam2r))
data <- melt(data, id = "Order")

ggplot(init_data, aes(y = Abundance, x = Order)) + 
  geom_point(size=3.5, color="black", fill= "white", shape=21, stroke=1) +
  theme_classic(base_size = 15) + theme(legend.position = c(0.2, 0.25)) + labs(fill = "Batch")+labs(colour = "Model:") +
  geom_line(data = data, aes(x = Order, y = value, colour = variable), alpha = 1, 
            linewidth = 1.2)+scale_color_d3() 

ggplot(init_data, aes(y = Abundance, x = Order)) + 
  geom_point(size=3.5, color="black", fill= "grey90", shape=21, stroke=1) + 
  theme_linedraw(base_size = 15) + theme(legend.position = "top", legend.title=element_blank(),
                                         legend.key.height  = unit(0.02, 'cm'), legend.text = element_text(size=15)) + 
  labs(fill = "Batch")+labs(colour = "Model:") +
  geom_line(data = data, aes(x = Order, y = value, color = variable), alpha = 1, linewidth = 1.5)+scale_color_jco() + 
  guides(color = guide_legend(keyheight = unit(0.1, 'cm')))

ggplot(init_data, aes(y = Abundance, x = Order)) + 
  geom_point(size=3.5, color="black", aes(fill= as.factor(Batch)), alpha=0.5, shape=21, stroke=0.5) + scale_fill_viridis_d(guide = 'none')+
  theme_bw(base_size = 15) + theme(legend.position = "none") + labs(fill = "Batch")+labs(colour = "Model:") +
  geom_line(data = data, aes(x = Order, y = value, colour = variable), alpha = 1, linewidth = 1.5)+scale_color_aaas()+
  facet_grid(~variable)

#................................................................