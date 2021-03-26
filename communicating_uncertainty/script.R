###communicating uncertainty in epidemic models

library(squire)
library(ggplot2)
library(dplyr)
library(matrixStats)
library(reshape2)
library(patchwork)
library(gridExtra)
library(pracma)


###fixed scenario
fixed <- readRDS("fixed_scenario.RDS")

# grab ICU data
imv_cols <- grep("^IMV", colnames(fixed$output))
patients_fixed <- do.call(
  rbind, 
  lapply(seq_len(dim(fixed$output)[3]), function(x) {
    data.frame(date = rownames(fixed$output[,,x]),
               patients = rowSums(fixed$output[,imv_cols,x]),
               replicate = x)
  }))

patients_fixed$date <- as.Date(patients_fixed$date)
###remove dummy replicate and trim dates
patients_fixed_trim <- patients_fixed %>% subset(date>="2020-09-01" & replicate!=0)

#numeric date sequence for new truncated date list
patients_fixed_trim$date_no <- rep(seq(from=0,by=1,length.out=nrow(patients_fixed_trim%>%subset(replicate==1))),100)

##gives peak and time for the max peak
patients_fixed_trim <- patients_fixed_trim %>% group_by(replicate) %>% mutate(peak = max(patients),timing = which.max(patients)-1)

###summary statistics 

first_peak_data_fixed <- data.frame("replicate"=0,"peak_first"=0,"timing_first"=0)

#extract first peak
for(i in 1:100){
  rep_subset <- patients_fixed_trim %>% subset(replicate==i)
  #smooth curve to be able to extract just the two 'main' peaks
  rep_subset_smooth <- predict(loess(rep_subset$patients~rep_subset$date_no))
  #extract fix peak
  first_peak <- findpeaks(rep_subset_smooth)[1,2]
  #search within +/- 10 days as smoothing can shift true date
  rep_subset_time <- rep_subset %>% subset(date_no >= first_peak-10&date_no<=first_peak+10)
  true_first_peak <- max(rep_subset_time$patients)
  true_timing <- min(rep_subset_time$date_no[which(rep_subset_time$patients==true_first_peak)])
  first_peak_data_fixed_it <- data.frame("replicate"=i,"peak_first"=true_first_peak,"timing_first"=true_timing)
  first_peak_data_fixed <- rbind(first_peak_data_fixed,first_peak_data_fixed_it)
}

#remove dummy first row
first_peak_data_fixed <- first_peak_data_fixed[-1,]

##merge max peak df and first peak df
patients_fixed_trim <- merge(patients_fixed_trim,first_peak_data_fixed,by="replicate",all.x=TRUE)

##assign variable that will be used to colour
patients_fixed_trim$rank <- match(rank(patients_fixed_trim$peak_first), sort(unique(rank(patients_fixed_trim$peak_first))))

##median for the line on the plot
patients_fixed_trim <- patients_fixed_trim %>% group_by(date_no) %>%  mutate(median = median(patients))

##remove the realisation that stops the axes from aligning
patients_fixed_trim <- patients_fixed_trim %>% subset(peak_first < 5000)

plot_a <- ggplot(patients_fixed_trim,aes(x=date_no,y=patients,group=replicate,col=rank))+
  geom_line(alpha=0.5)+
  geom_line(aes(x=date_no,y=median), inherit.aes = FALSE, lwd = 2) +
  theme_bw()+
  scale_color_viridis_c(end = 0.8) +
  labs(x="Time (Days)",y="ICU Demand",tag="a")+
  theme(legend.position = "NULL")+ylim(c(0,4700))

plot_c <- ggplot(patients_fixed_trim,aes(timing_first, peak_first, color = rank)) + 
  geom_point() + 
  labs(x="Day of First Peak",y="ICU Demand at First Peak",tag="c")+
  theme_bw() +
  scale_color_viridis_c(end = 0.8) +
  scale_size(range=c(0.5,3), breaks=quantile(patients_fixed_trim$ll, seq(0,1,0.2))) +
  theme(legend.position = "none")


###trigger scenario
trigger <- readRDS("trigger_scenario.RDS")

# grab ICU data
imv_cols <- grep("^IMV", colnames(trigger$output))
patients_trigger <- do.call(
  rbind, 
  lapply(seq_len(dim(trigger$output)[3]), function(x) {
    data.frame(date = rownames(trigger$output[,,x]),
               patients = rowSums(trigger$output[,imv_cols,x]),
               replicate = x)
  }))

patients_trigger$date <- as.Date(patients_trigger$date)
###remove dummy replicate and trim dates
patients_trigger_trim <- patients_trigger %>% subset(date>="2020-09-01" & replicate!=0)

#numeric date sequence for new truncated date list
patients_trigger_trim$date_no <- rep(seq(from=0,by=1,length.out=nrow(patients_trigger_trim%>%subset(replicate==1))),100)

##gives peak and time for the max peak
patients_trigger_trim <- patients_trigger_trim %>% group_by(replicate) %>% mutate(peak = max(patients),timing = which.max(patients)-1)

###summary statistics 

first_peak_data_trigger <- data.frame("replicate"=0,"peak_first"=0,"timing_first"=0)

#extract first peak
for(i in 1:100){
  rep_subset <- patients_trigger_trim %>% subset(replicate==i)
  #smooth curve to be able to extract just the two 'main' peaks
  rep_subset_smooth <- predict(loess(rep_subset$patients~rep_subset$date_no))
  #extract fix peak
  first_peak <- findpeaks(rep_subset_smooth)[1,2]
  #search within +/- 10 days as smoothing can shift true date
  rep_subset_time <- rep_subset %>% subset(date_no >= first_peak-10&date_no<=first_peak+10)
  true_first_peak <- max(rep_subset_time$patients)
  true_timing <- min(rep_subset_time$date_no[which(rep_subset_time$patients==true_first_peak)])
  first_peak_data_trigger_it <- data.frame("replicate"=i,"peak_first"=true_first_peak,"timing_first"=true_timing)
  first_peak_data_trigger <- rbind(first_peak_data_trigger,first_peak_data_trigger_it)
}

#remove dummy first row
first_peak_data_trigger <- first_peak_data_trigger[-1,]

##merge max peak df and first peak df
patients_trigger_trim <- merge(patients_trigger_trim,first_peak_data_trigger,by="replicate",all.x=TRUE)

##assign variable that will be used to colour
patients_trigger_trim$rank <- match(rank(patients_trigger_trim$peak_first), sort(unique(rank(patients_trigger_trim$peak_first))))

##median for the line on the plot
patients_trigger_trim <- patients_trigger_trim %>% group_by(date_no) %>%  mutate(median = median(patients))

##remove the realisation that stops the axes from aligning
patients_trigger_trim <- patients_trigger_trim %>% subset(peak_first < 5000)

plot_b <- ggplot(patients_trigger_trim,aes(x=date_no,y=patients,group=replicate,col=rank))+
  geom_line(alpha=0.5)+
  geom_line(aes(x=date_no,y=median), inherit.aes = FALSE, lwd = 2) +
  theme_bw()+
  scale_color_viridis_c(end = 0.8) +
  labs(x="Time (Days)",y="ICU Demand",tag="b")+
  theme(legend.position = "NULL")+ylim(c(0,4700))

plot_d <- ggplot(patients_trigger_trim,aes(timing_first, peak_first, color = rank)) + 
  geom_point() + 
  labs(x="Day of First Peak",y="ICU Demand at First Peak",tag="d")+
  theme_bw() +
  scale_color_viridis_c(end = 0.8) +
  theme(legend.position = "none")



###y axes of a and b the same
combine <- cowplot::plot_grid(plot_a, plot_b, plot_c, plot_d, ncol = 2)
ggsave(plot=combine,filename="Figure1_firstpeak_ab_ysame.pdf", width = 6, height = 5)
ggsave(plot=combine,filename="Figure1_firstpeak_ab_ysame.png", width = 6, height = 5)

###y axes of a, b, c, d the same
combine_ysame <- cowplot::plot_grid(plot_a, plot_b, plot_c+ylim(c(0,4700)), plot_d+ylim(c(0,4700)), ncol = 2)
ggsave(plot=combine_ysame,filename="Figure1_firstpeak_abcd_ysame.pdf", width = 6, height = 5)
ggsave(plot=combine_ysame,filename="Figure1_firstpeak_abcd_ysame.png", width = 6, height = 5)

##as above but c and d more similar x axis
combine_xysame <- cowplot::plot_grid(plot_a, plot_b, plot_c+ylim(c(0,4700))+xlim(c(60,136)), plot_d+ylim(c(0,4700)), ncol = 2)
ggsave(plot=combine_xysame,filename="Figure1_firstpeak_abcd_xysame.pdf", width = 6, height = 5)
ggsave(plot=combine_xysame,filename="Figure1_firstpeak_abcd_xysame.png", width = 6, height = 5)

##all axes the same
combine_xysame_all <- cowplot::plot_grid(plot_a, plot_b, plot_c+ylim(c(0,4700))+xlim(c(0,181)), plot_d+ylim(c(0,4700))+xlim(c(0,181)), ncol = 2)
ggsave(plot=combine_xysame,filename="Figure1_firstpeak_abcd_xysame_all.pdf", width = 6, height = 5)
ggsave(plot=combine_xysame,filename="Figure1_firstpeak_abcd_xysame_all.png", width = 6, height = 5)


