## communicating uncertainty in epidemic models

library(squire)
library(tidyverse)
library(matrixStats)
library(reshape2)
library(patchwork)
library(gridExtra)
library(pracma)
library(cowplot)


## figure 1 

cR0 <- c(3, 1.4, 2, 1.8, 1.1, 1.3, 2)
ctt <- c(0, 30, 40, 50, 60, 90, 110)

# deterministic model with different parameters
det_list <- lapply(1:10,
                   function(x) {
                     squire::run_deterministic_SEIR_model(
                       "Iran",
                       R0 = cR0 * runif(length(cR0), 0.95, 1.05),
                       tt_R0 = ctt * runif(length(ctt), 0.95, 1.05),
                       day_return = TRUE,
                       population = squire::get_population("Iran")$n/20,
                       contact_matrix_set = squire::get_mixing_matrix("Iran")
                     )
                     
                   })
# stochastic model with same parameters 
stoch <- squire::run_explicit_SEEIR_model("Iran",R0 = cR0,tt_R0 = ctt,day_return = TRUE,
                                          population = squire::get_population("Iran")$n/20,
                                          replicates = 10,
                                          contact_matrix_set = squire::get_mixing_matrix("Iran"))
# stochastic model with different parameters
stoch_list <- lapply(1:10,
                     function(x) {
                       squire::run_explicit_SEEIR_model(
                         "Iran",
                         R0 = cR0 * runif(length(cR0), 0.95, 1.05),
                         tt_R0 = ctt * runif(length(ctt), 0.95, 1.05),
                         day_return = TRUE,
                         replicates = 1,
                         population = squire::get_population("Iran")$n/20,
                         contact_matrix_set = squire::get_mixing_matrix("Iran")
                       )
                       
                     })

# plot stochastic model
stp <- plot(stoch, "deaths", replicates = TRUE, ci = FALSE,summarise = FALSE)
stp$layers[[1]]$aes_params$alpha <- 1
stp <- stp + scale_color_manual(values = "black")

# plot panels a, b, c
fig_1_a <- squire::projection_plotting(det_list, scenarios = letters[seq_along(det_list)], var_select = "deaths", add_parms_to_scenarios = FALSE) +
  theme(legend.position = "none") + ylab("Daily Deaths") + ylim(c(0,350)) + xlab("Time (Days)")
fig_1_b <- stp +
  theme(legend.position = "none", axis.title.y = element_blank()) + ylab("Daily Deaths") + ylim(c(0,350)) + xlab("Time (Days)")
fig_1_c <- squire::projection_plotting(stoch_list, scenarios = letters[seq_along(stoch_list)], var_select = "deaths", add_parms_to_scenarios = FALSE) +
  theme(legend.position = "none", axis.title.y = element_blank()) + ylab("Daily Deaths") + ylim(c(0,350)) + xlab("Time (Days)")

# panel a
fig_1_a$layers[[1]]$aes_params$alpha <- 0.2
fig_1_a$layers[[2]] <- NULL
fig_1_a <- fig_1_a + scale_color_manual(values = rep("black",10))

# panel b
fig_1_b$layers[[1]]$aes_params$alpha <- 0.2

# panel c
fig_1_c$layers[[1]]$aes_params$alpha <- 0.2
fig_1_c$layers[[2]] <- NULL
fig_1_c <- fig_1_c + scale_color_manual(values = rep("black",10))


### Figure 1: Uncertainty in epidemic models arising from different sources
cowplot::plot_grid(
  fig_1_a,
  fig_1_b,
  fig_1_c,
  ncol = 3, labels = "auto"
)


## For the rest of the figures, have to prepare data for the fixed and triggered scenarios and then plot at the end.

### fixed scenario

fixed <- readRDS("Data/fixed_scenario.RDS")

# grab ICU data
hosp_cols <- grep("^IOx", colnames(fixed$output))
hosp_fixed <- do.call(
  rbind, 
  lapply(seq_len(dim(fixed$output)[3]), function(x) {
    data.frame(date = as.Date(rownames(fixed$output[,,x])),
               patients = rowSums(fixed$output[,hosp_cols,x]),
               replicate = x,
               metric="Hospitalisations") %>% filter(date>="2020-09-01")
  }))

# grab ICU data
imv_cols <- grep("^IMV", colnames(fixed$output))
icu_fixed <- do.call(
  rbind, 
  lapply(seq_len(dim(fixed$output)[3]), function(x) {
    data.frame(date = as.Date(rownames(fixed$output[,,x])),
               patients = rowSums(fixed$output[,imv_cols,x]),
               replicate = x,
               metric="ICU") %>% filter(date>="2020-09-01")
  }))

# grab deaths data
death_cols <- 461:477
deaths_fixed <- data.frame(do.call(
  rbind,
  lapply(seq_len(dim(fixed$output)[3]), function(x) {
    data.frame(date = as.Date(rownames(fixed$output[,,x])),
               cumulative_deaths = rowSums(fixed$output[,death_cols,x]),
               replicate = x) %>% group_by(replicate) %>%
      mutate("patients"=c(NA,diff(cumulative_deaths)),
             "metric"="Death")
  }))) %>% select(date,patients,replicate,metric) %>% filter(date>="2020-09-01")


# bind three metrics together 
data_fixed <- rbind(hosp_fixed,icu_fixed,deaths_fixed)

# numeric date sequence for new truncated date list
data_fixed$date_no <- rep(rep(seq(from=0,by=1,length.out=nrow(deaths_fixed%>%filter(replicate==1))),100),3)

# gives peak and time for the max peak
data_fixed <- data_fixed %>% group_by(metric,replicate) %>% 
  mutate(peak = max(patients),timing = which.max(patients)-1)

# first peak of each trajectory for fixed scenarios

first_peak_data_fixed <- data.frame()

# extract first peak
for(i in 1:100){
  
  rep_hosp <- data_fixed %>% filter(replicate==i & metric=="Hospitalisations")
  rep_hosp_smooth <- predict(loess(rep_hosp$patients~rep_hosp$date_no))
  first_peak_hosp <- findpeaks(rep_hosp_smooth)[1,2]
  rep_time_hosp <- rep_hosp %>% subset(date_no >= first_peak_hosp-10&date_no<=first_peak_hosp+10)
  true_first_peak_hosp <- max(rep_time_hosp$patients)
  true_timing_hosp <- min(rep_time_hosp$date_no[which(rep_time_hosp$patients==true_first_peak_hosp)])
  
  rep_hosp_data <- data.frame("metric"="Hospitalisations","replicate"=i,
                              "peak_first"=true_first_peak_hosp,"timing_first"=true_timing_hosp)
  
  rep_icu <- data_fixed %>% filter(replicate==i & metric=="ICU")
  rep_icu_smooth <- predict(loess(rep_icu$patients~rep_icu$date_no))
  first_peak_icu <- findpeaks(rep_icu_smooth)[1,2]
  rep_time_icu <- rep_icu %>% subset(date_no >= first_peak_icu-10&date_no<=first_peak_icu+10)
  true_first_peak_icu <- max(rep_time_icu$patients)
  true_timing_icu <- min(rep_time_icu$date_no[which(rep_time_icu$patients==true_first_peak_icu)])
  
  rep_icu_data <- data.frame("metric"="ICU","replicate"=i,
                             "peak_first"=true_first_peak_icu,"timing_first"=true_timing_icu)
  
  
  rep_death <- data_fixed %>% filter(replicate==i & metric=="Death")
  rep_death_smooth <- predict(loess(rep_death$patients~rep_death$date_no))
  first_peak_death <- findpeaks(rep_death_smooth)[1,2]
  rep_time_death <- rep_death %>% subset(date_no >= first_peak_death-10&date_no<=first_peak_death+10)
  true_first_peak_death <- max(rep_time_death$patients)
  true_timing_death <- min(rep_time_death$date_no[which(rep_time_death$patients==true_first_peak_death)])
  
  rep_death_data <- data.frame("metric"="Death","replicate"=i,
                               "peak_first"=true_first_peak_death,"timing_first"=true_timing_death)
  
  ##want to rbind the individual reps then bind that with the bigger bind 
  first_peak_data_fixed_it <- rbind(rep_hosp_data,
                                    rep_icu_data,
                                    rep_death_data)
  
  first_peak_data_fixed <- rbind(first_peak_data_fixed,first_peak_data_fixed_it)
  
}

fixed_merge_check <- merge(data_fixed,first_peak_data_fixed,by=c("metric","replicate"),all.x=TRUE)

# manual fix for some replicates
icu_issues <- fixed_merge_check %>% filter(metric=="ICU") %>% select(replicate,peak,timing,peak_first,timing_first) %>% unique() %>% mutate("same"=peak-peak_first) %>% filter(same!=0) %>% select(replicate)

for(rep in icu_issues$replicate){
  first_peak_data_fixed$peak_first[which(first_peak_data_fixed$replicate==rep&
                                           first_peak_data_fixed$metric=="ICU")] <- fixed_merge_check %>% filter(replicate==rep&metric=="ICU"&date_no<100) %>% select(patients) %>% max()
  first_peak_data_fixed$timing_first[which(first_peak_data_fixed$replicate==rep&
                                             first_peak_data_fixed$metric=="ICU")] <- as.numeric(fixed_merge_check %>% filter(replicate==rep&metric=="ICU"&date_no<100) %>% mutate("true_first"=max(patients)) %>% filter(patients==max(patients)) %>% select(date_no) %>% min())
}

hosp_issues <- fixed_merge_check %>% filter(metric=="Hospitalisations") %>% select(replicate,peak,timing,peak_first,timing_first) %>% unique() %>% mutate("same"=peak-peak_first) %>% filter(same!=0) %>% select(replicate)

for(rep in hosp_issues$replicate){
  first_peak_data_fixed$peak_first[which(first_peak_data_fixed$replicate==rep&
                                           first_peak_data_fixed$metric=="Hospitalisations")] <- as.numeric(fixed_merge_check %>% filter(replicate==rep&metric=="Hospitalisations"&date_no<100) %>% select(patients) %>% max())
  first_peak_data_fixed$timing_first[which(first_peak_data_fixed$replicate==rep&
                                             first_peak_data_fixed$metric=="Hospitalisations")] <- as.numeric(fixed_merge_check %>% filter(replicate==rep&metric=="Hospitalisations"&date_no<100) %>% mutate("true_first"=max(patients)) %>% filter(patients==max(patients)) %>% select(date_no) %>% min())
}

death_issues <- fixed_merge_check %>% filter(metric=="Death") %>% select(replicate,peak,timing,peak_first,timing_first) %>% unique() %>% mutate("same"=peak-peak_first) %>% filter(same!=0) %>% select(replicate)

for(rep in death_issues$replicate){
  first_peak_data_fixed$peak_first[which(first_peak_data_fixed$replicate==rep&
                                           first_peak_data_fixed$metric=="Death")] <- as.numeric(fixed_merge_check %>% filter(replicate==rep&metric=="Death"&date_no<100) %>% select(patients) %>% max())
  first_peak_data_fixed$timing_first[which(first_peak_data_fixed$replicate==rep&
                                             first_peak_data_fixed$metric=="Death")] <- as.numeric(fixed_merge_check %>% filter(replicate==rep&metric=="Death"&date_no<100) %>% mutate("true_first"=max(patients)) %>% filter(patients==max(patients)) %>% select(date_no) %>% min())
}


# merge max peak df and first peak df
data_fixed <- merge(data_fixed,first_peak_data_fixed,by=c("metric","replicate"),all.x=TRUE)

# assign variable that will be used to colour

icu_rank_fixed <- data_fixed %>% filter(metric=="ICU")
icu_rank_fixed$rank <- match(rank(icu_rank_fixed$peak_first), sort(unique(rank(icu_rank_fixed$peak_first))))
hosp_rank_fixed <- data_fixed %>% filter(metric=="Hospitalisations")
hosp_rank_fixed$rank <- match(rank(hosp_rank_fixed$peak_first), sort(unique(rank(hosp_rank_fixed$peak_first))))
death_rank_fixed <- data_fixed %>% filter(metric=="Death")
death_rank_fixed$rank <- match(rank(death_rank_fixed$peak_first), sort(unique(rank(death_rank_fixed$peak_first))))

data_fixed <- rbind(icu_rank_fixed,hosp_rank_fixed,death_rank_fixed)

# median + quantiles 

data_fixed <- data.frame(data_fixed %>% group_by(date_no,metric) %>%  mutate(median = median(patients),
                                                                             q025 = quantile(patients,0.025),
                                                                             q25 = quantile(patients,0.25),
                                                                             q75 = quantile(patients,0.75),
                                                                             q975 = quantile(patients,0.975)))

# likelihood per simulation 
lls_fixed <- do.call(rbind,lapply(fixed$pmcmc_results$chains,  function(x) {x$results}))
lls_fixed <- lls_fixed$log_likelihood[match(rownames(fixed$replicate_parameters), rownames(lls_fixed))]

data_fixed$ll <- lls_fixed[match(data_fixed$replicate, 1:100)]

# set capacity limits
icu_capacity <- 1500
hosp_capacity <- 10000

fixed_capacity_calc_icu <- data.frame(data_fixed %>% filter(metric=="ICU") %>% group_by(replicate) %>% 
                                        mutate("capacity_breached"=any(patients>icu_capacity),
                                               "when_breached"=ifelse(capacity_breached,
                                                                      min(date_no[which(patients>icu_capacity)]),0),
                                               "how_long_breached"=ifelse(capacity_breached,
                                                                          length(which(patients>icu_capacity)),0)))

fixed_capacity_calc_hosp <- data.frame(data_fixed %>% filter(metric=="Hospitalisations") %>% group_by(replicate) %>% 
                                         mutate("capacity_breached"=any(patients>hosp_capacity),
                                                "when_breached"=ifelse(capacity_breached,
                                                                       min(date_no[which(patients>hosp_capacity)]),0),
                                                "how_long_breached"=ifelse(capacity_breached,
                                                                           length(which(patients>hosp_capacity)),0)))

fixed_capacity_calc <- rbind(fixed_capacity_calc_icu,fixed_capacity_calc_hosp)



### triggered scenario

triggered <- readRDS("Data/trigger_scenario.RDS")

# grab ICU data
hosp_cols <- grep("^IOx", colnames(triggered$output))
hosp_triggered <- do.call(
  rbind, 
  lapply(seq_len(dim(triggered$output)[3]), function(x) {
    data.frame(date = as.Date(rownames(triggered$output[,,x])),
               patients = rowSums(triggered$output[,hosp_cols,x]),
               replicate = x,
               metric="Hospitalisations") %>% filter(date>="2020-09-01")
  }))

# grab ICU data
imv_cols <- grep("^IMV", colnames(triggered$output))
icu_triggered <- do.call(
  rbind, 
  lapply(seq_len(dim(triggered$output)[3]), function(x) {
    data.frame(date = as.Date(rownames(triggered$output[,,x])),
               patients = rowSums(triggered$output[,imv_cols,x]),
               replicate = x,
               metric="ICU") %>% filter(date>="2020-09-01")
  }))

# grab deaths data
death_cols <- 461:477
deaths_triggered <- data.frame(do.call(
  rbind, 
  lapply(seq_len(dim(triggered$output)[3]), function(x) {
    data.frame(date = as.Date(rownames(triggered$output[,,x])),
               cumulative_deaths = rowSums(triggered$output[,death_cols,x]),
               replicate = x) %>% group_by(replicate) %>% 
      mutate("patients"=c(NA,diff(cumulative_deaths)),
             "metric"="Death") 
  }))) %>% select(date,patients,replicate,metric) %>% filter(date>="2020-09-01")


# bind three metrics together 
data_triggered <- rbind(hosp_triggered,icu_triggered,deaths_triggered)

# numeric date sequence for new truncated date list
data_triggered$date_no <- rep(rep(seq(from=0,by=1,length.out=nrow(deaths_triggered%>%filter(replicate==1))),100),3)

# gives peak and time for the max peak
data_triggered <- data_triggered %>% group_by(metric,replicate) %>% 
  mutate(peak = max(patients),timing = which.max(patients)-1)

# first peak of each trajectory for triggered scenarios

first_peak_data_triggered <- data.frame()

# extract first peak
for(i in 1:100){
  
  rep_hosp <- data_triggered %>% filter(replicate==i & metric=="Hospitalisations")
  rep_hosp_smooth <- predict(loess(rep_hosp$patients~rep_hosp$date_no))
  first_peak_hosp <- findpeaks(rep_hosp_smooth)[1,2]
  rep_time_hosp <- rep_hosp %>% subset(date_no >= first_peak_hosp-5&date_no<=first_peak_hosp+5)
  true_first_peak_hosp <- max(rep_time_hosp$patients)
  true_timing_hosp <- min(rep_time_hosp$date_no[which(rep_time_hosp$patients==true_first_peak_hosp)])
  
  rep_hosp_data <- data.frame("metric"="Hospitalisations","replicate"=i,
                              "peak_first"=true_first_peak_hosp,"timing_first"=true_timing_hosp)
  
  rep_icu <- data_triggered %>% filter(replicate==i & metric=="ICU")
  rep_icu_smooth <- predict(loess(rep_icu$patients~rep_icu$date_no))
  first_peak_icu <- findpeaks(rep_icu_smooth)[1,2]
  rep_time_icu <- rep_icu %>% subset(date_no >= first_peak_icu-5&date_no<=first_peak_icu+5)
  true_first_peak_icu <- max(rep_time_icu$patients)
  true_timing_icu <- min(rep_time_icu$date_no[which(rep_time_icu$patients==true_first_peak_icu)])
  
  rep_icu_data <- data.frame("metric"="ICU","replicate"=i,
                             "peak_first"=true_first_peak_icu,"timing_first"=true_timing_icu)
  
  
  rep_death <- data_triggered %>% filter(replicate==i & metric=="Death")
  rep_death_smooth <- predict(loess(rep_death$patients~rep_death$date_no))
  first_peak_death <- findpeaks(rep_death_smooth)[1,2]
  rep_time_death <- rep_death %>% subset(date_no >= first_peak_death-10&date_no<=first_peak_death+10)
  true_first_peak_death <- max(rep_time_death$patients)
  true_timing_death <- min(rep_time_death$date_no[which(rep_time_death$patients==true_first_peak_death)])
  
  rep_death_data <- data.frame("metric"="Death","replicate"=i,
                               "peak_first"=true_first_peak_death,"timing_first"=true_timing_death)
  
  ##want to rbind the individual reps then bind that with the bigger bind 
  first_peak_data_triggered_it <- rbind(rep_hosp_data,
                                        rep_icu_data,
                                        rep_death_data)
  
  first_peak_data_triggered <- rbind(first_peak_data_triggered,first_peak_data_triggered_it)
  
}

triggered_merge_check <- merge(data_triggered,first_peak_data_triggered,by=c("metric","replicate"),all.x=TRUE)

# manual fix for some replicates
icu_issues <- triggered_merge_check %>% filter(metric=="ICU") %>% select(replicate,peak,timing,peak_first,timing_first) %>% unique() %>% mutate("same"=peak-peak_first) %>% filter(same!=0) %>% select(replicate)

for(rep in icu_issues$replicate){
  first_peak_data_triggered$peak_first[which(first_peak_data_triggered$replicate==rep&
                                               first_peak_data_triggered$metric=="ICU")] <- triggered_merge_check %>% filter(replicate==rep&metric=="ICU"&date_no<100) %>% select(patients) %>% max()
  first_peak_data_triggered$timing_first[which(first_peak_data_triggered$replicate==rep&
                                                 first_peak_data_triggered$metric=="ICU")] <- as.numeric(triggered_merge_check %>% filter(replicate==rep&metric=="ICU"&date_no<100) %>% mutate("true_first"=max(patients)) %>% filter(patients==max(patients)) %>% select(date_no) %>% min())
}

hosp_issues <- triggered_merge_check %>% filter(metric=="Hospitalisations") %>% select(replicate,peak,timing,peak_first,timing_first) %>% unique() %>% mutate("same"=peak-peak_first) %>% filter(same!=0) %>% select(replicate)

for(rep in hosp_issues$replicate){
  first_peak_data_triggered$peak_first[which(first_peak_data_triggered$replicate==rep&
                                               first_peak_data_triggered$metric=="Hospitalisations")] <- as.numeric(triggered_merge_check %>% filter(replicate==rep&metric=="Hospitalisations"&date_no<100) %>% select(patients) %>% max())
  first_peak_data_triggered$timing_first[which(first_peak_data_triggered$replicate==rep&
                                                 first_peak_data_triggered$metric=="Hospitalisations")] <- as.numeric(triggered_merge_check %>% filter(replicate==rep&metric=="Hospitalisations"&date_no<100) %>% mutate("true_first"=max(patients)) %>% filter(patients==max(patients)) %>% select(date_no) %>% min())
}

death_issues <- triggered_merge_check %>% filter(metric=="Death") %>% select(replicate,peak,timing,peak_first,timing_first) %>% unique() %>% mutate("same"=peak-peak_first) %>% filter(same!=0) %>% select(replicate)

for(rep in death_issues$replicate){
  first_peak_data_triggered$peak_first[which(first_peak_data_triggered$replicate==rep&
                                               first_peak_data_triggered$metric=="Death")] <- as.numeric(triggered_merge_check %>% filter(replicate==rep&metric=="Death"&date_no<100) %>% select(patients) %>% max())
  first_peak_data_triggered$timing_first[which(first_peak_data_triggered$replicate==rep&
                                                 first_peak_data_triggered$metric=="Death")] <- as.numeric(triggered_merge_check %>% filter(replicate==rep&metric=="Death"&date_no<100) %>% mutate("true_first"=max(patients)) %>% filter(patients==max(patients)) %>% select(date_no) %>% min())
}



# merge max peak df and first peak df
data_triggered <- merge(data_triggered,first_peak_data_triggered,by=c("metric","replicate"),all.x=TRUE)

# assign variable that will be used to colour

icu_rank_triggered <- data_triggered %>% filter(metric=="ICU")
icu_rank_triggered$rank <- match(rank(icu_rank_triggered$peak_first), sort(unique(rank(icu_rank_triggered$peak_first))))
hosp_rank_triggered <- data_triggered %>% filter(metric=="Hospitalisations")
hosp_rank_triggered$rank <- match(rank(hosp_rank_triggered$peak), sort(unique(rank(hosp_rank_triggered$peak))))
death_rank_triggered <- data_triggered %>% filter(metric=="Death")
death_rank_triggered$rank <- match(rank(death_rank_triggered$peak_first), 
                                   sort(unique(rank(death_rank_triggered$peak_first))))

data_triggered <- rbind(icu_rank_triggered,hosp_rank_triggered,death_rank_triggered)

# median + quantiles 

data_triggered <- data.frame(data_triggered %>% group_by(date_no,metric) %>%  mutate(median = median(patients),
                                                                                     q025 = quantile(patients,0.025),
                                                                                     q25 = quantile(patients,0.25),
                                                                                     q75 = quantile(patients,0.75),
                                                                                     q975 = quantile(patients,0.975)))



# likelihood per simulation 
lls_triggered <- do.call(rbind,lapply(triggered$pmcmc_results$chains,  function(x) {x$results}))
lls_triggered <- lls_triggered$log_likelihood[match(rownames(triggered$replicate_parameters), rownames(lls_triggered))]

data_triggered$ll <- lls_triggered[match(data_triggered$replicate, 1:100)]


# arbitrary capacity limits set above

triggered_capacity_calc_icu <- data.frame(data_triggered %>% filter(metric=="ICU") %>% group_by(replicate) %>% 
                                            mutate("capacity_breached"=any(patients>icu_capacity),
                                                   "when_breached"=ifelse(capacity_breached,
                                                                          min(date_no[which(patients>icu_capacity)]),0),
                                                   "how_long_breached"=ifelse(capacity_breached,
                                                                              length(which(patients>icu_capacity)),0)))

triggered_capacity_calc_hosp <- data.frame(data_triggered %>% filter(metric=="Hospitalisations") %>% group_by(replicate) %>% 
                                             mutate("capacity_breached"=any(patients>hosp_capacity),
                                                    "when_breached"=ifelse(capacity_breached,
                                                                           min(date_no[which(patients>hosp_capacity)]),0),
                                                    "how_long_breached"=ifelse(capacity_breached,
                                                                               length(which(patients>hosp_capacity)),0)))

triggered_capacity_calc <- rbind(triggered_capacity_calc_icu,triggered_capacity_calc_hosp)



# combine main data frame 
data_fixed$scenario <- fixed_capacity_calc$scenario <- "Fixed"
data_triggered$scenario <- triggered_capacity_calc$scenario <- "Triggered"


data_comb <- rbind(data_fixed,data_triggered)
data_comb$metric_label <- ifelse(data_comb$metric=="ICU","ICU Demand",
                                 ifelse(data_comb$metric=="Death","Daily Deaths","Hospital Demand"))

# combine capacity data frame

capacity_calc <- rbind(fixed_capacity_calc,triggered_capacity_calc) %>% filter(metric!="Death") %>%
  select(metric,replicate,rank,ll,capacity_breached,when_breached,how_long_breached,scenario) %>% 
  arrange(capacity_breached) %>% unique()

capacity_calc$breached_label <- factor(ifelse(capacity_calc$capacity_breached,"Breached","Not breached"),levels=c("Not breached","Breached"))

capacity_calc_melt <- melt(capacity_calc,measure.vars = c("when_breached","how_long_breached"))
capacity_calc_melt$variable_label <- factor(ifelse(capacity_calc_melt$variable=="when_breached","Day of capacity breach","Duration of capacity breach"))


# summarise the key metrics 
data_summary <- data.frame(data_comb %>% group_by(scenario,metric_label) 
                           %>% mutate("median_peak"=median(peak_first),
                                      "lower_peak"=quantile(peak_first,0.025),
                                      "upper_peak"=quantile(peak_first,0.975),
                                      "median_time"=median(timing_first),
                                      "lower_time"=quantile(timing_first,0.025),
                                      "upper_time"=quantile(timing_first,0.975))
                           %>% select(scenario,metric_label,median_peak,lower_peak,upper_peak,
                                      median_time,lower_time,upper_time) %>% unique())


data_summary$metric_label <- ifelse(data_summary$metric_label=="Daily Deaths","Deaths",data_summary$metric_label)
data_summary$height <- rep(c(1500/286,7977/286,1),2)#from fixed peaks

data_summary$upper_time[which(data_summary$scenario=="Fixed"&data_summary$metric_label=="Hospital Demand")] <- 73.8
data_summary$scenario_label <- ifelse(data_summary$scenario=="Fixed","Scheduled suppression measures","Reactive suppression measures")
data_summary$scenario_label <- factor(data_summary$scenario_label,
                                      levels=c("Scheduled suppression measures","Reactive suppression measures"))

#### now ready to plot

# Figure 2: Traditional presentation of epidemic forecasts under two suppression strategies

fig_2_a <- ggplot(data_comb %>% filter(metric=="ICU"&scenario=="Fixed"),aes(x=date_no,y=median))+
  geom_ribbon(aes(x=date_no,ymax=q975,ymin=q75,fill="95% quantile"),alpha=0.2)+
  geom_ribbon(aes(x=date_no,ymax=q75,ymin=q25,fill="50% quantile"),alpha=0.5)+
  geom_ribbon(aes(x=date_no,ymax=q25,ymin=q025,fill="95% quantile"),alpha=0.2)+
  scale_fill_manual(values=c(alpha("#2A788EFF",0.2),alpha("#2A788EFF",0.5),alpha("#2A788EFF",0.2)))+
  geom_line(lwd=1,aes(col="median"))+
  scale_color_manual(values=c("#2A788EFF"))+
  theme_bw()+
  labs(x="Time (Days)",y="ICU Demand",fill="",col="",tag="a")+
  ylim(c(0,4000))+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", color = "black"))

fig_2_b <- ggplot(data_comb %>% filter(metric=="ICU"&scenario=="Triggered"),aes(x=date_no,y=median))+
  geom_ribbon(aes(x=date_no,ymax=q975,ymin=q75,fill="95% quantile"),alpha=0.2)+
  geom_ribbon(aes(x=date_no,ymax=q75,ymin=q25,fill="50% quantile"),alpha=0.5)+
  geom_ribbon(aes(x=date_no,ymax=q25,ymin=q025,fill="95% quantile"),alpha=0.2)+
  scale_fill_manual(values=c(alpha("#440154FF",0.1),alpha("#440154FF",0.5),alpha("#440154FF",0.1)))+
  geom_line(lwd=1,aes(col="median"))+
  scale_color_manual(values=c("#440154FF"))+
  theme_bw()+
  labs(x="Time (Days)",y="ICU Demand",fill="",col="",tag="b")+
  ylim(c(0,4000))+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", color = "black"))


fig_2_c <- ggplot(data_summary,aes(x=median_time,y=median_peak,col=scenario_label))+
  geom_point()+
  geom_errorbarh(aes(xmin=lower_time,xmax=upper_time,height=height*20))+
  geom_errorbar(aes(x=median_time,ymin=lower_peak,ymax=upper_peak,width=1.5))+
  theme_bw()+
  scale_colour_manual(values=c("#2A788EFF","#440154FF"))+
  facet_wrap(~metric_label,scales="free",strip.position="top")+
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "white", color = "black"))+
  labs(x="Time (Days)",y="Estimated Daily Peak",col="",tag="c")

plot_grid(plot_grid(fig_2_a,fig_2_b),
          fig_2_c,nrow=2)


# Figure 3: Communicating uncertainty in epidemic forecasts of ICU demand under two suppression strategies

fig_3_a <- ggplot(data_comb %>% filter(metric=="ICU"&scenario=="Fixed"&peak_first<5000),
                  aes(x=date_no,y=patients,group=replicate,col=rank))+
  geom_line(alpha=0.5)+
  geom_line(aes(x=date_no,y=median), inherit.aes = FALSE, lwd = 2) +
  theme_bw()+
  scale_color_viridis_c(end = 0.8) +
  labs(x="Time (Days)",y="ICU Demand",tag="a")+
  theme(legend.position = "NULL")+ylim(c(0,4700))

fig_3_b <- ggplot(data_comb %>% filter(metric=="ICU"&scenario=="Triggered"),aes(x=date_no,y=patients,group=replicate,col=rank))+
  geom_line(alpha=0.5)+
  geom_line(aes(x=date_no,y=median), inherit.aes = FALSE, lwd = 2) +
  theme_bw()+
  scale_color_viridis_c(end = 0.8) +
  labs(x="Time (Days)",y="ICU Demand",tag="b")+
  theme(legend.position = "NULL")+ylim(c(0,4700))

fig_3_c <- ggplot(data_comb %>% filter(metric=="ICU"&scenario=="Fixed"&peak_first<5000),
                  aes(timing_first, peak_first, color = rank)) +
  geom_point() +
  labs(x="Day of First Peak",y="ICU Demand at First Peak",tag="c")+
  theme_bw() +
  scale_color_viridis_c(end = 0.8) +
  scale_size(range=c(0.5,3), breaks=quantile(data_comb$ll, seq(0,1,0.2))) +
  theme(legend.position = "none")

fig_3_d <- ggplot(data_comb %>% filter(metric=="ICU"&scenario=="Triggered"),aes(timing_first, peak_first, color = rank)) +
  geom_point() +
  labs(x="Day of First Peak",y="ICU Demand at First Peak",tag="d")+
  theme_bw() +
  scale_color_viridis_c(end = 0.8) +
  scale_size(range=c(0.5,3), breaks=quantile(data_comb$ll, seq(0,1,0.2))) +
  theme(legend.position = "none")

plot_grid(fig_3_a, fig_3_b, fig_3_c, fig_3_d,
          ncol = 2)


# Figure 4: Communicating uncertainty in epidemic forecasts of ICU demand in terms of ICU capacity.

fig_4_a <- ggplot(data_comb %>% filter(metric=="ICU"&scenario=="Fixed"&peak_first<5000),
                  aes(x=date_no,y=patients,group=replicate,col=rank))+
  geom_line(alpha=0.5)+
  geom_line(aes(x=date_no,y=median), inherit.aes = FALSE, lwd = 2) +
  theme_bw()+
  scale_color_viridis_c(end = 0.8) +
  labs(x="Time (Days)",y="ICU Demand",tag="a")+
  theme(legend.position = "NULL")+ylim(c(0,4700))+
  #annotate("text",x=100,y=4000,label="Capacity breached in \n50% of simulations")+
  geom_hline(yintercept=1500,lwd=1,linetype="dashed")


fig_4_b_data <- capacity_calc_melt %>% filter(scenario=="Fixed"&metric=="ICU")
fig_4_b_data$value[which(fig_4_b_data$value==0)] <- -10

fig_4_b <- ggplot(fig_4_b_data,
                  aes(x=variable_label,y=value,col=rank))+
  scale_color_viridis_c(end = 0.8)+
  geom_boxplot(data=filter(fig_4_b_data,value>0,metric=="ICU",scenario=="Fixed"),outlier.shape=NA,lwd=1)+
  geom_jitter(width=0.35)+
  coord_cartesian(ylim=c(1,100))+
  theme_bw()+
  labs(x=" ",y="Time (Days)",tag="b")+
  theme(legend.position = "NULL")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15))

fig_4_c <- ggplot(data_comb %>% filter(metric=="ICU"&scenario=="Triggered"),aes(x=date_no,y=patients,group=replicate,col=rank))+
  geom_line(alpha=0.5)+
  geom_line(aes(x=date_no,y=median), inherit.aes = FALSE, lwd = 2) +
  theme_bw()+
  scale_color_viridis_c(end = 0.8) +
  labs(x="Time (Days)",y="ICU Demand",tag="c")+
  theme(legend.position = "NULL")+ylim(c(0,4700))+
  #annotate("text",x=100,y=4000,label="Capacity breached in \n70% of simulations")+
  geom_hline(yintercept=1500,lwd=1,linetype="dashed")

fig_4_d_data <- capacity_calc_melt %>% filter(scenario=="Triggered"&metric=="ICU")
fig_4_d_data$value[which(fig_4_d_data$value==0)] <- -10

fig_4_d <- ggplot(fig_4_d_data,
                  aes(x=variable_label,y=value,col=rank))+
  scale_color_viridis_c(end = 0.8)+
  geom_boxplot(data=filter(fig_4_d_data,value>0,metric=="ICU",scenario=="Triggered"),outlier.shape=NA,lwd=1)+
  geom_jitter(width=0.35)+
  coord_cartesian(ylim=c(1,100))+
  theme_bw()+
  labs(x=" ",y="Time (Days)",tag="d")+
  theme(legend.position = "NULL")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15))

plot_grid(fig_4_a,fig_4_b,fig_4_c,fig_4_d,
          axis="bt",align="hv")


#Supplementary Figure 1: Traditional presentation of epidemic forecasts under two suppression strategies

supp_1_a <- ggplot(data_comb %>% filter(metric=="Hospitalisations"&scenario=="Fixed"),aes(x=date_no,y=median))+
  geom_ribbon(aes(x=date_no,ymax=q975,ymin=q75,fill="95% quantile"),alpha=0.2)+
  geom_ribbon(aes(x=date_no,ymax=q75,ymin=q25,fill="50% quantile"),alpha=0.5)+
  geom_ribbon(aes(x=date_no,ymax=q25,ymin=q025,fill="95% quantile"),alpha=0.2)+
  scale_fill_manual(values=c(alpha("#2A788EFF",0.2),alpha("#2A788EFF",0.5),alpha("#2A788EFF",0.2)))+
  geom_line(lwd=1,aes(col="median"))+
  scale_color_manual(values=c("#2A788EFF"))+
  theme_bw()+
  ylim(c(0,21000))+
  labs(x="Time (Days)",y="Hospital Demand",fill="",col="",tag="a")+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", color = "black"))

supp_1_b <- ggplot(data_comb %>% filter(metric=="Hospitalisations"&scenario=="Triggered"),aes(x=date_no,y=median))+
  geom_ribbon(aes(x=date_no,ymax=q975,ymin=q75,fill="95% quantile"),alpha=0.2)+
  geom_ribbon(aes(x=date_no,ymax=q75,ymin=q25,fill="50% quantile"),alpha=0.5)+
  geom_ribbon(aes(x=date_no,ymax=q25,ymin=q025,fill="95% quantile"),alpha=0.2)+
  scale_fill_manual(values=c(alpha("#440154FF",0.1),alpha("#440154FF",0.5),alpha("#440154FF",0.1)))+
  geom_line(lwd=1,aes(col="median"))+
  scale_color_manual(values=c("#440154FF"))+
  theme_bw()+
  ylim(c(0,21000))+
  labs(x="Time (Days)",y="Hospital Demand",fill="",col="",tag="b")+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", color = "black"))

supp_1_c <- ggplot(data_comb %>% filter(metric=="Death"&scenario=="Fixed"),aes(x=date_no,y=median))+
  geom_ribbon(aes(x=date_no,ymax=q975,ymin=q75,fill="95% quantile"),alpha=0.2)+
  geom_ribbon(aes(x=date_no,ymax=q75,ymin=q25,fill="50% quantile"),alpha=0.5)+
  geom_ribbon(aes(x=date_no,ymax=q25,ymin=q025,fill="95% quantile"),alpha=0.2)+
  scale_fill_manual(values=c(alpha("#2A788EFF",0.2),alpha("#2A788EFF",0.5),alpha("#2A788EFF",0.2)))+
  geom_line(lwd=1,aes(col="median"))+
  scale_color_manual(values=c("#2A788EFF"))+
  theme_bw()+
  ylim(c(0,750))+
  labs(x="Time (Days)",y="Deaths",fill="",col="",tag="c")+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", color = "black"))

supp_1_d <- ggplot(data_comb %>% filter(metric=="Death"&scenario=="Triggered"),aes(x=date_no,y=median))+
  geom_ribbon(aes(x=date_no,ymax=q975,ymin=q75,fill="95% quantile"),alpha=0.2)+
  geom_ribbon(aes(x=date_no,ymax=q75,ymin=q25,fill="50% quantile"),alpha=0.5)+
  geom_ribbon(aes(x=date_no,ymax=q25,ymin=q025,fill="95% quantile"),alpha=0.2)+
  scale_fill_manual(values=c(alpha("#440154FF",0.1),alpha("#440154FF",0.5),alpha("#440154FF",0.1)))+
  geom_line(lwd=1,aes(col="median"))+
  scale_color_manual(values=c("#440154FF"))+
  theme_bw()+
  ylim(c(0,750))+
  labs(x="Time (Days)",y="Deaths",fill="",col="",tag="d")+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", color = "black"))


plot_grid(supp_1_a,supp_1_b,supp_1_c,supp_1_d,
          nrow=2)



# Supplementary Figure 2: Communicating uncertainty in epidemic forecasts of hospital demand under two suppression strategies

supp_2_a <- ggplot(data_comb %>% filter(metric=="Hospitalisations"&scenario=="Fixed"&peak_first<30000),
                   aes(x=date_no,y=patients,group=replicate,col=rank))+
  geom_line(alpha=0.5)+
  geom_line(aes(x=date_no,y=median), inherit.aes = FALSE, lwd = 2) +
  theme_bw()+
  scale_color_viridis_c(end = 0.8) +
  labs(x="Time (Days)",y="Hospital Demand",tag="a")+
  theme(legend.position = "NULL")

supp_2_b <- ggplot(data_comb %>% filter(metric=="Hospitalisations"&scenario=="Triggered"),
                   aes(x=date_no,y=patients,group=replicate,col=rank))+
  geom_line(alpha=0.5)+
  geom_line(aes(x=date_no,y=median), inherit.aes = FALSE, lwd = 2) +
  theme_bw()+
  scale_color_viridis_c(end = 0.8) +
  labs(x="Time (Days)",y="Hospital Demand",tag="b")+
  theme(legend.position = "NULL")+ylim(c(0,25000))

supp_2_c <- ggplot(data_comb %>% filter(metric=="Hospitalisations"&scenario=="Fixed"&peak_first<30000),
                   aes(timing_first, peak_first, color = rank)) +
  geom_point() +
  labs(x="Day of First Peak",y="Hospital Demand at First Peak",tag="c")+
  theme_bw() +
  scale_color_viridis_c(end = 0.8) +
  scale_size(range=c(0.5,3), breaks=quantile(data_comb$ll, seq(0,1,0.2))) +
  theme(legend.position = "none")

supp_2_d <- ggplot(data_comb %>% filter(metric=="Hospitalisations"&scenario=="Triggered"),
                   aes(timing_first, peak_first, color = rank)) +
  geom_point() +
  labs(x="Day of First Peak",y="Hospital Demand at First Peak",tag="d")+
  theme_bw() +
  scale_color_viridis_c(end = 0.8) +
  scale_size(range=c(0.5,3), breaks=quantile(data_comb$ll, seq(0,1,0.2))) +
  theme(legend.position = "none")

plot_grid(supp_2_a, supp_2_b, supp_2_c, supp_2_d, 
          nrow = 2)


#Supplementary Figure 3: Linking epidemic forecasts of hospital demand under two suppression strategies to potential capacity breaches

supp_3_a <- ggplot(data_comb %>% filter(metric=="Hospitalisations"&scenario=="Fixed"&peak_first<30000),
                   aes(x=date_no,y=patients,group=replicate,col=rank))+
  geom_line(alpha=0.5)+
  geom_line(aes(x=date_no,y=median), inherit.aes = FALSE, lwd = 2) +
  theme_bw()+
  scale_color_viridis_c(end = 0.8) +
  labs(x="Time (Days)",y="Hospital Demand",tag="a")+
  theme(legend.position = "NULL")+
  #annotate("text",x=100,y=4000,label="Capacity breached in \n50% of simulations")+
  geom_hline(yintercept=10000,lwd=1,linetype="dashed")


supp_3_b_data <- capacity_calc_melt %>% filter(scenario=="Fixed"&metric=="Hospitalisations")
supp_3_b_data$value[which(supp_3_b_data$value==0)] <- -10

supp_3_b <- ggplot(supp_3_b_data,
                   aes(x=variable_label,y=value,col=rank))+
  scale_color_viridis_c(end = 0.8)+
  geom_boxplot(data=filter(supp_3_b_data,value>0,metric=="Hospitalisations",scenario=="Fixed"),outlier.shape=NA,lwd=1)+
  geom_jitter(width=0.35)+
  coord_cartesian(ylim=c(1,100))+
  theme_bw()+
  labs(x=" ",y="Time (Days)",tag="b")+
  theme(legend.position = "NULL")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15))

supp_3_c <- ggplot(data_comb %>% filter(metric=="Hospitalisations"&scenario=="Triggered"),
                   aes(x=date_no,y=patients,group=replicate,col=rank))+
  geom_line(alpha=0.5)+
  geom_line(aes(x=date_no,y=median), inherit.aes = FALSE, lwd = 2) +
  theme_bw()+
  scale_color_viridis_c(end = 0.8) +
  labs(x="Time (Days)",y="Hospital Demand",tag="c")+
  theme(legend.position = "NULL")+ylim(c(0,25000))+
  #annotate("text",x=100,y=4000,label="Capacity breached in \n70% of simulations")+
  geom_hline(yintercept=10000,lwd=1,linetype="dashed")

supp_3_d_data <- capacity_calc_melt %>% filter(scenario=="Triggered"&metric=="Hospitalisations")
supp_3_d_data$value[which(supp_3_d_data$value==0)] <- -10

supp_3_d <- ggplot(supp_3_d_data,
                   aes(x=variable_label,y=value,col=rank))+
  scale_color_viridis_c(end = 0.8)+
  geom_boxplot(data=filter(supp_3_d_data,value>0,metric=="Hospitalisations",scenario=="Triggered"),outlier.shape=NA,lwd=1)+
  geom_jitter(width=0.35)+
  coord_cartesian(ylim=c(1,100))+
  theme_bw()+
  labs(x=" ",y="Time (Days)",tag="d")+
  theme(legend.position = "NULL")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15))

plot_grid(supp_3_a,supp_3_b,supp_3_c,supp_3_d,
          axis="bt",align="hv")

# text for figures 

capacity_calc_melt %>% filter(metric=="ICU"&variable=="when_breached") %>% select(scenario,breached_label) %>% table()
capacity_calc_melt %>% filter(metric=="Hospitalisations"&variable=="when_breached") %>% 
  select(scenario,breached_label) %>% table()

capacity_calc_melt %>% filter(value!=0) %>% group_by(metric,scenario,variable) %>% mutate("median"=median(value)) %>% select(metric,scenario,variable_label,median) %>% unique() %>% arrange(scenario)

