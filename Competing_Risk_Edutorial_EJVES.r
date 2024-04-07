#### Part 1 - Simulate time-to-event data ####
# Set working directory
setwd("~/Desktop/") # Adjust the working directory as desired


# Install and load necessary packages if not already installed
if(!require(survival)) install.versions("survival", version = '3.5-7')
if(!require(cmprsk)) install.versions("cmprsk", version = '2.2-11')
if(!require(survminer)) install.versions("survminer", version = '0.4.9')

# Set seed for reproducibility of simulations
set.seed(333)

# Number of observations
n_obs <- 800

# Create group variable A and B 
group_variable <- rep(c("GroupA", "GroupB"), each = n_obs/2)

# Simulate survival times from a Weibull distribution
shape_parameter <- 2
scale_parameter <- 40
survival_times <- rweibull(n_obs, shape = shape_parameter, scale = scale_parameter)

# Create a status variable indicating event occurrence (1 for Occlusion, 0 for Censoring)
event_status <- sample(0:1, n_obs, replace = TRUE, prob = c(0.7, 0.3))

# Create a data frame with the simulated data
simulated_data <- data.frame(
  ID = 1:n_obs,
  Time = survival_times,
  Event = event_status,
  Group = group_variable
)

# Code Event variable as a factor
simulated_data$Event <- as.factor(simulated_data$Event)

# Assign levesl to variable
levels(simulated_data$Event) <- list("Occlusion"=1,
                                     "Censor"=0)

#### Part 2 - Kaplan Meier Estimator ####
# Calculate survival model for to plot KM-Estimators and alculate logrank for difference
fit_simulation <- survfit(Surv(simulated_data$Time, as.numeric(simulated_data$Event)) ~ Group, data = simulated_data)

# Obtain occlusion-free Survival at 5-years (=60 months) from the KM estimator
summary(fit_simulation, times=60) # Note: The occlusion rate = 1-"survival" as presented by the summary statistic.


# Plot and store KM eestimator with details
jpeg("KM_Estimators_simulated.jpeg", units="in", width=12, height=10, res=600)
ggsurvplot(fit_simulation, data = simulated_data, fun = "event", # Assign data; delete "fun = "event" to start curves from top left
           pval = TRUE, # plot p-value
           pval.method = TRUE, # plot test name
           conf.int = TRUE, # plots 95% confidence interval bands for each curve
           surv.scale = "percent", # can be changed to modified
           xlab = "Time in Months",
           ylab = "Patients with Occlusion",
           break.time.by = 6,
           xlim = c(0,72),# break X axis in time intervals by 30
           censor.size = 4, # customize theme with a grid for better readability 
           legend.title = "", # plots legend title
           legend.labs = c("GroupA", "GroupB"),
           legend = c(0.1,0.6),
           risk.table = "abs_pct",   # absolute number and percentage at risk
           cumevents = TRUE,
           cumcensor = TRUE,
           tables.height = 0.13,
           ncensor.plot.height = 0.5,
           tables.theme = theme_pubr(),
           palette = c("#212121", "#607D8B"),
           title = "Kaplan-Meier Estimators: Time to Occlusion by Group")
dev.off()

#### Part 3 - Competing Risk Analysis ####
# Lets assume, 120 patients were censored due to death and they were NOT evenly distributed between the groups

# Copy original simulation data to new dataframe
simulated_data2 <- simulated_data

# Assign number of observations to randomly change
n_to_change_A <- 100
n_to_change_B <- 20

# Randomly change the levels for n observations to become dead instead of censored for both groups
indices_to_change_A <- sample(which(simulated_data2$Event!="Occlusion" & 
                                      simulated_data2$Time < 100 & 
                                      simulated_data2$Group == "GroupA"), n_to_change_A)

indices_to_change_B <- sample(which(simulated_data2$Event!="Occlusion" & 
                                      simulated_data2$Time < 100 & 
                                      simulated_data2$Group == "GroupB"), n_to_change_B)

# Assign new level "Dead" to the Event variable
levels(simulated_data2$Event) <- c(levels(simulated_data2$Event), "Dead")
simulated_data2$Event[indices_to_change_A] <- "Dead"
simulated_data2$Event[indices_to_change_B] <- "Dead"

# Calculate competing risk model
fit_simulation_CumInc <- cuminc(simulated_data2$Time, as.factor(simulated_data2$Event), simulated_data2$Group, cencode="Censor", na.action=na.omit)

# Test difference in occlusion rate and mortlity between the groups
fit_simulation_CumInc[["Tests"]]

# 95% confidence intervals can be calculated from the provided estimates and the variance: 95%-CI = est +/- 1.96*sqrt(var)

# Get estimates at 5-years (=60 months)
timepoints(fit_simulation_CumInc,60)

# Plot cumulative incidence curves for occlusion and mortality per group
jpeg("Competing_Risk_Plot.jpeg", units="in", width=10, height=7, res=600)
ggcompetingrisks(fit_simulation_CumInc, palette = c("#212121","#D84315"),
                 title = "Freedom from Occlusion by Group with Death as Competing Risk",
                 legend = "top",
                 xlab="Time to Event (Months)",
                 ylab = "Cumulative Incidence of Events (Occlusion & Death)",
                 ylim=c(0,1),
                 conf.int = T,
                 xlim=c(0,72),
                 ggtheme = theme_pubr(),
                 multsumiple_panels = F)  +
  annotate("text", x=35, y = 0.7, label= "Occlusion p = .001; Mortality p < .001")
dev.off()


#### Part 4 - Extension to plot numbers at risk ####
# Additional packages required
if(!require(riskRegression)) install.versions("riskRegression", '2023.03.22')
if(!require(prodlim)) install.versions("prodlim", '2023.08.28')

CompRskAnalysis2 <- prodlim(Hist(Time, Event, cens.code="Censor") ~ Group, data = simulated_data2)

# Two separate figures are plotted (occlusion and death)
# Currently, there is no simple option available to plot number of events

jpeg("Competing_Risk_Plot_Occlusion_with_n_at_risk.jpeg", units="in", width=10, height=7, res=600)
plot(CompRskAnalysis2,
     cause = "Occlusion",
     xlim=c(0, 72),
     legend.x="topleft", # position of legend
     legend.cex=1.5, # font size of legend
     marktime = TRUE, # the curves are tick-marked at right censoring times by invoking the function markTime.
     legend.title="",
     atrisk.title="Number at Risk",
     axis2.at=seq(0,1,0.2),
     axis1.at=seq(0,72,6),
     background.horizontal=seq(0,1,0.2),
     axis2.las=2,                            # rotate labels of y-axis 
     percent = FALSE,
     confint = TRUE,
     atrisk.col="black",
     atrisk.at=seq(0,72,6),
     col = c("#212121", "#607D8B"),
     xlab="Time in Months",
     ylab="Cumulative Incidence Curve")
title(main = "Cumulative Incidence of Occlusion by Group")


text(25,0.85,adj=0,paste("Gray's test: p-value = ", round(fit_simulation_CumInc$Tests[1,2],3)), cex = 1.2)

dev.off()  # Close the jpeg device


jpeg("Competing_Risk_Plot_Deaths_with_n_at_risk.jpeg", units="in", width=10, height=7, res=600)
plot(CompRskAnalysis2,
     cause = "Dead",
     xlim=c(0, 72),
     legend.x="topleft", # position of legend
     legend.cex=1.5, # font size of legend
     marktime = TRUE, # the curves are tick-marked at right censoring times by invoking the function markTime.
     legend.title="", # no legend title
     atrisk.title="Number at Risk",
     axis2.at=seq(0,1,0.2),
     axis1.at=seq(0,72,6),
     background.horizontal=seq(0,1,0.2),
     axis2.las=2, # rotate labels of y-axis
     percent = FALSE,
     confint = TRUE,
     atrisk.col="black",
     atrisk.at=seq(0,72,6),
     col = c("#212121", "#607D8B"),
     xlab="Time in Months",
     ylab="Cumulative Incidence Curve")

text(25,0.85,adj=0,paste("Gray's test: p-value < .001"), cex = 1.2)
title(main = "Cumulative Incidence of Death by Group")

dev.off()  # Close the jpeg device