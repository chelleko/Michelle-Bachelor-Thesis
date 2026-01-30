# Libraries
library(dplyr)
library(DEoptim)
library(psych)
library(ggplot2)

# Citations
citation("dplyr")
citation("DEoptim")
citation("psych")
citation("ggplot2")

# Loading my data set
my_data <- read.csv("data.csv")
mean(my_data$age)

# Quick sanity check
min(my_data$valence)
max(my_data$valence)
# The range is from 0 to 1

min(my_data$outcome)
max(my_data$outcome)
# The range is from -0.5 to 0.5

############### Objective Function ################

objective_function <- function(parameter,yt,xt) {

  # Define parameters
  alpha <- parameter[1] # Baseline level
  beta <- parameter[2] # Scaling factor
  gamma <- parameter[3] # Forgetting factor (0<gamma<1)
  
  n <- length(yt)
  
  # Discounted sum
  
  S <- numeric(n) # Creates a vector of 140 zeros
  
  S[1] <- xt[1]
  for (t in 2:n) {
    S[t] <- xt[t] + gamma * S[t - 1]
  }
  
  # Predicted values of y
  yhat <- alpha + beta * S
  
  # SSE
  SSE  <- sum((yt - yhat)^2)
  
  return(SSE)
}

############ Fit model for one participant +##############

fit_participant <- function(df) {
  
  res <- DEoptim(
    fn    = objective_function,
    lower = c(-2, -5, 0),
    upper = c( 2,  5, 1),
    yt    = df$valence,
    xt    = df$outcome,
    control = DEoptim.control(itermax = 200, trace = FALSE)
  )
  
  tibble(
    alpha = res$optim$bestmem[1],
    beta  = res$optim$bestmem[2],
    gamma = res$optim$bestmem[3],
    SSE   = res$optim$bestval
  )
}

########## # Fit model for all participants #############

# Now with all participants
results_all <- my_data %>%
  group_by(participant) %>%
  group_modify(~ fit_participant(.x)) %>%
  ungroup()

# Inspect results
head(results_all)
summary(results_all$gamma)
hist(results_all$gamma, breaks = 30)

# Save results
write.csv(results_all, "discounting_parameters_all_participants.csv",
          row.names = FALSE)

############# Plot Exponential Discounting Model ############

gamma <- 0.3

df <- data.frame(
  lag = 0:15
)
df$weight <- gamma^(df$lag)

ggplot(df, aes(x = lag, y = weight)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = bquote("Exponential discounting (" * gamma == .(gamma) * ")"),
    x = "Lag j (Trials ago)",
    y = expression("Weight" ~ gamma^j)
  ) +
  theme_minimal(base_size = 13)

################# CES-D & Cronbach's Alpha ################

# Restructuring data for the CES-D

# Selecting only the columns participant, valence, and CES-D 1-20
my_data_subset = subset(my_data, select = -c(age,gender,education,format,continuity,selected_door,valence,positive_affect,negative_affect,total_time,outcome,total,sequence,trial_number,trial,ease_1,ease_2,ease_3,ease_4,ease_5,ease_6,ease_7,ease_8,ease_9,ease_10,ease_11,ease_12,ease_13,ease_14,ease_15,ease_16,ease_17) )

# Selecting only the first rows for each participant's CES-D
my_data_cesd <- my_data_subset %>%
  distinct(participant, .keep_all = TRUE)

cesd_items <- subset(my_data_cesd, select = -participant)

# Converting the scale from 1-4 to 0-3
cesd_only <- cesd_items - 1

# Reverse scoring the positively worded items (4,8,12,16)
reverse_items <- c("cesd_4", "cesd_8", "cesd_12", "cesd_16")
cesd_only[reverse_items] <- 3 - cesd_only[reverse_items]

# Quick sanity check
range(as.matrix(cesd_only), na.rm = TRUE) # Range is from 0-3

alpha_results <- psych::alpha(cesd_only)
alpha_results$total$raw_alpha
alpha_results

# Create participant-level CES-D sum score (after reverse scoring)
cesd_sum <- rowSums(cesd_only, na.rm = TRUE)

# Participant-level dataframe: participant id + CES-D sum
cesd_participant_level <- data.frame(
  participant = my_data_cesd$participant,
  cesd_sum    = cesd_sum
)

# Quick sanity checks
View(cesd_participant_level)
nrow(cesd_participant_level)          
summary(cesd_participant_level$cesd_sum)
head(cesd_participant_level)

mean(cesd_sum)
range(cesd_sum)
sd(cesd_sum)

# Save as cesd.csv
write.csv(cesd_participant_level, "cesd.csv", row.names = FALSE)

################### Regression #######################

# Load discounting parameters
params <- read.csv("discounting_parameters_all_participants.csv")
range(params$gamma)
mean(params$gamma)
sd(params$gamma)


# Load CES-D scores
cesd <- read.csv("cesd.csv")
range(cesd$cesd_sum)

analysis_df <- params %>%
  left_join(cesd, by = "participant")

# Linear regression
lm_linear <- lm(cesd_sum ~ gamma, data = analysis_df)
summary(lm_linear)

# Quadratic regression
lm_quadratic <- lm(cesd_sum ~ gamma + I(gamma^2), data = analysis_df)
summary(lm_quadratic)

#################### Plotting ########################

# Histogram of cesd scores
hist(cesd$cesd_sum,xlab="CES-D scores", ylab="Participants",main="Histogram of CES-D Scores")

# Histogram of gamma scores
hist(params$gamma,xlab=expression(gamma), ylab="Participants",main="Histogram of γ values")

# Plotting cesd and gamma
ggplot(analysis_df, aes(x = gamma, y = cesd_sum)) +
  geom_point(alpha = 0.4, color = "gray") +
  geom_smooth(aes(color = "Linear"), method = "lm", formula = y ~ x, se = FALSE) +
  geom_smooth(aes(color = "Quadratic"), method = "lm", formula = y ~ x + I(x^2), se = FALSE, linetype = "dashed") +
  scale_color_manual(
    name = "Model",
    values = c("Linear" = "blue", "Quadratic" = "red")
  ) +
  labs(
    x = expression(gamma),
    y = "CES-D scores",
    title = "Relationships between γ and CES-D Scores"
  ) +
  theme_classic(base_size = 12)
