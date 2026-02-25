pacman::p_load(haven, dplyr, survey, ggplot2, pROC, broom, stringr, mfx, arm, boot, purrr)
gss_data <- read_dta("gss7222_r4.dta")

data<-gss_data %>% 
  dplyr::select(finalter, hrs1, degree, childs, marital, sex, class, earnrs)

data <- data |>
  rename(
    worsening_situation = finalter,
    hours_worked = hrs1,
    children = childs,
    social_class = class,
    earners=earnrs
  )

# If worsening_situation is 2 (worse), then it's coded as 1, else it will be coded as 0 (better and stayed the same)
data$worsening_situation<-ifelse(data$worsening_situation == 2, 1, 0)

data$hours_worked <- ifelse(is.na(data$hours_worked), 0,  # Assign 1 for NA (not in work)
                            ifelse(data$hours_worked <= 34, 1,  # Assign 2 for 34 hours or less
                                   ifelse(data$hours_worked <= 49, 2,  # Assign 3 for 35 to 49 hours
                                          3)))  # Assign 4 for over 50 hours

# If no degree (junior and less), coded as 0, else coded as 1 (bachelors and graduate)
data$degree<- ifelse(data$degree <= 2, 0, 1)

# If married, coded as 1, else coded as 0 (widowed, divorced, separated, never married)
data$marital <- ifelse(data$marital == 1, 1, 0) 

# Male = 0, Female = 1 
data$sex <- data$sex - 1     

# Shifting values for social class so it starts at 0
data$social_class <- data$social_class - 1     
# Removing a the respondent who said 'no social class'
data$social_class[data$social_class == 4] <- NA

# Removing NA for complete case analysis
data_clean<- na.omit(data)

# Turning independent variables into factors
factor_columns <- function(data, columns) {
  data %>%
    mutate(across(all_of(columns), as.factor))
}

f_cols <- c("children", "social_class", "degree", "marital", "sex", "hours_worked")
data_clean <- factor_columns(data_clean, f_cols)

numeric_columns <- function(data, columns) {
  data %>%
    mutate(across(all_of(columns), as.numeric))
}

n_cols <- c("earners", "children")
data_clean <- numeric_columns(data_clean, n_cols)

# Creating test and training data to evaluate model performance
set.seed(1)

training.rows <- sample(nrow(data_clean),(nrow(data_clean)/2))
train.data <- data_clean[training.rows,]
test.data <- data_clean[-training.rows,]

# Creating the binary logistic models
fit_logit_model <- function(formula, data) {
  model <- glm(formula, family = binomial(link = "logit"), data = data)
}

model_1 <- fit_logit_model(worsening_situation ~ degree + social_class, train.data)
model_2 <- fit_logit_model(worsening_situation ~ hours_worked + children + earners + marital, train.data)
model_3 <- fit_logit_model(worsening_situation ~ degree + sex + marital + social_class + children, train.data)

# Creating confusion matrix and classification metrics
evaluate_logit_model <- function(model, test_data, threshold = 0.2) {
  
  # Predicted probabilities
  predicted_probs <- predict(model, newdata = test_data, type = "response")
  
  # Convert to binary class predictions
  predicted_classes <- ifelse(predicted_probs > threshold, 1, 0)
  
  # Confusion matrix
  actual <- test_data$worsening_situation
  cm <- table(predicted_classes, actual)
  
  # Metrics
  error_rate <- round((cm[1,2] + cm[2,1]) / length(actual) * 100, 2)
  sensitivity <- round(cm[2,2] / (cm[1,2] + cm[2,2]) * 100, 2)
  specificity <- round(cm[1,1] / (cm[1,1] + cm[2,1]) * 100, 2)
  
  # ROC curve and AUC
  roc_obj <- roc(actual, predicted_probs)
  auc_value <- roc_obj$auc
  plot(roc_obj, legacy.axes = TRUE, main = "ROC Curve: Evaluating Classification Performance")
  
  # Output
  return(list(
    confusion_matrix = cm,
    error_rate = error_rate,
    sensitivity = sensitivity,
    specificity = specificity,
    auc = auc_value
  ))
}

results_1 <- evaluate_logit_model(model_1, test.data)
results_2 <- evaluate_logit_model(model_2, test.data)
results_3 <- evaluate_logit_model(model_3, test.data)

# Examining the distribution of predicted probabilities 
# As most lie between 0.1 and 0.3, 0.2 was chosen as the threshold
predicted_probs <- predict(model_3, type = "response")
hist(predicted_probs, breaks = 100)

# Tidy the model and include confidence intervals
plotdata_1 <- tidy(model_3, conf.int = TRUE)

# Filter for relevant variables
plotdata_1 <- plotdata_1 %>%
  filter(str_detect(term, "degree|sex|marital|social_class|children"))

# Rename term labels for better readability in the plot
plotdata_1 <- plotdata_1 %>%
  mutate(term = recode(term,
                       "degree1" = "Has degree",
                       "sex1"= "Female",
                       "marital1" = "Married",
                       "social_class1" = "Working class",
                       "social_class2" = "Lower middle",
                       "social_class3" = "Upper middle",
                       "children" = "Children"
  ))

# Create the plot
ggplot(plotdata_1, aes(x = estimate, y = reorder(term, estimate))) +
  geom_point(color = "red", size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 color = "red", height = 0.2, size = 1) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  xlab("Estimated Coefficient") +
  ylab("") +
  ggtitle("Effect Sizes and Significance of Predictors") +
  theme(axis.text = element_text(size = 10),
        plot.title = element_text(size = 12))

#explain how much the variables matter in explaining incidence of diabetes.
ame <- logitmfx(
  formula = worsening_situation ~ degree + sex + marital +  social_class+ children ,
  data = data_clean,
  atmean = FALSE
)

# CIPP

# Simulate coefficients from model
sims <- sim(model_3, n.sims = 1000)
coefs <- coef(sims)

#check order of variables
names(coef(model_3))

#Creating profiles to compare predicted probabilities to 
#since middle class had a larger AME than upper class, it was chosen despite seeming counter-intuitive
low_risk <- c(1,  # Intercept
              1,           # degree = yes
              0,            #sex = male    
              1,                # married = yes
              0, 1, 0,         # social_class = middle class
              0) # children = 0



moderate_risk <- c(1,  # Intercept
                   0,           # degree = no
                   0,            #sex = male    
                   1,                # married = yes
                   1, 0, 0,         # social_class = working class
                   3) # children = 3


high_risk <- c(1,  # Intercept
               0,           # degree = no
               1,            #sex = female    
               0,                # married = no
               0, 0, 0,         # social_class = lower class
               8) # children = 8+

# Calculating the predicted probabilities for each profile
predict_risk <- function(profile, coefs) {
  return(inv.logit(profile %*% t(coefs)))
}

# Apply to your profiles
low_risk <- predict_risk(low_risk, coefs)
moderate_risk <- predict_risk(moderate_risk, coefs)
high_risk <- predict_risk(high_risk, coefs)

# Calculating the change in predicted probabilities between different profiles
low_to_mod <- moderate_risk - low_risk
mod_to_high <- high_risk - moderate_risk
low_to_high <- high_risk - low_risk

# Summary statistics
# Mean = average CIPP for each transition
# Confidence intervals for CIPP

summarise_cipp <- function(diff_vector) {
  list(
    mean = mean(diff_vector),
    ci_lower = quantile(diff_vector, 0.025),
    ci_upper = quantile(diff_vector, 0.975)
  )
}

low_mod_summary <- summarise_cipp(low_to_mod)
mod_high_summary <- summarise_cipp(mod_to_high)
low_high_summary <- summarise_cipp(low_to_high)

# Construct the data frame
results <- data.frame(
  "Profile" = as.factor(c("Low vs Moderate", "Moderate vs High", "Low vs High")),
  "CIPP" = c(mean(low_to_mod), mean(mod_to_high), mean(low_to_high)),
  "ciLB" = c(quantile(low_to_mod, 0.025), quantile(mod_to_high, 0.025), quantile(low_to_high, 0.025)),
  "ciUB" = c(quantile(low_to_mod, 0.975), quantile(mod_to_high, 0.975), quantile(low_to_high, 0.975))
)

# Plot
ggplot(results, aes(x = CIPP, y = Profile)) +
  geom_point(col = "red", size = 2) +
  geom_errorbar(aes(xmin = ciLB, xmax = ciUB), col = "red", width = 0.1, linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  xlab("Predicted Difference in Probability") +
  ylab("") +
  ggtitle("Predicted Differences in Probability of Worsening Finances Across Profiles") +
  theme(axis.text = element_text(size = 10),
        plot.title = element_text(size = 13))