# Data for the months and outbreaks
months <- c("April", "May", "June", "July", "August", "September", "October", "November", "December", "January", "February", "March")
number_of_outbreaks <- c(45, 35, 43, 84, 60, 80, 64, 64, 39, 48, 33, 43)

# Constants
N <- 638 # Total outbreaks
p <- 1/12 # Probability for each month

# Function to calculate p-value using binom.test
calculate_p_values <- function(n, N, p) {
  test <- binom.test(n, N, p, alternative = "two.sided")
  prob <- dbinom(n,N,p)
  p_value <- test$p.value
  significant <- ifelse(p_value < 0.05, "Sig.", "")
  return(c(as.numeric(prob), as.numeric(p_value), significant))
}

# Apply the function to each month's data
results <- t(sapply(number_of_outbreaks, calculate_p_values, N = N, p = p))

# Convert the numeric parts back to numeric explicitly
results_numeric <- as.data.frame(results)
results_numeric[, 1] <- as.numeric(results_numeric[, 1])  # Probability
results_numeric[, 2] <- as.numeric(results_numeric[, 2])  # P-value

# Combine results into a data frame
output_table <- data.frame(
  Cluster = paste0("T", 1:12),
  Months_of_Cluster = months,
  Number_of_Outbreaks_in_Cluster = number_of_outbreaks,
  Total_Number_of_Outbreaks = N,
  Probability = round(results_numeric[, 1], digits = 5),
  P_value = signif(results_numeric[, 2], digits = 3),
  Significant_Result = results_numeric[, 3]
)

# Add July - September and July - November clusters
july_to_september_outbreaks <- sum(number_of_outbreaks[4:6])
july_to_november_outbreaks <- sum(number_of_outbreaks[4:8])

july_to_september <- calculate_p_values(july_to_september_outbreaks, N, 3/12)
july_to_november <- calculate_p_values(july_to_november_outbreaks, N, 5/12)

output_table <- rbind(
  output_table,
  data.frame(
    Cluster = c("T13", "T14"),
    Months_of_Cluster = c("July - September", "July - November"),
    Number_of_Outbreaks_in_Cluster = c(july_to_september_outbreaks, july_to_november_outbreaks),
    Total_Number_of_Outbreaks = N,
    Probability = round(c(as.numeric(july_to_september[1]), as.numeric(july_to_november[1])), digits = 5),
    P_value = signif(c(as.numeric(july_to_september[2]), as.numeric(july_to_november[2])), digits = 3),
    Significant_Result = c(july_to_september[3], july_to_november[3])
  )
)

# Print the table
print(output_table)



n = 224
N = 638
p = 3 / 12

#Q1.
dbinom(n,N,p)

#Is it significant?
#Q2.
(n/p-N+1)*dbinom(n,N,p) + 2 * (1-pbinom(n,N,p))

binom.test(n,N,p, alternative="two.sided")
prop.test(n,N,p)
