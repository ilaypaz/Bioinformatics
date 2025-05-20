#Group members: Stephane Guillemot, Avery Allan-McKay, Ilay Paz

# Load necessary libraries
library(tidyverse)

# Problem 1: Define function to perform chi-square operations
chisq_operations <- function(observed, expected, alpha) {
  
  # Input validation
  if (!is.data.frame(observed) || !is.data.frame(expected)) {
    stop("Error: 'observed' and 'expected' must be data frames.")
  }
  if (ncol(observed) != ncol(expected) || nrow(observed) != nrow(expected)) {
    stop("Error: 'observed' and 'expected' must have the same dimensions.")
  }
  if (any(expected <= 0)) {
    stop("Error: 'expected' frequencies must be positive.")
  }
  if (any(observed < 0)) {
    stop("Error: 'observed' frequencies cannot be negative.")
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Error: 'alpha' must be a number between 0 and 1.")
  }
  
  # Calculate degrees of freedom and critical value
  df <- ncol(observed) - 1  # Degrees of freedom
  critv <- qchisq(1 - alpha, df = df)
  
  # Initialize results list
  results <- list()
  
  # Loop through each row in the observed data
  for (i in 1:nrow(observed)) {
    # Calculate chi-square value for this row
    chi_square_value <- sum((observed[i, ] - expected[i, ])^2 / expected[i, ])
    hypothesis_rejected <- chi_square_value > critv
    
    # Store the results for this iteration
    results[[i]] <- list(
      significance_level = alpha,
      chi_square_value = chi_square_value,
      degrees_of_freedom = df,
      hypothesis_rejected = hypothesis_rejected
    )
  }
  
  # Print results for each iteration
  for (i in 1:length(results)) {
    print(results[[i]])
  }
  
  return(results)
}

# Load population data from CSV
pop <- read.csv("pop_data_lab5_problem3.csv")


# Extract columns with observed frequencies for Problem 1
cols <- pop[, 3:(ncol(pop))]

# Set expected quantities based on the first row (assumed constant across iterations)
expda1 <- cols[2, "A1A1"]
expda2 <- cols[2, "A1A2"]
expdho <- cols[2, "A2A2"]

# Initialize observed and expected data
observed <- cols
expected <- data.frame(A1A1 = rep(expda1, nrow(observed)),
                       A1A2 = rep(expda2, nrow(observed)),
                       A2A2 = rep(expdho, nrow(observed)))

# Call chisq_operations for Problem 1
prob1result_chisq <- chisq_operations(observed, expected, 0.05)
print(prob1result_chisq)


# Problem 2: Define the function to simulate allele dynamics under positive selection of a recessive allele
# Define the function to simulate allele dynamics under positive selection of a recessive allele
simulate_recessive_selection <- function(init_q, s, alpha, generations, pop_size) {
  
  # Input validation
  if (!is.numeric(init_q) || init_q <= 0 || init_q >= 1) stop("Initial frequency of q must be between 0 and 1")
  if (!is.numeric(s) || s < 0) stop("Selection coefficient must be non-negative")
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("Alpha must be between 0 and 1")
  if (!is.numeric(generations) || generations <= 0) stop("Generations must be a positive integer")
  if (!is.numeric(pop_size) || pop_size <= 0) stop("Population size must be a positive integer")
  
  # Initialize values
  q <- init_q
  p <- 1 - q
  results <- data.frame(generation = 1:generations, p = NA, q = NA)
  equilibrium_gen <- NA
  
  # Loop through generations
  for (gen in 1:generations) {
    # Store frequencies
    results$q[gen] <- q
    results$p[gen] <- p
    
    # Calculate Δq
    delta_q <- (p * q^2 * s) / (1 + q^2 * s)
    q <- q + delta_q
    p <- 1 - q
    
    # Calculate expected and observed counts
    observed <- data.frame(A1A1 = round(p^2 * pop_size), A1A2 = round(2 * p * q * pop_size), A2A2 = round(q^2 * pop_size))
    expected <- data.frame(A1A1 = p^2 * pop_size, A1A2 = 2 * p * q * pop_size, A2A2 = q^2 * pop_size)
    
    # Check HW equilibrium
    chi_res <- chisq_operations(observed, expected, alpha)
    if (is.na(equilibrium_gen) && any(sapply(chi_res, `[[`, "hypothesis_rejected"))) {
      equilibrium_gen <- gen
    }
  }
  
  # Plot results
  plot <- ggplot(results, aes(x = generation)) +
    geom_line(aes(y = q, color = "q"), na.rm = TRUE) +
    geom_line(aes(y = p, color = "p"), na.rm = TRUE) +
    labs(title = "Allele Dynamics under Recessive Selection", y = "Frequency", x = "Generation") +
    theme_minimal()
  
  if (!is.na(equilibrium_gen)) {
    plot <- plot + geom_vline(xintercept = equilibrium_gen, linetype = "dashed", color = "red")
    print(paste("Equilibrium disrupted at generation:", equilibrium_gen))
  } else {
    print("Equilibrium not disrupted within given generations.")
  }
  
  print(plot)
  
  # Return results
  list(
    initial_q = init_q,
    selection_coefficient = s,
    population_size = pop_size,
    results = results,
    equilibrium_generation = equilibrium_gen,  # Added here for clarity
    plot = plot
  )
}

# Run the function with different population sizes
prob2result_100 <- simulate_recessive_selection(0.5,0.01, 0.05, 1000, 100)
prob2result_1000 <- simulate_recessive_selection(0.5,0.01, 0.05, 1000, 1000)
prob2result_10000 <- simulate_recessive_selection(0.5,0.01, 0.05, 1000, 10000)

# Print the generation of HW equilibrium disruption for each result
print(paste("For pop_size = 100, equilibrium disrupted at generation:", prob2result_100$equilibrium_generation))
print(paste("For pop_size = 1000, equilibrium disrupted at generation:", prob2result_1000$equilibrium_generation))
print(paste("For pop_size = 10000, equilibrium disrupted at generation:", prob2result_10000$equilibrium_generation))

# For us, the output plots did not change with population. We are unsure why. 

#Problem 3
# Function to simulate allele dynamics under positive selection for co-dominant alleles
simulate_codom_selection <- function(data, s = 0.01, alpha = 0.05) {
  
  # Input validation
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Error: 'data' must be a matrix or a dataframe.")
  }
  

  
  if (!is.numeric(s) || s < 0) stop("Selection coefficient must be non-negative")
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("Alpha must be between 0 and 1")
  
  # Initialize variables
  equilibrium_gen <- NA
  co_dominance <- FALSE
  chi_square_values <- list()
  
  # Loop through each generation
  for (i in 1:nrow(data)) {
    # Calculate allele frequencies
    total_pop <- sum(data[i, c("A1A1", "A1A2", "A2A2")])
    p <- (2 * data$A1A1[i] + data$A1A2[i]) / (2 * total_pop)
    q <- 1 - p
    
    # Calculate Δq using the formula for co-dominance
    delta_q <- (p * q * s) / (1 + 2 * q * s)
    q <- q + delta_q
    p <- 1 - q
    
    # Calculate observed and expected counts
    observed <- data.frame(
      A1A1 = data$A1A1[i],
      A1A2 = data$A1A2[i],
      A2A2 = data$A2A2[i]
    )
    expected <- data.frame(
      A1A1 = p^2 * total_pop,
      A1A2 = 2 * p * q * total_pop,
      A2A2 = q^2 * total_pop
    )
    
    # Perform chi-square test for HW equilibrium
    chi_res <- chisq_operations(observed, expected, alpha)
    chi_square_values[[i]] <- chi_res
    
    # Check if HW equilibrium is disrupted
    if (is.na(equilibrium_gen) && any(sapply(chi_res, `[[`, "hypothesis_rejected"))) {
      equilibrium_gen <- data$gen[i]
      co_dominance <- TRUE  # Assuming equilibrium disruption suggests co-dominance
      break
    }
  }
  
  # Print output
  if (!is.na(equilibrium_gen)) {
    print(paste("Equilibrium disrupted at generation:", equilibrium_gen))
    print("This is a co-dominance scenario.")
  } else {
    print("Population is in HW equilibrium and does not show co-dominance.")
  }
  
  # Return results
  list(
    input_data = data,
    selection_coefficient = s,
    equilibrium_generation = equilibrium_gen,
    chi_square_values = chi_square_values,
    co_dominance = co_dominance
  )
}

#Create a data frame with the file, and then run our function using it

readfile <- read.csv("pop_data_lab5_problem3.csv")
q3answer <- simulate_codom_selection(readfile)
