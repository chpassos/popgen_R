# Population Genetics 
## Concepts, functions, simulations
library(tidyverse)

# Hardy-Weinberg Equilibrium Lecture
# Function to generate a population
pop_generator <- function(n = 100, fa = 0.5){
  tibble(
    ind = 1:n,
    allele1 = sample(c("A", "B"), n, prob = c(fa, 1-fa), replace = TRUE),
    allele2 = sample(c("A", "B"), n, prob = c(fa, 1-fa), replace = TRUE)) |>
    unite(c(allele1, allele2), col = "geno", sep = "") |>
    mutate(geno = case_when(
      geno == "AB" ~ "AB",
      geno == "BA" ~ "AB",
      TRUE ~ as.character(geno)))
}

pop_teste <- pop_generator(n = 500, fa = 0.7)

# Genotype Frequencies
## Example: f(AA) = #NAA / Ntotal
geno_freq <- pop_teste |>
  group_by(geno) |>
  tally() |>
  mutate(prop = n/sum(n))
geno_freq

# Allele Frequencies
## Example: f(A) = fAA + 1/2*fAB
fa <- geno_freq[geno_freq$geno == "AA", ][["prop"]] + (geno_freq[geno_freq$geno == "AB", ][["prop"]] / 2)
fb <- geno_freq[geno_freq$geno == "BB", ][["prop"]] + (geno_freq[geno_freq$geno == "AB", ][["prop"]] / 2)
fa
fb

# Heterozygosity
## Example: H = 2*fA*fB ~ Equals the expected frequency of heterozygotes
H <- 2*fa*fb
H

# Hardy-Weinberg Equilibrium
## Premises:
### 1. No selection
### 2. No Migration
### 3. No Mutation
### 4. Infinite Pop. Size
### 5. Random Mating
## If then so, genotype frequencies can be determined by:
## p^2; 2*p*q; q^2 -> Those are the Expected Genotype Frequencies

tibble(
  genos = c("exp_AA", "exp_AB", "exp_BB",
            "obs_AA", "obs_AB", "obs_BB"),
  values = c(fa^2, 2*fa*fb, fb^2,
             geno_freq |> filter(geno == "AA") |> pull(prop),
             geno_freq |> filter(geno == "AB") |> pull(prop),
             geno_freq |> filter(geno == "BB") |> pull(prop)))

# Would you say that this population is on Hardy-Weinberg Equilibrium?


##########
# Genetic Drift Lecture

# Coin Toss function -> Proportion of "caras". 
## How many coins are tossed (coins)
## Times we toss coins = coins (replicates)
coin_toss <- function(coins = 10, replicates = 10000, coin_prob = 0.5){
  replicates |>
    rerun(sample(c("cara", "coroa"), coins, prob = c(coin_prob, 1 - coin_prob), replace = TRUE)) |>
    map_dbl(~sum(.x == "cara")/length(.x)) |>
    as_tibble()
}

# Testing out our function :)
coin_toss(coins = 100)
coin_toss(coins = 10) |> 
  ggplot(aes(value, y = after_stat(count / sum(count)))) + 
  geom_histogram() + 
  coord_cartesian(xlim = c(0,1), ylim = c(0, 0.4))
coin_toss(coins = 100) |> 
  ggplot(aes(value, y = after_stat(count / sum(count)))) + 
  geom_histogram() + 
  coord_cartesian(xlim = c(0,1), ylim = c(0, 0.4))
coin_toss(coins = 1000) |> 
  ggplot(aes(value, y = after_stat(count / sum(count)))) + 
  geom_histogram() + 
  coord_cartesian(xlim = c(0,1), ylim = c(0, 0.4))

################################################################################
# Wright-Fisher Model 
## Premises:
### 1. No Selection
### 2. No Mutation
### 3. No Migration
### 4. Finite population size N, with 2N gene copies
### 5. Discrete generations
## Then, following these rules, one generations will be composed of
## random sampling of the gene copies present in the previous generations

# Genetic Drift Simulation -> Appears to be working!
## gen -> how many generations
## n -> population size
## fa -> initial frequency of "A" allele
drift_sim <- function(gen = 100, n = 1000, fa = 0.5){
  drift <- vector(mode = "numeric", length = gen)
  drift[1] <- fa
  
  for(i in 2:gen){
    new_pop <- sample(c("A", "B"), n, prob = c(drift[i-1], 1 - drift[i-1]), replace = TRUE)
    drift[i] <- length(new_pop[new_pop == "A"])/length(new_pop)
    drift
  }
  drift
}

# Generating replicas of drift
## Adding a "times" argument, which is the number of replicas we're doing
reps_drift_sim <- function(times = 5, gen = 100, n = 100, fa = 0.5){
  reps_drift <- rerun(times, drift_sim(gen, n, fa)) |>
    flatten_dbl() |>
    as_tibble() |>
    mutate(generations = rep(1:gen, times = times)) |>
    mutate(replicas = rep(1:times, each = gen))
  
}

# Function to plot our replicates ~ this is Allele Frequency x generations
## Adding a "means" argument. Plot mean allele freq?
plot_drift <- function(seed = NULL, times = 5, gen = 100, n = 100, fa = 0.5, means = FALSE) {
  if(!is.null(seed)) set.seed(seed)
  
  reps_drift <- reps_drift_sim(times, gen, n, fa)
  
  if(means == FALSE){
    reps_drift |>
      ggplot(aes(generations, value)) +
      geom_line(aes(group = replicas), 
                color = "gray70",
                size = 0.7) +
      labs(x = "Generations", y = "Frequency of Allele A") +
      coord_cartesian(ylim = c(0,1)) +
      theme_bw()
    
  } else{
    reps_drift |>
      ggplot(aes(generations, value)) +
      geom_line(aes(group = replicas), 
                color = "gray70",
                size = 0.7) +
      geom_line(data = reps_drift |>
                  rename(freq = value) |>
                  group_by(generations) |>
                  summarise(mean_freq = mean(freq)),
                aes(generations, mean_freq), 
                color = "red",
                size = 0.9) +
      labs(x = "Generations", y = "Frequency of Allele A") +
      coord_cartesian(ylim = c(0,1)) +
      theme_bw()
    
  }
  
}

plot_drift(seed = 9, times = 5, gen = 100, n = 100, fa = 0.3)
plot_drift(seed = 9, times = 5, gen = 100, n = 100, fa = 0.3, means = TRUE)

# Heterozygosity over time
plot_heterozygosity_drift <- function(seed = NULL, times = 5, gen = 100, n = 100, fa = 0.5, means = FALSE){
  if(!is.null(seed)) set.seed(seed)
  reps_drift <- reps_drift_sim(times, gen, n, fa) 
  
  if(means == FALSE){
    reps_drift |>
      group_by(generations, replicas) |>
      summarise(h = 2 * value * (1 - value)) |>
      ungroup() |>
      ggplot(aes(generations, h)) +
      geom_line(aes(group = replicas),
                color = "gray70",
                size = 0.7) +
      labs(x = "Generations", y = "Heterozygosity") +
      coord_cartesian(ylim = c(0, 0.5)) +
      theme_bw()
    
  } else{
    reps_drift |>
      group_by(generations, replicas) |>
      summarise(h = 2 * value * (1 - value)) |>
      ungroup() |>
      ggplot(aes(generations, h)) +
      geom_line(aes(group = replicas),
                color = "gray70",
                size = 0.7) +
      geom_line(data = reps_drift |>
                  group_by(generations, replicas) |>
                  summarise(h = 2 * value * (1 - value)) |>
                  ungroup() |>
                  group_by(generations) |>
                  summarise(mean_h = mean(h)),
                aes(generations, mean_h), 
                color = "red",
                size = 0.9) +
      labs(x = "Generations", y = "Heterozygosity") +
      coord_cartesian(ylim = c(0, 0.5)) +
      theme_bw()
  }
  
}

plot_heterozygosity_drift(seed = 2)
plot_heterozygosity_drift(seed = 2, means = TRUE)


# Variance between populations over time
plot_var_drift <- function(seed = NULL, times = 5, gen = 100, n = 100, fa = 0.5){
  if(!is.null(seed)) set.seed(seed)
  
  reps_drift <- reps_drift_sim(times, gen, n, fa)
  
  reps_drift |>
    group_by(generations) |>
    summarise(pop_var = var(value)) |>
    ggplot(aes(generations, pop_var)) +
    geom_line(color = "gray70",
              size = 0.7) +
    labs(x = "Generations", y = "Variance") +
    theme_bw()
}

plot_var_drift()

################################################################################



################################################################################
# Wright-Fisher Model with Mutation
## Here we have the Wright-Fisher Model, but this time
## the allele a mutates into allele A with probability mu each generation.
## What is the expected allele frequency in generation t + 1 if 
## the allele in gen t is fA(t) ?

drift_mut_sim <- function(gen = 100, n = 100, fa = 0.2, muBA = 0.001, muAB = 0.001, both = FALSE){
  drift <- vector(mode = "numeric", length = gen)
  drift[1] <- fa
  
  if(both == FALSE){ # Only mutation B -> A
    for(i in 2:gen){
      new_pop <- sample(c("A", "B"),
                        n,
                        prob = c((drift[i-1] + (muBA * (1 - drift[i-1]))),
                                 1 - drift[i-1]),
                        replace = TRUE)
      drift[i] <- length(new_pop[new_pop == "A"])/length(new_pop)
      drift
    }
    
  } else { # Mutations going B -> A or A -> B
    for(i in 2:gen){
      new_pop <- sample(c("A", "B"),
                        n,
                        prob = c((drift[i-1] + (muBA * (1 - drift[i-1]))),
                                 (1 - drift[i-1] + (muAB * (drift[i-1])))),
                        replace = TRUE)
      drift[i] <- length(new_pop[new_pop == "A"])/length(new_pop)
      drift
    }
  }
  drift
}

# Generating replicas of drift with mutation
## Adding a "times" argument, which is the number of replicas we're doing
reps_drift_mut_sim <- function(times = 5, gen = 100, n = 100, fa = 0.2, muBA = 0.001, muAB = 0.001, both = FALSE){
  reps_drift <- rerun(times, drift_mut_sim(gen, n, fa, muBA, muAB, both)) |>
    flatten_dbl() |>
    as_tibble() |>
    mutate(generations = rep(1:gen, times = times)) |>
    mutate(replicas = rep(1:times, each = gen))
  
}

# Function to plot our replicates ~ this is Allele Frequency x generations
## Adding a "means" argument. Plot mean allele freq - yes or no?
plot_drift_mut <- function(seed = NULL, times = 5, gen = 100, n = 100, fa = 0.5, muBA = 0.001, muAB = 0.001, both = FALSE, means = FALSE) {
  if(!is.null(seed)) set.seed(seed)
  
  reps_drift <- reps_drift_mut_sim(times, gen, n, fa, muBA, muAB, both)
  
  if(means == FALSE){
    reps_drift |>
      ggplot(aes(generations, value)) +
      geom_line(aes(group = replicas), 
                color = "gray70",
                size = 0.7) +
      labs(x = "Generations", y = "Frequency of Allele A") +
      coord_cartesian(ylim = c(0,1)) +
      theme_bw()
    
  } else{
    reps_drift |>
      ggplot(aes(generations, value)) +
      geom_line(aes(group = replicas), 
                color = "gray70",
                size = 0.7) +
      geom_line(data = reps_drift |>
                  rename(freq = value) |>
                  group_by(generations) |>
                  summarise(mean_freq = mean(freq)),
                aes(generations, mean_freq), 
                color = "red",
                size = 0.9) +
      labs(x = "Generations", y = "Frequency of Allele A") +
      coord_cartesian(ylim = c(0,1)) +
      theme_bw()
    
  }
  
}

plot_drift_mut(seed = 2, times = 10, n = 100, fa = 0.1,
               gen = 100, muBA = 0.001, means = FALSE)
plot_drift_mut(seed = 2, times = 10, n = 100, fa = 0.1,
               gen = 100, muBA = 0.001, means = TRUE)
plot_drift_mut(seed = 2, times = 10, n = 100, fa = 0.1,
               gen = 100, muBA = 0.1, means = FALSE,
               muAB = 0.01, both = TRUE)
plot_drift_mut(seed = 2, times = 10, n = 100, fa = 0.1,
               gen = 100, muBA = 0.1, means = TRUE,
               muAB = 0.01, both = TRUE)
plot_drift_mut(seed = 2, times = 10, n = 100, fa = 0.1,
               gen = 100, muBA = 0.01, means = FALSE,
               muAB = 0.1, both = TRUE)
plot_drift_mut(seed = 2, times = 10, n = 100, fa = 0.1,
               gen = 100, muBA = 0.01, means = TRUE,
               muAB = 0.1, both = TRUE)

# Heterozygosity over time
plot_heterozygosity_drift_mut <- function(seed = NULL, times = 5, gen = 100, n = 100, fa = 0.5, muBA = 0.001, muAB = 0.001, both = FALSE, means = FALSE){
  if(!is.null(seed)) set.seed(seed)
  reps_drift <- reps_drift_mut_sim(times, gen, n, fa, muBA, muAB, both) 
  
  if(means == FALSE){
    reps_drift |>
      group_by(generations, replicas) |>
      summarise(h = 2 * value * (1 - value)) |>
      ungroup() |>
      ggplot(aes(generations, h)) +
      geom_line(aes(group = replicas),
                color = "gray70",
                size = 0.7) +
      labs(x = "Generations", y = "Heterozygosity") +
      coord_cartesian(ylim = c(0, 0.5)) +
      theme_bw()
    
  } else{
    reps_drift |>
      group_by(generations, replicas) |>
      summarise(h = 2 * value * (1 - value)) |>
      ungroup() |>
      ggplot(aes(generations, h)) +
      geom_line(aes(group = replicas),
                color = "gray70",
                size = 0.7) +
      geom_line(data = reps_drift |>
                  group_by(generations, replicas) |>
                  summarise(h = 2 * value * (1 - value)) |>
                  ungroup() |>
                  group_by(generations) |>
                  summarise(mean_h = mean(h)),
                aes(generations, mean_h), 
                color = "red",
                size = 0.9) +
      labs(x = "Generations", y = "Heterozygosity") +
      coord_cartesian(ylim = c(0, 0.5)) +
      theme_bw()
  }
  
}

plot_heterozygosity_drift_mut(seed = 2)
plot_heterozygosity_drift_mut(seed = 2, means = TRUE)


# Variance between populations over time
plot_var_drift_mut <- function(seed = NULL, times = 5, gen = 100, n = 100, fa = 0.5, muBA = 0.001, muAB = 0.001, both = FALSE, means = FALSE){
  if(!is.null(seed)) set.seed(seed)
  
  reps_drift <- reps_drift_mut_sim(times, gen, n, fa, muBA, muAB, both)
  
  reps_drift |>
    group_by(generations) |>
    summarise(pop_var = var(value)) |>
    ggplot(aes(generations, pop_var)) +
    geom_line(color = "gray70",
              size = 0.7) +
    labs(x = "Generations", y = "Variance") +
    theme_bw()
}

plot_var_drift_mut()
################################################################################




################################################################################
# Natural Selection
## Premises:
### 1. Variation in the population
### 2. This variation contributes to survival and differential reproduction
### 3. This variation is heritable
### 4. Infinite Pop size
### 5. Random Mating
### 6. No mutation; No migration

# Function that simulates natural selection
## gen -> How many generations for selection to work
## fA -> Initial frequency of the advantageous allele
## s -> selection coefficient (how strong selection is)
## type -> types of selection
### These types includes dominant, recessive, 
### aditive, advantage of heterozygous and disadvantage of heterozygous
selection_sim <- function(gen = 100, fA = 0.5, s, type){
  
  if (type == "dominant") {
    adapt_values <- c("WAA" = 1, "WAa" = 1, "Waa" = 1-s)
    selection <- vector(mode = "numeric", length = gen)
    selection[1] <- fA
    
    for(i in 2:gen){
      fAA <- selection[i-1]^2
      fAa <- 2 * selection[i-1] * (1 - selection[i-1])
      faa <- (1 - selection[i-1])^2
      
      fAA_WAA <- fAA * adapt_values["WAA"]
      fAa_WAa <- fAa * adapt_values["WAa"]
      faa_Waa <- faa * adapt_values["Waa"]
      
      W_bar <- fAA_WAA + fAa_WAa + faa_Waa
      
      fA_next <- (fAA_WAA + (fAa_WAa/2)) / W_bar
      selection[i] <- fA_next
      
    }
    
    selection |>
      as_tibble() |>
      mutate(row_n = row_number()) |>
      ggplot(aes(row_n, value)) +
      geom_line(size = 0.7) +
      labs(x = "Generations", y = "Frequency of Advantageous Allele",
           title = paste0("Selection with s = ", s, 
                          "; type of selection: ", type)) +
      coord_cartesian(ylim = c(0, 1)) +
      theme_bw() 
    
    
  } else if (type == "recessive"){
    adapt_values <- c("WAA" = 1, "WAa" = 1-s, "Waa" = 1-s)
    selection <- vector(mode = "numeric", length = gen)
    selection[1] <- fA
    
    for(i in 2:gen){
      fAA <- selection[i-1]^2
      fAa <- 2 * selection[i-1] * (1 - selection[i-1])
      faa <- (1 - selection[i-1])^2
      
      fAA_WAA <- fAA * adapt_values["WAA"]
      fAa_WAa <- fAa * adapt_values["WAa"]
      faa_Waa <- faa * adapt_values["Waa"]
      
      W_bar <- fAA_WAA + fAa_WAa + faa_Waa
      
      fA_next <- (fAA_WAA + (fAa_WAa/2)) / W_bar
      selection[i] <- fA_next
      
    }
    
    selection |>
      as_tibble() |>
      mutate(row_n = row_number()) |>
      ggplot(aes(row_n, value)) +
      geom_line(size = 0.7) +
      labs(x = "Generations", y = "Frequency of Advantageous Allele",
           title = paste0("Selection with s = ", s, 
                          "; type of selection: ", type)) +
      coord_cartesian(ylim = c(0, 1)) +
      theme_bw() 
    
    
  } else if (type == "aditive"){
    adapt_values <- c("WAA" = 1, "WAa" = 1-(s/2), "Waa" = 1-s)
    selection <- vector(mode = "numeric", length = gen)
    selection[1] <- fA
    
    for(i in 2:gen){
      fAA <- selection[i-1]^2
      fAa <- 2 * selection[i-1] * (1 - selection[i-1])
      faa <- (1 - selection[i-1])^2
      
      fAA_WAA <- fAA * adapt_values["WAA"]
      fAa_WAa <- fAa * adapt_values["WAa"]
      faa_Waa <- faa * adapt_values["Waa"]
      
      W_bar <- fAA_WAA + fAa_WAa + faa_Waa
      
      fA_next <- (fAA_WAA + (fAa_WAa/2)) / W_bar
      selection[i] <- fA_next
      
    }
    
    selection |>
      as_tibble() |>
      mutate(row_n = row_number()) |>
      ggplot(aes(row_n, value)) +
      geom_line(size = 0.7) +
      labs(x = "Generations", y = "Frequency of Advantageous Allele",
           title = paste0("Selection with s = ", s, 
                          "; type of selection: ", type)) +
      coord_cartesian(ylim = c(0, 1)) +
      theme_bw() 
    
    
  } else if (type == "advantage of heterozygous"){
    adapt_values <- c("WAA" = 1-s, "WAa" = 1, "Waa" = 1-s)
    selection <- vector(mode = "numeric", length = gen)
    selection[1] <- fA
    
    for(i in 2:gen){
      fAA <- selection[i-1]^2
      fAa <- 2 * selection[i-1] * (1 - selection[i-1])
      faa <- (1 - selection[i-1])^2
      
      fAA_WAA <- fAA * adapt_values["WAA"]
      fAa_WAa <- fAa * adapt_values["WAa"]
      faa_Waa <- faa * adapt_values["Waa"]
      
      W_bar <- fAA_WAA + fAa_WAa + faa_Waa
      
      fA_next <- (fAA_WAA + (fAa_WAa/2)) / W_bar
      selection[i] <- fA_next
      
    }
    
    selection |>
      as_tibble() |>
      mutate(row_n = row_number()) |>
      ggplot(aes(row_n, value)) +
      geom_line(size = 0.7) +
      labs(x = "Generations", y = "Frequency of Advantageous Allele",
           title = paste0("Selection with s = ", s, 
                          "; type of selection: ", type)) +
      coord_cartesian(ylim = c(0, 1)) +
      theme_bw() 
    
    
  } else if (type == "disadvantage of heterozygous"){
    adapt_values <- c("WAA" = 1, "WAa" = 1-s, "Waa" = 1)
    selection <- vector(mode = "numeric", length = gen)
    selection[1] <- fA
    
    for(i in 2:gen){
      fAA <- selection[i-1]^2
      fAa <- 2 * selection[i-1] * (1 - selection[i-1])
      faa <- (1 - selection[i-1])^2
      
      fAA_WAA <- fAA * adapt_values["WAA"]
      fAa_WAa <- fAa * adapt_values["WAa"]
      faa_Waa <- faa * adapt_values["Waa"]
      
      W_bar <- fAA_WAA + fAa_WAa + faa_Waa
      
      fA_next <- (fAA_WAA + (fAa_WAa/2)) / W_bar
      selection[i] <- fA_next
      
    }
    
    selection |>
      as_tibble() |>
      mutate(row_n = row_number()) |>
      ggplot(aes(row_n, value)) +
      geom_line(size = 0.7) +
      labs(x = "Generations", y = "Frequency of Advantageous Allele",
           title = paste0("Selection with s = ", s, 
                          "; type of selection: ", type)) +
      coord_cartesian(ylim = c(0, 1)) +
      theme_bw() 
    
    
  } else {
    print("I can't deal with this type of selection")
  }
  
}

selection_sim(gen = 100, fA = 0.1, s = 0.1, type = "recessive")

################################################################################
# Interaction between Selection and Drift


