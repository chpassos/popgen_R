# Admixing chromosomes in R
library(tidyverse)


# Building chromosomes
df <- data.frame("chromosome" = c(rep(c(1), 1000), rep(c(2), 1000)), 
                 "position" = c(1:1000, 1:1000),
                 "type" = c(rep(c("A"), 1000), rep(c("B"), 1000)))


types <- c("A", "B")
ancestry_choice <- sample(types, size = 1, prob = c(0.5, 0.5))


output <- data.frame(matrix(ncol=3, nrow=1000))
for(i in 1:1000){
  
  recombination <- sample(c(TRUE, FALSE), size = 1, prob = c(0.05, 1 - 0.05))
  
  if(recombination == FALSE){
    chr_pos <- filter(df, position == i & type %in% ancestry_choice)
    output[i, ] <- chr_pos
    output[i,]
    
  } else{
    ancestry_choice <- types[!types %in%ancestry_choice]
    chr_pos <- filter(df, position == i & type %in% ancestry_choice)
    output[i, ] <- chr_pos
    output[i,]
  }
  
}

output %>%
  select("position" = X2, "type" = X3) %>% 
  mutate("chromosome" = rep(1, 1000)) %>%
  ggplot(aes(chromosome, position, fill = type)) +
  geom_tile()
######################################################




## Aqui embaixo são minhas tentativas de fazer diploide, preciso pensar ainda
##############################
## EXAMPLE -- 
filter(df, type %in% ancestry_choice)
filter(df, !(type %in% ancestry_choice))
output <- data.frame(matrix(ncol=3, nrow=1000))
i <- 467
chr_pos <- filter(df, position == i & type == "B")
output[i, ] <- chr_pos
output[i,]
############################


### Diploid individual
# Parental individuals
df <- data.frame("chromosome" = c(rep(c("1"), 1000), rep(c("2"), 1000)), 
                 "position" = c(1:1000, 1:1000),
                 "cromatide1" = rep(1, 2000),
                 "cromatide2" = rep(2, 2000),
                 "type" = c(rep(c("A"), 1000), rep(c("B"), 1000))) %>%
  pivot_longer(c(cromatide1, cromatide2), names_to = "cromatide", values_to = "valor_cromatide")

df %>%
  ggplot(aes(chromosome, position, group = cromatide, fill = type)) +
  geom_tile(width=0.4, position = position_dodge(width=0.5))


# First generation
cromatides <- c("cromatide1", "cromatide2")

cromatide_choice <- sample(cromatides, size = 1, prob = c(0.5, 0.5))
A <- filter(df, cromatide %in% cromatide_choice & type == "A")

cromatide_choice <- sample(cromatides, size = 1, prob = c(0.5, 0.5))
B <- filter(df, cromatide %in% cromatide_choice & type == "B")

A %>% rbind(B) %>%
  select(position, type, cromatide) %>%
  mutate("chromosome" = rep(1, 2000))

df_ger1 <- df %>%
  add_rownames() %>% 
  select(rowname, chromosome, position, type) %>%
  mutate(rowname = as.numeric(rowname)) %>%
  mutate("cromatide" = case_when(
      rowname <= 1000 ~ "cromatide1",
      rowname > 1000 ~ "cromatide2"
      ))

df_ger1 %>% View()
  group_by(cromatide) %>%
  distinct(position, .keep_all = TRUE) %>% View()
  ggplot(aes(chromosome, position, group = cromatide, fill = type)) +
  geom_tile(position = position_dodge(width = 0.5))
