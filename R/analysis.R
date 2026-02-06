# analysis

library(here)
library(tidyverse)
library(broom)

load(here("data", "toxin_data.Rdata"))

# all products ------------------------------------------------------------

# number of toxins out of 34 in each of the 118 samples tested
number_toxins <- main_data |> 
  summarize(n = n(), n_toxins = sum(amount > 0), .by = ID)

# distribution of number of toxins
number_toxins |> 
  count(n_toxins)

range(number_toxins$n_toxins)
quantile(number_toxins$n_toxins, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
summary(number_toxins$n_toxins)

# Clopper-Pearson interval
n <- nrow(number_toxins)
# products with at least one toxin detected
x <- sum(number_toxins$n_toxins != 0)
x/n
binom.test(x, n) |> tidy()
# https://rpubs.com/brouwern/binomialCI2
# https://www.itl.nist.gov/div898/handbook/prc/section2/prc241.htm

# Wilson interval
epitools::binom.wilson(x, n)

# by product type ---------------------------------------------------------
by_food_type <- number_toxins |> 
  inner_join(products) |> 
  group_by(Food_type) |> 
  summarize(n = n(), any_toxin = sum(n_toxins > 0)) |> 
  mutate(no_toxin = n - any_toxin)
  # mutate(pct = 100 * any_toxin / n)

# Test for differences by type
# make contingency table as matrix
z <- by_food_type |> 
  filter(Food_type != "Miscellaneous") |> 
  select(any_toxin, no_toxin) |> 
  as.matrix(dimnames = list(by_food_type$Food_type, names(by_food_type)))
# row names not coming out for some reason; add manually
rownames(z) <- by_food_type$Food_type[-3]

# chi-square approximation questionable, use simulated p-value
chisq.test(z, simulate.p.value = TRUE)

# most prevalent toxins ---------------------------------------------------
toxin_tally <- main_data |> 
  filter(amount > 0) |> 
  count(toxin) |> 
  arrange(-n)

# toxins not occurring in any product
setdiff(unique(main_data$toxin), toxin_tally$toxin)

main_data |> 
  select(ID, Food_type) |> 
  distinct() |> 
  count(Food_type)

# RESUME HERE -------------------------------------------------------------

# most prevalent toxin by food type
food_type_list <- main_data |> 
  filter(amount > 0, Food_type != "Miscellaneous") |> 
  group_by(Food_type) |>
  count(toxin) |> 
  arrange(Food_type, -n) |> 
  group_split() 

walk(food_type_list, \(x) print(x, n = 5))
  
  
  
  inner_join(products) |> 
  summarize(n = n(), any_toxin = sum(n_toxins > 0)) |> 
  mutate(no_toxin = n - any_toxin)

# main_data |> 
#   group_by(ID, Food_type) |> 
#   count(detected = amount > 0) |> 
#   pivot_wider(id_cols = Food_type, names_from = detected, values_from = n) |>  
#   mutate(across(c(`FALSE`, `TRUE`), ~ replace_na(., 0)),
#          n = `FALSE` + `TRUE`) |> 
#   select(Food_type, toxins = `TRUE`, n)

# other analyses
#   at least one toxin by food type
#   specific toxins


# number_toxins <- main_data |> 
#   summarize(.by = ID, detected = amount > 0) |> 
#   ungroup() |> 
#   pivot_wider(id_cols = ID, names_from = detected, values_from = n) |>  
#   mutate(across(c(`FALSE`, `TRUE`), ~ replace_na(., 0)),
#          n = `FALSE` + `TRUE`) |> 
#   select(toxins = `TRUE`, n)

# number_toxins <- main_data |> 
#   group_by(ID) |> 
#   count(detected = amount > 0) |> 
#   ungroup() |> 
#   pivot_wider(id_cols = ID, names_from = detected, values_from = n) |>  
#   mutate(across(c(`FALSE`, `TRUE`), ~ replace_na(., 0)),
#          n = `FALSE` + `TRUE`) |> 
#   select(ID, toxins = `TRUE`, n)  

# # number of products tested

# main_data |> 
#   group_by(Food_type, ID) |> 
#   summarize(n = n(), n_toxins = sum(amount > 0))
