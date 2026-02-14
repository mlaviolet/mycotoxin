# analysis

library(here)
library(tidyverse)
library(broom)
library(flextable)

theme_set(theme_classic())

load(here("data", "toxin_data.Rdata"))

# # coerce toxin type and food groups to factors
# main_data <- main_data |> 
#   mutate(across(c(Food_type, toxin_type), factor)) |> 
#   # change "Miscellanous" level to "Other"
#   mutate(Food_type = fct_other(Food_type, drop = "Miscellaneous")) |> 
#   mutate(toxin_type = fct_collapse(toxin_type,
#     Trichothecene = c("Trichothecene Type A", "Trichothecene Type B"),
#     Difuranocoumarin = c("Difuranocoumarin",
#                          "Difuranocoumarin xanthone precursor to aflatoxin")
#     ))

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
  count(toxin_abb) |> 
  arrange(-n)

# toxins not occurring in any product
setdiff(unique(main_data$toxin_abb), toxin_tally$toxin_abb)

# distribution of food types
main_data |> 
  select(ID, Food_type) |> 
  distinct() |> 
  count(Food_type)

# most prevalent toxin by food type
food_type_list <- main_data |> 
  filter(amount > 0, Food_type != "Miscellaneous") |> 
  group_by(Food_type) |>
  count(toxin) |> 
  arrange(Food_type, -n) |> 
  group_split() 

walk(food_type_list, \(x) print(x, n = 5))

# most prevalent toxin types
toxin_grouped <- main_data |> 
  mutate(toxin_type = factor(toxin_type)) |> 
  filter_out(amount == 0) |> 
  count(toxin_type, .drop = FALSE) |> 
  # summarize(n = sum(n))
  arrange(-n) |> 
  mutate(pct = 100 * n / sum(n))

# Graph of toxin types ----------------------------------------------------

# graph of toxin_types, collapsing less frequent
work_data <- main_data |> 
  # use cases where a toxin was detected
  filter(amount > 0) |> 
  # collapse trichothecene and difuranocoumarin groups
  mutate(
    toxin_grp = 
      fct_collapse(
        toxin_type,
        Trichothecene = c("Trichothecene Type A", 
                          "Trichothecene Type B"),
        Difuranocoumarin = c("Difuranocoumarin",
                             "Difuranocoumarin xanthone precursor to aflatoxin")
             )) |> 
  # keep top 5 groups and collapse others
  mutate(toxin_grp = fct_lump_n(toxin_grp, 5)) |> 
  # make "Miscellaneous" food group "Other"
  mutate(Food_type = fct_other(Food_type, drop = "Miscellaneous"))

# labels of groups with individual toxins
toxin_grp_lbl <- rev(
  c("Cyclic hexadepsipeptide\n(Beau, Enn_A, Enn_A1, Enn_B, Enn_B1)",
    "Dibenzopyrone\n(AOH, AME)",
    "Difuranocoumarin\n(AFB1, AFB2, AFG1, AFG2)",
    "Polyketide\n(CIT, FB1, FB2, FB3, GRI)",
    "Trichothecene\n(3_Ace, 15_Ace, DAS, DOM, DON,\nDON_3_Glu, FUS-X, HT_2, NEO, NIV, T2)",
    "Other\n(\u03b1-ZEA, \u03b2-ZEA, CPA, OTA, Rocq, STC, ZEA)"))

# \u03b1 is Unicode for lower-case alpha; \u03b2 is lower-case beta

# put counts in decreasing order for graphing
work_data |> 
  count(toxin_grp) |> 
  # arrange(-n) |> 
  # set up horizontal bar graph--CAN THIS BE COMBINED WITH PREVIOUS LINE?
  ggplot(aes(x = reorder(toxin_grp, n), y = n)) +
  # bar graph
  geom_col(fill = "grey60") + 
  # add counts as labels on ends of bars
  geom_text(aes(label = n), hjust = -0.5) + 
  labs(y = "Number of occurrences (n = 4,012)", x = "Toxin group") +
  ylim(0, 235) +
  scale_x_discrete(labels = toxin_grp_lbl) + 
  theme(axis.title.y = element_blank()) +
  # make graph horizontal
  coord_flip() 
ggsave(here("output", "graph_toxin.png"))

# Table of toxin groups by food type --------------------------------------  
x <- work_data |> 
  filter_out(amount == 0) |> 
  count(Food_type, toxin_grp) |> 
  arrange(-n) |> 
  print(n = Inf)

ft <- xtabs(~ toxin_grp + Food_type, data = subset(work_data, amount > 0))|> 
  epitools::table.margins() |> 
  as_tibble(rownames = "Toxin") |> 
  flextable::flextable()  |> 
  add_header_row(values = "Food", colwidths = 9) |> 
  align(i = 1, j = NULL, align = "center", part = "header")
# https://rdrr.io/cran/flextable/man/align.html  
# use save_as_docx() to save to Word  
save_as_docx(ft, path = here("output", "table_food-by-toxin.docx"))

# OK TO HERE --------------------------------------------------------------

x <- work_data |> 
  filter_out(amount == 0) |> 
  count(Food_type, toxin_grp, .drop = FALSE) |> 
  arrange(-n) |> 
  pivot_wider(names_from = Food_type, values_from = n)

# inner_join(products) |> 
# summarize(n = n(), any_toxin = sum(n_toxins > 0)) |> 
# mutate(no_toxin = n - any_toxin)

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
