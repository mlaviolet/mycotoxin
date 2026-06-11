# analysis

library(here)
library(tidyverse)
library(broom)
library(flextable)
# library(widyr)

theme_set(theme_classic())

# load(here("data", "toxin_data.Rdata"))
load(here("data", "toxin_data_2026-06-06.Rdata"))

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
  summarize(n = n(), 
            n_toxins = sum(amount > 0), 
            .by = ID) |> 
  arrange(ID)

# distribution of number of toxins
number_toxins |> 
  count(n_toxins)

# 45 of 118 have two toxins or fewer

range(number_toxins$n_toxins)
quantile(number_toxins$n_toxins, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
summary(number_toxins$n_toxins)
# MEDIAN IS NOW 3, MEAN IS 3.6
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   1.000   3.000   3.627   6.000  13.000 

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
# x   n proportion    lower     upper conf.level
# 1 103 118  0.8728814 0.800821 0.9214291       0.95

# by product type ---------------------------------------------------------
by_food_type <- number_toxins |> 
  inner_join(products) |> 
  group_by(Food_type) |> 
  summarize(n = n(), any_toxin = sum(n_toxins > 0)) |> 
  mutate(no_toxin = n - any_toxin)
  # mutate(pct = 100 * any_toxin / n)

# Test for differences by type
# make contingency table as matrix
food_by_any_toxin <- by_food_type |> 
  # filter(Food_type != "Miscellaneous") |> 
  select(any_toxin, no_toxin) |> 
  as.matrix(dimnames = list(by_food_type$Food_type, names(by_food_type)))
# row names not coming out for some reason; add manually
rownames(food_by_any_toxin) <- by_food_type$Food_type

# chi-square approximation questionable, use simulated p-value
set.seed(42)
chisq.test(food_by_any_toxin, simulate.p.value = TRUE)
# P = 0.66

# most prevalent toxins ---------------------------------------------------
toxin_tally <- main_data |> 
  filter(amount > 0) |> 
  count(toxin_abb) |> 
  arrange(-n)

# toxins not occurring in any product
setdiff(unique(main_data$toxin_abb), toxin_tally$toxin_abb)
# [1] "Rocq"  "AFG2"  "DAS"   "FUS-X" "NEO"   "CIT"   "a-ZEA" "b-ZEA"

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
  mutate(toxin_grp = fct_lump_n(toxin_grp, 5)) 

# labels of groups with individual toxins
toxin_grp_lbl <-
  c("Cyclic hexadepsipeptide\n(BEA, EnnA, EnnA1, EnnB, EnnB1)",
    "Dibenzopyrone\n(AOH, AME)",
    "Difuranocoumarin\n(AFB1, AFB2, AFG1, AFG2)",
    "Polyketide\n(CIT, FB1, FB2, FB3, GRI)",
    "Trichothecene\n(3-Ace, 15-Ace, DAS, DOM, DON,\nDON-3-Glu, FUS-X, HT-2, NEO, NIV, T2)",
    "Other\n(\u03b1-ZEA, \u03b2-ZEA, CPA, OTA, Rocq, STC, ZEA)") |> 
  rev()

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
  labs(y = "Number of positive samples", 
       caption = "Number of products tested: 118") +
  ylim(0, 235) +
  scale_x_discrete(labels = toxin_grp_lbl) + 
  theme(axis.title.y = element_blank(),
        plot.caption = element_text(hjust = 0.5)) +
  # make graph horizontal
  coord_flip() 

ggsave(here("output",
            paste0("graph_toxin_", as.character(today()), ".png")))

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
# save_as_docx(ft, path = here("output", "table_food-by-toxin.docx"))

# crosstab of toxin group by food type
work_data |> 
  filter_out(amount == 0) |> 
  count(Food_type, toxin_grp, .drop = FALSE) |> 
  arrange(-n) |> 
  pivot_wider(names_from = Food_type, values_from = n)

# Table 4: Summary of samples contaminated --------------------------
# max and median contamination level 
# reproduce Table 4 in manuscript
step1 <- work_data |> 
  summarize(n = n(),
            max_amt = max(amount),
            med_amt = median(amount),
            .by = toxin_abb) |> 
  right_join(toxins, by = "toxin_abb") |> 
  mutate(n = replace_na(n, 0)) 

# Ace_3 has two with max 100; choose one at random
set.seed(42)
step2 <- work_data |> 
  filter(amount ==  max(amount), .by = toxin_abb) |> 
  # filter(toxin_abb == "Ace_3") |> 
  slice_sample(n = 1, by = toxin_abb) |> 
  select(toxin_abb, toxin_grp, amount, Food_type) 
  
table_4 <- left_join(step1, step2, by = "toxin_abb") |>
  select(starts_with("toxin"), LOQ, n, max_amt, Food_type, med_amt) 
levels(table_4$Food_type)[7] <- "Baby apple juice"
  
table_4 <- table_4 |> 
  arrange(toxin_grp, toxin)

writexl::write_xlsx(
  table_4, 
  here("output", 
       paste0("Table4_", as.character(today()), ".xlsx")
  ))


# OK TO HERE --------------------------------------------------------------


  # all "Other" are apple juice; rename
  # arrange(toxin_abb)


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

# number_toxins <- main_data |> n_levels <- nlevels(my_factor)

# Calculate combinations (n choose 2)
# choose(n_levels, 2)
#   ungroup() |> 
#   pivot_wider(id_cols = ID, names_from = detected, values_from = n) |>  
#   mutate(across(c(`FALSE`, `TRUE`), ~ replace_na(., 0)),
#          n = `FALSE` + `TRUE`) |> 
#   select(ID, toxins = `TRUE`, n)  

# # number of products tested

# main_data |> 
#   group_by(Food_type, ID) |> 
#   summarize(n = n(), n_toxins = sum(amount > 0))

# references for heatmaps
# https://davetang.org/muse/2010/12/06/making-a-heatmap-with-r/
# https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html

# NOT WORKING
unique_vals <- unique(work_data$toxin_grp)
df <- work_data |> 
  group_by(ID)
  mutate(x = combn(unique_vals, 2))

choose(length(unique_vals), 2)
x <- all_pairs <- combn(unique_vals, 2)
class(all_pairs) 

df <- work_data |> 
  widyr::pairwise_count(item = toxin_grp, feature = ID, diag = FALSE) |> 
  print(n = Inf)
  
df2 <- work_data %>%
  # Join data onto itself by the group identifier
  inner_join(work_data, by = "ID", relationship = "many-to-many") %>%
  # Filter to avoid matching an item with itself or counting pairs twice
  filter(toxin_grp.x < toxin_grp.y) %>%
  # Count the unique pairs per group
  count(ID, toxin_grp.x, toxin_grp.y, name = "pair_count")

 # NEED TO BUILD SIMPLE TEST DATA SET


