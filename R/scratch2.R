# debug tallying pairs of toxin groups

library(here)
library(tidyverse)

# load data
load(here("data", "toxin_data_2026-06-06.Rdata"))

toxin_grp_lbl <-
  c("Cyclic hexadepsipeptide\nBEA, EnnA, EnnA1, EnnB, EnnB1",
    "Dibenzopyrone\nAOH, AME",
    "Difuranocoumarin\nAFB1, AFB2, AFG1, AFG2",
    "Polyketide\nCIT, FB1, FB2, FB3, GRI",
    "Trichothecene\n3-Ace, 15-Ace, DAS, DOM, DON,\nDON-3-Glu, FUS-X, HT-2, NEO, NIV, T2",
    "Other\n\u03b1-ZEA, \u03b2-ZEA, CPA, OTA, Rocq, STC, ZEA") |> 
  rev()

# bar graph of detected and quantified
bar_data <- main_data |> 
  select(toxin_type, amount, detected) |> 
  mutate(quantified = if_else(amount > 0, "Quantified", "Not quantified")) |> 
  select(-amount) |> 
  filter_out(detected == "No") |> 
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
  select(-toxin_type) |> 
  pivot_longer(-toxin_grp) |> 
  filter_out(value == "Not quantified") |> 
  select(-name) |> 
  mutate(value = factor(value, levels = c("Detected", "Quantified")))
  
bar_data |> 
  ggplot() +
  aes(x = toxin_grp, fill = value) +
  labs(y = "Count", x = NULL) +
  geom_bar(position = position_dodge(reverse = TRUE)) +
  scale_fill_discrete(palette = c("#ef8a62", "#67a9cf")) +
  
  geom_text(
    stat = "count", # Tells ggplot to calculate counts for the text labels
    aes(label = after_stat(count)), # Pulls the calculated counts into the label
    position = position_dodge(width = 0.9, reverse = TRUE),
    hjust = -0.5) +
  scale_x_discrete(limits = rev, labels = toxin_grp_lbl) + 
  ylim(0, 300) +
  coord_flip() +
  theme(legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.1))


# THIS IS IT--TWEAK -------------------------------------------------------

bar_data |> 
  count(quantified)

  pivot_longer()
  
  count(toxin_type, detected, quantified)

bar_data |> 
  summarize(q = sum(n), .by = detected)
bar_data |> 
  summarize(q = sum(n), .by = qu)

            
            bar_data_q <- main_data |> 
  select(ID, toxin_type, amount) |> 
  mutate(quantified = if_else(amount > 0, "Quantified", "Not quantified")) |> 
  filter(quantified == "Quantified") |> 
  select(-amount)
# 426
  
bar_data_d <- main_data |> 
  select(ID, toxin_type, detected) |> 
  filter(detected == "Detected")
# 706

# function to get pairs of toxins
get_combn <- function(df) {
  unique_vals <- unique(as.character(df$toxin_grp))
  t(combn(unique_vals, 2))
  }

# data to extract pairs
df <- main_data |> 
  select(ID, toxin_type, amount) |> 
  filter_out(amount == 0) |> 
  # distinct(ID) # 103 unique IDs
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
  select(-toxin_type, -amount) |> 
  # keep top 5 groups and collapse others
  mutate(toxin_grp = fct_lump_n(toxin_grp, 5)) |> 
  # distinct(ID) # still 103 unique IDs
  group_by(ID) |> 
  group_split() |> 
  # length is 103
  # remove data frames with only one row
  keep(~ nrow(.x) > 1) |> 
  # length is 84
  keep(~ length(unique(.x$toxin_grp)) > 1)
  # length is 67

df2 <- map(df, get_combn) 
df3 <- do.call(rbind, df2) |> 
  as.data.frame()

df4 <- df3 |> 
  count(V1, V2)



# df3 <- map_dfr(df2)
# 
# df3 <- do.call(rbind, df2) |> # this works
#   set_names(c("grp1", "grp2"))
# 
# df3


ggplot(mpg, aes(x = class, fill = drv)) +
  geom_bar(position = "dodge")


bind_rows(df2) # not working, why?



# test function S251551306-237716
# NEED TO REMOVE DUPLICATES?

df[[84]] |> get_combn()

df |> get_combn()

# CHECK OUT [[23]]
# ONLY ONE GROUP

df[[23]] |> unique() |> length()
  

