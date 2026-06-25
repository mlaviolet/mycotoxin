# debug tallying pairs of toxin groups

library(here)
library(tidyverse)

# load data
load(here("data", "toxin_data_2026-06-06.Rdata"))

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
  pivot_longer(-toxin_grp)

bar_data |> 
  filter_out(value == "Not quantified") |> 
  ggplot(aes(x = toxin_grp, fill = value)) +
  geom_bar(position = "dodge") +
  coord_flip()
# THIS IS IT--TWEAK

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
  

