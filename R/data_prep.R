# script to test import

# questions:
# how to handle foods with same labels?

library(tidyverse)
library(readxl)
library(here)

# need dates in YYYY-MM-DD
# column C is blank
# still have columns beginning with numbers: "15_Ace", "3_Ace";
#   suggest "Ace_15" and "Ace_3"--will change to full names later
# "corn flour" file--all names alike; in general, how to treat duplicated names?

# linking table of food and classification
# only need "raw data" file

df <- read_xlsx("data-raw/Toddler study raw data 2025-12-05.xlsx",
                na = c("nd", "np", "dn"),
                # cells E4 and M40 are apparently typos
                col_types = c(rep("text", 3), rep("numeric", 34)),
                n_max = 119) |> 
  select(- `...3`) |> 
  rename(Ace_15 = `15_Ace`, Ace_3 = `3_Ace`) |> 
  mutate(across(everything(), ~ replace_na(., 0)))

toxins <- read_xlsx("data-raw/2025-12-05 data legend.xlsx") |> 
  mutate(Abbreviation = str_replace(Abbreviation, "15_Ace", "Ace_15"),
         Abbreviation = str_replace(Abbreviation, "3_Ace", "Ace_3"))

df_long <- df |> 
  pivot_longer(cols = 3:36, names_to = "toxin", values_to = "ug_kg")


# # list of files to import
# x <- list.files(path = "data-raw", pattern = "\\.xlsx$", full.names = TRUE) 
# # remove "data legend" file
# x <- x[str_detect(x, pattern = "legend", negate = TRUE)]

# x <- readxl::read_excel(here("data-raw", "breakfast cereals.xlsx")) |> 
#   # drop empty column 3
#   select(-`...3`) |> 
#   # change space in column name to underscore
#   rename_with(~ str_replace(.x, " ", "_")) |> 
#   # change "nd" to numeric 0
#   # one value shown as missing because cell is blank, probably should be 0
#   mutate(across(3:36, ~ str_replace(.x, "nd", "0"))) |> 
#   mutate(across(3:36, as.numeric)) |> 
#   # this part does the reshaping
#   pivot_longer(cols = 3:36, names_to = "toxin", values_to = "amount") |> 
#   arrange(Food_tested, Number)
#   
# # number of detectable toxins per product  
#   
# x |> 
#   # group_by(Food_tested)
#   summarize(.by = Food_tested, n = sum(amount > 0, na.rm = TRUE)) |> 
#   print(n = Inf)
# 
# writexl::write_xlsx(x, here("data", "breakfast-cereals-tidy.xlsx"))
  