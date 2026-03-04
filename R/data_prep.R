# script to import and process data

library(tidyverse)
library(here)
library(readxl)
library(writexl)

# FIXED
# need dates in YYYY-MM-DD
# column C is blank
# still have columns beginning with numbers: "15_Ace", "3_Ace";
#   suggest "Ace_15" and "Ace_3"--will change to full names later
# only need "raw data" file

# QUESTIONS
# "corn flour" file--all names alike; in general, how to treat duplicated names?
# Coconut_Flour	S251511351-237971 *2 (why the "*2"?) FIXED
# in raw data spreadsheet cells E4 and M40 are apparently typos, confirm?
#   rows 25 and 29 have duplicate ID S250551125-223922 FIXED

# linking table of food and classification
products <- read_xlsx("data-raw/product-key.xlsx") |> 
  # cleanup
  mutate(ID = str_remove(ID, " \\*2")) |> 
  filter(Food_tested != "Organic_fruit_and_veggie_bars") |> 
  mutate(Food_type = fct_other(Food_type, drop = "Miscellaneous"))

# full name of toxin for reporting, and type of toxin
toxins <- read_xlsx("data-raw/2026-03-04 data legend.xlsx") |> 
  mutate(Abbreviation = str_replace(Abbreviation, "15_Ace", "Ace_15"),
         Abbreviation = str_replace(Abbreviation, "3_Ace", "Ace_3")) |> 
  rename(toxin_abb = Abbreviation, toxin = Mycotoxin,
         toxin_type = `Mycotoxin Type`) |> 
  mutate(toxin_type = factor(toxin_type)) |> 
  rename(LOQ = `LOQ (ppm)`)

# import data from Excel file
main_data <- read_xlsx("data-raw/Toddler study raw data 2025-12-05.xlsx",
                       na = c("nd", "np", "dn"),
                       # cells E4 and M40 are apparently typos
                       col_types = c(rep("text", 3), rep("numeric", 34)),
                       n_max = 119) |> 
  # delete empty column
  select(- `...3`) |> 
  # change names to syntactic
  rename(Ace_15 = `15_Ace`, Ace_3 = `3_Ace`, BEA = Beau) |> 
  # change names to match updated linking table of 2026-02-07
  rename(# BETA = Beau,
         `FUS-X` = FX,
         NEO = NEOS,
         ZEA = ZONE,
         STC = Sterig,
         GRI = GRIS,
         `a-ZEA` = Alpha,
         `b-ZEA` = Beta,
         AOH =  ALT) |> 
  # change NA's to 0 (no toxin detected)
  mutate(across(everything(), ~ replace_na(., 0))) |>
  # remove duplicated sample
  filter(Food_tested != "Organic_fruit_and_veggie_bars") |> 
  # remove extraneous characters from sample ID
  mutate(Number = str_remove(Number, " \\*2")) |> 
  # tidy to long format
  pivot_longer(cols = 3:36, names_to = "toxin_abb", values_to = "amount") |> 
  # join with table of food categories
  rename(ID = Number) |> 
  inner_join(products, by = c("ID", "Food_tested")) |> 
  # JOIN WITH TOXIN CATEGORY WHEN AVAILABLE
  inner_join(toxins, by = "toxin_abb") |> 
  select(ID, Food_tested, toxin, amount, Food_type, toxin_abb, toxin_type,
         LOQ) |> 
  mutate(across(c(Food_type, toxin_type), factor)) |> 
  mutate(Food_type = fct_other(Food_type, drop = "Miscellaneous"))

# save data in .Rdata and .xlsx formats -----------------------------------
save(main_data, products, toxins, 
     file = here("data", "toxin_data.Rdata"))
writexl::write_xlsx(
  main_data, 
  here("data", 
       paste0("working_data-", as.character(today()), ".xlsx")
       ))

# OK TO HERE --------------------------------------------------------------


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
  