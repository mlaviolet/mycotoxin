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
# "corn flour" file--all names alike; combine with specimen ID so unique

# QUESTIONS
# Coconut_Flour	S251511351-237971 *2 (why the "*2"?) FIXED
# in raw data spreadsheet cells E4 and M40 are apparently typos, confirm? YES
#   rows 25 and 29 have duplicate ID S250551125-223922 FIXED

# linking table of food and classification
products <- read_xlsx("data-raw/product-key.xlsx") |> 
  # cleanup
  mutate(ID = str_remove(ID, " \\*2")) |> 
  filter(Food_tested != "Organic_fruit_and_veggie_bars") |> 
  mutate(Food_type = fct_other(Food_type, drop = "Miscellaneous"))

# full name of toxin and type of toxin for reporting
toxins <- read_xlsx("data-raw/2026-03-04 data legend.xlsx") |> 
  mutate(Abbreviation = str_replace(Abbreviation, "15_Ace", "Ace_15"),
         Abbreviation = str_replace(Abbreviation, "3_Ace", "Ace_3")) |> 
  rename(toxin_abb = Abbreviation, toxin = Mycotoxin,
         toxin_type = `Mycotoxin Type`) |> 
  mutate(toxin_type = factor(toxin_type)) |> 
  rename(LOQ = `LOQ (ppm)`)

# import data from Excel file
# THIS IS DATA FOR LOQ--IMPORT DATA FROM LOD AND MERGE
# DATA SOURCE FOR LOD DATA IS
#   "Toddler study raw with LOD data saved 2026-06-15 with DS.xlsx"
main_data <- 
  read_xlsx(here("data-raw", 
                 # "data-raw/Toddler raw data 2026-06-05.xlsx"),
                 # "Toddler study raw with LOD data saved 2026-06-11.xlsx"),
                 "Toddler raw data with just LOQ 2026-06-15 with DS.xlsx"),
                       na = c("nd", "np", "dn", "< 0"),
                       # cells E4 and M40 are apparently typos
                       col_types = c(rep("text", 3), rep("numeric", 35)),
                       n_max = 118) |> 
  rename(ID = Number) |> 
  # delete empty column and count
  select(- `...3`, -Count) |> 
  # change names to syntactic; change ID field to "ID"
  rename(Ace_15 = `15_Ace`, 
         Ace_3 = `3_Ace`, 
         BEA = Beau) |> 
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
  ## remove duplicated sample
  # filter(Food_tested != "Organic_fruit_and_veggie_bars") |> 
  # remove extraneous characters from sample ID
  mutate(ID = str_remove(ID, " \\*2")) |> 
  # correct ID on "Rice Krispies"
  mutate(ID = str_replace(ID,
                          "S242960827-202086", 
                          "S242960827-202068")) |> 
  # tidy to long format
  pivot_longer(cols = 3:36, names_to = "toxin_abb", values_to = "amount") |> 
  # join with table of food categories
  inner_join(products, by = c("ID", "Food_tested")) |> 
  # join with table of toxin names and abbreviations
  inner_join(toxins, by = "toxin_abb") |> 
  select(ID, Food_tested, toxin, amount, Food_type, toxin_abb, toxin_type) |> 
  # change Food_type and toxin_type to factors
  mutate(across(c(Food_type, toxin_type), factor)) |> 
  # change "Miscellaneous" food type to "Other" so it appears last
  mutate(Food_type = fct_other(Food_type, drop = "Miscellaneous")) |> 
  arrange(ID)

detect_df <- 
  readxl::read_xlsx(
    here("data-raw", "Toddler study raw with LOD data saved 2026-06-11.xlsx"),
    na = "nd") |> 
  select(-`...3`, -Food_tested) |> 
  rename(ID = `...2`) |> 
  mutate(across(!ID, ~ if_else(is.na(.x), "No", "Detected"))) |> 
  pivot_longer(cols = !ID, names_to = "toxin_abb", values_to = "detected") |> 
  # getting 120 specimens instead of 118
  filter_out(is.na(ID)) |> 
  mutate(ID = str_remove(ID, " \\*2")) |> 
  # summarize(n_detect = sum(value == "Detected"), .by = toxin_abb) |> 
  # make abbreviations match
  mutate(toxin_abb = str_replace(toxin_abb, "Alpha", "a-ZEA"),
         toxin_abb = str_replace(toxin_abb, "Beta", "b-ZEA"),
         toxin_abb = str_replace(toxin_abb, "15_Ace", "Ace_15"),
         toxin_abb = str_replace(toxin_abb, "3_Ace", "Ace_3"),
         toxin_abb = str_replace(toxin_abb, "Beau", "BEA"),
         toxin_abb = str_replace(toxin_abb, "FX", "FUS-X"),
         toxin_abb = str_replace(toxin_abb, "GRIS", "GRI"),
         toxin_abb = str_replace(toxin_abb, "Sterig", "STC"),
         toxin_abb = str_replace(toxin_abb, "ALT", "AOH"),
         toxin_abb = str_replace(toxin_abb, "NEOS", "NEO"),
         toxin_abb = str_replace(toxin_abb, "ZONE", "ZEA")
  ) |> 
  mutate(detected = factor(detected))

main_data <- main_data |> 
  inner_join(detect_df, by = c("ID", "toxin_abb")) |> 
  arrange(ID, toxin_abb) |> 
  select(c(1,2,5,3,8,4,7,6))
rm(detect_df)

# save data in .Rdata and .xlsx formats -----------------------------------
# .Rdata format
# save(main_data, products, toxins,
#      file = here("data",
#                  paste0("toxin_data_", as.character(today()), ".Rdata")))
# # Excel format
# writexl::write_xlsx(
#   main_data,
#   path = here("data",
#               paste0("working_data_", as.character(today()), ".xlsx")
#               ))
# .rds format
# saveRDS(main_data,
#      file = here("data",
#                  paste0("toxin_data_", as.character(today()), ".rds")))

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
  