# generate pairs of toxins
library(here)
library(tidyverse)

load(here("data", "toxin_data_2026-06-06.Rdata"))

# data with collapsed categories
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

# find number of toxin groups by specimen ID
df_groups <- work_data |> 
  select(ID, toxin_grp) |> 
  distinct() |> 
  count(ID) |> 
  arrange(ID)

count(df_groups, n)

# build subset data set
#   choose one ID from 2, 3, and 4 groups
set.seed(42) 
df_ids <- df_groups |> 
# remove IDs with only one toxin group
  filter_out(n == 1) |> 
  # choose one each from 2, 3 , 4, 5, 6 groups
  slice_sample(n = 1, by = n) |> 
  arrange(n)

df_sub <- work_data |> 
  select(ID, toxin_grp) |> 
  filter(ID %in% df_ids$ID)

# looks good
# work_data |> 
#   filter(ID == "S250221430-217975")

# # wrap this into a function
# unique_vals <- unique(as.character(work_data$toxin_grp))
# t(combn(unique_vals, 2))

get_combn <- function(df) {
  unique_vals <- unique(as.character(df$toxin_grp))
  t(combn(unique_vals, 2))
  }

# test function
work_data |> 
    select(ID, toxin_grp) |> 
    filter(ID == "S242960827-202082") |> 
    get_combn()


# search "purrr map apply function by group"
# build dataframe with list-colmun
nested_df <- df_sub %>% 
    group_by(ID) %>% 
    nest() |> 
  mutate(matrix_col = map(data, ~ get_combn(.x))) 
  # mutate(x = map(data, ~ get_combn(.x))) 

# dataframe with pairs of toxin groups by ID
result_df <- nested_df |> 
  mutate(matrix_col = map(matrix_col, as.data.frame)) %>%
  # unnest_wider(matrix_col)
  # USE THIS INSTEAD
  unnest(matrix_col) |> 
  select(-data) |> 
  ungroup()

# count occurrences of each pair
counts <- result_df |> 
  select(-ID) |> 
  count(V1, V2) |> 
  arrange(desc(n))

# all unique?
counts |> 
  select(-n) |> 
  distinct() |> 
  nrow()

# "Polyketide" and "Other" appear in rows 13 and 15, in reverse order
# Remove "Other" from data? -- possibly not meaningful anyway
# TRY WITH ENTIRE DATA SET

# OK TO HERE --------------------------------------------------------------
summary(df_sub$toxin_grp)
levels(df_sub$toxin_grp) |> combn(2)


# THIS LOOKS LIKE IT--CONFIRM
# TRY WITH WHOLE DATA AND RESHAPE AS MATRIX

result_df <- nested_df %>% 
  mutate(x = map(data, ~ get_combn(.x))) |> 
  unnest(x) |> 
  select(-data) 
  # set_names(c("ID", "grp1", "grp2"))

# show matrix from list column
result_df$x[[1]]


# RENAME COLUMNS
# looks good so far -- next, tally pairs 

# r extract all columns of matrix from list-column
# Extract all columns dynamically

  



# EXAMPLE: extract all columns from matrix --------------------------------
# Create a sample data frame with a matrix in a list-column
df <- tibble(
  id = 1:2,
  matrix_col = list(
    matrix(1:4, nrow = 2, dimnames = list(NULL, c("A", "B"))),
    matrix(5:8, nrow = 2, dimnames = list(NULL, c("A", "B")))
    )
  )

# Extract all columns dynamically
df_expanded <- df %>%
  mutate(matrix_col = map(matrix_col, as.data.frame)) %>%
  # unnest_wider(matrix_col)
  # USE THIS INSTEAD
  unnest(matrix_col)


  
df_test$ID

# example of each number of toxins
# 0: S242960827-202068
# 1: S242781609-199274 
# 2: S240241332-165026
# 3: S242960827-202077
# 4: S242960827-202090 
# 5: S242960827-202080
# 6: S242781609-199271

# test data 
df <- work_data |> 
  select(ID, toxin_grp) |> 
  filter(ID %in% c("S242960827-202068",
                   "S242781609-199274",
                   "S240241332-165026",
                   "S242960827-202077",
                   "S242960827-202090",
                   "S242960827-202080",
                   "S242781609-199271"))
count(df, ID) |> 
  arrange(n)



df_test <- work_data |> 
  filter(ID == "S242781609-199271")
get_combn(df_test)  



test1 <- map()
  group_by(ID)


# use following for testing; has six toxins
# S242961615-201994


unique_vals <- unique(as.character(df$toxin_grp))

# use head(number_toxins, 10) for testing


# is this what I want? DOESN'T SEEM SO
df <- pairwise_count(work_data, toxin_grp, ID, sort = TRUE, diag = FALSE,
                     upper = FALSE) |> 
  print(n = Inf)
sum(df$n)

# TRY USING LEVELS AS ARGUMENT TO as.character
# get_combn <- function(df) {
#   unique_vals <- as.character(levels(df$toxin_grp))
#   t(combn(unique_vals, 2))
#   }


# result_df |> 
#   mutate(x = map(x, as.data.frame)) |> 
#   unnest(x) |> 
#   select(-data)

# df_sub |> 
#   distinct()
#   split(df_groups$ID) |