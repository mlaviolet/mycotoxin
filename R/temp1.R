new_order <- c("Wheat flour", "Corn flour", "Non-grain", "Breakfast cereals",   
               "Pasta", "Snack foods", "Other")

hmap_data <- work_data |> 
  mutate(Food_type = fct_relevel(Food_type, new_order),
         Food_type = fct_recode(Food_type, "First foods" = "Other"  )) |> 
  count(Food_type, toxin_grp, .drop = FALSE) |> 
  group_by(Food_type) |> 
  mutate(norm_count = scale(n)) |> 
  ungroup() 
  
hmap_data |>   
  ggplot(aes(x = Food_type, y = toxin_grp, fill = norm_count)) +
  geom_tile(color = "grey") +   # borders between tiles
  # scale_fill_viridis_c(option = "magma") +
  scale_fill_distiller(palette = "RdYlBu") +
  scale_y_discrete(limits = rev) +
  # labs(fill = "Standardized\ncount", x = "Food type", y = "Toxin class") +
  labs(x = "Food type", y = "Toxin class", 
       caption = "Colors are based on normalized counts") +
  geom_text(aes(label = n), color = "black") +
  theme(legend.position = "none")
  # theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
  #       legend.position = "none")
  # geom_tile(
  #   data = subset(work_data, norm_count == 0),
  #   fill = NA,
  #   color = "black",
  #   linewidth = 1) 

by_raw_count <- hmap_data |> 
  select(-norm_count) |> 
  pivot_wider(names_from = Food_type, values_from = n) |> 
  writexl::write_xlsx(
    here("output",
         paste0("by_raw_count_", as.character(today()), ".xlsx")
    ))

by_std_count <- hmap_data |> 
  select(-n) |> 
  pivot_wider(names_from = Food_type, values_from = norm_count) |> 
  writexl::write_xlsx(
    here("output",
         paste0("by_std_count_", as.character(today()), ".xlsx")
    ))

  


# reconstructed Table 3

tbl3 <- 
  readxl::read_xlsx(
    here("data-raw", "Toddler study raw with LOD data saved 2026-06-11.xlsx"),
    na = "nd") |> 
  select(-`...3`, -Food_tested) |> 
  rename(ID = `...2`) |> 
  mutate(across(!ID, ~ if_else(is.na(.x), "No", "Detected"))) |> 
  pivot_longer(cols = !ID, names_to = "toxin_abb") |> 
  # getting 120 specimens instead of 118
  filter_out(is.na(ID)) |> 
  summarize(n_detect = sum(value == "Detected"), .by = toxin_abb) |> 
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
         )


df <- table_4 |> 
  full_join(tbl3, by = "toxin_abb")

# x <- sort(tbl3$toxin_abb)
# y <- sort(table_4$toxin_abb)
# x == y
# setdiff(y, x)
  