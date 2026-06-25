# construct heatmap of counts by food type and toxin class

# reordered factor levels; put flours adjacent
new_order <- c("Wheat flour", "Corn flour", "Non-grain", "Breakfast cereals",   
               "Pasta", "Snack foods", "Other")

# use tidied data
dfx <- main_data |> 
  # want to count detecteds, remove others
  filter_out(detected == "No") |> 
  # collapse similar toxin groups to get manageable number of groups
  mutate(
    toxin_grp = 
      fct_collapse(
        toxin_type,
        Trichothecene = c("Trichothecene Type A", 
                          "Trichothecene Type B"),
        Difuranocoumarin = c("Difuranocoumarin",
                             "Difuranocoumarin xanthone precursor to aflatoxin")
      )) |> 
  # keep top 5 toxin groups and collapse others
  mutate(toxin_grp = fct_lump_n(toxin_grp, 5)) |> 
  # reorder levels to have flours together
  # change "Other" food category to "First foods"
  mutate(Food_type = fct_relevel(Food_type, new_order),
         Food_type = fct_recode(Food_type, "First foods" = "Other")) |> 
  # include zeros in counts
  count(Food_type, toxin_grp, .drop = FALSE) |> 
  # compute normalized toxin counts by food groups, shows up better on map
  group_by(Food_type) |> 
  mutate(z = scale(n)) |> 
  ungroup()
  
# make the heatmap with food types on horizontal and toxin classes on vertical
dfx |> 
  ggplot(aes(x = Food_type, y = toxin_grp, fill = z)) +
  geom_tile(color = "grey") +
  scale_fill_distiller(palette = "RdYlBu") +
  scale_y_discrete(limits = rev) +
  labs(x = "Food type", y = "Toxin class", 
       caption = "Numbers are counts of detections; colors are based on normalized counts") +
  geom_text(aes(label = n), color = "black") +
  theme(legend.position = "none")
  
  