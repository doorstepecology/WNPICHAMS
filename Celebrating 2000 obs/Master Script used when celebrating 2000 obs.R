#Loading...----
##the packages needed ----
library(tidyverse)
library(vegan)
library(geosphere)
library(ggrepel)
## the latest version of the data ----
raw <- read_csv(file.choose())

# Cleaning ----
## extract the useful columns  ----
rawfilteredwithduplicates <- raw %>%
  select(
    id,
    observed_on,
    user_login,
    quality_grade,
    num_identification_agreements,
    latitude,
    longitude,
    positional_accuracy,
    private_latitude,
    private_longitude,
    public_positional_accuracy,
    scientific_name,
    url,
    `field:where is the plant growing?`,
    `field:modelled pre-1750 pct`,
    `field:optional detail: "native" to where?`,
    `field:optional detail: could the soil be mostly the original, local soil?`,
    `field:optional detail: types of "other" softscapes and hardscapes`,
    `field:tree canopy aerial imagery dates`,
    `field:types of wall surfaces - wnpichams optional field`,
    `field:what vegetation management affected this plant?`
  ) %>% 
  filter(str_count(scientific_name, "\\S+") >= 2)
## (but only when the obs are ID'd to species)

## remove likely repeats of same individual plant ----

exclusion_distance <- 10   # duplicate radius in metres — change this to tweak sensitivity

dedup_group <- function(df, radius) {
  df <- df %>% arrange(observed_on)
  n <- nrow(df)
  keep <- rep(TRUE, n)
  
  if (n > 1) {
    coords <- cbind(df$longitude, df$latitude)
    d <- distm(coords, fun = distHaversine)
    
    for (i in seq_len(n)) {
      if (keep[i]) {
        same_user_and_date <- (df$user_login == df$user_login[i]) & (df$observed_on == df$observed_on[i])
        later_within_radius <- (seq_len(n) > i) & (d[i, ] <= radius) & keep & !same_user_and_date
        keep[later_within_radius] <- FALSE
      }
    }
  }
  df$keep <- keep
  df
}

dedup_result <- rawfilteredwithduplicates %>%
  group_by(scientific_name, across(starts_with("field:"))) %>%
  group_modify(~ dedup_group(.x, radius = exclusion_distance)) %>%
  ungroup()

rawfiltered <- dedup_result %>%
  filter(keep) %>%
  select(-keep)

kept_and_removed <- dedup_result %>%
  group_by(scientific_name, across(starts_with("field:"))) %>%
  filter(any(!keep)) %>%
  ungroup() %>%
  mutate(status = if_else(keep, "kept", "removed")) %>%
  select(id, scientific_name, observed_on, url, status, starts_with("field:"))

#Initial analysis and processing ----
## Check to see how many obs, columns and species ----
dim(rawfiltered)
n_distinct(rawfiltered$scientific_name)
## on the version seen on 6/8, it's 1834 obs, 239 species, and 20 columns
## with deduplication at exclusion_distance = 4, it's 1818 obs

## See the grand species accumulation curve ----
## create the datapoints to plot
species_over_time <- rawfiltered %>%
  arrange(observed_on) %>%
  mutate(
    obs_number = row_number(),
    new_species = !duplicated(scientific_name),
    cumulative_species = cumsum(new_species)
  )
## create the plot itself
ggplot(species_over_time, aes(x = obs_number, y = cumulative_species)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  labs(x = "Observations (chronological)", y = "Species",
       title = "Species accumulation over time") +
  theme_minimal()

## assign habitat groupings ----
## this creates the three column sheet "obs", where the data is just an ID, a species, and a habitat
obs <- rawfiltered %>%
  mutate(
    habitat = case_when(
      `field:where is the plant growing?` == "a footpath" ~ "footpath",
      `field:where is the plant growing?` == "a garden bed" ~ "garden bed",
      `field:where is the plant growing?` == "a road" ~ "road",
      `field:where is the plant growing?` == "on another plant (epiphytes and aerial hemiparasites)" ~ "epiphytes and mistletoes",
      
      `field:where is the plant growing?` == "a lawn or any other regularly mown or slashed area" &
        `field:optional detail: could the soil be mostly the original, local soil?` == "Commonly cultivated running-grasses (e.g. Stenotaphrum) are present; the substrate seems to have been affected by turf-laying" ~ "lawn - turf",
      
      `field:where is the plant growing?` == "a lawn or any other regularly mown or slashed area" &
        `field:optional detail: could the soil be mostly the original, local soil?` == "maybe yes or maybe no, but I'm sure the original soil (if present) has been very impacted by construction works and/or earthworks" ~ "lawn - footprint",
      
      `field:where is the plant growing?` == "a lawn or any other regularly mown or slashed area" ~ "lawn - possible original soil",
      
      `field:where is the plant growing?` == "a wall" &
        `field:types of wall surfaces - wnpichams optional field` == "brick" ~ "wall - brick",
      
      `field:where is the plant growing?` == "a wall" &
        `field:types of wall surfaces - wnpichams optional field` %in% c(
          "dry stone walls", "earthworks cutting", "gabion",
          "natural exposed rock-faces in built-up areas", "natural-stone masonry"
        ) ~ "wall - natural stone",
      
      `field:where is the plant growing?` == "a wall" ~ "wall - other",
      
      `field:optional detail: types of "other" softscapes and hardscapes` == "drainage infrastructure" ~ "drainage infrastructure",
      
      `field:optional detail: types of "other" softscapes and hardscapes` == "softscapes with no clear vegetation management, which are in the immediate vicinity of a constructed hardscape" ~ "unmanged in hardscape",
      
      `field:where is the plant growing?` == "other constructed hardscape" ~ "other",
      `field:where is the plant growing?` == "other softscape" ~ "other",
      
      TRUE ~ NA_character_
    )
  ) %>%
  select(id, scientific_name, habitat)

## create the community table ----
## this creates a table where every row is an observation, and every column is a species. 
## (as an aside, this is where the logic of analysing the iNat data as though it is community 
## data fails a little; ideally each row is a replicate of a sample of a community,not a lone 
## point. In other words, all those zeroes are not as meaningful as they should be. Still,
## it's the best way to analyse things I can think of right now)
comm <- obs %>%
  count(habitat, scientific_name) %>%
  pivot_wider(names_from = scientific_name, values_from = n, values_fill = 0) %>%
  column_to_rownames("habitat")

## create a Hurlbert's rarefaction curve for final habitats ----
rare_data <- rarecurve(as.matrix(comm), step = 5, tidy = TRUE)

ggplot(rare_data, aes(x = Sample, y = Species, color = Site)) +
  geom_line(linewidth = 0.8) +
  labs(x = "Number of observations", y = "Species richness", color = "Habitat") +
  theme_minimal()

## so, still got a lot of sampling to do. Lets pretend otherwise for fun
## rarefy ----
##I'd be a little bit more comfy rarefying to n=100, when there is enough data
##rarefaction parameters ----
min_n <- 50
set.seed(1312)
n_reps <- 1000
rarefied_list <- replicate(n_reps, rrarefy(comm, min_n), simplify = FALSE)

comm_rarefied_avg <- Reduce(`+`, rarefied_list) / n_reps

rarefied_array <- simplify2array(rarefied_list)
comm_rarefied_sd <- apply(rarefied_array, c(1,2), sd)
## I'm not quite sure how to handle the variance associated with our rarefaction

## rarefied richness bar chart ----
richness_per_rep <- sapply(rarefied_list, function(mat) rowSums(mat > 0))

richness_summary <- data.frame(
  Habitat = rownames(comm_rarefied_avg),
  Mean_Richness = rowMeans(richness_per_rep),
  SD_Richness = apply(richness_per_rep, 1, sd)
)

ggplot(richness_summary, aes(x = reorder(Habitat, Mean_Richness), y = Mean_Richness)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = Mean_Richness - SD_Richness, ymax = Mean_Richness + SD_Richness), width = 0.2) +
  coord_flip() +
  labs(x = "Habitat", y = "Species richness (rarefied to n = 50, 1 s.d.)", title = "Rarefied species richness by habitat") +
  theme_minimal()

##Calculate ecological distances ----
bc_dist <- vegdist(comm_rarefied_avg, method = "bray")

hc <- hclust(bc_dist, method = "average")

##Make a dendrogram ----
library(ggdendro)
ggdendrogram(hc, rotate = TRUE) +
  labs(title = "Habitat dissimilarity (Bray-Curtis, UPGMA)")

comm_final <- comm_rarefied_avg[!rownames(comm_rarefied_avg) %in% c("epiphytes and mistletoes", "other", "unmanged in hardscape"), ]

bc_dist_final <- vegdist(comm_final, method = "bray")

##Make an ordination plot ----
set.seed(1213)
nmds_final <- metaMDS(bc_dist_final, k = 2, trymax = 100)

nmds_final_scores <- as.data.frame(scores(nmds_final))
nmds_final_scores$Habitat <- rownames(nmds_final_scores)

ggplot(nmds_final_scores, aes(x = NMDS1, y = NMDS2, label = Habitat)) +
  geom_point(size = 3, color = "steelblue") +
  geom_text_repel(size = 3) +
  theme_minimal() +
  labs(title = "nMDS of major habitat communities",
       subtitle = paste("Stress =", round(nmds_final$stress, 3)))