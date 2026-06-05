# Load packages
library(tidyverse)
library(cowplot)
library(sf) # for handling ecoregions

# Import dataset
# From Arik's "signatures" folder on Google Drive
df <- read.csv("../data/julia_signatures_wPettitt-trends.csv")

# Read the column names
colnames(df)

# Make a correlation matrix of the median values
cor_med_plot <-
  df %>% 
  select(contains("_median")) %>%
  drop_na() %>%
  cor() %>% 
  # Convert to long format, keeping only unique pairs (upper triangle)
  as.data.frame() %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "correlation") %>%
  filter(var1 < var2) %>%  # keeps only unique pairs, drops diagonal & duplicates
  mutate(var1 = str_remove(var1, "_median"),
         var2 = str_remove(var2, "_median")) %>% 
  # Plot
  ggplot(aes(x = var1, y = var2, fill = correlation)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = RColorBrewer::brewer.pal(11, "PuOr"),
    limits = c(-1, 1)
  ) +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0))  +
  labs(x = NULL, y = NULL, fill = "Correlation", 
       title = "Pairwise Correlations of Median Values")
cor_med_plot
save_plot(plot = cor_med_plot,
          filename = "figures/signatures_median_corr_matrix.png",
          base_width = 18, base_height = 14)
  
# Make a correlation matrix of the mean values
cor_mean_plot <-
  df %>% 
  select(contains("_mean")) %>%
  select(!contains("pettitt")) %>% 
  drop_na() %>%
  cor() %>% 
  # Convert to long format, keeping only unique pairs (upper triangle)
  as.data.frame() %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "correlation") %>%
  filter(var1 < var2) %>%  # keeps only unique pairs, drops diagonal & duplicates
  mutate(var1 = str_remove(var1, "_mean"),
         var2 = str_remove(var2, "_mean")) %>% 
  # Plot
  ggplot(aes(x = var1, y = var2, fill = correlation)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = RColorBrewer::brewer.pal(11, "PuOr"),
    limits = c(-1, 1)
  ) +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0))  +
  labs(x = NULL, y = NULL, fill = "Correlation", 
       title = "Pairwise Correlations of Mean Values")
cor_mean_plot
save_plot(plot = cor_mean_plot,
          filename = "figures/signatures_mean_corr_matrix.png",
          base_width = 18, base_height = 14)

# Make a correlation matrix of the trends
cor_senslope_plot <-
  df %>% 
  select(contains("_senn_slp")) %>%
  drop_na() %>%
  cor() %>% 
  # Convert to long format, keeping only unique pairs (upper triangle)
  as.data.frame() %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "correlation") %>%
  filter(var1 < var2) %>%  # keeps only unique pairs, drops diagonal & duplicates
  mutate(var1 = str_remove(var1, "_senn_slp"),
         var2 = str_remove(var2, "_senn_slp")) %>% 
  # Plot
  ggplot(aes(x = var1, y = var2, fill = correlation)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = RColorBrewer::brewer.pal(11, "PuOr"),
    limits = c(-1, 1)
  ) +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0)) +
  labs(x = NULL, y = NULL, fill = "Correlation", 
       title = "Pairwise Correlations of Sen's Slopes")
cor_senslope_plot
save_plot(plot = cor_senslope_plot,
          filename = "figures/signatures_senslope_corr_matrix.png",
          base_width = 18, base_height = 14)

# Make table to determine how frequent each metric trend is significant
sig_thresh = 0.05
sig_table <- df %>% 
  select(contains("_mk_pval")) %>% 
  pivot_longer(cols = everything(), names_to = "metric", values_to = "value") %>% 
  drop_na() %>%
  group_by(metric) %>% 
  summarize(n_total = n(),
            n_sig = sum(value <= sig_thresh),
            pct_sig = (n_sig/n_total) * 100)

# Make a table of sig/non-sig trend values
trend_table <- left_join(
  df %>% 
    select(gage_id, contains("_mk_pval")) %>% 
    pivot_longer(cols = -gage_id, names_to = "metric", values_to = "mk_pval") %>% 
    mutate(metric = str_remove(metric, "_mk_pval")),
  df %>% 
    select(gage_id, contains("_senn_slp")) %>% 
    pivot_longer(cols = -gage_id, names_to = "metric", values_to = "sen_slope") %>% 
    mutate(metric = str_remove(metric, "_senn_slp")),
  by = c("gage_id", "metric")
) %>% left_join(
  df %>% 
    select(gage_id, contains("_median")) %>% 
    pivot_longer(cols = -gage_id, names_to = "metric", values_to = "median") %>% 
    mutate(metric = str_remove(metric, "_median")),
  by = c("gage_id", "metric")
) %>% 
  mutate(sig = ifelse(mk_pval <= sig_thresh, "yes", "no")) %>% 
  group_by(metric, sig) %>% 
  drop_na() %>% 
  summarize(n_total = n(),
            median_value = median(median),
            median_trend = median(sen_slope),
            median_trend_pct_median_value = (median_trend/median_value) * 100
            )


# Import N. Am. ecoregion data 
ecoregions <- st_read("../data/na_cec_eco_l2/NA_CEC_Eco_Level2.shp")

# Convert df to sf object using lat/lon columns
df_sf <- df %>%
  select(gage_id:longitude) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Reproject ecoregions to WGS84 if needed (EPA shapefiles are often in Albers)
ecoregions <- st_transform(ecoregions, crs = 4326) %>% 
  st_make_valid()

# Spatial join to extract ecoregion for each point
df_with_ecoregion <- df_sf %>%
  st_join(ecoregions, join = st_within)

# Bind to main df for summarization
df <- df %>% 
  left_join(df_with_ecoregion %>% 
              select(gage_id, eco_level2 = NA_L2NAME, eco_level1 = NA_L1NAME),
            by = "gage_id")

# Summarize sig trends by ecoregion
sig_table_eco <- df %>% 
  select(contains("_mk_pval"), eco_level2) %>% 
  pivot_longer(cols = -eco_level2, names_to = "metric", values_to = "value") %>% 
  drop_na() %>%
  group_by(metric, eco_level2) %>% 
  summarize(n_total = n(),
            n_sig = sum(value <= sig_thresh),
            pct_sig = (n_sig/n_total) * 100)

# Plot significance by ecoregion
sig_by_eco_plot <- 
  sig_table_eco %>%
  filter(!grepl('pettitt', metric)) %>% 
  mutate(metric = str_remove(metric, "_mk_pval")) %>% 
  ggplot(aes(x = metric, y = eco_level2, fill = pct_sig)) +
  geom_tile() +
  scale_fill_viridis_c() +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0))  +
  labs(x = NULL, y = NULL, fill = "% Significant", 
       title = "Metric Significance by Level 2 Ecoregion")
sig_by_eco_plot
save_plot(plot = sig_by_eco_plot,
          filename = "figures/significance_by_ecoregion.png",
          base_width = 18, base_height = 14)


# Summarize significant trend magnitudes by ecoregion
trend_table_eco <- left_join(
  df %>% 
    select(gage_id, eco_level2, contains("_mk_pval")) %>% 
    pivot_longer(cols = -c(gage_id, eco_level2), names_to = "metric", values_to = "mk_pval") %>% 
    mutate(metric = str_remove(metric, "_mk_pval")),
  df %>% 
    select(gage_id, eco_level2, contains("_senn_slp")) %>% 
    pivot_longer(cols = -c(gage_id, eco_level2), names_to = "metric", values_to = "sen_slope") %>% 
    mutate(metric = str_remove(metric, "_senn_slp")),
  by = c("gage_id", "metric", "eco_level2")
) %>% 
  filter(mk_pval <= sig_thresh) %>% 
  drop_na() %>% 
  group_by(metric) %>% 
  mutate(trend_max = max(sen_slope), 
         trend_min = min(sen_slope),
         # trend_nrml = (sen_slope - trend_min) / 
         #   (trend_max - trend_min),
         trend_nrml = (sen_slope) / max(abs(trend_max), abs(trend_min))) %>% 
  ungroup() %>% 
  group_by(metric, eco_level2) %>% 
  summarize(trend_mean = mean(sen_slope),
            trend_median = median(sen_slope),
            trend_mean_nrml = mean(trend_nrml),
            trend_median_nrml = median(trend_nrml))
  
# Plot median normalized trend by ecoregion
trend_median_nrml_by_eco_plot <- 
  trend_table_eco %>%
  filter(!grepl('pettitt', metric)) %>% 
  ggplot(aes(x = metric, y = eco_level2, fill = trend_median_nrml)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = RColorBrewer::brewer.pal(11, "PuOr"),
    limits = c(-1, 1)
  ) +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0))  +
  labs(x = NULL, y = NULL, fill = "Normalized Sen Slope", 
       title = "Median Normalized Trend Magnitude by Level 2 Ecoregion")
trend_median_nrml_by_eco_plot
save_plot(plot = trend_median_nrml_by_eco_plot,
          filename = "figures/median_trend_normalized_by_ecoregion.png",
          base_width = 18, base_height = 14)

# Plot mean normalized trend by ecoregion
trend_mean_nrml_by_eco_plot <- 
  trend_table_eco %>%
  filter(!grepl('pettitt', metric)) %>% 
  ggplot(aes(x = metric, y = eco_level2, fill = trend_mean_nrml)) +
  geom_tile() +
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(
    colors = RColorBrewer::brewer.pal(11, "PuOr"),
    limits = c(-1, 1)
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0))  +
  labs(x = NULL, y = NULL, fill = "Normalized Sen Slope", 
       title = "Mean Normalized Trend Magnitude by Level 2 Ecoregion")
trend_mean_nrml_by_eco_plot
save_plot(plot = trend_mean_nrml_by_eco_plot,
          filename = "figures/mean_trend_normalized_by_ecoregion.png",
          base_width = 18, base_height = 14)


