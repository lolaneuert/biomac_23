### ITALIAN FOREST DATA FILTERING
library(sf)
library(rgdal)
library(tidyverse)
library(rgeos)

# UPLOAD DATA
path_forest <- "/home/nicola/OneDrive/forest_ita/sampling_for/R_sampling_for/data"

# forest <- readRDS(file.path(path_forest, paste("preferential_accurate_data", ".rds", sep=""))) # list of data.frame(s)

forest_veg <- readRDS(file.path(path_forest, paste("pref_3gr", ".rds", sep="")))

forest_plot <- readRDS(file.path(path_forest, paste("preferential_accurate_data", ".rds", sep=""))) %>% 
        purrr::pluck("plots") %>% 
        pivot_longer(!"sp", names_to = "id", values_to = "cover") %>% 
        pivot_wider(names_from = "sp", values_from = "cover") %>% 
        mutate(id = as.numeric(id))
forest_header <- readRDS(file.path(path_forest, paste("preferential_accurate_data", ".rds", sep=""))) %>% 
        pluck("header") 

forest_env <- readRDS(file.path(path_forest, paste("preferential_accurate_data", ".rds", sep=""))) %>% 
  pluck("env") 

forest_point <- forest_veg %>% ### Preferential sites
        mutate(id = as.numeric(id)) %>% 
        left_join(forest_header,
                  dplyr::select(c("id", "lon", "lat")),
                  by = "id") %>% 
        SpatialPointsDataFrame(coords = .[, c("lon", "lat")],
                               data = .,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

grid_16km <- readOGR(dsn = "/home/nicola/OneDrive/gis/shapes/", 
                  layer = "grid16km") %>% 
        spTransform(., CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

forest_site <- raster::extract(grid_16km, forest_point, df = TRUE) %>% 
        cbind(forest_point$id) %>% 
        dplyr::select(c(8, 2)) %>%
        rename(id = 1, site = 2)

### ITALIAN FOREST: Preferential
forest_ita <- forest_veg %>% 
        rename(veg = nc_clus) %>% 
        mutate(id = as.numeric(id)) %>% 
        mutate(veg = replace(veg, veg == "N", "azon")) %>% 
        mutate(veg = replace(veg, veg == "s.5", "warm")) %>% 
        mutate(veg = replace(veg, veg == "s.3", "cool")) %>% 
        mutate(veg = replace(veg, veg == "s.4", "cold")) %>% 
        mutate(veg = factor(veg, levels = c("azon", "warm", "cool", "cold"))) %>% 
        left_join(forest_header %>%
                          dplyr::select(c("id", "lon", "lat", "size", "metric_inaccuracy", "year")),
                  by = "id") %>% 
        left_join(forest_site) %>% 
        relocate(veg, .after = size) %>% 
        relocate(site, .before = size) %>% 
        left_join(forest_plot, by = "id")

rm(forest_header, forest_plot, forest_point, forest_site, forest_veg, grid_16km)
# forest_ita %>% # relative size of vegetation types
#         pull(veg) %>% 
#         table() %>% 
#         as.data.frame() %>% 
#         rename(veg = ".",
#                n = Freq) %>% 
#         mutate(perc = round((n / sum(n)) * 100, 0)) 


### FILTERING
forest_sample <- forest_ita %>% 
        filter(size >= 100 | size <= 300) %>% # plot size of 100 m²
        group_by(veg) %>% 
        sample_frac(0.01) %>%
        ungroup

forest_sample_plot <- forest_sample %>% 
        dplyr::select(9:ncol(forest_sample)) %>% 
        dplyr::select(which(!colSums(., na.rm=TRUE) %in% 0))
        
forest_sample <- forest_sample %>% 
        dplyr::select(1:8) %>% 
        cbind(forest_sample_plot)
str(forest_sample_env)

forest_sample_env <- forest_env %>% 
  dplyr::filter(id %in% forest_sample[, "id"]) %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(data.frame(id = forest_sample$id), # order it according plots
            forest_sample_env, 
            by = "id") %>% 
  dplyr::select("L", "N", "altitude", "AMT", "MeTCoQ", "AP", "PS") %>% 
  mutate(AMT = round(AMT / 10, 1),
         MeTCoQ = round(MeTCoQ / 10, 1))

forest_sample <- list(plot = forest_sample,
                      env = forest_sample_env)
 

# forest_sample %>% # relative size of vegetation types
#         pull(veg) %>% 
#         table() %>% 
#         as.data.frame() %>% 
#         rename(veg = ".",
#                n = Freq) %>% 
#         mutate(perc = round((n / sum(n)) * 100, 1))
# 
# forest_sample %>% 
#         filter(veg == "azon") %>% 
#         pull(site) %>% 
#         unique() %>% 
#         length()


## ICP data set
icp <- readRDS(file.path(path_forest, paste("probabilistic_filtered", ".rds", sep="")))

icp_sample_warm <- icp %>% 
  filter(veg == "warm") %>% 
  group_by(site) %>%
  sample_n(2) %>%
  nest() %>% 
  ungroup() %>% 
  sample_n(6) %>% 
  unnest()

icp_sample_cool <- icp %>% 
  filter(veg == "cool") %>% 
  group_by(site) %>%
  sample_n(2) %>%
  nest() %>% 
  ungroup() %>% 
  sample_n(3) %>% 
  unnest()

icp_sample_cold <- icp %>% 
  filter(veg == "cold") %>% 
  group_by(site) %>%
  sample_n(2) %>%
  nest() %>% 
  ungroup() %>% 
  sample_n(2) %>% 
  unnest()

icp_sample <- icp_sample_warm %>% 
  rbind(icp_sample_cool) %>% 
  rbind(icp_sample_cold)

tree_sp <- readxl::read_xlsx(file.path(path_forest, paste("sp_traits", ".xlsx", sep="")),
                             na = "NA", sheet = 1) %>%
  filter(Form_1 == "Tree" |
           Form_1 == "Shrub") %>%
  dplyr::select(sp = 2) %>% 
  as.data.frame()

icp_head <- icp_sample %>% 
  dplyr::select(c("id", "lon", "lat", "site", "size", "veg"))
icp_tree_sp <- icp_sample %>% 
  dplyr::select(7:1105) %>% 
  rownames_to_column("id") %>% 
  pivot_longer(!id, names_to = "sp", values_to = "cover") %>% 
  left_join(tree_sp %>% 
              mutate(check = 1), 
            by = "sp") %>% 
  dplyr::filter(check == 1) %>% 
  dplyr::select(!check) %>% 
  pivot_wider(names_from = "sp",
              values_from = "cover") %>% 
  column_to_rownames("id")

sp_names = data.frame(sp = colnames(icp_tree_sp),
                      code = paste("sp.", seq(1:ncol(icp_tree_sp)), sep = "")) 
colnames(icp_tree_sp) <- sp_names$code


probabilistic_sample <- icp_head %>% 
  cbind(icp_tree_sp) %>% 
  mutate(metric_inaccuracy = 3,
         year = 2007) %>% 
  relocate(metric_inaccuracy, .after = veg) %>% 
  relocate(year, .after = metric_inaccuracy)
icp_env <- readRDS(file.path(path_forest, paste("probabilistic_plot_data", ".rds", sep=""))) 
icp_env <- icp_env$env

probabilistic_sample_env <- icp_env %>% 
  dplyr::filter(id %in% probabilistic_sample$id) %>% 
  mutate(id = as.numeric(id)) 

probabilistic_sample_env <- probabilistic_sample_env %>% 
  left_join(data.frame(id = probabilistic_sample$id), # order it according plots
            probabilistic_sample_env, 
            by = "id") %>% 
  dplyr::select("L", "N", "altitude", "AMT", "MeTCoQ", "AP", "PS") %>% 
  mutate(AMT = round(AMT / 10, 1),
         MeTCoQ = round(MeTCoQ / 10, 1))
probabilistic_sample <- list(plot = probabilistic_sample,
                             env = probabilistic_sample_env)

saveRDS(probabilistic_sample, "data/probabilistic_sample.rds")
saveRDS(sp_names, "data/sp_names_probabilistic_sample.rds")
saveRDS(forest_sample, "data/forest_sample.rds")
