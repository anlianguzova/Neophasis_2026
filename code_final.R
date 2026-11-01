# PT = post testicular length
# OS = oral sucker
# VS = ventral sucker
# Ph = pharynx
# Oes = oesophagus
# IB-VS = intestinal bifurcation - ventral sucker
# Vit-VS = vitellarium - ventral sucker
# VS-Ov = ventral sucker - ovary
# Ov = ovary
# Tes (A) = testis (anterior)
# Tes (B) = testis (posterior)
# Tes overlap = testis overlap
# Lat tes overlap = lateral testis overlap
# C-PE - C (or caeca/ceca) refers to the blind-ending intestinal tubes of the fluke. PE stands for the posterior extremity (the very end of the worm's tail/body).

library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(vegan)
library(car)
library(multcomp)
library(lme4)
library(lmerTest)
library(emmeans)
library(sf)
library(rnaturalearth)

# Data ####

maritae <- read_xlsx('data/Dimensions_250926_final.xlsx', sheet = 1)
maritae$Haplogroup <- as.factor(maritae$Haplogroup)
maritae$Location <- as.factor(maritae$Location)
maritae$Region <-  as.factor(maritae$Region)
maritae[, -c(1:4)] <- lapply(maritae[, -c(1:4)], as.numeric)

eggs <- read_xlsx('data/Dimensions_250926_final.xlsx', sheet = 2)

loc_levels <- c("Barents Sea", "White Sea", "Pechora Sea", "East Kamchatka")
loc_colors <- c("East Kamchatka" = '#E51A4B',
                "White Sea" = '#5CA595',
                "Pechora Sea" = '#00A0E3',
                "Barents Sea" = '#A79ECD')

# PCA analysis ####

## Feature selection ####

selection_ref <- c(haplogroup = "Haplogroup",
                   id = "Parasite ID",
                   region = "Region",
                   location = "Location",
                   body_l = "Body, l",
                   body_w = "Body, w (mid VS)",
                   forebody_l = "Forebody, l",
                   osuc_l = "OS, l",
                   osuc_w = "OS, w",
                   vsuc_l = "VS, l",
                   vsuc_w = "VS, w",
                   ph_l = "Ph, l",
                   ph_w = "Ph, w (mid)")

mar_ref <- maritae |> dplyr::select(all_of(selection_ref)) |> as.data.frame()
rownames(mar_ref) <- mar_ref$id

## Transformation ####

center <- function(x){
  x - mean(x, na.rm = TRUE)
}

mar_trans <- as.data.frame(lapply(mar_ref[, c(5:13)], log1p))
mar_cent <- t(apply(mar_trans, 1, center))

## PCA ####

pca_mar <- rda(mar_cent, scale = FALSE)
summary(pca_mar)

pca_scores   <- as.data.frame(scores(pca_mar, display = "sites",   choices = 1:2, scaling = 1))
pca_loadings <- as.data.frame(scores(pca_mar, display = "species", choices = 1:2, scaling = 2))

var_exp <- eigenvals(pca_mar)[1:2] / sum(eigenvals(pca_mar)) * 100

pca_scores$haplogroup <- mar_ref$haplogroup
pca_scores$location   <- mar_ref$location
pca_scores$region     <- mar_ref$region
pca_scores$id         <- mar_ref$id

centroids <- pca_scores |>
  group_by(region) |>
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

## Plots ####

p_load <- ggplot(pca_loadings, aes(x = PC1, y = PC2)) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "indianred") +
  geom_text(aes(label = row.names(pca_loadings)), color = "indianred", size = 3, vjust = -0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  labs(x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
       y = paste0("PC2 (", round(var_exp[2], 1), "%)")) +
  coord_equal() +
  theme_bw()

p_scores <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = location, shape = haplogroup)) +
  stat_ellipse(aes(group = haplogroup), linetype = "dashed", linewidth = 0.6) +
  geom_point(size = 4, alpha = 0.8) +
  geom_point(data = centroids, aes(x = PC1, y = PC2),
             color = "black", shape = 4, size = 5, stroke = 1.2,
             inherit.aes = FALSE) +
  geom_text_repel(data = centroids, aes(x = PC1, y = PC2, label = region),
                  color = "black", size = 4, fontface = "bold",
                  nudge_y = 0.1,
                  point.padding = 0.1,
                  box.padding = 0.1,
                  segment.color = NA,
                  inherit.aes = FALSE) +
  geom_text_repel(aes(label = id), size = 2.8, color = "black",
                  max.overlaps = 20, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  scale_color_manual(values = loc_colors, breaks = loc_levels) +
  labs(x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
       y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
       color = "Location", shape = "Haplogroup") +
  coord_equal() +
  theme_bw()

p_pca <- p_scores + p_load + plot_annotation(tag_levels = "A")
p_pca

# PERMANOVA ####

dist_mar <- dist(mar_cent, method = "euclidean")

adonis2(dist_mar ~ location + haplogroup, data = mar_ref,
        permutations = 9999, by = "margin")

## Dispersion ####

bd <- betadisper(dist_mar, mar_ref$location)
permutest(bd, permutations = 9999)
bd$group.distances
table(mar_ref$location)
boxplot(bd)

# Egg size analysis ####

## Data ####

eggs_raw <- eggs |>
  rename(length = 4, width = 5) |>
  mutate(Location = factor(Location, levels = loc_levels))

eggs_loc <- eggs_raw |>
  group_by(Haplogroup, `Parasite ID`, Location) |>
  summarise(mean_length = mean(length, na.rm = TRUE),
            mean_width  = mean(width,  na.rm = TRUE),
            n_eggs      = n(),
            .groups = "drop")

## Mixed model ####

m_mix <- lmer(length ~ Location + (1 | `Parasite ID`), data = eggs_raw)

## Model diagnostics ####

mix_diag <- eggs_raw |>
  mutate(.fitted    = fitted(m_mix),
         .std.resid = resid(m_mix, scaled = TRUE))

ggplot(data = mix_diag, aes(x = .fitted, y = .std.resid)) +
  geom_jitter(width = 0.3, height = 0) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "grey50") +
  labs(x = "Fitted values", y = "Standardized residuals")

ggplot(data = mix_diag, aes(x = Location, y = .std.resid)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  labs(y = "Standardized residuals")

qqPlot(resid(m_mix, scaled = TRUE), id = FALSE)
qqPlot(ranef(m_mix)$`Parasite ID`[, 1], id = FALSE)

leveneTest(.std.resid ~ Location, data = mix_diag)

infl_worm <- influence(m_mix, groups = "Parasite ID")
cd_worm <- as.matrix(cooks.distance(infl_worm))

ggplot(data.frame(worm = rownames(cd_worm), cd = cd_worm[, 1]),
       aes(x = worm, y = cd)) +
  geom_col() +
  geom_hline(yintercept = 4 / nrow(cd_worm), linetype = "dashed", color = "red") +
  labs(x = "Parasite ID", y = "Cook's distance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

mix_diag[abs(mix_diag$.std.resid) > 3, ]

## Model results ####

anova(m_mix)
summary(m_mix)

m_mix_full <- lmer(length ~ Location + Haplogroup + (1 | `Parasite ID`), data = eggs_raw)
anova(m_mix_full)
confint(m_mix_full, parm = "HaplogroupHg2")

## Post-hoc ####

emm_loc <- emmeans(m_mix, ~ Location, lmer.df = "satterthwaite")
pairs(emm_loc, adjust = "tukey")

cld_mix <- cld(emm_loc, Letters = letters, adjust = "tukey")

## Plot ####

MyData <- as.data.frame(emm_loc) |>
  dplyr::select(Location, fit = emmean, lwr = lower.CL, upr = upper.CL) |>
  left_join(as.data.frame(cld_mix) |>
              dplyr::select(Location, .group) |>
              mutate(.group = trimws(.group)),
            by = "Location")

y_top <- max(eggs_loc$mean_length) + 3

pos_j <- position_jitter(width = 0.08, height = 0, seed = 42)

gg_means <- ggplot(eggs_loc, aes(x = Location, y = mean_length, color = Location)) +
  geom_point(aes(shape = Haplogroup),
             position = pos_j, size = 2.2, alpha = 0.7, stroke = 1) +
  geom_text_repel(aes(label = `Parasite ID`),
                  position = pos_j, size = 2.5, max.overlaps = Inf,
                  min.segment.length = 0, segment.size = 0.2,
                  seed = 42, show.legend = FALSE) +
  geom_pointrange(data = MyData, aes(y = fit, ymin = lwr, ymax = upr),
                  position = position_nudge(x = 0.25),
                  size = 0.6, linewidth = 1) +
  geom_text(data = MyData, aes(y = y_top, label = .group),
            color = "black", size = 5) +
  scale_color_manual(values = loc_colors, guide = "none") +
  scale_shape_manual(values = c("Hg1" = 16, "Hg2" = 2)) +
  labs(x = NULL, y = "mean egg length (µm)", shape = "Haplogroup") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        legend.position = "top")

gg_means
# Combined figure ####

p_combined <- p_scores +
  (p_load /
     (gg_means +
        scale_shape_manual(values = c("Hg1" = 16, "Hg2" = 17)) +
        theme(legend.position = "none"))) +
  plot_layout(widths = c(1.3, 1)) +
  plot_annotation(tag_levels = "A")

p_combined

ggsave("images/mar_pca_eggs.svg", p_combined, width = 14, height = 8, units = "in")

ggsave("images/mar_pca_eggs.eps", p_combined,
       device = cairo_ps, width = 14, height = 8,
       units = "in", fallback_resolution = 600)

ggsave("images/mar_pca_eggs.pdf", p_combined, width = 14, height = 8, units = "in")

# Mantel test ####

## Distance matrices ####

pops <- c("East_Siberian_Sea", "East_Kamchatka", "WS",
          "Pechora_Sea", "Barents_Sea", "Sea_of_Japan")

# builds a symmetric matrix from a row-wise lower triangle incl. diagonal
tri_to_sym <- function(vals, labels) {
  n <- length(labels)
  stopifnot(length(vals) == n * (n + 1) / 2)
  m <- matrix(0, n, n, dimnames = list(labels, labels))
  m[upper.tri(m, diag = TRUE)] <- vals
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}

fst_matrix <- tri_to_sym(c(
  0.00000,
  0.00140, 0.00000,
  0.49371, 0.49244, 0.00000,
  0.23264, 0.23341, 0.36251, 0.00000,
  0.46377, 0.46667, 0.70653, 0.12698, 0.00000,
  0.19048, 0.22533, 0.51267, 0.33054, 0.58915, 0.00000
), pops)

geo_matrix <- tri_to_sym(c(
  0,
  2116, 0,
  4390, 5860, 0,
  4081, 5572, 1180, 0,
  4413, 5640,  592,  868, 0,
  3086, 1889, 6631, 6307, 6459, 0
), pops)

fst_dist <- as.dist(fst_matrix)
geo_dist <- as.dist(geo_matrix)

## Test ####

mantel_result <- mantel(fst_dist, geo_dist, method = "pearson", permutations = 10000)
mantel_result

## Plot ####

plot_df <- data.frame(geographic_distance_km = as.vector(geo_dist),
                      genetic_distance_FST   = as.vector(fst_dist))

p_ibd <- ggplot(plot_df, aes(x = geographic_distance_km, y = genetic_distance_FST)) +
  geom_point(size = 3, color = "blue", alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink") +
  labs(title = "Isolation by Distance within Hg1+2 (N. oculata s.s.)",
       x = "Geographic distance (km)",
       y = expression(pairwise ~ F[ST])) +
  annotate("text",
           x = max(plot_df$geographic_distance_km) * 0.7,
           y = max(plot_df$genetic_distance_FST) * 0.9,
           label = paste("Mantel r =", round(mantel_result$statistic, 3),
                         "\np =", round(mantel_result$signif, 3)),
           hjust = 0, size = 4) +
  theme_bw()

p_ibd

# Map ####

## Settings ####
lon_min <- 0;  lon_max <- 205   # degrees E; 205 = 155 W, so Chukchi and Bering Seas are included
lat_min <- 30; lat_max <- 80    # degrees N
land_col <- "grey30"
sea_col  <- "white"

crs_map <- "+proj=lcc +lat_1=40 +lat_2=70 +lat_0=60 +lon_0=100 +datum=WGS84 +units=m +no_defs" # central meridian moved east to balance the wider extent

## Adding  land and lakes ####
land <- ne_countries(scale = "medium", returnclass = "sf")
land <- land[land$continent != "Antarctica", ]
lakes <- ne_download(scale = 50, type = "lakes", category = "physical", returnclass = "sf")  # needs internet once

## Setting map limits ####
lon_seq <- seq(lon_min, lon_max, length.out = 200)
lat_seq <- seq(lat_min, lat_max, length.out = 200)
ring <- rbind(cbind(lon_seq, lat_min), cbind(lon_max, lat_seq),
              cbind(rev(lon_seq), lat_max), cbind(lon_min, rev(lat_seq)))
ring <- rbind(ring, ring[1, ])
frame <- st_transform(st_sfc(st_polygon(list(ring)), crs = 4326), crs_map)
lims  <- st_bbox(frame)

## Graticule (longitudes above 180 are allowed; PROJ wraps them) ####
parallels <- lapply(seq(30, 80, 10), function(la)
  st_linestring(cbind(seq(-60, 260, length.out = 600), la)))
meridians <- lapply(seq(-30, 210, 30), function(lo)
  st_linestring(cbind(lo, seq(20, 89, length.out = 200))))
grat <- st_sfc(c(parallels, meridians), crs = 4326)

## Plot ####
p <- ggplot() +
  geom_sf(data = land,  fill = land_col, colour = "grey60", linewidth = 0.2) +
  geom_sf(data = lakes, fill = sea_col,  colour = "grey60", linewidth = 0.1) +
  geom_sf(data = grat,  colour = "grey70", linewidth = 0.25) +
  coord_sf(crs = crs_map,
           xlim = unname(lims[c("xmin", "xmax")]),
           ylim = unname(lims[c("ymin", "ymax")]),
           expand = FALSE) +
  theme_void() +
  theme(panel.background = element_rect(fill = sea_col, colour = NA),
        panel.border     = element_rect(fill = NA, colour = "black", linewidth = 0.6))
p

ggsave("images/overview_map.svg", p, width = 11, height = 7, device = svglite::svglite)
ggsave("images/overview_map.pdf", p, width = 11, height = 7, device = cairo_pdf)

