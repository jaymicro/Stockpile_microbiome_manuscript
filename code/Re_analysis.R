
# Import packages ---------------------------------------------------------


library(qiime2R)
library(tidyverse)
library(GUniFrac)
library(vegan)
library(picante)
library(ape)
library(janitor)
library(ggpubr)
library(ranger)
library(paletteer)
library(pairwiseAdonis)
library(ggtext)




# Bacterial analyses ------------------------------------------------------



## Data import and wrangling -----------------------------------------------





bact_count <- read_tsv("../data/bact/bact_count.tsv") %>% 
  column_to_rownames(., var = "hash_id") %>% 
  t(.) %>% 
  as.data.frame(.)%>% 
  filter(rowSums(.)>1000) %>% 
  select_if(colSums(.)>100) %>% 
  adorn_percentages(,, `e8a524e024cfebb795310da5dd143dfe`:`bbd53bf9c1c05fe1cdfae19c7e947ec0`)



samp_name_bact <- rownames(bact_count)

bact_metadata <- readxl::read_xlsx("../data/bact/bacteria_meta.xlsx") %>% 
  filter(Sample_id %in% samp_name_bact) %>% 
  mutate(depth_group = case_when( sampling_depth == "reference" ~ "reference",
                                  site == "QR_mill" & sampling_depth %in% c("10","20", "60")  ~ "0 - 60",
                                  #  site == "QR_mill" & sampling_depth ==  "20" ~ "11-20",
                                  #  site == "QR_mill" & sampling_depth %in% c( "75", "90","105", "108") ~ "75-120",
                                  site == "QR_mill" & sampling_depth %in% c("75", "90","105", "108", "200", "210","250") ~ "75 - 260",
                                  site == "QR_mill" & sampling_depth %in% c("350", "365", "380") ~ "350 - 390",
                                  site == "QR_mill" & sampling_depth %in% c("500", "550", "565") ~ "500 - 575",
                                  site == "NAF" & sampling_depth %in% c("30", "60")  ~ "0 - 60",
                                  site == "NAF" & sampling_depth %in% c("90", "120", "150")  ~ "61 - 150",
                                  site == "NAF" & sampling_depth %in% c("152", "300","460","610"  )  ~ "152 - 610",
                                  site == "NAF" & sampling_depth %in% c("760", "910", "1070", "1220","1370","610"  )  ~ "611 - 1372")) %>% 
 # drop_na() %>% 
  column_to_rownames(., var = "Sample_id") %>% 
  select(-sample_name)


bact <- merge(bact_count, bact_metadata,
              by = 'row.names', all = F) %>% 
  column_to_rownames(., var = "Row.names")

bact_asv <- data.frame(hash_id = colnames(bact)) %>% 
  mutate(asv_id = case_when(hash_id == "site" ~ "site",
                            hash_id == "sampling_depth" ~ "sampling_depth",
                            hash_id == "depth_group" ~ "depth_group",
                            TRUE  ~  paste("ASV", row_number(), sep = "_"))
  )

bact_rf <- bact
colnames(bact_rf) <- bact_asv$asv_id


tree_tips_bact <- bact %>% 
  select(!c(site, sampling_depth, depth_group)) %>% 
  colnames()

tree_bact <- read_qza("../data/bact/ashley-bacterial.rooted-tree.qza")[["data"]] %>% 
  keep.tip(., tip = tree_tips_bact )


bact_tax <- read_qza("../data/bact/ashley-bacterial.taxonomy.qza")


tax_bact <- bact_tax[["data"]] %>% 
  as.data.frame() %>% 
  filter(Feature.ID %in% tree_tips_bact) %>% 
  # column_to_rownames(., var = "Feature.ID") %>% 
  select(!Confidence) %>% 
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") 


bact_naf_abund <- bact %>%
  filter(site == "NAF") %>% 
  select(!c( site, sampling_depth, depth_group)) %>% 
  select_if(colSums(.)>0) 


bact_qr_abund <- bact %>%
  filter(site == "QR_mill") %>% 
  select(!c( site, sampling_depth, depth_group))%>% 
  select_if(colSums(.)>0) 

asv_naf_bact <- colnames(bact_naf_abund)
asv_qr_bact <- colnames(bact_qr_abund)

tree_naf_bact <- tree_bact %>% 
  keep.tip(., tip = asv_naf_bact)

tree_qr_bact <- tree_bact %>% 
  keep.tip(., tip = asv_qr_bact)


## Alpha diversity analysis ------------------------------------------------

#faith's PD analysis

faiths_pd_naf_bact <- pd(bact_naf_abund, tree = tree_naf_bact) %>% 
  merge(., bact_metadata, by = 'row.names', all = F) %>% 
  #mutate(sampling_depth = as.numeric(sampling_depth)) %>% 
  filter(depth_group != "reference") %>% 
  mutate(depth_group = factor(depth_group, levels = c("reference", "0 - 60", "61 - 150", "152 - 610", "611 - 1372")))

naf_pd_sr_ref_bact <- pd(bact_naf_abund, tree = tree_naf_bact) %>% 
  merge(., bact_metadata, by = 'row.names', all = F) %>% 
  filter(depth_group == "reference") %>% 
  group_by(site) %>% 
  summarise(across(c(PD,SR), mean))

mod_pd_naf <- lm(PD ~ depth_group, faiths_pd_naf_bact)
summary.aov(mod_pd_naf)

mod_SR_naf <- lm(SR ~ depth_group, faiths_pd_naf_bact)
summary.aov(mod_SR_naf)


naf_bact_pd <- ggplot(faiths_pd_naf_bact, aes(x = depth_group, y = PD, fill = depth_group)) +
  geom_hline(yintercept = naf_pd_sr_ref_bact$PD, color = "steelblue" , linetype = "dotdash") +
  geom_violin(alpha = 0.3, show.legend = F) +
  geom_boxplot(width=0.2)+
  geom_pwc(ref.group = "0 - 60",
           method = "t.test",
           p.adjust.method = "BH",
           label = "p.adj.signif") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_fill_viridis_d(begin = 0.1, end = 0.9) +
  labs(x = "Stockpile depth (cm)",
       y = "Faith's phylogenetic diveristy",
       title = "16S rRNA gene")
naf_bact_pd

naf_bact_SR <- ggplot(faiths_pd_naf_bact, aes(x = depth_group, y = SR, fill = depth_group)) +
  geom_hline(yintercept = naf_pd_sr_ref_bact$SR, color = "steelblue" , linetype = "dotdash") +
  geom_violin(alpha = 0.3, show.legend = F) +
  geom_boxplot(width=0.2)+
  geom_pwc(ref.group = "0 - 60",
           method = "t.test",
           p.adjust.method = "BH",
           label = "p.adj.signif") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_fill_viridis_d(begin = 0.1, end = 0.9) +
  labs(x = "Stockpile depth (cm)",
       y = "Observed species richness",
       title = "16S rRNA gene")
naf_bact_SR


bact_faiths_pd_qr <- pd(bact_qr_abund, tree = tree_qr_bact)%>% 
  merge(., bact_metadata, by = 'row.names', all = F)  %>% 
  filter(depth_group != "reference") %>% 
  mutate(depth_group = factor(depth_group, levels = c("reference", "0 - 60","75 - 260","350 - 390","500 - 575" )))

bact_qr_pd_sr_ref <- pd(bact_qr_abund, tree = tree_qr_bact) %>% 
  merge(., bact_metadata, by = 'row.names', all = F) %>% 
  filter(depth_group == "reference") %>% 
  group_by(site) %>% 
  summarise(across(c(PD,SR), mean))

mod_pd_qr_bact <- lm(PD ~ depth_group, bact_faiths_pd_qr)
summary.aov(mod_pd_qr_bact)

mod_SR_qr_bact <- lm(SR ~ depth_group, bact_faiths_pd_qr)
summary.aov(mod_pd_qr_bact)


qr_bact_pd <- ggplot(bact_faiths_pd_qr, aes(x = depth_group, y = PD, fill = depth_group)) +
  geom_hline(yintercept = bact_qr_pd_sr_ref$PD, color = "steelblue" , linetype = "dotdash") +
  geom_violin(alpha = 0.3, show.legend = F) +
  geom_boxplot(width=0.1)+
  geom_pwc(ref.group = "0 - 60",
           method = "t.test",
           p.adjust.method = "BH",
           label = "p.adj.signif") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_fill_viridis_d(begin = 0.1, end = 0.9) +
  labs(x = "Stockpile depth (cm)",
       y = "Faith's phylogenetic diveristy",
       title = "16S rRNA gene")
qr_bact_pd

qr_bact_SR <- ggplot(bact_faiths_pd_qr, aes(x = depth_group, y = SR, fill = depth_group)) +
  geom_hline(yintercept = bact_qr_pd_sr_ref$SR, color = "steelblue" , linetype = "dotdash") +
  geom_violin(alpha = 0.3, show.legend = F) +
  geom_boxplot(width=0.1)+
  geom_pwc(ref.group = "0 - 60",
           method = "t.test",
           p.adjust.method = "BH",
           label = "p.adj.signif") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_fill_viridis_d(begin = 0.1, end = 0.9) +
  labs(x = "Stockpile depth (cm)",
       y = "Observed species richness",
       title = "16S rRNA gene")
qr_bact_SR



## NMDS on unweighed unifrac -----------------------------------------------




du_naf_bac <- GUniFrac(bact_naf_abund, tree = tree_naf_bact)$unifracs
unweighed_unifracs_naf_bac <- du_naf_bac[, , "d_UW"]

set.seed(1)
bact_nmds_naf_unifrac <- metaMDS(unweighed_unifracs_naf_bac, trymax = 999)
plot(bact_nmds_naf_unifrac)
text(bact_nmds_naf_unifrac)

plot_bact_naf_nmds <- data.frame(scores(bact_nmds_naf_unifrac)) %>% 
  merge(., bact_metadata, by = 'row.names', all = F)



du_qr_bact <- GUniFrac(bact_qr_abund, tree = tree_qr_bact)$unifracs
unweighed_unifracs_qr_bact <- du_qr_bact[, , "d_UW"]

set.seed(2)
bact_nmds_qr_unifrac <- metaMDS(unweighed_unifracs_qr_bact, trymax = 999)
plot(bact_nmds_qr_unifrac)
text(bact_nmds_qr_unifrac)

plot_bact_qr_nmds <- data.frame(scores(bact_nmds_qr_unifrac)) %>% 
  merge(., bact_metadata, by = 'row.names', all = F)

## Permanova ---------------------------------------------------------------

naf_meta_bact <- merge(bact_naf_abund, bact_metadata,
                       by = 'row.names')

naf_bray_bact <- naf_meta_bact %>% 
  select(!c(Row.names, site, sampling_depth,depth_group )) %>% 
  decostand(., method = "hellinger") %>% 
  vegdist(.)

set.seed(111)
mod_naf_bact <- adonis2(naf_bray_bact ~ naf_meta_bact$depth_group)
set.seed(990)
pairwise.adonis2(naf_bray_bact ~ depth_group, data = naf_meta_bact )
mod_naf_bact
plot(betadisper(naf_bray_bact, naf_meta_bact$depth_group))

qr_meta_bact <- merge(bact_qr_abund, bact_metadata,
                      by = 'row.names')

qr_bray_bact <- qr_meta_bact %>% 
  select(!c(Row.names, site, sampling_depth, depth_group)) %>% 
  decostand(., method = "hellinger") %>% 
  vegdist(.)

set.seed(111)
mod_qr_bact <- adonis2(qr_bray_bact ~ depth_group, qr_meta_bact)
mod_qr_bact
set.seed(2311)
pairwise.adonis2(qr_bray_bact ~ depth_group, data = qr_meta_bact )
plot(betadisper(qr_bray_bact, qr_meta_bact$depth_group))




## Random forest model -----------------------------------------------------


##new afton
bact_naf_rf_df <- bact_rf %>%
  filter(site == "NAF") %>% 
  select(!c( site, sampling_depth)) %>% 
  mutate(depth_group = as.factor(depth_group))

bact_rf_naf_names <-  bact_naf_rf_df %>% 
  select(-depth_group) %>% 
  select_if(colSums(.)>0) %>% 
  colnames()

bact_naf_rf <- bact_naf_rf_df %>% 
  select(depth_group, all_of(bact_rf_naf_names))

set.seed(112)
bact_mod_naf_rf <- ranger(depth_group ~., bact_naf_rf, classification = T, importance = "impurity_corrected")

bact_naf_vip <- vip::vip(bact_mod_naf_rf) +
  theme_bw(base_size = 18) +
  labs(title = "Relative importance of of 10 16S ASV at New Afton")

set.seed(54654654)
p_val_naf <- importance_pvalues(
  bact_mod_naf_rf,
  method = "janitza", formula = depth_group ~., bact_naf_rf,
  num.permutations = 999)

naf_imp_bact_asv <- p_val_naf %>%
  as.data.frame() %>% 
  # filter(pvalue < 0.05) %>%
  rownames_to_column(var = "asv_id")  %>% 
  arrange(desc(importance)) %>% 
  slice(1:10) %>% 
  #  rownames_to_column(., var = "asv_id") %>% 
  inner_join(., bact_asv, by = "asv_id") %>% 
  pull(hash_id)

bact_naf_imp <- bact %>% 
  select(all_of(naf_imp_bact_asv),depth_group) %>% 
  pivot_longer(-depth_group) %>% 
  inner_join(bact_asv, by = c("name" = "hash_id")) %>% 
  mutate(depth_group = factor(depth_group, levels = c("reference", "0 - 60", "61 - 150", "152 - 610", "611 - 1372"))) %>% 
  drop_na()

bact_naf_imp_asv_plot <- ggplot(bact_naf_imp, aes(x = depth_group, y = value, fill = depth_group)) +
  geom_boxplot() +
  facet_wrap(~asv_id, scales = "free") +
  scale_y_continuous(labels = scales::percent) +
   geom_pwc(ref.group = "reference",
            label = "p.signif",
            hide.ns = T) +
  theme_bw(base_size = 18) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.15),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        panel.grid = element_blank()) +
  scale_fill_viridis_d(name = "Stockpile depth (cm)") +
  labs(x = NULL,
       y = "Relative abundance of important \n 16S ASVs at New Afton")

bact_naf_imp_asv_plot
# 
# bact_asv %>% 
#   filter(hash_id %in% naf_imp_bact_asv) %>% 
#   inner_join(., tax_bact, by = c("hash_id" = "Feature.ID")) %>% 
#   view()

##QR mill
bact_qr_rf_df <- bact_rf %>%
  filter(site == "QR_mill") %>% 
  select(!c( site, sampling_depth)) %>% 
  mutate(depth_group = as.factor(depth_group))

bact_rf_qr_names <- bact_qr_rf_df %>% 
  select(-depth_group) %>% 
  select_if(colSums(.)>0) %>% 
  colnames()

bact_qr_rf <- bact_qr_rf_df %>% 
  select(depth_group, all_of(bact_rf_qr_names))

set.seed(1122)
bact_mod_qr_rf <- ranger(depth_group ~., bact_qr_rf, classification = T, importance = "impurity_corrected")

bact_qr_vip <- vip::vip(bact_mod_qr_rf) +
  theme_bw(base_size = 18)+
  labs(title = "Relative importance of top 10 16S ASV at QR mill")


set.seed(5465)
p_val_qr <- importance_pvalues(
  bact_mod_qr_rf,
  method = "janitza", formula = depth_group ~., bact_qr_rf,
  num.permutations = 999)

qr_imp_bact_asv <- p_val_qr %>%
  as.data.frame() %>% 
  # filter(pvalue < 0.05) %>%
  rownames_to_column(var = "asv_id")  %>% 
  arrange(desc(importance)) %>% 
  slice(1:10) %>% 
  #  rownames_to_column(., var = "asv_id") %>% 
  inner_join(., bact_asv, by = "asv_id") %>% 
  pull(hash_id)

bact_qr_imp <- bact %>% 
  select(all_of(qr_imp_bact_asv),depth_group) %>% 
  pivot_longer(-depth_group) %>% 
  inner_join(bact_asv, by = c("name" = "hash_id")) %>% 
  mutate(depth_group = factor(depth_group,
                              levels = c("reference", "0 - 60","75 - 260","350 - 390","500 - 575" ))) %>% 
  drop_na()

bact_qr_imp_asv_plot <- ggplot(bact_qr_imp, aes(x = depth_group, y = value, fill = depth_group)) +
  geom_boxplot() +
  facet_wrap(~asv_id, scales = "free") +
  scale_y_continuous(labels = scales::percent) +
  geom_pwc(ref.group = "reference",
           label = "p.signif",
           hide.ns = T) +
  theme_bw(base_size = 18) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.15),
        axis.text.x = element_text(angle = 45, vjust = 0.7),
        panel.grid = element_blank()) +
  scale_fill_viridis_d(name = "Stockpile depth (cm)") +
  labs(x = NULL,
       y = "Relative abundance of important \n16S ASVs at QR mill")

bact_qr_imp_asv_plot

bact_asv %>% 
  filter(hash_id %in% qr_imp_bact_asv) %>% 
  inner_join(., tax_bact, by = c("hash_id" = "Feature.ID"))

## Constrained ordination --------------------------------------------------

## New Afton
bact_chem_naf <- readxl::read_xlsx("../data/bact/na_chem_data.xlsx") %>% 
  drop_na() %>% 
  select(-depth_group) %>% 
  column_to_rownames(., var  = "sample_id") %>% 
  decostand(., method = "normalize", MARGIN = 2) %>% 
  #dplyr::select(c(Zn:NH4)) %>% 
  select(-c(Mn,P,NH4,C.N,B))


bact_sample_name_naf <- rownames(bact_chem_naf)

bact_naf_rda <- bact_naf_abund %>% 
  filter(row.names(.) %in% bact_sample_name_naf) %>% 
  decostand(., method = "hellinger")

bac_rda_name <- rownames(bact_naf_rda)

bact_chem_naf <- bact_chem_naf %>% 
  filter(row.names(.) %in% bac_rda_name)
#correlation::correlation(chem_naf) %>% summary() %>% plot()

set.seed(111)

mod_rdaa_bact_naf <- rda(bact_naf_rda ~.,bact_chem_naf)
naf_bact_constrain <- anova.cca(mod_rdaa_bact_naf, by = "term")
naf_bact_varname <- naf_bact_constrain %>% 
  as.data.frame() %>% 
  filter(row.names(.) != "Residual")
#filter(`Pr(>F)` < 0.05) %>% 
#row.names()

plot(mod_rdaa_bact_naf)



naf_contrain_ord_bact <- mod_rdaa_bact_naf$CCA$biplot[,1:2] %>% 
  as.data.frame() %>% 
  #  filter(row.names(.) %in% naf_bact_varname ) %>% 
  merge(., naf_bact_varname, by = "row.names") %>% 
  # rownames_to_column(., var = "var_name") %>% 
  mutate(line_type = `Pr(>F)` > 0.05  )






## QR mill
chem_qr_bact <- read.csv("../data/bact/qr_chem_data.csv") %>% 
  select(!c( H)) %>% 
  column_to_rownames(., var  = "sample_id")%>% 
  decostand(., method = "normalize", MARGIN = 2) 


sample_name_qr_bact <- rownames(chem_qr_bact)

bact_qr_rda <- bact_qr_abund %>% 
  filter(row.names(.) %in% sample_name_qr_bact) %>% 
  decostand(., method = "hellinger")

bact_qr_name <-row.names(bact_qr_rda)

chem_qr_bact <- chem_qr_bact %>% 
  filter(row.names(.) %in% bact_qr_name)


set.seed(1211)

mod_rdaa_bact_qr <- rda(bact_qr_rda ~.,chem_qr_bact )
qr_bact_constrain <- anova.cca(mod_rdaa_bact_qr, by = "term")
qr_bact_varname <- qr_bact_constrain %>% 
  as.data.frame() %>% 
  filter(row.names(.) != "Residual")

#plot(mod_rdaa_bact_qr)



bact_qr_contrain_ord <- mod_rdaa_bact_qr$CCA$biplot[,1:2] %>% 
  as.data.frame() %>% 
  merge(., qr_bact_varname, by = "row.names")%>% 
  mutate(line_type = `Pr(>F)` > 0.05  ) %>% 
  mutate(Row.names = str_remove(Row.names, pattern = ".1"))



## Betadiversity plot ------------------------------------------------------

bact_naf_beta <- plot_bact_naf_nmds %>% 
  mutate(depth_group = factor(depth_group, levels = c("reference", "0 - 60", "61 - 150", "152 - 610", "611 - 1372"))) %>% 
  ggplot(., aes(x = NMDS1, y = NMDS2, fill = depth_group)) +
  geom_segment(data = naf_contrain_ord_bact, aes(x = 0, xend = RDA1,
                                                 y = 0, yend = RDA2,
                                                 linetype = line_type), inherit.aes = F,
               color = "gray70",
               arrow = arrow(length = unit(0.03, "npc")),
               linewidth = 1,
               show.legend = F) +
  ggrepel::geom_text_repel(data = naf_contrain_ord_bact, aes(x = RDA1, y = RDA2, label = Row.names),
                           size = 5,
                           # fill = NA,
                           inherit.aes = F) +
  geom_point(pch=21, size = 5) +
  # stat_ellipse(geom = "polygon", alpha =0.3, show.legend = F) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "bottom") +
  scale_fill_viridis_d(name = "stockpile\ndepth (cm)",
                       begin = 0.1, end = 0.9) +
  annotate(geom = "text", label = glue::glue("R\U00B2 = {round(mod_naf_bact$R2[1],2)},  p < {mod_naf_bact$`Pr(>F)`[1]}"),
           x = -0.38, y = 0.35, size = 5, fontface = "italic") +
    labs(title = "16S rRNA gene")

bact_naf_beta


bact_qr_beta <- plot_bact_qr_nmds %>% 
  mutate(depth_group = factor(depth_group, levels = c("reference", "0 - 60","75 - 260","350 - 390","500 - 575" ))) %>% 
  ggplot(., aes(x = NMDS1, y = NMDS2, fill = depth_group)) +
  geom_segment(data = bact_qr_contrain_ord, aes(x = 0, xend = RDA1,
                                                y = 0, yend = RDA2,
                                                linetype = line_type), inherit.aes = F,
               color = "gray70",
               arrow = arrow(length = unit(0.03, "npc")),
               linewidth = 1,,
               show.legend = F) +
  geom_text(data = bact_qr_contrain_ord, aes(x = RDA1, y = RDA2, label = Row.names),
            size = 5,
            # fill = NA,
            inherit.aes = F) +
  geom_point(pch=21, size = 5) +
  # stat_ellipse(geom = "polygon", alpha =0.2, show.legend = F) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "right") +
  scale_fill_viridis_d(name = "Stockpile\ndepth (cm)",
                       begin = 0.1, end = 0.9) +
  annotate(geom = "text", label = glue::glue("R\U00B2 = {round(mod_qr_bact$R2[1],2)},  p < {mod_qr_bact$`Pr(>F)`[1]}"),
           x = -0.35, y = 0.40, size = 5, fontface = "italic") +
  labs(title = "16S rRNA gene")
 
bact_qr_beta

## Relative abundance ------------------------------------------------------


## New Afton

naf_bact_rel <- naf_meta_bact %>% 
  pivot_longer(-c(Row.names, site, sampling_depth, depth_group)) %>% 
  inner_join(., tax_bact, by = c("name" = "Feature.ID")) %>% 
  select(Row.names, site,  value, depth_group, Phylum) %>% 
  drop_na(Phylum) %>% 
  group_by(Row.names, site,depth_group, Phylum) %>% 
  summarise(val_sum = sum(value)) %>% 
  ungroup() %>% 
  group_by(site,depth_group, Phylum) %>% 
  summarise(mean_abund = mean(val_sum)) %>% 
  #mutate(depth_group = as.numeric(depth_group)) %>% 
  filter(Phylum != "p__bact_phy_Incertae_sedis",
         mean_abund >0.01) %>% 
  drop_na()

phylum_others_bact_naf <-  naf_bact_rel %>% 
  group_by(depth_group) %>% 
  summarize(mean_abund = sum(mean_abund)) %>% 
  mutate(mean_abund = 1- mean_abund,
         site = "NAF",
         Phylum = " Others")


naf_bact_rel_all <- bind_rows(naf_bact_rel, phylum_others_bact_naf) %>% 
  mutate(Phylum = str_remove(Phylum, pattern = "p__")) %>% 
  mutate(Phylum = str_remove_all(Phylum, pattern = "_.+")) %>% 
  filter(mean_abund >0) %>% 
  mutate(
    Phylum = factor(Phylum, levels = c(" Acidobacteriota"," Actinobacteriota"," Bacteroidota",
                                       " Chloroflexota", " Gemmatimonadota", " Proteobacteria", 
                                       " Firmicutes"," Thermoproteota"," Myxococcota", " Patescibacteria",
                                       " Others")),
    depth_group = factor(depth_group, levels = c("reference", "0 - 60","61 - 150", "152 - 610", "611 - 1372"))) %>% 
  mutate(phylum = if_else(Phylum != " Others", glue::glue("<i>{Phylum}</i>"), "Others"))



bact_naf <- ggplot(naf_bact_rel_all, aes(x = depth_group, y = mean_abund, fill = phylum)) +
  #geom_area(stat = "identity", color = "black", linewidth = 0.001,alpha = 0.7)+
  geom_bar(position="fill", stat="identity", color = "black", alpha = 0.7) +
  scale_fill_paletteer_d(name = "Phylum",`"miscpalettes::pastel"`)  +
  theme_bw(base_size = 22) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.key.height =  unit(5, 'mm') ,
        plot.title = element_text(hjust = 0.5),
        legend.text = element_markdown(size = 24),
        legend.title = element_text(size = 24) ) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Stockpile depth (cm)",
       y = "Relative abundance",
       title = "16S rRNA gene")

bact_naf

## QR mills

qr_bact_rel <- qr_meta_bact %>% 
  pivot_longer(-c(Row.names, site, sampling_depth, depth_group)) %>% 
  inner_join(., tax_bact, by = c("name" = "Feature.ID")) %>% 
  select(Row.names, site,  value, depth_group, Phylum) %>% 
  drop_na(Phylum) %>% 
  mutate(Phylum = str_remove(Phylum, pattern = "p__")) %>% 
  separate(Phylum, into = c("Phylum"), sep = "_") %>% 
  group_by(Row.names, site, depth_group, Phylum) %>% 
  summarise(val_sum = sum(value)) %>% 
  ungroup() %>% 
  group_by(site,depth_group, Phylum) %>% 
  summarise(mean_abund = mean(val_sum)) %>% 
  #mutate(depth_group = as.numeric(depth_group)) %>% 
  filter(mean_abund >0.01) 

phylum_others_bact <-  qr_bact_rel %>% 
  group_by(depth_group) %>% 
  summarize(mean_abund = sum(mean_abund)) %>% 
  mutate(mean_abund = 1- mean_abund,
         site = "QR_mill",
         Phylum = " Others")


qr_bact_rel_all <- bind_rows(qr_bact_rel, phylum_others_bact) %>% 
  mutate(Phylum = str_remove(Phylum, pattern = "p__"),
         # Phylum = str_remove(Phylum, pattern = "_.+")
  ) %>% 
  mutate(
    Phylum = factor(Phylum, levels = c(" Acidobacteriota"," Actinobacteriota"," Bacteroidota",
                                       " Chloroflexota"," Gemmatimonadota"," Proteobacteria",
                                       " Desulfobacterota"," Firmicutes"," Halobacteriota",
                                       " Myxococcota"," Spirochaetota", " Fibrobacterota",
                                       " Dormibacterota"," Others")),
    depth_group = factor(depth_group, levels = c("reference", "0 - 60","75 - 260","350 - 390","500 - 575" ))) %>% 
  mutate(phylum = if_else(Phylum != " Others", glue::glue("<i>{Phylum}</i>"), "Others"))



bact_qr <- ggplot(qr_bact_rel_all, aes(x = depth_group, y = mean_abund, fill = phylum)) +
  #geom_area(stat = "identity", color = "black", linewidth = 0.001,alpha = 0.7)+
  geom_bar(position="fill", stat="identity", color = "black", alpha = 0.7) +
  scale_fill_paletteer_d(name = "Phylum",`"miscpalettes::pastel"`) +
  #scale_fill_viridis_d()+
  theme_bw(base_size = 22) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.key.height =  unit(5, 'mm') ,
        plot.title = element_text(hjust = 0.5),
        legend.text = element_markdown(size = 24),
        legend.title = element_text(size = 24)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Stockpile depth (cm)",
       y = "Relative abundance",
       title = "16S rRNA gene")
bact_qr




# Fungal analyses ---------------------------------------------------------

## Data import and wrangling -----------------------------------------------



fungi_qza <- read_qza("../data/fungi/ashley-fungal.featuretable.qza")

fungi_count <- fungi_qza[["data"]] %>% 
  t(.) %>% 
  as.data.frame(.)%>% 
  filter(rowSums(.)>1000) %>% 
  select_if(colSums(.)>50) %>% 
  adorn_percentages(,, `bceb0fc1e2e85f4ba9a0253dd5811da0`:`6590a3fe2ac7a6c502d552c10fd2c66a`) %>% 
  # adorn_percentages(,, `bceb0fc1e2e85f4ba9a0253dd5811da0`:`91322b4090368abf952873c2c493c3ec`) %>% 
  filter(!row.names(.) %in% c("SC4-01.5a", "NA-R3"))




samp_name_fungi <- rownames(fungi_count)

fungi_metadata <- readxl::read_xlsx("../data/fungi/fungi_metadata.xlsx") %>% 
  mutate(depth_group = case_when( sampling_depth == "reference" ~ "reference",
                                  site == "QR_mill" & sampling_depth %in% c("10","20", "60")  ~ "0 - 60",
                                  #  site == "QR_mill" & sampling_depth ==  "20" ~ "11-20",
                                  #  site == "QR_mill" & sampling_depth %in% c( "75", "90","105", "108") ~ "75-120",
                                  site == "QR_mill" & sampling_depth %in% c("75", "90","105", "108", "200", "210","250") ~ "75 - 260",
                                  site == "QR_mill" & sampling_depth %in% c("350", "365", "380") ~ "350 - 390",
                                  site == "QR_mill" & sampling_depth %in% c("500", "550", "565") ~ "500 - 575",
                                  site == "NAF" & sampling_depth %in% c("30", "60")  ~ "0 - 60",
                                  site == "NAF" & sampling_depth %in% c("90", "120", "150")  ~ "61 - 150",
                                  site == "NAF" & sampling_depth %in% c("152", "300","460","610"  )  ~ "152 - 610",
                                  site == "NAF" & sampling_depth %in% c("760", "910", "1070", "1220","1370","610"  )  ~ "611 - 1372")) %>% 
  drop_na() %>% 
  column_to_rownames(., var = "Sample_id")


fungi <- merge(fungi_count, fungi_metadata,
               by = 'row.names', all = F) %>% 
  column_to_rownames(., var = "Row.names")

fungi_asv <- data.frame(hash_id = colnames(fungi)) %>% 
  mutate(asv_id = case_when(hash_id == "site" ~ "site",
                            hash_id == "sampling_depth" ~ "sampling_depth",
                            hash_id == "depth_group" ~ "depth_group",
                            TRUE  ~  paste("ASV", row_number(), sep = "_"))
  )

fungi_rf <- fungi
colnames(fungi_rf) <- fungi_asv$asv_id


tree_tips_fungi <- fungi %>% 
  select(!c(site, sampling_depth, depth_group)) %>% 
  colnames()

tree_fungi <- read_qza("../data/fungi/ashley-fungal.rooted-tree.qza")[["data"]] %>% 
  keep.tip(., tip = tree_tips_fungi )


fungi_tax <- read_qza("../data/fungi/ashley-fungal.taxonomy.qza")


tax_fungi <- fungi_tax[["data"]] %>% 
  as.data.frame() %>% 
  filter(Feature.ID %in% tree_tips_fungi) %>% 
  # column_to_rownames(., var = "Feature.ID") %>% 
  select(!Confidence) %>% 
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") 


fungi_naf_abund <- fungi %>%
  filter(site == "NAF") %>% 
  select(!c( site, sampling_depth, depth_group)) %>% 
  select_if(colSums(.)>0) 


fungi_qr_abund <- fungi %>%
  filter(site == "QR_mill") %>% 
  select(!c( site, sampling_depth, depth_group))%>% 
  select_if(colSums(.)>0) 

fungi_asv_naf <- colnames(fungi_naf_abund)
fungi_asv_qr <- colnames(fungi_qr_abund)


tree_naf_fungi <- tree_fungi %>% 
  keep.tip(., tip = fungi_asv_naf)

tree_qr_fungi <- tree_fungi %>% 
  keep.tip(., tip = fungi_asv_qr)


## Alpha diversity analysis ------------------------------------------------

#faith's PD analysis

fungi_faiths_pd_naf <- pd(fungi_naf_abund, tree = tree_naf_fungi) %>% 
  merge(., fungi_metadata, by = 'row.names', all = F) %>% 
  #mutate(sampling_depth = as.numeric(sampling_depth)) %>% 
  filter(depth_group != "reference") %>% 
  mutate(depth_group = factor(depth_group, levels = c("reference", "0 - 60", "61 - 150", "152 - 610", "611 - 1372")))

fungi_naf_pd_sr_ref <- pd(fungi_naf_abund, tree = tree_naf_fungi) %>% 
  merge(., fungi_metadata, by = 'row.names', all = F) %>% 
  filter(depth_group == "reference") %>% 
  group_by(site) %>% 
  summarise(across(c(PD,SR), mean))

fungi_mod_pd_naf <- lm(PD ~ depth_group, fungi_faiths_pd_naf)
summary.aov(fungi_mod_pd_naf)

fungi_mod_SR_naf <- lm(SR ~ depth_group, fungi_faiths_pd_naf)
summary.aov(fungi_mod_SR_naf)


naf_fungal_pd <- ggplot(fungi_faiths_pd_naf, aes(x = depth_group, y = PD, fill = depth_group)) +
  geom_hline(yintercept = fungi_naf_pd_sr_ref$PD, color = "steelblue" , linetype = "dotdash") +
  geom_violin(alpha = 0.3, show.legend = F) +
  geom_boxplot(width=0.2)+
  geom_pwc(ref.group = "0 - 60",
           method = "t.test",
           p.adjust.method = "BH",
           label = "p.adj.signif") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_fill_viridis_d(begin = 0.1, end = 0.9) +
  labs(x = "Stockpile depth (cm)",
       y = "Faith's phylogenetic diveristy",
       title = "ITS gene")
naf_fungal_pd

naf_fungal_SR <- ggplot(fungi_faiths_pd_naf, aes(x = depth_group, y = SR, fill = depth_group)) +
  geom_hline(yintercept = fungi_naf_pd_sr_ref$SR, color = "steelblue" , linetype = "dotdash") +
  geom_violin(alpha = 0.3, show.legend = F) +
  geom_boxplot(width=0.2)+
  geom_pwc(ref.group = "0 - 60",
           method = "t.test",
           p.adjust.method = "BH",
           label = "p.adj.signif") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_fill_viridis_d(begin = 0.1, end = 0.9) +
  labs(x = "Stockpile depth (cm)",
       y = "Observed species richness",
       title = "ITS gene")
naf_fungal_SR


fungi_faiths_pd_qr <- pd(fungi_qr_abund, tree = tree_qr_fungi)%>% 
  merge(., fungi_metadata, by = 'row.names', all = F)  %>% 
  filter(depth_group != "reference") %>% 
  mutate(depth_group = factor(depth_group, levels = c("reference", "0 - 60","75 - 260","350 - 390","500 - 575" )))

fungi_qr_pd_sr_ref <- pd(fungi_qr_abund, tree = tree_qr_fungi) %>% 
  merge(., fungi_metadata, by = 'row.names', all = F) %>% 
  filter(depth_group == "reference") %>% 
  group_by(site) %>% 
  summarise(across(c(PD,SR), mean))

# mod_pd_qr <- lm(PD ~ depth_group, fungi_faiths_pd_qr)
# summary.aov(mod_pd_qr)
# 
# mod_SR_qr <- lm(SR ~ depth_group, fungi_faiths_pd_qr)
# summary.aov(mod_SR_qr)


qr_fungal_pd <- ggplot(fungi_faiths_pd_qr, aes(x = depth_group, y = PD, fill = depth_group)) +
  geom_hline(yintercept = fungi_qr_pd_sr_ref$PD, color = "steelblue" , linetype = "dotdash") +
  geom_violin(alpha = 0.3, show.legend = F) +
  geom_boxplot(width=0.1)+
  geom_pwc(ref.group = "0 - 60",
           method = "t.test",
           p.adjust.method = "BH",
           label = "p.adj.signif") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_fill_viridis_d(begin = 0.1, end = 0.9) +
  labs(x = "Stockpile depth (cm)",
       y = "Faith's phylogenetic diveristy",
       title = "ITS gene")  +
  ylim(c(0, 170))
qr_fungal_pd

qr_fungal_SR <- ggplot(fungi_faiths_pd_qr, aes(x = depth_group, y = SR, fill = depth_group)) +
  geom_hline(yintercept = fungi_qr_pd_sr_ref$SR, color = "steelblue" , linetype = "dotdash") +
  geom_violin(alpha = 0.3, show.legend = F) +
  geom_boxplot(width=0.1)+
  geom_pwc(ref.group = "0 - 60",
           method = "t.test",
           p.adjust.method = "BH",
           label = "p.adj.signif") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_fill_viridis_d(begin = 0.1, end = 0.9) +
  labs(x = "Stockpile depth (cm)",
       y = "Observed species richness",
       title = "ITS gene") +
  ylim(c(0, 1100))
qr_fungal_SR



## NMDS on unweighed unifrac -----------------------------------------------




du_naf_fungi <- GUniFrac(fungi_naf_abund, tree = tree_naf_fungi)$unifracs
unweighed_unifracs_naf_fungi <- du_naf_fungi[, , "d_UW"]

set.seed(1)
fungi_nmds_naf_unifrac <- metaMDS(unweighed_unifracs_naf_fungi, trymax = 9999)
plot(fungi_nmds_naf_unifrac)
text(fungi_nmds_naf_unifrac)

plot_fungi_naf_nmds <- data.frame(scores(fungi_nmds_naf_unifrac)) %>% 
  merge(., fungi_metadata, by = 'row.names', all = F)



du_qr_fungi <- GUniFrac(fungi_qr_abund, tree = tree_qr_fungi)$unifracs
unweighed_unifracs_qr_fungi <- du_qr_fungi[, , "d_UW"]

set.seed(2)
fungi_nmds_qr_unifrac <- metaMDS(unweighed_unifracs_qr_fungi, trymax = 999)
plot(fungi_nmds_qr_unifrac)
text(fungi_nmds_qr_unifrac)

plot_fungi_qr_nmds <- data.frame(scores(fungi_nmds_qr_unifrac)) %>% 
  merge(., fungi_metadata, by = 'row.names', all = F)

## Permanova ---------------------------------------------------------------

naf_meta_fungi <- merge(fungi_naf_abund, fungi_metadata,
                  by = 'row.names')

naf_bray_fungi <- naf_meta_fungi %>% 
  select(!c(Row.names, site, sampling_depth,depth_group )) %>% 
  decostand(., method = "hellinger") %>% 
  vegdist(.)

set.seed(111)
mod_naf_fungi <- adonis2(naf_bray_fungi ~ depth_group, naf_meta_fungi)
mod_naf_fungi
set.seed(231)
pairwise.adonis2(naf_bray_fungi ~ depth_group, naf_meta_fungi)
plot(betadisper(naf_bray_fungi, naf_meta_fungi$depth_group))

qr_meta_fungi <- merge(fungi_qr_abund, fungi_metadata,
                 by = 'row.names')

qr_bray_fungi <- qr_meta_fungi %>% 
  select(!c(Row.names, site, sampling_depth, depth_group)) %>% 
  decostand(., method = "hellinger") %>% 
  vegdist(.)

set.seed(111)
mod_qr_fungi <- adonis2(qr_bray_fungi ~ depth_group, qr_meta_fungi)
mod_qr_fungi
set.seed(12221)
pairwise.adonis2(qr_bray_fungi ~ depth_group, qr_meta_fungi)
plot(betadisper(qr_bray_fungi, qr_meta_fungi$depth_group))




## Random forest model -----------------------------------------------------


##new afton
fungi_naf_rf_df <- fungi_rf %>%
  filter(site == "NAF") %>% 
  select(!c( site, sampling_depth)) %>% 
  mutate(depth_group = as.factor(depth_group))

rf_naf_names_fungi <- fungi_naf_rf_df %>% 
  select(-depth_group) %>% 
  select_if(colSums(.)>0) %>% 
  colnames()

fungi_naf_rf <- fungi_naf_rf_df %>% 
  select(depth_group, all_of(rf_naf_names_fungi))

set.seed(11223)
mod_naf_rf_fungi <- ranger(depth_group ~., fungi_naf_rf, classification = T, importance = "impurity_corrected")

fungi_naf_vip <- vip::vip(mod_naf_rf_fungi) +
  theme_bw(base_size = 18)+
  labs(title = "Relative importance of top 10 ITS ASV at New Afton")

set.seed(54654654)
p_val_naf_fungi <- importance_pvalues(
  mod_naf_rf_fungi,
  method = "janitza", formula = depth_group ~., fungi_naf_rf,
  num.permutations = 999)

naf_imp_fungal_asv <- p_val_naf_fungi %>%
  as.data.frame() %>% 
  # filter(pvalue < 0.05) %>%
  rownames_to_column(var = "asv_id")  %>% 
  arrange(desc(importance)) %>% 
  slice(1:10) %>% 
  #  rownames_to_column(., var = "asv_id") %>% 
  inner_join(., fungi_asv, by = "asv_id") %>% 
  pull(hash_id)

fungi_naf_imp <- fungi %>% 
  select(all_of(naf_imp_fungal_asv),depth_group) %>% 
  pivot_longer(-depth_group) %>% 
  inner_join(fungi_asv, by = c("name" = "hash_id")) %>% 
  mutate(depth_group = factor(depth_group, levels = c("reference", "0 - 60", "61 - 150", "152 - 610", "611 - 1372"))) %>% 
  drop_na()

fungi_naf_imp_asv_plot <- ggplot(fungi_naf_imp, aes(x = depth_group, y = value, fill = depth_group)) +
  geom_boxplot() +
  facet_wrap(~asv_id, scales = "free") +
  scale_y_continuous(labels = scales::percent) +
   geom_pwc(ref.group = "reference",
            label = "p.signif",
            hide.ns = T) +
  theme_bw(base_size = 18) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.15),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.7)) +
  scale_fill_viridis_d(name = "Stockpile depth (cm)") +
  labs(x = NULL,
       y = "Relative abundance of important \nITS ASVs at New Afton")

fungi_naf_imp_asv_plot

fungi_asv %>% 
  filter(hash_id %in% naf_imp_fungal_asv) %>% 
  inner_join(., tax_fungi, by = c("hash_id" = "Feature.ID"))

##QR mill
fungi_qr_rf_df <- fungi_rf %>%
  filter(site == "QR_mill") %>% 
  select(!c( site, sampling_depth)) %>% 
  mutate(depth_group = as.factor(depth_group))

fungi_rf_qr_names <- fungi_qr_rf_df %>% 
  select(-depth_group) %>% 
  select_if(colSums(.)>0) %>% 
  colnames()

fungi_qr_rf <- fungi_qr_rf_df %>% 
  select(depth_group, all_of(fungi_rf_qr_names))

set.seed(1122334)
fungi_mod_qr_rf <- ranger(depth_group ~., fungi_qr_rf, classification = T, importance = "impurity_corrected")
fungi_qr_vip <- vip::vip(fungi_mod_qr_rf) +
  theme_bw(base_size = 18) +
  labs(title = "Relative importance of top 10 ITS ASV at QR")

set.seed(5465)
p_val_qr_fungi <- importance_pvalues(
  fungi_mod_qr_rf,
  method = "janitza", formula = depth_group ~., fungi_qr_rf,
  num.permutations = 999)

qr_imp_fungal_asv <- p_val_qr_fungi %>%
  as.data.frame() %>% 
  # filter(pvalue < 0.05) %>%
  rownames_to_column(var = "asv_id")  %>% 
  arrange(desc(importance)) %>% 
  slice(1:10) %>% 
  #  rownames_to_column(., var = "asv_id") %>% 
  inner_join(., fungi_asv, by = "asv_id") %>% 
  pull(hash_id)

fungi_qr_imp <- fungi %>% 
  select(all_of(qr_imp_fungal_asv),depth_group) %>% 
  pivot_longer(-depth_group) %>% 
  inner_join(fungi_asv, by = c("name" = "hash_id")) %>% 
  mutate(depth_group = factor(depth_group,
                              levels = c("reference", "0 - 60","75 - 260","350 - 390","500 - 575" ))) %>% 
  drop_na()

fungi_qr_imp_asv_plot <- ggplot(fungi_qr_imp, aes(x = depth_group, y = value, fill = depth_group)) +
  geom_boxplot() +
  facet_wrap(~asv_id, scales = "free") +
  scale_y_continuous(labels = scales::percent) +
  geom_pwc(ref.group = "reference",
           label = "p.signif",
           hide.ns = T) +
  theme_bw(base_size = 18) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.15),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.7)) +
  scale_fill_viridis_d(name = "Stockpile depth (cm)") +
  labs(x = NULL,
       y = "Relative abundance of important \n ITS ASVs at QR mills")

fungi_qr_imp_asv_plot

fungi_asv %>% 
  filter(hash_id %in% qr_imp_fungal_asv) %>% 
  inner_join(., tax_fungi, by = c("hash_id" = "Feature.ID"))

## Constrained ordination --------------------------------------------------

## New Afton
chem_naf_fungi <- readxl::read_xlsx("../data/fungi/na_chem_data.xlsx") %>% 
  select(-depth_group) %>% 
  column_to_rownames(., var  = "sample_id") %>% 
  drop_na() %>% 
   decostand(., method = "normalize", MARGIN = 2) %>% 
  select(-NH4.1, B)



sample_name_naf_fungi <- rownames(chem_naf_fungi)

fungi_naf_rda <- fungi_naf_abund %>% 
  filter(row.names(.) %in% sample_name_naf_fungi) %>% 
  decostand(., method = "hellinger")


set.seed(111)

mod_rdaa_fungi_naf <- rda(fungi_naf_rda ~.,chem_naf_fungi )
naf_fungi_constrain <- anova.cca(mod_rdaa_fungi_naf, by = "term")
naf_fungi_varname <- naf_fungi_constrain %>% 
  as.data.frame() %>% 
  filter(row.names(.) != "Residual")

#plot(naf_fungi_constrain)



fungi_naf_contrain_ord <- mod_rdaa_fungi_naf$CCA$biplot[,1:2] %>% 
  as.data.frame() %>% 
  merge(., naf_fungi_varname, by = "row.names") %>% 
  mutate(Row.names = str_remove(Row.names, pattern = ".1"))%>% 
  mutate(line_type = `Pr(>F)` > 0.05  )





## QR mill
  chem_qr_fungi <- read.csv("../data/fungi/qr_chem_data.csv") %>% 
  select(!c(S.1, H)) %>% 
  column_to_rownames(., var  = "sample_id")%>% 
  decostand(., method = "normalize", MARGIN = 2) 


sample_name_qr_fungi <- rownames(chem_qr_fungi)

fungi_qr_rda <- fungi_qr_abund %>% 
  filter(row.names(.) %in% sample_name_qr_fungi) %>% 
  decostand(., method = "hellinger")

mod_rdaa_fungi_qr <- rda(fungi_qr_rda ~.,chem_qr_fungi )
qr_fungi_constrain <- anova.cca(mod_rdaa_fungi_qr, by = "term")
qr_fungi_varname <- qr_fungi_constrain %>% 
  as.data.frame() %>% 
  filter(row.names(.) != "Residual")

#plot(mod_rdaa_fungi_qr)



fungi_qr_contrain_ord <- mod_rdaa_fungi_qr$CCA$biplot[,1:2] %>% 
  as.data.frame() %>% 
  merge(., qr_fungi_varname, by = "row.names")%>% 
  mutate(Row.names = str_remove(Row.names, pattern = ".1")) %>% 
  mutate(line_type = `Pr(>F)` > 0.05  )



## Betadiversity plot ------------------------------------------------------

fungi_naf_beta <- plot_fungi_naf_nmds %>% 
  mutate(depth_group = factor(depth_group, levels = c("reference", "0 - 60", "61 - 150", "152 - 610", "611 - 1372"))) %>% 
  ggplot(., aes(x = NMDS1, y = NMDS2, fill = depth_group)) +
  geom_segment(data = fungi_naf_contrain_ord, aes(x = 0, xend = RDA1,
                                                  y = 0, yend = RDA2,
                                                  linetype = line_type), inherit.aes = F,
               color = "gray70",
               arrow = arrow(length = unit(0.03, "npc")),
               linewidth = 1,
               show.legend = F) +
  ggrepel::geom_text_repel(data = fungi_naf_contrain_ord, aes(x = RDA1, y = RDA2, label = Row.names),
                           size = 5,
                           # fill = NA,
                           inherit.aes = F,
                           max.overlaps = Inf) +
  geom_point(pch=21, size = 5) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "bottom") +
  scale_fill_viridis_d(name = "stockpile\ndepth (cm)",
                       begin = 0.1, end = 0.9) +
  annotate(geom = "text", label = glue::glue("R\U00B2 = {round(mod_naf_fungi$R2[1],2)},  p < {mod_naf_fungi$`Pr(>F)`[1]}"),
           x = 0.25, y = 0.40, size = 5, fontface = "italic") +
  
  labs(title = "ITS gene") +
  ylim(c(NA, 0.4))

fungi_naf_beta


fungi_qr_beta <- plot_fungi_qr_nmds %>% 
  mutate(depth_group = factor(depth_group, levels = c("reference", "0 - 60","75 - 260","350 - 390","500 - 575" ))) %>% 
  ggplot(., aes(x = NMDS1, y = NMDS2, fill = depth_group)) +
  geom_segment(data = fungi_qr_contrain_ord, aes(x = 0, xend = RDA1,
                                                 y = 0, yend = RDA2,
                                                 linetype = line_type), inherit.aes = F,
               color = "gray70",
               arrow = arrow(length = unit(0.03, "npc")),
               linewidth = 1,
               show.legend = F) +
  ggrepel::geom_text_repel(data = fungi_qr_contrain_ord, aes(x = RDA1+0.03, y = RDA2, label = Row.names),
                           size = 5,
                           inherit.aes = F) +
   geom_point(pch=21, size = 5) +
  # stat_ellipse(geom = "polygon", alpha =0.2, show.legend = F) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "bottom") +
  scale_fill_viridis_d(name = "Stockpile\ndepth (cm)",
                       begin = 0.1, end = 0.9) +
  annotate(geom = "text", label = glue::glue("R\U00B2 = {round(mod_qr_fungi$R2[1],2)},  p < {mod_qr_fungi$`Pr(>F)`[1]}"),
           x = 0.45, y = 0.75, size = 5, fontface = "italic") +
    labs(title = "ITS gene")

fungi_qr_beta

## Relative abundance ------------------------------------------------------


## New Afton

naf_fungi_rel <- naf_meta_fungi %>% 
  pivot_longer(-c(Row.names, site, sampling_depth, depth_group)) %>% 
  inner_join(., tax_fungi, by = c("name" = "Feature.ID")) %>% 
  select(Row.names, site,  value, depth_group, Phylum) %>% 
  drop_na(Phylum) %>% 
  group_by(Row.names, site,depth_group, Phylum) %>% 
  summarise(val_sum = sum(value)) %>% 
  ungroup() %>% 
  group_by(site,depth_group, Phylum) %>% 
  summarise(mean_abund = mean(val_sum)) %>% 
  #mutate(depth_group = as.numeric(depth_group)) %>% 
  filter(Phylum != "p__Fungi_phy_Incertae_sedis",
         mean_abund >0.01) 

fungi_phylum_others <-  naf_fungi_rel %>% 
  group_by(depth_group) %>% 
  summarize(mean_abund = sum(mean_abund)) %>% 
  mutate(mean_abund = 1- mean_abund,
         site = "NAF",
         Phylum = "Others")

naf_fungi_rel_all <- bind_rows(naf_fungi_rel, fungi_phylum_others) %>% 
  mutate(Phylum = str_remove(Phylum, pattern = "p__")) %>% 
  filter(mean_abund >0) %>% 
  mutate(Phylum = factor(Phylum, levels = c("Ascomycota",
                                            "Basidiomycota",
                                            "Mortierellomycota",
                                            "Rozellomycota",
                                            "Others")),
         depth_group = factor(depth_group, levels = c("reference", "0 - 60","61 - 150", "152 - 610", "611 - 1372"))) %>% 
  mutate(phylum = if_else(Phylum != "Others", glue::glue("<i>{Phylum}</i>"), "Others"))



fungal_naf_rel <- ggplot(naf_fungi_rel_all, aes(x = depth_group, y = mean_abund, fill = phylum)) +
  #geom_area(stat = "identity", color = "black", linewidth = 0.001,alpha = 0.7)+
  geom_bar(position="fill", stat="identity", color = "black", alpha = 0.7) +
  scale_fill_paletteer_d(name = "Phylum",`"miscpalettes::pastel"`) +
  theme_bw(base_size = 22) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.key.height =  unit(5, 'mm') ,
        plot.title = element_text(hjust = 0.5),
        legend.text = element_markdown(size = 24),
        legend.title = element_text(size = 24)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Stockpile depth (cm)",
       y = "Relative abundance",
       title = "ITS gene")

fungal_naf_rel
## QR mills

qr_fungi_rel <- qr_meta_fungi %>% 
  pivot_longer(-c(Row.names, site, sampling_depth, depth_group)) %>% 
  inner_join(., tax_fungi, by = c("name" = "Feature.ID")) %>% 
  select(Row.names, site,  value, depth_group, Phylum) %>% 
  drop_na(Phylum) %>% 
  group_by(Row.names, site,depth_group, Phylum) %>% 
  summarise(val_sum = sum(value)) %>% 
  ungroup() %>% 
  group_by(site,depth_group, Phylum) %>% 
  summarise(mean_abund = mean(val_sum)) %>% 
  #mutate(depth_group = as.numeric(depth_group)) %>% 
  filter(Phylum != "p__Fungi_phy_Incertae_sedis",
         mean_abund >0.01) %>% 
  mutate(phylum = glue::glue("<i>{Phylum},</i>"))

fungi_phylum_others <-  qr_fungi_rel %>% 
  group_by(depth_group) %>% 
  summarize(mean_abund = sum(mean_abund)) %>% 
  mutate(mean_abund = 1- mean_abund,
         site = "QR_mill",
         Phylum = "Others")


qr_fungi_rel_all <- bind_rows(qr_fungi_rel, fungi_phylum_others) %>% 
  mutate(Phylum = str_remove(Phylum, pattern = "p__")) %>% 
  mutate(Phylum = factor(Phylum, levels = c("Ascomycota",
                                            "Basidiomycota",
                                            "Mortierellomycota",
                                            "Rozellomycota",
                                            "Others")),
         depth_group = factor(depth_group, levels = c("reference", "0 - 60","75 - 260","350 - 390","500 - 575" ))) %>% 
  mutate(phylum = if_else(Phylum != "Others", glue::glue("<i>{Phylum}</i>"), "Others"))



fungal_qr <- ggplot(qr_fungi_rel_all, aes(x = depth_group, y = mean_abund, fill = phylum)) +
  #geom_area(stat = "identity", color = "black", linewidth = 0.001,alpha = 0.7)+
  geom_bar(position="fill", stat="identity", color = "black", alpha = 0.7) +
  scale_fill_paletteer_d(name = "Phylum",`"miscpalettes::pastel"`) +
  theme_bw(base_size = 22) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.key.height =  unit(5, 'mm') ,
        plot.title = element_text(hjust = 0.5),
        legend.text = element_markdown(size = 24),
        legend.title = element_text(size = 24)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Stockpile depth (cm)",
       y = "Relative abundance",
       title = "ITS gene")

fungal_qr






# Save figures ------------------------------------------------------------

library(ggpubr)

fig1 <- ggarrange(bact_naf, fungal_naf_rel, bact_qr, fungal_qr,  labels = "AUTO", font.label = list(size = 20))
ggsave(fig1, file = "../plots/fig1.jpg", height = 18, width = 24, units = "in", dpi = 800)



fig2 <- ggarrange(naf_bact_pd, naf_fungal_pd, naf_bact_SR, naf_fungal_SR,
                  ncol = 2, nrow = 2, labels = "AUTO", font.label = list(size = 20))
ggsave(fig2, file = "../plots/fig2_naf.jpg", height = 9, width = 12, units = "in", dpi = 800)



fig3 <- ggarrange(qr_bact_pd, qr_fungal_pd, qr_bact_SR, qr_fungal_SR,
                  ncol = 2, nrow = 2, labels = "AUTO", font.label = list(size = 20))
ggsave(fig3, file = "../plots/fig3_qr.jpg", height = 9, width = 12, units = "in", dpi = 800)



fig4 <- ggarrange(bact_naf_beta, fungi_naf_beta, labels = "AUTO", font.label = list(size = 16),
                  common.legend = T, legend = "bottom")
ggsave(fig4, file = "../plots/fig4_naf.jpg", height = 6, width = 12, units = "in", dpi = 800)                 



fig5 <- ggarrange(bact_qr_beta, fungi_qr_beta, labels = "AUTO", font.label = list(size = 16),
                  common.legend = T, legend = "bottom")
ggsave(fig5, file = "../plots/fig5_qr.jpg", height = 6, width = 12, units = "in", dpi = 800)     



# Supplementary figures ---------------------------------------------------




ggsave(bact_naf_imp_asv_plot, file = "../plots/Fig_SI1.jpg", height = 14, width = 14, units = "in", dpi = 800) 



ggsave(bact_qr_imp_asv_plot, file = "../plots/Fig_SI2.jpg", height = 14, width = 14, units = "in", dpi = 800) 



ggsave(fungi_naf_imp_asv_plot, file = "../plots/Fig_SI3.jpg", height = 14, width = 14, units = "in", dpi = 800) 



ggsave(fungi_qr_imp_asv_plot, file = "../plots/Fig_SI4.jpg", height = 14, width = 14, units = "in", dpi = 800) 



fig_si5 <- ggarrange(bact_naf_vip, bact_qr_vip, fungi_naf_vip, fungi_qr_vip,  labels = "AUTO", font.label = list(size = 20))
ggsave(fig_si5, file = "../plots/Fig_SI5.jpg", height = 18, width = 24, units = "in", dpi = 800)



fig_SI6 <- naf_meta_bact %>% 
  pivot_longer(-c(Row.names, site, sampling_depth, depth_group)) %>% 
  inner_join(., tax_bact, by = c("name" = "Feature.ID")) %>% 
  select(Row.names, site,  value, depth_group, Phylum) %>% 
  drop_na(Phylum) %>% 
  group_by(Row.names, site,depth_group, Phylum) %>% 
  summarise(val_sum = sum(value)) %>% 
  # filter(val_sum > 0.01) %>% 
  mutate(Phylum = str_remove(Phylum, pattern = "p__")) %>% 
  separate(Phylum, into = c("Phylum"), sep = "_") %>% 
  filter(Phylum %in% c(" Chloroflexota",  " Firmicutes")) %>%  
  mutate(depth_group = factor(depth_group, levels = c("reference", "0 - 60", "61 - 150", "152 - 610", "611 - 1372"))) %>% 
  # filter(depth_group != "reference") %>% 
  #filter(Phylum == " Firmicutes") %>% 
  ggplot(., aes(x = depth_group, y = val_sum, fill = depth_group)) +
  geom_boxplot()+
  geom_pwc(ref.group = "reference",
           label = "p.signif") +
  facet_wrap(~ Phylum, scales = "free") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "bottom") +
  scale_fill_viridis_d(name = "stockpile\ndepth (cm)",
                       begin = 0.1, end = 0.9) +
  labs(x = "Stockpile depth (cm)",
       y = "Relative abundance") +
  scale_y_continuous(labels = scales::percent)

ggsave(fig_SI6, file = "../plots/fig_SI6.jpg", height = 6, width = 12, units = "in", dpi = 800)


Fig_SI7 <- qr_meta_bact %>% 
  pivot_longer(-c(Row.names, site, sampling_depth, depth_group)) %>% 
  inner_join(., tax_bact, by = c("name" = "Feature.ID")) %>% 
  select(Row.names, site,  value, depth_group, Phylum) %>% 
  drop_na(Phylum) %>% 
  group_by(Row.names, site,depth_group, Phylum) %>% 
  summarise(val_sum = sum(value)) %>% 
  # filter(val_sum > 0.01) %>% 
  mutate(Phylum = str_remove(Phylum, pattern = "p__")) %>% 
  separate(Phylum, into = c("Phylum"), sep = "_") %>% 
  filter(Phylum %in% c(" Acidobacteriota",  " Actinobacteriota"," Desulfobacterota",
                       " Halobacteriota"," Verrucomicrobiota", " Proteobacteria" )) %>%  
  mutate(depth_group = factor(depth_group, levels = c("reference", "0 - 60","75 - 260","350 - 390","500 - 575" ))) %>% 
  # filter(depth_group != "reference") %>% 
  #filter(Phylum == " Firmicutes") %>% 
  ggplot(., aes(x = depth_group, y = val_sum, fill = depth_group)) +
  geom_boxplot()+
  geom_pwc(ref.group = "reference",
           label = "p.signif") +
  facet_wrap(~ Phylum, scales = "free") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "bottom") +
  scale_fill_viridis_d(name = "stockpile\ndepth (cm)",
                       begin = 0.1, end = 0.9) +
  labs(x = "Stockpile depth (cm)",
       y = "Relative abundance") +
  scale_y_continuous(labels = scales::percent)

ggsave(Fig_SI7, file = "../plots/Fig_SI7.jpg", height = 10, width =17, units = "in", dpi = 800)


fig_SI8 <- naf_meta_fungi %>% 
  pivot_longer(-c(Row.names, site, sampling_depth, depth_group)) %>% 
  inner_join(., tax_fungi, by = c("name" = "Feature.ID")) %>% 
  select(Row.names, site,  value, depth_group, Phylum) %>% 
  drop_na(Phylum) %>% 
  group_by(Row.names, site,depth_group, Phylum) %>% 
  summarise(val_sum = sum(value)) %>% 
  # filter(val_sum > 0.01) %>% 
  mutate(Phylum = str_remove(Phylum, pattern = "p__")) %>% 
  filter(Phylum %in% c("Glomeromycota",  "Kickxellomycota")) %>%  
  mutate(depth_group = factor(depth_group, levels = c("reference", "0 - 60", "61 - 150", "152 - 610", "611 - 1372"))) %>% 
  ggplot(., aes(x = depth_group, y = val_sum, fill = depth_group)) +
  geom_boxplot()+
  geom_pwc(ref.group = "reference",
           label = "p.signif") +
  facet_wrap(~ Phylum, scales = "free") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "bottom") +
  scale_fill_viridis_d(name = "stockpile\ndepth (cm)",
                       begin = 0.1, end = 0.9) +
  labs(x = "Stockpile depth (cm)",
       y = "Relative abundance") +
  scale_y_continuous(labels = scales::percent)

ggsave(fig_SI8, file = "../plots/fig_SI8.jpg", height = 6, width = 12, units = "in", dpi = 800)

