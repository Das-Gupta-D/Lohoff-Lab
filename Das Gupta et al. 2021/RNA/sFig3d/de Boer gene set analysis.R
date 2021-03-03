# prerequisits 
library(tidyverse)
library(readxl)
library(pheatmap)
library(fuzzyjoin)

# create funtcion to negate %in%
`%notin%` <- Negate(`%in%`)

#load data
Gene_Set <- read_excel("data/Gene Set de Boer_refined.xlsx")
TUMvBM <- read_excel("data/GSEA_scenarioTumvBM.xlsx")
NEJM <- read_excel("data/NEJM Genes Ph like.xlsx")
Gene_Set_human_refined <- read_excel("data/Gene Set de Boer human_refined.xlsx")

# create table of logFC of gene set genes for "TUMvBM" 
TUMvBM_diff <- 
  TUMvBM %>%
  pivot_longer(c(BM1:TUM3), names_to = "sample") %>%
  mutate(
    state = if_else(sample %in% c("BM1", "BM2"), "contr", "tumour")
  ) %>%
  filter(
    NAME %in% as_vector(Gene_Set)
  ) %>%
  group_by(NAME, state) %>%
  summarise(mean = mean(value)) %>%
  ungroup() %>%
  group_by(NAME) %>%
  mutate(
    entry = c(1:2)
  ) %>%
  ungroup() %>%
  mutate(
    foldchange = ifelse(entry==2, mean/lag(mean) , 0)
  ) %>%
  filter(
    entry == 2
  ) %>%
  mutate(
    log_FC = log10(foldchange)
  ) %>%
  select(NAME, log_FC) %>%
  arrange(log_FC) %>%
  rename(GENE = NAME)

# create table of logFC of gene set genes for "NEJM" 
NEJM_diff <- 
  NEJM %>%
  filter(
    GENE %in% as_vector(Gene_Set_human_refined)
  ) %>%
  arrange(`Ph-like_v_others_logFC`)

Gene_Set_human_refined_2 <- 
  Gene_Set_human_refined %>%
  rename(GENE = NAME)

refine_hum <- 
    anti_join(Gene_Set_human_2, NEJM_diff)

# combining the two tables (TUMvBM_diff and NEJM_diff)
pooled <- regex_inner_join(TUMvBM_diff, NEJM_diff, by="GENE", ignore_case = TRUE) %>%
  select(-GENE.x) %>%
  pivot_longer(c("log_FC", "Ph-like_v_others_logFC"), names_to = "set") %>%
  filter(
    GENE.y %notin% c("BIRC7", "TCFL5")  
  ) %>%
  group_by() %>%
  arrange(value)

# plotting the gene regulation patterns
pooled %>%
  ggplot(aes(x = reorder(GENE.y, -value), y = value, fill = set)) +
  geom_col(position = "dodge") +
  labs(
    title = "Comparison of C.Mulligan Ph-like genes to Irf4-ko leukemia",
    subtitle = "Geneset from de Boer"
  ) +
  ylab("log10(fold change)") +
  xlab(element_blank()) +
  scale_fill_discrete(labels = c("IRF4-ko Leukemia vs BM", "NEJM Ph-like ALL vs. Other ALL"))
