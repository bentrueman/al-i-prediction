
# setup -------------------------------------------------------------------

library("tidyverse")
library("visdat")

# read/rearrange ----------------------------------------------------------

data_train <- read_csv("data-clean/data-al-i-train.csv")

data_long <- data_train %>%
  rename_if(is.numeric, \(x) str_replace_all(x, "(^(?!cens))", "value_\\1")) %>% 
  pivot_longer(
    -c(site:value_longitude),
    names_to = c(".value", "param"),
    names_pattern = "([^_]*?)_(.+)"
  ) 

# missings ----------------------------------------------------------------

data_train %>% 
  select(-starts_with("cens_")) %>% 
  vis_dat()

# proportion missing: 

hi_pmiss <- data_long %>% 
  group_by(param) %>% 
  summarize(pmiss = mean(is.na(value))) %>% 
  ungroup() %>% 
  arrange(desc(pmiss)) %>% 
  filter(pmiss > .5)

# proportion non-detect ---------------------------------------------------

hi_pcens <- data_long %>% 
  group_by(param) %>% 
  summarize(pcens = mean(cens == "left")) %>% 
  ungroup() %>% 
  arrange(desc(pcens)) %>% 
  filter(pcens > .5)

# correlations ------------------------------------------------------------

data_train %>% 
  select(ph_sonde:avg_colour_colour_units) %>% 
  pivot_longer(-ali_ug_l, names_to = "param") %>% 
  group_by(param) %>% 
  summarize(
    cor = cor(ali_ug_l, value, method = "pearson", use = "complete")
  ) %>% 
  ungroup() %>% 
  arrange(desc(abs(cor))) %>% 
  anti_join(hi_pcens, by = "param") %>% 
  anti_join(hi_pmiss, by = "param") %>% 
  filter(abs(cor) >= .1)
