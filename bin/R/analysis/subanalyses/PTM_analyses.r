GRAPHS <- list()
TABLES <- list()
#' PTM analyses
has_mods <- run$first %>% filter(!is.na(Mods))
UNIQUE_MODS <- run$first$Mods %>%
  discard(is.na) %>%
  map(\(x) str_split_1(x, "\\|")) %>%
  unlist() %>%
  map_chr(\(x) gsub(" [1-9]+", "", x)) %>%
  unique() %>%
  discard(\(x) grepl("0$", x))

ptms_first <- modMetrics(run$first)
ptms_sec <- modMetrics(run$sec)
ptm_percent_diff <- abs(ptms_sec$percentages - ptms_first$percentages)


# Try to identify qualitatively which PTMs tend to occur together
# Biased towards more commonly found PTMs?
dist <- ptms_first[[1]] |>
  t() |>
  vegan::vegdist(method = "bray")
pcoa <- vegan::wcmdscale(dist, eig = TRUE)
plot <- pcoa$points %>% ggplot(aes(x = Dim1, y = Dim2, color = rownames(.))) +
  geom_point(size = 4) +
  xlab("V1") +
  ylab("V2") +
  labs(title = "PCOA of PTMs (Bray-Curtis distance)", color = "Modification")
GRAPHS$ptm_pcoa <- plot

# Test if modifications are associated using chi square

# Transform tb so that proteins are in rows and mods in cols
# A protein has TRUE in the mod col if it has that mod
tb <- run$first
chi <- tb %>%
  filter(category != "other") %>%
  select(ProteinId, category) %>%
  left_join(., df2Tb(ptms_first$count_df, "ProteinId")) %>%
  mutate(across(is.double, \(x) ifelse(is.na(x), FALSE, TRUE)))

no_mods <- chi %>%
  select(-c(ProteinId, category)) %>%
  apply(1, \(x) ifelse(any(x), FALSE, TRUE)) %>%
  unlist()
chi$none <- no_mods

# Create 2 x 2 contingency table with mods as columns and categories as rows
# But since mods are not mutually exclusive, will need to test each
# modification individually against the categories
# TODO: Need to format this nicely
chosen <- "Met_Oxidation"
tab <- table(chi$category, chi[[chosen]])
ptm_tests <- list()
for (mod in UNIQUE_MODS) {
  ptm_tests[[mod]]$test <- chisq.test(chi$category, chi[[mod]])
  ptm_tests[[mod]]$table <- table(chi$category, chi[[mod]])
  print(ptm_tests[[mod]])
}
