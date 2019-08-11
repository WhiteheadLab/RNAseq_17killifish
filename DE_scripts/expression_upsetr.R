library(UpSetR)

dim(sig_clade_physiology_interaction)
dim(sig_salinity_physiology_interaction)
dim(sig_salinity_clade_interaction)
dim(sig_threeway)
dim(sig_main_clade)
dim(sig_main_physiology)
dim(sig_main_salinity)

# Dataset
input <- c(
  Clade = length(sig_main_clade),
  Physiology = length(sig_main_physiology),
  Salinity = length(sig_main_salinity),
  `Physiology&Clade` = length(sig_clade_physiology_interaction),
  `Physiology&Salinity` = length(sig_salinity_physiology_interaction),
  `Salinity&Clade` = length(sig_salinity_clade_interaction),
  `Salinity&Clade&Physiology` = length(sig_threeway)
  )

# Plot
upset(fromExpression(input), order.by = "freq", decreasing = T)
