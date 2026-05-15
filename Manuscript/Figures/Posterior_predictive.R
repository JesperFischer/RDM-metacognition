pacman:p_load(cmdstanr, posterior)

###################################
# Read RDS
####################################

fit = readRDS("C:\\Users\\User\\OneDrive\\Documenten\\psychologie 1e Master\\internship\\Hierarchical_nolapse_50_clamped_standard.rds")

fit$summary()

