#### calculate-K2P.R ####
## Kevin Labrador
## 2023-12-29

#### INTRODUCTION ####
# This script uses R to calculate the K2P genetic distances across taxonomic hierarchies - species, genus, and family. This was done to go over the limitations of BOLD systems workbench on classifying families:

# 1. Scaridae and Labridae are treated as separate families when in fact, Scaridae should nested within Labridae as the subfamily, Scarini.
# 2. BOLD assigns all groupers and allied species as Serranidae. However, recent taxonomic developments divided Serranidae into four families: (1) Serranidae, (2) Epinephelidae, (3) Anthiadidae, and (4) Liopropomatidae.

#### INITIALIZE ####

#### Housekeeping ####

# Set-up the working directory in the source file location: 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

wd <- getwd() # Do this after the working directory was set to source file location

# Assign directories to objects in the global environment
indir <- paste0("../data/")
outdir <- paste0("../results/")
scripts <- paste0("../scripts/")

#### SCRIPTS ####

# Coerce upper triangle and diagonal to NA
coerce_upper_triangle_to_na <- function(mat) {
  mat[upper.tri(mat, diag=T)] <- NA
  return(mat)
}

# Calculate pairwise K80 distances within a taxonomic hierarchy.
calculate_pairwise_k80 <- function (seq_list, model, level) {
  
  dist.df <- map (seq_list,
                  
                  # Calculate K80 distances per group
                  ~ dist.dna(.x, 
                             model = model,
                             pairwise.deletion = T,
                             as.matrix = T) %>% 
                    # Coerce upper triangle and diagonal to NA
                    coerce_upper_triangle_to_na() %>% 
                    
                    # Convert to data frame
                    as.data.frame() %>% 
                    
                    # Assign rownames to column
                    rownames_to_column("taxa_1") %>%
                    
                    # Pivot the data
                    pivot_longer(cols = contains("MINDA"), 
                                 names_to = "taxa_2",
                                 values_to = "distance") %>% 
                    
                    # Rearrange columns
                    select (taxa_1, taxa_2, distance) %>% 
                    
                    # Remove NAs
                    na.omit() %>% 
                    
                    # multiply distance by 100
                    mutate (distance = distance * 100)
                  
  ) %>% plyr::ldply(., .id = "taxa") %>% 
    mutate (level = level)
  
  return (dist.df)
}


#### PAIRWISE COMPARISONS ####
# Isolate the observations needed for pairwise comparison. Use the output from BOLD to identify these observations based on their process_id.
df <- data.bkp %>% 
  mutate (label = paste (scientific_name, process_id, sep="|"))

# WITHIN SPECIES
seq_list_species <- df %>% 
  group_by (scientific_name) %>% 
  
  # exclude species with only one sequence
  filter (n() > 1) %>% 
  split (., .$scientific_name) %>% 
  map (~.x %>% as_DNAbin(labels = label, sequences = sequence))

k80.species.df <- calculate_pairwise_k80(seq_list_species, model="K80", level = "species") %>% 
  separate (col = taxa_1, into = c("species_1", "process_id_1"), remove = T, sep = "\\|") %>% 
  separate (col = taxa_2, into = c("species_2", "process_id_2"), remove = T, sep ="\\|")


# WITHIN GENUS
seq_list_genus <- df %>% 
  
  # exclude genus with only one representative  
  group_by (genus) %>% 
  filter (n() > 1) %>% 
  split (., .$genus) %>% 
  map (~.x %>% as_DNAbin(labels = label, sequences = sequence))

k80.genus.df <- calculate_pairwise_k80(seq_list_genus, model = "K80", level = "genus") %>% 
  
  # Remove intraspecific comparisons
  separate (col = taxa_1, into = c("species_1", "process_id_1"), remove = T, sep = "\\|") %>%   
  separate (col = taxa_2, into = c("species_2", "process_id_2"), remove = T, sep ="\\|") %>% 
  filter (species_1 != species_2)

                                             
# WITHIN FAMILY
seq_list_family <- df %>% 
  group_by (family) %>% 
  
  # exclude families with only one genus
  filter (n() > 1) %>% 
  split (., .$family) %>% 
  map (~.x %>% as_DNAbin(labels = label, sequences = sequence))

k80.family.df <- calculate_pairwise_k80(seq_list_family, model = "K80", level = "family") %>% 
  
  # Remove intraspecific comparisons
  separate (col = taxa_1, into = c("species_1", "process_id_1"), remove = T, sep = "\\|") %>% 
  separate (col = taxa_2, into = c("species_2", "process_id_2"), remove = T, sep ="\\|") %>% 
  filter (species_1 != species_2) %>% 
  
  # Remove intra-genera comparisons
  separate (col = species_1, into = c("genus_1", NA), remove = F, sep = "\\ ") %>% 
  separate (col = species_2, into = c("genus_2", NA), remove = F, sep = "\\ ") %>%
  filter (genus_1 != genus_2) %>% 
  select (-c(genus_1, genus_2))

dist.df <- rbind (k80.species.df,
                  k80.genus.df,
                  k80.family.df) %>% 
  mutate (level = fct_relevel(level, "species", "genus", "family"))

write.csv(dist.df, paste0(outdir, "pairwise-k2p-distance.csv"))
