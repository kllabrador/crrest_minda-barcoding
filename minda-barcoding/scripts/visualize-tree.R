#### visualize-tree.R ####
## Kevin Labrador
## 2023-12-29

#### INTRODUCTION ####
# This script uses R to visualize the phylogenetic tree generated from RAxML. The newick file (*.tre) is used as input data.
# Tree node labels were initially renamed outside of R using Notepad++. These can also be done inside of R.

#### INITIALIZE ####

#### Housekeeping ####

# Set-up the working directory in the source file location: 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

wd <- getwd() # Do this after the working directory was set to source file location

# Assign directories to objects in the global environment
indir <- paste0("../data/")
outdir <- paste0("../results/")
scripts <- paste0("../scripts/")

# Load dataset. Export tree as a newick file. Make sure to include branch lengths and bootstrap support values.
infile <- paste0(outdir, "tree/K80/RAxML_bipartitionsBranchLabels.20240115_phylo-aln_renamed.tre")


#### PREPARE TREE FILE ####
NJ <- read.tree(infile) 

# Root tree with outgroup
outgroup <- NJ$tip.label[grep("Neotrygon", NJ$tip.label)]
NJ.rooted <- root(NJ, outgroup)


#### PREPARE BASE TREE ####
# Draw the base tree
tree.base <- ggtree (NJ
                     , layout = "rectangular"
                     , ladderize = T
)

# Extract data from the tree
meta.df <- tree.base$data

# Extract bootstrap values
bootstrap <- meta.df %>% 
  filter (isTip == F) %>% 
  mutate (label = as.numeric(label)) %>% 
  filter (label > 0.5)


# Create species column for aesthetics
tipLabel <- meta.df %>% 
  filter (isTip == T) %>% 
  full_join(., seq.metadata) %>% 
  mutate (flag = case_when (grepl("flag", Final.ID) ~ "Y"
                            , TRUE ~ "N")
  )

t <- paste(tipLabel$SampleID_alt2
           , paste (substr(tipLabel$Initial.ID, 1,3), str_remove(tipLabel$Initial.ID, ".*\\ "), sep=".")
           , sep = "_"
)

tipLabel <- tipLabel %>% 
  mutate (tipLabel = t) %>% 
  mutate (tipLabel = case_when (grepl ("Neoceratodus", label) ~ "Outgroup_Neo.forsteri"
                                , TRUE ~ tipLabel)
  )


tree <- tree.base + 
  #geom_hilight(node = 553, fill ="steelblue", alpha = 0.30, extend = 0.025) +   
  #geom_text(data=bootstrap, aes(label=node), col='red', size = 2) +
  geom_point(data=bootstrap, col='red', size = 1) +
  geom_point(data=tipLabel, aes(x = x, fill=flag), pch=21, col="black", size = 2, alpha=0.75) +
  #geom_highlight(node = 302)+
  geom_tiplab2(data=tipLabel, aes(label=Family, col=Family), hjust=-0.050, size = 1.25, align=T) +  
  guides(color = "none"
         , fill = "none"
  )+
  #theme(legend.position = "bottom"
  #      , legend.text = element_text(size = 5)
  scale_fill_manual(values = c("steelblue", "coral")) 

tree

ggsave(tree
       , filename = paste0(subdir, "FD_COI_NJ3.jpg")
       , width = 7.5, height = 7.5, units = "in", dpi = 300)


```