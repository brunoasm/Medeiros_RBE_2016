## This script is used to generate tree for the 
## beetle species that were collected here.
## Species identified only to family/subfamily
## I will use a modified version of phytools function add.species.to.genus
## If more than one species in one genus, will




library(ape)
library(phytools)
tree <- read.nexus(file = "FinalBEAST_forTreeEdit.tre")
#remove outgroups
tree <- drop.tip(tree,tree$tip.label[grepl('OG|Neu',tree$tip.label)])
#read my species data table
data_table <- read.csv("../all_data.csv")
species <- apply(unique(data_table[c('higher_taxon','species')]),1,function(x){paste(x, collapse='_')})
#read table with tip taxonomy
taxonomy <- read.csv("taxa_molecular_phylogeny.csv")


translate_names <- function(x){
  if (grepl('Bruchinae$', x)){
    return('Bruchinae')
  }
  else if (taxonomy$Subfamily[taxonomy$name_in_phylogeny == x] == 'Scolytinae'){
    return('Scolytinae')
  }
  else if (taxonomy$Subfamily[taxonomy$name_in_phylogeny == x] == 'Scydmaeninae'){
    return('Scydmaeninae')
  }
  else if (taxonomy$Subfamily[taxonomy$name_in_phylogeny == x] == 'Sphaeridiinae'){
    return('Sphaeridiinae')
  }
  else{
    return(as.character(taxonomy$Family[taxonomy$name_in_phylogeny == x]))
  }
}


tree$tip.label = paste(sapply(tree$tip.label,translate_names),paste('sp',1:length(tree$tip.label), sep = ''), sep = '_')

tree_list = list()
class(tree_list) = "multiPhylo"

for (j in 1:100){
  cat (paste('TREE',j))
  tree_list[[j]] = tree #start with backbone tree
  for (i in sample(1:length(species),size = length(species),replace = F)){ #add species randomly
    tree_list[[j]] = add.species.to.genus(tree_list[[j]],species[i],where='random')
  }
  tree_list[[j]] = drop.tip(tree_list[[j]],grep('_sp',tree_list[[j]]$tip.label)) #remove species not in this study
  write.tree(tree_list, 'final_100_trees.tre')
}
