#######
# Earth Hologenome Initiative taxonomy colour palette generation
# by Antton Alberdi (antton.alberdi@sund.ku.dk)
# version=1.0
# 17/06/2023
#######

## LOAD LIBRARIES
library(ape)
library(tidyverse)
library(grid)

# GET INPUT DATA FROM GTDB
options(timeout=1000)

tree_URL="https://data.gtdb.ecogenomic.org/releases/release214/214.0/bac120_r214.tree"
download.file(tree_URL, "bac120_r214.tree")
tree <- read.tree("bac120_r214.tree")

taxonomy_URL="https://data.gtdb.ecogenomic.org/releases/release214/214.0/bac120_taxonomy_r214.tsv"
download.file(taxonomy_URL, "bac120_taxonomy_r214.tsv")
taxonomy <- read.table("bac120_taxonomy_r214.tsv")

# RUN ANALYSIS

#Create one (first) representative genome per phylum
phylum_genomes <- taxonomy %>%
  rename(genome=V1) %>%
  separate(V2, into = c("division","phylum","class","order","family","genus","species"), sep = ";", remove = FALSE) %>%
  select(genome,phylum) %>%
  filter(genome %in% tree$tip.label) %>%
  distinct(phylum, .keep_all = TRUE)

#Prune the phylogenetic tree with one representative per phylum
phylum_tree <- keep.tip(tree,phylum_genomes[,1])

#Sort phylum table by tree topology
phylum_genomes_sorted <- phylum_genomes %>%
  arrange(match(genome, phylum_tree$tip.label))

#Replace genome names in tips with phylum names
phylum_tree$tip.label <- phylum_genomes_sorted[,2]

#Declare rainbow palette of 255 colors
cols <- c("#5B0A76", "#63098B", "#7007AB", "#7C07CA", "#8207DF", "#8007EA","#7807EE", "#6E07EE", "#6407EF", "#5807EF", "#4907EF", "#3607EF","#2208EE", "#1208EC", "#0B08E9", "#0808E4", "#0808E0", "#0808DD","#0808D9", "#0808D4", "#0808CF", "#0808C9", "#0808C4", "#0808BF","#0808B9", "#0808B3", "#0808AD", "#0808A6", "#08089F", "#080899","#080B93", "#08108B", "#081782", "#081E79", "#082672", "#082E6D","#08366A", "#083E68", "#084668", "#084E68", "#08566A", "#085D6B","#08626D", "#086470", "#086372", "#086275", "#08617A", "#086282","#08658B", "#086893", "#086C98", "#08719B", "#08769D", "#087B9E","#08809F", "#08859F", "#0888A0", "#088CA2", "#0891A3", "#0896A5","#089CA6", "#08A1A7", "#08A6A9", "#08ACAC", "#08B1B1", "#08B6B6","#08BBBB", "#08C0C0", "#08C5C5", "#08C8C8", "#08CCCC", "#08D1D1","#08D5D5", "#08D8D8", "#08DBDB", "#08DEDE", "#08E1E1", "#08E3E3","#08E6E6", "#08E9E9", "#08EBEB", "#08EDEC", "#08EEEA", "#08EEE5","#08EEDF", "#08EDD8", "#08EBD1", "#08EACA", "#08E8C0", "#08E7B2","#08E6A2", "#08E593", "#08E385", "#08E27B", "#08E074", "#08DD72","#08DB72", "#08D873", "#08D575", "#08D176", "#08CC76", "#08C775","#08C173", "#08BC70", "#08B86B", "#08B563", "#08B25B", "#08AF54","#08AC51", "#08A950", "#08A450", "#089F50", "#089950", "#08944F","#088F4F", "#088A4F", "#08864E", "#08844B", "#088646", "#08893F","#088C38", "#088C30", "#088B25", "#098A18", "#0E8B0E", "#188D09","#259208", "#2F9708", "#349C08", "#37A108", "#39A608", "#3CAC08","#41B108", "#46B608", "#4BBB08", "#50C008", "#55C408", "#59C708","#5FC908", "#67CC08", "#6ECE08", "#76D108", "#7ED308", "#86D608","#8ED908", "#97DB08", "#A0DE08", "#A9E108", "#B3E308", "#BBE508","#C3E608", "#CDE608", "#D7E608", "#DFE508", "#E2E308", "#E1E008","#E0DD08", "#E0D908", "#E0D308", "#E0CD08", "#E0C608", "#E0BF08","#E0B908", "#E0B408", "#E0AE08", "#E0A608", "#E09F08", "#E09708","#E08F08", "#E08708", "#E07F08", "#E07708", "#E06F08", "#E06708","#E05F08", "#E05908", "#E05408", "#E04E08", "#DF4608", "#DF3F08","#DF3708", "#DE2E08", "#DD2508", "#DB1C08", "#D81308", "#D50C08","#D10908", "#CC0808", "#C80808", "#C40808", "#BE0808", "#B60808","#AF0808", "#A70808", "#9F0808", "#970808", "#8F0808", "#870808","#7F0808", "#780808", "#730808", "#6E0808", "#680808", "#6C0D0D","#701313", "#721616", "#731A1A", "#761E1E", "#772221", "#792626","#7B2929", "#7D2D2D", "#7F3131", "#803435", "#823838", "#843C3C","#86403F", "#884343", "#8A4746", "#8C4B4A", "#8E4F4E", "#8F5252","#915555", "#93595A", "#955D5D", "#966161", "#986464", "#9A6868","#9C6B6C", "#9E706F", "#A07373", "#A27777", "#A47B7B", "#A57E7E","#A88282", "#A98686", "#AA8A89", "#AD8D8D", "#AF9191", "#B09595","#B29999", "#B59C9C", "#B69F9F", "#B7A3A3", "#B9A7A7", "#BCAAAB","#BEAEAE", "#BFB2B2", "#C1B6B6", "#C3B9B9", "#C5BDBD", "#C7C0C1","#C8C5C5", "#CAC8C9", "#CCCCCC")

#Subset colors to the number of phyla
set.seed(100)
cols_subset <- cols[sort(sample(length(cols),nrow(phylum_genomes_sorted),replace=FALSE))]

#Visualise color
#grid.raster(cols, interpolate = FALSE)
#grid.raster(cols_subset, interpolate = FALSE)

#Add colors to phylum table
phylum_genomes <- phylum_genomes %>%
  mutate(colors=cols_subset)

# GENERATE OUTPUT

#Phylum-color table
write.table(phylum_genomes[,-1],"ehi_phylum_colors.tsv",col.names=T,row.names=F,quote=F,sep="\t")

#Phylum tree
write.tree(phylum_tree,"phylum_tree.tree")

#Phylum tree images
pdf("phylum_tree.pdf",width=4,height=25)
plot(phylum_tree, tip.color = cols_subset,cex=0.75)
dev.off()

png("phylum_tree.png",width=600,height=1800)
plot(phylum_tree, tip.color = cols_subset,cex=0.75)
dev.off()
