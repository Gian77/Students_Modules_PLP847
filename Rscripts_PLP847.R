# BASIC R COMMUNTIY ANALYSIS ---------
# Michigan State University
# written by Gian MN Benucci
# Dec 20, 2018

# load packages ---------------------
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(indicspecies)
library(vegan)

# create a phyloseq object ----------
otus_ITS <- read.delim("otu_table_ITS_UPARSE.txt",row.names=1) 
otus_phy_ITS <-otu_table(otus_ITS,taxa_are_rows = TRUE)
metadata_ITS <-read.delim("mapping_ITS_new.txt",row.names=1)
metadata_phy_ITS <-sample_data(metadata_ITS)
taxonomy_ITS <-read.delim("consensus_taxonomy_ITS.txt", header=TRUE, row.names=1)
taxonomy_phy_ITS <- tax_table(as.matrix(taxonomy_ITS))
otus_seq_ITS <- readDNAStringSet("otus_ITS.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

physeq_obj_ITS <- phyloseq(otus_phy_ITS,metadata_phy_ITS,taxonomy_phy_ITS,otus_seq_ITS)
tax_table(physeq_obj_ITS)[tax_table(physeq_obj_ITS)==""]<- "NA"


physeq_obj_ITS
head(tax_table(physeq_obj_ITS))
head(otu_table(physeq_obj_ITS))
head(sample_data(physeq_obj_ITS))


# filtering dataset --------------------
tax_table(physeq_obj_ITS)[tax_table(physeq_obj_ITS)==""]<- "Unclassified"

any(tax_table(physeq_object) == "Unclassified")
any(tax_table(physeq_object) == "Protista")

physeq_object <- subset_taxa(physeq_object, Kingdom == "Fungi")
physeq_object

tax_table(physeq_object)[tax_table(physeq_object)==""]<- "Unclassified"
tax_table(physeq_object)


# >>> Filtering out OTUs ----------------------------------
# Oliver et al. 2015, PCR errors and tag switching
# Lindhal et al. 2013, tag switching - that's a good  one!
# Barberan et al. 2012, removing OTUs that appear in less than x samples

physeq_object -> physeq_object_filt
otu_table(physeq_object_filt)[otu_table(physeq_object_filt) <= 4] <- 0 ### tag switching
otu_table(physeq_object_filt) <- otu_table(physeq_object_filt)[which(rowSums(otu_table(physeq_object_filt)) >= 10),] ### PCR Errors 
physeq_object_filt

# OTUs that are found in at least 5% of samples 
library(metagMisc)
physeq_object_filt <- phyloseq_filter_prevalence(physeq_object_filt, prev.trh = 0.05, abund.trh = NULL)
physeq_object_filt

...
