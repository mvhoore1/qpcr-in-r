# library
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))

# two factors experimental design + balanced number of biol replicates and genes
qpcr = read.delim("qpcr-mock-data-balanced-two-factors.tsv", header = TRUE,stringsAsFactors = F)

# first, let's make a unique identifier by concatenating the xp_factor with the biol_rep columns
qpcr$unique_id = paste(qpcr$xp_factor_1,
                       qpcr$xp_factor_2,
                       qpcr$biol_rep,
                       sep = "_")

# Transform ct_vals by doing 2^-X
# this "2" should be able to change at some point because different primers will have different effciencies. However this is only important when you compare different primers/genes with each other in absolute terms
qpcr = mutate(qpcr,exp_ct_vals=2^-ct_vals)

#######################################################
# Finding possible outliers among biological replicates 
#######################################################

# calculate the spread among technical replicates within one sample

# if one technical replicate is far away from the others, there's an issue with that sample
# calculate the standard deviation based on all technical replicates
# if one is too far away (e.g. +-2SD or +-3SD) then it should be removed
# remove that technical replicate 

# plot the variation among technical replicates per sample
qpcr %>% 
  ggplot(aes(unique_id,exp_ct_vals,fill=xp_factor_1)) + 
  geom_boxplot() +
  geom_point() +
  facet_wrap(~ gene,scales = "free") +
  scale_color_brewer(palette = "Set2") +
  labs(x="Biological replicate",y="Expression (2^-Ct)",title = "Variation among technical replicates") +
  theme(axis.text.x = element_text(angle=90))

##########################################################
# average the technical replicates of the reference genes
#########################################################
reference_values = qpcr %>% 
  filter(gene_type == "reference") %>% 
  group_by(xp_factor_1,xp_factor_2,biol_rep,gene) %>%
  summarise(avg_reference = mean(exp_ct_vals))

# assigning a unique identifier for left_join
# we add the sample unique identifier to be able to combine the reference gene Ct values with the target genes Ct values
# making sure that reference genes Ct values are next to 
reference_values$unique_id = paste(
  reference_values$xp_factor_1,
  reference_values$xp_factor_2,
  reference_values$biol_rep,
  sep="_") 

# simplifying the reference_values dataframe
reference_values = reference_values[,c("avg_reference","unique_id")]

#########################################
# Left join to add the reference Ct value
#########################################
# adding reference gene 2^-Ct values to every 2^-Ct value for every sample (identified by xp-factor & biological rep)

qpcr.with.ref.values = dplyr::left_join(x = qpcr,y = reference_values,by="unique_id")

qpcr.with.ref.values = mutate(qpcr.with.ref.values,delta.ct=exp_ct_vals / avg_reference)
head(qpcr.with.ref.values)

################################################
# Average the technical replicates for all genes
################################################
qpcr.with.ref.values.avg = qpcr.with.ref.values %>% 
  group_by(xp_factor_1,xp_factor_2,biol_rep,gene) %>%
  summarise(avg_tech = mean(delta.ct))

# plot per xp_factor and facet per gene
qpcr.with.ref.values.avg %>% 
  ggplot(aes(x=xp_factor_1,y=avg_tech,fill=xp_factor_2)) +
  geom_boxplot() +
  # geom_point() +
  scale_y_log10() +
  facet_wrap(~ gene) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Experimental condition",y="Normalised expression (AU)") +
  theme(axis.text= element_text(color="black"))


