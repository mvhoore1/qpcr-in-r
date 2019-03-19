# a script to generate a mock qPCR file

# a one-factor qPCR design file will have 7 columns, the last one being the Ct values

# make each column individually
wells = seq(1:96)
sample_ids = paste("sample",seq(1:96),sep = "")
xp_factor = c(
  rep("wild-type",48),
  rep("mutant",48))
biological_replicates = rep(c(1,1,2,2,3,3),times = 96/6)
genes = rep(c("ACT7","DREB","SAND"),times = 96/6)
gene_type = c(
  rep('reference',times=6),
  rep('target',times=96-6)
)
ct_values = rnorm(n = 96,sd = 10,mean = 25)


# making the dataframe
qpcr_mock_data = data.frame(
  wells = wells,
  sampleids = sample_ids,
  xp_factor = xp_factor,
  biol_rep = biological_replicates,
  gene = genes,
  gene_type = gene_type,
  ct_vals = ct_values
)

write.table(qpcr_mock_data,file = "qpcr-mock-data.tsv",quote = F,sep = "\t",row.names = F)