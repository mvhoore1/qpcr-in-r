# In the tutorial they added the samples-names by hand, I made it now that it will just take the first row from the excel sheet and use this as for samples
qPCR = read.csv("qPCR.csv", header = TRUE)
samples <- as.character(qPCR[,1])
data <- qPCR[-c(1)]
res <- pcr_analyze(data,
                  group_var = samples,
                  reference_gene = 'SAND',
                  reference_group = 'WTC')
