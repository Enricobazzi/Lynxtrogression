# load libraries
library(tidyverse)
library(RColorBrewer)

# input file from pop_vcf_missing_gts.sh
filename <- "../../data/variant_fitlering/lp_ll_introgression_perpop_missing_gts.csv"

# load table
miss_table <- read_csv(filename, progress = T, show_col_types = FALSE)

# extract each population's counts of missing values and calculate missing proportion
# by dividing by number of missing samples (Var1) by total number of samples in population:
freq_miss_table <- data.frame()

for (pop in c("lpa", "wel", "eel", "sel")) {
  print(pop)  
  pop_miss_table <- data.frame(table(miss_table[,pop]))
  pop_miss_table$Var1 <- as.numeric(levels(pop_miss_table$Var1))/max(as.numeric(levels(pop_miss_table$Var1)))
  pop_miss_table$Freq_prop <- (pop_miss_table$Freq)/NROW(miss_table)
  pop_miss_table$Freq_cumsum <- cumsum(pop_miss_table$Freq_prop)
  pop_miss_table$population <- pop
  freq_miss_table <- rbind(freq_miss_table, pop_miss_table)
}

# extract each population into its own table for easier plot
lpa <- freq_miss_table %>% filter(population == "lpa")
wel <- freq_miss_table %>% filter(population == "wel")
eel <- freq_miss_table %>% filter(population == "eel")
sel <- freq_miss_table %>% filter(population == "sel")

# cumulative proportion of SNPs included vs decreasing missingness filter strictness
cummiss <- ggplot() +
  geom_line(data=lpa, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="blue", size=1.5)+
  geom_point(data=lpa, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="blue", size=4, shape=15)+
  geom_line(data=wel, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="red", size=1.5)+
  geom_point(data=wel, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="red", size=4, shape=15)+
  geom_line(data=eel, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="black", size=1.5)+
  geom_point(data=eel, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="black", size=4, shape=15)+
  geom_line(data=sel, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="green", size=1.5)+
  geom_point(data=sel, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="green", size=4, shape=15)+
  scale_x_continuous(limits = c(0, 0.4))+
  xlab("Proportion of missing data")+
  ylab("Proportion of SNPs included")+
  theme_minimal()

# write summary table
write.table(x = freq_miss_table,
            file = paste0("../../data/variant_fitlering/freq_miss_table.tsv"),
            quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")

# save plot
pdf(file = paste0("../../plots/variant_filtering/cumulative_miss.pdf"), width = 8, height = 4)
cummiss
dev.off()

# get bed file of snps to filter in each population and write it

lpa_bed <- miss_table[,c(1,2,3)] %>% filter(lpa >= 3)
prop_filtered <- NROW(lpa_bed)/NROW(miss_table)
print(paste("proportion filtered for LPA:", prop_filtered))
lpa_bed$start <- (lpa_bed$position)-1
lpa_bed <- lpa_bed[,c(1,4,2)]
write.table(x = lpa_bed,
            file = paste0("../../data/variant_fitlering/lpa_miss_filter.bed"),
            quote=FALSE,  col.names = F, row.names = FALSE, sep= "\t")


wel_bed <- miss_table[,c(1,2,4)] %>% filter(wel >= 3)
prop_filtered <- NROW(wel_bed)/NROW(miss_table)
print(paste("proportion filtered for WEL:", prop_filtered))
wel_bed$start <- (wel_bed$position)-1
wel_bed <- wel_bed[,c(1,4,2)]
write.table(x = wel_bed,
            file = paste0("../../data/variant_fitlering/wel_miss_filter.bed"),
            quote=FALSE,  col.names = F, row.names = FALSE, sep= "\t")

eel_bed <- miss_table[,c(1,2,5)] %>% filter(eel >= 3)
prop_filtered <- NROW(eel_bed)/NROW(miss_table)
print(paste("proportion filtered for EEL:", prop_filtered))
eel_bed$start <- (eel_bed$position)-1
eel_bed <- eel_bed[,c(1,4,2)]
write.table(x = eel_bed,
            file = paste0("../../data/variant_fitlering/eel_miss_filter.bed"),
            quote=FALSE,  col.names = F, row.names = FALSE, sep= "\t")

sel_bed <- miss_table[,c(1,2,6)] %>% filter(sel >= 2)
prop_filtered <- NROW(sel_bed)/NROW(miss_table)
print(paste("proportion filtered for SEL:", prop_filtered))
sel_bed$start <- (sel_bed$position)-1
sel_bed <- sel_bed[,c(1,4,2)]
write.table(x = sel_bed,
            file = paste0("../../data/variant_fitlering/sel_miss_filter.bed"),
            quote=FALSE,  col.names = F, row.names = FALSE, sep= "\t")
