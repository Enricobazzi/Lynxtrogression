# load libraries
library(tidyverse)
library(RColorBrewer)
options(scipen=500)

# input data folder
input_folder <- "../../data/variant_filtering/rd_beds/"

# input suffix
suffix <- "_10kb_depth.bed"

# loop through populations
for (pop in c("lpa", "wel", "eel", "sel")){

 # list of chromosomes
 chrs <- system(paste0("ls ", input_folder, " | ", "grep ", pop, " | ",
                       "sed 's/", pop, "_//g'", " | ",
                       "sed 's/", suffix, "//g'", " | ",
                       "grep -v Super_Scaffold_10"), intern=T)
 
 # prepare bed
 pop_rd_bed <- data.frame()
 # start counter for color
 c=1
 
 # import each chromosome in a loop
 for (chr in chrs){
   chr_bed <- read_tsv(file=paste0(input_folder,pop,"_",chr,suffix), col_names = F,
                       progress = T, show_col_types = FALSE)
   
   if ((c %% 2) == 0){
     chr_bed$color <- "#BEBEBE"
   } else {
     chr_bed$color <- "#000000"
   }
   
   pop_rd_bed <- rbind(pop_rd_bed, chr_bed)
 
   c = c + 1
 }
 
 # add window number column
 pop_rd_bed$win_num <- c(1:NROW(pop_rd_bed))
 # rename columns
 colnames(pop_rd_bed) <- c("chromosome", "start", "end", "mean_rd", "color", "win_num")
 
 # calculate maximum and minimum depth based on mean and sd of mean rd values
 maxdepth <- mean(pop_rd_bed$mean_rd) + 0.5*sd(pop_rd_bed$mean_rd)

 # see how many windows are excluded based on max and min depth
 filtered_pop_rd_bed <- pop_rd_bed %>% filter(mean_rd < maxdepth)
 
 n_rows_filtered <- NROW(pop_rd_bed) - NROW(filtered_pop_rd_bed)
 print(paste(pop,":", "number of windows filtered:", n_rows_filtered))

 # plot
 rd_plot <- ggplot()+
   #geom_point(data=pop_rd_bed, aes(x=win_num, y=mean_rd), alpha=0.5, color=pop_rd_bed$color, size=4, shape=15)+
   geom_histogram(data=pop_rd_bed, aes(x = mean_rd), bins = 200)+
   scale_x_continuous(breaks = 0:4000*50, limits = c(0, mean(pop_rd_bed$mean_rd)*2)) +
   geom_vline(xintercept=maxdepth, linetype="dashed", colour="red", size=0.3) +
   theme_classic()
 
 # save plot
 ggsave(filename = paste0("../../plots/variant_filtering/rd_filters/", pop, "_rd_limits.pdf"),
        plot = rd_plot, 
        width = 4,
        height = 4,
        units = "in")
 
 # save bed of windows to be filtered
 bad_pop_rd_bed <- pop_rd_bed %>% filter(mean_rd > maxdepth)
 write.table(x = bad_pop_rd_bed,
             file = paste0("../../data/variant_filtering/", pop, "_rd_filter.bed"),
             quote=FALSE,  col.names = F, row.names = FALSE, sep= "\t")
 

}
