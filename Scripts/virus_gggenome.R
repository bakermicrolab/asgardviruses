#Purpose
#Create linear genome representations of virus genomes from Rambo et al.

library(tidyverse)
library("gggenomes")
library(RColorBrewer)

#Set this to the Data directory path in the GitHub repo
setwd("../Data")

virus_anno_df <- read.table("DataS5_GeneMap.txt", header = TRUE, sep = "\t", quote = "") %>% as_tibble()

virus_seqs <- read.csv("virus_seq_info.csv", header = TRUE) %>% as_tibble()

virus_anno_df$Proposed_Name <- gsub("`", "", virus_anno_df$Proposed_Name)

protospacer <- read.table("protospacer_coords.tsv", header = TRUE, sep = "\t", quote = "")

protospacer <- protospacer %>% rename(feat_id = query_id) %>% select(feat_id, seq_id, mismatch, start, end) %>%
  mutate(seq_id = gsub("_scaffold_2kb_scaffold_[0-9]+", "", seq_id),
         mismatch = as.factor(mismatch),
         feat_id = gsub("_____", "-", feat_id),
         feat_id = gsub("_[0-9]+_[0-9]+$", "", feat_id)) %>%
  filter(!grepl("Nidhogg", seq_id))


va_gene_df <- virus_anno_df %>% dplyr::mutate(Gene_Coordinates = gsub("\\(|\\)", "", Gene_Coordinates),
                                            Gene_Coordinates = gsub("\\.\\.", "-", Gene_Coordinates)) %>%
    tidyr::separate(Gene_Coordinates, c("start", "end"), sep = "-", remove = TRUE) %>%
    dplyr::mutate(start = as.numeric(start),
                  end = as.numeric(end)) %>%
    dplyr::select(Proposed_Name, Gene_ID,
           start, end, strand,
           gene_name) %>%
    dplyr::mutate(width = abs(start - end) + 1,
           type = "CDS",
           geom_id = Proposed_Name) %>%
    rename(seq_id = Proposed_Name,
           feat_id = Gene_ID,
           Note = gene_name) %>%
  dplyr::distinct()


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]

qual_col_vec <- unlist(mapply(brewer.pal,
                              qual_col_pals$maxcolors, rownames(qual_col_pals)))

#Create a dataframe of gene colors
gene_note_uniq <- unique(va_gene_df$Note)

note_color <- data.frame(gene_note_uniq, qual_col_vec[1:length(gene_note_uniq)])

colnames(note_color) <- c("Note", "color")

note_color$color[which(note_color$Note == "hypothetical")] <- toupper("#d6d6d4")

va_gene_df <- va_gene_df

#Remove hypothetical annotations to avoid clutter
va_gene_df$Note[which(grepl("[Hh]ypo", va_gene_df$Note))] <- ""

col <- as.character(note_color$color)
names(col) <- as.character(note_color$Note)
names(col)[2] <- ""

gggenomes(va_gene_df %>% filter(!grepl("Nidhogg", seq_id)), virus_seqs %>% filter(!grepl("Nidhogg", seq_id))) %>%
  add_feats(protospacer=protospacer) +
  geom_seq() +
  geom_bin_label() +
  geom_feat(data=feats(protospacer), aes(color = mismatch), position = "identity", size = 12) +
  geom_gene(aes(fill=Note), show.legend = FALSE) +
  geom_gene_tag(aes(label=Note), nudge_y=0.08,
                check_overlap = FALSE,
                size = 2) +
  scale_fill_manual(values = col)
