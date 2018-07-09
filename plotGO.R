target_gene_id <- unique(read.delim("miRNA-gene interactions.txt")$EntrezID)
# BiocInstaller::biocLite("clusterProfiler")
# BiocInstaller::biocLite("org.Hs.eg.db")

display_number = c(15, 10, 15)
## GO enrichment with clusterProfiler
library(clusterProfiler)
ego_MF <- enrichGO(OrgDb="org.Hs.eg.db",
             gene = target_gene_id,
             pvalueCutoff = 0.05,
             ont = "MF",
             readable=TRUE)
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[1], ]
# ego_result_MF <- ego_result_MF[order(ego_result_MF$Count),]

ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id,
                   pvalueCutoff = 0.05,
                   ont = "CC",
                   readable=TRUE)
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
# ego_result_CC <- ego_result_CC[order(ego_result_CC$Count),]

ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id,
                   pvalueCutoff = 0.05,
                   ont = "BP",
                   readable=TRUE)
ego_result_BP <- na.omit(as.data.frame(ego_BP)[1:display_number[3], ])
# ego_result_BP <- ego_result_BP[order(ego_result_BP$Count),]

go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                                   Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                                   GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                                   type=factor(c(rep("biological process", display_number[1]), rep("cellular component", display_number[2]),
                                          rep("molecular function", display_number[3])), levels=c("molecular function", "cellular component", "biological process")))

## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                       collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}

labels=(sapply(
  levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)],
  shorten_names))
names(labels) = rev(1:nrow(go_enrich_df))

## colors for bar // green, blue, orange
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
library(ggplot2)
p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_manual(values = CPCOLS) + theme_bw() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")

p

pdf("go_enrichment_of_miRNA_targets.pdf")
p
dev.off()

svg("go_enrichment_of_miRNA_targets.svg")
p
dev.off()
