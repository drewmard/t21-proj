library(data.table)
library(scales)
library(ggplot2)
df1.full.mg = fread("~/Documents/Research/t21-proj/out/full/DEG_list/liver_v_femur.txt",data.table = F,stringsAsFactors = F)
genes_of_interest <- c("TOP2A", "CDC20",  "BIRC5", "CXCR4" , "IL6R",
                       "GATA1","TAL1","FOXO3","KLF1","AHSP",
                       "ITGA2B","SERPINE2","FOXO3","APOE","DUSP1")
df1.full.mg$geneLabel <- NA
df1.full.mg$geneLabel[df1.full.mg$names %in% genes_of_interest] <- df1.full.mg$names[df1.full.mg$names %in% genes_of_interest] 

tmp1 <- df1.full.mg[,c("logFC.h","P.Value.h","class",'names')]
tmp2 <- df1.full.mg[,c("logFC.t21","P.Value.t21","class",'names')]
colnames(tmp1) <- c("LFC","P","class",'names')
colnames(tmp2) <- c("LFC","P","class",'names')
tmp1$P <- -log10(tmp1$P)
tmp2$P <- log10(tmp2$P)
tmp <- rbind(tmp1,tmp2)
rng <- seq(floor(min(tmp$P)),ceiling(max(tmp$P)),1)
rng <- rng[rng%%2 == 0]

tmp$geneLabel <- NA
tmp$geneLabel[tmp$names %in% genes_of_interest] <- tmp$names[tmp$names %in% genes_of_interest]

g <- ggplot(tmp,aes(x=LFC,y=P,col=class)) + 
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid = element_blank()) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs(x="LFC (Healthy)",y="LFC (T21)") +
  # scale_color_brewer(palette="Set3",
  labs(col="",x="LFC",y=expression(-log[10](italic(P)))) +
  scale_y_continuous(breaks=rng,labels=abs(rng)) +
  scale_x_continuous(breaks=pretty_breaks())

pdf(paste0("~/Documents/Research/t21-proj/out/full/figures/volcano.liver_v_femur.pdf"),width=8,height=5)
gout <- g +
  geom_point(aes(alpha=ifelse(class=="none",0.5,1))) +
  scale_color_manual(
    labels=c("Environment-driven","Not significant","T21-induced","T21-reverted","T21-induced or\nenvironment-driven"),
    values=c("#8DD3C7" ,"#FFFFB3", "#FB8072" ,"#BEBADA","#80B1D3")) +
  guides(alpha="none")
print(gout)
dev.off()

pdf(paste0("~/Documents/Research/t21-proj/out/full/figures/volcano.liver_v_femur.labels.pdf"),width=8,height=5)
gout <- g +
  geom_point(aes(alpha=ifelse(class=="none",0.5,1))) +
  scale_color_manual(
    labels=c("Environment-driven","Not significant","T21-induced","T21-reverted","T21-induced or\nenvironment-driven"),
    values=c("#8DD3C7" ,"#FFFFB3", "#FB8072" ,"#BEBADA","#80B1D3")) +
  guides(alpha="none") +
  geom_text(aes(label=geneLabel),vjust=2)
print(gout)
dev.off()


j=2
system("mkdir -p ~/Documents/Research/t21-proj/out/full/figures/single_gene")
# for (j in 1:2) {
for (j in 1:length(genes_of_interest)) {
  # gout <- g + 
  gout <- ggplot(tmp[order(as.numeric(tmp$names==genes_of_interest[j])),],aes(x=LFC,y=P,col=class)) + 
    theme_bw() + 
    theme(panel.border = element_blank(),
          panel.grid = element_blank()) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    labs(x="LFC (Healthy)",y="LFC (T21)") +
    # scale_color_brewer(palette="Set3",
    labs(col="",x="LFC",y=expression(-log[10](italic(P)))) +
    scale_y_continuous(breaks=rng,labels=abs(rng)) +
    scale_x_continuous(breaks=pretty_breaks()) +
    geom_point(
    aes(
      alpha=ifelse(names!=genes_of_interest[j],0,1),
      col=(names!=genes_of_interest[j])
    )
  ) +
    scale_color_manual(values=c("black","orange")) +
    guides(alpha="none") +
    labs(title=genes_of_interest[j]) #+
    # geom_text(label=tmp$geneLabel,vjust=2,col='black')
  pdf(paste0("~/Documents/Research/t21-proj/out/full/figures/single_gene/volcano.liver_v_femur.",genes_of_interest[j],".pdf"),width=8,height=5)
  print(gout)
  dev.off()
}

# +




# g <- ggplot(df1.full.mg,aes(x=logFC.h,y=-log10(P.Value.h),col=class)) + 
#   geom_point(aes(alpha=ifelse(class=="none",0.5,1))) + 
#   theme_bw() + 
#   theme(panel.border = element_blank(),
#         panel.grid = element_blank()) +
#   geom_vline(xintercept = 0) +
#   geom_hline(yintercept = 0) +
#   labs(x="LFC (Healthy)",y="LFC (T21)") +
#   # scale_color_brewer(palette="Set3",
#   scale_color_manual(
#     labels=c("Environment-driven","Not significant","T21-induced","T21-reverted","T21-induced or\nenvironment-driven"),
#     values=c("#8DD3C7" ,"#FFFFB3", "#FB8072" ,"#BEBADA","#80B1D3")) +
#   guides(alpha="none") +
#   labs(col="") +
#   scale_x_continuous(breaks=pretty_breaks()) +
#   scale_y_continuous(breaks=pretty_breaks())
# g2 <- ggplot(df1.full.mg,aes(x=logFC.t21,y=-log10(P.Value.t21),col=class)) + 
#   geom_point(aes(alpha=ifelse(class=="none",0.5,1))) + 
#   theme_bw() + 
#   theme(panel.border = element_blank(),
#         panel.grid = element_blank()) +
#   geom_vline(xintercept = 0) +
#   geom_hline(yintercept = 0) +
#   labs(x="LFC (Healthy)",y="LFC (T21)") +
#   # scale_color_brewer(palette="Set3",
#   scale_color_manual(
#     labels=c("Environment-driven","Not significant","T21-induced","T21-reverted","T21-induced or\nenvironment-driven"),
#     values=c("#8DD3C7" ,"#FFFFB3", "#FB8072" ,"#BEBADA","#80B1D3")) +
#   guides(alpha="none") +
#   labs(col="") +
#   scale_x_continuous(breaks=pretty_breaks()) +
#   scale_y_continuous(breaks=pretty_breaks())
# library(cowplot)
# plot_grid(g,g2,nrow=2)
