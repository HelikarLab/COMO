library(FactoMineR)
library(factoextra)



res.mca <- MCA(df, graph=FALSE, axes = c(1,2)) # mca analysis
options(ggrepel.max.overlaps = Inf)
plotname <- paste(c(mca_fig_dir,"/MCA_", technique, "_", c, "_", ".png"),collapse="")
print(pal)
print(grp)
mplot <- fviz_mca_ind(res.mca, repel = TRUE, ggtheme = theme_minimal(), labelsize=2, col.ind=grp)
ggsave(plotname, mplot, scale=5)