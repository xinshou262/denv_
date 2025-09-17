setwd("D:\\DENV-3\\drug-degcount\\4genes\\ssGSEA免疫浸润")
####5.免疫浸润ssGSEA分析####
library(pheatmap)
# 读取表达矩阵数据
GSE39582_series <- read.table("TPM.txt", check.names = FALSE, header = TRUE, sep = "\t")

# Step 2: 去掉重复的基因名，只保留第一次出现的基因名
GSE39582_series <- GSE39582_series[!duplicated(GSE39582_series[, 1]), ]

# Step 3: 将第一列（基因名）设置为行名
rownames(GSE39582_series) <- GSE39582_series[, 1]

# 删除原来的基因名列（如果需要）
GSE39582_series <- GSE39582_series[, -1]
GSE39582_series=log2(GSE39582_series+1)
group<- read.table("group.txt", check.names = FALSE, header = TRUE,  sep = "\t")
# 创建一个新的因子列，以确保对照组和治疗组按顺序排列
# 假设分组信息是 con 和 treat
group$group <- factor(group$group, levels = c("con", "treat"), labels = c("Control", "DENV"))

rownames(group) <- group[, 1]
# 排序分组信息，根据分组顺序重新排列样本（列名）
sorted_samples <- group[order(group$group), ]$Sample  # 假设 group 里有一个 "Sample" 列表示样本名称

#GSE39582_series_tumor = GSE39582_series[,1:19]

library(tidyverse)
cellMarker1 <- read.delim("CellReports.txt", header = F, sep = "\t") #
geneSet<-split(cellMarker1$V1, cellMarker1$V2)


library(GSVA)
GSE39582_series_tumor=GSE39582_series
GSE39582_series_tumor=as.matrix(GSE39582_series_tumor)
# 使用 match() 对 GSE39582_series_tumor 进行重新排序
GSE39582_series_tumor <- GSE39582_series_tumor[, match(sorted_samples, colnames(GSE39582_series_tumor))]
gsva_matrix1 <- gsva(GSE39582_series_tumor, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
# Estimating ssGSEA scores for 28 gene sets.
# [1] "Calculating ranks..."
# [1] "Calculating absolute values from ranks..."
# |===========================================================================================| 100%
# 
# [1] "Normalizing..."
# 图1
res <- gsva_matrix1
#pheatmap(res, show_colnames = F)
# 确保样本顺序已按分组排列（Control在前，DENV在后）
 GSE39582_series_tumor <- GSE39582_series_tumor[, sorted_samples]

 # 创建列注释数据框（样本分组信息）
  annotation_col <- data.frame(Group = group[sorted_samples, "group"])  # 按排序后的样本提取分组
 rownames(annotation_col) <- sorted_samples  # 行名为样本名
 
   # 定义分组颜色（可选）
  annotation_colors <- list(
         Group = c(Control = "#1B9E77", DENV = "#D95F02")  # 使用颜色编码或名称
     )
  png("gsva_heatmap.png", width = 800, height = 600)  # 可以调整图像宽度和高度（像素）
   # 绘制热图并添加注释
   hotmap<-pheatmap(
        gsva_matrix1,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_colnames = FALSE,          # 不显示列名
         cluster_cols = FALSE,            # 不聚类列以保持原始顺序
         main = "GSVA Scores Heatmap",    # 添加标题
         color = colorRampPalette(c("#6238ff","#ffffff","#ff220e"))(100)  # 自定义颜色梯度
     )
   
   dev.off()
##加上样本名####
#group_list <- colnames(GSE39582_series_tumor)
#annotation <- data.frame(group_list)
#rownames(annotation) <- colnames(res)
#pheatmap(res,
         #show_colnames = F,
         #annotation_col = annotation,
         #cellwidth = 15, # 每个单元格的宽度
         #cellheight = 15, # 每个单元格的高度
         #fontsize = 8, # 字体大小
#)
#加上分组信息


# 图2
# 免疫细胞区分
# annotation1 <- data.frame(num1 = colnames(GSE39582_series_tumor),
#                           num2 = colnames(GSE39582_series_tumor)
# )
# library(pheatmap)
# gsva_matrix1<- t(scale(t(gsva_matrix)))
# gsva_matrix1[gsva_matrix1< -2] <- -2
# gsva_matrix1[gsva_matrix1>2] <- 2
# anti_tumor <- c('Activated CD4 T cell', 'Activated CD8 T cell', 'Central memory CD4 T cell', 'Central memory CD8 T cell', 'Effector memeory CD4 T cell', 'Effector memeory CD8 T cell', 'Type 1 T helper cell', 'Type 17 T helper cell', 'Activated dendritic cell', 'CD56bright natural killer cell', 'Natural killer cell', 'Natural killer T cell')
# pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'CD56dim natural killer cell', 'Immature dendritic cell', 'Macrophage', 'MDSC', 'Neutrophil', 'Plasmacytoid dendritic cell')
# anti<- gsub('^ ','',rownames(gsva_matrix1))%in%anti_tumor
# pro<- gsub('^ ','',rownames(gsva_matrix1))%in%pro_tumor
# non <- !(anti|pro)
# gsva_matrix1<- rbind(gsva_matrix1[anti,],gsva_matrix1[pro,],gsva_matrix1[non,])
# normalization<-function(x){
#   return((x-min(x))/(max(x)-min(x)))}
# nor_gsva_matrix1 <- normalization(gsva_matrix1)
# annotation_col = data.frame(patient=annotation1$num2)
# rownames(annotation_col)<-colnames(GSE39582_series_tumor)
# bk = unique(c(seq(0,1, length=100)))
# plot = pheatmap(
#   nor_gsva_matrix1, # 输入数据矩阵
#   show_colnames = F, # 是否显示列名
#   cluster_rows = F, # 是否对行进行聚类
#   cluster_cols = F, # 是否对列进行聚类
#   annotation_col = annotation_col, # 列注释信息
#   breaks = bk, # 颜色分段的断点
#   cellwidth = 20, # 每个单元格的宽度
#   cellheight = 20, # 每个单元格的高度
#   fontsize = 8, # 字体大小
#   gaps_row = c(12, 20) # 行之间的间隔大小
# )
# 
# #filename = 'ssgsea.pdf',width = 8)
# 
# # 图3
# # score加和后，ggplot2进行绘图
# anti_tumor <- c('Activated CD4 T cell', 'Activated CD8 T cell', 'Central memory CD4 T cell', 'Central memory CD8 T cell', 'Effector memeory CD4 T cell', 'Effector memeory CD8 T cell', 'Type 1 T helper cell', 'Type 17 T helper cell', 'Activated dendritic cell', 'CD56bright natural killer cell', 'Natural killer cell', 'Natural killer T cell')
# pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'CD56dim natural killer cell', 'Immature dendritic cell', 'Macrophage', 'MDSC', 'Neutrophil', 'Plasmacytoid dendritic cell')
# # load('score.Rdata')
# anti<- as.data.frame(gsva_matrix1[gsub('^ ','',rownames(gsva_matrix1))%in%anti_tumor,])
# pro<- as.data.frame(gsva_matrix1[gsub('^ ','',rownames(gsva_matrix1))%in%pro_tumor,])
# anti_n<- apply(anti,2,sum)
# pro_n<- apply(pro,2,sum)
# patient <- annotation1$num2[match(colnames(gsva_matrix1),annotation1$num1)]
# library(ggplot2)
# data <- data.frame(anti=anti_n,pro=pro_n,patient=patient)
# anti_pro<- cor.test(anti_n,pro_n,method='pearson')
# gg<- ggplot(data,aes(x = anti, y = pro),color=patient) + 
#   xlim(-20,15)+ylim(-15,10)+
#   labs(x="Anti-tumor immunity", y="Pro-tumor suppression") +
#   geom_point(aes(color=patient),size=3)+geom_smooth(method='lm')+
#   annotate("text", x = -5, y =7.5,label=paste0('R=',round(anti_pro$estimate,4),'\n','p<0.001'))
# # ggsave(gg,filename = 'cor.pdf', height = 6, width = 8)
# 
# ###箱线图####
# # 转换数据格式为适合ggplot的长格式
 library(tidyr)
library(ggpubr)
# 
# # 转置并转换为数据框（行名为样本，列名为基因集）
 a <- as.data.frame(t(gsva_matrix1))

# # 添加分组信息（确保样本顺序与分组信息一致）
a$group <- group[sorted_samples, "group"]
# 
# # 添加样本名（可选）
 a$sample <- rownames(a)
# 
# # 转换为长格式
 ggsea <- pivot_longer(
   a,
   cols = -c(group, sample),
   names_to = "GeneSet",
   values_to = "Expression"
 )
 
# 绘制箱线图
ggplot(ggsea, aes(x = GeneSet, y = Expression, fill = group)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    position = position_dodge(0.8)
  ) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  labs(y = "Expression", x = NULL, fill = "Group") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  ) +
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,
    label.y = max(ggsea$Expression) * 1.05  # 调整标签位置
  )

#绘制marker基因和细胞之间的相关性热图####

library(Hmisc)  # 确保已安装并加载Hmisc包
library(pheatmap)
library(dplyr)
library(grid)

# 定义目标基因（注意检查基因名大小写是否与表达矩阵一致）
mygene <- c("CXCL10", "EPHB2", "EZH2")

# 假设：
# gsva_matrix1 是GSVA结果矩阵（行：基因集，列：样本）
# GSE39582_series_tumor 是表达矩阵（行：基因，列：样本）

# Step 1: 确保样本顺序一致（重要！）
common_samples <- intersect(colnames(gsva_matrix1), colnames(GSE39582_series_tumor))
gsva_common <- gsva_matrix1[, common_samples]
exp_common <- GSE39582_series_tumor[, common_samples]

# Step 2: 提取目标基因表达数据（检查基因是否存在）
missing_genes <- setdiff(mygene, rownames(exp_common))
if (length(missing_genes) > 0) {
  warning("以下基因在表达矩阵中不存在: ", paste(missing_genes, collapse = ", "))
  mygene <- intersect(mygene, rownames(exp_common))
}

gene_exp <- exp_common[mygene, ]

# Step 3: 合并GSVA评分与基因表达数据（样本为列，合并后转置）
combined_data <- t(rbind(gsva_common, gene_exp))  # 转置后行是样本，列是变量

# Step 4: 计算相关性矩阵
cor_result <- rcorr(combined_data)
# 定义关键索引
num_genesets <- nrow(gsva_common)
num_genes <- length(mygene)
total_variables <- num_genesets + num_genes


# 提取相关性矩阵和P值
m <- cor_result$r[1:num_genesets, (num_genesets + 1):total_variables]
p <- cor_result$P[1:num_genesets, (num_genesets + 1):total_variables]

# Step 5: 创建显著性标记矩阵
tmp <- matrix(
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "***",
           ifelse(p < 0.01, "**",
                  ifelse(p < 0.05, "*", "")
           )
    )
  ),
  nrow = nrow(p),
  ncol = ncol(p)
)

## Step 6: 绘制热图####
p <- pheatmap(
  t(m),
  display_numbers = t(tmp),
  angle_col = 45,
  color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
  border_color = "grey60",
  cellwidth = 20,
  cellheight = 20,
  treeheight_col = 0,
  treeheight_row = 0,
  main = "",  # 隐藏默认标题
  silent = TRUE
)

# 绘制热图
grid.newpage()
grid.draw(p$gtable)

# 在右上方添加标题
pushViewport(viewport(
  x = unit(1.08, "npc"),  # 右侧位置（0.95表示距离左边95%画布宽度）
  y = unit(0.72, "npc"),  # 顶部位置
  width = unit(0.4, "npc"),  # 标题区域宽度
  just = c("right")   # 右对齐+顶部对齐
))

grid.text(
  "Correlation",
  gp = gpar(fontsize = 12, fontface = "bold"),
  hjust = 1  # 文本右对齐
)

popViewport()
#
grid.text(
      "*p<0.05",
     gp = gpar(fontsize = 9),  # 默认字体，不加粗
       hjust = 1,  # 文本右对齐
       x = unit(0.86, "npc"),  # 调整 x 位置（0.95表示靠近右边）
       y = unit(0.75, "npc")  # 调整 y 位置（垂直位置）
   )
grid.text(
  "**p<0.01",
  gp = gpar(fontsize = 9),  # 默认字体，不加粗
  hjust = 1,  # 文本右对齐
  x = unit(0.86, "npc"),  # 调整 x 位置（0.95表示靠近右边）
  y = unit(0.79, "npc")  # 调整 y 位置（垂直位置）
)
grid.text(
  "***p<0.001",
  gp = gpar(fontsize = 9),  # 默认字体，不加粗
  hjust = 1,  # 文本右对齐
  x = unit(0.86, "npc"),  # 调整 x 位置（0.95表示靠近右边）
  y = unit(0.83, "npc")  # 调整 y 位置（垂直位置）
)
#### 新增：免疫细胞相关性热图 ####
# 计算不同免疫细胞/基因集之间的相关系数矩阵
cor_matrix <- cor(t(gsva_matrix1))  # 转置后计算行（基因集）之间的相关性
pdf("免疫细胞相关性热图.pdf", width = 12, height = 8)  # 宽度和高度单位是英寸
# 绘制相关性热图
pheatmap(
  cor_matrix,
  color = colorRampPalette(c("#0033FF", "white", "#FF9900"))(100),
  border_color = "grey80",
  clustering_method = "complete",   # 聚类方法
  treeheight_row = 30,             # 行聚类树高度
  treeheight_col = 30,             # 列聚类树高度
  show_colnames = TRUE,            # 显示列名
  show_rownames = TRUE,            # 显示行名
  fontsize_row = 8,                # 行名字体大小
  fontsize_col = 8,                # 列名字体大小
  main = "",
  annotation_colors = if(exists("annotation_colors")) annotation_colors else list(), # 使用之前定义的配色
  cellwidth = 12,
  cellheight = 12
)
dev.off()
