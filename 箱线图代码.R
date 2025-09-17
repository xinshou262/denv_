rm(list = ls())
setwd("D:\\DENV-3\\drug-degcount\\4genes\\验证组TPM数据\\ROC曲线")
library(ggplot2)
library(ggpubr)
library(reshape2)
expFile = "完整箱线图数据.txt"
# 读取数据
data <- read.table(expFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE,fileEncoding = "UTF-8")

data$group <- factor(data$group, levels = c("con", "treat"), labels = c("Control", "DENV"))

# 循环绘制每个基因的箱线图并保存,从第三列开始是基因表达量
for (gene in colnames(data)[3:ncol(data)]) {
  
  # 筛选出当前基因的数据
  gene_data <- data[, c("group", gene)]  # 提取group和当前基因的列
  # 删除含有NA的行
  gene_data <- na.omit(gene_data)
  
  
  # 创建箱线图，设置填充色（fill）和边框颜色（color）
  p <- ggboxplot(gene_data, x = "group", y = gene,
                 fill = "group",       # 设置箱体填充色
                 color = "black",      # 设置边框颜色
                 palette = c("#FFD699", "#CCCCFF"))  # 自定义填充色
  
  # 添加比较（p值）
  com <- list(c("Control", "DENV"))
  p <- p + stat_compare_means(comparisons = com, aes(label = ..p.format..))
  
  # 加粗字体
  p <- p + theme(
    axis.title = element_text(face = "bold", size = 14),  # 坐标轴标题加粗
    axis.text = element_text(face = "bold", size = 12),   # 坐标轴刻度加粗
    legend.title = element_text(face = "bold", size = 12), # 图例标题加粗
    legend.text = element_text(face = "bold", size = 12),  # 图例文本加粗
    plot.title = element_text(face = "bold", size = 16),   # 图形标题加粗
    strip.text = element_text(face = "bold", size = 14)    # Facet标签加粗
  )
  
  # 保存图形，文件名基于基因名
  ggsave(paste0("boxplot_", gene, ".png"), plot = p, width = 8, height = 6)
}

