library(ggVolcano)
library(openxlsx)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

setwd("D:\\芯片数据集医院数据\\GSE6631训练集\\差异表达基因")

# 读取文件
degdata <- read.xlsx("差异logfc.xlsx", sheet = 1)
# 设定差异基因阈值
degdata$status <- 'stable'
degdata$status[degdata$logFC>0.585 & degdata$adj.P.Val<0.05] <- 'up'
degdata$status[degdata$logFC< -0.585 & degdata$adj.P.Val<0.05] <- 'down'

# 挑选出差异基因子集，用于后续添加label
labeldata <- subset(degdata, abs(logFC)>0.585 & adj.P.Val<0.05)
id <- order(-log10(labeldata$adj.P.Val),decreasing = T)
labeldata <- labeldata[id,]
head(labeldata)
# 计算Log2FC的最大绝对值
max_abs_logFC <- max(abs(degdata$logFC), na.rm = TRUE)
# 设置对称范围（例如取整或略大于最大值）
x_limit <- c(-max_abs_logFC, max_abs_logFC)
# 若需要更美观的整数范围，可手动调整：
# x_limit <- c(-3, 3)  # 根据实际情况调整数值
# (3)彩色渐变版####
ggplot(data = degdata,
       mapping = aes(
         x = logFC,
         y = -log10(adj.P.Val)
       )) +
  # 绘制散点（保留颜色和大小梯度）
  geom_point(aes(
    color = -log10(adj.P.Val),
    size = -log10(adj.P.Val)
  ),
  alpha = 1) +
  
  # 绘制垂直线（FC阈值）也就是log2（x）=0.5
  geom_vline(
    xintercept = c(log2(0.7071068), log2(1.414214)),
    linetype = "dashed",
    color = "#363636"
  ) +
  
  # 绘制水平线（显著性阈值）
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "#363636"
  ) +
  
  theme_bw() +
  
  # 设置颜色渐变（RdYlBu调色板反向）
  scale_color_gradientn(colours = brewer.pal(11, "RdYlBu") %>% rev()) +
  scale_x_continuous(limits = x_limit)+
  # 设置点的大小范围
  scale_size_continuous(range = c(0.3, 3)) +
  
  # 删除以下注释的标签图层：
  # geom_text_repel(
  #   data = labeldata[1:10, ],
  #   mapping = aes(label = label, color = -log10(adj.P.Val)),
  #   size = 2.5
  # ) +
  
  # 坐标轴标签
  xlab("Log2FC") +
  ylab("-Log10(adj.Pvalue)") +
  
  # 图例设置
  guides(
    size = FALSE,
    color = guide_colorbar(title = "-Log10(adjPvalue)")
  )

ggsave("plot3center.pdf", width = 7, height = 5.5)


             