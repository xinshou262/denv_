library(pheatmap)
library(grid)
library(RColorBrewer)
setwd("D:\\DENV-3\\drug-degcount\\4genes\\对接打分的热图")
# 读取数据
retu <- read.csv("对接打分.csv", row.names = 1, check.names = FALSE)
retu_matrix <- as.matrix(retu)

# 创建显示文本矩阵（0显示为--）
display_num <- ifelse(retu_matrix == 0, "--", sprintf("%.2f", retu_matrix))

# 将0值替换为NA以实现白色背景
retu_na <- retu_matrix
retu_na[retu_na == 0] <- NA

# 生成热图
pdf(file="新热图111.pdf", width=10, height=12)
p <- pheatmap(retu_na,
              border_color = "white",
              cluster_cols = FALSE,
              cluster_rows = FALSE,
              fontsize_row = 12,
              fontsize_col = 16,
              cellwidth = 40,
              cellheight = 20,
              display_numbers = display_num,
              fontsize_number = 10,
              number_color = "red",
              na_col = "#E6E6FA",  # 0值改为浅紫色
              color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))

# 调整列标签角度
g <- p$gtable
col_names_index <- which(g$layout$name == "col_names")
g$grobs[[col_names_index]]$rot <- 45
g$grobs[[col_names_index]]$hjust <- 1

grid.newpage()
grid.draw(g)
dev.off()
