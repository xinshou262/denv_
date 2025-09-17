#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af
rm(list = ls())
setwd("D:\\DENV-3\\drug-degcount\\4genes\\单基因GSEA\\高通量单基因")
library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
library(DESeq2)
#输入文件
expFile="有基因名的所有患病样本count.txt"         
sgene="EPHB2"       #输入进行单基因GSEA的基因名称

rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)# 检查第一列是否有重复名
rownames(rt)=rt[,1]
rt <- rt[, -1]
# 转换为数值矩阵
data <- as.matrix(rt)

# 将矩阵转换为字符型（处理特殊符号）
data_char <- matrix(gsub("[^0-9.-]", "", data), nrow = nrow(data))

# 再转换为数值矩阵
data_numeric <- matrix(as.numeric(data_char), nrow = nrow(data))

# 保留行列名
rownames(data_numeric) <- rownames(data)
colnames(data_numeric) <- colnames(data)

# 再转换为数值矩阵
data_numeric <- matrix(as.numeric(data_char), nrow = nrow(data))

# 保留行列名
rownames(data_numeric) <- rownames(data)
colnames(data_numeric) <- colnames(data)

# 检查矩阵类型
class(data_numeric)  # 应为 "matrix" "array"
typeof(data_numeric) # 应为 "double"

# 查看目标基因的表达值（应无引号/空格）
data_numeric[sgene, ]

data <- data_numeric
# 根据目标基因（CXCL10）表达分组
sgene <- "EPHB2"
group <- ifelse(data[sgene,] > median(data[sgene,]), "High", "Low")   
group <- factor(group,levels = c("High","Low"))

# 创建正式的colData对象（关键步骤）
coldata <- data.frame(
  condition = factor(group, levels = c("High", "Low")),
  row.names = colnames(data)  # 必须与data矩阵的列名（样本名）一致
)

# 构建DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData = coldata,
  design = ~ condition
)

# 检查重复的基因名
dup_genes <- rownames(data)[duplicated(rownames(data))]
print(dup_genes)  # 输出重复的基因名（示例：[1] "CXCL10" "TP53"）


# 执行DESeq2标准化和差异分析
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","High","Low"))
write.csv(as.data.frame(res), file = "DESeq单基因中值分组差异分析结果.csv")
# 重命名列以匹配 limma 代码
res_df <- as.data.frame(res) %>%
  dplyr::rename(
    logFC = log2FoldChange,  # 对应 limma 的 logFC
    P.Value = pvalue,        # 对应 limma 的 P.Value
    adj.P.Val = padj         # 对应 limma 的 adj.P.Val
  )

# 添加基因符号列（假设行名为基因名）
res_df$symbol <- rownames(res_df)
library(clusterProfiler)
library(org.Hs.eg.db)

# 转换基因符号到 Entrez ID
df <- bitr(unique(res_df$symbol),
           fromType = "SYMBOL",
           toType = c("ENTREZID"),
           OrgDb = org.Hs.eg.db)

# 合并 Entrez ID 到 DESeq2 结果
DEG <- merge(res_df, df, by.x = "symbol", by.y = "SYMBOL")

# 按 logFC 降序排序
data_all_sort <- DEG %>% 
  dplyr::arrange(desc(logFC))




# 生成geneList,#把foldchange按照从大到小提取出来
geneList = data_all_sort$logFC
names(geneList) <- data_all_sort$ENTREZID

# 显式降序排序（关键！）
geneList <- sort(geneList, decreasing = TRUE)

# 过滤缺失的ENTREZID和重复基因
geneList <- geneList[!is.na(names(geneList))]
geneList <- geneList[!duplicated(names(geneList))]

# 运行gseKEGG
kk2 <- gseKEGG(
  geneList = geneList,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 200,
  pvalueCutoff = 0.05,
  pAdjustMethod = "none"
)
class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af=as.data.frame(kk2@result)
write.table(af,file=paste0("2.","all_GSEA的结果.xls"),sep="\t",quote=F,col.names=T)


#排序后分别取GSEA结果的前5个和后5个
num=5
pdf(paste0("EPHB2.","1down_GSEA.pdf"),width = 8,height = 8)
gseaplot2(kk2, 
          geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore),num)],
          base_size = 20)  # 增大所有文本大小
dev.off()

pdf(paste0("EPHB2.","1up_GSEA.pdf"),width = 8,height = 8)
gseaplot2(kk2, 
          geneSetID = rownames(kk2@result)[tail(order(kk2@result$enrichmentScore),num)],
          base_size = 20)  # 增大所有文本大小
dev.off()

#排序后取前5个和后5个一起展示
num=5
pdf(paste0("EPHB2.","1all_GSEA.pdf"),width = 14,height = 12)  # 增大图像尺寸以容纳更大文本
gseaplot2(kk2, 
          geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))],
          base_size = 20)  # 增大所有文本大小
dev.off()




