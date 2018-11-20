################################################
### 作者：徐洲更                             ###
### 更新时间：2018-11-07                     ###
################################################

## 向量篇
### 内置向量
letters # 小写字母
LETTERS # 大写字母
month.abb  # 月份缩写
month.name # 月份全称

### 创建向量
c("xu", "zhou", "geng")
1:10
seq(1,10, 2) # 每隔2个数字
?seq
rep(c("a","b","c"), 2) # 重复

### 向量索引
letters[10]
letters[1:10] 
letters[seq(1,10,2)] 
letters[seq(1,10,2) > 2]


### 向量运算
1:10 + 11:20 
11:20 - 1:10 
11:20 %% 2 == 0

### 变量赋值
myname <- c("xu", "zhou", "geng")
c("xu", "zhou", "geng") -> yourname

## 矩阵篇
library(Matrix)

m <- matrix(1:20, nrow = 5, ncol =4 ) # 创建矩阵
 
Matrix::rowSums(m)  #行和
Matrix::rowMeans(m) #行均值
Matrix::colSums(m)  #列和
Matrix::colMeans(m) #列均值


## 设置项目路径
getwd()
setwd("/Users/zgxu/Desktop/R语言基础/")

## 文件读取
### 数据存储无非两种方式: 二进制, 文本
### 核心函数 read.table
### 其他函数 read.xxx
options(stringsAsFactors = FALSE)
expr_df <- read.csv("exprData.csv", header = FALSE
                          )
expr_df[1:10,1:6]                          
meta_data <- read.table("metadata.txt",sep = "\t", header = TRUE)
head(meta_data)

gene_names <- readLines(con = "gene_name.txt")
head(gene_names)

## 数据框篇

is.data.frame(expr_df) # 判断是否是数据框
dim(expr_df) # 行列数


### 变量重命名
row.names(expr_df) <- gene_names # 行名
colnames(expr_df) <- meta_data$Run #列名
head(expr_df)

#### 根据dex进行命名 
new_name <- paste(meta_data$dex,rep(c(1,2,3,4),each=2), sep = "_") # 合并字符串
colnames(expr_df) <- new_name
head(expr_df)

### 创建新变量: 根据现有值进行变换
#### 形式: 变量名 <- 表达式

#### 案例1: 获取一些关于基因和样本的元信息
total_count_per_gene <- rowSums(expr_df) #每个基因的表达量
total_count_per_sample <- colSums(expr_df) # 每个样本的表达量


### 案例2: 根据深度深度进行标准化
new_df <- expr_df / rep(total_count_per_sample, each=nrow(expr_df)) * 1000000

### 类型转换
#### isis..xxxx 用于判断, as.xxxx用于转换
is.matrix(new_df)
expr_matrix <- as.matrix(new_df)
is.matrix(expr_matrix)

### 数据排序

#### 案例: 按照基因的表达量总数对表达矩阵排序

df_ordered <- expr_df[order(total_count_per_gene),] # 升序
df_ordered <- expr_df[order(total_count_per_gene, decreasing = T),] # 降序

### 数据集合并

#### 案例：增加样本元信息

cbind(meta_data, total_count_per_sample) # 增加列, column
meta_data$total_count_per_sample <- total_count_per_sample  # 增加列

#### 案例: 增加每个样本的总数行
df_tmp <- rbind(expr_df, total_count_per_sample) # 增加行, row
tail(df_tmp)

### 数据集取子集

### 案例： 提取表达量前10的基因
top_10 <- df_ordered[1:10, ]

#### 案例：过滤低表达的基因
gene_to_filter <- total_count_per_gene > 10 
df_filterd <- expr_df[gene_to_filter, ]
dim(df_filterd)

#### 案例: 随机抽样
set.seed(19930831)
df_sample <- expr_df[sample(1:nrow(expr_df), 1000, replace = FALSE ), ]
dim(df_sample)
