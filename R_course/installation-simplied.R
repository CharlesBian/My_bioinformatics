################################################
### 作者：果子
### 更新时间：2018-11-7
### 联系方式: guoshipeng2008@126.com
### 微信：guo2sky

### 一定要认真阅读：遇到不确定的情况退回去看视频。 
### 注意哦！！ 这里面有的有可能和视频里面文字顺序不一样，不要紧的。
### 假如打开后是乱码，也不要紧，点击左上方file，选择Reopen with encodings
### 选择UTF-8

#bioconductor
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
#CRAN
options(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")


# 检查是否设定完毕
getOption("BioC_mirror")
getOption("CRAN")

## 查看R包安装的位置
.libPaths()
#如果只出现一个地点，也是正常的。

###########################
# 上面完成不了也不要紧，从这里开始我们的旅程
# 配置bioconductor安装环境,
# 视频中说的删除，然后tab键可以补齐，如果一不小心删除错了，不知道如何处理。
# 此时，不确定就不要保存这个文件，退出重新进入。
install.packages("BiocManager")
library("BiocManager")

##安装bioconductor
BiocManager::install("DESeq2")
library(DESeq2)

##安装CRAN上的包
install.packages("dplyr")
library("dplyr")
install.packages("ggplot2")
library("ggplot2")
