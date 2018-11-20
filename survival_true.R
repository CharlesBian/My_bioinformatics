library(openxlsx)

table <- read.xlsx(choose.files())
Svival <- read.xlsx(choose.files(),2)

library(survival)
library(survminer)
#"55693" : KDM4D
#"1499"  : CTNNB1

onegene_survival(special_geneid="55693")

classify_suvival(special_geneid="55693",
                 classify_geneid= "1499",
                 param= "Down")


onegene_survival <- function(special_geneid)
  {
  special_express <- table[table$gene_id==special_geneid,3:433]
  E <- as.matrix(special_express)
  E <- as.numeric(E)
  median_express <- median(E)
  Svival$special_express <- as.numeric(special_express)  #add to Svival

#create a survival formula
  mysurvival <- survfit(Surv(Svival$start.time,Svival$days,Svival$vital_status)~ Svival$special_express>=median_express,
                        conf.type = "logit")
#caculate the p_value
  p_val <- survdiff(Surv(Svival$days,Svival$vital_status)~ Svival$special_express>=median_express,data = Svival)
  p.value <- 1-pchisq(p_val$chisq,length(p_val$n)-1)
  
#plot the kaplan-Meier curve
  plot(mysurvival,
       main=paste("Kaplan-Meier estimate",
                  "(",table[table$gene_id==special_geneid,1],")",sep = ""),
       xlab = "Time(Days)", 
       ylab = "Survival probability",
       col = c('blue','red'),
       sub = paste("p.value=",p.value))
  legend("topright", 
         legend=c(paste(table[table$gene_id==special_geneid,1]," >=Median","(N=",sum(Svival$special_express>=median_express),")",sep=""),
                  paste(table[table$gene_id==special_geneid,1]," < Median","(N=",sum(Svival$special_express<median_express),")",sep="")),
         col=c("red", "blue"),
         lwd=2)
  }


#classify_gene up and down
classify_suvival <- function(special_geneid, classify_geneid, param)
  {
  special_express <- table[table$gene_id==special_geneid,3:433]
  E <- as.matrix(special_express)
  E <- as.numeric(E)
  median_express <- median(E)
  Svival$special_express <- as.numeric(special_express)   #add to Svival
  
  classify_express <- table[table$gene_id==classify_geneid,3:433]
  C <- as.matrix(classify_express)
  C <- as.numeric(C)
  C_median_express <- median(C)
#subset
  Svival$threshold <- as.factor(ifelse(classify_express>=C_median_express,"Up","Down"))
  Svival <- subset(Svival, threshold == param)   #notice : reset the paramter
  E <- as.matrix(Svival$special_express)
  E <- as.numeric(E)
  median_express <- median(E)

#caculate the p_value
  p_val <- survdiff(Surv(Svival$days,Svival$vital_status)~ Svival$special_express>=median_express,data = Svival)
  p.value <- 1-pchisq(p_val$chisq,length(p_val$n)-1)

  #create a survival formula
  mysurvival <- survfit(Surv(Svival$start.time,Svival$days,Svival$vital_status)~ Svival$special_express>=median_express,
                      conf.type = "logit")

#plot the kaplan-Meier curve
  plot(mysurvival,
     main=paste(table[table$gene_id==classify_geneid,1],"-",Svival$threshold[1]," -- Kaplan-Meier estimate",
                "(",table[table$gene_id==special_geneid,1],")",sep = ""),
     xlab = "Time(Days)", 
     ylab = "Survival probability",
     col = c('blue','red'),
     sub = paste("p.value=",p.value))
  legend("topright", 
       legend=c(paste(table[table$gene_id==special_geneid,1]," >=Median","(N=",sum(Svival$special_express>=median_express),")",sep=""),
                paste(table[table$gene_id==special_geneid,1]," < Median","(N=",sum(Svival$special_express<median_express),")",sep="")),
       col=c("red", "blue"),
       lwd=2)
  }


#create a Cox formula
mycox <- survfit(Surv(Svival$start.time,Svival$days,Svival$vital_status)~ E>=median_express,
                 conf.type = "logit") 
survfitcoxph.fit()

#apply function
p <- matrix()
p <- sapply(table$gene_id,p_value_calculate2)


for(i in 1:20531)
{
  special_express <- table[i,3:433]
  E <- as.matrix(special_express)
  E <- as.numeric(E)
  median_express <- median(E)
  Svival$special_express <- as.numeric(special_express)
  ifelse(
    min(E) >= median_express,
    p.value <- 1,
    p.value <- p_value_calculate())
  p[i]<-p.value
}

classify_geneid= "55693"
param= "Up"

p_calculate <- function()
{
  p_val <- survdiff(Surv(a$days,a$vital_status)~ a$special_express>=median_express,data = a)
  p.value <- 1-pchisq(p_val$chisq,length(p_val$n)-1)
  p.value
}

for(i in 1:20531)
{
  special_express <- table[i,3:433]
  E <- as.matrix(special_express)
  E <- as.numeric(E)
  median_express <- median(E)
  Svival$special_express <- as.numeric(special_express)
  classify_express <- table[table$gene_id==classify_geneid,3:433]
  C <- as.matrix(classify_express)
  C <- as.numeric(C)
  C_median_express <- median(C)
  #subset
  Svival$threshold <- as.factor(ifelse(classify_express>=C_median_express,"Up","Down"))
  a <- subset(Svival, threshold == param)   #notice : reset the paramter
  
  E <- as.matrix(a$special_express)
  E <- as.numeric(E)
  median_express <- median(E)
  ifelse(
    min(E) >= median_express,
    p.value <- 1,          
    p.value <- p_calculate())
  p[i]<-p.value
}
write.xlsx(p,choose.files(),append = FALSE)


for(i in 1:20531)
{
  p.value<-p_value_calculate2(table$gene_id[i])
  p[i]<-p.value
}

table$p.value_survival<-p

#p_value calculate and p_value calculate2 code
p_value_calculate2(1677)

p_value_calculate <- function()
  {
  p_val <- survdiff(Surv(Svival$days,Svival$vital_status)~ Svival$special_express>=median_express,data = Svival)
  p.value <- 1-pchisq(p_val$chisq,length(p_val$n)-1)
  p.value
}

special_geneid<-1

p_value_calculate2 <- function(special_geneid)
{
  special_express <- table[table$gene_id==special_geneid,3:433]
  E <- as.matrix(special_express)
  E <- as.numeric(E)
  median_express <- median(E)
  Svival$special_express <- as.numeric(special_express)
  ifelse(
    min(E) >= median_express,
    p.value <- 1,
    p.value <- p_value_calculate())
  p.value
}

p_value_calculate3()

p_value_calculate3 <- function(special_geneid, classify_geneid, param)
  {
  special_express <- table[table$gene_id==special_geneid,3:433]
  E <- as.matrix(special_express)
  E <- as.numeric(E)
  median_express <- median(E)
  Svival$special_express <- as.numeric(special_express)   #add to Svival
  classify_express <- table[table$gene_id==classify_geneid,3:433]
  C <- as.matrix(classify_express)
  C <- as.numeric(C)
  C_median_express <- median(C)
  #subset
  Svival$threshold <- as.factor(ifelse(classify_express>=C_median_express,"Up","Down"))
  Svival <- subset(Svival, threshold == param)   #notice : reset the paramter
  E <- as.matrix(Svival$special_express)
  E <- as.numeric(E)
  median_express <- median(E)
  
  Svival$special_express <- as.numeric(special_express)
  ifelse(
    min(E) >= median_express,
    p.value <- 1,
    p.value <- p_value_calculate())
  p.value
  }

#survival first classify



for(i in 1:20531)
{
  special_express <- table[i,3:433]
  E <- as.matrix(special_express)
  E <- as.numeric(E)
  median_express <- median(E)
  Svival$special_express <- as.numeric(special_express)
  
  classify_express <- table[table$gene_id==classify_geneid,3:433]
  C <- as.matrix(classify_express)
  C <- as.numeric(C)
  C_median_express <- median(C)
  #subset
  Svival$threshold <- as.factor(ifelse(classify_express>=C_median_express,"Up","Down"))
  a <- subset(Svival, threshold == param)   #notice : reset the paramter
  
  E <- as.matrix(a$special_express)
  E <- as.numeric(E)
  median_express <- median(E)
  ifelse(
    min(E) >= median_express,
    p.value <- 1,          
    p.value <- p_calculate())
  p[i]<-p.value
}

special_geneid="55693"
special_express <- table[table$gene_id==special_geneid,3:433]
E <- as.matrix(special_express)
E <- as.numeric(E)
median_express <- median(E)
quan <- quantile(special_express)
Svival$special_express <- as.numeric(special_express)  #add to Svival
Svival$status <- as.factor(ifelse(Svival$special_express<=quan$`75%`,
                                   ifelse(Svival$special_express<=quan$`25%`,"down","median"),"up"))

Svival$status <- as.factor(ifelse(Svival$special_express<=quan$`75%`,"Middle&down","up"))

#create a survival formula
mysurvival <- survfit(Surv(Svival$start.time,Svival$days,Svival$vital_status)~ Svival$status,
                        conf.type = "logit")
#caculate the p_value
p_val <- survdiff(Surv(Svival$days,Svival$vital_status)~ Svival$status,data = Svival)
p.value <- 1-pchisq(p_val$chisq,length(p_val$n)-1)
  
#plot the kaplan-Meier curve
plot(mysurvival,
      main=paste("Kaplan-Meier estimate",
                  "(",table[table$gene_id==special_geneid,1],")",sep = ""),
      xlab = "Time(Days)", 
      ylab = "Survival probability",
      col = c('blue','red'),
      sub = paste("p.value=",p.value))
legend("topright", 
        legend=c(paste(table[table$gene_id==special_geneid,1]," >=75%","(N=",sum(Svival$status=="up"),")",sep=""),
   #             paste(table[table$gene_id==special_geneid,1]," = middle","(N=",sum(Svival$status=="median"),")",sep=""),
                paste(table[table$gene_id==special_geneid,1]," <=25%","(N=",sum(Svival$status=="down"),")",sep="")),
       col=c("red","green","blue"),
       lwd=2)





















