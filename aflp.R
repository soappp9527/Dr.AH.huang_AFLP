library(data.table)
library(poppr)
library(phangorn)
aflp_data <- fread("aflp_data.csv", stringsAsFactors = FALSE)
aflp_data$bp <- as.numeric(gsub(",", "", aflp_data$bp))
aflp_data <- aflp_data[bp < 2000 & conc > 2,]
aflp_data <- aflp_data[!(bp == 104 & sample == "sp001"),]
#
aflp_data[bp == 39,]$bp <- 40
aflp_data[bp %in% 86:92,]$bp <- 90
aflp_data[bp %in% 100:110,]$bp <- 105
aflp_data[bp %in% 231:240,]$bp <- 235
aflp_data[bp %in% 243:256,]$bp <- 250
aflp_data[bp %in% 261:268,]$bp <- 265
aflp_data[bp %in% 271:283,]$bp <- 275
aflp_data[bp %in% 296:308,]$bp <- 300
aflp_data[bp %in% 342:347,]$bp <- 345
aflp_data[bp %in% 371:378,]$bp <- 375
aflp_data[bp %in% 417:426,]$bp <- 420
aflp_data[bp %in% 427:438,]$bp <- 430
aflp_data[bp %in% 599,]$bp <- 602
aflp_data[bp %in% 618,]$bp <- 620
aflp_data[bp %in% 912:919,]$bp <- 915
aflp_data[bp %in% 923:926,]$bp <- 925
aflp_data[bp %in% 929:934,]$bp <- 930
aflp_data[bp %in% 1097:1098,]$bp <- 1098
aflp_data[bp %in% 1120:1126,]$bp <- 1125
aflp_data[bp %in% 1134:1139,]$bp <- 1135
aflp_data[bp %in% 1165:1168,]$bp <- 1165
#
aflp_matrix <- data.frame(dcast(aflp_data, sample~bp))
aflp_matrix <- data.frame(aflp_matrix[,-1], row.names = aflp_matrix[,1])
aflp_matrix[is.na(aflp_matrix)] <- 0
aflp_matrix[aflp_matrix != 0] <- 1
#
aflp_nei <- nei.dist(as.matrix(aflp_matrix))
aflp_upgma <- upgma(aflp_nei)
#
gd_matrix <- data.frame(t(combn(rownames(aflp_matrix),2)), as.numeric(aflp_nei))
names(gd_matrix) <- c("sample1", "sample2", "genetic distance")
write.csv(gd_matrix, file = "gd_matrix.csv")

#
png(file = "phylogram.png", height = 1500, width = 1000)
plot(aflp_upgma, type="phylogram", cex = 2)
dev.off()

png(file = "fan.png", height = 1000, width = 1000)
plot(aflp_upgma, type="fan",cex = 2)
dev.off()
