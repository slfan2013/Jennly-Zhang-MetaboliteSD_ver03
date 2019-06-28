label_data = fread("Metabolite Label_KOMP correlation.csv", header = FALSE)


# 6 in null, the correlation phenotype metabolite (correlation male).csv

pacman::p_load(data.table, corrplot)
cor_data_male= fread("6 in null, the correlation phenotype metabolite (correlation male).csv")


cor_data_male = merge(cor_data_male, label_data, by.x = "V1", by.y = "V2", sort = FALSE)

cor_data_male_label = cor_data_male[["V3"]]
cor_male = data.matrix(cor_data_male[,-c(1, ncol(cor_data_male), ncol(cor_data_male)-1),with=FALSE])
# rownames(cor_male) = colnames(cor_male) = NULL
rownames(cor_male) = cor_data_male_label
cor_p_data_male= fread("6 in null, the correlation phenotype metabolite (correlation p value male).csv")
cor_p_male = data.matrix(cor_p_data_male[,-1,with=FALSE])

# 6 in null, the correlation phenotype metabolite (correlation female).csv

cor_data_female= fread("6 in null, the correlation phenotype metabolite (correlation female).csv")


cor_data_female = merge(cor_data_female, label_data, by.x = "V1", by.y = "V2", sort = FALSE)

cor_data_female_label = cor_data_female[["V3"]]
cor_female = data.matrix(cor_data_female[,-c(1, ncol(cor_data_female), ncol(cor_data_female)-1),with=FALSE])
# rownames(cor_female) = colnames(cor_female) = NULL
rownames(cor_female) = cor_data_female_label
cor_p_data_female= fread("6 in null, the correlation phenotype metabolite (correlation p value female).csv")
cor_p_female = data.matrix(cor_p_data_female[,-1,with=FALSE])








################## DIFFERENT DIRECTION ##################
different_direction_index = which(cor_male * cor_female < 0 & cor_p_female < 0.05 & cor_p_male < 0.05, arr.ind = T)
different_direction_index_row_unique = unique(different_direction_index[,1])
different_direction_index_col_unique = unique(different_direction_index[,2])

p_data = cor_p_male[different_direction_index_row_unique,different_direction_index_col_unique]

par(mfrow=c(1,1))
data = cor_female[different_direction_index_row_unique,different_direction_index_col_unique] - cor_male[different_direction_index_row_unique,different_direction_index_col_unique]

data_row_order =  hclust(dist(data))$order
data_col_order =  hclust(dist(t(data)))$order



fake_p = matrix(0, nrow = nrow(cor_female), ncol = ncol(cor_female))
for(i in 1:nrow(different_direction_index)){
  fake_p[different_direction_index[i,1],different_direction_index[i,2]] = 1
}

fake_p = fake_p[different_direction_index_row_unique, different_direction_index_col_unique]

 

corrplot(data[data_row_order,data_col_order], is.corr = FALSE, insig = c("pch"), p.mat = fake_p[data_row_order,data_col_order], pch = 22, pch.col = "black", pch.cex = 2)


par(mfrow=c(1,1))
corrplot(cor_male[different_direction_index_row_unique,different_direction_index_col_unique][data_row_order,data_col_order], p.mat = fake_p[data_row_order,data_col_order], method="circle", pch = 22, pch.col = "black", pch.cex = 2)

corrplot(cor_female[different_direction_index_row_unique,different_direction_index_col_unique][data_row_order,data_col_order], , p.mat = fake_p[data_row_order,data_col_order], method="circle", pch = 22, pch.col = "black", pch.cex = 2)


################## DIFFERENT SIZE ################## 
different_size_index = which(cor_male * cor_female > 0 & cor_p_female < 0.05 & cor_p_male < 0.05, arr.ind = T)
different_size_index_row_unique = unique(different_size_index[,1])
different_size_index_col_unique = unique(different_size_index[,2])

p_data = cor_p_male[different_size_index_row_unique,different_size_index_col_unique]

par(mfrow=c(1,1))
data = cor_female[different_size_index_row_unique,different_size_index_col_unique] - cor_male[different_size_index_row_unique,different_size_index_col_unique]

data_row_order =  hclust(dist(data))$order
data_col_order =  hclust(dist(t(data)))$order

fake_p = matrix(0, nrow = nrow(cor_female), ncol = ncol(cor_female))
for(i in 1:nrow(different_size_index)){
  fake_p[different_size_index[i,1],different_size_index[i,2]] = 1
}

fake_p = fake_p[different_size_index_row_unique, different_size_index_col_unique]


corrplot(data[data_row_order,data_col_order], is.corr = FALSE, insig = c("pch"), p.mat = fake_p[data_row_order,data_col_order], pch = 22, pch.col = "black", pch.cex = 2)


par(mfrow=c(1,1))
corrplot(cor_male[different_size_index_row_unique,different_size_index_col_unique][data_row_order,data_col_order], p.mat = fake_p[data_row_order,data_col_order], method="circle", pch = 22, pch.col = "black", pch.cex = 5)

corrplot(cor_female[different_size_index_row_unique,different_size_index_col_unique][data_row_order,data_col_order], p.mat = fake_p[data_row_order,data_col_order], method="circle", pch = 22, pch.col = "black", pch.cex = 5)


################## SIGNIFICANT IN ONE GENDER FEMALE ################## 
one_gender_female_index = which(cor_p_female < 0.05 & cor_p_male > 0.05, arr.ind = T)
one_gender_female_index_row_unique = unique(one_gender_female_index[,1])
one_gender_female_index_col_unique = unique(one_gender_female_index[,2])

p_data = cor_p_male[one_gender_female_index_row_unique,one_gender_female_index_col_unique]

par(mfrow=c(1,1))
data = cor_female[one_gender_female_index_row_unique,one_gender_female_index_col_unique] - cor_male[one_gender_female_index_row_unique,one_gender_female_index_col_unique]

data_row_order =  hclust(dist(data))$order
data_col_order =  hclust(dist(t(data)))$order

fake_p = matrix(0, nrow = nrow(cor_female), ncol = ncol(cor_female))
for(i in 1:nrow(one_gender_female_index)){
  fake_p[one_gender_female_index[i,1],one_gender_female_index[i,2]] = 1
}

fake_p = fake_p[one_gender_female_index_row_unique, one_gender_female_index_col_unique]


corrplot(data[data_row_order,data_col_order], is.corr = FALSE, insig = c("pch"), p.mat = fake_p[data_row_order,data_col_order], pch = 22, pch.col = "black", pch.cex = 2)


par(mfrow=c(1,1))
corrplot(cor_female[one_gender_female_index_row_unique,one_gender_female_index_col_unique][data_row_order,data_col_order], p.mat = fake_p[data_row_order,data_col_order], method="circle", pch = 22, pch.col = "red", pch.cex = 5, tl.cex = 2)


################## SIGNIFICANT IN ONE GENDER MALE ################## 
one_gender_male_index = which(cor_p_female > 0.05 & cor_p_male < 0.05, arr.ind = T)
one_gender_male_index_row_unique = unique(one_gender_male_index[,1])
one_gender_male_index_col_unique = unique(one_gender_male_index[,2])

p_data = cor_p_male[one_gender_male_index_row_unique,one_gender_male_index_col_unique]

par(mfrow=c(1,1))
data = cor_female[one_gender_male_index_row_unique,one_gender_male_index_col_unique] - cor_male[one_gender_male_index_row_unique,one_gender_male_index_col_unique]

data_row_order =  hclust(dist(data))$order
data_col_order =  hclust(dist(t(data)))$order

fake_p = matrix(0, nrow = nrow(cor_female), ncol = ncol(cor_female))
for(i in 1:nrow(one_gender_male_index)){
  fake_p[one_gender_male_index[i,1],one_gender_male_index[i,2]] = 1
}

fake_p = fake_p[one_gender_male_index_row_unique, one_gender_male_index_col_unique]


corrplot(data[data_row_order,data_col_order], is.corr = FALSE, insig = c("pch"), p.mat = fake_p[data_row_order,data_col_order], pch = 22, pch.col = "black", pch.cex = 2)


par(mfrow=c(1,1))
corrplot(cor_male[one_gender_male_index_row_unique,one_gender_male_index_col_unique][data_row_order,data_col_order], p.mat = fake_p[data_row_order,data_col_order], method="circle", pch = 22, pch.col = "red", pch.cex = 5, tl.cex = 2)



























































