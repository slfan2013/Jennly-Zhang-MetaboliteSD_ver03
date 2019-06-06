Null_Metabolomics = wcmc::read_data("Null_Metabolomics.xlsx")
Null_Metabolomics_p = Null_Metabolomics$p
Null_Metabolomics_f = Null_Metabolomics$f
Null_Metabolomics_e = Null_Metabolomics$e_matrix





## 1 wild - sex
### deal with missing value. 
#### filter missing value > 80%


num_missing_male = apply(Null_Metabolomics_e,1,function(x){
  sum(is.na(x[Null_Metabolomics_p$Gender %in% "Male"]))
})
num_missing_female = apply(Null_Metabolomics_e,1,function(x){
  sum(is.na(x[Null_Metabolomics_p$Gender %in% "Female"]))
})
names(num_missing_male) = names(num_missing_female) = Null_Metabolomics_f$label
missing_index = ((num_missing_male/sum(Null_Metabolomics_p$Gender %in% "Male")) > 0.8)|((num_missing_female/sum(Null_Metabolomics_p$Gender %in% "Female")) > 0.8)
Null_Metabolomics_e = Null_Metabolomics_e[!missing_index,]
Null_Metabolomics_f = Null_Metabolomics_f[!missing_index,]




#### filter missing value > 80%
for(i in 1:nrow(Null_Metabolomics_e)){
  Null_Metabolomics_e[i,is.na(Null_Metabolomics_e[i,])] = 0.5 * min(Null_Metabolomics_e[i,], na.rm = TRUE)
}
rownames(Null_Metabolomics_e) = Null_Metabolomics_f$label

## statistics
p_val1 = apply(Null_Metabolomics_e,1,function(x){
  summary(lm(x ~ Null_Metabolomics_p$Gender))$coefficients[2,4]
})
fc1 = apply(Null_Metabolomics_e,1,function(x){
  mean(x[Null_Metabolomics_p$Gender %in% "Male"])/mean(x[Null_Metabolomics_p$Gender %in% "Female"])
})
p_val_adj1 = p.adjust(p_val1, method = "fdr")

sum(p_val1<0.05)/length(p_val1)

fwrite(data.table(label = Null_Metabolomics_f$label, p_value = p_val1, fold_change = fc1, adjusted_p_value = p_val_adj1),"1 wild - sex.csv")





## 2 adjust by body weight
Phenotype_UCDavisSelect = wcmc::read_data("Phenotype_UCDavisSelect.xlsx")
Phenotype_UCDavisSelect_p = Phenotype_UCDavisSelect$p
Phenotype_UCDavisSelect_f = Phenotype_UCDavisSelect$f
Phenotype_UCDavisSelect_e = Phenotype_UCDavisSelect$e_matrix

body_weight_all = Phenotype_UCDavisSelect_e[,Phenotype_UCDavisSelect_p[["phenotype"]] %in% "Body weight"]
names(body_weight_all) = Phenotype_UCDavisSelect_f$label
body_weight = body_weight_all[Null_Metabolomics_p$label]

p_val2 = apply(Null_Metabolomics_e,1,function(x){
  summary(lm(x ~ Null_Metabolomics_p$Gender + body_weight))$coefficients[2,4]
})
fc2 = apply(Null_Metabolomics_e,1,function(x){
  summary(lm(x ~ Null_Metabolomics_p$Gender + body_weight))$coefficients[2,1]
})
p_val_adj2 = p.adjust(p_val2,"fdr")

sum(p_val2<0.05)/length(p_val2)
fwrite(data.table(label = rownames(Null_Metabolomics_e), p_value = p_val2, fold_change = fc2, adjusted_p_value = p_val_adj2),"2 adjust by body weight.csv")


## 3 phenotype -  sex
Raw_Phenotype_UCDavis = wcmc::read_data("Raw_Phenotype_UCDavis.xlsx")
Raw_Phenotype_UCDavis_p = Raw_Phenotype_UCDavis$p
Raw_Phenotype_UCDavis_f = Raw_Phenotype_UCDavis$f
Raw_Phenotype_UCDavis_e = Raw_Phenotype_UCDavis$e_cat_matrix



Raw_Phenotype_UCDavis_f$label_procedure = sapply(strsplit(Raw_Phenotype_UCDavis_f$label,"___"),function(x){x[1]})
Raw_Phenotype_UCDavis_f$label_phenotype = sapply(strsplit(Raw_Phenotype_UCDavis_f$label,"___"),function(x){x[2]})

rownames(Raw_Phenotype_UCDavis_e) = Raw_Phenotype_UCDavis_f$label

# combine the phenotypes because one phenotype may have multiple procedures.
unique_phenotype = unique(Raw_Phenotype_UCDavis_f$label_phenotype)
for(u in 1:length(unique_phenotype)){
  if(sum(Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u])>1){
    if(length(table(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],]))>1){
      duplicated_value_index = which(apply(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],],2,function(x){
        if(length(unique(x)) == 1){
          if(!unique(x)=="NA"){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }else{
          return(FALSE)
        }
      }))
      if(length(duplicated_value_index)>0){
        print(unique_phenotype[u])
        print(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],1:10])
        Sys.sleep(20)
      }
    }
  }
}




Raw_Phenotype_UCDavis_e_merge = matrix(NA, nrow = length(unique(Raw_Phenotype_UCDavis_f$label_phenotype)), ncol = ncol(Raw_Phenotype_UCDavis_e))
rownames(Raw_Phenotype_UCDavis_e_merge) = unique(Raw_Phenotype_UCDavis_f$label_phenotype)

i = 1
for(u in 1:length(unique_phenotype)){
  
  if(sum(Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u])>1){
    # stop(u)
    for(j in 1:ncol(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],])){
      the_j_th_col = Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],j]
      if(!length(unique(the_j_th_col)) == 1){
        Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],][the_j_th_col == "NA",j] = Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],][!the_j_th_col == "NA",j] 
      } 
    }
    Raw_Phenotype_UCDavis_e_merge[i,] = unique(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],])
    i = i + 1
  }else{
    Raw_Phenotype_UCDavis_e_merge[i,] = Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],]
    i = i + 1
  }
}



# 
# for(i in 1:nrow(Raw_Phenotype_UCDavis_e_merge)){
#   
#   if(sum(is.na(as.numeric(Raw_Phenotype_UCDavis_e_merge[i,]))) > nrow(Raw_Phenotype_UCDavis_p)/2){
#     print(i)
#     print(table(Raw_Phenotype_UCDavis_e_merge[i,]))
#     cat("\n")
#     Sys.sleep(2)
#   }
#   
# }
# 
# 
Raw_Phenotype_UCDavis_e_merge[Raw_Phenotype_UCDavis_e_merge == "NA"] = NA


Raw_Phenotype_UCDavis_e_null = Raw_Phenotype_UCDavis_e_merge[,Raw_Phenotype_UCDavis_p$Genotype %in% "(null)"]
Raw_Phenotype_UCDavis_p_null = Raw_Phenotype_UCDavis_p[Raw_Phenotype_UCDavis_p$Genotype %in% "(null)",]

p_val3 = fc3 = c()
continuous_index = c()
for(i in 1:nrow(Raw_Phenotype_UCDavis_e_merge)){
  
  continuous = sum(is.na(as.numeric(Raw_Phenotype_UCDavis_e_merge[i,!is.na(Raw_Phenotype_UCDavis_e_merge[i,])]))) == 0
  
  
  if(continuous){
    continuous_index[i] = TRUE
    # p_val3[length(p_val3)+1] = tryCatch(summary(lm(as.numeric(Raw_Phenotype_UCDavis_e_null[i,])~ Raw_Phenotype_UCDavis_p_null$Gender))$coef[2,4],error = function(er){
    #   return(NA)
    # })
    
    
    p_val3[length(p_val3)+1] = tryCatch(
      summary(lm(as.numeric(Raw_Phenotype_UCDavis_e_null[i,])~ Raw_Phenotype_UCDavis_p_null$Gender))$coef[2,4]
      ,error = function(er){
      return(NA)
    })
    
    fc3[length(fc3)+1] = tryCatch(
      summary(lm(as.numeric(Raw_Phenotype_UCDavis_e_null[i,]) ~ Raw_Phenotype_UCDavis_p_null$Gender))$coef[2,1]
      ,error = function(er){
        return(NA)
      })
    
    
  }else{
    continuous_index[i] = FALSE
    
    
    p_val3[length(p_val3)+1] = summary(xtabs(~Raw_Phenotype_UCDavis_e_null[i,] + Raw_Phenotype_UCDavis_p_null$Gender))$p.value
    fc3[length(fc3)+1] = summary(xtabs(~Raw_Phenotype_UCDavis_e_null[i,] + Raw_Phenotype_UCDavis_p_null$Gender))$statistic
  }
  
}



sum(p_val3[continuous_index]<0.05,na.rm = TRUE)/sum(continuous_index)




fwrite(data.table(label = rownames(Raw_Phenotype_UCDavis_e_merge), p_value = p_val3, fold_change = fc3, adjusted_p_value = p.adjust(p_val3,'fdr')),"3 phenotype -  sex.csv")




## 4 phenotype -  sex adjusted by body weight
# body_weight_all
# sum(Raw_Phenotype_UCDavis_p_null$label%in%names(body_weight_all))/length(Raw_Phenotype_UCDavis_p_null$label)

body_weight = body_weight_all[Raw_Phenotype_UCDavis_p_null$label]


p_val4 = fc4 = c()
continuous_index = c()
for(i in 1:nrow(Raw_Phenotype_UCDavis_e_merge)){
  
  continuous = sum(is.na(as.numeric(Raw_Phenotype_UCDavis_e_merge[i,!is.na(Raw_Phenotype_UCDavis_e_merge[i,])]))) == 0
  
  
  if(continuous){
    continuous_index[i] = TRUE
    # p_val4[length(p_val4)+1] = tryCatch(summary(lm(as.numeric(Raw_Phenotype_UCDavis_e_null[i,])~ Raw_Phenotype_UCDavis_p_null$Gender))$coef[2,4],error = function(er){
    #   return(NA)
    # })
    
    
    p_val4[length(p_val4)+1] = tryCatch(
      summary(lm(as.numeric(Raw_Phenotype_UCDavis_e_null[i,])~ Raw_Phenotype_UCDavis_p_null$Gender+body_weight))$coef[2,4]
      ,error = function(er){
        return(NA)
      })
    
    fc4[length(fc4)+1] = tryCatch(
      summary(lm(as.numeric(Raw_Phenotype_UCDavis_e_null[i,])~ Raw_Phenotype_UCDavis_p_null$Gender+body_weight))$coef[2,1]
      ,error = function(er){
        return(NA)
      })
    
  }else{
    continuous_index[i] = FALSE
    
    
    p_val4[length(p_val4)+1] = summary(xtabs(~Raw_Phenotype_UCDavis_e_null[i,] + Raw_Phenotype_UCDavis_p_null$Gender + body_weight))$p.value
    fc4[length(fc4)+1] = summary(xtabs(~Raw_Phenotype_UCDavis_e_null[i,] + Raw_Phenotype_UCDavis_p_null$Gender + body_weight))$statistic
  }
  
}

sum(p_val4[continuous_index]<0.05,na.rm = TRUE)/sum(continuous_index)
fwrite(data.table(label = rownames(Raw_Phenotype_UCDavis_e_merge), p_value = p_val4, fold_change = fc4, adjusted_p_value = p.adjust(p_val4,'fdr')),"4 phenotype -  sex adjusted by body weight.csv")





# 5 for each genotype, check which metabolite are associated with gene * gender, gene, gender, or not significant at all.

All_Metabolomics = wcmc::read_data("All_Metabolomics.xlsx")
All_Metabolomics_p = All_Metabolomics$p
All_Metabolomics_f = All_Metabolomics$f
All_Metabolomics_e = All_Metabolomics$e_matrix

unique_genes = unique(All_Metabolomics_p$Genotype)


All_Metabolomics_e_no_mising = All_Metabolomics_e

All_Metabolomics_e[!is.finite(All_Metabolomics_e)] = NA
# All_Metabolomics_e[!is.finite(All_Metabolomics_e)] = 0



for(g in 1:length(unique_genes)){
  
  current_gene = unique_genes[g]
  current_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c(current_gene)]
  current_e = All_Metabolomics_e[, current_label]
  
  # deal with missing value here.
  for(i in 1:nrow(current_e)){
    
    if(sum(is.na(current_e[i,]))>(ncol(current_e)-2)){
      current_e[i,] = NA
    }else{
      current_e[i,is.na(current_e[i,])] = runif(sum(is.na(current_e[i,])),min=0.49,max=0.51) * min(current_e[i,!is.na(current_e[i,])])
    }
    
  }
  
  
  All_Metabolomics_e_no_mising[,current_label] = current_e
  
}
null_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c('null')]
null_e = All_Metabolomics_e[, null_label]
null_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c('null') & All_Metabolomics_p$Gender %in% "Female"]
null_e_female = All_Metabolomics_e[, null_label]
null_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c('null') & All_Metabolomics_p$Gender %in% "Male"]
null_e_male = All_Metabolomics_e[, null_label]


p_values_female = list()
fold_changes_female = list()
adjusted_p_values_female = list()

for(g in 2:length(unique_genes)){
  current_gene = unique_genes[g]
  current_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c(current_gene) & All_Metabolomics_p$Gender %in% "Female"]
  
  current_e = All_Metabolomics_e_no_mising[, current_label]
  
  p_values_female[[current_gene]] = fold_changes_female[[current_gene]] = adjusted_p_values_female[[current_gene]] = c()
  for(i in 1:nrow(All_Metabolomics_e_no_mising)){
    
    if((sum(is.na(current_e[i,])) == ncol(current_e))|(sum(is.na(null_e_female[i,])) == ncol(null_e_female))){
      p_values_female[[current_gene]][i] = fold_changes_female[[current_gene]][i] = NA
    }else{
      dta = data.table(y = scale(c(current_e[i,], null_e_female[i,])), group = rep(c("Gene","Control"), c(ncol(current_e), ncol(null_e_female))))
      colnames(dta) = c("y","group")
      
      p_values_female[[current_gene]][i] = summary(lm(y~group,data = dta))$coef[2,4]
      fold_changes_female[[current_gene]][i] = summary(lm(y~group,data = dta))$coef[2,1]
    }
  }
  
  adjusted_p_values_female[[current_gene]] = p.adjust(p_values_female[[current_gene]],'fdr')
  
  
  
}
sapply(p_values_female,function(x){sum(x<0.05,na.rm = TRUE)})

for(i in 1:length(names(p_values_female))){
  
  fwrite(data.table(label = All_Metabolomics_f$label, p_value = p_values_female[[i]], fold_change = fold_changes_female[[i]], adjusted_p_value = adjusted_p_values_female[[i]]),paste0("5 FEMALE for each genotype (",names(p_values_female)[i],"), check which metabolite are associated with gene.gender, gene, gender, or not significant at all.csv"))
  
  
}



p_values_male = list()
fold_changes_male = list()
adjusted_p_values_male = list()

for(g in 2:length(unique_genes)){
  current_gene = unique_genes[g]
  current_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c(current_gene) & All_Metabolomics_p$Gender %in% "Male"]
  
  current_e = All_Metabolomics_e_no_mising[, current_label]
  
  p_values_male[[current_gene]] = fold_changes_male[[current_gene]] = adjusted_p_values_male[[current_gene]] = c()
  for(i in 1:nrow(All_Metabolomics_e_no_mising)){
    
    if((sum(is.na(current_e[i,])) == ncol(current_e))|(sum(is.na(null_e_male[i,])) == ncol(null_e_male))){
      p_values_male[[current_gene]][i] = fold_changes_male[[current_gene]][i] = NA
    }else{
      dta = data.table(y = scale(c(current_e[i,], null_e_male[i,])), group = rep(c("Gene","Control"), c(ncol(current_e), ncol(null_e_male))))
      colnames(dta) = c("y","group")
      
      p_values_male[[current_gene]][i] = summary(lm(y~group,data = dta))$coef[2,4]
      fold_changes_male[[current_gene]][i] = summary(lm(y~group,data = dta))$coef[2,1]
    }
  }
  
  adjusted_p_values_male[[current_gene]] = p.adjust(p_values_male[[current_gene]],'fdr')
  
  
  
}
sapply(p_values_male,function(x){sum(x<0.05,na.rm = TRUE)})


for(i in 1:length(names(p_values_male))){
  
  fwrite(data.table(label = All_Metabolomics_f$label, p_value = p_values_male[[i]], fold_change = fold_changes_male[[i]], adjusted_p_value = adjusted_p_values_male[[i]]),paste0("5 MALE for each genotype (",names(p_values_male)[i],"), check which metabolite are associated with gene.gender, gene, gender, or not significant at all.csv"))
  
  
}











p_values = list()
# fold_changes = list()
adjusted_p_values = list()

for(g in 2:length(unique_genes)){
  current_gene = unique_genes[g]
  current_label = All_Metabolomics_p$label[All_Metabolomics_p$Genotype %in% c(current_gene)]
  
  current_e = All_Metabolomics_e_no_mising[, current_label]
  
  current_p = All_Metabolomics_p[All_Metabolomics_p$label %in% current_label,]
  null_p = All_Metabolomics_p[All_Metabolomics_p$Genotype %in% 'null',]
  
  p_values[[current_gene]] = adjusted_p_values[[current_gene]] = c()
  # fold_changes[[current_gene]] = c()
  for(i in 1:nrow(All_Metabolomics_e_no_mising)){
    
    if((sum(is.na(current_e[i,])) == ncol(current_e))|(sum(is.na(null_e[i,])) == ncol(null_e))){
      p_values[[current_gene]][i] = NA
      # fold_changes[[current_gene]][i] = NA
    }else{
      dta = data.table(y = scale(c(current_e[i,], null_e[i,])), group = rep(c("Gene","Control"), c(ncol(current_e), ncol(null_e))), sex=c(current_p$Gender,null_p$Gender))
      colnames(dta) = c("y","group",'sex')
      
      summary = summary(lm(y~group*sex,data = dta))$coef
      
      if(nrow(summary) == 4){
        p_values[[current_gene]][i] =summary[4,4]
      }else{
        p_values[[current_gene]][i]= NA
      }
      
      
      # fold_changes[[current_gene]][i] = summary(lm(y~group,data = dta))$coef[2,1]
    }
  }
  
  adjusted_p_values[[current_gene]] = p.adjust(p_values[[current_gene]],'fdr')
  
  
  
}
sapply(p_values,function(x){sum(x<0.05,na.rm = TRUE)})


for(i in 1:length(names(p_values_male))){
  
  fwrite(data.table(label = All_Metabolomics_f$label, p_value = p_values[[i]], fold_change = NA, adjusted_p_value = adjusted_p_values[[i]]),paste0("5 INTERACTION for each genotype (",names(p_values_male)[i],"), check which metabolite are associated with gene.gender, gene, gender, or not significant at all.csv"))
  
  
}





# 6 in null, the correlation between each phenotype and each metabolite.

Null_Metabolomics = wcmc::read_data("Null_Metabolomics.xlsx")
Null_Metabolomics_p = Null_Metabolomics$p
Null_Metabolomics_f = Null_Metabolomics$f
Null_Metabolomics_e = Null_Metabolomics$e_matrix


num_missing_male = apply(Null_Metabolomics_e,1,function(x){
  sum(is.na(x[Null_Metabolomics_p$Gender %in% "Male"]))
})
num_missing_female = apply(Null_Metabolomics_e,1,function(x){
  sum(is.na(x[Null_Metabolomics_p$Gender %in% "Female"]))
})
names(num_missing_male) = names(num_missing_female) = Null_Metabolomics_f$label
missing_index = ((num_missing_male/sum(Null_Metabolomics_p$Gender %in% "Male")) > 0.8)|((num_missing_female/sum(Null_Metabolomics_p$Gender %in% "Female")) > 0.8)
Null_Metabolomics_e = Null_Metabolomics_e[!missing_index,]
Null_Metabolomics_f = Null_Metabolomics_f[!missing_index,]




#### filter missing value > 80%
for(i in 1:nrow(Null_Metabolomics_e)){
  Null_Metabolomics_e[i,is.na(Null_Metabolomics_e[i,])] = 0.5 * min(Null_Metabolomics_e[i,], na.rm = TRUE)
}
rownames(Null_Metabolomics_e) = Null_Metabolomics_f$label




Raw_Phenotype_UCDavis = wcmc::read_data("Raw_Phenotype_UCDavis.xlsx")
Raw_Phenotype_UCDavis_p = Raw_Phenotype_UCDavis$p
Raw_Phenotype_UCDavis_f = Raw_Phenotype_UCDavis$f
Raw_Phenotype_UCDavis_e = Raw_Phenotype_UCDavis$e_cat_matrix



Raw_Phenotype_UCDavis_f$label_procedure = sapply(strsplit(Raw_Phenotype_UCDavis_f$label,"___"),function(x){x[1]})
Raw_Phenotype_UCDavis_f$label_phenotype = sapply(strsplit(Raw_Phenotype_UCDavis_f$label,"___"),function(x){x[2]})

rownames(Raw_Phenotype_UCDavis_e) = Raw_Phenotype_UCDavis_f$label

# combine the phenotypes because one phenotype may have multiple procedures.
unique_phenotype = unique(Raw_Phenotype_UCDavis_f$label_phenotype)
for(u in 1:length(unique_phenotype)){
  if(sum(Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u])>1){
    if(length(table(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],]))>1){
      duplicated_value_index = which(apply(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],],2,function(x){
        if(length(unique(x)) == 1){
          if(!unique(x)=="NA"){
            return(TRUE)
          }else{
            return(FALSE)
          }
        }else{
          return(FALSE)
        }
      }))
      if(length(duplicated_value_index)>0){
        print(unique_phenotype[u])
        print(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],1:10])
        Sys.sleep(20)
      }
    }
  }
}




Raw_Phenotype_UCDavis_e_merge = matrix(NA, nrow = length(unique(Raw_Phenotype_UCDavis_f$label_phenotype)), ncol = ncol(Raw_Phenotype_UCDavis_e))

colnames(Raw_Phenotype_UCDavis_e_merge) = colnames(Raw_Phenotype_UCDavis_e)

rownames(Raw_Phenotype_UCDavis_e_merge) = unique(Raw_Phenotype_UCDavis_f$label_phenotype)

i = 1
for(u in 1:length(unique_phenotype)){
  
  if(sum(Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u])>1){
    # stop(u)
    for(j in 1:ncol(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],])){
      the_j_th_col = Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],j]
      if(!length(unique(the_j_th_col)) == 1){
        Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],][the_j_th_col == "NA",j] = Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],][!the_j_th_col == "NA",j] 
      } 
    }
    Raw_Phenotype_UCDavis_e_merge[i,] = unique(Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],])
    i = i + 1
  }else{
    Raw_Phenotype_UCDavis_e_merge[i,] = Raw_Phenotype_UCDavis_e[Raw_Phenotype_UCDavis_f$label_phenotype %in% unique_phenotype[u],]
    i = i + 1
  }
}



# 
# for(i in 1:nrow(Raw_Phenotype_UCDavis_e_merge)){
#   
#   if(sum(is.na(as.numeric(Raw_Phenotype_UCDavis_e_merge[i,]))) > nrow(Raw_Phenotype_UCDavis_p)/2){
#     print(i)
#     print(table(Raw_Phenotype_UCDavis_e_merge[i,]))
#     cat("\n")
#     Sys.sleep(2)
#   }
#   
# }
# 
# 
Raw_Phenotype_UCDavis_e_merge[Raw_Phenotype_UCDavis_e_merge == "NA"] = NA


Raw_Phenotype_UCDavis_e_null = Raw_Phenotype_UCDavis_e_merge[,Raw_Phenotype_UCDavis_p$Genotype %in% "(null)"]
Raw_Phenotype_UCDavis_p_null = Raw_Phenotype_UCDavis_p[Raw_Phenotype_UCDavis_p$Genotype %in% "(null)",]
rownames(Raw_Phenotype_UCDavis_e_null) = rownames(Raw_Phenotype_UCDavis_e_merge)


Raw_Phenotype_UCDavis_e_null = Raw_Phenotype_UCDavis_e_null[,colnames(Null_Metabolomics_e)]
Raw_Phenotype_UCDavis_e_null = apply(Raw_Phenotype_UCDavis_e_null,2,as.numeric)
rownames(Raw_Phenotype_UCDavis_e_null) = rownames(Raw_Phenotype_UCDavis_e_merge)

# all gender
cor = cor(t(Null_Metabolomics_e), t(Raw_Phenotype_UCDavis_e_null),use = "pairwise.complete.obs", method = "spearman")

rownames(cor) = rownames(Null_Metabolomics_e)
colnames(cor) = rownames(Raw_Phenotype_UCDavis_e_null)

cor_test = cor
for(i in 1:nrow(cor_test)){
  print(i)
  for(j in 1:ncol(cor_test)){
    if(is.na(cor[i,j])){
      cor_test[i,j] = NA
    }else{
      cor_test[i,j] = cor.test(Null_Metabolomics_e[i,], Raw_Phenotype_UCDavis_e_null[j,], method = "spearman")$p.value
    }
  }
}



# female only
female_labels = Raw_Phenotype_UCDavis_p$label[Raw_Phenotype_UCDavis_p$Gender %in% "F"]
cor_female = cor(t(Null_Metabolomics_e[,colnames(Null_Metabolomics_e) %in% female_labels]), t(Raw_Phenotype_UCDavis_e_null[,colnames(Raw_Phenotype_UCDavis_e_null) %in% female_labels]),use = "pairwise.complete.obs", method = "spearman")

rownames(cor_female) = rownames(Null_Metabolomics_e[,colnames(Null_Metabolomics_e) %in% female_labels])
colnames(cor_female) = rownames(Raw_Phenotype_UCDavis_e_null[,colnames(Raw_Phenotype_UCDavis_e_null) %in% female_labels])

cor_test_female = cor_female
for(i in 1:nrow(cor_test_female)){
  print(i)
  for(j in 1:ncol(cor_test_female)){
    if(is.na(cor[i,j])){
      cor_test_female[i,j] = NA
    }else{
      cor_test_female[i,j] = cor.test(Null_Metabolomics_e[i,], Raw_Phenotype_UCDavis_e_null[j,], method = "spearman")$p.value
    }
  }
}
write.csv(cor_female,paste0("6 in null, the correlation between each phenotype and each metabolite (correlation female).csv"))
write.csv(cor_test_female,paste0("6 in null, the correlation between each phenotype and each metabolite (correlation p value female).csv"))


# male only
male_labels = Raw_Phenotype_UCDavis_p$label[Raw_Phenotype_UCDavis_p$Gender %in% "M"]
cor_male = cor(t(Null_Metabolomics_e[,colnames(Null_Metabolomics_e) %in% male_labels]), t(Raw_Phenotype_UCDavis_e_null[,colnames(Raw_Phenotype_UCDavis_e_null) %in% male_labels]),use = "pairwise.complete.obs", method = "spearman")

rownames(cor_male) = rownames(Null_Metabolomics_e[,colnames(Null_Metabolomics_e) %in% male_labels])
colnames(cor_male) = rownames(Raw_Phenotype_UCDavis_e_null[,colnames(Raw_Phenotype_UCDavis_e_null) %in% male_labels])

cor_test_male = cor_male
for(i in 1:nrow(cor_test_male)){
  print(i)
  for(j in 1:ncol(cor_test_male)){
    if(is.na(cor[i,j])){
      cor_test_male[i,j] = NA
    }else{
      cor_test_male[i,j] = cor.test(Null_Metabolomics_e[i,], Raw_Phenotype_UCDavis_e_null[j,], method = "spearman")$p.value
    }
  }
}

write.csv(cor_male,paste0("6 in null, the correlation between each phenotype and each metabolite (correlation male).csv"))
write.csv(cor_test_male,paste0("6 in null, the correlation between each phenotype and each metabolite (correlation p value male).csv"))























































