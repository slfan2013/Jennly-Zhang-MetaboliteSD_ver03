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

p_val3 = c()
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
    
    
  }else{
    continuous_index[i] = FALSE
    
    
    p_val3[length(p_val3)+1] = summary(xtabs(~Raw_Phenotype_UCDavis_e_null[i,] + Raw_Phenotype_UCDavis_p_null$Gender))$p.value
    
  }
  
}



sum(p_val3[continuous_index]<0.05,na.rm = TRUE)/sum(continuous_index)

## 4 phenotype -  sex adjusted by body weight
# body_weight_all
# sum(Raw_Phenotype_UCDavis_p_null$label%in%names(body_weight_all))/length(Raw_Phenotype_UCDavis_p_null$label)

body_weight = body_weight_all[Raw_Phenotype_UCDavis_p_null$label]


p_val4 = c()
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
    
    
  }else{
    continuous_index[i] = FALSE
    
    
    p_val4[length(p_val4)+1] = summary(xtabs(~Raw_Phenotype_UCDavis_e_null[i,] + Raw_Phenotype_UCDavis_p_null$Gender + body_weight))$p.value
    
  }
  
}

sum(p_val4[continuous_index]<0.05,na.rm = TRUE)/sum(continuous_index)



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


cor = cor(t(Null_Metabolomics_e), t(Raw_Phenotype_UCDavis_e_null),use = "pairwise.complete.obs")

rownames(cor) = rownames(Null_Metabolomics_e)
colnames(cor) = rownames(Raw_Phenotype_UCDavis_e_null)


