# 6 in null, the correlation phenotype metabolite (correlation male).csv

pacman::p_load(data.table, corrplot)
cor_data= fread("6 in null, the correlation phenotype metabolite (correlation male).csv")
cor = data.matrix(cor_data[,-1,with=FALSE])





cor_p_data= fread("6 in null, the correlation phenotype metabolite (correlation p value male).csv")
cor_p = data.matrix(cor_p_data[,-1,with=FALSE])

met_sig = apply(cor_p,1,function(x){sum(x<0.05, na.rm = TRUE)})
phe_sig = apply(cor_p,2,function(x){sum(x<0.05, na.rm = TRUE)})





cor = cor[order(met_sig, decreasing = TRUE), order(phe_sig, decreasing = TRUE)]
cor_p = cor_p[order(met_sig, decreasing = TRUE), order(phe_sig, decreasing = TRUE)]
met_sig = met_sig[order(met_sig, decreasing = TRUE)]
phe_sig = phe_sig[order(phe_sig, decreasing = TRUE)]




cor = cor[!met_sig == 0,!phe_sig == 0]
cor_p = cor_p[!met_sig == 0,!phe_sig == 0]
met_sig = met_sig[!met_sig==0]


cor = cor[met_sig > 10,]
cor_p = cor_p[met_sig > 10,]


met_order = hclust(dist(cor))$order
phe_order = hclust(dist(t(cor)))$order


cor = cor[met_order,phe_order]
cor_p = cor_p[met_order,phe_order]

colnames(cor) = NULL


corrplot(cor, p.mat = cor_p, method="circle", insig = "blank")

