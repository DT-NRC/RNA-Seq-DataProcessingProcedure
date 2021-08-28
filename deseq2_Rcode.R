library("DESeq2")

genome="iwgsc_v1" 

dir_count_table=paste("PATH/TO/DATA", genome, ".txt", sep="")
dir_conditions=paste("PATH/TO/conditionFile", genome, ".txt", sep ="")
dir_lib_types=paste("PATH/TO/libTypeFile", genome, ".txt", sep="")
conditions=read.table(dir_conditions)
count_table=read.table(dir_count_table, header=TRUE, row.names=1 )
lib_types=read.table(dir_lib_types)
head(count_table)

# Wild Type samples and FHB samples
conditions=c("WT","WT","WT","WT","WT", "FHB","FHB","FHB","FHB","FHB", "ABA","ABA","ABA","ABA","ABA")

# treatment samples
expts=c("FHB","FHB","FHB","FHB","FHB", "ABA","ABA","ABA","ABA","ABA")

# control samples
ctrls=c("WT","WT","WT","WT","WT")

#1. FHB vs. mock 
expts[1]=c("FHB","FHB","FHB","FHB","FHB")
ctrls[1]=c("WT","WT","WT","WT","WT")

#2. ABA+FHB vs. mock
expts[2]=c("ABA","ABA","ABA","ABA","ABA") 
ctrls[2]=c("WT","WT","WT","WT","WT")

# more comparison pairs can be added ......

meta_data=data.frame(condition=conditions, libType=lib_types, row.names=colnames(count_table))

for (i in 1:length(expts))
{
ctrl=ctrls[i]
expt=expts[i]

dir_save_sorted=paste("PATH/TO/SaveOutputDATA", genome, "_", expt, "_VS_", ctrl, ".txt", sep="")
dir_save_unsorted=paste("PATH/TO/SaveOutputDATA", genome, "_", expt, "_vs_", ctrl, "_unsorted_test.txt", sep="")

#build meta-data
meta_data=data.frame(condition=conditions, libType=lib_types, row.names=colnames(count_table))
colnames(meta_data)=c("condition","libType")
meta_data

deseq_dataset=DESeqDataSetFromMatrix(countData = count_table, colData = meta_data, design = ~ condition)
deseq_dataset
deseq_dataset$condition = relevel(deseq_dataset$condition, ctrl)
deseq_dataset

print("Reads Counts before normalization:")
head(counts(deseq_dataset),n=20)

# call DESeq
deg=DESeq(deseq_dataset) 
# note: the counts in deg is unchanged after calling DESeq
print("Reads Counts after normalization:")
head(round(counts(deg,normalized=TRUE),4),n=20)

deg_results=results(deg, addMLE=TRUE, alpha=0.01,contrast=c("condition",expt,ctrl))
deg_results
mcols(deg_results)$description

# ranking
ind=order(deg_results$padj)
deg_results_ordered=data.frame(deg_results[ind,])
deg_results=data.frame(deg_results)
ncol_result=ncol(deg_results_ordered)
deg_results_unordered=data.frame( round(deg_results[,1:(ncol_result-2)],4), format(deg_results[,(ncol_result-1):ncol_result],digits=4, scientific=TRUE) ) 
deg_results_ordered=data.frame( round(deg_results_ordered[,1:(ncol_result-2)],4), format(deg_results_ordered[,(ncol_result-1):ncol_result],digits=4, scientific=TRUE) ) 
count_table_original_unordered=counts(deg,normalized=FALSE)
count_table_original_ordered=counts(deg,normalized=FALSE)[ind,]
count_table_normalized_unordered=round(counts(deg,normalized=TRUE),4)
count_table_normalized_ordered=round(counts(deg,normalized=TRUE),4)[ind,]

summary(deg_results)
count_table_deg_result_unordered=data.frame(count_table_original_unordered, count_table_normalized_unordered, deg_results_unordered)
head(count_table_deg_result_unordered,n=20)

# save results
write.table(count_table_deg_result_unordered, file=dir_save_unsorted, sep="\t", row.names=TRUE, col.name=NA, quote=FALSE)
write.table(count_table_deg_result_ordered, file=dir_save_sorted, sep="\t", row.names=TRUE, col.name=NA, quote=FALSE)
} # end of for
