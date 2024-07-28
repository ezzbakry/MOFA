if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

install.packages("TCGAbiolinks")

#Start

library("limma")
library("vegan")
library("cluster")
library("factoextra")
library("gridExtra")
library("PerformanceAnalytics")
library("corrplot")
library("Hmisc")
library("RColorBrewer")
library("impute")
library("pathview")
library("glmnet")
library("reshape2")
library("readr")
library("magrittr")
library("ggrepel")
library("tidyverse")
library("ggfortify")
library("ggplot2")
library("ggrepel")
library("EnhancedVolcano")
library("TCGAbiolinks")
library("MOFA2")
library("MOFAdata")
library("data.table")
library("magrittr") 
library("dplyr") 
library("GGally")
library("reticulate")
reticulate::py_config()
mirna <- read.csv("D:\\مشروع_التخرج\\data\\Normalized mirna_data.csv")
colnames(mirna)
row.names(mirna) = mirna$X
mirna <- mirna[-1]
mirna <- as.data.frame(t(mirna))

rna <- read.csv("D:\\مشروع_التخرج\\data\\Normalized rnaseq_data.csv")
row.names(rna) = rna$X
rna <- rna[-1]
rna <- as.data.frame(t(rna))


methyl <- read.csv("D:\\مشروع_التخرج\\data\\Normalized methylation_data.csv")
row.names(methyl) = methyl$X
methyl <- methyl[-1]
methyl <- as.data.frame(t(methyl))

mutation <- read.csv("D:\\مشروع_التخرج\\data\\Normalized mutation_data.csv")
row.names(mutation) = mutation$X
mutation <- mutation[-1]
mutation <- as.data.frame(t(mutation))


#bring all samples names(columns) common between all omics types
common_cols=intersect(colnames(mirna),colnames(rna))
common_cols=intersect(common_cols,colnames(methyl))
common_cols=intersect(common_cols,colnames(mutation))

mirna_upd=mirna[,common_cols]
rna_upd=rna[,common_cols]
methyl_upd=methyl[,common_cols]
mutation_upd=mutation[,common_cols]

mirna_mat=as.matrix(mirna_upd)
rna_mat=as.matrix(rna_upd)
methyl_mat=as.matrix(methyl_upd)
mutation_mat=as.matrix(mutation_upd)


clinical_data=read.table("D:\\مشروع_التخرج\\data\\Human__TCGA_LAML__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi", header=TRUE, na.strings=c("NA","NaN", ""),row.names=1, sep="\t")
class(clinical_data)
#Filter it on common samples (columns)
dim(clinical_data)
clinical_data=t(clinical_data)
clinical_data=data.frame(clinical_data)
clinical_data$race <- replace(clinical_data$race, is.na(clinical_data$race), "white")
clinical_data <- subset(clinical_data, select = -overallsurvival)

dim(clinical_data)
clinical_data=t(clinical_data)
clinical_data1=clinical_data[,common_cols]
clinical_data1_t=data.frame(t(clinical_data1))
dim(clinical_data1_t)
clinical_data1_t$sample=rownames(clinical_data1_t)
dim(clinical_data1_t)


data=list(mirna_mat,rna_mat,methyl_mat,mutation_mat)
lapply(data, dim) 
lapply(data, class) 
MOFAobject<-create_mofa(data)

views_names(MOFAobject)=c('mirna','rna_seq','methylation','mutation')
MOFAobject
plot_data_overview(MOFAobject)


data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE
head(data_opts,1)

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 7
head(model_opts,1)

train_opts <- get_default_training_options(MOFAobject)
head(train_opts,1)

Nsamples = sum(MOFAobject@dimensions$N)


# samples_metadata(MOFAobject) <- clinical_data1_t
head(MOFAobject@samples_metadata, n=3)

MOFAobject <- prepare_mofa(object = MOFAobject ,
                           data_options = data_opts , 
                           model_options = model_opts , 
                           training_options = train_opts
                           )


outfile = paste0("D:\\model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject,outfile, use_basilisk=FALSE)
model <- load_model(outfile)
plot_data_overview(model)

head(model@samples_metadata, n=3)

head(model@cache$variance_explained$r2_total[[1]]) 
head(model@cache$variance_explained$r2_per_factor[[1]]) 
plot_variance_explained(model, x="view", y="factor")
plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]

Nsamples = sum(model@dimensions$N)
#We set sample metadata to the clinical data
sample_metadata=clinical_data1_t
samples_metadata(model) <- sample_metadata
dim(samples_metadata(model))
head(sample_metadata)
#column_to_rownames(samples_metadata)
head(model@samples_metadata, n=4)
head(model@cache$variance_explained$r2_total[[1]])
print(model)

# gender ethnicity  status   race  

#mirna     rna_seq     methylation


plot_factor(model, 
            factor = 1:4,
            color_by  = "gender")
plot_factors(model, 
            factor = 1:4,
            color_by  = "race")

plot_factor(model, 
                 factors = c(1,2,3,4),
                 color_by  = "status",
                 dot_size = 3,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = F,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)


# gender ethnicity  status   race  

plot_data_heatmap(model,
                  view = "mutation",         # view of interest
                  factor = 1,             # factor of interest
                  features = 20,          # number of features to plot (they are selected by weight)
                  cluster_rows = TRUE, cluster_cols = TRUE,
show_rownames = TRUE, show_colnames = TRUE,annotation_samples = "race")

plot_data_scatter(model,
                  view = "mutation",         # view of interest
                  factor = 1,             # factor of interest
                  features = 5,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,          # add linear regression
                  color_by = "race"
)





plot_weights(model,
             view = "mirna",
             factor = 1,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_weights(model,
             view = "rna_seq",
             factor = 1,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)
plot_weights(model,
             view = "methylation",
             factor = 1,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_weights(model,
             view = "mutation",
             factor = 1,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

###############################################
plot_top_weights(model,
                 view = "mirna",
                 factor = 1,
                 nfeatures = 10
)
