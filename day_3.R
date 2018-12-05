




mat = read.table("ALLcancerdata.txt", check.names = FALSE)
mat = as.matrix(mat)
anno = read.table("ALLannotation.txt")

head(mat)
head(anno)

## pure numbers as sample names is very bad habit
head(mat[, 39, drop = FALSE])
head(mat[, "39", drop = FALSE])

colnames(mat) = paste0("sample_", colnames(mat))
rownames(anno) = paste0("sample_", anno$Samples)
anno = anno[, "ALL.AML", drop = FALSE]

# double check the consistency of samples in the two objects
identical(colnames(mat), rownames(anno))

quantile(mat)

mat[mat < 100] = 100
mat[mat > 16000] = 16000



l = apply(mat, 1, function(x) {
	max(x)/min(x) <= 5 || max(x) - min(x) <= 500
})

mat[!l, ]


mat = mat[!l, ]

quantile(mat)



gene_max = apply(mat, 1, max)
gene_min = apply(mat, 1, min)

# max/min <= 5 or max - min <= 500
l = gene_max/gene_min <= 5 | gene_max - gene_min <= 500


mat = log10(mat)


## PCA analysis
fit = prcomp(t(mat))  # why do we need to transpose the matrix?
loc = fit$x[, 1:2]  # the first two PCs
plot(loc, col = ifelse(anno[, 1] == "AML", "red", "blue"))


# proportion of variance in each PC
all_var = apply(fit$x, 2, var)
p = all_var/sum(all_var)


# MDS analysis
loc = cmdscale(dist(t(mat)), k = 2)
plot(loc, col = ifelse(anno[, 1] == "AML", "red", "blue"))

# what if we change to a different distance measurement
cor_mat = cor(mat)
cor_dist = as.dist(1 - cor_mat)
loc = cmdscale(cor_dist, k = 2)
plot(loc, col = ifelse(anno[, 1] == "AML", "red", "blue"))

## kmeans
set.seed(123)
km = kmeans(t(mat), centers = 2)
km_cluster = km$cluster
table(km_cluster, anno[, 1])

# repeat kmeans multiple times


# perfrom k-means with random centers 100 times
all_km = NULL
for(i in 1:100) {
	km = kmeans(t(mat), centers = 2)$cluster
	group_mean = tapply(colMeans(mat), km, mean)
	if(group_mean[1] > group_mean[2]) {
		km = c(2, 1)[km]
	}
	all_km = rbind(all_km, km)
}

library(ComplexHeatmap)
# visualize how the 100 k-means clusters look like
Heatmap(all_km,
	top_annotation = HeatmapAnnotation(df = anno),
	row_title = "100 k-means")

# calculate the consistency of clustering from 100 k-means
freq = apply(all_km, 2, function(x) {
	c("1" = sum(x == 1),
	  "2" = sum(x == 2))
})

# the proportion or the probability to belong to each cluster
prop = apply(freq, 2, function(x) x/sum(x))

# add the probability to heatmap
Heatmap(all_km,
	top_annotation = HeatmapAnnotation(df = anno,
		prop = t(prop)),
	row_title = "100 k-means")

## heatmap
Heatmap(mat, 
	top_annotation = HeatmapAnnotation(df = anno), 
	show_row_names = FALSE,
	show_column_names = FALSE)

## top 1000 genes with highest variance
gene_var = apply(mat, 1, var)
top_genes = order(gene_var, decreasing = TRUE)[1:1000]

mat_top = mat[top_genes, ]

fit = prcomp(t(mat_top))  # why do we need to transpose the matrix?
loc = fit$x[, 1:2]  # the first two PCs
plot(loc, col = ifelse(anno[, 1] == "AML", "red", "blue"))

Heatmap(mat_top, 
	top_annotation = HeatmapAnnotation(df = anno), 
	show_row_names = FALSE,
	show_column_names = FALSE, 
	column_title = "top 1000 genes with higest variance")

