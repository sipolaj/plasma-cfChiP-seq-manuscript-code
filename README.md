# plasma-cfChiP-seq-manuscript-code
This repository contains the code used to carry out the analyses in our article "Plasma cell-free DNA histone methylation enables phenotypic and clinical segmentation of metastatic prostate cancer".

### Scaling fragment counts between samples and correcting GC bias differences
Two Julia functions, `normalise_counts(counts, gc_frac)` and `normalise_counts_with_control_peaks(counts, gc_frac)`, in the file `data_normalisation.jl` were used to normalise count matrices. 

`normalise_counts(counts, gc_frac)` takes in a fragment count matrix `counts`, where columns are samples and rows are genomic features (e.g. transcription start sites), and a GC fraction vector `gc_frac`, which includes GC fractions for each genomic feature. The function scales the samples so that their means are equal, while simultaneously correcting for differences in their GC biases.

`normalise_counts_with_control_peaks(counts, gc_frac)` takes in a fragment count matrix `counts`, where columns are samples and rows are 100bp genomic windows, and a GC fraction vector `gc_frac`, which includes GC fractions for each genomic window. The function first calls 20,000 consensus peaks, then selects 2000 peaks with the least between-sample variance, and finally scales each sample so that the mean of the 2000 peaks is 1.0, while simultaneously correcting for their GC biases.


### Fitting linear models to count data
Below is an example Julia code of a how we fit linear models to fragment count and local ctDNA fraction data.

Load normalised fragment counts and local ctDNA fractions from Supplementary Table 8.
```
using DelimitedFiles, Statistics, GLM, DataFrames, CategoricalArrays
d = readdlm("supplementary_table_8.tsv", '\t')
samples = d[1, 2:end]
genes = d[2:end, 1]
d = d[2:end, 2:end]
counts = [parse(Float64, split(x, ";")[1]) for x=d]
cancer_frac = [parse(Float64, split(x, ";")[2]) for x=d]
```

Extract counts and local ctDNA fractions of a single gene, HOXB13
```
gene_idx = only(findall(genes .== "HOXB13"))
gene_counts = counts[gene_idx, :]
gene_cancer_frac = cancer_frac[gene_idx, :]
```

Calculate correlation between counts and cancer fraction
```
corr_cancer_frac = cor(gene_counts, gene_cancer_frac)
```

Fit a simple linear model `β0 + β1F`, where `F` is cancer fraction
```
data = DataFrame(y=gene_counts, x1=gene_cancer_frac)
model0 = GLM.lm(@formula(y ~ x1), data)
β0, β1 = coef(model0)
normal = β0 # Model count for sample with 0% cancer
cancer = β0 + β1 # Model count for sample with 100% cancer
```

Fit a bivariate linear model `β0 + β1F + β2FG`, where `F` is cancer fraction and `G` is sample group
```
group = zeros(length(gene_counts))
group[1:5] .= 1 # First 5 samples are assigned to group 1, the rest to group 0
group = categorical(group)
data = DataFrame(y=gene_counts, x1=gene_cancer_frac, x2=group)
model1 = GLM.lm(@formula(y ~ x1 + x1 & x2), data)
β0, β1, β2 = coef(model0)
normal = β0 # Model count for sample with 0% cancer
cancer0 = β0 + β1 # Model count for group 0 sample with 100% cancer
cancer1 = β0 + β1 +β2 # Model count for group 1 sample with 100% cancer
```

Calculate F-test p-value for the group term `β2FG` in the bivariate model
```
f = ftest(model0.model, model1.model)
p = f.pval[2]
```
