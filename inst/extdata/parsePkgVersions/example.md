---
title: "Example for parsePkgVersions()"
output: 
  html_document:
    theme: sandstone
    toc: true
    toc_float: true
editor_options: 
  chunk_output_type: console
---



# Load packages


```r
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(umap)
  library(DropletUtils)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(swissknife)
  library(stringr)
  library(sctransform)
  library(BiocParallel)
})
```


# Session info


```r
date()
```

```
## [1] "Fri Aug 16 17:35:03 2019"
```

```r
sessionInfo()
```

```
## R version 3.6.0 (2019-04-26)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /usr/lib64/libopenblas-r0.3.3.so
## 
## locale:
## [1] C
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] sctransform_0.2.0           stringr_1.4.0              
##  [3] swissknife_0.12             tidyr_0.8.3                
##  [5] dplyr_0.8.3                 tibble_2.1.3               
##  [7] DropletUtils_1.4.2          umap_0.2.2.0               
##  [9] scran_1.12.1                scater_1.12.2              
## [11] ggplot2_3.2.0               SingleCellExperiment_1.6.0 
## [13] SummarizedExperiment_1.14.0 DelayedArray_0.10.0        
## [15] BiocParallel_1.18.0         matrixStats_0.54.0         
## [17] Biobase_2.44.0              GenomicRanges_1.36.0       
## [19] GenomeInfoDb_1.20.0         IRanges_2.18.1             
## [21] S4Vectors_0.22.0            BiocGenerics_0.30.0        
## [23] RColorBrewer_1.1-2         
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6             RcppAnnoy_0.0.12        
##  [3] dynamicTreeCut_1.63-1    tools_3.6.0             
##  [5] R6_2.4.0                 irlba_2.3.3             
##  [7] KernSmooth_2.23-15       HDF5Array_1.12.1        
##  [9] vipor_0.4.5              uwot_0.1.3              
## [11] lazyeval_0.2.2           colorspace_1.4-1        
## [13] withr_2.1.2              tidyselect_0.2.5        
## [15] gridExtra_2.3            compiler_3.6.0          
## [17] BiocNeighbors_1.2.0      labeling_0.3            
## [19] scales_1.0.0             digest_0.6.20           
## [21] Rsamtools_2.0.0          rmarkdown_1.14          
## [23] R.utils_2.9.0            XVector_0.24.0          
## [25] pkgconfig_2.0.2          htmltools_0.3.6         
## [27] limma_3.40.4             rlang_0.4.0             
## [29] FNN_1.1.3                DelayedMatrixStats_1.6.0
## [31] jsonlite_1.6             R.oo_1.22.0             
## [33] RCurl_1.95-4.12          magrittr_1.5            
## [35] BiocSingular_1.0.0       GenomeInfoDbData_1.2.1  
## [37] Matrix_1.2-17            Rcpp_1.0.1              
## [39] ggbeeswarm_0.6.0         munsell_0.5.0           
## [41] Rhdf5lib_1.6.0           reticulate_1.12         
## [43] viridis_0.5.1            R.methodsS3_1.7.1       
## [45] stringi_1.4.3            yaml_2.2.0              
## [47] edgeR_3.26.5             MASS_7.3-51.4           
## [49] zlibbioc_1.30.0          Rtsne_0.15              
## [51] rhdf5_2.28.0             plyr_1.8.4              
## [53] grid_3.6.0               listenv_0.7.0           
## [55] dqrng_0.2.1              crayon_1.3.4            
## [57] lattice_0.20-38          Biostrings_2.52.0       
## [59] locfit_1.5-9.1           knitr_1.23              
## [61] pillar_1.4.2             igraph_1.2.4.1          
## [63] codetools_0.2-16         future.apply_1.3.0      
## [65] reshape2_1.4.3           glue_1.3.1              
## [67] evaluate_0.14            RcppParallel_4.4.3      
## [69] png_0.1-7                gtable_0.3.0            
## [71] purrr_0.3.2              future_1.14.0           
## [73] assertthat_0.2.1         xfun_0.8                
## [75] rsvd_1.0.1               RSpectra_0.15-0         
## [77] viridisLite_0.3.0        GenomicAlignments_1.20.1
## [79] beeswarm_0.2.3           globals_0.12.4          
## [81] statmod_1.4.32
```
