Illustrate error with stan\_rdump
================
Jacqueline Buros
12/14/2018

Review data
-----------

Let's load some example data from the `bladderbatch` package.

``` r
data(bladderdata, package = 'bladderbatch')
# Get the expression data
edata = exprs(bladderEset)
# Get the pheno data
pdata = pData(bladderEset)
```

Next we convert to a long format for input to brms

``` r
expr_data <- edata %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gene') %>%
  tidyr::gather(sample_desc, normalized_expression, -gene) %>%
  dplyr::left_join(pdata %>%
                     dplyr::mutate(sample_desc = colnames(edata)),
                   by = 'sample_desc')
```

Note that this dataset comprises expression for 57 samples

``` r
assertthat::assert_that(
  expr_data %>% dplyr::group_by(sample) %>%
    dplyr::n_groups()
  == 57
  )
```

    ## [1] TRUE

and 22,283 genes

``` r
assertthat::assert_that(
  expr_data %>% dplyr::group_by(gene) %>% 
    dplyr::n_groups()
  == 22283
)
```

    ## [1] TRUE

In this example we will look at only the "Cancer" samples for simplicity

``` r
training_data <- expr_data %>% 
  dplyr::filter(cancer == 'Cancer') 
```

among these samples, we have 3 outcomes

``` r
training_data %>%
  dplyr::distinct(outcome)
```

    ##    outcome
    ## 1 sTCC+CIS
    ## 2 sTCC-CIS
    ## 3     mTCC

Let's imagine we want to evaluate how gene expression differs in sTCC vs mTCC

``` r
training_data <- training_data %>%
  dplyr::mutate(is_mtcc = ifelse(stringr::str_detect(outcome, pattern = 'mTCC'), 1, 0))
training_data %>% dplyr::distinct(is_mtcc, outcome)
```

    ##    outcome is_mtcc
    ## 1 sTCC+CIS       0
    ## 2 sTCC-CIS       0
    ## 3     mTCC       1

We are left with data for 40 samples. This is too many samples when fitting with brms, so we select a subset for illustration.

``` r
assertthat::assert_that(
  training_data %>% dplyr::group_by(sample) %>% n_groups() 
  == 40)
```

    ## [1] TRUE

``` r
training_sample <- training_data %>%
  dplyr::group_by(batch) %>%
  dplyr::distinct(sample_desc) %>%
  dplyr::sample_n(3) %>%
  dplyr::ungroup()
```

Prepare data for `stan_rdump`
-----------------------------

``` r
# brms inputs (hypothetically)
brms_inputs <- list(formula = brms::bf(normalized_expression ~ 1 + is_mtcc + is_mtcc:gene + (1 | batch) + (1 | gene)),
                    data = training_data %>% dplyr::semi_join(training_sample, by = c('sample_desc', 'batch')),
                    sparse = T,
                    priors = brms::set_prior(horseshoe(df = 3, par_ratio = 0.1), class = 'b'),
                    family = hurdle_gamma)
brms_stan_code <- purrr::invoke(brms::make_stancode, brms_inputs)
brms_stan_data <- purrr::invoke(brms::make_standata, brms_inputs)
```

Write out to file using `stan_rdump`

``` r
## errors
try(
  rstan::stan_rdump(names(brms_stan_data), file = 'test_rdump.Rds', envir = list2env(brms_stan_data))
  )
```

Session info
------------

``` r
sessionInfo()
```

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Debian GNU/Linux 9 (stretch)
    ## 
    ## Matrix products: default
    ## BLAS: /usr/lib/openblas-base/libblas.so.3
    ## LAPACK: /usr/lib/libopenblasp-r0.2.19.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] bindrcpp_0.2.2      brms_2.6.0          Rcpp_1.0.0         
    ##  [4] forcats_0.3.0       stringr_1.3.1       dplyr_0.7.8        
    ##  [7] purrr_0.2.5         readr_1.3.0         tidyr_0.8.2        
    ## [10] tibble_1.4.2        ggplot2_3.1.0       tidyverse_1.2.1    
    ## [13] bladderbatch_1.20.0 Biobase_2.42.0      BiocGenerics_0.28.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-137         matrixStats_0.54.0   xts_0.11-2          
    ##  [4] lubridate_1.7.4      threejs_0.3.1        httr_1.4.0          
    ##  [7] rstan_2.18.2         tools_3.5.1          backports_1.1.2     
    ## [10] R6_2.3.0             DT_0.5               lazyeval_0.2.1      
    ## [13] colorspace_1.3-2     withr_2.1.2          tidyselect_0.2.5    
    ## [16] gridExtra_2.3        prettyunits_1.0.2    processx_3.2.0      
    ## [19] Brobdingnag_1.2-6    compiler_3.5.1       cli_1.0.1           
    ## [22] rvest_0.3.2          shinyjs_1.0          xml2_1.2.0          
    ## [25] colourpicker_1.0     scales_1.0.0         dygraphs_1.1.1.6    
    ## [28] mvtnorm_1.0-8        ggridges_0.5.1       callr_3.0.0         
    ## [31] digest_0.6.18        StanHeaders_2.18.0   rmarkdown_1.11      
    ## [34] base64enc_0.1-3      pkgconfig_2.0.2      htmltools_0.3.6     
    ## [37] htmlwidgets_1.3      rlang_0.3.0.1        readxl_1.1.0        
    ## [40] rstudioapi_0.8       shiny_1.2.0          bindr_0.1.1         
    ## [43] generics_0.0.2       zoo_1.8-4            jsonlite_1.6        
    ## [46] gtools_3.8.1         crosstalk_1.0.0      inline_0.3.15       
    ## [49] magrittr_1.5         loo_2.0.0            bayesplot_1.6.0     
    ## [52] Matrix_1.2-14        munsell_0.5.0        abind_1.4-5         
    ## [55] yaml_2.2.0           stringi_1.2.4        pkgbuild_1.0.2      
    ## [58] plyr_1.8.4           grid_3.5.1           promises_1.0.1      
    ## [61] crayon_1.3.4         miniUI_0.1.1.1       lattice_0.20-35     
    ## [64] haven_2.0.0          hms_0.4.2            knitr_1.21          
    ## [67] ps_1.2.1             pillar_1.3.0         igraph_1.2.2        
    ## [70] markdown_0.9         shinystan_2.5.0      reshape2_1.4.3      
    ## [73] stats4_3.5.1         rstantools_1.5.1     glue_1.3.0          
    ## [76] evaluate_0.12        BiocManager_1.30.4   modelr_0.1.2        
    ## [79] httpuv_1.4.5         cellranger_1.1.0     gtable_0.2.0        
    ## [82] assertthat_0.2.0     xfun_0.4             mime_0.6            
    ## [85] xtable_1.8-3         broom_0.5.1          coda_0.19-2         
    ## [88] later_0.7.5          rsconnect_0.8.12     shinythemes_1.1.2   
    ## [91] bridgesampling_0.6-0
