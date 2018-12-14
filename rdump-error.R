if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("bladderbatch", quietly = TRUE))
  BiocManager::install("bladderbatch")
library(bladderbatch)
library(tidyverse)
library(brms)

data(bladderdata, package = 'bladderbatch')
# Get the expression data
edata = exprs(bladderEset)
# Get the pheno data
pdata = pData(bladderEset)

# construct denormalized dataset for input to brms
expr_data <- edata %>%
  as.data.frame() %>%
  tibble::rownames_to_column('gene') %>%
  tidyr::gather(sample_desc, normalized_expression, -gene) %>%
  dplyr::left_join(pdata %>%
                     dplyr::mutate(sample_desc = colnames(edata)),
                   by = 'sample_desc')

# this data comprises expression for 57 samples
assertthat::assert_that(
  expr_data %>% dplyr::group_by(sample) %>%
    dplyr::n_groups()
  == 57
  )

# and 22,283 genes
assertthat::assert_that(
  expr_data %>% dplyr::group_by(gene) %>% 
    dplyr::n_groups()
  == 22283
)

# We will look at only the "Cancer" samples for simplicity
training_data <- expr_data %>% 
  dplyr::filter(cancer == 'Cancer') 

# among these samples, we have 3 outcomes
training_data %>%
  dplyr::distinct(outcome)

# Let's imagine we want to evaluate how gene expression differs in sTCC vs mTCC
training_data <- training_data %>%
  dplyr::mutate(is_mtcc = ifelse(stringr::str_detect(outcome, pattern = 'mTCC'), 1, 0))
training_data %>% dplyr::distinct(is_mtcc, outcome)

# we are left with data for 40 samples
assertthat::assert_that(
  training_data %>% dplyr::group_by(sample) %>% n_groups() 
  == 40)

# This is too many for us to work with when preparing data.
# let's sample 3 per bath
training_sample <- training_data %>%
  dplyr::group_by(batch) %>%
  dplyr::distinct(sample_desc) %>%
  dplyr::sample_n(3) %>%
  dplyr::ungroup()

# brms inputs (hypothetically)
brms_inputs <- list(formula = brms::bf(normalized_expression ~ 1 + is_mtcc + is_mtcc:gene + (1 | batch) + (1 | gene)),
                    data = training_data %>% dplyr::semi_join(training_sample, by = c('sample_desc', 'batch')),
                    sparse = T,
                    priors = brms::set_prior(horseshoe(df = 3, par_ratio = 0.1), class = 'b'),
                    family = hurdle_gamma)
brms_stan_code <- purrr::invoke(brms::make_stancode, brms_inputs)
brms_stan_data <- purrr::invoke(brms::make_standata, brms_inputs)
rstan::stan_rdump(names(brms_stan_data), file = 'test_rdump.Rds', envir = list2env(brms_stan_data))

