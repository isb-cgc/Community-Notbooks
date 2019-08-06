One-Way ANOVA Test in BigQuery
================

In this notebook, we will cover how to do a one-way ANOVA test in
BigQuery. This statistical test can be used to determine whether there
is a statistically significant difference between the means of two or
more independent groups. Although in this example, I’m only looking at
two groups, it would not be difficult to extend this to any number of
groups, assuming there is a reasonable number of samples within each
group.

Consider the model \(y_{ij} = m + a_i + e_{ij}\), \(y_{ij}\) is a
continuous variable over samples \(j\), in groups \(i\), and \(a_i\) is
a constant for each group \(i\), and \(e_{ij}\) is a gaussian error term
with mean 0.

Using this model, we are describing the data as being sampled from
groups, with each group having a mean value equal to \(m + a_i\). The
null hypothesis is that each of the group means is the same (ie that the
ai terms are zero), while the alternative hypothesis is that at least
one of the ai terms is not zero.

We use the F-test to compare these two hypotheses. To compute the test
statistic, we compute the within-group variation and the between-group
variation. Recall that sample variance is defined as the sum of squared
differences between observations and the mean, divided by the number of
samples (normalized).

First, we will need to import the required libraries and create a client
variable. For more information see [‘Quick Start Guide to
ISB-CGC’](INSERT%20LINK) and alternative authentication methods can be
found
[here](https://googleapis.github.io/google-cloud-python/latest/core/auth.html).

``` r
library(bigrquery)
library(dplyr)
```

Then let us set up the billing variables we will be using in this
notebook:

``` r
billing <- 'isb-cgc-02-0001' # Insert your project ID in the ''
```

Let’s look at the query:

``` r
anova_query = "
WITH
  -- using standard SQL,
  -- we'll select our cohort and gene expression
  --
  cohortExpr AS (
  SELECT
    sample_barcode,
    LOG10(normalized_count) AS expr
  FROM
    `isb-cgc.TCGA_hg19_data_v0.RNAseq_Gene_Expression_UNC_RSEM`
  WHERE
    project_short_name = 'TCGA-BRCA'
    AND HGNC_gene_symbol = 'TP53'
    AND normalized_count IS NOT NULL
    AND normalized_count > 0),
  --
  -- And we'll select the variant data for our cohort,
  -- we're going to be comparing variant types (SNP, DEL, etc)
  --
  cohortVar AS (
  SELECT
    Variant_Type,
    sample_barcode_tumor AS sample_barcode
  FROM
    `isb-cgc.TCGA_hg19_data_v0.Somatic_Mutation_MC3`
  WHERE
    SYMBOL = 'TP53' ),
  --
  -- then we join the expression and variant data using sample barcodes
  --
  cohort AS (
  SELECT
    cohortExpr.sample_barcode AS sample_barcode,
    Variant_Type AS group_name,
    expr
  FROM
    cohortExpr
  JOIN
    cohortVar
  ON
    cohortExpr.sample_barcode = cohortVar.sample_barcode ),
  --
  -- First part of the calculation, the grand mean (over everything)
  --
  grandMeanTable AS (
  SELECT
    AVG(expr) AS grand_mean
  FROM
    cohort ),
  --
  -- Then we need a mean per group, and we can get a count of samples
  -- per group.
  --
  groupMeansTable AS (
  SELECT
    AVG(expr) AS group_mean,
    group_name,
    COUNT(sample_barcode) AS n
  FROM
    cohort
  GROUP BY
    group_name),
  --
  -- To get the between-group variance
  -- we take the difference between the grand mean
  -- and the means for each group and sum over all samples
  -- ... a short cut being taking the product with n.
  -- Later we'll sum over the n_sq_diff
  --
  ssBetween AS (
  SELECT
    group_name,
    group_mean,
    grand_mean,
    n,
    n*POW(group_mean - grand_mean,2) AS n_diff_sq
  FROM
    groupMeansTable
  CROSS JOIN
    grandMeanTable ),
  --
  -- Then, to get the variance within each group
  -- we have to build a table matching up the group mean
  -- with the values for each group. So we join the group
  -- means to the values on group name. We are going to
  -- sum over this table just like ssBetween
  --
  ssWithin AS (
  SELECT
    a.group_name AS group_name,
    expr,
    group_mean,
    b.n AS n,
    POW(expr - group_mean, 2) AS s2
  FROM
    cohort a
  JOIN
    ssBetween b
  ON
    a.group_name = b.group_name ),
  --
  -- The F stat comes from a ratio, the numerator is
  -- calculated using the between group variance, and
  -- dividing by the number of groups (k) minus 1.
  --
  numerator AS (
  SELECT
    'dummy' AS dummy,
    SUM(n_diff_sq) / (count(group_name) - 1) AS mean_sq_between
  FROM
    ssBetween ),
  --
  -- The denominator of the F stat ratio is found using the
  -- variance within groups. We divide the sum of the within
  -- group variance and divide it by (n-k).
  --
  denominator AS (
  SELECT
    'dummy' AS dummy,
    COUNT(distinct(group_name)) AS k,
    COUNT(group_name) AS n,
    SUM(s2)/(COUNT(group_name)-COUNT(distinct(group_name))) AS mean_sq_within
  FROM
    ssWithin),
  --
  -- Now we're ready to calculate F!
  --
  Ftable AS (
  SELECT
    n,
    k,
    mean_sq_between,
    mean_sq_within,
    mean_sq_between / mean_sq_within AS F
  FROM
    numerator
  JOIN
    denominator
  ON
    numerator.dummy = denominator.dummy)

SELECT
  *
FROM
  Ftable
"
```

``` r
# To see the R console output with query processing information, turn queit to FALSE
anova_result <- bq_project_query(billing, anova_query, quiet = TRUE)
# Transform the query result into a tibble
anova_result <- bq_table_download(anova_result, quiet = TRUE)
anova_result
```

    ## # A tibble: 1 x 5
    ##       n     k mean_sq_between mean_sq_within     F
    ##   <int> <int>           <dbl>          <dbl> <dbl>
    ## 1   375     3            6.23         0.0956  65.1

OK, so let’s check our work. Using the BRCA cohort and TP53 as our gene,
we have 375 samples with a variant in this gene. We’re going to look at
whether the type of variant is related to the gene expression we
observe. If we just pull down the data using the ‘cohort’ subtable (as
above), we can get a small data frame, which let’s us do the standard F
stat table in R.

``` r
dat_query = "
WITH
  -- using standard SQL,
  -- we'll select our cohort and gene expression
  --
  cohortExpr AS (
  SELECT
    sample_barcode,
    LOG10(normalized_count) AS expr
  FROM
    `isb-cgc.TCGA_hg19_data_v0.RNAseq_Gene_Expression_UNC_RSEM`
  WHERE
    project_short_name = 'TCGA-BRCA'
    AND HGNC_gene_symbol = 'TP53'
    AND normalized_count IS NOT NULL
    AND normalized_count > 0),
  --
  -- And we'll select the variant data for our cohort,
  -- we're going to be comparing variant types (SNP, DEL, etc)
  --
  cohortVar AS (
  SELECT
    Variant_Type,
    sample_barcode_tumor AS sample_barcode
  FROM
    `isb-cgc.TCGA_hg19_data_v0.Somatic_Mutation_MC3`
  WHERE
    SYMBOL = 'TP53' ),
  --
  -- then we join the expression and variant data using sample barcodes
  --
  cohort AS (
  SELECT
    cohortExpr.sample_barcode AS sample_barcode,
    Variant_Type AS group_name,
    expr
  FROM
    cohortExpr
  JOIN
    cohortVar
  ON
    cohortExpr.sample_barcode = cohortVar.sample_barcode ),
  --
  -- First part of the calculation, the grand mean (over everything)
  --
  grandMeanTable AS (
  SELECT
    AVG(expr) AS grand_mean
  FROM
    cohort )
SELECT sample_barcode, group_name, expr
FROM cohort"
```

``` r
# To see the R console output with query processing information, turn queit to FALSE
dat <- bq_project_query(billing, dat_query, quiet = TRUE)
# Transform the query result into a tibble
dat <- bq_table_download(dat, quiet = TRUE)
dat
```

    ## # A tibble: 375 x 3
    ##    sample_barcode   group_name  expr
    ##    <chr>            <chr>      <dbl>
    ##  1 TCGA-A8-A084-01A INS         2.49
    ##  2 TCGA-E2-A108-01A SNP         3.17
    ##  3 TCGA-EW-A6SD-01A SNP         3.45
    ##  4 TCGA-AR-A1AR-01A DEL         2.84
    ##  5 TCGA-D8-A1J9-01A SNP         3.11
    ##  6 TCGA-AO-A03V-01A SNP         3.23
    ##  7 TCGA-D8-A13Y-01A DEL         2.64
    ##  8 TCGA-D8-A27N-01A SNP         3.04
    ##  9 TCGA-BH-A0EE-01A DEL         2.42
    ## 10 TCGA-AC-A8OQ-01A SNP         2.86
    ## # ... with 365 more rows

``` r
dat %>% group_by(group_name) %>% summarize(mean=mean(expr), sd=sd(expr))
```

    ## # A tibble: 3 x 3
    ##   group_name  mean    sd
    ##   <chr>      <dbl> <dbl>
    ## 1 DEL         2.79 0.322
    ## 2 INS         2.64 0.116
    ## 3 SNP         3.22 0.313

``` r
anova(lm(data=dat, expr~group_name))
```

    ## Analysis of Variance Table
    ## 
    ## Response: expr
    ##             Df Sum Sq Mean Sq F value    Pr(>F)    
    ## group_name   2 12.460  6.2302  65.147 < 2.2e-16 ***
    ## Residuals  372 35.576  0.0956                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
