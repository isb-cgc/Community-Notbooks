Gene Set Scoring using hg38 mutation tables
================

### Created By: David L Gibbs (david.gibbs ( \~ at \~ ) systemsbiology ( \~ dot \~ ) org)

#### Updated By: Lauren Hagen (lauren.hagen ( \~ at \~ ) systemsbiology ( \~ dot \~ ) org)

#### Original Text from [Query of the Month January 2018](https://isb-cancer-genomics-cloud.readthedocs.io/en/latest/sections/QueryOfTheMonthClub.html#january2018)

In this notebook, we’re going to implement a common bioinformatics task:
gene set scoring. In this procedure, we will compare the *joint*
expression of a set of genes between two groups.

Gene sets frequently result from experiments in which, for example, an
expression is compared between two groups (e.g. control and treatment),
and the genes that are significantly differentially expressed between
the two groups form the “gene set”. Given a gene set, numerous
approaches have been proposed to assign functional interpretations to a
particular list of genes.

A related approach is to consider pathways as gene sets: each gene set
will, therefore, be a canonical representation of a biological process
compiled by domain experts. We will use the
[WikiPathways](https://www.wikipathways.org/index.php/Download_Pathways)
gene sets that have been previously assembled and stored in BigQuery. In
total there were 381 pathways downloaded from WikiPathways. In the BQ
table, each row contains a pathway and a gene associated with that
pathway.

Our task will be to determine which genesets are differentially
expressed when comparing two groups of samples. We will define groups of
samples using the new hg38 Somatic Mutations table from the NCI Genomic
Data Commons.

The BigQuery table at ISB-CGC is found at:
isb-cgc.TCGA\_hg38\_data\_v0.Somatic\_Mutation\_DR10

The BigQuery pathway table is found at:
isb-cgc.QotM.WikiPathways\_20170425\_Annotated

We will implement a method found in a paper titled “Gene set enrichment
analysis made simple” (Rafael A Irrazary et al, PMID 20048385). They
propose a simple solution that outperforms a popular and more complex
method known as GSEA.

In short, (I’ll explain later), the gene set score comes from an average
of T-tests, where the T-statistics come from testing each gene for
differential expression between the two groups. The statistic is then
weighted by the square root of the sample size (number of genes in the
set), so that with larger gene sets, the ‘significant’ effect size can
get pretty small.

At the end of it, the full query processes 29.8GB in 10.1s and costs
\~$0.15.

First, we will need to import the required libraries and create a client
variable. For more information see [‘Quick Start Guide to
ISB-CGC’](INSERT%20LINK) and alternative authentication methods can be
found
[here](https://googleapis.github.io/google-cloud-python/latest/core/auth.html).

``` r
library(bigrquery)
library(dplyr)
library(stringr)
library(ggplot2)
```

Then let us set up the billing variables we will be using in this
notebook:

``` r
billing <- 'isb-cgc-02-0001' # Insert your project ID in the ''
```

Let’s get started on creating the full query. First, which tissue type
should we focus on? Let’s choose PARP1 as an interesting gene – it
encodes an enzyme involved in DNA damage repair and is also the target
of some therapeutic drugs – and use the Somatic Mutation table to choose
a cancer type:

``` r
PARP1_query <- "SELECT
  project_short_name,
  COUNT(DISTINCT(sample_barcode_tumor)) AS n
FROM
  `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation_DR10`
WHERE
  Hugo_Symbol = 'PARP1'
GROUP BY
  project_short_name
ORDER BY
  n DESC"
# To see the R console output with query processing information, turn queit to FALSE
PARP1_result <- bq_project_query(billing, PARP1_query, quiet = TRUE)
# Transform the query result into a tibble
PARP1_result <- bq_table_download(PARP1_result, quiet = TRUE)
PARP1_result
```

    ## # A tibble: 22 x 2
    ##    project_short_name     n
    ##    <chr>              <int>
    ##  1 TCGA-UCEC             46
    ##  2 TCGA-COAD             22
    ##  3 TCGA-STAD             22
    ##  4 TCGA-BRCA             17
    ##  5 TCGA-SKCM             16
    ##  6 TCGA-BLCA             15
    ##  7 TCGA-LUSC             12
    ##  8 TCGA-CESC              9
    ##  9 TCGA-GBM               8
    ## 10 TCGA-LUAD              6
    ## # ... with 12 more rows

The result of that query shows that 46 tumor samples have PARP1
mutations in UCEC, followed by COAD and STAD with 22 each. That’s a big
lead by UCEC, so let’s focus our work there.

Here’s where our main query will begin. Since this is standard SQL,
we’ll be naming each subtable, and the full query can be constructed
by concatenating each of the following sub-queries. Many of the sub
query code chunks that display the sub query as a table are commented
out, so that the extra queries will be run and therefore are not be
charged. To quickly remove the comments in these code chucks, select the
code within the code chuck and press ctrl+shift+c.

``` r
s1 <- "WITH
s1 AS (
  SELECT
    sample_barcode_tumor AS sample_barcode
  FROM
    `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation_DR10`
  WHERE
    project_short_name = 'TCGA-UCEC'
  GROUP BY
    1
  )"
query1 <- str_c(s1,"SELECT * FROM s1")
```

``` r
# To see the R console output with query processing information, turn queit to FALSE
query1_result <- bq_project_query(billing, query1, quiet = TRUE) 
# Transform the query result into a tibble
query1_result <- bq_table_download(query1_result, quiet = TRUE)
query1_result
```

    ## # A tibble: 530 x 1
    ##    sample_barcode  
    ##    <chr>           
    ##  1 TCGA-A5-A0VP-01A
    ##  2 TCGA-AJ-A5DW-01A
    ##  3 TCGA-AP-A1DK-01A
    ##  4 TCGA-AX-A2HC-01A
    ##  5 TCGA-B5-A1MR-01A
    ##  6 TCGA-E6-A1LX-01A
    ##  7 TCGA-EO-A22X-01A
    ##  8 TCGA-QF-A5YS-01A
    ##  9 TCGA-2E-A9G8-01A
    ## 10 TCGA-A5-A0G1-01A
    ## # ... with 520 more rows

This query returns 530 tumor sample barcodes with at least one known
somatic mutation. (The TCGA Biospecimen table includes information for a
total of 553 UCEC tumor samples, but some may have not been sequenced or
may have no somatic mutations – the former being more likely than the
latter.) Recall that somatic mutations are variants in the DNA that are
found when comparing the tumor sequence to the ‘matched normal’ sequence
(typically from a blood sample, but sometimes from adjacent tissue).
Next, for all these samples, we’ll want to restrict our analysis to
samples for which we also have mRNA expression data:

``` r
sampleGroup <- "sampleGroup AS (
  SELECT
    sample_barcode
  FROM
    `isb-cgc.TCGA_hg38_data_v0.RNAseq_Gene_Expression`
  WHERE
    project_short_name = 'TCGA-UCEC'
    AND sample_barcode IN
    (select sample_barcode from s1)
  GROUP BY
    1
  )"
query2 <- str_c(s1,",", sampleGroup, "SELECT * FROM sampleGroup")
```

``` r
# # To see the R console output with query processing information, turn queit to FALSE
# query2_result <- bq_project_query(billing, query2, quiet = TRUE) 
# # Transform the query result into a tibble
# query2_result <- bq_table_download(query2_result, quiet = TRUE)
# query2_result
```

Now we have 526 samples for which we have gene expression and somatic
mutation calls. We are interested in partitioning this group into two
parts: one with a mutation of interest, and one without. So let’s gather
barcodes for tumors with non-silent mutations in PARP1.

``` r
groups <- "
--
-- The first group has non-synonymous mutations in PARP1
--
grp1 AS (
  SELECT
    sample_barcode_tumor AS sample_barcode
  FROM
    `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation_DR10`
  WHERE
    Hugo_Symbol = 'PARP1'
    AND One_Consequence <> 'synonymous_variant'
    AND sample_barcode_tumor IN (
      SELECT
        sample_barcode
      FROM
        sampleGroup )
    GROUP BY sample_barcode
  ),
--
-- group 2 is the rest of the samples
--
grp2 AS (
  SELECT
    sample_barcode
  FROM
    sampleGroup
  WHERE
    sample_barcode NOT IN (
      SELECT
        sample_barcode
      FROM
        grp1)
)
"
query3_grp1 <- str_c(s1, ",", sampleGroup, ",", groups, "SELECT * FROM grp1")
query3_grp2 <- str_c(s1, ",", sampleGroup, ",", groups, "SELECT * FROM grp2")
```

``` r
# # To see the R console output with query processing information, turn queit to FALSE
# query3_grp1_result <- bq_project_query(billing, query3_grp1, quiet = TRUE) 
# # Transform the query result into a tibble
# query3_grp1_result <- bq_table_download(query3_grp1_result, quiet = TRUE)
# query3_grp1_result
# length(query3_grp1_result)
```

``` r
# # To see the R console output with query processing information, turn queit to FALSE
# query3_gp2_result <- bq_project_query(billing, query3_grp2, quiet = TRUE) 
# # Transform the query result into a tibble
# query3_grp2_result <- bq_table_download(query3_grp2_result, quiet = TRUE)
# query3_grp2_result
# length(query3_grp2_result)
```

This results in 41 tumor samples with non-synonymous PARP1 variants and
485 samples without.

Next we’re going to summarize the gene expression within each of these
groups. This will be used for calulating T-statistics in the following
portion of the query we are constructing. For each gene, we’ll take the
mean, variance, and count of samples.

``` r
summaries = "
-- Summaries for Group 1 (with mutation)
--
summaryGrp1 AS (
  select
    gene_name as symbol,
    AVG(LOG10( HTSeq__FPKM_UQ +1)) as genemean,
    VAR_SAMP(LOG10( HTSeq__FPKM_UQ +1)) as genevar,
    count(sample_barcode) as genen
  FROM
    `isb-cgc.TCGA_hg38_data_v0.RNAseq_Gene_Expression`
  WHERE
    sample_barcode IN (select sample_barcode FROM grp1)
    AND gene_name IN (
      SELECT
        Symbol as gene_name
      FROM
        `isb-cgc.QotM.WikiPathways_20170425_Annotated`
    )
  GROUP BY
    gene_name
),
--
-- Summaries for Group 2 (without mutation)
--
summaryGrp2 AS (
  select
    gene_name as symbol,
    AVG(LOG10( HTSeq__FPKM_UQ +1)) as genemean,
    VAR_SAMP(LOG10( HTSeq__FPKM_UQ +1)) as genevar,
    count(sample_barcode) as genen
  FROM
    `isb-cgc.TCGA_hg38_data_v0.RNAseq_Gene_Expression`
  WHERE
    sample_barcode IN (select sample_barcode FROM grp2)
    AND gene_name IN (
      SELECT
        Symbol as gene_name
      FROM
        `isb-cgc.QotM.WikiPathways_20170425_Annotated`
    )
  GROUP BY
    gene_name
)
--
"
query4_grp1 <- str_c(s1, ",", sampleGroup, ",", groups, ",", summaries, "SELECT * FROM summaryGrp1")
query4_grp2 <- str_c(s1, ",", sampleGroup, ",", groups, ",", summaries, "SELECT * FROM summaryGrp2")
```

``` r
# # To see the R console output with query processing information, turn queit to FALSE
# query4_grp1_result <- bq_project_query(billing, query4_grp1, quiet = TRUE) 
# # Transform the query result into a tibble
# query4_grp1_result <- bq_table_download(query4_grp1_result, quiet = TRUE)
# query4_grp1_result
# length(query4_grp1_result)
```

``` r
# # To see the R console output with query processing information, turn queit to FALSE
# query4_gp2_result <- bq_project_query(billing, query4_grp2, quiet = TRUE)
# # Transform the query result into a tibble
# query4_grp2_result <- bq_table_download(query4_grp2_result, quiet = TRUE)
# query4_grp2_result
# length(query4_grp2_result)
```

This results in two sets of summaries for 4,822 genes. (There are 4,962
unique gene symbols in the WikiPathways table, but 140 of them do not
match any of the symbols in the TCGA hg38 expression table.) With this,
we are ready to calculate T-statistics. Here we’re going to use a two
sample T-test assuming independent variance (and that we have enough
samples to assume that). The T-statistic is found by taking the
difference in means (of gene expression between our two groups), and
normalizing it by measures of variance and sample size. Here, we want to
keep T-statistics that are zero, which might come from having zero
variance, because having a T-stat for each gene is important in the gene
set score, even if it’s a zero. To do that, you’ll see the use of an IF
statement below.

``` r
tStatsPerGene = "
tStatsPerGene AS (
SELECT
  grp1.symbol as symbol,
  grp1.genen as grp1_n,
  grp2.genen AS grp2_n,
  grp1.genemean AS grp1_mean,
  grp2.genemean AS grp2_mean,
  grp1.genemean - grp2.genemean as meandiff,
  IF ((grp1.genevar > 0
       AND grp2.genevar > 0
       AND grp1.genen > 0
       AND grp2.genen > 0),
    (grp1.genemean - grp2.genemean) / SQRT( (POW(grp1.genevar,2)/grp1.genen)+ (POW(grp2.genevar,2)/grp2.genen) ),
    0.0) AS tstat
FROM
  summaryGrp1 as grp1
  JOIN
  summaryGrp2 AS grp2
  ON
  grp1.symbol = grp2.symbol
GROUP BY
  grp1.symbol,
  grp1.genemean,
  grp2.genemean,
  grp1.genevar,
  grp2.genevar,
  grp1.genen,
  grp2.genen
)
"
query5 <- str_c(s1, ",", sampleGroup, ",", groups, ",", summaries, ",", tStatsPerGene, "SELECT * FROM tStatsPerGene")
```

``` r
# To see the R console output with query processing information, turn queit to FALSE
tStatsPerGene_results <- bq_project_query(billing, query5, quiet = TRUE) 
# Transform the query result into a tibble
tStatsPerGene_results <- bq_table_download(tStatsPerGene_results, quiet = TRUE)
head(tStatsPerGene_results)
```

    ## # A tibble: 6 x 7
    ##   symbol   grp1_n grp2_n grp1_mean grp2_mean meandiff  tstat
    ##   <chr>     <int>  <int>     <dbl>     <dbl>    <dbl>  <dbl>
    ## 1 CYP2E1       41    489      3.48      3.65  -0.174  -2.73 
    ## 2 GPR174       41    489      3.52      2.99   0.532   2.56 
    ## 3 UBE2D2       41    489      6.03      5.99   0.0434  8.15 
    ## 4 HLA-DRB1     41    489      6.59      6.63  -0.0402 -0.947
    ## 5 DTX3L        41    489      5.34      5.51  -0.167  -8.43 
    ## 6 KIF20B       41    489      4.73      4.70   0.0387  3.51

We are now going to plot the t-statistic using ggplot2 library.

``` r
ggplot(data=tStatsPerGene_results) +
  geom_point(mapping = aes(x = tstat, y = abs(meandiff), color=abs(tstat)))
```

![](Gene_Set_Scoring_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggplot(data=tStatsPerGene_results) +
  geom_density(mapping = aes(x=tstat))
```

![](Gene_Set_Scoring_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

OK\! We have a distribution of T statistics. The question is whether
there’s some hidden structure to these values. Are there gene sets with
unusually high T-statistics? And do these gene sets make any sort of
biological sense? Let’s find out\!

Now we are going to integrate our gene set table. This is as easy as
doing a table join.

``` r
geneSetTable = "geneSetTable AS (
SELECT
  gs.pathway,
  gs.wikiID,
  gs.Symbol,
  st.grp1_n,
  st.grp2_n,
  st.grp1_mean,
  st.grp2_mean,
  st.meandiff,
  st.tstat
FROM
  `isb-cgc.QotM.WikiPathways_20170425_Annotated` as gs
JOIN
  tStatsPerGene as st
ON
  st.symbol = gs.symbol
GROUP BY
  gs.pathway,
  gs.wikiID,
  gs.Symbol,
  st.grp1_n,
  st.grp2_n,
  st.grp1_mean,
  st.grp2_mean,
  st.meandiff,
  st.tstat
)
"
query6 = str_c(s1, ",", sampleGroup, ",", groups, ",", summaries, ",", tStatsPerGene, ",", geneSetTable, "SELECT * FROM geneSetTable")
```

``` r
# # To see the R console output with query processing information, turn queit to FALSE
# geneSetTable_results <- bq_project_query(billing, query6, quiet = TRUE) 
# # Transform the query result into a tibble
# geneSetTable_results <- bq_table_download(geneSetTable_results, quiet = TRUE)
# head(geneSetTable_results)
```

That’s it\! For each gene in the pathways (gene sets) table, we have
joined in the T-statistic comparing our two groups. Now for the gene set
score\! To get this, we’re going to simply average over the T’s within
each pathway, and scale the result by the square root of the number of
genes. When the number of genes gets large (reasonably so), the value
approximates a Z-score. In this way, using R for example, we could get a
p-value and perform multiple testing correction in order to control the
false discovery rate.

``` r
geneSetScores = "geneSetScores AS (
SELECT
  pathway,
  wikiID,
  COUNT(symbol) AS n_genes,
  AVG(ABS(meandiff)) AS avgAbsDiff,
  (SQRT(COUNT(symbol))/COUNT(symbol)) * SUM(tstat) AS score
FROM
  geneSetTable
GROUP BY
  pathway,
  wikiID )
--
--
SELECT
  *
FROM
  geneSetScores
ORDER BY
  score DESC "
query7 = str_c(s1, ",", sampleGroup, ",", groups, ",", summaries, ",", tStatsPerGene, ",", geneSetTable, ",", geneSetScores)
```

``` r
# To see the R console output with query processing information, turn queit to FALSE
geneSetScores_results <- bq_project_query(billing, query7, quiet = TRUE)
# Transform the query result into a tibble
geneSetScores_results <- bq_table_download(geneSetScores_results, quiet = TRUE)
head(geneSetScores_results)
```

    ## # A tibble: 6 x 5
    ##   pathway                         wikiID           n_genes avgAbsDiff score
    ##   <chr>                           <chr>              <int>      <dbl> <dbl>
    ## 1 Retinoblastoma (RB) in Cancer   WikiPathways_20~      85     0.0907 101. 
    ## 2 DNA Replication                 WikiPathways_20~      41     0.0964  87.3
    ## 3 Cell Cycle                      WikiPathways_20~      98     0.0808  63.0
    ## 4 Eukaryotic Transcription Initi~ WikiPathways_20~      37     0.0467  55.0
    ## 5 Glycolysis and Gluconeogenesis  WikiPathways_20~      48     0.131   54.4
    ## 6 miRNA Regulation of DNA Damage~ WikiPathways_20~      70     0.0845  53.2

``` r
ggplot(geneSetScores_results,aes(x=avgAbsDiff, y=score, color=n_genes, size=n_genes, label=pathway)) +
  geom_point() +
  geom_text(data=subset(geneSetScores_results, score >= score[3]),hjust=0,nudge_x = 0.005)
```

![](Gene_Set_Scoring_files/figure-gfm/unnamed-chunk-4-1.png)<!-- --> So,
we see that ‘Retinoblastoma (RB) in Cancer’ is in the top spot with a
score way above the \#2 position. Why might that be? Well, PARP1 is
involved in DNA damage repair, specifically through the non-homologous
endjoining (NHEJ) mechanism. Samples that are deficient in PARP1 are
going to have a hard time repairing DNA breaks, which makes cancer more
likely. So, RB1 might need to take up the slack, and indeed it’s known
as a ‘tumor suppressor protein’: when DNA is damaged, the cell cycle
needs to freeze, which happens to be one of RB1’s special tricks, and
also probably why we see the next two top ranked pathways ‘DNA
Replication’ and ‘Cell Cycle’.

For more on this topic see:

  - Retinoblastoma (RB) in Cancer (Homo sapiens) ([Wiki pathway
    WP2446](https://www.wikipathways.org/index.php/Pathway:WP2446))
  - RB1 gene ([Wikipedia
    entry](https://en.wikipedia.org/wiki/Retinoblastoma_protein))
  - Direct involvement of retinoblastoma family proteins in DNA repair
    by non-homologous end-joining. ([Cook et
    al, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25818292))

The full query is written below and called upon with the BigQuery python
magic command to ease the readablity of the SQL. It has been commented
out to save money by not running the same query a second time.

``` r
full_query <- "WITH
s1 AS (
  SELECT
    sample_barcode_tumor AS sample_barcode
  FROM
    `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation_DR10`
  WHERE
    project_short_name = 'TCGA-UCEC'
  GROUP BY
    1
),
---
sampleGroup AS (
SELECT
  sample_barcode
FROM
  `isb-cgc.TCGA_hg38_data_v0.RNAseq_Gene_Expression`
WHERE
  project_short_name = 'TCGA-UCEC'
  AND sample_barcode IN
  (select sample_barcode from s1)
GROUP BY
  1
),
--
-- The first group has non-synonymous mutations in PARP1
--
grp1 AS (
SELECT
  sample_barcode_tumor AS sample_barcode
FROM
  `isb-cgc.TCGA_hg38_data_v0.Somatic_Mutation_DR10`
WHERE
  Hugo_Symbol = 'PARP1'
  AND One_Consequence <> 'synonymous_variant'
  AND sample_barcode_tumor IN (
    SELECT
      sample_barcode
    FROM
      sampleGroup )
  GROUP BY sample_barcode
),
--
-- group 2 is the rest of the samples
--
grp2 AS (
SELECT
  sample_barcode
FROM
  sampleGroup
WHERE
  sample_barcode NOT IN (
    SELECT
      sample_barcode
    FROM
      grp1)
),
---
-- Summaries for Group 1 (with mutation)
--
summaryGrp1 AS (
  select
    gene_name as symbol,
    AVG(LOG10( HTSeq__FPKM_UQ +1)) as genemean,
    VAR_SAMP(LOG10( HTSeq__FPKM_UQ +1)) as genevar,
    count(sample_barcode) as genen
  FROM
    `isb-cgc.TCGA_hg38_data_v0.RNAseq_Gene_Expression`
  WHERE
    sample_barcode IN (select sample_barcode FROM grp1)
    AND gene_name IN (
      SELECT
        Symbol as gene_name
      FROM
        `isb-cgc.QotM.WikiPathways_20170425_Annotated`
    )
  GROUP BY
    gene_name
),
--
-- Summaries for Group 2 (without mutation)
--
summaryGrp2 AS (
  select
    gene_name as symbol,
    AVG(LOG10( HTSeq__FPKM_UQ +1)) as genemean,
    VAR_SAMP(LOG10( HTSeq__FPKM_UQ +1)) as genevar,
    count(sample_barcode) as genen
  FROM
    `isb-cgc.TCGA_hg38_data_v0.RNAseq_Gene_Expression`
  WHERE
    sample_barcode IN (select sample_barcode FROM grp2)
    AND gene_name IN (
      SELECT
        Symbol as gene_name
      FROM
        `isb-cgc.QotM.WikiPathways_20170425_Annotated`
    )
  GROUP BY
    gene_name
),
---
tStatsPerGene AS (
SELECT
  grp1.symbol as symbol,
  grp1.genen as grp1_n,
  grp2.genen AS grp2_n,
  grp1.genemean AS grp1_mean,
  grp2.genemean AS grp2_mean,
  grp1.genemean - grp2.genemean as meandiff,
  IF ((grp1.genevar > 0
       AND grp2.genevar > 0
       AND grp1.genen > 0
       AND grp2.genen > 0),
    (grp1.genemean - grp2.genemean) / SQRT( (POW(grp1.genevar,2)/grp1.genen)+ (POW(grp2.genevar,2)/grp2.genen) ),
    0.0) AS tstat
FROM
  summaryGrp1 as grp1
  JOIN
  summaryGrp2 AS grp2
  ON
  grp1.symbol = grp2.symbol
GROUP BY
  grp1.symbol,
  grp1.genemean,
  grp2.genemean,
  grp1.genevar,
  grp2.genevar,
  grp1.genen,
  grp2.genen
),
---
geneSetTable AS (
SELECT
  gs.pathway,
  gs.wikiID,
  gs.Symbol,
  st.grp1_n,
  st.grp2_n,
  st.grp1_mean,
  st.grp2_mean,
  st.meandiff,
  st.tstat
FROM
  `isb-cgc.QotM.WikiPathways_20170425_Annotated` as gs
JOIN
  tStatsPerGene as st
ON
  st.symbol = gs.symbol
GROUP BY
  gs.pathway,
  gs.wikiID,
  gs.Symbol,
  st.grp1_n,
  st.grp2_n,
  st.grp1_mean,
  st.grp2_mean,
  st.meandiff,
  st.tstat
),
---
geneSetScores AS (
SELECT
  pathway,
  wikiID,
  COUNT(symbol) AS n_genes,
  AVG(ABS(meandiff)) AS avgAbsDiff,
  (SQRT(COUNT(symbol))/COUNT(symbol)) * SUM(tstat) AS score
FROM
  geneSetTable
GROUP BY
  pathway,
  wikiID )
--
--
SELECT
  *
FROM
  geneSetScores
ORDER BY
  score DESC"
```

``` r
# # To see the R console output with query processing information, turn queit to FALSE
# full_query_results <- bq_project_query(billing, full_query, quiet = TRUE)
# # Transform the query result into a tibble
# full_query_results <- bq_table_download(full_query_results, quiet = TRUE)
# head(full_query_results)
```
