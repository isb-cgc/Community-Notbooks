How to Create Cohorts
================

### Created By: David L Gibbs (david.gibbs ( \~ at \~ ) systemsbiology ( \~ dot \~ ) org)

#### Updated By: Lauren Hagen (lauren.hagen ( \~ at \~ ) systemsbiology ( \~ dot \~ ) org)

#### Original Text from [Create Cohorts – part 1](https://github.com/isb-cgc/examples-Python/blob/master/notebooks/Creating%20TCGA%20cohorts%20--%20part%201.ipynb)

# Creating TCGA cohorts

This notebook will show you how to create a TCGA cohort using the
publicly available TCGA BigQuery tables that the
[ISB-CGC](http://isb-cgc.org) project has produced based on the
open-access [TCGA](http://cancergenome.nih.gov/) data available at the
[Data Portal](https://tcga-data.nci.nih.gov/tcga/). You will need to
have access to a Google Cloud Platform (GCP) project in order to use
BigQuery. If you don’t already have one, you can sign up for a
[free-trial](https://cloud.google.com/free-trial/).

We are not attempting to provide a thorough BigQuery or IPython tutorial
here, as a wealth of such information already exists. Here are some
links to some resources that you might find useful: \*
[BigQuery](https://cloud.google.com/bigquery/what-is-bigquery) \* the
BigQuery [web UI](https://console.cloud.google.com/bigquery) \* where
you can run queries interactively \* [Jupyter
Notebooks](http://jupyter.org/) \* [Google Cloud
Datalab](https://cloud.google.com/datalab/) \* interactive cloud-based
platform for analyzing data built on the Jupyter Notebooks \* [Google
Colaboratory](https://colab.research.google.com/) \* Free Jupyter
Notebook environment that runs in your browser

There are also many tutorials and samples available on github (see, in
particular, the
[datalab](https://github.com/GoogleCloudPlatform/datalab) repo, the
[Google Genomics](https://github.com/googlegenomics) project), and our
own [Community
Notebooks](https://github.com/isb-cgc/Community-Notebooks).

OK then, let’s get started\! In order to work with BigQuery, the first
thing you need to do is load the libraries required for using BigQuery
in R:

``` r
library(bigrquery)
library(dplyr)
library(dbplyr)
```

    ## 
    ## Attaching package: 'dbplyr'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     ident, sql

Then let us set up some of the variables we will be using in this
notebook:

``` r
billing <- 'isb-cgc-02-0001' # Insert your project ID in the ''
theTable <- "isb-cgc.TCGA_bioclin_v0.Clinical" # The convention for calling a table is project.dataset.table
```

The next thing you need to know is how to access the specific tables you
are interested in. BigQuery tables are organized into datasets, and
datasets are owned by a specific GCP project. The tables we will be
working with in this notebook are in a dataset called
**`TCGA_bioclin_v0`**, owned by the **`isb-cgc`** project. A full table
identifier is of the form `<project_id>.<dataset_id>.<table_id>`. Let’s
start by getting some basic information about the tables in this
dataset:

``` r
# Let us look which tables are in the TCGA_bioclin_v0 dataset
tables<-list_tables("isb-cgc", "TCGA_bioclin_v0") # the convention is project name then dataset
tables
```

    ## [1] "Annotations" "Biospecimen" "Clinical"    "clinical_v1"

In this tutorial, we are going to look at a few different ways that we
can use the information in these tables to create cohorts. Now, you
maybe asking what we mean by “cohort” and why you might be interested in
*creating* one, or maybe what it even means to “create” a cohort. The
TCGA dataset includes clinical, biospecimen, and molecular data from
over 10,000 cancer patients who agreed to be a part of this landmark
research project to build [The Cancer Genome
Atlas](http://cancergenome.nih.gov/). This large dataset was originally
organized and studied according to [cancer
type](http://cancergenome.nih.gov/cancersselected) but now that this
multi-year project is nearing completion, with over 30 types of cancer
and over 10,000 tumors analyzed, **you** have the opportunity to look at
this dataset from whichever angle most interests you. Maybe you are
particularly interested in early-onset cancers, or gastro-intestinal
cancers, or a specific type of genetic mutation. This is where the idea
of a “cohort” comes in. The original TCGA “cohorts” were based on cancer
type (aka “study”), but now you can define a cohort based on virtually
any clinical or molecular feature by querying these BigQuery tables. A
cohort is simply a list of samples, using the [TCGA
barcode](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/)
system. Once you have created a cohort you can use it in any number of
ways: you could further explore the data available for one cohort, or
compare one cohort to another, for example.

## Exploring the Clinical data table

Let’s start by looking at the clinical data table. The TCGA dataset
contains a few very basic clinical data elements for almost all
patients, and contains additional information for some tumor types only.
For example smoking history information is generally available only for
lung cancer patients, and BMI (body mass index) is only available for
tumor types where that is a known significant risk factor. Let’s take a
look at the clinical data table and see how many different pieces of
information are available to us:

``` r
# Create the SQL query
sql_query1 <- "SELECT
                column_name
              FROM
                `isb-cgc.TCGA_bioclin_v0.INFORMATION_SCHEMA.COLUMNS`
              WHERE table_name = 'Clinical'"
# Use BigQuery to run the SQL query on the Clincal table from the TCGA_bioclin_v0 dataset
# To see the R console output with query processing information, turn queit to FALSE
result <- bq_project_query(billing, sql_query1, quiet = TRUE) 
# Transform the query result into a tibble
result <- bq_table_download(result, quiet = TRUE)
result
```

    ## # A tibble: 73 x 1
    ##    column_name                   
    ##    <chr>                         
    ##  1 program_name                  
    ##  2 case_barcode                  
    ##  3 case_gdc_id                   
    ##  4 program_dbgap_accession_number
    ##  5 project_short_name            
    ##  6 project_name                  
    ##  7 disease_code                  
    ##  8 gender                        
    ##  9 vital_status                  
    ## 10 race                          
    ## # ... with 63 more rows

That’s a lot of fields\! We can also get at the schema through dplyr:

``` r
# Connect to BigQuery
con <- dbConnect(
  bigrquery::bigquery(),
  project = "isb-cgc",
  dataset = "TCGA_bioclin_v0",
  billing = billing
)
con
```

    ## <BigQueryConnection>
    ##   Dataset: isb-cgc.TCGA_bioclin_v0
    ##   Billing: isb-cgc-02-0001

``` r
# Create a table from the Clinical table in the TCGA_bioclin_v0 dataset
tcga_clinical <- tbl(con, "Clinical")
# Get the columns names
columns <- colnames(tcga_clinical)
# Print interesting information about the columns
cat("The first 5 columns names are: ")
```

    ## The first 5 columns names are:

``` r
cat(columns[1:5], sep=", ")
```

    ## program_name, case_barcode, case_gdc_id, program_dbgap_accession_number, project_short_name

``` r
cat("\nThere are", length(columns), " columns in the TCGA_bioclin_v0 table.")
```

    ## 
    ## There are 73  columns in the TCGA_bioclin_v0 table.

Let’s look at these fields and see which ones might be the most
“interesting”, by looking at how many times they are filled-in (not
NULL), or how much variation exists in the values. If we wanted to look
at just a single field, “tobacco\_smoking\_history” for example, we
could use a very simple query to get a basic summary:

``` r
tobacco_query1 <- "SELECT tobacco_smoking_history,
                COUNT(*) AS n
               FROM `isb-cgc.TCGA_bioclin_v0.Clinical`
               GROUP BY tobacco_smoking_history
               ORDER BY n DESC"
# Run the query
tobacco1 <- bq_project_query(billing, tobacco_query1, quiet = TRUE) 
# Create a dataframe with the results from the query
tobacco1 <- bq_table_download(tobacco1, quiet = TRUE)
# Show the dataframe
tobacco1
```

    ## # A tibble: 6 x 2
    ##   tobacco_smoking_history     n
    ##                     <int> <int>
    ## 1                      NA  8316
    ## 2                       1   865
    ## 3                       4   799
    ## 4                       2   710
    ## 5                       3   568
    ## 6                       5    57

``` r
# Using pipes and dbplyr to Query:
tobacco_query2 <- tcga_clinical %>%
  count(tobacco_smoking_history) %>%
  arrange(desc(n))
# Show the SQL Query
tobacco_query2 %>%
  show_query()
```

    ## <SQL>
    ## SELECT `tobacco_smoking_history`, count(*) AS `n`
    ## FROM `Clinical`
    ## GROUP BY `tobacco_smoking_history`
    ## ORDER BY `n` DESC

``` r
# Create a dataframe with the results from the query 
tobacco2 <- tobacco_query2 %>% 
  collect()
```

    ## 
    Running job 'isb-cgc-02-0001.job_w4NHkIrJkl1iznBiJGjvkxQfjrOO.US' [-]  1s
    Running job 'isb-cgc-02-0001.job_w4NHkIrJkl1iznBiJGjvkxQfjrOO.US' [\]  2s
    Running job 'isb-cgc-02-0001.job_w4NHkIrJkl1iznBiJGjvkxQfjrOO.US' [|]  2s
    ## Complete
    ## Billed: 10.49 MB
    ## Downloading 6 rows in 1 pages.

``` r
# Show the resulting dataframe
tobacco2
```

    ## # A tibble: 6 x 2
    ##   tobacco_smoking_history     n
    ##                     <int> <int>
    ## 1                      NA  8316
    ## 2                       1   865
    ## 3                       4   799
    ## 4                       2   710
    ## 5                       3   568
    ## 6                       5    57

For more information on dbplyr visit [dbplyer
Overview](https://dbplyr.tidyverse.org/) and [Writing SQL with
dbplyr](https://dbplyr.tidyverse.org/articles/sql.html). Another useful
resourse for using dplyr, dbplyr and tidyverse is “R for Data Science”
by Garrett Grolemund and Hadley Wickham which can be accessed for free
[here](https://www.tidyverse.org/learn/).

But if we want to loop over all fields and get a sense of which fields
might provide us with useful criteria for specifying a cohort, we’ll
want to automate that. We’ll put a threshold on the minimum number of
patients that we expect information for, and the maximum number of
unique values (since fields such as the “ParticipantBarcode” will be
unique for every patient and, although we will need that field later,
it’s probably not useful for defining a
cohort).

``` r
clinicalTable <- "isb-cgc.TCGA_bioclin_v0.Clinical" # The convention for calling a table is project.dataset.table
# Create the SQL query
sql <- "SELECT
          COUNT(program_name)
        FROM
          `isb-cgc.TCGA_bioclin_v0.Clinical`"
# Use BigQuery to run the SQL query on the Clincal table from the TCGA_bioclin_v0 dataset
# and get the number of Patients in the dataset
numPatients <- bq_project_query(billing, sql, quiet = TRUE)
numPatients <- bq_table_download(numPatients, quiet = TRUE)
numPatients <- as.integer(numPatients) # Convert to integer
# Print the total number of patients
cat("The Clinical table describes a total of", numPatients, "patients")
```

    ## The Clinical table describes a total of 11315 patients

``` r
# let's set a threshold for the minimum number of values that a field should have,
# the maximum number of unique values, and either the highest cancer type or
# the mean and sigma of the row.
minNumPatients <- numPatients*0.80
maxNumValues <- 50

# Create a variable to be filled in by the for loop with the number
# interesting features
numInteresting <- 0

# Create a list to hold the results from the loop below
iList <- c()

# Loop over the fields and find the number of values with the number of unique
# values and the
for (field in columns) {
  query <- paste("SELECT",field,"FROM `isb-cgc.TCGA_bioclin_v0.Clinical`", sep = " ")
  tb <- bq_project_query(billing, query, quiet = TRUE)
  df <- bq_table_download(tb, quiet = TRUE)
  type <- class(df[[1]])
  if(type=="character") {
    freq <- as.data.frame(table(df))
    order_freq <- freq[order(freq[,2], decreasing = TRUE),]
    topFrac <- order_freq$Freq[1]/sum(freq[,2])
    if (sum(freq[,2]) >= minNumPatients) {
      if (length(freq[,1]) <= maxNumValues && length(freq[,1]) > 1){
        if ( topFrac < 0.90 ) {
          numInteresting <- numInteresting + 1
          iList <- append(iList, field)
          cat("\n     > ", field, " has ", sum(freq[,2]), " values with ", length(freq[,1]), " unique (",
              as.character(order_freq[1,1]), " occurs ", order_freq[1,2], " times)")
        }
      }
    }
  } else {
    if ( length(which(is.na(df[,1]) == FALSE)) >= minNumPatients) {
      iSd <- round(sd(df[[1]], na.rm = TRUE))
      iMean <- round(mean(df[[1]], na.rm = TRUE))
      if ( iSd > 0.1 ) {
        numInteresting <- numInteresting + 1
        iList <- append(iList, field)
        cat("\n     > ", field, " has ", length(which(is.na(df[,1]) == FALSE)), "value(S) (mean = ", iMean, ", sigma =", iSd, ")")
      }
    }
  }
}
```

    ## 
    ##      >  project_short_name  has  11315  values with  33  unique ( TCGA-BRCA  occurs  1098  times)
    ##      >  project_name  has  11315  values with  33  unique ( Breast Invasive Carcinoma  occurs  1098  times)
    ##      >  disease_code  has  11315  values with  33  unique ( BRCA  occurs  1098  times)
    ##      >  gender  has  11160  values with  2  unique ( FEMALE  occurs  5815  times)
    ##      >  vital_status  has  11156  values with  2  unique ( Alive  occurs  7534  times)
    ##      >  race  has  9835  values with  5  unique ( WHITE  occurs  8186  times)
    ##      >  age_at_diagnosis  has  11109 value(S) (mean =  59 , sigma = 14 )
    ##      >  days_to_birth  has  11041 value(S) (mean =  -21763 , sigma = 5266 )
    ##      >  days_to_last_known_alive  has  11102 value(S) (mean =  1037 , sigma = 1041 )
    ##      >  year_of_initial_pathologic_diagnosis  has  11030 value(S) (mean =  2008 , sigma = 4 )
    ##      >  person_neoplasm_cancer_status  has  10236  values with  2  unique ( TUMOR FREE  occurs  6507  times)
    ##      >  batch_number  has  11160 value(S) (mean =  203 , sigma = 135 )

``` r
cat("\n Found ", numInteresting, "potentially interesting features: \n", iList)
```

    ## 
    ##  Found  12 potentially interesting features: 
    ##  project_short_name project_name disease_code gender vital_status race age_at_diagnosis days_to_birth days_to_last_known_alive year_of_initial_pathologic_diagnosis person_neoplasm_cancer_status batch_number

The above helps us narrow down on which fields are likely to be the most
useful, but if you have a specific interest, for example in menopause or
HPV status, you can still look at those in more detail very easily:

``` r
# Using pipes and dbplyr to Query:
menopause_stat_query <- tcga_clinical %>%
  filter(!is.na(menopause_status)) %>%
  count(menopause_status) %>%
  arrange(desc(n))
# Show the SQL Query
menopause_stat_query %>%
  show_query()
```

    ## <SQL>
    ## SELECT `menopause_status`, count(*) AS `n`
    ## FROM `Clinical`
    ## WHERE (NOT(((`menopause_status`) IS NULL)))
    ## GROUP BY `menopause_status`
    ## ORDER BY `n` DESC

``` r
# Show the results of the Query
menopause_stat <- menopause_stat_query %>% 
  collect()
```

    ## Complete

    ## Billed: 10.49 MB

    ## Downloading 4 rows in 1 pages.

``` r
menopause_stat
```

    ## # A tibble: 4 x 2
    ##   menopause_status                                                        n
    ##   <chr>                                                               <int>
    ## 1 Post (prior bilateral ovariectomy OR >12 mo since LMP with no prio~  1291
    ## 2 Pre (<6 months since LMP AND no prior bilateral ovariectomy AND no~   389
    ## 3 Peri (6-12 months since last menstrual period)                         82
    ## 4 Indeterminate (neither Pre or Postmenopausal)                          54

We might wonder which specific tumor types have menopause information:

``` r
# Using pipes and dbplyr to Query:
menopause_type_query <- tcga_clinical %>%
  filter(!is.na(menopause_status)) %>%
  count(project_short_name) %>%
  arrange(desc(n))
# Show the SQL Query
menopause_type_query %>%
  show_query()
```

    ## <SQL>
    ## SELECT `project_short_name`, count(*) AS `n`
    ## FROM `Clinical`
    ## WHERE (NOT(((`menopause_status`) IS NULL)))
    ## GROUP BY `project_short_name`
    ## ORDER BY `n` DESC

``` r
# Show the results of the Query
menopause_type <- menopause_type_query %>% 
  collect()
```

    ## 
    Running job 'isb-cgc-02-0001.job_M6fooM6OGuejdF1SKdadR7ZZuwt3.US' [-]  1s
    Running job 'isb-cgc-02-0001.job_M6fooM6OGuejdF1SKdadR7ZZuwt3.US' [\]  1s
    Running job 'isb-cgc-02-0001.job_M6fooM6OGuejdF1SKdadR7ZZuwt3.US' [|]  2s
    Running job 'isb-cgc-02-0001.job_M6fooM6OGuejdF1SKdadR7ZZuwt3.US' [/]  2s
    Running job 'isb-cgc-02-0001.job_M6fooM6OGuejdF1SKdadR7ZZuwt3.US' [-]  2s
    Running job 'isb-cgc-02-0001.job_M6fooM6OGuejdF1SKdadR7ZZuwt3.US' [\]  3s
    Running job 'isb-cgc-02-0001.job_M6fooM6OGuejdF1SKdadR7ZZuwt3.US' [|]  3s
    ## Complete
    ## Billed: 10.49 MB
    ## Downloading 4 rows in 1 pages.

``` r
menopause_type
```

    ## # A tibble: 4 x 2
    ##   project_short_name     n
    ##   <chr>              <int>
    ## 1 TCGA-BRCA           1007
    ## 2 TCGA-UCEC            517
    ## 3 TCGA-CESC            237
    ## 4 TCGA-UCS              55

``` r
# Using pipes and dbplyr to Query:
hpv_stat_query <- tcga_clinical %>%
  filter(!is.na(hpv_status)) %>%
  count(hpv_status, hpv_calls) %>%
  filter( n > 20) %>%
  arrange(desc(n))

# Show the SQL Query
hpv_stat_query %>%
  show_query()
```

    ## <SQL>
    ## SELECT *
    ## FROM (SELECT `hpv_status`, `hpv_calls`, count(*) AS `n`
    ## FROM (SELECT *
    ## FROM `Clinical`
    ## WHERE (NOT(((`hpv_status`) IS NULL)))) `dbplyr_001`
    ## GROUP BY `hpv_status`, `hpv_calls`) `dbplyr_002`
    ## WHERE (`n` > 20.0)
    ## ORDER BY `n` DESC

``` r
# Show the results of the Query
hpv_stat <- hpv_stat_query %>% 
  collect()
```

    ## 
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [-]  1s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [\]  1s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [|]  2s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [/]  2s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [-]  3s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [\]  3s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [|]  3s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [/]  3s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [-]  4s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [\]  4s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [|]  4s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [/]  5s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [-]  5s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [\]  5s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [|]  6s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [/]  6s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [-]  6s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [\]  6s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [|]  7s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [/]  7s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [-]  8s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [\]  8s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [|]  8s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [/]  8s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [-]  9s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [\]  9s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [|]  9s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [/]  9s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [-] 10s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [\] 10s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [|] 11s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [/] 11s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [-] 11s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [\] 11s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [|] 12s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [/] 12s
    Running job 'isb-cgc-02-0001.job_JylvDjtcmi3wMktrc08h3lg5Icj3.US' [-] 12s
    ## Complete
    ## Billed: 10.49 MB
    ## Downloading 5 rows in 1 pages.

``` r
hpv_stat
```

    ## # A tibble: 5 x 3
    ## # Groups:   hpv_status [2]
    ##   hpv_status hpv_calls     n
    ##   <chr>      <chr>     <int>
    ## 1 Negative   <NA>        664
    ## 2 Positive   HPV16       238
    ## 3 Positive   HPV18        41
    ## 4 Positive   HPV33        25
    ## 5 Positive   HPV45        24

An additional factor to consider, when creating a cohort is that there
may be additional information that might lead one to exclude a
particular patient from a cohort. In certain instances, patients have
been redacted or excluded from analyses for reasons such as prior
treatment, etc, but since different researchers may have different
criteria for using or excluding certain patients or certain samples from
their analyses, an overview of the annoations can be found
[here](https://docs.gdc.cancer.gov/Encyclopedia/pages/Annotations_TCGA/).
These annotations have also been uploaded into a BigQuery table and can
be used in conjuction with the other BigQuery tables.

# Create a Cohort from Two Tables

Now that we have a better idea of what types of information is available
in the Clinical data table, let’s create a cohort consisting of female
breast-cancer patients, diagnosed at the age of 50 or younger.

In this next code sections, we define several queries with dplyr pipes
which will then allow us to use them in a final query. We will then save
the query to a tibble to allow it to be analyzed later. \* the first
query, called **`select_on_annotations`**, finds all patients in the
Annotations table which have either been ‘redacted’ or had ‘unacceptable
prior treatment’;  
\* the second query, **`select_on_clinical`** selects all female
breast-cancer patients who were diagnosed at age 50 or younger, while
also pulling out a few additional fields that might be of interest; and
\* the final query joins these two together and returns just those
patients that meet the clinical-criteria and do **not** meet the
exclusion-criteria.

``` r
tcga_annotations <- tbl(con, "Annotations")
```

``` r
select_on_annotations <- tcga_annotations %>%
  # Find all patients which have either been 'redacted' or had 'unacceptable prior treatment'
  filter(entity_type=="Patient" && (category=="History of unacceptable prior treatment related to a prior/other malignancy" | classification=="Redaction" )) %>%
  # Group by the case_barcode, category, classification
  group_by(case_barcode, category, classification) %>%
  # Summarise to have the group by clause added to the query
  summarise()

# Show the SQL Query but not run the query
select_on_annotations %>%
  show_query()
```

    ## <SQL>
    ## SELECT `case_barcode`, `category`, `classification`
    ## FROM `Annotations`
    ## WHERE (`entity_type` = 'Patient' AND (`category` = 'History of unacceptable prior treatment related to a prior/other malignancy' OR `classification` = 'Redaction'))
    ## GROUP BY `case_barcode`, `category`, `classification`

``` r
select_on_clinical <- tcga_clinical %>%
  # Find all of the female, under 50, breast cancer patients
  filter(disease_code == "BRCA" && age_at_diagnosis<=50 && gender=="FEMALE") %>%
  # Select interesting fields from the table
  select(case_barcode, vital_status, days_to_last_known_alive, ethnicity, histological_type, menopause_status)

# Show the SQL Query but not run the query
select_on_clinical %>%
  show_query()
```

    ## <SQL>
    ## SELECT `case_barcode`, `vital_status`, `days_to_last_known_alive`, `ethnicity`, `histological_type`, `menopause_status`
    ## FROM `Clinical`
    ## WHERE (`disease_code` = 'BRCA' AND `age_at_diagnosis` <= 50.0 AND `gender` = 'FEMALE')

``` r
early_onset_breast_cancer_query <- select_on_annotations %>%
  # Full join the two tables
  full_join(select_on_clinical, by = "case_barcode") %>%
  # patients that meet the clinical-criteria and do **not** meet the exclusion-criteria
  filter((is.na(category) | is.na(classification) ) && !is.na(case_barcode)) %>%
  # Only select the case barcode column
  select(case_barcode)

# Run the query and save to a tibble
early_onset_breast_cancer <- early_onset_breast_cancer_query %>%
  collect()
```

    ## 
    Running job 'isb-cgc-02-0001.job_Lbt8AICMkqRCZHkIWZgNLoFV2SCG.US' [-]  1s
    Running job 'isb-cgc-02-0001.job_Lbt8AICMkqRCZHkIWZgNLoFV2SCG.US' [\]  1s
    Running job 'isb-cgc-02-0001.job_Lbt8AICMkqRCZHkIWZgNLoFV2SCG.US' [|]  2s
    Running job 'isb-cgc-02-0001.job_Lbt8AICMkqRCZHkIWZgNLoFV2SCG.US' [/]  2s
    Running job 'isb-cgc-02-0001.job_Lbt8AICMkqRCZHkIWZgNLoFV2SCG.US' [-]  2s
    Running job 'isb-cgc-02-0001.job_Lbt8AICMkqRCZHkIWZgNLoFV2SCG.US' [\]  3s
    Running job 'isb-cgc-02-0001.job_Lbt8AICMkqRCZHkIWZgNLoFV2SCG.US' [|]  3s
    ## Complete
    ## Billed: 20.97 MB
    ## Downloading 327 rows in 1 pages.

``` r
head(early_onset_breast_cancer, 5)
```

    ## # A tibble: 5 x 2
    ## # Groups:   case_barcode, category [5]
    ##   category case_barcode
    ##   <chr>    <chr>       
    ## 1 <NA>     TCGA-A2-A0D3
    ## 2 <NA>     TCGA-A8-A09A
    ## 3 <NA>     TCGA-AC-A23E
    ## 4 <NA>     TCGA-AR-A24U
    ## 5 <NA>     TCGA-AR-A2LJ
