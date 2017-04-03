Pharmacy Visit Pattern Mining
================
azure31
April 3, 2017

I have an anonymized transaction database of customer visits to different pharmacies for purchase of drug X and I want to investigate the pharmacy visit sequence for patients who are frequent shoppers. Since drug X has potential for abuse/overuse, I would like to see if these customers exhibit any patterns that can be used to identify such behaviour apriori.

``` r
#Load dataset
pattern= read.csv("pattern_mining_hash.csv")
pattern$MBR_SYS_ID= as.character(pattern$MBR_SYS_ID)
pattern$NABP_NBR= as.character(pattern$NABP_NBR)
head(pattern)
```

    ##   MBR_SYS_ID FILL_DT_SYS_ID NABP_NBR
    ## 1   e492c38f           7291 8149c4ee
    ## 2   e492c38f           7298 5dd25528
    ## 3   e492c38f           7312 08d19759
    ## 4   e492c38f           7323 8149c4ee
    ## 5   e492c38f           7357 4f82c384
    ## 6   e492c38f           7357 4f82c384

Here, we have three variables, member\_id, fill\_dt\_sys\_id encodes the date of pharmacy visit and the hashed pharmacy NABP number.

``` r
library(data.table)
pattern= data.table(pattern)
pattern=pattern[order(MBR_SYS_ID, FILL_DT_SYS_ID, NABP_NBR)]

#remove rows where multiple pharmacies visited on same day
pattern= pattern[, head(.SD, 1), by=c("MBR_SYS_ID", "FILL_DT_SYS_ID")]
```

Here, we are concerned with the visit pattern and not the pharmacies themselves so two members with visit pattern ababab and cdcdcd should be considered similar in terms of visit behavior. In order to incorporate this into the sequence we encode it in the following manner (assuming window size 10 \*) -

|0 | member visits pharmacy previously not visited |1 | members visits pharmacy of lag 1 (ie previously visited pharmacy) |2 | member visits pharmacy of lag 2 |3 | member visits pharmacy of lag3 |4 | member visits pharmacy of lag4 etc

Therefore, sequences ababab and cdcdcd are both encoded to 002222. Using this encoding, we apply sequential pattern mining to find frequent sub-sequences with support &gt;=0.2 (ie present in atleast 20% of all visit sequences)

\*In the data set I have, the average number of pharmacy visits for a customer over a year is around 10.

We begin by encoding our transaction dataset such that we capture the visit sequence in the way described above. Here the seq variable is used to encode the visit sequence.

``` r
#Function to convert pharmacy visit pattern to sequence

pattern$seq= rep(0, nrow(pattern))
sequence <-function(id){
  
  df= pattern[pattern$MBR_SYS_ID==id,]
  
  df$lag1= ifelse(shift(df$NABP_NBR, 1, "lag")==df$NABP_NBR, 1, 0)
  df$lag2= ifelse(shift(df$NABP_NBR, 2, "lag")==df$NABP_NBR, 1, 0)
  df$lag3= ifelse(shift(df$NABP_NBR, 3, "lag")==df$NABP_NBR, 1, 0)
  df$lag4= ifelse(shift(df$NABP_NBR, 4, "lag")==df$NABP_NBR, 1, 0)
  df$lag5= ifelse(shift(df$NABP_NBR, 5, "lag")==df$NABP_NBR, 1, 0)
  df$lag6= ifelse(shift(df$NABP_NBR, 6, "lag")==df$NABP_NBR, 1, 0)
  df$lag7= ifelse(shift(df$NABP_NBR, 7, "lag")==df$NABP_NBR, 1, 0)
  df$lag8= ifelse(shift(df$NABP_NBR, 8, "lag")==df$NABP_NBR, 1, 0)
  df$lag9= ifelse(shift(df$NABP_NBR, 9, "lag")==df$NABP_NBR, 1, 0)
  df$lag10= ifelse(shift(df$NABP_NBR, 10, "lag")==df$NABP_NBR, 1, 0)
  df$lag0=as.numeric(with(df, !(lag1 | lag2| lag3 | lag4| lag5| lag6 | lag7 | lag8 | lag9 | lag10)))
  df1= df[, c(15, 5:14)]
  df1$seq= sapply(1:nrow(df1), function(x){which.max(df1[x])})
  df1$seq=df1$seq-1
  return(df1$seq)
  
}


ids= unique(pattern$MBR_SYS_ID)
seq= lapply(ids, sequence)

#apply function to our data
for( i in ids){
  pattern$seq[pattern$MBR_SYS_ID==i]= sequence(i)
}

pattern=pattern[order(MBR_SYS_ID, FILL_DT_SYS_ID)]
head(pattern)
```

    ##    MBR_SYS_ID FILL_DT_SYS_ID NABP_NBR seq
    ## 1:   0001d1d5           7489 17163ddb   0
    ## 2:   0001d1d5           7502 a23c85b7   0
    ## 3:   0001d1d5           7509 a23c85b7   1
    ## 4:   0001d1d5           7515 172a84a8   0
    ## 5:   0001d1d5           7523 a23c85b7   2
    ## 6:   0001d1d5           7529 4de1491c   0

We will use the arulesSequence package in R to perform the sequential pattern mining. As a pre-processing step, we will convert our dataset to conform with the data format required by the package- Each patient's visit sequence will be captured by the tuple (sequenceID, eventID, items) where sequenceID corresponds to the patientID, eventID to the visit\_dt and items to the basket of pharmacies visited. Now we will convert data to conform to arulesSequences data format

``` r
MBR_SYS_ID= sort(unique(pattern$MBR_SYS_ID))
sequenceID= 1:length(MBR_SYS_ID)
df_seqid= data.frame(MBR_SYS_ID, sequenceID)

pattern= merge(pattern, df_seqid, by= "MBR_SYS_ID", all.x=T )

library(dplyr)
```

    ## -------------------------------------------------------------------------

    ## data.table + dplyr code now lives in dtplyr.
    ## Please library(dtplyr)!

    ## -------------------------------------------------------------------------

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:data.table':
    ## 
    ##     between, first, last

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
pattern= pattern %>% 
  group_by(MBR_SYS_ID) %>% 
  mutate(eventID = dense_rank(FILL_DT_SYS_ID))

pattern$MBR_SYS_ID=NULL
pattern$FILL_DT_SYS_ID=NULL
pattern$NABP_NBR=NULL
pattern$size= rep(1, nrow(pattern))
pattern= pattern[, c(2, 3, 4, 1)]

head(pattern)
```

    ## # A tibble: 6 × 4
    ##   sequenceID eventID  size   seq
    ##        <int>   <int> <dbl> <dbl>
    ## 1          1       1     1     0
    ## 2          1       2     1     0
    ## 3          1       3     1     1
    ## 4          1       4     1     0
    ## 5          1       5     1     2
    ## 6          1       6     1     0

``` r
#save down results
write.table(pattern, "pattern_shoppers.csv", 
            row.names = FALSE, 
            col.names = FALSE,
            sep = ' ',
            quote = FALSE)
```

We will use the cspade algorithm in arulesSequences to perform pattern mining. Using support 20% and adjusting maxgap=1, we want to consider only consecutive series.

``` r
library(arulesSequences)
```

    ## Loading required package: arules

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'arules'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     recode

    ## The following objects are masked from 'package:base':
    ## 
    ##     abbreviate, write

``` r
pattern_tran <- read_baskets("pattern_shoppers.csv", info = c("sequenceID","eventID","size"))
head(as(pattern_tran, "data.frame"))
```

    ##   items sequenceID eventID size
    ## 1   {0}          1       1    1
    ## 2   {0}          1       2    1
    ## 3   {1}          1       3    1
    ## 4   {0}          1       4    1
    ## 5   {2}          1       5    1
    ## 6   {0}          1       6    1

``` r
s1 <- cspade(pattern_tran, 
             parameter = list(support = 0.4, maxlen = 10, maxgap = 1), 
             control = list(verbose = TRUE, summary = TRUE, bfstype = TRUE))
```

    ## 
    ## parameter specification:
    ## support : 0.4
    ## maxsize :  10
    ## maxlen  :  10
    ## maxgap  :   1
    ## 
    ## algorithmic control:
    ## bfstype  :  TRUE
    ## verbose  :  TRUE
    ## summary  :  TRUE
    ## tidLists : FALSE
    ## 
    ## preprocessing ... 1 partition(s), 1.51 MB [0.68s]
    ## mining transactions ... 0 MB [0.53s]
    ## reading sequences ... [0.02s]
    ## 
    ## total elapsed time: 1.23s

``` r
summary(s1)
```

    ## set of 9 sequences with
    ## 
    ## most frequent items:
    ##       0       1       2 (Other) 
    ##       6       4       2       2 
    ## 
    ## most frequent elements:
    ##     {0}     {1}     {2} (Other) 
    ##       6       4       2       2 
    ## 
    ## element (sequence) size distribution:
    ## sizes
    ## 1 2 3 
    ## 3 5 1 
    ## 
    ## sequence length distribution:
    ## lengths
    ## 1 2 3 
    ## 3 5 1 
    ## 
    ## summary of quality measures:
    ##     support      
    ##  Min.   :0.4696  
    ##  1st Qu.:0.5328  
    ##  Median :0.6340  
    ##  Mean   :0.6729  
    ##  3rd Qu.:0.7797  
    ##  Max.   :1.0000  
    ## 
    ## includes transaction ID lists: FALSE 
    ## 
    ## mining info:
    ##          data ntransactions nsequences support
    ##  pattern_tran         68207       7091     0.4

``` r
rules=as(s1, "data.frame")

#Print all sequences
print(rules$sequence)
```

    ## [1] <{0}>         <{1}>         <{2}>         <{0},{2}>     <{0},{1}>    
    ## [6] <{1},{1}>     <{0},{0}>     <{1},{0}>     <{0},{0},{0}>
    ## 9 Levels: <{0},{0},{0}> <{0},{0}> <{0},{1}> <{0},{2}> <{0}> ... <{2}>
