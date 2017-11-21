# OAFilter
Trim and filter reads by the Overall accuracy



### Overall Accuracy

Phred Score **Q** is a famous term used for evaluating data quality:

$$Q=-10\log_{10}P$$

Where $P$ represent the error discover rate.

Since we are not dealt with assembly task, we focus on the whole quality of a continuously framgent. So we use $Overall Accuray$ to assess it:

$$OA=\prod_{i=1}^{L} (1-P_i)$$

Where $L$ represent the fragment end position, $i$ is the specific base position and $P_i$ is the error discover rate. 

Natually, $OA$ is the Conditional probability of a whole fragment with **all bases are corrected identified**.

In case we gererally obtained the $Q$ score, so we can use the following formula:

$$OA= $\prod_{i=1}^{L} (1-10^{-Q_i/10})$$

The relation between $Q$ and $OA$ can be show like below table:

If each base have the same quality, we get:

| Phrd Q | error P | Accuracy | 30bp_OA  | 50bp_OA  | 100bp_OA    | note                                     |
| ------ | ------- | -------- | -------- | -------- | ----------- | ---------------------------------------- |
| 50     | 0.00001 | 0.99999  | 0.9997   | 0.9995   | 0.999000495 |                                          |
| 40     | 0.0001  | 0.9999   | 0.997004 | 0.995012 | 0.990049339 | the Hiseq 2000 ceiling                   |
| 37     | 0.0002  | 0.9998   | 0.994031 | 0.990072 | 0.980243162 | the ZEBRA-500 ceiing(mid-2017)           |
| 30     | 0.001   | 0.999    | 0.970431 | 0.951206 | 0.904792147 | Often regarded as 'High quality fragment' |
| 20     | 0.01    | 0.99     | 0.7397   | 0.605006 | 0.366032341 |                                          |
| 10     | 0.1     | 0.9      | 0.042391 | 0.005154 | 2.65614E-05 |                                          |

To summary, the OA is effected by fragment length, and much sensitive to base quality. So this algrathm will try to trim both tails of a read to obtain the fragment area with highest OA score.

# E.g

#### 1. Overall Accuracy algrithm :
```{r, echo=TRUE}
OA0_fun <- function(x){
  x <- as.numeric(x)
  oa <- rep(1,length(x))
  acc <- 1-10^(-x[1]/10)
  oa[1] <- acc
  for(i in 2:length(x)){
    acc <- 1-10^(-x[i]/10)
    oa[i] <- oa[i-1] * acc
  }
  return(oa)
}
```

#### 2. Overall Accuracy with one worst quality base ignored
```{r, echo=TRUE}
OA1_fun <- function(x){
  x <- as.numeric(x)
  oa <- rep(1,length(x))
  acc <- 1-10^(-x[1]/10)
  ma <- acc
  for(i in 2:length(x)){
    acc <- 1-10^(-x[i]/10)
    if(acc < ma){
      oa[i] <- oa[i-1] * ma
      ma <- acc
    }else{
      oa[i] <- oa[i-1] * acc
    }
  }
  return(oa)
}
```
#### 3. Overall Accuracy with two worst quality bases ignored
```{r, echo=TRUE}
OA2_fun <- function(x){
  x <- as.numeric(x)
  oa <- rep(1,length(x))
  ma <- c(1-10^(-x[1]/10),1-10^(-x[2]/10))
  for(i in 3:length(x)){
    acc <- 1-10^(-x[i]/10)
    ma <- sort(c(ma[1:2],acc))
    oa[i] <- oa[i-1] * ma[3]
  }
  return(oa)
}
```