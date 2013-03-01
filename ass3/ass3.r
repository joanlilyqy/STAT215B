#### Ying Qiao
#### SID: 21412301
#### STAT 215B, Assignment 3, 3/1/2013

rm(list=ls())

#### Detective Work
#### N : total subject counts for various subgroups/margins
unsub <- read.table("part6_907.txt", header=T)
head(unsub)
with(unsub, table(act))

# --- from codebook | for reference here ---
# sempl:  0-unemployed;  1-employed;
# act:  treatment assigned | treatment delivered
#         1-arr|arr;     2-non|non;    3-non|arr;   4-arr|non;

# follow "intension-to-treat" to regroup "act"
# trass: 0-no arrest;    1-arrest;
unsub$trass <- as.numeric(unsub$act == 1 | unsub$act == 4)

# create N table, margins can be derived
N.tbl <- with(unsub, table(sempl, trass))
N.tbl

# report the rate of unemployment among the subjects
rowSums(N.tbl)[1] / sum(sum(N.tbl)) * 100   # PH gives 29% for unemployment rate


#### n : number of re-abusing subjects
# --- from PH Figure 1.---
PH.per <- cbind(c(7.1,12.3),c(16.7,6.2)) / 100 # following the same 0/1 in N.tbl

# Similarly, create n table, margins can be derived
n.tbl <- round(N.tbl * PH.per)
n.tbl

# report the rate of re-abuse among non/-arrestees
colSums(n.tbl)[1] / colSums(N.tbl)[1] * 100 # PH gives 10.6% for non-arrestees
colSums(n.tbl)[2] / colSums(N.tbl)[2] * 100 # PH gives 9.0%  for arrestees


#### Statistical Work
# a - Fisher's exact test
# b - two sample test of equal binomial proportions
#     z = (p1-p2)/sqrt(p(1-p)(1/n1+1/n2))
#     where p=(n1p1+n2p2)/(n1+n2)
two.binom.test <- function (tbl) {
  # 2-by-2 table similar to Fisher test
  p <- tbl[1, ] / colSums(tbl)
  p.hat <- rowSums(tbl)[1] / sum(sum(tbl))
  del.hat <- p[1] - p[2]
  var.hat <- p.hat*(1-p.hat)*sum(1/colSums(tbl))
  z.hat <- (del.hat) / sqrt(var.hat)
  return (z.hat)
}


# 1. employed; re-abuse ~ arrest
# among employed subjects
# H0: the occurrence of a subsequent assault is the same between the two groups
# HA: the re-abuse rate is larger in the non-arrest group than the treatment group
cl1 <- matrix(c(n.tbl[2, ], (N.tbl-n.tbl)[2, ]), nrow=2, byrow=T,
              dimnames = list(c("reabuse","non-recid"), c("non-arrest","arrest")))
cl1
fisher.test(cl1, alternative="greater")
1 - pnorm(two.binom.test(cl1))


# 2. unemployed; re-abuse ~ arrest
# among unemployed subjects
# H0: the occurrence of a subsequent assault is the same between the two groups
# HA: the re-abuse rate is smaller in the non-arrest group than the treatment group
cl2 <- matrix(c(n.tbl[1, ], (N.tbl-n.tbl)[1, ]), nrow=2, byrow=T,
              dimnames = list(c("reabuse","non-recid"), c("non-arrest","arrest")))
cl2
fisher.test(cl2, alternative="less")
pnorm(two.binom.test(cl2))


# 3. all; re-abuse ~ arrest
# among all subjects
# H0: the occurrence of a subsequent assault is the same between the two groups
# HA: the re-abuse rate is different in the non-arrest group from the treatment group
cl3 <- cl1 + cl2
cl3
fisher.test(cl3, alternative="two.sided")
pnorm(two.binom.test(cl3)) / 2



