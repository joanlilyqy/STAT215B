#### Ying Qiao
#### SID: 21412301
#### STAT 215B, Assignment 3, 3/1/2013

rm(list=ls())

#### Detective Work
## N : total subject counts for various subgroups/margins
unsub <- read.table("part6_907.txt", header=T)
head(unsub)

# --- from codebook (for reference here) ---
# sempl:  0-unemployed;  1-employed;
# act (treatment assigned | treatment delivered):
#         1-arr|arr;     2-non|non;    3-non|arr;   4-arr|non;

# initial look
with(unsub, table(sempl)) #table(unsub$sempl)["0"]
table(unsub$act)

# follow "intension-to-treat" to regroup "act"
# trass: 0-no arrest;    1-arrest;
unsub$trass <- as.numeric(unsub$act == 1 | unsub$act == 4)
with(unsub, table(trass))

# create N table and margins
tbl <- with(unsub, table(sempl, trass))
N.unempl <- sum(tbl[1, ])
N.empl   <- sum(tbl[2, ])
N.noarr  <- sum(tbl[ ,1])
N.arr    <- sum(tbl[ ,2])
N <- sum(sum(tbl))

# report the rate of unemployment among the subjects
N.unempl / N  # PH gives 29%


## n : number of re-abusing subjects
# --- from PH Figure 1.---
PH.per <- cbind(c(7.1,12.3),c(16.7,6.2)) / 100 # following the same 0-1 in tbl

# create n table and margins
recid <- round(tbl*PH.per)
n.unempl <- sum(recid[1, ])
n.empl   <- sum(recid[2, ])
n.noarr  <- sum(recid[ ,1])
n.arr    <- sum(recid[ ,2])
n <- sum(sum(recid))

n.noarr / N.noarr # PH gives 10.6%
n.arr / N.arr     # PH gives 9.0%


#### Statistical Work
# Fisher's exact test


# 1. employed; recid ~ arrest



# 2. unemployed; recid ~ arrest



# 3. recid ~ arrest

