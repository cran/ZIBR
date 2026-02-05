## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ZIBR)
library(dplyr)
library(nlme)
library(tibble)

set.seed(19683)

## ----echo=FALSE, out.width="75%", fig.cap="Details of the statistical model"----
knitr::include_graphics("zibr_stat_model.png")

## -----------------------------------------------------------------------------
sim <- simulate_zero_inflated_beta_random_effect_data(
  subject_n = 100, time_n = 5,
  X = as.matrix(c(rep(0, 50 * 5), rep(1, 50 * 5))),
  Z = as.matrix(c(rep(0, 50 * 5), rep(1, 50 * 5))),
  alpha = as.matrix(c(-0.5, 1)),
  beta = as.matrix(c(-0.5, 0.5)),
  s1 = 1, s2 = 0.8,
  v = 5,
  sim_seed = 100
)
names(sim)

## -----------------------------------------------------------------------------
zibr_fit <- zibr(
  logistic_cov = sim$X, beta_cov = sim$Z, Y = sim$Y,
  subject_ind = sim$subject_ind, time_ind = sim$time_ind
)
zibr_fit

## -----------------------------------------------------------------------------
zibr_fit <- zibr(
  logistic_cov = sim$X, beta_cov = sim$X, Y = sim$Y,
  subject_ind = sim$subject_ind, time_ind = sim$time_ind
)
zibr_fit

## -----------------------------------------------------------------------------
head(ibd, n = 20)

## -----------------------------------------------------------------------------
zibr_fit <- zibr(
  logistic_cov = data.frame(Treatment = ibd$Treatment),
  beta_cov = data.frame(Treatment = ibd$Treatment),
  Y = ibd$Abundance, subject_ind = ibd$Subject,
  time_ind = ibd$Time
)
zibr_fit

## -----------------------------------------------------------------------------
### Load the raw data
PLEASE.file <- "https://raw.githubusercontent.com/chvlyl/PLEASE/master/1_Data/Raw_Data/MetaPhlAn/PLEASE/G_Remove_unclassfied_Renormalized_Merge_Rel_MetaPhlAn_Result.xls" # nolint
PLEASE.raw <- read.table(PLEASE.file,
  sep = "\t", header = TRUE, row.names = 1,
  check.names = FALSE, stringsAsFactors = FALSE
)
taxa.raw <- t(PLEASE.raw)

### Make sure you load the data correctly
cat("samples", "taxa", dim(taxa.raw), "\n")
taxa.raw[1:3, 1:3]

## -----------------------------------------------------------------------------
### Load total non-human read counts
human.read.file <- "https://raw.githubusercontent.com/chvlyl/PLEASE/master/1_Data/Raw_Data/MetaPhlAn/Human_Reads/please_combo_human_reads.xls" # nolint
human.read <- read.table(human.read.file,
  sep = "\t", header = TRUE,
  row.names = 1, stringsAsFactors = FALSE
)

### Filter low depth samples (low non human reads)
low.depth.samples <- subset(human.read, NonHumanReads < 10000)
low.depth.samples[, 1:5]

### Delete these samples from PLEASE data.
rownames(taxa.raw)[which(rownames(taxa.raw) %in% rownames(low.depth.samples))]

### Before deletion
dim(taxa.raw)

### After deletion
taxa.raw <- taxa.raw[-which(rownames(taxa.raw) %in% rownames(low.depth.samples)), ]
dim(taxa.raw)


### Filter low abundant bacterial data
filter.index1 <- apply(taxa.raw, 2, function(X) {
  sum(X > 0) > 0.4 * length(X)
})
filter.index2 <- apply(taxa.raw, 2, function(X) {
  quantile(X, 0.9) > 1
})
taxa.filter <- taxa.raw[, filter.index1 & filter.index2]
taxa.filter <- 100 * sweep(taxa.filter, 1, rowSums(taxa.filter), FUN = "/")
cat("after filter:", "samples", "taxa", dim(taxa.filter), "\n")
cat(colnames(taxa.filter), "\n")
head(rowSums(taxa.filter))

###
taxa.data <- taxa.filter
dim(taxa.data)

## -----------------------------------------------------------------------------
#### Load sample information
sample.info.file <- "https://raw.githubusercontent.com/chvlyl/PLEASE/master/1_Data/Processed_Data/Sample_Information/2015_02_13_Processed_Sample_Information.csv" # nolint
sample.info <- read.csv(sample.info.file, row.names = 1)

## -----------------------------------------------------------------------------
#### create covariates,
#### Time, antiTNF+EEN
reg.cov <-
  data.frame(Sample = rownames(taxa.data), stringsAsFactors = FALSE) %>%
  left_join(rownames_to_column(sample.info, var = "Sample"), by = "Sample") %>%
  dplyr::filter(Treatment.Specific != "PEN") %>%
  dplyr::select(Sample, Time, Subject, Response, Treatment.Specific) %>%
  group_by(Subject) %>%
  summarise(count = n()) %>%
  dplyr::filter(count == 4) %>%
  dplyr::select(Subject) %>%
  left_join(rownames_to_column(sample.info, var = "Sample"), by = "Subject") %>%
  mutate(Treat = ifelse(Treatment.Specific == "antiTNF", 1, 0)) %>%
  dplyr::select(Sample, Subject, Time, Response, Treat) %>%
  dplyr::mutate(Subject = paste("S", Subject, sep = "")) %>%
  dplyr::mutate(Time = ifelse(Time == "1", 0,
                              ifelse(Time == "2", 1,
                                     ifelse(Time == "3", 4,
                                            ifelse(Time == "4", 8, NA))))) %>%
  dplyr::mutate(Time.X.Treatment = Time * Treat) %>%
  as.data.frame()

### take out first time point
reg.cov.t1 <- subset(reg.cov, Time == 0)
rownames(reg.cov.t1) <- reg.cov.t1$Subject
reg.cov.t234 <- subset(reg.cov, Time != 0)
reg.cov.t234 <- data.frame(
  baseline.sample = reg.cov.t1[reg.cov.t234$Subject, "Sample"],
  baseline.subject = reg.cov.t1[reg.cov.t234$Subject, "Subject"],
  reg.cov.t234,
  stringsAsFactors = FALSE
)

#### Fit ZIBR model on the real data
spe.all <- colnames(taxa.data)
p.species.list.zibr <- list()
p.species.list.lme <- list()
for (spe in spe.all) {
  # spe = 'g__Collinsella'
  ###### create covariates
  X <- data.frame(
    Baseline = taxa.data[reg.cov.t234$baseline.sample, spe] / 100,
    # reg.cov.t234[,c('log.days','Delivery','Delivery.X.log.days')]
    reg.cov.t234[, c("Time", "Treat")]
  )
  rownames(X) <- reg.cov.t234$Sample
  Z <- X
  subject.ind <- reg.cov.t234$Subject
  time.ind <- reg.cov.t234$Time
  ####
  # cat(spe,'\n')
  Y <- taxa.data[reg.cov.t234$Sample, spe] / 100
  # cat('Zeros/All',sum(Y==0),'/',length(Y),'\n')
  ####################
  ## fit linear random effect model with arcsin square transformation on Y
  tdata <- data.frame(Y.tran = asin(sqrt(Y)), X, SID = subject.ind)
  lme.fit <- nlme::lme(Y.tran ~ Baseline + Time + Treat, random = ~ 1 | SID, data = tdata)
  coef.mat <- summary(lme.fit)$tTable[-1, c(1, 5)]
  p.species.list.lme[[spe]] <- coef.mat[, 2]
  ###################
  if (sum(Y > 0) < 10 || sum(Y == 0) < 10 || max(Y) < 0.01) {
    p.species.list.zibr[[spe]] <- p.species.list.lme[[spe]]
    next
  } else {
    est <- zibr(
      logistic_cov = X, beta_cov = Z, Y = Y,
      subject_ind = subject.ind,
      time_ind = time.ind,
      quad_n = 30, verbose = TRUE
    )
    p.species.list.zibr[[spe]] <- est$joint.p
  }
  # break
}


#### adjust p values
p.species.zibr <- t(as.data.frame(p.species.list.zibr))
p.species.zibr.adj <-
  rownames_to_column(as.data.frame(p.species.zibr), var = "Species") %>% mutate(across(-Species, ~p.adjust(., "fdr")))

# write.csv(p.species.zibr.adj,file=paste('4_Results/Real_Data_Estimation_Results_antiTNF_EEN_ZIBR.csv',sep=''))


p.species.lme <- t(as.data.frame(p.species.list.lme))
p.species.lme.adj <-
  rownames_to_column(as.data.frame(p.species.lme), var = "Species") %>% mutate(across(-Species, ~p.adjust(., "fdr")))

# write.csv(p.species.lme.adj,file=paste('4_Results/Real_Data_Estimation_Results_antiTNF_EEN_LME.csv',sep=''))

####
intersect(
  p.species.lme.adj[p.species.lme.adj$Treat < 0.05, "Species"],
  p.species.zibr.adj[p.species.zibr.adj$Treat < 0.05, "Species"]
)
setdiff(
  p.species.lme.adj[p.species.lme.adj$Treat < 0.05, "Species"],
  p.species.zibr.adj[p.species.zibr.adj$Treat < 0.05, "Species"]
)
setdiff(
  p.species.zibr.adj[p.species.zibr.adj$Treat < 0.05, "Species"],
  p.species.lme.adj[p.species.lme.adj$Treat < 0.05, "Species"]
)

intersect(
  p.species.lme.adj[p.species.lme.adj$Time < 0.05, "Species"],
  p.species.zibr.adj[p.species.zibr.adj$Time < 0.05, "Species"]
)
setdiff(
  p.species.lme.adj[p.species.lme.adj$Time < 0.05, "Species"],
  p.species.zibr.adj[p.species.zibr.adj$Time < 0.05, "Species"]
)
setdiff(
  p.species.lme.adj[p.species.lme.adj$Time < 0.05, "Species"],
  p.species.zibr.adj[p.species.zibr.adj$Time < 0.05, "Species"]
)

