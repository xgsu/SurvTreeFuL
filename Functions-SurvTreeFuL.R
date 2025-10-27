
# ###########################################
# FUNCTIONS FOR SURVIVAL TREEFUL
# ###########################################


options(warn=-1)  # TURN OFF WARNINGS GLODBALLY
library(tidyverse)
library(survival)
# library(survminer) 
require(rpart)
require(MASS)
# library(maxstat)     # ORDINAL VARIABLES ONLY
library(rpart.plot)
library("viridis")
library(glmnet)
library(survAUC)
library(gmp)
# library(plyr)
# remove.packages("plyr")
select <- dplyr::select
rowSums <- base::rowSums; colSums <- base::colSums




factor2numeric <- function(x) {as.numeric(levels(x))[x]} 
expit <- function(x) (tanh(x/2)+1)/2  	# CORRECT & VERIFIED
remove.rows.AllZeros <- function(M, threshold=0) return(M[rowSums(M)>threshold,])
remove.columns.AllZeros <- function(M, threshold=0) return(M[,colSums(M)>threshold]) 



# CONTROL OPTION FOR rpart TO OBTAIN ONE SINGLE SPLIT
ctr0.rpart <- rpart.control(minsplit = 3, minbucket=1, cp = -1, 
                          maxcompete = 0, maxsurrogate = 0, xval = 1,
                          surrogatestyle = 0, maxdepth = 1)

# CONTROL OPTIONS FOR FUNCTION best.cut()
ctr0.bestcut <- list(
  method0 = "default", # OPTION TO ENFORCE GREEDY SEARCH
  min.ndeath=5,      # MINIMUM # DEATHS FOR GS
  d= 0.03,            # SET SEARCH RANGE FOR SSS
  a0 = 50, epi = .Machine$double.eps^0.25, # PARAMETERS IN SSS
  threshold.ncuts=20) # THRESHOLD # CUTS TO CHOOSE BT GS & SSS

# =============================================================
# FUNCTION Replace() REPACES A SET OF VALUES WITH ANOTHER SET
# =============================================================

Replace <- function(x, old.values, new.values){
  if (length(old.values) != length(new.values)) stop("The legnth of old.values is not the same as that of new.values.")
  replace(
    x,
    x %in% old.values,
    new.values[na.omit(match(x, old.values))]
  )
}





# =============================================================
# FUNCTION rpart.split1rule() EXTRACTS SPLIT FROM rpart OBJECT
# =============================================================
# THE SPLITTING VARIABLE AND CUTPOINT FOR THE FIRST SPLIT
rpart.split1rule <- function(split.rpart){
  xsplit.rpart <- split.rpart$frame[1, 1]; # xsplit.rpart
  cutpoint <- rpart.plot::rpart.rules(split.rpart, roundint =FALSE)[1, 5]
  list(xsplit=xsplit.rpart, cutpoint=cutpoint)
}


# ------------------------------------------------------------------
# FUNCTION int2binary CONVERTS N (NATURAL NUMBER) INTO BINARY
# -------------------------------------------------------------------
N2binary <- function(x) {
  y <- rev(as.integer(intToBits(x))) 
  id.rm <- seq_len(match(1,y,length(y))-1)
  x <- paste(y[-id.rm], collapse = "")
  return(x)
}
# N2binary(123)



# =======================================================
# FUNCTION Obj.LogRank() - THE OBJECTIVE FUNCTION IN SSS
# =======================================================
# LOGRANK IS THE SCORE TEST IN coxph. NO NEED TO FIT ALTERNATIVE MODEL. 
ctr0.coxph <-  coxph.control(eps = 1e-02, toler.chol = 1e-3,
                             iter.max = 0, toler.inf = 1e-2, outer.max =1, timefix=TRUE)
Obj.LogRank <- function(c, a= 50, Time, status, x) 
{
  z <- expit(a*(x-c));
  # z <- sign(x < c)
  logrank <- coxph(Surv(Time, status) ~ z, y=FALSE, control=ctr0.coxph)$score
  return(logrank)
}

# ===============================================================
# FUNCTION LogRank() - SELF-WRITTEN FOR LOGRANK TEST STATISTIC
# ===============================================================
# EVEN FASTER THAN survdiff()
LogRank <- function(Time, status, z, min.event.childnode=2) {
  logrank <- NA
  n.L1 <- sum(z==0 & status==1); n.R1 <- sum(z==1 & status==1);
  # print(cbind(n=length(Time), n.R=sum(z==1), n.1=sum(status==1), n.L1, n.R1))
  if (min(n.L1, n.R1) >= min.event.childnode) {
    # Time.death <- Time[status==1]
    Risk <- outer(Time, Time, FUN=">=")
    R11 <- colSums(Risk[z==1, (z==1) & (status==1)])
    R12 <- colSums(Risk[z==1, (z==0) & (status==1)])
    R21 <- colSums(Risk[z==0, (z==1) & (status==1)])
    R22 <- colSums(Risk[z==0, (z==0) & (status==1)])
    Diff <- sum(R21/(R11+R21))-sum(R12/(R12+R22))
    Var <- sum(R11*R21/(R11+R21)^2) + sum(R12*R22/(R12+R22)^2)
    logrank <- Diff^2/Var
    # p.value <- pchisq(logrank, df=1, lower.tail =FALSE)
  }
  return(logrank)
}


# --------------------------------
# FUNCTION rdat() GENERATES DATA 
# --------------------------------
rdat <- function(n=500, beta=c(-2, 0, 0), cutoff=.5, digits=2, 
                 model=c("linear", "MARS1", "tree"), 
                 p0=2, add.nominal=TRUE, n.LETTERS = 5, 
                 icensor=1, details=FALSE) 
{
  # GENERATE COVARIATES
  x1 <- sample(x=0:1, size=n, replace=TRUE)
  x2 <- round(runif(n=n, min=0, max=1), digits =digits)
  X <- cbind(matrix(sample(x=0:1, size=n*p0, replace=TRUE), nrow = n, ncol=p0),
       matrix(round(runif(n=n*p0, min=0, max=1), digits = digits), nrow = n, ncol=p0))
  if (add.nominal) x3 <- sample(LETTERS[1:n.LETTERS], size=n, replace=TRUE)
  if (model=="linear") rate <- exp(beta[1] + beta[2]*x2 + beta[3]*X[, p0+1])   
  else if (model=="tree") rate <- exp(beta[1] + beta[2]*x1 + beta[3]*sign(x2<=cutoff))
  else if (model=="tree1") rate <- exp(beta[1] + beta[2]*sign(sin(x2*6*pi) >=0))
  else if (model=="tree2") rate <- exp(beta[1] + beta[2]*x1*sign(0.25 <= x2 & x2 <= 0.75))
  # else if (model=="tree3") rate <- exp(beta[1] + beta[2]*sign(is.element(x3, LETTERS[1:3]))*sign(x2 < 0.5))
  else if (model=="KAN") rate <- exp(beta[1] + beta[2]*sin(5*pi*x2^2) + beta[3]*sin(5*pi*X[, p0+1]*X[, p0+1]))
  # rate <- exp(-1 + 2*sin(2*pi*x2^2) + 2*sin(2*pi*X[, p0+1]*X[, p0+1]))
  else print("Error in specification of model=. Check!")
  xobs <- rexp(n, rate=rate); 
  if (details) print(mean(xobs))
  #### Generate observed failure time and status
  if (icensor==0) status <- 1
  else {     
    censor <- rexp(n, rate=rate)
    status <- sign(xobs <= censor*icensor)
    xobs <- pmin(xobs, censor*icensor) 
  }
  if (details) print(mean(status))
  ##### Output
  dat <- data.frame(id=1:n, time=xobs, status=status, x1=x1, x2=x2, X)
  if (add.nominal) dat <- data.frame(dat, x3=x3)
  dat <- dat %>% mutate_if(is.character, as.factor)   
  names(dat) <- c("id", "time", "status", paste("x", 1:(NCOL(dat)-3), sep=""))
  return(dat)
}  

# NON-PH TREE MODEL
rdat.AFT <- function(n=500, beta=c(-1, 1, 3), cutoff=.5, digits=2, 
                 p0=2, add.nominal=TRUE, n.LETTERS = 5, 
                 icensor=1, details=FALSE) 
{
  # GENERATE COVARIATES
  x1 <- sample(x=0:1, size=n, replace=TRUE)
  x2 <- round(runif(n=n, min=0, max=1), digits =digits)
  X <- cbind(matrix(sample(x=0:1, size=n*p0, replace=TRUE), nrow = n, ncol=p0),
             matrix(round(runif(n=n*p0, min=0, max=1), digits = digits), nrow = n, ncol=p0))
  if (add.nominal) x3 <- sample(LETTERS[1:n.LETTERS], size=n, replace=TRUE)
  xobs <- exp(beta[1] + beta[2]*x1 + beta[3]*sign(x2<=cutoff) + rlogis(n, location = 0, scale = 1))
  if (icensor==0) status <- 1
  else {     
    censor <- exp(beta[1] + beta[2]*x1 + beta[3]*sign(x2<=cutoff) + rlogis(n, location = 0, scale = 1))
    status <- sign(xobs <= censor*icensor)
    xobs <- pmin(xobs, censor*icensor) 
  }
  if (details) print(mean(status))
  ##### Output
  dat <- data.frame(id=1:n, time=xobs, status=status, x1=x1, x2=x2, X)
  if (add.nominal) dat <- data.frame(dat, x3=x3)
  dat <- dat %>% mutate_if(is.character, as.factor)   
  names(dat) <- c("id", "time", "status", paste("x", 1:(NCOL(dat)-3), sep=""))
  return(dat)
}  






# =========================================================================
# FUNCTION best.cut FINDS THE BEST CUTOFF POINT OF A CONTINUOUS VARIABLE
# =========================================================================

best.cut <- function(time, status, x, 
                     control.bestcut=list(
                      method0 = "default", # OPTION TO ENFORCE GREEDY SEARCH
                      min.ndeath=3,      # MINIMUM # DEATHS FOR GS
                      d= 0.03,            # SET SEARCH RANGE FOR SSS
                      a0 = 50, epi = .Machine$double.eps^0.25, # PARAMETERS IN SSS
                      threshold.ncuts=20) # THRESHOLD # CUTS TO CHOOSE BT GS & SSS
                     ) 
{
  bcut <- NA; score.max <- -1e10; 
  # print(x)
  temp <- sort(unique(x)); 
  if (is.character(x)) x <- as.factor(x)
  if (is.factor(x) && length(temp) >1)  {     # HANDLE NOMINAL VARIABLES
    split1.rpart <- rpart(Surv(time, status)~x, control=ctr0.rpart)
    bcut <- as.character(temp[split1.rpart$csplit==1])
    z <- sign(is.element(x, bcut)) + 0
    score.max <- LogRank(time, status, z=z)
    if (length(bcut) >1)  bcut <- paste(bcut, collapse=",")
  } else {
    # print(control.bestcut$method0); print(temp)
    if (control.bestcut$method0 != "GreedySearch" && length(temp) > control.bestcut$threshold.ncuts) {  # USE SSS
      # STANDARDIZE x IN ORDER FOR a TO BE FIXED
      x.scale <- scale(x, center = TRUE, scale = TRUE)
      LB <- quantile(x.scale, probs = control.bestcut$d); 
      UB <- quantile(x.scale, probs =1-control.bestcut$d);
      # LB <- min(x); UB <- max(x)
      bcut0 <- optimize(Obj.LogRank, lower=LB, upper=UB, maximum=TRUE, tol=control.bestcut$epi,
                           a=control.bestcut$a0, Time=time, status=status, x=x.scale)$maximum
      bcut <- attr(x.scale,"scaled:center") + bcut0*attr(x.scale,"scaled:scale")
      z <- sign(x < bcut) + 0
      score.max <- LogRank(time, status, z=z)
    } else {
      # zcut <- temp[-1]
      zcut <- temp[-length(temp)] + diff(temp)/2
      for(j in zcut) {
        score <- NA; 
        z <- sign(x<j) + 0
        n.L <- sum(z==0); n.R <- n.R1 <- sum(z==1)
        n.L1 <- sum(z==0 & status==1); n.R1 <- sum(z==1 & status==1);
        if (min(n.L, n.R) >= 5 && min(n.L1, n.R1) >= control.bestcut$min.ndeath) {
          score <- LogRank(time, status, z)
          # score <- survdiff(Surv(time, status) ~ z, rho = 0)$chisq
        }
        if (!is.na(score) && score >= score.max) {
          coin <- TRUE
          if (score == score.max) coin <- sample(x=c(TRUE, FALSE), size=1, replace=FALSE) 
          if (coin) {score.max <- score; bcut <- j }
          } 
        }
    }
  }
  # print(score.max)
  if (!is.na(score.max) && score.max== -1e10) {bcut <- score.max <- NA}
  return(list(bcut=bcut, score.max=score.max))
}


# =========================================================================
# FUNCTION z.create() CREATES z VECTOR GIVEN SPLITTING VARIABLE AND CUTOFF
# =========================================================================
z.create <- function(vname, cutoff, dat) {
  if (is.na(vname)||is.na(cutoff)) stop("vname or cutoff is missing.")
  if (!is.element(vname, names(dat))) stop("Variable vname is not a variable in dat.")
  z <- NA
  x <- dat[, vname]
  # print(c(vname=vname, cutoff=cutoff))
  if (!is.factor(x)) {
    # print("x is not a factor")
    if (is.character(cutoff))  cutoff <- as.numeric(cutoff)
    # if (is.na(cutoff)) stop("Something is wrong with cutoff.")
    z <- sign(x< cutoff)  
  } else {
    # print("x is a factor")
    cutoff <- strsplit(cutoff, ",")[[1]]
    z <- sign(is.element(x, cutoff))  
  }
  return(z+0)
}




# ==================================================
# FUNCTION GUIDE() FINDS THE BEST SPLIT VIA GUIDE
# ==================================================
# rho=[0,1] REFERS TO DIFFERENT WEIGHTING SCHEMES IN LOGRANK TEST (rho=0)
GUIDE <- function(time, status, X, use.survdiff=TRUE, rho=0) {
  logworth.max <- -1e10 
  split.var.GUIDE <-  NA;  
  p <- NCOL(X); 
  for (j in 1:p) {
    x <- X[,j]; 
    logworth <- NULL
    if (length(unique(x)) > 1) {
      if (!is.factor(x)) {
        if (is.character(x) || length(unique(x)) <= 4) x <- as.factor(x)
        else {
          # print(x); print(quantile(x, probs = seq(0, 1, 0.25)))
          x <- cut(x, breaks=unique(quantile(x, probs = seq(0, 1, 0.25))), 
                   include.lowest =TRUE, ordered_result = FALSE)
        }
      }
      if (!anyNA(table(x)) && min(table(x), na.rm = TRUE)  >= 5) {
        if (use.survdiff) {
          fit <- survdiff(Surv(time, status)~x, rho = rho);  # LOGRANK TEST
          logworth <- -log10(fit$pvalue); 
        } else {
          fit <- coxph(Surv(time, status)~x, ties = "efron")
          # fit <- coxph(Surv(time, status)~x, ties = "breslow") 
          # fit <- coxph(Surv(time, status)~x, ties = "exact") 
          logworth <- - log10(summary(fit)$logtest[3]);    # LRT 
          # summary(fit)$sctest[3]; summary(fit)$waldtest[3]   # SCORE OR WALD TEST
        }
      }
    }
    # print(cbind(j, names(X)[j], cat=is.factor(x), value=x[1], logworth=logworth))
    if (!is.null(logworth) && !is.na(logworth) && logworth >= logworth.max) {
      coin <- TRUE
      if (logworth == logworth.max) coin <- sample(x=c(TRUE, FALSE), size=1, replace=FALSE) 
      if (coin) {
        logworth.max <- logworth; split.var.GUIDE <- j
      }
    }
  }
  if (length(split.var.GUIDE) >1) split.var.GUIDE <- split.var.GUIDE[1]
  x.split <- ifelse(is.na(split.var.GUIDE), NA, names(X)[split.var.GUIDE])
  return(x.split)
}


# ================================================
# FUNCTION partition() PROVIDES ONE SPLIT OF DATA
# ================================================
partition <- function(dat, test, name="1", cols.x, mtry=length(cols.x), 
                      max.depth=20, min.nodesize=20, min.event.childnode=3, 
                      method.split=c("default", "rpart", "GUIDE", "combined"),
                      stratify.status=TRUE, use.rpart=TRUE,       
                      stop.for.testdata=TRUE,
                      control.bestcut=ctr0.bestcut, details=FALSE) 
{
  call <- match.call(); out <- match.call(expand = F)
  out$info <- out$name.l <- out$name.r <- out$left <- out$right <- out$... <- NULL
  n <- NROW(dat); n.events <- sum(dat$status==1); 
  median.surv <- median(survfit(Surv(time, status)~ 1, data=dat))
  var <- vname <- NA; cutoff <- NA; score <- NA;    
  depth <- nchar(name) # CONTROL THE MAX TREE DEPTH
  vnames <- colnames(dat)
  if (is.null(method.split)) method.split <- "default"
  if (!is.element("id", vnames)) stop("Error in id! Please create an id variable with values 1:n.")
  if (!is.null(test) && depth < max.depth && 
      n >= min.nodesize && n.events >= 2*min.event.childnode) {
    # CHECK AND APPLY mtry
    if (!is.null(mtry)) {
      if (!is.integer(mtry) | mtry <1 | mtry > length(cols.x)) 
        stop("Error in mtry! It must be an integer between 1 and p=length(cols.x)") 
      else cols.x <- sample(cols.x, size=mtry, replace=F)
    }
    form0 <- as.formula(paste("Surv(time, status) ~ ", paste(vnames[cols.x], collapse =" + "), sep=""))
    # DIFFERENT WAYS OF SPLITTING
    if (method.split == "rpart") {
      split.rpart <- rpart.split1rule(rpart(formula=form0, data=dat, y=FALSE, control=ctr0.rpart))
      vname <- split.rpart$xsplit; # SPLITTING VARIABLE 
      # cutoff <- split.rpart$cutpoint
    } else if (method.split == "GUIDE") {
      vname <- GUIDE(time=dat$time, status=dat$status, X=dat[,cols.x], use.survdiff=TRUE, rho=0)
    } else if (method.split == "combined") {
      if (stratify.status)  {
        dat.b <- dat %>% group_by(status) %>%
			# sample_frac(2/3, replace = FALSE) %>%   # TEST SAMPLE
			sample_frac(1, replace = TRUE) %>% # BOOTSTRAP
			# slice_sample(n=trunc(n/2), replace =TRUE) %>%        # BOOTSTRAP    # NOT RIGHT
			# slice_sample(n=trunc(n/3), replace =FALSE) %>%         # TEST SAMPLE
			as.data.frame() # %>% mutate(time=jitter(time))  # JITTERING TIME TO REDUCE TIES
        dat.oob <- dat[-dat.b$id, ]
      } else  {
        id.b <- sample(x=1:n, size=n, replace=TRUE)          # BOOTSTRAP
        # id.b <- sample(x=1:n, size=n/3, replace=FALSE)   # TEST SAMPLE
        dat.b <- dat[id.b, ]; dat.oob <- dat[-id.b, ]
      }
      # GUIDE
      xsplit.b.GUIDE <- GUIDE(time=dat.b$time, status=dat.b$status, X=dat.b[,cols.x])
      bcut.b.GUIDE <- NA
      # print(xsplit.b.GUIDE)
      if (!is.na(xsplit.b.GUIDE)) {
        tmp <- best.cut(time=dat.b$time, status=dat.b$status, x=unlist(dat.b[xsplit.b.GUIDE]), 
                      control.bestcut=control.bestcut)
        bcut.b.GUIDE <- tmp$bcut
      } 
      if (is.na(bcut.b.GUIDE)) xsplit.b.GUIDE <- NA
      # print(cbind(xsplit.b.GUIDE, bcut.b.GUIDE))
      
      # GREEDY SEARCH
      if (use.rpart) {
        split.b.GS <- rpart.split1rule(rpart(formula=form0, data=dat.b, y=FALSE, control=ctr0.rpart))
        xsplit.b.GS <- split.b.GS$xsplit; # SPLITTING VARIABLE 
        bcut.b.GS <- split.b.GS$cutpoint
      } else {
        xsplit.b.GS <- bcut.b.GS <- NA
        max.score.b <- -1e10
        if (!is.na(xsplit.b.GUIDE) && !is.na(tmp$score.max)) {
          xsplit.b.GS <- xsplit.b.GUIDE; 
          bcut.b.GS <- bcut.b.GUIDE
          max.score.b <- tmp$score.max
        }
        for (j in cols.x) {
          if (!is.na(xsplit.b.GUIDE) && vnames[j] != xsplit.b.GUIDE) {
            split.j <- best.cut(dat.b$time, dat.b$status, x=dat.b[,j], 
                                control.bestcut=control.bestcut)
            score.j <- split.j$score.max
            # print(cbind(score.j, max.score.b))
            if (!is.na(score.j) && score.j >= max.score.b) {
              max.score.b <- score.j; xsplit.b.GS <- vnames[j]; bcut.b.GS <- split.j$bcut
            }
          }
        }
      }
      if (is.na(bcut.b.GS)) xsplit.b.GS <- NA
      
      # COMPARE GS & GUIDE VIA CROSS VALIDATION
      if (is.na(xsplit.b.GUIDE) && is.na(xsplit.b.GS)) vname <- NA
      else if (is.na(xsplit.b.GUIDE) && !is.na(xsplit.b.GS)) vname <- xsplit.b.GS
      else if (!is.na(xsplit.b.GUIDE) && is.na(xsplit.b.GS)) vname <- xsplit.b.GUIDE
      else { 
        if (xsplit.b.GS==xsplit.b.GUIDE) {
          vname <- xsplit.b.GS 
          cutoff <- bcut.b.GS
        } else {
          z.GS <- z.create(vname=xsplit.b.GS, cutoff=bcut.b.GS, dat=dat.oob)
          score.GS <- LogRank(dat.oob$time, dat.oob$status, z=z.GS, min.event.childnode=min.event.childnode)
          z.GUIDE <- z.create(vname=xsplit.b.GUIDE, cutoff=bcut.b.GUIDE, dat=dat.oob)
          score.GUIDE <- LogRank(dat.oob$time, dat.oob$status, z=z.GUIDE, min.event.childnode=min.event.childnode)
          if (is.na(score.GUIDE) && is.na(score.GS)) vname <- sample(c(xsplit.b.GS, xsplit.b.GUIDE), size=1)   # RANDOMLY SELECT ONE
          else if (is.na(score.GUIDE) && !is.na(score.GS)) vname <- xsplit.b.GS
          else if (!is.na(score.GUIDE) && is.na(score.GS)) vname <- xsplit.b.GUIDE 
          else vname <- ifelse(score.GS > score.GUIDE, xsplit.b.GS, xsplit.b.GUIDE)
        }
      }
    } else {  # THE DEFAULT - INTERSECTED VALIDATION (IV)
      ########  PREPARE DATA FOR INTERTSECTED VALIDATION W/ OVERLAPPED BAGS
      #### THE TRAINING SET
      dat <- dat %>% mutate(group = sample(1:3, size=n, replace=TRUE))
      id.grp2 <- dat$id[dat$group==2]
      dat.b.tmp <- dat %>% filter(group!=3) %>% 
        # group_by(status) %>%                          # STRATIFIED BY CENSORING STATUS
        slice_sample(prop=1.0, replace =TRUE, by=status) 
      dat.b <- rbind(dat %>% filter(group==1),  dat.b.tmp) %>%        
        as.data.frame() # %>% mutate(time=jitter(time))  
      #### THE VALIDATION SET
      id2.oob <- id.grp2[which(is.element(id.grp2, dat.b.tmp$id))]
      dat.grp2.oob <- dat %>% filter(is.element(id, id2.oob))
      dat.ind <- rbind(dat %>% filter(group==3),  dat.grp2.oob)
      n.remaining <- n - NROW(dat.ind)
      dat.oob.tmp <- dat %>% filter(group!=1) %>% 
        # group_by(status) %>% 
        slice_sample(prop=n.remaining/n*3/2, replace =TRUE, by=status)       
      dat.oob <- rbind(dat.ind, dat.oob.tmp) %>% as.data.frame() # %>% mutate(time=jitter(time))  # JITTERING TIME TO REDUCE TIES      
      dat <- dat %>% dplyr::select(-group)
      # print(cbind(NROW(dat.b), NROW(dat.oob)))
      
      score.max <- -1e10
      xsplit.BV <-  NA; 
      # x.split <- xcut <- NA
      min.ndeath <- 3
      for (j in 1:length(cols.x)) {
        x <- dat.b[, cols.x[j]]
        xname <- names(dat)[cols.x[j]]
        if (is.character(x)) x <- as.factor(x)
        if (is.factor(x) && length(levels(x)) > 2 ) {
          split.rpart <- rpart(Surv(time, status)~x, data=dat.b, control=ctr0.rpart)
          # split.var <- split.rpart$frame[1, 1]
          # split.bcut <- split.rpart$splits[4]; split.bcut
          bcut <- sort(unique(x))[split.rpart$csplit==1]
          z.test <- sign(is.element(dat.oob[, cols.x[j]], bcut)) + 0
        } else {
          bcut <- best.cut(dat.b$time, dat.b$status, x)$bcut
          z.test <- sign(dat.oob[, cols.x[j]] < bcut) + 0
        }
        # print(c(xname, bcut))
        if (!anyNA(z.test)) {
          score <- LogRank(Time=dat.oob$time, status=dat.oob$status, z=z.test)
          # print(cbind(j, xname, bcut, score))
          # score <- survdiff(Surv(time, status) ~ z.test, data=dat.oob, rho = 0)$chisq
          if (!is.na(score) && score >= score.max) {
            coin <- TRUE
            if (score == score.max) coin <- sample(x=c(TRUE, FALSE), size=1, replace=FALSE) 
            if (coin) {
              score.max <- score; xsplit.BV <- xname
              # x.split <-  cols.x[j]   # COLUMN NUMBER
              # xcut <- bcut            # BEST CUTOFF POINT
              }
            }
        }
      }
      vname <- xsplit.BV
    }
  }
  # print(vname); 
  
  # OBTAIN INFO ON THE BEST SPLITTING VARIABL 
  if (!(is.na(vname))) {
    var <- which(vnames==vname)     # COLUMN NUMBER OF BEST SPLITTING VARIABLE
    split.tmp <- best.cut(dat$time, dat$status, x=dat[,vname], 
                          control.bestcut=control.bestcut)
    cutoff <- split.tmp$bcut
    score <- split.tmp$score.max
    if (is.na(score)) var <- vname <- cutoff <- NA
    else {
      out$name.l <- paste(name, 0, sep=""); out$name.r <- paste(name, 1, sep="")
      z.dat <- z.create(vname=vname, cutoff=cutoff, dat=dat)
      out$left  <- dat[z.dat==1,]
      out$right <- dat[z.dat==0, ]  
    }
  }
  
  # HANDLE THE TEST SAMPLE PART
  ntest <- ntest.events <- median.surv.test <- score.test <- NA
  out$left.test <- out$right.test <- NULL
  if (!is.null(test)) {
    if (base::identical(dat, test)) {
      ntest <- n; ntest.events <- n.events; median.surv.test <- median.surv
      if (!is.na(score)) {score.test <- score; out$left.test <-  out$left; out$right.test <- out$right} 
    } else {
      ntest <- NROW(test); ntest.events <- sum(test$status==1); 
      if (!is.na(score) & ntest > 3 & ntest.events > 1)  {
        median.surv.test <- median(survfit(Surv(time, status)~ 1, data=test))
        z.test <- z.create(vname=vname, cutoff=cutoff, dat=test)
        score.test <- LogRank(test$time, test$status, z=z.test, min.event.childnode=min.event.childnode) 
        if (!is.na(score.test)) {
          out$left.test <- test[z.test==1,  ]
          out$right.test <- test[z.test==0,  ]
        } else if(stop.for.testdata) {
          out$left  <- out$right <- NULL;
          vname <- var <- score <- cutoff <- NA;
          } 
      }
    } 
  }
  # OUTPUT
  out$info <- data.frame(node=name, size = n, n.events=n.events, median.surv=median.surv, var = var, vname=vname, 
                         cutoff= cutoff, score=score, size.test=ntest, 
                         ntest.events=ntest.events, median.surv.test=median.surv.test, score.test=score.test)
  return(out) 
}


# ===================================
# FUNCTION grow() GROWS A LARGE TREE
# ===================================

grow <- function(data, test, cols.x, mtry=length(cols.x), 
                 max.depth=20, min.nodesize=20, min.event.childnode=5, 
                 method.split="default-combined", stratify.status=TRUE, use.rpart=TRUE,
                 stop.for.testdata=TRUE, control.bestcut=ctr0.bestcut)
{
  out <- list.nd <- list.test <- temp.list <- temp.test <- temp.name <- NULL
  list.nd <- list(data); list.test <- list(test)
  name <- 1
  while (length(list.nd)!=0) {    
    for (i in 1:length(list.nd)){
      if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 1){ 
        split <- partition(dat=list.nd[[i]], test=list.test[[i]], name=name[i],
                          cols.x=cols.x, mtry=mtry, max.depth=max.depth, 
                          min.nodesize=min.nodesize, min.event.childnode=min.event.childnode, 
                          method.split=method.split, stratify.status=stratify.status, use.rpart=use.rpart,       # COMBING GS & GUIDE 
                          stop.for.testdata=stop.for.testdata, control.bestcut=control.bestcut) 
        # print(split$info)
        out <- rbind(out, split$info)
        if (stop.for.testdata) cond0 <- (!is.null(split$left) && !is.null(split$left.test))
        else cond0 <- (!is.null(split$left))
        if (cond0){
          temp.list <- c(temp.list, list(split$left, split$right))
          temp.test <- c(temp.test, list(split$left.test, split$right.test))
          temp.name <- c(temp.name, split$name.l, split$name.r)
        }}}
    list.nd <- temp.list; list.test <- temp.test; name <- temp.name
    temp.list <- temp.test <- temp.name <- NULL
  }   
  out$node <- as.character(out$node)
  tre <- out[order(out$node), ] %>% as.data.frame()
  
  # QUICK CLEARING 
  # -----------------
  names.2NA <- c("var", "vname", "cutoff", "score", "score.test")
  nodes.problematic <- nodes.problematic <- tre %>% 
    mutate(size=as.numeric(size), n.events=as.numeric(size)) %>% 
    filter(size < 5 & n.events <= min.event.childnode) %>% 
    select(node) %>% unlist() %>% unname()
  if (length(nodes.problematic) >0) {
    # print(tre)
    for (h in nodes.problematic) {
      nod.rm <- substring(h, first=1, last=nchar(h)-1)
      if (is.element(nod.rm, tre$node)) {
        tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
        tre[match(nod.rm, tre$node), names.2NA] <- NA
      }
    }
    # print(tre)
  }
  return(tre) 
}



# ===================================================================
# PLOTTING IT TREE STRUCTURE, MODIFIED FROM PETER CALHOUN'S CODES
# ===================================================================

plot.tree <- function(tree, cols.nominal=NULL, color0=NULL, 
                      leaf.label=FALSE, pch.leaf=21, plot.score=TRUE, 
                      textDepth=10, lines="rectangle", digits=2)
{
  depth<-max(nchar(tree[,1]))
  par(xaxs='i')
  par(mar=c(1,1,1,1))
  par(xpd=TRUE)
  plot(1, type="n", xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), axes=FALSE,xaxs="i",yaxs="i")
  nodes <- tree$node
  # nodesBin <- gsub("1", "0", nodes)
  # nodesBin <- gsub("2", "1", nodesBin)
  nodesBin <- nodes
  lastObs <- nchar(nodesBin)
  nodesBin <- substr(nodesBin,2,lastObs)
  var <- tree$vname 
  col.var <- tree$var   
  cut <- as.character(tree$cutoff) 
  size <- tree$size  
  m <- tree$median.surv
  if (plot.score) score <- round(as.numeric(as.character(tree$score)), digits = digits)
  if (is.null(color0)) {
    if (!is.null(tree$color))  color0 <- tree$color
    else color0 <- rep("gray60", length(nodesBin))
  } else if (length(color0) != NROW(tree)) 
    stop("color0 must have the same length as nrow(tree).")
  
  for(i in 1:length(nodesBin)){
    nChar <- nchar(nodesBin[i])
    # print(cbind(i, nChar))
    if(!is.na(var[i])){
      j <- ifelse(i==1, i+1, i)
      if(lines=="rectangle"){
        # print(cbind(i, nodesBin[i], strtoi(nodesBin[i])))
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+1),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+1),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+1)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+3),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
      } else if(lines=="triangle"){
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+2),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+1)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+2),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
      }         
      if(nChar <= textDepth){ 
        # print(cbind(i=i, var[i], cut[i]))
        if (plot.score) {
          score.i <- as.character(score[i])
          text(x=(1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+2), 
                                       y=(depth-nChar)/(depth+1)+1/(depth+20)-0.03, 
                             labels=substitute(paste(bold(score.i))), cex=1, col="red")
        }
        # HANDLE SPLITS WITH NOMINAL OR ORDINAL X
        if (is.na(as.numeric(as.character(cut[i]))) | 
            (!is.null(cols.nominal) && is.element(col.var[i], cols.nominal))) {
          cutpoint <- unlist(strsplit(as.character(cut[i]), split=" "))
          cutpoint <- paste("{", paste(cutpoint, collapse=" "), "}", sep="")
          text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+2), 
               (depth-nChar)/(depth+1)+1/(depth+20), 
               bquote(.(as.character(var[i]))%in%.(cutpoint)),cex=1, col="blue")
        } else {
          cutpoint <- round(as.numeric(cut[i]), digits=digits)
          text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[j], base = 2)+2), 
               (depth-nChar)/(depth+1)+1/(depth+20), 
               bquote(.(as.character(var[i]))<.(cutpoint)),cex=1, col="blue")
        }
      }
    } else {
      if(nChar <= textDepth){
        median.i <- round(as.numeric(m[i]), digits=digits)
        points((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)+0.01, 
               pch=pch.leaf, col=color0[i], bg=color0[i], cex=3)
        # points((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)+0.01, 
        #       pch=21, col="gray65", bg="gray65", cex=3)
        
        if (leaf.label) {
          text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)-0.02, 
               nodes[i],cex=1, offset=1, col=color0[i])
        } else {
          text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)-0.02, 
             paste("n=", size[i], sep=""),cex=1, offset=1, col=color0[i])
          text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2), (depth-nChar)/(depth+1)-0.05, 
             paste("m=", median.i, sep=""), cex=1, offset=1, col=color0[i])
        }
      }
    }
  }
}



# ==========================================================
# FUNCTION de() FINDS ALL THE DESCENDENTS OF NODE x IN tree
# DOWN UP TO depth
# ==========================================================
de <- function(x, tree, depth=NULL)
{
  if(length(x) != 1) {
	# print(x)
	stop("The length of x in function de must be 1.") 
	}
  if (is.null(depth)) depth <- max(nchar(tree$node)) + 10
  y <- tree$node;  de <- NA
  if(sum(match(x, y), na.rm = T) != 0) {
    temp <- 1:length(y)
    start <- match(x, y) + 1    
    end <- length(y)
    if(start <= length(y) & nchar(y[start]) > nchar(x)) {
      temp1 <- temp[temp >= start & nchar(y) <= nchar(x)][1] - 1
      if(!is.na(temp1)) end <- temp1
      de <- y[start:end]
    }}
  de <- de[(nchar(de)-nchar(x)) <= depth]
  de
}



# =================================================
# FUNCTION senddown() SENDS DATA DOWN TO A TREE
# =================================================

senddown <- function(data, tre){
  if (is.null(tre)) stop("tre is NULL. Cann't send dat down.")
  dat <- cbind(data, node=1)
  if (NROW(tre) > 1) {
    for (i in 1:NROW(tre)){
      if (!is.na(tre$var[i])){
        # print(cbind(i, tre$node[i], tre$vname[i], tre$cutoff[i]))
        x.split <-  dat[, tre$vname[i]];
        cutoff <- tre$cutoff[i]
        in.node <- (dat$node)==(tre$node[i]); 
        # print(sum(in.node))
        if (is.factor(x.split)) {
          cutoff <- strsplit(cutoff, ",")[[1]]
          l.nd <- dat$node[in.node & is.element(x.split, cutoff)] 
          dat$node[in.node & is.element(x.split, cutoff)] <- paste(l.nd, 0, sep="")  
          r.nd <- dat$node[in.node & !is.element(x.split, cutoff)] 
          dat$node[in.node & !is.element(x.split, cutoff)] <- paste(r.nd, 1, sep="") 
        } else {
          cutoff <- as.numeric(as.character(cutoff))
          l.nd <- dat$node[in.node & x.split < cutoff] 
          dat$node[in.node & x.split < cutoff] <- paste(l.nd, 0, sep="")  
          r.nd <- dat$node[in.node & x.split >= cutoff] 
          dat$node[in.node & x.split >=  cutoff] <- paste(r.nd, 1, sep="") 
        }
      }
    }
  }
  dat$node <- as.factor(dat$node) 
  if (!is.null(tre$grp)) {
    dat$grp <- as.numeric(Replace(as.character(dat$node), tre$node, tre$grp))  
  }
  return(dat) 
}



# ====================================================================
# Pruning and Size Selection Based on LeBlanc and Crowley (JASA, 1992)
# ====================================================================
pruning <- function(tre, ntmnl.min=1, print.it=FALSE)
{
  if(is.null(dim(tre))) stop("No Need to Prune Further.")
  n.tmnl <- sum(is.na(tre$vname)); 
  result <- NULL; 
  if (n.tmnl > 1) {
    subtree <- 1            
    names.2NA <- c("var", "vname", "cutoff", "score", "score.test")
    while (n.tmnl > max(ntmnl.min, 2)) {
      # print(tre)
      internal <- tre$node[!is.na(tre$vname)]; n.links <- length(internal); 
      r.value <- 1:n.links
      for(i in 1:n.links) {
        branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
        score <- as.numeric(as.vector(branch$score))
        r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
      }
      alpha <- min(r.value)
      nod.rm <- internal[r.value == alpha]; 
      if (length(nod.rm)>1) print("Multiple Nodes will be pruned. Check!")
      G <- sum(as.numeric(as.vector(tre$score.test)), na.rm=T); 
      result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tre), 
                                    size.tmnl=nrow(tre)-n.links, alpha=alpha, G=G))
      tre0 <- tre
      tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
      tre[match(nod.rm, tre$node), names.2NA] <- NA
      n.tmnl <- sum(is.na(tre$vname))
      if (n.tmnl < ntmnl.min)  tre <- tre0
      else subtree <- subtree + 1          
    }
    # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
    if (ntmnl.min==1) result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre), 
                                                    size.tmnl=1, alpha=9999, G=0))     
    result <- as.data.frame(result)
    if (print.it) print(result)
  }
  return(tre)
}
# tre.pruned <- pruning(tre0, ntmnl.min=10, print.it=TRUE)
# tre.pruned



# ==========================================================================
# FUNCTION prune.validate() IMPLEMTNS SPLIT-COMPLEXITY PRUNING PLUS VALIDATION
# ==========================================================================

prune.validate <- function(dat, test, tre)
{
  # SPLIT-COMPLEXITY PRUNING AND VALIDATION
  #  tre0.dev <- tre0.bic <- tre
  y.test <- Surv(test$time, test$status)
  n.tmnl <- sum(is.na(tre$vname)); 
  result <- NULL; 
  # min.bic <- min.dev <- 1e10;  
  subtree <- 1            
  names.2NA <- c("var", "vname", "cutoff", "score", "score.test")
  while (n.tmnl >= 1) {
    # COMPUTE (PREDICTED) DEVIANCE 
    dat.tmp <- senddown(data=dat, tre=tre) %>% dplyr::select(time, status, node)
    if (n.tmnl ==1) fit.tmp <- coxph(Surv(time, status)~ 1, data=dat.tmp)
    else fit.tmp <- coxph(Surv(time, status)~ node, data=dat.tmp)
    pred.tmp <- predict(fit.tmp, newdata=senddown(data=test, tre=tre), type="lp")
    dev.tmp <- glmnet::coxnet.deviance(pred = pred.tmp, y = y.test)
    #     bic.tmp <- dev.tmp + log(sum(test$status))*n.tmnl
    #    if (dev.tmp <= min.dev) {tre0.dev <- tre; min.dev <- dev.tmp}
    #    if (bic.tmp <= min.bic) {tre0.bic <- tre; min.bic <- bic.tmp}
    if (n.tmnl >=2) {
      # print(tre)
      internal <- tre$node[!is.na(tre$vname)]; n.links <- length(internal); 
      r.value <- 1:n.links
      for(h in 1:n.links) {
        branch <- tre[is.element(tre$node, c(internal[h], de(internal[h], tree=tre))),]
        score <- as.numeric(as.vector(branch$score))
        r.value[h] <- sum(score, na.rm=T) / sum(!is.na(score))
      }
      # print("-------------------------------------")
      # print(cbind(internal, r.value))
      alpha <- min(r.value)
      nod.rm <- internal[r.value == alpha]; 
      if (length(nod.rm)>1) {
        print("Multiple Nodes will be pruned. Check!")
        nod.rm <- nod.rm[1]
      }
      G <- sum(as.numeric(as.vector(tre$score.test)), na.rm=T); 
      result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tmnl=n.tmnl, 
                                    alpha=alpha, G=G, dev=dev.tmp))
      tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
      tre[match(nod.rm, tre$node), names.2NA] <- NA
      n.tmnl <- sum(is.na(tre$vname))
    } else {
      result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tmnl=1, 
                                    alpha=9999, G=0, dev=dev.tmp)) 
      n.tmnl <- 0
    } 
    subtree <- subtree + 1         
  }
  result <- result %>% as.data.frame %>%  
    mutate(G=as.numeric(as.character(G)), size.tmnl=as.numeric(as.character(size.tmnl)),
           G.2=G-2*size.tmnl, dev=as.numeric(as.character(dev)), 
           BIC=dev+log(NROW(test))*size.tmnl) 
  return(result)
}





# =====================================================
# FUNCTION obtain.tree() REMOVES SOME NODES FROM TREE
# =====================================================
obtain.tree <- function(tre, nodes.rm=NULL) {
  if (!is.null(nodes.rm) && !(sum(is.na(nodes.rm)) > 0)) {
    names.2NA <- c("var", "vname", "cutoff", "score", "score.test")
    for (h in nodes.rm) {
      desc.h <- de(h, tre)
      tre <- tre %>% filter(!is.element(node, desc.h))
      tre[match(h, tre$node), names.2NA] <- NA
    }
  }
  return(tre)
}




# ========================================================================
# FUNCTION clean.up() TRUNCATES NOISY BRANCHES OF INITIAL TREE 
# =====================================================================
clean.up <- function(tre, print.it=FALSE, sig.level=0.05, 
                     use.score.test=TRUE, 
                     remove.missing.median.survival=FALSE) {
  cols.change <- c("var", "vname", "cutoff", "score", "score.test")
  # STEP I:  REMOVE NODES WITH MISSING MEDIAN SURVIVAL
  Nodes0 <- NULL
  if (remove.missing.median.survival) {
    # if (use.score.test) 
      Nodes0 <- tre$node[is.na(tre$median.surv) | is.na(tre$median.surv.test)]
    # else Nodes0 <- tre$node[is.na(tre$median.surv)]
  }
  if (use.score.test) Nodes0 <- c(Nodes0, tre$node[is.na(tre$score.test) & !is.na(tre$score)])
  # print(Nodes)
  if (!is.null(Nodes0)) {
    Nodes0 <- Nodes0[order(nchar(Nodes0))]
    n.Nodes0 <- length(Nodes0)
    if ( n.Nodes0 > 0) {
        for (i in 1:n.Nodes0) {
          node.i <- Nodes0[i]    
          if (is.element(node.i, tre$node)) {
            mother.i <- substr(node.i, start=1, stop=nchar(node.i)-1)
            nodes.rm <- de(mother.i, tre)
            tre <- tre[!is.element(tre$node, nodes.rm), ]
            tre[tre$node==mother.i, cols.change] <- NA
          }
        }
    }
  }
  
  # STEP II: REMOVE INTERNAL NODES WITH NO SIGNIFICANT SPLIT IN DESCEDANTS
  nodes.internal <- tre$node[!is.na(tre$score)] 
  nodes.internal <- nodes.internal[order(nchar(nodes.internal))]
  out <- NULL
  critival.value <- qchisq(p=sig.level, df=2, lower.tail = FALSE)  # 3.84 or 6
  critival.value.test <- qchisq(p=sig.level*2, df=1, lower.tail = FALSE) # MORE LIBERAL
  for (h in nodes.internal) {
    score <- as.numeric(as.character(tre$score)) 
    score.test <- as.numeric(as.character(tre$score.test)) 
    # print(h)
    score.h <- tre$score[tre$node==h]
    score.test.h <- tre$score.test[tre$node==h]
    if (use.score.test) cond1 <- (!is.na(score.h) && !is.na(score.test.h))
    else cond1 <- (!is.na(score.h))
    # print(cbind(h, score.h, score.test.h, cond1))
    if (is.element(h, tre$node) & cond1) { 
      des.h <- c(h, de(h, tre, depth=NULL))
      max.score.h <- max(score[is.element(tre$node, des.h)], na.rm =TRUE)
      # print(cbind(h=h, max.score=max.score.h))
      if (use.score.test) {
        max.score.test.h <- max(score.test[is.element(tre$node, des.h)], na.rm =TRUE)
        cond2 <- (max.score.h < critival.value | max.score.test.h < critival.value.test)
      } else cond2 <- (max.score.h < critival.value)
      if (cond2) {
        out <- rbind(out, c(node.cut=h, max.score=max.score.h))
        tre <- tre[!is.element(tre$node, de(h, tre)), ]
        tre[tre$node==h, cols.change] <- NA
      }
    }
  }
  
  out <- out |> as.data.frame()
  if (print.it) {print(out); print(tre)} 
  return(tre)
}








# =================================================
# FUNCTION hazard.sorting() SORTS NODES WRT HAZARD
# =================================================
# IT ALSO PREPARES INPUTS FOR SUBSEQUENT glmnet FITTING
hazard.sorting <- function(dat, tre=NULL, according=NULL) {
  if (!is.element("node", names(dat))) stop("dat doesnot have the node variable.")
  if (!is.factor(dat$node)) dat$node <- as.factor(dat$node)
  leaves <- levels(dat$node)
  K <- length(leaves)  
  if (!is.null(according) && according=="median.survival") {
    if (is.null(tre)) stop("With this option, you need a tre to extract median survival times.")
    leaves.tre <- tre$node[is.na(tre$vname)]
    # print(sort(leaves)); print(sort(leaves.tre))
    if (!identical(sort(leaves), sort(leaves.tre))) stop("Leaves in tre and dat don't match.")
    median.surv <- as.numeric(as.character(tre$median.surv[is.na(tre$vname)]))
    measure.sorted <- sort(median.surv)
    leaves.sorted <- leaves.tre[order(median.surv)]
    } else {
      #node.base0 <- leaves[1]
      fit0 <- coxph(Surv(time, status)~node, data=dat) 
      Beta0 <- c(0, fit0$coefficients); names(Beta0) <- leaves
      if (sum(is.na(Beta0))) {
        print("This is a bit weird, but you got a NA coeffecticy estimate in coxph. Let's check:")
        print(Beta0); print(table(dat$node)); 
        print(tre); plot.tree(tre) 
      }
      measure.sorted <- sort(Beta0, na.last=FALSE)   # THE na.ast=TRUE OPTION FOR NA BETA ESTIMATES 
      leaves.sorted <- names(measure.sorted)
    }
  dat <- dat %>% 
    mutate(node=factor(node, ordered = TRUE, levels=leaves.sorted), 
           node.rank=as.numeric(node)) # NO NEED OF RECODING
  contrasts(dat$node) <- contr.treatment(n=K, base = 1)
  fit1 <- coxph(Surv(time, status)~node, data=dat, x=TRUE) 
  # print(coef(fit1))
  key <- data.frame(node=leaves.sorted, measure.sorted=measure.sorted, rank=1:K, 
                    beta.coxph=c(0, fit1$coefficients)) 
  X1 <- as.matrix(fit1$x)
  # attr(X1, "contrasts")$node
  B.inverse <- matrix(1, nrow=K-1, ncol=K-1) 
  B.inverse[upper.tri(B.inverse)] <- 0
  X <- X1%*%B.inverse
  colnames(X) <- leaves.sorted[-1]  # NOT EXACTLY RIGHT, BUT USEFUL LATER
  y <- Surv(time=dat$time, event=dat$status)
  # penalty.factor.glmnet <- c(0, rep(1, K-2))   
  # penalty.factor.glmnet <- rep(1, K-1) # TO ALLOW FOR 0 ESTIAMTE OF beta1
  penalty.factor.glmnet <- 1/abs(diff(key$beta.coxph))   # ADAPTIVE LASSO
  return(list(dat=dat, K=K, leaves.sorted=leaves.sorted, key=key, 
              B.inverse=B.inverse, X=X, y=y, 
              penalty.factor.glmnet=penalty.factor.glmnet))
}





# ==================================================
# FUNCTION shear() TRUNCATES TREES VIA MERGED GROUP
# ==================================================

shear <- function(tre, group.name=NULL) {
  nodes.internal <- tre$node[!is.na(tre$vname)] 
  nodes.internal <- nodes.internal[order(nchar(nodes.internal))]
  cols.change <- c("var", "vname", "cutoff", "score", "score.test")
  if (is.null(group.name)) group.name <- "grp.BIC"
  for (h in nodes.internal) {
    grp <- tre[, group.name]
    grp.h <- grp[is.element(tre$node, de(h, tre))] |> na.omit() |> unique()
    if (length(grp.h) == 1) {
      tre <- tre[!is.element(tre$node, de(h, tre)), ]
      tre[tre$node==h, cols.change] <- NA
      tre[tre$node==h, group.name] <- grp.h 
    }
  }
  grp <- tre[, group.name] %>% as.numeric()
  n.grp <- length(unique(grp[!is.na(grp)])); # print(n.grp) 
  palet <- viridis_pal(alpha = 1, begin = 0.2, end = 1, option="D")(n.grp)
  # print(tre)
  tre$color <- palet[grp] %>% replace_na("gray50")
  # tre$color <- tre[, group.name]
  return(tre)
}




# ====================================================================================
# FUNCTION SurvTreeFuL.TestSample() CONSTRUCTS SURVIVAL TREEFUL MODELS VIA TEST SAMPLE
# ====================================================================================

RowSD <- function(x) {sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))}

SurvTreeFuL.TestSample <- function(dat, test, cols.x, mtry=length(cols.x), 
                                    plot.it=FALSE, details=FALSE, relaxed=FALSE, 
                                    max.depth=20, min.nodesize=20, min.event.childnode=min.nodesize/4, 
                                    remove.missing.median.survival=FALSE,
                                    method.split="default-combined", stratify.status=TRUE, 
                                    use.rpart=TRUE, starting.tree=NULL, 
                                    control.bestcut=ctr0.bestcut) 
{
  name.metrics <- c("dev", "BIC")
  if(is.null(test)) stop("In the test sample method, you need a test data set.")
  
  # FIRST OBTAIN CART TREES
  # ***************************************************
  form0 <- base::paste("Surv(time, status) ~ ",  base::paste(names(dat)[cols.x], collapse =" + ")) |> 
    as.formula()
  fit.cart <- cart.test(formula=form0, data=dat, test=test, n0=min.nodesize, plot.it=FALSE)
  size.cart <- fit.cart$tree.size 
  cart.initial <- fit.cart$tre0
  cart.bic <- fit.cart$BTREES$btree.BIC
  cart.dev <- fit.cart$BTREES$btree.dev
  if (plot.it) {plot.cart(cart.bic); title(main="CART via BIC")}
  
  # GROW THE LARGE INITIAL TREE WITH ENHANCEMENTS
  # ***************************************************
  tre0.initial <- grow(data=dat, test=test, cols.x=cols.x, mtry=mtry, 
                       max.depth=max.depth, min.nodesize=min.nodesize, min.event.childnode=min.event.childnode, 
                       method.split=method.split, stratify.status=stratify.status, use.rpart=use.rpart,
                       control.bestcut=control.bestcut)
  if (plot.it) plot.tree(tre0.initial, cols.nominal=NULL, textDepth=10, leaf.label=TRUE,
                         plot.score=TRUE, lines="rectangle", digits=3)
  tre0.clean <- clean.up(tre0.initial, print.it=FALSE, sig.level=0.05, use.score.test=TRUE, 
                         remove.missing.median.survival=remove.missing.median.survival)
  if (plot.it) plot.tree(tre0.clean, cols.nominal=NULL, textDepth=10, leaf.label=TRUE,
                         plot.score=TRUE, lines="rectangle", digits=3)
  
  # COMBINE LEBLANC AND CROWLEY (2021) WITH CV-ED DEVIANCE - TEST SAMPLE
  # ---------------------------------------------------------------------
  tre <- tre0.clean
  tre0.dev <- tre0.bic <- tre
  y.test <- Surv(test$time, test$status)
  n.tmnl <- sum(is.na(tre$vname)); 
  result <- NULL; min.bic <- min.dev <- 1e10;  subtree <- 1            
  names.2NA <- c("var", "vname", "cutoff", "score", "score.test")
  while (n.tmnl >= 1) {
    # COMPUTE (PREDICTED) DEVIANCE 
    dat.tmp <- senddown(data=dat, tre=tre) %>% dplyr::select(time, status, node)
    if (n.tmnl ==1) fit.tmp <- coxph(Surv(time, status)~ 1, data=dat.tmp)
    else fit.tmp <- coxph(Surv(time, status)~ node, data=dat.tmp)
    pred.tmp <- predict(fit.tmp, newdata=senddown(data=test, tre=tre), type="lp")
    dev.tmp <- glmnet::coxnet.deviance(pred = pred.tmp, y = y.test)
    bic.tmp <- dev.tmp + log(sum(test$status))*n.tmnl
    if (dev.tmp <= min.dev) {tre0.dev <- tre; min.dev <- dev.tmp}
    if (bic.tmp <= min.bic) {tre0.bic <- tre; min.bic <- bic.tmp}
    if (n.tmnl >=2) {
      # print(tre)
      internal <- tre$node[!is.na(tre$vname)]; n.links <- length(internal); 
      r.value <- 1:n.links
      for(h in 1:n.links) {
        branch <- tre[is.element(tre$node,c(internal[h], de(internal[h], tree=tre))),]
        score <- as.numeric(as.vector(branch$score))
        r.value[h] <- sum(score, na.rm=T) / sum(!is.na(score))
      }
      # print("-------------------------------------")
      # print(cbind(internal, r.value))
      alpha <- min(r.value)
      nod.rm <- internal[r.value == alpha]; 
      if (length(nod.rm)>1) {
        print("Multiple Nodes will be pruned. Check!")
        nod.rm <- nod.rm[1]
      }
      G <- sum(as.numeric(as.vector(tre$score.test)), na.rm=T); 
      result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tmnl=n.tmnl, 
                                    alpha=alpha, G=G, dev=dev.tmp))
      tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
      tre[match(nod.rm, tre$node), names.2NA] <- NA
      n.tmnl <- sum(is.na(tre$vname))
    } else {
      result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tmnl=1, 
                                    alpha=9999, G=0, dev=dev.tmp)) 
      n.tmnl <- 0
    } 
    subtree <- subtree + 1         
  }
  result <- result %>% as.data.frame %>%  
    mutate(G=as.numeric(as.character(G)), size.tmnl=as.numeric(as.character(size.tmnl)),
           G.2=G+2*size.tmnl, dev=as.numeric(as.character(dev)), 
           BIC=dev+log(NROW(test))*size.tmnl) 
  if (details) print(result)
  
  # THE STARTING TREE FOR FUSED REGULARIZATION
  # --------------------------------------------
  if(is.null(starting.tree)) starting.tree <- "tree.dev"
  tre0 <- switch(starting.tree,
                 "cart.initial" = cart.initial,
                 "cart.dev" = cart.dev,
                 "cart.bic" = cart.bic,
                 "tre0.initial" = tre0.initial,
                 "tre0.clean" = tre0.clean,
                 "tre0.dev" = tre0.dev,
                 "tre0.bic" = tre0.bic,
                 tre0.dev)    # DEFAULT
  if (is.element(starting.tree, c("cart.initial", "cart.dev", "cart.bic"))) 
    tre0 <- rpart.to.TreeFuL(tre0, dat=dat, test=test) 
  if (plot.it) {
    is.nominal <- dat %>% 
      dplyr::summarise(across(names(dat)[cols.x], ~ is.factor(.x))) %>% 
      unlist() %>% as.vector()
    cols.nominal <- cols.x[is.nominal]
    plot.tree(tre0, cols.nominal=cols.nominal, textDepth=10, leaf.label=TRUE,
              plot.score=TRUE, lines="rectangle", digits=3)
  }
  TREES <- vector("list", 8 + length(name.metrics))  # INITIALIZE LIST TREES
  TREES[1:8] <- list(cart.initial, cart.dev, cart.bic, tre0.initial, tre0.clean, 
                     tre0.dev, tre0.bic, tre0) 
  sizes <- c(sum(is.na(tre0.initial$vname)), sum(is.na(tre0.clean$vname)), 
             sum(is.na(tre0.dev$vname)), sum(is.na(tre0.bic$vname)), 
             sum(is.na(tre0$vname)))
  names(TREES) <- c("cart.initial", "cart.dev", "cart.bic", "tre0.initial", "tre0.clean", 
                    "tre0.dev", "tre0.bic", "tre0", paste("besttree", name.metrics, sep=".")) 
  
  # FUSED (ADAPTIVE LASSO) REGULARIZATION USING glmnet WITH REPARAMETERIZATION 
  # ****************************************************************************
  if (NROW(tre0) <= 3) {
    key <- lambda.selected <- metrics <- NA;
    tre0$grp <- 1; tre0$color <- "#7AD151FF"
    if (NROW(tre0) == 3) {tre0$grp <- c(NA, 1, 2); tre0$color <- c("gray50", "#7AD151FF", "#414487FF")}
    tree.size <- rep(sum(is.na(tre0$vname)), length(name.metrics))
    btree.dev <- btree.bic <- tre0
  } else {
    dat0 <- senddown(data=dat, tre=tre0) %>% dplyr::select(time, status, node)
    # SORT THE TERMINAL NODES AND PREPARE FOR glmnet
    out.hazard.sorting <- hazard.sorting(dat0)  
    key <- key.hazard.sorting <- out.hazard.sorting$key
    B.inverse <- out.hazard.sorting$B.inverse
    X <- out.hazard.sorting$X; y <- out.hazard.sorting$y 
    penalty.factor.glmnet <- out.hazard.sorting$penalty.factor.glmnet
    # print(out.hazard.sorting$K)
    
    # FIT FUSED LASSO WITH glmnet AFTER REPARAMETERIZATION 
    fit0 <- glmnet(x=X, y=y, family = "cox", relax=relaxed, standardize = FALSE,
                   penalty.factor = penalty.factor.glmnet)
    Lambda <- fit0$lambda
    n.Lambda <- length(Lambda)
    Gamma <- fit0$beta
    Beta <- B.inverse%*%Gamma |> as.matrix() |> as.data.frame() 
    complexity <- Beta %>% as.matrix() %>%  as.data.frame() %>% 
      dplyr::summarise(across(.cols=colnames(Beta), ~ length(unique(.x)))) %>% 
      unlist() %>% as.vector()
    # print(dim(Beta)); print(complexity)
    
    # VALIDATION VIA TEST SAMPLE 
    # ------------------------------
    if (identical(dat, test)) stop("No test data, but you wanna use the test sample method. R u kidding? :)")
    y.test <- Surv(time=test$time, event=test$status)
    test0 <- senddown(data=test, tre=tre0) %>% 
      dplyr::select(time, status, node) %>% 
      mutate(node=factor(node, ordered = TRUE, levels=out.hazard.sorting$leaves.sorted)) 
    contrasts(test0$node) <- contr.treatment(n=out.hazard.sorting$K, base = 1)
    X.test <- model.matrix(~node, data=test0) |> # tail(n = NROW(test)) %>% 
      as.data.frame() %>% dplyr::select(-1) %>% as.matrix()
    X.test <- X.test%*%B.inverse
    pred <- predict(fit0, newx=X.test, type = "link")
    # print(NROW(test)); print(dim(pred))
    
    #  PERFORMANCE METRICS
    # ---------------------------
    dev <- coxnet.deviance(pred = pred, y = y.test)
    # dev0 <- coxnet.deviance(x=X.test, y=y.test,  beta = Gamma)  # DEVIANCE 
    # print(cbind(dev, dev0))    # SAME 
    bic <- dev + log(sum(test$status)) * complexity 
    metrics <- data.frame(dev=dev, bic=bic) %>%  
      setNames(name.metrics)
    # print(metrics)
    lambda.selected <- tree.size <- rep(0, length(name.metrics))
    if (plot.it) par(mfrow=c(1, 2), mar=rep(5,4), xpd = FALSE)
    for (i in 1:length(name.metrics)) {
      tre.i <- tre0
      metric <- metrics[, i]
      if (is.element(i, c(1, 2))) lambda.best <- Lambda[which.min(metric)]
      else lambda.best <- Lambda[which.max(metric)]
      Beta.i <- c(0, Beta[, Lambda==lambda.best])
      grp.i <- as.numeric(as.factor(Beta.i))
      Main <- paste("(", letters[i], ")  ", name.metrics[i], sep="")
      if (plot.it) {
        plot(Lambda, metric, type="l", col="blue", lwd=1.5, main=Main)
        abline(v=lambda.best, col="orange", lty=2)
        # legend(x=quantile(Lambda, probs=0.95), y = quantile(dev, probs=0.90), 
        #       legend=c("deviance", "BIC"), fill = c("blue", "orange"))
      }
      lambda.selected[i] <- lambda.best
      tree.size[i] <- length(unique(Beta.i))   # #GROUPS
      key <- cbind(key, Beta.i, grp.i)
      
      # OBTAIN THE FINAL TREE
      recoding <- c(grp.i, NA); names(recoding) <- c(key.hazard.sorting$node, "999") 
      grp.tre <- tre.i$node; grp.tre[!is.na(tre.i$vname)] <- "999"
      tre.i$grp <- dplyr::recode(grp.tre, !!!recoding)
      tre.i <- shear(tre.i, group.name="grp")
      TREES[[i+8]] <- tre.i
    }
    if (plot.it) par(mfrow=c(1, 1))
    metrics <- cbind(Lambda, metrics) %>% as.data.frame() %>% 
      setNames(nm=c("lambda", name.metrics))
    key <- key %>% setNames(nm=c("node", "measure.sorted", "rank",  "beta.coxph", 
                                 paste(c("Beta", "grp"), rep(name.metrics, rep(2, length(name.metrics))), sep="."))) 
  }
  
  # PREPARE FOR PREDICTION
  # ******************************
  FITS <- vector("list", length(name.metrics))
  for (i in 1:length(name.metrics)) {
    tre.i <- TREES[[i+8]]
    # print(tre.i)
    if (is.null(tre.i) | NROW(tre.i) <= 1) fit.i <- coxph(Surv(time, status)~1, data=dat)
    else {
      recoding <- tre.i$grp[is.na(tre.i$vname)] 
      # print(tre.i)
      names(recoding) <- tre.i$node[is.na(tre.i$vname)]; 
      dat.i <- senddown(data=dat, tre=tre.i) %>% 
        dplyr::select(time, status, node) %>% 
        mutate(group = as.factor(dplyr::recode(as.character(node), !!!recoding))) 
      fit.i <- coxph(Surv(time, status)~group, data=dat.i)
    }
    FITS[[i]] <- fit.i
  }
  names(FITS) <- name.metrics
  tree.size <- c(size.cart, sizes, tree.size); names(tree.size) <- names(TREES) 
  list(key=key, lambda.selected=lambda.selected, tree.size=tree.size, 
       dat=dat, starting.tree=starting.tree, tree0 = tre0, 
	   TREES=TREES, metrics=metrics, FITS=FITS, 
	   formula=form0, min.nodesize=min.nodesize, cols.x=cols.x)
}




# ======================================================================
# FUNCTION SurvTreeFuL.CV() CONSTRUCTS SURVIVAL TREEFUL TREES VIA CV
# ======================================================================

SurvTreeFuL.CV <- function(dat, cols.x, mtry=length(cols.x), 
                            n.folds=5, starting.tree=NULL, plot.it=FALSE, 
                            cleaning.up=TRUE, remove.missing.median.survival=FALSE, 
                            relaxed=FALSE, details=FALSE,
                            max.depth=10, min.nodesize=20, min.event.childnode=min.nodesize/4, 
                            method.split="default-combined", stratify.status=TRUE, 
                            use.rpart=TRUE, control.bestcut=ctr0.bestcut, ...) 
{
  name.metrics <- c("dev1", "bic1", "dev2", "bic2", "cvm")
  TREES <- NULL
  
  # OBTAIN rpart TREES
  # -----------------------
  form0 <- paste("Surv(time, status) ~ ", paste(names(dat)[cols.x], collapse = " + "), sep="") |> as.formula()
  fit.rpart <- cart.CV(formula=form0, data=dat, V=n.folds) 
  rpart.initial <- fit.rpart$TREE[[1]]
  rpart.0SE <- fit.rpart$TREES[[2]]; rpart.1SE <- fit.rpart$TREES[[3]];
  TREES <- c(TREES, list(rpart.initial=rpart.initial, 
                rpart.0SE=rpart.0SE, rpart.1SE=rpart.1SE))
  if (plot.it) {
    par(mfrow=c(1, 2), mar=rep(4, 4))
    plot.cart(rpart.0SE); plot.cart(rpart.1SE)
  }
  
  # SurvTreeFuL STEP1: OBTAIN THE LARGE INITIAL TREE
  # -------------------------------------------------
  id.cv <- sample(x=1:n.folds, size=NROW(dat), replace=TRUE); 
  out.SplitCompPrune <- NA
  dat0 <- dat
  if (is.null(starting.tree)) starting.tree <- "rpart.1SE"
  if (!is.element(starting.tree, c("rpart.initial", "rpart.0SE", "rpart.1SE"))) {
    tre0 <- grow(data=dat, test=dat, cols.x=cols.x, mtry=mtry, 
                 max.depth=max.depth, min.nodesize=min.nodesize, min.event.childnode=min.event.childnode, 
                 method.split=method.split, stratify.status=stratify.status, use.rpart=use.rpart,
                 control.bestcut=control.bestcut)
    if (plot.it) {
      par(mfrow=c(1, 1), mar=rep(4, 4))
      plot.tree(tre0, cols.nominal=NULL, textDepth=10, leaf.label=TRUE,
                plot.score=TRUE, lines="rectangle", digits=3)
      title(main ="Initial Large Tree")
    }
    TREES <- c(TREES, list(initial.tree=tre0))
    if (cleaning.up) {
            tre0 <- clean.up(tre0, remove.missing.median.survival=remove.missing.median.survival); 
            TREES <- c(TREES, list(initial.tree.clean=tre0))}  # OPTIONAL CLEANING UP
    if (plot.it) {
      par(mfrow=c(1, 1), mar=rep(4, 4))
      plot.tree(tre0, cols.nominal=NULL, textDepth=10, leaf.label=TRUE,
                           plot.score=TRUE, lines="rectangle", digits=3)
      title(main ="Initial Large Tree after Clearning Up")
    }
    
    #  SurvTreeFuL STEP2: SPLIT-COMPLEXITY PRUNING PLUS VALIDATION
    # ---------------------------------------------------------------
    out0.prune <- prune.validate(dat=dat, test=dat, tre=tre0)
    alpha0 <- c(0, out0.prune$alpha) |> as.numeric()
    # Out.prune.CV <- as.list(n.folds)
    Dev.cv <- matrix(0, nrow=NROW(out0.prune), ncol=n.folds)
    for (v in 1:n.folds) {
      if (details) print(paste("============ V FOLD with v = ", v, "==================", sep=""))
      dat.v <- dat[id.cv!=v, ]; test.v <- dat[id.cv==v,]
      tre.v <- grow(data=dat.v, test=dat.v, cols.x=cols.x, mtry=mtry, 
                    max.depth=max.depth, min.nodesize=min.nodesize, min.event.childnode=min.event.childnode, 
                    method.split=method.split, stratify.status=stratify.status, use.rpart=use.rpart,
                    control.bestcut=control.bestcut)
      if (cleaning.up) tre.v <- clean.up(tre.v, 
                              remove.missing.median.survival=remove.missing.median.survival)
      outv.prune <- prune.validate(dat=dat.v, test=test.v, tre=tre.v)
      if (details) print(outv.prune)
      # Out.prune.CV[[v]] <- outv.prune
      alpha.v <- c(0, outv.prune$alpha) |> as.numeric()
      for (i0 in 1:NROW(out0.prune)) {
        if (i0==1) Dev.cv[i0, v] <- outv.prune$dev[1]
        else if (i0==NROW(out0.prune)) Dev.cv[i0, v]  <- tail(outv.prune$dev, n=1)
        else {
          alpha.i <- (alpha0[i0] + alpha0[i0+1])/2 
          pos.i <- rank(c(alpha.i, alpha.v))[1]-1  
          Dev.cv[i0, v] <- outv.prune$dev[pos.i]
        }
      }
    }
    if (details) print(Dev.cv)
    out.SplitCompPrune <- out0.prune %>% 
      mutate(dev.sum.cv=apply(Dev.cv, 1, sum),
             bic.cv=dev.sum.cv + log(sum(dat$status)) * size.tmnl,
             dev.mean.cv=apply(Dev.cv, 1, mean),
             dev.se.cv=apply(Dev.cv, 1, sd)/sqrt(n.folds)
      )
    if (details) print(out.SplitCompPrune) 
    # 0SE
    pos.btre <- which.min(out.SplitCompPrune$dev.mean.cv)
    nodes.rm <- out.SplitCompPrune$node.rm[1:(pos.btre-1)]
    # print(tre0); print(nodes.rm) 
    tre0.0SE <- obtain.tree(tre0, nodes.rm)
    size.tre0.0SE <- sum(!is.na(tre0.0SE$vname))
    # 1SE
    threshold.1SE <- min(out.SplitCompPrune$dev.mean.cv) + 
      out.SplitCompPrune$dev.se.cv[which.min(out.SplitCompPrune$dev.mean.cv)]
    pos.btre <- tail(which(out.SplitCompPrune$dev.mean.cv <= threshold.1SE), n=1) 
    nodes.rm <- out.SplitCompPrune$node.rm[1:(pos.btre-1)]
    tre0.1SE <- obtain.tree(tre0, nodes.rm) 
    size.tre0.1SE <- sum(!is.na(tre0.1SE$vname))
    # BIC
    pos.btre <- which.min(out.SplitCompPrune$bic.cv)
    nodes.rm <- out.SplitCompPrune$node.rm[1:(pos.btre-1)]
    tre0.BIC <- obtain.tree(tre0, nodes.rm)
    size.tre0.BIC <- sum(!is.na(tre0.BIC$vname))
    TREES <- c(TREES, list(tre0.0SE=tre0.0SE, tre0.1SE=tre0.1SE, 
                           tre0.BIC=tre0.BIC))
  } 
  
  # SurvTreeFuL STEP3: FUSED REGULARIZATION
  # ----------------------------------------
  if(is.null(starting.tree)) starting.tree <- "rpart.0SE"
  tre0 <- switch(starting.tree,
                 "rpart.initial" = rpart.initial,
                 "rpart.0SE" = rpart.0SE,
                 "rpart.1SE" = rpart.1SE,
                 "tre0.initial" = tre0,
                 "tre0.0SE" = tre0.0SE,
                 "tre0.1SE" = tre0.1SE,
                 "tre0.BIC" = tre0.BIC,
                 rpart.0SE)    # DEFAULT
  if (is.element(starting.tree, c("rpart.initial", "rpart.0SE", "rpart.1SE"))) 
    tre0 <- rpart.to.TreeFuL(tre0, dat=dat, test=dat) 
  if (plot.it) {
    is.nominal <- dat %>% 
      dplyr::summarise(across(names(dat)[cols.x], ~ is.factor(.x))) %>% 
      unlist() %>% as.vector()
    cols.nominal <- cols.x[is.nominal]
    plot.tree(tre0, cols.nominal=cols.nominal, textDepth=10, leaf.label=TRUE,
              plot.score=TRUE, lines="rectangle", digits=3)
  }
  
  if (NROW(tre0) <= 3) {
    key <- lambda.selected <- metrics <- NA;
    tre0$grp <- 1; tre0$color <- "#7AD151FF"
    if (NROW(tre0) == 3) {
      # USE rpart TO HANDLE THIS PECULIAR CASE 其实这一小步没必要
      form.0 <- paste("Surv(time, status) ~ ", tre0$vname[1], sep="") |> as.formula()
      # form.0 <- paste("Surv(time, status) ~ ", paste(names(dat)[cols.x], collapse = " + "), sep="") |> as.formula()
      size0 <- cart.CV(formula=form.0, data=dat, control=rpart.control(maxdepth=1, cp=0), V=10,
                       plot.it=FALSE)$tree.size; # print(size0)
      size0 <- size0[2]  # 1SE 
      if(size0 == 1) {
        # stop("The final tree is the null tree.")
        tre0 <- tre0 %>%  head(n=1) %>% 
          mutate(var=NA, vname=NA, cutoff=NA, score=NA, score.test=NA)
      } else {
        tre0$grp <- c(NA, 1, 2); tre0$color <- c("gray50", "#7AD151FF", "#414487FF") 
      }
    }
    tree.size <- rep(sum(is.na(tre0$vname)), length(name.metrics))
    TREES <- c(TREES, rep(list(tre0), length(name.metrics)))
  } else {
    # FUSED REGULARIZATION
    dat <- senddown(data=dat, tre=tre0)
    # print(table(dat$node)); print(tre0$node[is.na(tre0$vname)])
    # SORT THE TERMINAL NODES AND PREPARE FOR glmnet
    # out.hazard.sorting <- hazard.sorting(dat, tre=tre0, according="median.survival")      
    out.hazard.sorting <- hazard.sorting(dat, tre=tre0)      
    
    key <- key.hazard.sorting <- out.hazard.sorting$key
    B.inverse <- out.hazard.sorting$B.inverse
    X <- out.hazard.sorting$X; y <- out.hazard.sorting$y 
    penalty.factor.glmnet <- out.hazard.sorting$penalty.factor.glmnet
    # print(out.hazard.sorting$K)
    
    # FIT FUSED LASSO WITH glmnet AFTER REPARAMETERIZATION 
    fit0 <- glmnet(x=X, y=y, family = "cox", 
                   relax = TRUE,      # RELAXED LASSO 
                   penalty.factor = penalty.factor.glmnet, 
                   standardize = FALSE)
    # RELAXED LASSO
    # if (FALSE) fit0 <- relax.glmnet(fit=fit0, x=X, y=y, family = "cox", path = TRUE,
    #                    standardize = FALSE, penalty.factor = penalty.factor.glmnet)
    Lambda <- fit0$lambda
    n.Lambda <- length(Lambda)
    Gamma <- fit0$beta
    Beta <- B.inverse%*%Gamma |> as.matrix() |> as.data.frame() 
    complexity <- Beta %>% as.matrix() %>%  as.data.frame() %>% 
      dplyr::summarise(across(.cols=colnames(Beta), ~ length(unique(.x)))) %>% 
      unlist() %>% as.vector()
    
    # V-FOLD CV
    time.cv <- status.cv <- pred <- NULL
    # dat0 <- dat %>% dplyr::select(-node)
    dev <- matrix(0, nrow=n.Lambda, ncol=n.folds)
    for (v in 1:n.folds) {
      if (details) print(paste("============ V FOLD with v = ", v, "==================", sep=""))
      dat.v <- dat0[id.cv!=v, ]; test.v <- dat0[id.cv==v,]
      time.v <- test.v$time; status.v <- test.v$status
      ytest.v <- Surv(time=time.v, event=status.v)
      time.cv <- c(time.cv, time.v); status.cv <- c(status.cv, status.v)
      tre.v <- grow(data=dat.v, test=dat.v, cols.x=cols.x, mtry=mtry, 
                    max.depth=max.depth, min.nodesize=min.nodesize, min.event.childnode=min.event.childnode, 
                    method.split=method.split, stratify.status=stratify.status, use.rpart=use.rpart,
                    control.bestcut=control.bestcut)
      if (cleaning.up) tre.v <- clean.up(tre.v, remove.missing.median.survival=remove.missing.median.survival);
      tre.v <- pruning(tre.v, ntmnl.min=sum(is.na(tre0$vname)), print.it=FALSE) 
      # print(v); print(tre.v)
      
      # HANLDE SCENARIOS WHERE tre.v IS OF SIZE 1 OR 2
      if (NROW(tre.v) <= 3) {
        fit0.v <- coxph(Surv(time, status)~1, data=dat.v)
        pred0.v <- predict(fit0.v, newdata=test.v, type="lp")
        if (NROW(tre.v) == 1) pred.v <- matrix(rep(pred0.v, n.Lambda), ncol=n.Lambda, byrow = FALSE)
        else {
          dat.v <- senddown(data=dat.v, tre=tre.v)
          test.v <- senddown(data=test.v, tre=tre.v)
          fit1.v <- coxph(Surv(time, status)~ node, data=dat.v)
          pred1.v <- predict(fit1.v, newdata=test.v, type="lp")
          n.nulltree <- which.max(which(complexity == 1)) 
          pred.v <- matrix(c(rep(pred0.v, n.nulltree), rep(pred1.v, n.Lambda-n.nulltree)), 
                           ncol=n.Lambda, byrow = FALSE)
        }
      } else {
        dat.v <- senddown(data=dat.v, tre=tre.v)
        # print(tre.v); print(table(dat.v$node))
        out.v <- hazard.sorting(dat.v, tre=tre.v) 
        # out.v <- hazard.sorting(dat.v, tre=tre.v, according="median.survival") 
        # print(out.v$X)
        fit.v <- glmnet(x=out.v$X, y=out.v$y, family = "cox", lambda =Lambda,
                        penalty.factor = out.v$penalty.factor.glmnet, 
                        standardize = FALSE, relax = relaxed)
        # if (relaxed) fit.v <- relax.glmnet(fit=fit.v, x=out.v$X, y=out.v$y, family = "cox", 
        #                               lambda =Lambda, path = TRUE, standardize = FALSE, 
        #                               penalty.factor = out.v$penalty.factor.glmnet)
        test0.v <- senddown(data=test.v, tre=tre.v) %>% 
          dplyr::select(time, status, node) %>% 
          mutate(node=factor(node, ordered = TRUE, levels=out.v$leaves.sorted)) 
        contrasts(test0.v$node) <- contr.treatment(n=out.v$K, base = 1)
        X.test.v <- model.matrix(~node, data=test0.v) |> # tail(n = NROW(test)) %>% 
          as.data.frame() %>% dplyr::select(-1) %>% as.matrix()
        X.test.v <- X.test.v%*% (out.v$B.inverse)
        pred.v <- predict(fit.v, newx=X.test.v, type = "link")
      }
      pred <- rbind(pred, pred.v)
      dev[, v] <- coxnet.deviance(pred = pred.v, y = ytest.v)
    }
    #  PERFORMANCE METRICS
    # ----------------------
    # WAY I:  CVed DEVIANCE AND 1SE 
    dev1 <- rowSums(dev)
    bic1 <- dev1 + log(sum(dat$status)) * complexity
    # bic <- dev1 + log(NROW(dat)) * complexity
    # print(cbind(dev, bic))
    cvm <- rowMeans(dev)
    cvsd <- RowSD(dev)/n.folds 
    aSE <- min(cvm) + 1*cvsd[which.min(cvm)]          # 1SE
    lambda.best.1SE <- Lambda[cvm  <= aSE] |> head(n=1)
    # WAY II: COMPUTING DEVIANCE
    y.test <- Surv(time=time.cv, event=status.cv)
    dev2 <- coxnet.deviance(pred = pred, y = y.test)
    bic2 <- dev2 + log(sum(dat$status)) * complexity
    metrics <- data.frame(dev1=dev1, bic1=bic1, dev2=dev2, bic2=bic2, 
                          cvm=cvm) %>%  
      setNames(name.metrics)
    # print(metrics)
    lambda.selected <- tree.size <- rep(0, length(name.metrics))
    if (plot.it) par(mfrow=c(1, length(name.metrics)), mar=rep(5,4), xpd = FALSE)
    for (i in 1:length(name.metrics)) {
      tre.i <- tre0
      metric <- metrics[, i]
      lambda.best <- Lambda[which.min(metric)] 
      if (name.metrics[i] == "cvm") lambda.best <- lambda.best.1SE   # 1SE
      Beta.i <- c(0, Beta[, Lambda==lambda.best])
      grp.i <- as.numeric(as.factor(Beta.i))
      Main <- paste("(", letters[i], ")  ", name.metrics[i], sep="")
      if (plot.it) {
        plot(Lambda, metric, type="l", col="blue", lwd=1.5, main=Main)
        abline(v=lambda.best, col="orange", lty=2)
      }
      lambda.selected[i] <- lambda.best
      tree.size[i] <- length(unique(Beta.i))   # #GROUPS
      key <- cbind(key, Beta.i, grp.i)
      
      # OBTAIN THE FINAL TREE
      recoding <- c(grp.i, NA); names(recoding) <- c(key.hazard.sorting$node, "999") 
      grp.tre <- tre.i$node; grp.tre[!is.na(tre.i$vname)] <- "999"
      tre.i$grp <- dplyr::recode(grp.tre, !!!recoding)
      tre.i <- shear(tre.i, group.name="grp")
      TREES <- c(TREES, list(tre.i))
    }
    if (plot.it) par(mfrow=c(1, 1))
    metrics <- cbind(Lambda, metrics) %>% as.data.frame() %>% 
      setNames(nm=c("lambda", name.metrics))
    key <- key %>% setNames(nm=c("node", "measure.sorted", "rank",  "beta.coxph", 
                      paste(c("Beta", "grp"), rep(name.metrics, rep(2, length(name.metrics))), sep="."))) 
  }
  
  # PREPARE FOR PREDICTION
  # ----------------------------
  FITS <- vector("list", length(name.metrics))
  for (i in 1:length(name.metrics)) {
    i0 <- length(TREES) - length(name.metrics) + i
    tre.i <- TREES[[i0]]
    # print(tre.i)
    n.grp <- length(unique(na.omit(tre.i$grp)))  
    if (is.null(tre.i) | NROW(tre.i) <= 1 | n.grp <= 1) fit.i <- coxph(Surv(time, status)~1, data=dat0)
    else {
      recoding <- tre.i$grp[is.na(tre.i$vname)] 
      # print(tre.i)
      names(recoding) <- tre.i$node[is.na(tre.i$vname)]; 
      dat.i <- senddown(data=dat0, tre=tre.i) %>% 
        dplyr::select(time, status, node) %>% 
        mutate(group = as.factor(dplyr::recode(as.character(node), !!!recoding))) 
      fit.i <- coxph(Surv(time, status)~group, data=dat.i)
    }
    FITS[[i]] <- fit.i
  }
  names(TREES)[(length(TREES)-length(name.metrics)+1):length(TREES)] <- paste("besttree", name.metrics, sep=".") 
  names(tree.size) <- names(FITS) <- name.metrics
  list(out.SplitCompPrune=out.SplitCompPrune, key=key, lambda.selected=lambda.selected, tree.size=tree.size, 
       dat=dat, starting.tree=starting.tree, tree0 = tre0, 
	   TREES=TREES, metrics=metrics, FITS=FITS, 
	   formula=form0, min.nodesize=min.nodesize, cols.x=cols.x)
}










# ========================================
# FUNCTION predict.SurvTreeFuL() PREDICTS
# ========================================

predict.SurvTreeFuL <- function(fit, newdata, type="BIC", type.coxph="lp") {
  # OBTAIN TRE AND PREPARE DATA
  tre <- fit$TREES[[paste("besttree.", type, sep="")]] 
  recoding <- tre$grp[is.na(tre$vname)] 
  names(recoding) <- tre$node[is.na(tre$vname)]; 
  if (is.null(tre) && NROW(tre) <=1) {
    dat.test <- newdata %>% 
      mutate(node=1, group=as.factor(1))
  } else {
    dat.test <- senddown(data=newdata, tre=tre) %>% 
      dplyr::select(time, status, node) %>% 
      mutate(group = as.factor(dplyr::recode(as.character(node), !!!recoding))) 
  }
  # OBTAIN THE FITTED COXPH MODEL
  fit.coxph <- fit$FITS[[type]] 
  pred <- predict(fit.coxph, newdata=dat.test, type=type.coxph)
  return(list(pred=pred, Newdata=dat.test))
}




# ================================================
# FUNCTION group.summary() SUMMARIZES EACH GROUP
# ================================================

group.summary <- function(object.SurvTreeFuL, dat, best.tree = "besttree.BIC",  
                          bootstrap.bias.correction=TRUE, B=20, plot.it=FALSE, seed=NULL,
                          ...) {
  # EXTRACT INFO FROM object.SurvTreeFuL
  starting.tree <- object.SurvTreeFuL$starting.tree
  min.nodesize <-  object.SurvTreeFuL$min.nodesize
  cols.x <- object.SurvTreeFuL$cols.x
  TREES <- object.SurvTreeFuL$TREES; # names(TREES)
  btre <- TREES[[best.tree]]  # BEST TREE
  dat0 <- senddown(data=dat, tre=btre); # head(dat0)
  
  # OBTAIN SUMMARY INFO ON MEDIAN SURVIVAL
  # ----------------------------------------
  fit.0 <- survfit(Surv(time, status) ~ grp, data = dat0, 
                   stype=1, ctype=1, se.fit=TRUE, conf.int=0.95)
  # names(fit.0); names(summary(fit.0))
  P0 <- 1-(1-fit.0$"conf.int")/2   # 0.975 for 95% CI
  Z0 <- qnorm(p=P0)
  result <- survival:::survmean(fit.0,rmean="none")$matrix %>% 
    as.data.frame() %>%
    dplyr::select(-c(n.max, n.start)) %>% 
    setNames(nm=c("n", "n.events", "median", "LB", "UB")) # %>% 
  # mutate(Median = (LB+UB)/2, SE = (UB-LB)/(2*Z0),
  #       lb = median - Z0*SE, ub = median + Z0*SE)  %>%    # SLIGHT DIFFERENCES
  # dplyr::select(n, n.events, median, LB, UB)
  
  # OBTAIN SUMMARY INFO ON HAZARD RATIO
  # ------------------------------------
  fit0 <- coxph(Surv(time, status) ~ factor(grp), data = dat0)
  result0 <- coef(summary(fit0)) %>% as.data.frame() %>% 
    set_names(nm=c("beta", "HR", "SE.beta", "Z", "p.value")) %>% 
    mutate(sd.n=sqrt(NROW(dat0))*SE.beta, 
           sd.n0 = sqrt(sum(dat0$status==1))*SE.beta) %>% 
    dplyr::select(-c(sd.n0)) 
  
  # BOOTSTRAP BIAS CORRECTION TO THE sd.n FOR BETA ESTIMATES
  # ---------------------------------------------------------
  if (bootstrap.bias.correction) {
    size.btre <-  sum(is.na(btre$vname))
    depth.btre <- max(nchar(btre$node))
    size.group <- length(base::unique(na.omit(btre$grp)))
    depth0 <- depth.btre+2
    bias <- matrix(0, ncol=2, nrow=size.group-1)
    if (!is.null(seed)) set.seed(seed=seed)
    print("Starting bootstrap bias correction, b =")
    b <- 1
    while (b <= B) {
      print(b)
      # STRATIFIED BOOTSTRAP 
      dat.b <- dat %>% group_by(status) %>%
        sample_frac(1, replace = TRUE) %>% 
        as.data.frame()
      
      # CONSTRUCT A TREE FROM BOOTSTRAP SAMPLE OF LIMITED DEPTH & SIZE 
      if (is.null(starting.tree)) starting.tree <- "tree.dev"
      if (is.element(starting.tree, c("cart.initial", "cart.dev", "cart.bic"))) {
        form <- fit.SurvTreeFuL$formula
        tre0.b <- rpart(formula=form, data=dat.b, control=rpart.control(minsplit=min.nodesize, minbucket=round(min.nodesize/3), 
                                          maxdepth=depth0, cp=0, maxcompete=0, xval=0,
                                          maxsurrogate=0, usesurrogate=2, surrogatestyle=0), 
                        model=TRUE);
        tre1.b <- prune.cart(tre=tre0.b, size=size.btre)
        tre1.b <- rpart.to.TreeFuL(tre1.b, dat=dat.b) 
      } else {
        tre0.b <- grow(data=dat.b, test=dat.b, cols.x=cols.x, max.depth=depth0, 
            min.nodesize=min.nodesize, min.event.childnode=min.nodesize/4)
        tre1.b <- pruning(tre0.b, ntmnl.min=size.btre, print.it=FALSE)
      }
      
      # FUSION WITH glmnet()
      dat.b.tmp <- senddown(data=dat.b, tre=tre1.b)
      out.hazard.sorting <- hazard.sorting(dat.b.tmp, tre=tre1.b)  
      key <- out.hazard.sorting$key
      fit.glmnet <- glmnet(x=out.hazard.sorting$X, y=out.hazard.sorting$y , family = "cox", 
                           relax = TRUE,      # RELAXED LASSO 
                           penalty.factor = out.hazard.sorting$penalty.factor.glmnet, 
                           standardize = FALSE)
      Lambda <- fit.glmnet$lambda
      Beta <- (out.hazard.sorting$B.inverse)%*%(fit.glmnet$beta) |> as.matrix() |> as.data.frame() 
      complexity <- Beta %>% as.matrix() %>%  as.data.frame() %>% 
        dplyr::summarise(across(.cols=colnames(Beta), ~ length(unique(.x)))) %>% 
        unlist() %>% as.vector()
      # cbind(Lambda, complexity)
      # PICK UP THE GROUPING STRUCTURE VIA GROUP SIZE FROM btre
      lambda.best <- Lambda[which.max(complexity >= size.group)]
      Beta.b <- c(0, Beta[, Lambda==lambda.best])
      grp.b <- as.numeric(as.factor(Beta.b))
      key <- cbind(key, Beta.b, grp.b)
      
      grp <- tre1.b$node; 
      grp[!(grp %in% (key$node))] <- NA
      tre1.b$grp <- Replace(grp, old.values = key$node, new.values = key$grp.b)
      tre.b <- shear(tre1.b, group.name="grp")  # FINAL TREE FOR bTH BOOTSTRAP SAMPLE
      # print(tre.b)
      if (NROW(tre.b) >= 3) {     # PREVENT NULL TREE
        if (plot.it) {
          par(mfrow=c(1, 2))
          plot.tree(tre1.b); plot.tree(tre.b)
          par(mfrow=c(1, 1))
        }
        
        # COMPUTE BIAS 
        dat.b <- senddown(data=dat.b, tre=tre.b) %>% 
          dplyr::select(time, status, grp); # head(dat.b)
        fit.b <- coxph(Surv(time, status) ~ factor(grp), data = dat.b)
        out.b <- coef(summary(fit.b)) %>% as.data.frame() %>% 
          set_names(nm=c("beta", "HR", "SE.beta", "Z", "p.value")) %>% 
          mutate(sd.n=sqrt(NROW(dat.b))*SE.beta, 
                 sd.n0 = sqrt(sum(dat.b$status==1))*SE.beta) %>% 
          dplyr::select(-c(HR, SE.beta, Z, p.value, sd.n0)) 
        
        dat0.b <- senddown(data=dat, tre=tre.b)  %>% 
          dplyr::select(time, status, grp); 
        fit0.b <- coxph(Surv(time, status) ~ factor(grp), data = dat0.b)
        out0.b <- coef(summary(fit0.b)) %>% as.data.frame() %>% 
          set_names(nm=c("beta", "HR", "SE.beta", "Z", "p.value")) %>% 
          mutate(sd.n=sqrt(NROW(dat0.b))*SE.beta, 
                 sd.n0 = sqrt(sum(dat0.b$status==1))*SE.beta) %>% 
          dplyr::select(-c(HR, SE.beta, Z, p.value, sd.n0)) 
        
        tab <- table(dat0$grp, dat0.b$grp)
        M.prop <- prop.table(tab, margin =1)[-1, -1] 
        bias.b <- as.matrix(out0.b - out.b) 
        bias.b <- M.prop%*% bias.b
        # print(bias.b) 
        bias <- bias + bias.b
        b <- b+1
      }  
    }
    bias <- bias/B
    result0$bias.beta <- bias[,1]
    result0$bias.sd <- bias[,2]
    result0$beta.debiased <- result0$beta + bias[, 1]
    result0$sd.debaised <- result0$sd.n + bias[, 2] 
    result0 <- result0 %>% 
      add_row(beta=0, HR=1, SE.beta=NA, Z=NA, p.value=NA, sd.n=NA,
              bias.beta=NA, bias.sd=NA, beta.debiased=0, sd.debaised=NA, 
              .before = 1)
  } else {
    result0 <- result0 %>% 
      add_row(beta=0, HR=1, SE.beta=NA, Z=NA, p.value=NA, sd.n=NA,
              .before = 1)
  }
  result <- cbind(result, result0) 
  return(result)
}






# =======================================================
# WRAP-UP FUNCTION FOR TEST-SAMPLE BASED CART VIA rpart
# =======================================================
# vignette("longintro", package = "rpart")

cart.test <- function(formula, data, test, n0=NROW(data)/30, plot.it=TRUE){
  control.cart <- rpart.control(minsplit=n0, minbucket=round(n0/3), 
                                maxdepth=10, cp=0, maxcompete=0,   # NUMBER OF COMPETITIVE SPLITS
                                maxsurrogate=0, usesurrogate=2, surrogatestyle=0,  	# SURROGATE SPLITS FOR MISSING DATA
                                xval=0)  # NO V-FOLD CV NEEDED   
  tre0 <- rpart(formula=formula, data=data, control=control.cart, model=TRUE);
  CP <- tre0$cptable %>% as.data.frame %>% dplyr::select(CP) %>% unlist() %>% as.numeric()
  nsplit <- tre0$cptable %>% as.data.frame %>% dplyr::select(nsplit) %>% unlist() %>% as.numeric()
  nsubtree <- length(CP)
  dev <- rep(0, nsubtree)
  y.test <- Surv(test$time, test$status)
  for (i in 1:nsubtree) {
    tree.i <- rpart::prune(tre0, cp=CP[i])
    # dev[i] <- dev.pred.rpart(tre=tree.i, dat=data, test=test)$dev    # Bug in rpart
    pred.i <- predict(tree.i, newdata=test)    # WOULD GET DIFFERET RESULTS
    dev[i] <- glmnet::coxnet.deviance(pred = pred.i, y = y.test)
  }
  BIC <- dev + log(sum(test$status))*nsplit
  if (plot.it) {
    par.original <- par() 
    par(mfrow=c(1,2), mar=rep(4,4))
    plot(nsplit, dev, type="b", col="red", xlab="# splits")
    plot(nsplit, BIC, type="b", col="blue", xlab="# splits")
    par(par.original)
  }
  i.star <- c(which.min(dev), which.min(BIC))
  Best.cp <- CP[i.star]
  n.splits <- nsplit[i.star];
  # TO HANDLE ERROR WITH LIST OF rpart TREES
  btree.dev <- rpart::prune(tre0, cp=Best.cp[1])  
  btree.BIC <- rpart::prune(tre0, cp=Best.cp[2])  
  names(n.splits) <- c("Deviance", "BIC")
  return(list(tre0=tre0, tree.size=c(initial=sum(tre0$frame$var=="<leaf>"), n.splits+1), 
              BTREES=list(btree.dev=btree.dev, btree.BIC=btree.BIC)))
}

# =================================================================
# FUNCTION dev.pred.rpart() COMPUTED VALIDTED DEVIANCE FROM rpart
# =================================================================
dev.pred.rpart <- function(tre, dat, test) {
  dat0 <- dat %>% 
    dplyr::select(time, status) %>%
    mutate(node = as.factor(rpart.predict(tre, newdata = dat, nn=TRUE)$nn)) 
  # mutate(node=as.factor(btre.cart$where))   # SOMETIMES CAN BE WRONG
  # print(table(dat0$node))
  if (length(unique(dat0$node)) == 1) fit0 <- coxph(Surv(time, status)~1, data=dat0) 
  else fit0 <- coxph(Surv(time, status)~node, data=dat0) 
  test0 <- test %>% 
    dplyr::select(time, status) %>% 
    mutate(node = as.factor(rpart.predict(tre, newdata = test, nn=TRUE)$nn)) 
  if (!all(is.element(unique(test0$node), unique(dat0$node)))) {   # Bug in rpart !!
    print(unique(tre$where)); print(unique(dat0$node)); print(unique(test0$node)) 
  }
  pred <- predict(fit0, test0, type="lp")
  y.test <- Surv(test$time, test$status)
  dev <- glmnet::coxnet.deviance(pred = pred, y = y.test)
  # print(fit0); print(table(test0$node))
  Cstat <- survival::concordance(object=fit0, newdata=test0)$concordance
  # Cstat <- survAUC::BeggC(Surv.rsp=Surv(dat$time, dat$status), 
                               # Surv.rsp.new=y.test,
                               # lp=fit0$linear.predictors,
                               # lpnew=pred)
  # Uno's C-STAT CAN BE LESS THAN 0.5, WHICH IS HARD TO INTERPRET
  # Cstat <- survAUC::UnoC(Surv.rsp=Surv(dat$time, dat$status), 
  #                       Surv.rsp.new=y.test, lpnew=pred)
  return(list(pred=pred, dev=dev, Cstat=Cstat))
}



# ==========================================================================
# WRAP-UP FUNCTION cart.CV() FOR 0-SE & 1-SE V-FOLD CV BASED CART VIA rpart
# ==========================================================================
cart.CV <- function(formula, data, control=NULL, V=10, plot.it=FALSE){
  TREES <- as.list(1:3); tree.size <- 1:3
  if (is.null(control)) control <- rpart.control(minsplit=20, minbucket=5, 
                                                 maxdepth=10, cp=0, maxcompete=0,  # NUMBER OF COMPETITIVE SPLITS
                                                 maxsurrogate=0, usesurrogate=2, surrogatestyle=0,  	# SURROGATE SPLITS FOR MISSING DATA
                                                 xval=V)  
  tre0 <- rpart(formula=formula, data=data, control=control);
  tree.size[1] <- sum(tre0$frame$var=="<leaf>")
  TREES[[1]] <- tre0
  if (plot.it) plot.cart(tre0, color="gray80")
  if (plot.it) plotcp(tre0, minline = TRUE) # 1SE
  # 0-SE SELECTION OF TREE SIZE
  # ---------------------------
  opt <- which.min(tre0$cptable[,"xerror"])
  best.cp <- tre0$cptable[opt, "CP"]; # print(cp.best) 
  best.tree <- prune(tre0, cp = best.cp)
  TREES[[2]] <- best.tree
  btree.size <- sum(best.tree$frame$var=="<leaf>")
  tree.size[2] <-  btree.size
  if (plot.it) plot.cart(best.tree, color="gray80", main="Final 0-SE Tree") # 1SE
  # 1-SE SELECTION OF TREE SIZE
  # ---------------------------
  cv.error <- (tre0$cptable)[,4]
  SE1 <- min(cv.error) + ((tre0$cptable)[,5])[which.min(cv.error)]      # 1SE; CAN BE EASILY MODIFIED AS aSE FOR SOME a
  position <- min(which(cv.error <= SE1)); # print(position)
  # best.cp <-  sqrt(tre0$cptable[position,1]*tre0$cptable[(position-1),1]); 
  best.cp <- tre0$cptable[position,1]; 
  best.tree <- prune(tre0, cp=best.cp)
  TREES[[3]] <- best.tree
  btree.size <- sum(best.tree$frame$var=="<leaf>")
  tree.size[3] <-  btree.size
  if (plot.it) plot.cart(best.tree, color="gray80", main="Final 1-SE Tree") # 1SE
  names(TREES) <- names(tree.size) <- c("TreeInitial", "0-SE", "1-SE")
  return(list(TREES=TREES, tree.size=tree.size))
}

# ==============================================================
# FUNCTION prune.cart() PRUNES AN rpart TREE TO A NEEDED SIZE
# ==============================================================

prune.cart <- function(tre, size) {
  cp.tbl <- tre$cptable %>% as.data.frame()
  cp0 <- cp.tbl$CP[which.max(cp.tbl$nsplit>= (size-1))]
  tre <- prune.rpart(tre, cp=cp0)
  return(tre)
}



# ============================================
# PLOT CART (rpart) TREE WITH COLORED LEAVES
# ============================================
plot.cart <- function(tree, color="gray80", main="CART Tree") {
  prp(tree, extra = 1, faclen=0, nn = T, roundint=FALSE,
      fallen.leaves=FALSE, type=0, main=main,
      branch.col="darkgray",
      split.col="gray30", nn.col="gray45",
      # split.box.col = "lightgray", 
      # split.border.col = "black",
      # shadow.col = "gray", col="red",
      box.col=color) 
}


# ====================================================================
# FUNCTION predict.rpart.leaf() PREDICTS LEAD WHERE EACH OBS FALLS 
# ====================================================================
predict.rpart.leaf <- function(fit.rpart, newdata) {
  pred <- predict(fit.rpart, newdata=newdata, type = "vector") |> 
    as.numeric() |>
    round(digits = 6) |> 
    as.character()
  recoding <- row.names(fit.rpart$frame);  
  names(recoding) <- as.character(round(fit.rpart$frame$yval, digits = 6))
  # print(recoding)
  leaf <- dplyr::recode(pred, !!!recoding) |> as.integer()
  return(leaf)
}



# ===============================================================
# FUNCTION rpart.to.TreeFuL() CONVERTS rpart OBJECT INTO TreeFuL 
# ===============================================================
rpart.to.TreeFuL <- function(btre.rpart, dat, test=NULL, min.event.childnode=5){
  # PREPARE THE rpart TREE MODEL
  tre.print.info <- capture.output(print(btre.rpart), type = "output")[-(1:5)] %>% 
    gsub("< ", "<", ., fixed = TRUE)     # REMOVE SPACE AFTER <
  tre.cart <- read.table(textConnection(tre.print.info), skip=0, sep = "",
                         fill = TRUE, row.names=NULL, header=FALSE, 
                         col.names=c("node.number", "split", "n", "deviance", "yval", "leaf"))
  node.num.rpart <- tre.cart$node.number %>%  gsub(")", "", ., fixed = TRUE) %>% as.numeric()
  tre.cart$node00 <- tre.cart$node0 <- tre.cart$node <- sapply(node.num.rpart, FUN=N2binary)
  tre.cart0 <- tre.cart
  
  cond0 <- grepl("<", tre.cart$split, fixed = TRUE)
  l0 <- nchar(tre.cart$node)
  substr(tre.cart$node[cond0], l0[cond0], l0[cond0]) <- "0"
  cond0 <- grepl(">=", tre.cart$split, fixed = TRUE)
  l0 <- nchar(tre.cart$node)
  substr(tre.cart$node[cond0], l0[cond0], l0[cond0]) <- "1"
  for (i in 1:NROW(tre.cart)) {
    node0.i <- tre.cart$node0[i]
    des0.i <- de(node0.i, tre.cart0)
    if (!is.na(des0.i[1])) {
      node.i <- tre.cart$node[i]
      des.i <- tre.cart$node[is.element(tre.cart$node0, des0.i)]
      substr(des.i, 1, nchar(node.i)) <- node.i
      tre.cart$node[is.element(tre.cart$node0, des0.i)] <- des.i
    }
  }
  tre.cart <- tre.cart[order(tre.cart$node),] %>% 
    mutate(node0=node00, node00=NULL)
  # print(btre.rpart); print(tre.cart)
  
  # OBTAIN TREEFUL MODEL
  if (is.null(test)) test <- dat
  tre <- matrix(NA, nrow=NROW(tre.cart), ncol=12) 
  dat0 <- dat %>% mutate(leaf=predict.rpart.leaf(btre.rpart, newdata=dat))
  test0 <- test %>% mutate(leaf=predict.rpart.leaf(btre.rpart, newdata=test))
  tre.cart0 <- tre.cart %>% mutate(node=node0, node0=NULL)  
  for (i in 1:NROW(tre.cart)) {
    node0.i <- tre.cart$node0[i]
    leaves.i <-  tre.cart %>% 
      filter(is.element(node0, c(node0.i, de(node0.i, tre.cart0))) & leaf=="*") %>% 
      dplyr::select(node0) %>% unlist() %>% gmp::as.bigz() %>% base::strtoi(base = 2) 
    dat.i <- dat0 %>% filter(is.element(leaf, leaves.i))
    test.i <- test0 %>% filter(is.element(leaf, leaves.i))
    n <- NROW(dat.i); n.events <- sum(dat.i$status==1);
    # print(cbind(i, node0.i, tre.cart$node[i], n))
    # print(leaves.i); print(unique(dat0$leaf))
    median.surv <- median(survfit(Surv(time, status)~ 1, data=dat.i))
    ntest <- NROW(test.i); ntest.events <- sum(test.i$status==1);
    median.surv.test <- NA
    if (!is.null(ntest) && !is.na(ntest) && ntest > 5) {median.surv.test <- median(survfit(Surv(time, status)~ 1, data=test.i))}
    xsplit <- x.col <- cutoff <- score <- score.test <- NA 
    if ((tre.cart$leaf)[i]!="*") {
      child.left <- paste(tre.cart$node0[i], "0", sep = "")
      split.i <- tre.cart$split[tre.cart$node0==child.left]
      pattern.tmp <- str_extract(split.i, pattern=c(">", "<", "=")) |> na.omit() |> as.vector()
      tmp.strsplit <- base::strsplit(split.i, split=pattern.tmp)
      xsplit <- tmp.strsplit[[1]][1]  # EXTRACT SPLITTING VARIABLE vname
      x.col <- which(names(dat.i)==xsplit)  # COLUMN NUMBER OF xsplit, var
      cutoff <- gsub(pattern="=", replacement="", tmp.strsplit[[1]][2]) 
      cutoff <- ifelse(!is.na(as.numeric(cutoff)), as.numeric(cutoff), cutoff)
      # print(split.i); print(pattern.tmp); print(cbind(xsplit, cutoff))
      z <- z.create(vname=xsplit, cutoff=cutoff, dat=dat.i)
      score <- LogRank(dat.i$time, dat.i$status, z=z, min.event.childnode=min.event.childnode) 
      if (base::identical(dat.i, test.i)) score.test <- score 
      else {
        z.test <- z.create(vname=xsplit, cutoff=cutoff, dat=test.i)
        score.test <- LogRank(test.i$time, test.i$status, z=z.test, min.event.childnode=min.event.childnode) 
      }
    } 
    tre[i, ]  <- c(node=tre.cart$node[i], size=n, n.events=n.events, median.surv=median.surv, 
                   var=x.col, vname=xsplit,   cutoff=cutoff, score=score, 
                   size.test=ntest, ntest.events=ntest.events, 
                   median.surv.test=median.surv.test, score.test=score.test)  
  }
  tre <- tre %>%  as.data.frame() %>% 
    setNames(nm=c("node", "size", "n.events", "median.surv", "var", "vname", 
                  "cutoff", "score", "size.test", "ntest.events", 
                  "median.surv.test", "score.test")) 
  return(tre)
}









# ====================================================================
# FUNCTION rdat.cutpoint() GENERATES DATA FOR ONE SPLIT INVESTIGATION 
# VIA SMOOTH SIGMOID SURROGATE (SSS) VS. GREEDY SEARCH (GS)
# ====================================================================

rdat.cutpoint <- function(n=300,  beta=c(1, -0.5),
                          uniform = "discrete", K=10, digits=5, 
                          cut=.5, icensor=1, details=FALSE) 
{
  # GENERATE COVARIATES
  if (uniform == "discrete") x <- sample(x=0:K, size=n, replace=TRUE)/K
  else if (uniform == "continuous") x <- round(runif(n=n, min=0, max=1), digits =digits)
  else stop("Wrong options for the uniform= argument. Plz check!")
  
  rate <- exp(beta[1] + beta[2]*sign(x<=cut))
  xobs <- rexp(n, rate=rate); 
  if (details) print(mean(xobs))
  #### Generate observed failure time and status
  if (icensor==0) status <- 1
  else {     
    censor <- rexp(n, rate)
    status <- sign(xobs <= censor*icensor)
    xobs <- pmin(xobs, censor*icensor) 
  }
  if (details) print(mean(status))
  ##### Output
  data.frame(id=1:n, time=xobs, status=status, x=x)
}  

# set.seed(123)
# dat <- rdat.cutpoint(n=300, K=10, beta=c(1, -2), 
#                                  cut=.5, icensor=1, details=TRUE) 



# ======================================================================
# FUNCTION rdat.bias() GENERATES DATA FOR UNBIASED VARIABLE SELECTION
# ======================================================================

rdat.bias <- function(n=500, beta=c(-2, 0, 0), cutoff=.5, digits=10, 
                      add.nomial=TRUE, n.LETTERS = 5, 
                      icensor=1, details=FALSE, add.variables="binary", p0=10) 
{
  # GENERATE COVARIATES
  x1 <- sample(x=0:1, size=n, replace=TRUE)
  x2 <- round(runif(n=n, min=0, max=1), digits =digits)
  if (add.variables=="binary") X <- matrix(sample(x=0:1, size=n*p0, replace=TRUE), nrow = n, ncol=p0)
  else X <- matrix(round(runif(n=n*p0, min=0, max=1), digits = digits), nrow = n, ncol=p0)
  
  rate <- exp(beta[1] + beta[2]*x1 + beta[3]*sign(x2<=cutoff))
  xobs <- rexp(n, rate=rate); 
  if (details) print(mean(xobs))
  #### Generate observed failure time and status
  if (icensor==0) status <- 1
  else {     
    censor <- rexp(n, rate=rate)
    status <- sign(xobs <= censor*icensor)
    xobs <- pmin(xobs, censor*icensor) 
  }
  if (details) print(mean(status))
  ##### Output
  dat <- data.frame(id=1:n, time=xobs, status=status, x1=x1, x2=x2, X)
  if (add.nomial) dat <- data.frame(dat, 
                                    x3=sample(LETTERS[1:n.LETTERS], size=n, replace=TRUE))
  dat <- dat %>% mutate_if(is.character, as.factor)   
  names(dat) <- c("id", "time", "status", paste("x", 1:(NCOL(dat)-3), sep=""))
  return(dat)
}  


# =================================================================
# FUNCTION rdat.NonNull() GENERATES DATA FROM A NULL/NONNULL SETTING
# =================================================================

rdat.NonNull <- function(n=500, beta=1, icensor=1) 
{
  # GENERATE COVARIATES
  K2 <- 10; K3 <- 50; 
  x1 <- sample(x=0:1, size=n, replace=TRUE)
  x2 <- sample(x=1:K2, size=n, replace=TRUE)/K2
  x3 <- sample(x=1:K3, size=n, replace=TRUE)/K3
  x4 <- round(runif(n=n, min=0, max=1), digits =5)
  x5 <- sample(LETTERS[1:10], size=n, replace=TRUE)
  
  # CREATE DUMMY VARIABLES
  z2 <- sign(x2<=0.5); z3 <- sign(x3<=0.5); z4 <- sign(x4<=0.5); 
  z5 <- is.element(x5, LETTERS[1:5]) + 0

  rate <- exp(-beta + beta*x1 - beta*z2 + beta*z3 - beta*z4 + beta*z5)
  xobs <- rexp(n, rate=rate); 
  
  #### Generate observed failure time and status
  if (icensor==0) status <- 1
  else {     
    censor <- rexp(n, rate=rate)
    status <- sign(xobs <= censor*icensor)
    xobs <- pmin(xobs, censor*icensor) 
  }
  ##### Output
  dat <- data.frame(id=1:n, time=xobs, status=status, 
                    x1=x1, x2=x2, x3=x3, x4=x4, x5=x5)
  dat <- dat %>% mutate_if(is.character, as.factor)   
  names(dat) <- c("id", "time", "status", paste("x", 1:(NCOL(dat)-3), sep=""))
  return(dat)
}  


# ==============================================================================
# FUNCTION rdat.binary.signal() GENERATES DATA FROM A NULL SETTING WITH BINARY SIGNAL
# ==============================================================================

rdat.binary.signal <- function(n=500, beta=1, icensor=1) 
{
  # GENERATE COVARIATES
  K2 <- 10; K3 <- 50; 
  x1 <- sample(x=0:1, size=n, replace=TRUE)
  x2 <- sample(x=1:K2, size=n, replace=TRUE)/K2
  x3 <- sample(x=1:K3, size=n, replace=TRUE)/K3
  x4 <- round(runif(n=n, min=0, max=1), digits =5)
  x5 <- sample(LETTERS[1:10], size=n, replace=TRUE)

  #### Generate observed failure time and status
  rate <- exp(-beta + beta*x1)
  xobs <- rexp(n, rate=rate);   
  if (icensor==0) status <- 1
  else {     
    censor <- rexp(n, rate=rate)
    status <- sign(xobs <= censor*icensor)
    xobs <- pmin(xobs, censor*icensor) 
  }
  ##### Output
  dat <- data.frame(id=1:n, time=xobs, status=status, 
                    x1=x1, x2=x2, x3=x3, x4=x4, x5=x5)
  dat <- dat %>% mutate_if(is.character, as.factor)   
  names(dat) <- c("id", "time", "status", paste("x", 1:(NCOL(dat)-3), sep=""))
  return(dat)
}  


# ============================================================================================
# FUNCTION rdat.continuous.signal() GENERATES DATA FROM A NULL SETTING WITH CONTINUOUS SIGNAL
# ============================================================================================
rdat.continuous.signal <- function(n=500, beta=1, icensor=1) 
{
  # GENERATE COVARIATES
  x1 <- round(runif(n=n, min=0, max=1), digits = 5)
  p0 <- 9
  X <- matrix(sample(x=0:1, size=n*p0, replace=TRUE), nrow = n, ncol=p0)
  #### Generate observed failure time and status
  rate <- exp(-beta + beta*sign(x1<=0.5))
  xobs <- rexp(n, rate=rate); 
  if (icensor==0) status <- 1
  else {     
    censor <- rexp(n, rate=rate)
    status <- sign(xobs <= censor*icensor)
    xobs <- pmin(xobs, censor*icensor) 
  }
  ##### Output
  dat <- data.frame(id=1:n, time=xobs, status=status, x1=x1, X)
  dat <- dat %>% mutate_if(is.character, as.factor)   
  names(dat) <- c("id", "time", "status", paste("x", 1:(NCOL(dat)-3), sep=""))
  return(dat)
}
