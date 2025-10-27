
rm(list=ls(all=TRUE), envir = environment())
source("Functions-SurvTreeFuL.R") 

# SIMULATE DATA FROM TREE2 MODEL
set.seed(123)
beta0 <- c(-1, 3) 
p0 <- 1
dat <- rdat(n=600, beta=beta0, model="tree2", p0=p0)
head(dat); dim(dat); names(dat)


# FIT SURVTREEFUL WITH V-FOLD CV 
# --------------------------------
# SPECIFY COLS OF PREDICTORS 
cols.x <- 4:7
fit <- SurvTreeFuL.CV(dat=dat, cols.x=cols.x, n.folds=5, 
                      # starting.tree = "rpart.1SE", 
                      starting.tree = "tre0.0SE",
                      max.depth=6, min.nodesize=20, min.event.childnode=5)
TREES <- fit$TREES; 
SIZE <- fit$tree.size; SIZE
# names(TREES); names(FITS)

# CART TREES 
# --------------
plot.cart(TREES$rpart.initial)  # LARGE INITIAL rpart TREE
plot.cart(TREES$rpart.1SE)  # BEST rpart TREE BY 1-SE


# K-M SURVIVAL CURVES
btre.CART <- TREES$rpart.1S
dat1 <- dat %>% 
  mutate(leaf.CART=predict.rpart.leaf(btre.CART, dat))

library(ggsurvfit)
survfit2(Surv(time, status) ~ leaf.CART, data = dat1) %>% 
  ggsurvfit() +
  labs(x = "Time", y = "Survival") +
  ggtitle("K-M Survival Curves") +
  theme(plot.title = element_text(hjust = 0.5))



# ONE BEST SurvTreeFuL MODEL
# -----------------------------

name.btre <- "besttree.bic1"
btre <- TREES[[name.btre]] 
plot.tree(btre)
plot.tree(btre, plot.score = FALSE)

# K-M SURVIVAL CURVES
dat2 <- senddown(dat, btre); names(dat2) 
par(mfrow=c(1, 2))
plot(survfit(Surv(time, status) ~ node, data = dat2), 
     main="KM Curves for Nodes")
plot(survfit(Surv(time, status) ~ grp, data = dat2), 
     main="KM Curves for Groups")

# USING ggplot2
library(ggsurvfit)
fig1 <- survfit2(Surv(time, status) ~ node, data = dat2) %>% 
  ggsurvfit() +
  labs(x = "Time", y = "Survival") +
  ggtitle("K-M Survival Curves") +
  theme(plot.title = element_text(hjust = 0.5))

fig2 <- survfit2(Surv(time, status) ~ grp, data = dat2) %>% 
  ggsurvfit() +
  labs(x = "Time", y = "Survival") +
  ggtitle("K-M Survival Curves") +
  theme(plot.title = element_text(hjust = 0.5))

library(gridExtra)
grid.arrange(fig1, fig2, nrow = 1)


# BOOTSTRAP BIAS CORRECTION FOR GROUP SUMMARY
# -----------------------------------------------

B <- 25
out <- group.summary(object.SurvTreeFuL=fit, dat, best.tree = name.btre,
              bootstrap.bias.correction=TRUE, B=B, seed=NULL) 
out
out %>% 
  # dplyr::select(beta, sd.n, beta.debiased, sd.debaised) %>% 
  # slice(-1)   %>%  # REMOVE THE FIRST ROW
  mutate(HR.debiased = exp(beta.debiased),
         SE.debiased=sd.debaised/sqrt(n), 
         z.debiased=beta.debiased/SE.debiased, 
         pvalue.debiased = pchisq(z.debiased^2, df=1, lower.tail = FALSE))







#
