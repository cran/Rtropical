## ----knitr-setup, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(Rtropical)
library(ape)

## ----trop-svm-setup-----------------------------------------------------------
set.seed(101)
data(sim_trees)
treevecs = do.call("rbind", lapply(sim_trees, as.vector))
labels = as.factor(rep(c(1, 2), each = nrow(treevecs)/2))
# generate training data set
trn_ind = sample(1: nrow(treevecs), nrow(treevecs)*0.8)
x = treevecs[trn_ind, ]
y = labels[trn_ind]

# generate testing data set
newx = treevecs[-trn_ind, ]
newy = labels[-trn_ind]

## ----trop-svm-----------------------------------------------------------------
# run tropical svm
start = Sys.time()
trop_fit <- tropsvm(x, y, auto.assignment = TRUE)
end = Sys.time()
# predict for testing data
trop_pred <- predict(trop_fit, newx)
# compute classification accuracy
sum(as.vector(trop_pred) == newy)/length(newy)
print(paste("The running time is: ", round(end - start, digits = 3), "s", sep = ""))

## ----cv-trop-svm--------------------------------------------------------------
# tropical svm with cross-validation
start = Sys.time()
cv_trop_fit <- cv.tropsvm(x, y, nassignment = 100, parallel = TRUE)
end = Sys.time()
cv_trop_pred <- predict(cv_trop_fit, newx)
# compute classification accuracy for testing data
sum(cv_trop_pred == newy)/length(newy)
print(paste("The running time is: ", round(end - start, digits = 3), "min", sep = ""))

## ----svm-e1071----------------------------------------------------------------
svm_fit <- e1071::svm(x, y)
svm_pred <- predict(svm_fit, newx)
sum(svm_pred == newy)/length(newy)

## ----trop-pca-----------------------------------------------------------------
data(apicomplexa)
treevecs <- as.matrix(apicomplexa, parallel = TRUE)
pca_fit = troppca.poly(treevecs)

## ----plot-trop-pca------------------------------------------------------------
plot(pca_fit, fw = T)

## ----trop-fw------------------------------------------------------------------
tropFW(treevecs)

