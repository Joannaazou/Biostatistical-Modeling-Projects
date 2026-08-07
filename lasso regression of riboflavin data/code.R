library(ScaleSpikeSlab)
library(glmnet)
data(riboflavin)

#learn about the data first -> the reason why the riboflavin production takes on negative value is that they are in log scale
help(riboflavin)

#check the structure of riboflavin
str(riboflavin)

#I found a "dinames" symbol in the console, so I check the colnames and the rownames of this data frame.
colnames(riboflavin)
rownames(riboflavin)
colnames(riboflavin$x) 
rownames(riboflavin$x) #Now I see what "dinames" mean
colnames(riboflavin$y) #It seems like riboflavin$y does not have dinames
rownames(riboflavin$y)

#convert the "data frame"'s first column to the response variable and the second column into the predictor variable matrix
y <- as.vector(riboflavin$y)
x <- as.matrix(riboflavin$x)

#Here I check the class and dimension of them again
length(y)
dim(x)
class(y)
class(x)

#visualize the data first
hist(y, main = "Riboflavin Production (log scale)", xlab = "y", breaks = 20) #The distribution of Y is skewed to the left. But that's ok 'cause neither LASSO nor OLS requires normal distribution of error or Y
summary(y)
x[1:5, 1:5]

#perform cross-validation to select the optimal lambda using cv.glmnet
#first set the random seed in order to make this program replicable
set.seed(42)
cv_fit <- cv.glmnet(x, y, alpha = 1, family = "gaussian", nfolds = 10) #alpha = 1 indicates the lasso method, and gaussian indicates that this is a linear regression(not logistic) with continuous variable Y

#choose lambda(visualization)
plot(cv_fit)

#choose lambda value and see how many genes they have chosen
cv_fit$lambda.min
cv_fit$lambda.1se
cv_fit #I chose 1se 'cause it includes lesser gene than min (27<40)

#extract the name of the selected genes
coef_1se <- coef(cv_fit, s = "lambda.1se")
selected_genes_1se <- rownames(coef_1se)[which(coef_1se != 0)]

#see the genes and the number of genes
selected_genes_1se
length(selected_genes_1se) - 1

#evaluate the lasso model
#1. look at the cv_mse ad compare it with the null model mse
cv_mse <- cv_fit$cvm[cv_fit$lambda == cv_fit$lambda.1se]
null_mse <- var(y)
null_mse
cv_mse
cv_mse/null_mse

#2.calculate the R^2 value
pred_1se <- predict(cv_fit, newx = x, s = "lambda.1se")
r_squared <- 1 - sum((y - pred_1se)^2)/sum((y - mean(y))^2)
r_squared #R^2 >0.5 which is good


##Note: According to further research, not all selected genes belong to the riboflavin biosynthetic pathway. These genes should therefore be considered predictive biomarkers and candidate regulators rather than direct causal genes.
##Mediatior analysis could probably be the next step. 

####################################################################
############LASSO REGRESSION WITH TRAINING&TEST SETS################
####################################################################
set.seed(221)
n <- nrow(x)
train_idx <- sample(1:n, size = round(0.75 * n))  # 85% for training set

x_train <- x[train_idx, ]
y_train <- y[train_idx]
x_test  <- x[-train_idx, ]
y_test  <- y[-train_idx]

#check sample size
length(train_idx)   # sample size in the training set
n - length(train_idx)  # sample size in the test set

#visualize the data of the training set first
hist(y_train, main = "Riboflavin Production of the Training Set(log scale)", xlab = "y", breaks = 20) 
summary(y_train)
x_train[1:5, 1:5]

#run the cross validation lasso regression model(repeat the steps above)
set.seed(232)
cv_fit_train <- cv.glmnet(x_train, y_train, alpha = 1, family = "gaussian", nfolds = 10)

plot(cv_fit_train)
cv_fit_train
best_lambda_train <- cv_fit$lambda.1se #only 15 genes selected!

# predict on the test set and evaluate
pred_test <- predict(cv_fit_train, s = best_lambda_train, newx = x_test)
pred_test
mse_test <- mean((y_test - pred_test)^2)
sst_test <- sum((y_test - mean(y_test))^2)
sse_test <- sum((y_test - pred_test)^2)
r2_test <- 1 - sse_test/sst_test

mse_test
r2_test #r2>0.5 which is fair
