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
hist(y, main = "Riboflavin Production (log scale)", xlab = "y") #The distribution of Y is skewed to the left. But that's ok 'cause neither LASSO nor OLS requires normal distribution of error or Y
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

##NOTE: In the current version of code I do not divide training set and test set. Neither did I check the distribution of errors. These were expected to added later.
##Also according to further research, not all selected genes belong to the riboflavin biosynthetic pathway. These genes should therefore be considered predictive biomarkers and candidate regulators rather than direct causal genes.
##Mediatior analysis could probably be the next step.
