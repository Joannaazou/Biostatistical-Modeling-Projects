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

#FURTHER VISUALIZATION
#PLot1. Bar chart - coefficient sizes of selected genes
# Extract coefficients at lambda.1se
coef_1se_train <- coef(cv_fit_train, s = "lambda.1se")

#Convert to data frame, drop intercept, drop zero coefficients
coef_df_train <- data.frame(
  gene = rownames(coef_1se_train)[-1],           # remove Intercept row
  coefficient = as.vector(coef_1se_train)[-1]
)

# Keep only non-zero coefficients
coef_df_train <- coef_df_train[coef_df_train$coefficient != 0, ]

# Sort by coefficient size for a cleaner plot
coef_df_train <- coef_df_train[order(coef_df_train$coefficient), ]
coef_df_train$gene <- factor(coef_df_train$gene, levels = coef_df_train$gene)

# Plot
library(ggplot2)
ggplot(coef_df_train, aes(x = gene, y = coefficient, fill = coefficient > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +   # horizontal bars so gene names don't overlap
  scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
                    labels = c("Negative", "Positive"), name = "") +
  labs(title = "Genes Selected by Lasso and Their Coefficients (lambda.1se)",
       x = "Gene", y = "Coefficient", 
       caption = "**Positive coefficient value indicates contribution to riboflavin production") +
  theme_minimal()

## A lot of coefficients are negative? Why? This may reflect competition for metabolic resources, negative feedback regulation, or co-expression with the true causal genes.
## Which means statistics can not directly give the genes that "cause" the production of riboflavin
## Also, I did not conduct stability checks (e.g., repeated data splitting or bootstrap selection frequency, as mentioned earlier). So the risk of false positives among these negatively-associated genes is not negligible.

#Plot2. Heatmap
library(pheatmap) 

# Names of selected genes
selected_genes_train <- coef_df_train$gene

# Extract expression data for just these genes from the original matrix
# x is samples (rows) x genes (columns), so transpose since heatmaps put genes in rows
expr_selected <- t(x[, selected_genes_train])

# Optional: add y (riboflavin production) as a top annotation to see if expression
# pattern lines up with production level
annotation_col <- data.frame(riboflavin_y = y)
rownames(annotation_col) <- rownames(x)

pheatmap(expr_selected,
         scale = "row",
         annotation_col = annotation_col,
         show_colnames = FALSE,
         main = "Expression Heatmap of Lasso-Selected Genes",
         color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
         fontsize = 11,                    
         fontsize_row = 10,                
         border_color = NA,                
         treeheight_row = 40,              
         treeheight_col = 30,              
         cellwidth = 6,                    
         cellheight = 16,                 
         angle_col = 45)   

## This is messy, not that clear

#Plot3. Coefficient path plot — shows the full lasso selection process
fit_train <- glmnet(x_train, y_train, alpha = 1)

plot(fit_train, xvar = "lambda", label = TRUE)
abline(v = -log(cv_fit_train$lambda.min), lty = 2, col = "gray40")
abline(v = -log(cv_fit_train$lambda.1se), lty = 2, col = "gray40")
legend("topright", legend = c("lambda.min", "lambda.1se"), lty = 2, col = "gray40")
