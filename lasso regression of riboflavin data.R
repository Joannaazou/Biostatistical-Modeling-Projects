library(ggplot2)
iris <- iris
ggplot(iris, aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point(size = 2.5, alpha = 0.7) + 
  geom_smooth(method = "lm", se = TRUE) + 
  labs(title = "Sepal Length vs. Petal Length by Species", 
       x = "Sepal Length", 
       y = "Petal Length") +
  theme_minimal()

setosa_iris <- subset(iris, iris$Species == "setosa")
ggplot(setosa_iris, aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point(size = 2.5, alpha = 0.7) + 
  geom_smooth(method = "lm", se = TRUE) + 
  labs(title = "Sepal Length vs. Petal Length of Setosa", 
       x = "Sepal Length", 
       y = "Petal Length") +
  theme_minimal()

setosa_iris_lm <- lm(Petal.Length ~ Sepal.Length, data = setosa_iris)
