
c <- read.csv(file = 'outs/corr_gene2.csv',header=FALSE)
c <- as.data.frame(c)
df2 <- as.matrix(c)

corrplot(df2)
corrplot(df2, method = 'square')