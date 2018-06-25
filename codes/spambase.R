setwd("D:/Bi/utc/sy09/Project 2")
#install.packages("tree")
#install.packages(factoextra)
#install.packages("FactoMineR")
#install.packages("rpart")
#install.packages("partykit")
#install.packages("ISLR")
library(ISLR)
library(partykit)
library(rpart)
library(corrplot)
library(MASS)
library(cluster)
library("FactoMineR")
library(factoextra)
library(tree)
source("./fonctions/separ1.R")
source("./fonctions/separ2.R")
source("./fonctions/logistic.R")
source("./fonctions/anadisc.R")
source("./fonctions/mvdnorm.R")
source("./fonctions/prob.ad.R")
source("./fonctions/prob.log.R")
source("./fonctions/prob.log2.R")
source("./fonctions/interv_conf.R")

spambase <- read.csv("./donnees/spambase.csv", header = TRUE)
str(spambase)
head(spambase)
summary(spambase)
Z <- spambase[,"Z"]
table(Z)
spambase <- spambase[,c(-1,-59)]

#Analyse descritive 
head(spambase)
str(spambase)
summary(spambase)
boxplot(spambase)
correlation <- cor(spambase)
corrplot(correlation)
##X1 corrole avec tous les autres
covariance <- cov(spambase)
#ACP
res.pca <- princomp(spambase, cor = TRUE)
str(res.pca)
eig <- res.pca$sdev^2/sum(res.pca$sdev^2)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50), ncp = 15)
fviz_pca_biplot(res.pca, axes = c(1,2), habillage = Z, addEllipses = TRUE, ellipse.level = 0.95)
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
fviz_pca_var(res.pca, col.var = "black")
fviz_pca_ind(res.pca, axes = c(1,2), col.ind = "black", habillage = Z)
fviz_pca_ind(res.pca, axes = c(1,3), col.ind = "black", habillage = Z)
fviz_pca_ind(res.pca, axes = c(2,3), col.ind = "black", habillage = Z)

X <- spambase
z <- Z
for (i in 1:dim(X)[2])
{
        X[,i] <- X[,i]/max(X[,i])
}
        #Naif bayesien
taux_app <- rep(0, 100)
taux_test <- rep(0,100)
for (i in 1:100)
{
        donn.sep <- separ1(X,z)
        Xapp <- donn.sep$Xapp
        zapp <- donn.sep$zapp
        Xtst <- donn.sep$Xtst
        ztst <- donn.sep$ztst
        #apprentissage
        param <- nba.app(Xapp,zapp)
        res <- ad.val(param, Xapp)
        count <- 0
        for (j in 1:dim(Xapp)[1])
        {
                if ((!is.na(res$pred[j])) && (res$pred[j] != zapp[j]))
                        count <- count + 1
        }
        taux_app[i] <- count / dim(Xapp)[1]
        #test
        count <- 0
        #param <- nba.app(Xtst,ztst)
        res <- ad.val(param, Xtst)
        for (j in 1:dim(Xtst)[1])
        {
                if ((!is.na(res$pred[j])) && (res$pred[j] != ztst[j]))
                        count <- count + 1
        }
        taux_test[i] <- count / dim(Xtst)[1]
}
mean(taux_app)
mean(taux_test)

        #Analyse discriminante quadratique
X<-spambase
for (i in 1:dim(X)[2])
{
        X[,i] <- X[,i]/max(X[,i])
}
taux_app <- rep(0, 100)
taux_test <- rep(0,100)
for (i in 1:100)
{
        donn.sep <- separ1(X,z)
        Xapp <- donn.sep$Xapp
        zapp <- donn.sep$zapp
        Xtst <- donn.sep$Xtst
        ztst <- donn.sep$ztst
        #apprentissage
        param <- adq.app(Xapp,zapp)
        res <- ad.val(param, Xapp)
        count <- 0
        for (j in 1:dim(Xapp)[1])
        {
                if ((!is.na(res$pred[j])) && (res$pred[j] != zapp[j]))
                        count <- count + 1
        }
        taux_app[i] <- count / dim(Xapp)[1]
        #test
        count <- 0
        param <- adq.app(Xtst,ztst)
        res <- ad.val(param, Xtst)
        for (j in 1:dim(Xtst)[1])
        {
                if ((!is.na(res$pred[j])) && (res$pred[j] != ztst[j]))
                        count <- count + 1
        }
        taux_test[i] <- count / dim(Xtst)[1]
}
mean(taux_app)
mean(taux_test)
        
        #ERREUR : exp(-infi) = 0 donc f_k = 0 => ne pouvoir pas choisir la class 
        #analyse discriminante quadratique par LDA
X <- spambase
for (i in 1:dim(X)[2])
{
        X[,i] <- X[,i]/max(X[,i])
}
Xnew <- X[,c(-32,-34,-36,-40)]
z<-Z
Xx <- cbind(Xnew,z)
taux_app <- rep(0, 100)
taux_test <- rep(0,100)
formule <- paste(names(Xx)[1:53], collapse = '+')
formule <- paste(c("z", formule), collapse = '~')
formule <- as.formula(formule)
for (i in 1:100)
{
        print(i)
        donn.sep <- separ1(Xx,z)
        Xapp <- donn.sep$Xapp
        zapp <- donn.sep$zapp
        Xtst <- donn.sep$Xtst
        ztst <- donn.sep$ztst
        qda.lm <- qda(formule, data=Xapp)
        
        #apprentissage
        zpred <- predict(qda.lm, Xapp)
        for (j in 1:dim(Xapp)[1])
        {
                if (zpred$class[j] != zapp[j])
                        count <- count + 1
        }
        taux_app[i] <- count / dim(Xapp)[1]
        print(i)
        #test
        zpred <- predict(qda.lm, Xtst)
        count <- 0
        for (j in 1:dim(Xtst)[1])
        {
                if (zpred$class[j] != ztst[j])
                        count <- count + 1
        }
        taux_test[i] <- count / dim(Xtst)[1]
}
mean(taux_app)
mean(taux_test)


        #Analyse discriminante linaire
taux_app <- rep(0, 100)
taux_test <- rep(0,100)
for (i in 1:100)
{
        donn.sep <- separ1(X,z)
        Xapp <- donn.sep$Xapp
        zapp <- donn.sep$zapp
        Xtst <- donn.sep$Xtst
        ztst <- donn.sep$ztst
        #apprentissage
        param <- adl.app(Xapp,zapp)
        res <- ad.val(param, Xapp)
        count <- 0
        for (j in 1:dim(Xapp)[1])
        {
                if ((!is.na(res$pred[j])) && (res$pred[j] != zapp[j]))
                        count <- count + 1
        }
        taux_app[i] <- count / dim(Xapp)[1]
        #test
        count <- 0
        param <- adl.app(Xtst,ztst)
        res <- ad.val(param, Xtst)
        for (j in 1:dim(Xtst)[1])
        {
                if ((!is.na(res$pred[j])) && (res$pred[j] != ztst[j]))
                        count <- count + 1
        }
        taux_test[i] <- count / dim(Xtst)[1]
}
mean(taux_app)
mean(taux_test)


        #ERREUR : exp(-infi) = 0 donc f_k = 0 => ne pouvoir pas choisir la class 
        #analyse discriminante lineaire par LDA
X <- spambase
for (i in 1:dim(X)[2])
{
        X[,i] <- X[,i]/max(X[,i])
}
Xx <- cbind(X,z)
taux_app <- rep(0, 100)
taux_test <- rep(0,100)
for (i in 1:100)
{
        donn.sep <- separ1(Xx,z)
        Xapp <- donn.sep$Xapp
        zapp <- donn.sep$zapp
        Xtst <- donn.sep$Xtst
        ztst <- donn.sep$ztst
        lda.lm <- lda(z~., data=Xapp)
        
        #apprentissage
        zpred <- predict(lda.lm, Xapp)
        for (j in 1:dim(Xapp)[1])
        {
                if (zpred$class[j] != zapp[j])
                        count <- count + 1
        }
        taux_app[i] <- count / dim(Xapp)[1]
        
        #test
        zpred <- predict(lda.lm, Xtst)
        count <- 0
        for (j in 1:dim(Xtst)[1])
        {
                if (zpred$class[j] != ztst[j])
                        count <- count + 1
        }
        taux_test[i] <- count / dim(Xtst)[1]
}
mean(taux_app)
mean(taux_test)

        #analyse logistic avec library ISLR
X<-spambase
for (i in 1:dim(X)[2])
{
        X[,i] <- X[,i]/max(X[,i])
}
zprime <- Z
zprime[which(zprime==2)] <- 0
Xx <- cbind(X,zprime)
Xx <- Xx[,-32]
taux_app <- rep(0, 100)
taux_test <- rep(0,100)
formule <- paste(names(Xx)[1:56], collapse = '+')
formule <- paste(c("zprime", formule), collapse = '~')
formule <- as.formula(formule)
for (i in 1:100)
{
        print(i)
        donn.sep <- separ1(Xx,z)
        Xapp <- donn.sep$Xapp
        zapp <- donn.sep$zapp
        Xtst <- donn.sep$Xtst
        ztst <- donn.sep$ztst
        log <- glm(formule, data=Xapp, family = binomial())
        
        #apprentissage
        count<-0
        zpred <- predict(log, Xapp, type = "response")
        zpred <- ifelse(zpred>0.5,1,2)
        for (j in 1:dim(Xapp)[1])
        {
                if (zpred[j] != zapp[j])
                        count <- count + 1
        }
        taux_app[i] <- count / dim(Xapp)[1]
        print(i)
        #test
        zpred <- predict(log, Xtst, type = "response")
        zpred <- ifelse(zpred>0.5,1,2)
        
        count <- 0
        for (j in 1:dim(Xtst)[1])
        {
                if (zpred[j] != ztst[j])
                        count <- count + 1
        }
        taux_test[i] <- count / dim(Xtst)[1]
}
mean(taux_app)
mean(taux_test)



#Analyse logistic avec intercept = 0
taux_test <- rep(0,100)
taux_app <- rep(0,100)
for (i in 1:100)
{
        donn.sep <- separ1(X,z)
        Xapp <- donn.sep$Xapp
        zapp <- donn.sep$zapp
        Xtst <- donn.sep$Xtst
        ztst <- donn.sep$ztst
        #apprentissage
        param <- log.app(Xapp, zapp, intr = FALSE, epsi = 1e-5)
        count <- 0
        res <- log.val(param$beta, Xapp)
        for (j in 1:dim(Xapp)[1])
        {
                if (res$pred[j] != zapp[j])
                        count <- count + 1
        }
        taux_app[i] <- count /dim(Xapp)[1]
        #test
        count <- 0
        res <- log.val(param$beta, Xtst)
        for (j in 1:dim(Xtst)[1])
        {
                if (res$pred[j] != ztst[j])
                        count <- count + 1
        }
        taux_test[i] <- count /dim(Xtst)[1]
}
mean(taux_app)
mean(taux_test)

        #Arbre de decision 
X <- spambase
for (i in 1:dim(X)[2])
{
        X[,i] <- X[,i]/max(X[,i])
}
taux_app <- rep(0,100)
taux_test <- rep(0,100)
for ( i in 1:100) 
{
        donn.sep <- separ1(X,Z)
        Xapp <- donn.sep$Xapp
        zapp <- donn.sep$zapp
        Xtst <- donn.sep$Xtst
        ztst <- donn.sep$ztst;
        zpred <- rep(0,length(ztst))
        zapp <- as.factor(zapp)
        
        rpart.tree <- rpart(zapp ~ ., data=Xapp)
        prune_tree <- prune(rpart.tree, cp=1e-5)
        
        pred <- predict(prune_tree, Xapp)
        zpred <- max.col(pred)
        
        count <- 0
        for (j in 1:dim(Xapp)[1])
        {
                if (zpred[j] != zapp[j])
                        count <- count + 1
        }
        taux_app[i] <- count /dim(Xapp)[1]
        
        
        pred <- predict(prune_tree, Xtst)
        zpred <- max.col(pred)
        
        count <- 0
        for (j in 1:dim(Xtst)[1])
        {
                if (zpred[j] != ztst[j])
                        count <- count + 1
        }
        taux_test[i] <- count /dim(Xtst)[1]
}
mean(taux_app)
mean(taux_test)

#tracer tree
rparty.tree <- as.party(prune_tree)
plot(rparty.tree)
