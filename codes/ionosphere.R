setwd("D:/Bi/utc/sy09/Project 2")
#install.packages("tree")
#install.packages(factoextra)
#install.packages("FactoMineR")
#install.packages("rpart")
#install.packages("partykit")

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

ionosphere <- read.csv("./donnees/ionosphere.csv", header = TRUE)
str(ionosphere)
head(ionosphere)
summary(ionosphere)
Z <- ionosphere[,"Z"]
table(Z)
ionosphere <- ionosphere[,c(-1,-36)]

#Analyse descritive 
head(ionosphere)
str(ionosphere)
summary(ionosphere)
        ## X2 toujours 0 => supprimer
ionosphere <- ionosphere[,-2]
boxplot(ionosphere)
correlation <- cor(ionosphere)
corrplot(correlation)
        ##X1 corrole avec tous les autres
covariance <- cov(ionosphere)
#ACP
res.pca <- princomp(ionosphere, cor = TRUE)
str(res.pca)
eig<-res.pca$sdev^2/sum(res.pca$sdev^2)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50), ncp = 20)
fviz_pca_biplot(res.pca, axes = c(1,2), habillage = Z, addEllipses = TRUE, ellipse.level = 0.95)
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
fviz_pca_var(res.pca, col.var = "black", axes = c(2,3))
fviz_pca_ind(res.pca, axes = c(1,2), col.ind = "black", habillage = Z)
fviz_pca_ind(res.pca, axes = c(1,3), col.ind = "black", habillage = Z)
fviz_pca_ind(res.pca, axes = c(2,3), col.ind = "black", habillage = Z)

#Hiearachie cluster
rec_scaled <- scale(as.matrix(ionosphere))
dissimilarite <- dist(rec_scaled, diag = TRUE)
breast.hc <- hclust(dissimilarite, method = "ward.D2")
plot(breast.hc, main = "HC")

#kmeans cluster
rec_scaled <- scale(as.matrix(ionosphere))
dissimilarite <- dist(rec_scaled, diag = TRUE)
breast_kmeans <- kmeans(dissimilarite, 2)
clusplot(ionosphere, breast_kmeans$cluster)
breast_kmeans <- kmeans(dissimilarite, 3)
clusplot(ionosphere, breast_kmeans$cluster)
breast_kmeans <- kmeans(dissimilarite, 4)
clusplot(ionosphere, breast_kmeans$cluster)
breast_kmeans <- kmeans(dissimilarite, 5)
clusplot(ionosphere, breast_kmeans$cluster)
breast_kmeans <- kmeans(dissimilarite, 6)
clusplot(ionosphere, breast_kmeans$cluster)
breast_kmeans <- kmeans(dissimilarite, 7)
clusplot(ionosphere, breast_kmeans$cluster)
breast_kmeans <- kmeans(dissimilarite, 8)
clusplot(ionosphere, breast_kmeans$cluster)
breast_kmeans <- kmeans(dissimilarite, 9)
clusplot(ionosphere, breast_kmeans$cluster)

X <- ionosphere
z <- Z
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
                if (res$pred[j] != zapp[j])
                        count <- count + 1
        }
        taux_app[i] <- count / dim(Xapp)[1]
        #test
        count <- 0
        param <- nba.app(Xtst,ztst)
        res <- ad.val(param, Xtst)
        for (j in 1:dim(Xtst)[1])
        {
                if (res$pred[j] != ztst[j])
                        count <- count + 1
        }
        taux_test[i] <- count / dim(Xtst)[1]
}
mean(taux_app)
mean(taux_test)
        #Analyse discriminante quadratique
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
                if (res$pred[j] != zapp[j])
                        count <- count + 1
        }
        taux_app[i] <- count / dim(Xapp)[1]
        #test
        count <- 0
        param <- adq.app(Xtst,ztst)
        res <- ad.val(param, Xtst)
        for (j in 1:dim(Xtst)[1])
        {
                if (res$pred[j] != ztst[j])
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
                if (res$pred[j] != zapp[j])
                        count <- count + 1
        }
        taux_app[i] <- count / dim(Xapp)[1]
        #test
        count <- 0
        param <- adl.app(Xtst,ztst)
        res <- ad.val(param, Xtst)
        for (j in 1:dim(Xtst)[1])
        {
                if (res$pred[j] != ztst[j])
                        count <- count + 1
        }
        taux_test[i] <- count / dim(Xtst)[1]
}
mean(taux_app)
mean(taux_test)

        #Analyse logistic avec intercept = 0
Xprime <- res.pca$scores[,1:25]

taux_test <- rep(0,100)
taux_app <- rep(0,100)
for (i in 1:100)
{
        donn.sep <- separ1(Xprime,z)
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

        #Analyse logistique avec intercept = 1
Xprime <- res.pca$scores[,1:25]
taux_app <- rep(0,100)
taux_test <- rep(0,100)
for (i in 1:100)
{
        donn.sep <- separ1(Xprime,z)
        Xapp <- donn.sep$Xapp
        zapp <- donn.sep$zapp
        Xtst <- donn.sep$Xtst
        ztst <- donn.sep$ztst
        #apprentissage
        param <- log.app(Xapp, zapp, intr = TRUE, epsi = 1e-5)
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
        #Analyse logistique quadratique  avec intercept = 0
Xprime <- res.pca$scores[,1:6]
Xnew <- quadratique(Xprime)
taux_app <- rep(0,100)
taux_test <- rep(0,100)
for (i in 1:100)
{
        donn.sep <- separ1(Xnew,z)
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
        taux_test[i] <- count /dim(Xapp)[1]
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
mean(taux_test)

        #Abre de decision


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

        # graphique
donn.sep <- separ1(X,Z)
Xapp <- donn.sep$Xapp
zapp <- donn.sep$zapp
Xtst <- donn.sep$Xtst
ztst <- donn.sep$ztst;
zpred <- rep(0,length(ztst))
zapp <- as.factor(zapp)
rpart.tree <- rpart(zapp ~ ., data=Xapp)
prune_tree <- prune(rpart.tree, cp=0.02)
rparty.tree <- as.party(prune_tree)
plot(rparty.tree)

