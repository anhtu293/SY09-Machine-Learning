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
source("./fonctions/interv_conf.R")

sonar <- read.csv("./donnees/sonar.csv", header = TRUE)
str(sonar)
head(sonar)
summary(sonar)
Z <- sonar[,"Z"]
table(Z)
sonar <- sonar[,c(-1,-62)]

#Analyse descritive 
head(sonar)
str(sonar)
summary(sonar)
boxplot(sonar)
correlation <- cor(sonar)
corrplot(correlation)
covariance <- cov(sonar)
        #ACP
res.pca <- princomp(sonar, cor = TRUE)
str(res.pca)
eig <- res.pca$sdev^2/sum(res.pca$sdev^2)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50), ncp = 20)
fviz_pca_biplot(res.pca, axes = c(1,3), habillage = Z, addEllipses = TRUE, ellipse.level = 0.95)
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
fviz_pca_var(res.pca, axes = c(1,3), col.var = "black")
fviz_pca_ind(res.pca, axes = c(1,2), col.ind = "black", habillage = Z)
fviz_pca_ind(res.pca, axes = c(1,3), col.ind = "black", habillage = Z)
fviz_pca_ind(res.pca, axes = c(2,4), col.ind = "black", habillage = Z)

#Hiearachie cluster
rec_scaled <- scale(as.matrix(sonar))
dissimilarite <- dist(rec_scaled, diag = TRUE)
breast.hc <- hclust(dissimilarite, method = "ward.D2")
plot(breast.hc, main = "HC")

#kmeans cluster
rec_scaled <- scale(as.matrix(sonar))
dissimilarite <- dist(rec_scaled, diag = TRUE)
breast_kmeans <- kmeans(dissimilarite, 2)
clusplot(sonar, breast_kmeans$cluster)
breast_kmeans <- kmeans(dissimilarite, 3)
clusplot(sonar, breast_kmeans$cluster)
breast_kmeans <- kmeans(dissimilarite, 4)
clusplot(sonar, breast_kmeans$cluster)
breast_kmeans <- kmeans(dissimilarite, 5)
clusplot(sonar, breast_kmeans$cluster)
breast_kmeans <- kmeans(dissimilarite, 6)
clusplot(sonar, breast_kmeans$cluster)
breast_kmeans <- kmeans(dissimilarite, 7)
clusplot(sonar, breast_kmeans$cluster)
breast_kmeans <- kmeans(dissimilarite, 8)
clusplot(sonar, breast_kmeans$cluster)
breast_kmeans <- kmeans(dissimilarite, 9)
clusplot(sonar, breast_kmeans$cluster)

X <- sonar
z <- Z
        #Naif bayesien
Xnew <- res.pca$scores[,1:25]
taux_app <- rep(0, 100)
taux_test <- rep(0,100)
for (i in 1:100)
{
        donn.sep <- separ1(Xnew,z)
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
        #param <- nba.app(Xtst,ztst)
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
Xnew <- res.pca$scores[,1:25]
taux_app <- rep(0, 100)
taux_test <- rep(0,100)
for (i in 1:100)
{
        donn.sep <- separ1(Xnew,z)
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
        #param <- adq.app(Xtst,ztst)
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

#trop de parametre (g=60) => overfitting. utiliser pca, choisir 30 composants => conserver 95% informations
# => bon resultat malgre aux mauvais resultats de naif(trop simple pour apprentissage)
#ADL marche les donnees itinitales et donne le resultat aceptable car assez complique et pas trop simple
        #Analyse discriminante linaire
taux_app <- rep(0, 100)
taux_test <- rep(0,100)
for (i in 1:100)
{
        donn.sep <- separ1(Xnew,z)
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
        #param <- adl.app(Xtst,ztst)
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

        #essayer de diminuer les parametresm, un peu mieux mais n'est pas remarque

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
        ##ne pas tres marquer

        #Analyse logistique quadratique  avec intercept = 0
Xprime <- res.pca$scores[,1:6]
Xnew <- quadratique(X)
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

#Analyse logistique quadratique avec intercept = 1
Xprime <- res.pca$scores[,1:4]
Xnew <- quadratique(X)
taux_test <- rep(0,100)
for (i in 1:100)
{
        donn.sep <- separ1(Xnew,z)
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
        donn.sep <- separ1(as.data.frame(X),Z)
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
