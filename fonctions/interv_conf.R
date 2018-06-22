interv_conf <- function(vector){
        interv <- c(0,0,0);
        mu <- mean(vector);
        s <- sd(vector);
        n <- length(vector);
        interv[1]<- mu
        interv[2]<-mu-1.96* (s/sqrt(n))
        interv[3]<-mu+1.96* (s/sqrt(n))
        
        interv
}# renvoie moyenne, borne inf, borne sup.
