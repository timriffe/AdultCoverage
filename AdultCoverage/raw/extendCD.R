# TODO: make a data object that we can work with...
mx <- cdmltw()$nmx
extend <- apply(cdmltw()$nmx[,c("65","70","75","80","85","90","95")],1,function(y){
    x        <- seq(65, 95, by = 5)
    y2       <- predict(lm(log(y) ~ x), newdata = list(x = c(95, 100, 105, 110)))
    y2[2:4]  <- y2[1] + cumsum(diff(y2) / c(1.5, 2, 3)) 
    exp(y2[2:4])
  })
rownames(extend) <- c(100,105,110)
graduate <- cbind(mx[,ncol(mx)],t(extend))
y        <- graduate[1,]
mxgrad <- apply(graduate,1,function(y){
    exp(approx(c(95,100,105,110),log(y),xout = 95:114,rule=2)$y)
  })

ax <- mxgrad * 0 + .5

qx <- mxgrad / (1 + (1- ax) * mxgrad)
qx[nrow(qx),] <- 1
lx <- rbind(1, apply(1-qx,2,cumprod))
Lx <- (lx[1:(nrow(lx)-1),] + lx[2:nrow(lx),]) / 2

e95 <- colSums(Lx)
e100 <- colSums(Lx[6:nrow(Lx),]) / lx[6,]
e105 <- colSums(Lx[11:nrow(Lx),]) / lx[11,]
e110 <- colSums(Lx[16:nrow(Lx),]) / lx[16,]

library(demogR)
install.packages("http://cran.r-project.org/src/contrib/Archive/demogR/demogR_0.4.2.tar.gz", repos=NULL)
ex <- cdmltw()$ex
exnew <- cbind(ex,e100,e105,e110)
colnames(exnew) <- c(0,1,seq(5,110,by=5))
ex <- exnew
cdmltw110 <- function(){

}
















