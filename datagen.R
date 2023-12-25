rm(list=ls());gc();dev.off()
library(sna)
set.seed(4395834L)

n <- 100

w1<-rgraph(n)               #Draw the AR matrix
w2<-w1               #Draw the MA matrix
x<-matrix(rnorm(n*5),n,5) #Draw some covariates
r1<-0.2                       #Set the model parameters
r2<-0.1
sigma<-0.1
beta<-rnorm(5)
#Assemble y from its components:
nu<-rnorm(n,0,sigma)          #Draw the disturbances
e<-qr.solve(diag(n)-r2*w2,nu) #Draw the effective errors
y<-qr.solve(diag(n)-r1*w1,x%*%beta+e)  #Compute y

# write info
write.csv(file="./data/w1.csv", x=w1, row.names=FALSE)
write.csv(file="./data/w2.csv", x=w2, row.names=FALSE)
write.csv(file="./data/x.csv", x=x, row.names=FALSE)
write.csv(file="./data/beta.csv", x=beta, row.names=FALSE)
write.csv(file="./data/nu.csv", x=nu, row.names=FALSE)
write.csv(file="./data/e.csv", x=e, row.names=FALSE)
write.csv(file="./data/y.csv", x=y, row.names=FALSE)

#Now, fit the autocorrelation model:
fit<-lnam(y = y,x = x,W1 = w1,W2 = w2)

# spatial fit
library(spatialreg)
library(spdep)

w1_sp <- spdep::mat2listw(x = w1)
w2_sp <- w1_sp

# method 1: MLE
fit_sp <- spatialreg::sacsarlm(
  formula = y ~ x - 1,
  listw = w1_sp,
  listw2 = w2_sp,
  Durbin = FALSE
)

# method 2: GMM
fit_gmm <- spatialreg::gstsls(
  formula = y ~ x - 1,
  listw = w1_sp,
  listw2 = w2_sp
)

# # spsur
# library(spsur)
# 
# fit_spsur <- spsur::spsurml(
#   formula = y ~ x - 1,
#   data = data.frame(y=y,x=y),
#   listw = w1,
#   type = "sarar"
# )

# compare
results <- c(beta,r1,r2,sigma)
results <- rbind(
  results,
  c(as.vector(fit$beta),fit$rho1,fit$rho2,fit$sigma),
  c(as.vector(fit_sp$coefficients),fit_sp$rho,fit_sp$lambda,sqrt(fit_sp$s2)),
  c(fit_gmm$coefficients[2:6],fit_gmm$coefficients[1],fit_gmm$lambda,sqrt(fit_gmm$s2))
)
results <- as.data.frame(results)
rownames(results) <- c("data","sna","sp MLE","sp GMM")
colnames(results) <- c(paste0("beta",1:5),"rho1","rho2","sigma")
results
