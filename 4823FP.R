library(survival)
library(tidyverse)
library(e1071)
library(GGally)
library(ggfortify)
library(ggplot2)
library(flexsurv)
library(kableExtra)
library(stargazer)
library(cowplot)
library(flexsurvcure)
data <- colon[seq(2,1858,2),]
data <- na.omit(data)
data$sex      <- factor(data$sex)
data$obstruct <- factor(data$obstruct)
data$perfor   <- factor(data$perfor)
data$adhere   <- factor(data$adhere)
data$differ   <- factor(data$differ)
data$extent   <- factor(data$extent)
data$surg     <- factor(data$surg)


nonpara.model <- survfit(Surv(time,status)~rx, data = data, type = "kaplan-meier")
par(mgp=c(2.5,1,0))
plot(nonpara.model, 
     main="Kaplan Meier Survival Probability Estimates",
     ylab =c(expression(paste(hat(S),"(t)"))),
     xlab="Time",
     col=c("black","red","blue"))
legend("topright", 
       legend=c("Obs", "Lev", "Lev+5FU"),
       col=c("black", "red", "blue"),
       lty=1, 
       cex=0.8,
       lwd=2,
       bty="n")

#survdiff on rx
survdiff(Surv(time,status)~rx, data = data)


#AFT MODEL
model.AFT <-survreg(Surv(time,status)~rx + sex + age + obstruct + perfor + adhere +
                      nodes + differ + extent + surg + node4,
                    data = data,
                    dist = "weibull")
summary(model.AFT)

#Weibull Para PH model
model.PH=flexsurvreg(Surv(time,status)~rx + sex + age + obstruct + perfor + adhere +
                       nodes + differ + extent + surg + node4,
                     data = data,
                     dist = "weibullPH")
summary(model.PH)

#Reduced AFT
model.AFT2 <-survreg(Surv(time,status)~ rx + obstruct +  adhere +
                       nodes + extent + surg + node4,
                     data = data,
                     dist = "weibull")
summary(model.AFT2)

#Full PH
model.PH2 <-coxph(Surv(time,status)~rx + sex + age + obstruct + perfor + adhere +
                    nodes + differ + extent + surg + node4,
                  data = data)
summary(model.PH2)

#reduced model
model.PH3 <-coxph(Surv(time,status)~rx + obstruct + adhere +
                    nodes + extent + surg + node4,
                  data = data)

summary(model.PH3)

#Reduced AFT covariarnace matrix
cv.mx <- model.AFT2$var
cv.mx[,12]<-cv.mx[,12]*model.AFT2$scale
cv.mx[12,]<-cv.mx[12,]*model.AFT2$scale 

psi.hat <- (log(600)-(9.555+0.736-0.259-0.3052-0.064*(5)-1.695-0.823))/model.AFT2$scale
var.psi <- 1/model.AFT2$scale^2*c(1,0,1,1,1,5,0,0,1,0,1,psi.hat)%*%cv.mx%*%t(t(c(1,0,1,1,1,5,0,0,1,0,1,psi.hat)))

psi.hat+qnorm(0.975)*sqrt(0.07171199)
psi.hat-qnorm(0.975)*sqrt(0.07171199)


#Model Checking AFT
model.AFT2.res<-exp(-model.AFT2$linear.predictor/model.AFT2$scale)* (Surv(data$time, data$status)[,1])^(1/model.AFT2$scale)
np.fit.res<-survfit(Surv(model.AFT2.res,data$status)~1)
par(mgp=c(2.5,1,0))
par(las=1)
#par(mar=c(5,6,4,2)-1)
plot(np.fit.res$time,np.fit.res$surv,
     type="s", 
     xlab="Residual", 
     ylab=c(expression(paste(hat(S),"(t)"))), 
     main="KM Estimate vs Weibull Estimate of Survival Function")
x <- seq(min(model.AFT2.res), max(model.AFT2.res), length.out = 400) #adding a smooth plot of S(ehat)
y <- exp(- exp(x))
lines(x, y, col = "red", lwd = 2)
legend("topright", 
       c("KM Estimate", "Weibull Estimate"), 
       lty=1,
       col=c("black","red"), 
       bty = "n")

#Mode AFT model checking
np.fit.res<-survfit(Surv(model.AFT2.res,data$status)~1)
plot(log(np.fit.res$time), 
     log(-log(np.fit.res$surv)), 
     type="s",xlab="log(t)", 
     ylab="log(-log(S(t)))", 
     main = "log(-log(S(t))) vs log(t)") 
abline(2,2)


#AFT residual plots (not in analysis)
de.res<-residuals(model.AFT2, type="deviance")

plot(model.AFT2$linear.predictor, 
     de.res, 
     main = "Residuals vs Linear Predictors", 
     xlab = "Linear Predictor",
     ylab="Residual")

plot(Surv(data$time, data$status)[,1], 
     de.res, 
     main = "Residuals vs Time", 
     xlab = "Time (days)",
     ylab="Residual")


#PH model checking
cox.zph(model.PH3)
plot(data$time[data$status==1], 
     residuals(model.PH3, "schoenfeld")[,10], 
     main = "Schoenfeld residuals over time of node4", 
     xlab = "Time (of observed events)",
     ylab = "Schoenfeld Resid")
lines(lowess(data$time[data$status==1], 
             residuals(model.PH3, type="schoenfeld")[,10]))

#attempting to fix time varying node4
data$ind.500 <-factor(as.numeric(data$time>700))

model.PH4 <- coxph(Surv(time,status)~rx +obstruct + adhere +
                     nodes + extent + surg + +node4:ind.500+node4,
                   data = data)

cox.zph(model.PH4)

#back to residual checking in reduced model
#no pattern mean 0.
plot(residuals(model.PH3, "martingale"),
     model.PH3$time,  
     main = "Martingale residuals over time of rx", 
     xlab="Time (of observed events)",
     ylab="Martingale Resid")

plot(data$time, 
     residuals(model.PH3, "score")[,1], 
     main = "Score residuals over time of rx", 
     xlab="Time (of observed events)",
     ylab="Score Resid")
lines(lowess(data$time, residuals(model.PH3, type="score")[,1]))

plot(data$time, residuals(model.PH3, "deviance"), 
     main = "Deviance residuals over time of rx", 
     xlab="Time (of observed events)",
     ylab="Deviance Resid")
lines(lowess(data$time, residuals(model.PH3, type="deviance")))

plot(futime[fustat==1], 
     residuals(cox.fit, type="scaledsch"), 
     xlab="Time",
     ylab="Residual",
     main="scaledsch")
lines(lowess(futime[fustat==1], residuals(cox.fit, type="scaledsch"))) 

#KM vs COXPH
plot(survfit(Surv(time,status)~1,data=data)$time, 
     survfit(Surv(time,status)~1,data=data)$surv,
     type="s",
     xlab="Time", 
     ylab=c(expression(paste(hat(S),"(t)"))), 
     main="KM Estimate vs Cox PH Estimate of Survival Function")

lines(survfit(model.PH3,type="aalen"), col="red")
legend("topright", 
       c("KM Estimate", "Cox PH Estimate"), 
       lty=1,
       col=c("black","red"), 
       bty = "n")


#CURE MODEL
cure_model2 <- flexsurvcure(Surv(time,status)~rx + sex + age + obstruct + perfor + adhere +
                              nodes + differ + extent + surg + node4, data=data, 
                            link="logistic", 
                            dist="weibull", 
                            mixture=T)
#cure_model2

plot(cure_model2, 
     main = "AFT Cure Model vs KM Estimated Model",
     xlab="Time", 
     ylab=c(expression(paste(hat(S),"(t)"))))
legend("topright", 
       legend=c("KM Estimates", "AFT Mixture"),
       col=c("black", "red"), 
       lty=1, 
       cex=0.8,
       lwd=2,
       bty="n")

#Making 2 new data sets, with only 1 type of rx each
data.simple.fu <- data[data$rx=="Lev+5FU",]
data.simple.obs <- data[data$rx=="Obs",]

data.simple.fu$rx <- factor(data.simple.fu$rx)
data.simple.obs$rx <- factor(data.simple.obs$rx)

#vs 1 data set with 2 types of rx, for a simple AFT model
data.simple.reg <- data[data$rx!="Lev",]
data.simple.reg$rx <- factor(data.simple.reg$rx)
model.simple <- survreg(Surv(time,status)~rx,data=data.simple.reg)

km.estimates <- survfit(Surv(time,status)~rx,data=data.simple.reg)

#SIMPLE MIXTURE MODEL
#mixture log likelihood
Mix.Lik.p<-function(t,param,status){
  p=param[1]
  u=param[2]
  b=param[3]
  -sum(log((p/b*exp((log(t)-u)/b)*exp(-exp((log(t)-u)/b)))^(status)*(p*exp(-exp((log(t)-u)/b))+(1-p))^(1-status)))
}

param.est.fu=nlminb(start=c(0.5,5,0.5),
                    Mix.Lik.p,
                    status=data.simple.fu$status,
                    t=data.simple.fu$time)$par

param.est.obs=nlminb(start=c(0.5,5,0.5),
                     Mix.Lik.p,
                     status=data.simple.obs$status,
                     t=data.simple.obs$time)$par

#survival probabilities using MLE
st=param.est.fu[1]*exp(-exp((log(data.simple.fu$time)-param.est.fu[2])/param.est.fu[3]))+(1-param.est.fu[1])
st2=param.est.obs[1]*exp(-exp((log(data.simple.obs$time)-param.est.obs[2])/param.est.obs[3]))+(1-param.est.obs[1])

df.test <- data.frame(time=data.simple.fu$time, surv=st)
df.test <- df.test[order(df.test$time),,drop=F]

df.test2 <- data.frame(time=data.simple.obs$time, surv=st2)
df.test2 <- df.test2[order(df.test2$time),,drop=F]

par(las=1)
plot(km.estimates, main="Mixture Model vs. Traditional Weibull AFT",
     ylab =c(expression(paste(hat(S),"(t)"))),xlab="Time")
lines(predict(model.simple, 
              newdata=list(rx="Obs"),
              type="quantile",
              p=seq(0,1,by=.01)),
      seq(1,0,by=-.01),
      col="darkgreen")

lines(predict(model.simple,
              newdata=list(rx="Lev+5FU"),
              type="quantile",
              p=seq(0,1,by=.01)),
      seq(1,0,by=-.01),
      col="darkgreen")
lines(df.test$time,df.test$surv,type="l",col="red")
lines(df.test2$time,df.test2$surv,type="l",col="red")
legend("topright", legend=c("KM Estimates", "AFT", "Mixture"),
       col=c("black", "darkgreen", "red"), lty=1, cex=0.8,lwd=2,bty="n")



######ADDITIONAL CODE TO PLOT KM ESTIMATES FOR EACH SINGLE COVARIATE MODEL
plot1<-autoplot(test.model11,
                conf.int = FALSE,
                censor.size = 2,
                surv.size = 1.2, 
                ylab = c(expression(paste(hat(S),"(t)"))),
                xlab = "Time (days)",
                main = "KM",
                legend=F)+
  labs(color="Trt") +
  theme(legend.key.size = unit(0.5, "cm"))+
  theme(axis.title.y = element_text(angle=0,vjust=0.5))
plot3<-autoplot(test.model3,
                conf.int = FALSE,
                censor.size = 2,
                surv.size = 1.2, 
                ylab = "",
                xlab = "Time (days)",
                main = "KM",
                legend=F)+
  labs(color="Trt") +
  theme(legend.key.size = unit(0.5, "cm"))+
  theme(axis.title.y = element_text(angle=0,vjust=0.5))

plot4<-autoplot(test.model4,
                conf.int = FALSE,
                censor.size = 2,
                surv.size = 1.2, 
                ylab = "",
                xlab = "Time (days)",
                main = "KM",
                legend=F)+
  labs(color="Trt") +
  theme(legend.key.size = unit(0.5, "cm"))+
  theme(axis.title.y = element_text(angle=0,vjust=0.5))
plot5<-autoplot(test.model5,
                conf.int = FALSE,
                censor.size = 2,
                surv.size = 1.2, 
                ylab = c(expression(paste(hat(S),"(t)"))),
                xlab = "Time (days)",
                main = "KM",
                legend=F)+
  labs(color="Trt") +
  theme(legend.key.size = unit(0.5, "cm"))+
  theme(axis.title.y = element_text(angle=0,vjust=0.5))
plot7<-autoplot(test.model7,
                conf.int = FALSE,
                censor.size = 2,
                surv.size = 1.2, 
                ylab = "",
                xlab = "Time (days)",
                main = "KM",
                legend=F)+
  labs(color="Trt") +
  theme(legend.key.size = unit(0.5, "cm"))+
  theme(axis.title.y = element_text(angle=0,vjust=0.5))
plot8<-autoplot(test.model8,
                conf.int = FALSE,
                censor.size = 2,
                surv.size = 1.2, 
                ylab = "",
                xlab = "Time (days)",
                main = "KM",
                legend=0)+
  labs(color="Trt") +
  theme(legend.key.size = unit(0.5, "cm"))+
  theme(axis.title.y = element_text(angle=0,vjust=0.5))
plot9<-autoplot(test.model9,
                conf.int = FALSE,
                censor.size = 2,
                surv.size = 1.2, 
                ylab = c(expression(paste(hat(S),"(t)"))),
                xlab = "Time (days)",
                main = "KM",
                legend=F)+
  labs(color="Trt") +
  theme(legend.key.size = unit(0.5, "cm"))+
  theme(axis.title.y = element_text(angle=0,vjust=0.5))
plot10<-autoplot(test.model10,
                 conf.int = FALSE,
                 censor.size = 2,
                 surv.size = 1.2, 
                 ylab = c(expression(paste(hat(S),"(t)"))),
                 xlab = "Time (days)",
                 main = "KM",
                 legend=F)+
  labs(color="Trt") +
  theme(legend.key.size = unit(0.5, "cm"))+
  theme(axis.title.y = element_text(angle=0,vjust=0.5))

plot_grid(plot1,plot3,plot4,plot5,plot7,plot8,plot9,plot10, ncol=3,nrow=3)

