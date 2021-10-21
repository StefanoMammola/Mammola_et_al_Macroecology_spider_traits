library(bbgdm)


model <-readRDS("Models/Null_BBGDM.rds")
m1 <- model$Trial_1$Btotal$model
resids <- diagnostics(m1$model)
par(mfrow=c(2,2))
plot(resids)
plot.bbgdm(m1)


model$Trial_2$Btotal$beta_metric
