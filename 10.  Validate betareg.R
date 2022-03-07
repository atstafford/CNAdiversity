# VALIDATION

beta <- beta_52backrem

# Transpose as input table requires cols=bins
input <- t(testcnBinned[,-c(1:4)]) 

# Convert matrix data to character then factor
input <- data.frame(apply(input, 2, as.character),check.names = FALSE) 
input <- data.frame(lapply(input, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)

# Dummy the relevant cols
input <- input[ ,which(colnames(input) %in% hg19predictors$bin)]
names(input) = paste("bin_", names(input), sep="")
input <- dummy_cols(input, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
rownames(input) <- colnames(testcnBinned[,-c(1:4)])

# Predict
predictedIth_test <- data.frame(predicted = predict(beta, input))

# Enter actual and predicted ITH into a comparative dataframe
predictedIth_test <- data.frame(cbind(actualIth_test, predictedIth_test))
predictedIth_test$type <- 'test'
reg <- lm(actual~predicted, data=predictedIth_test)
summary(reg)

# Add two labels for plot
ggplot(data = predictedIth_test, aes(x = actual, y = predicted, label=sample)) +
  geom_line(aes(group = patient), size=0.2, colour = "#003366") +
  geom_point(fill = "#003366", colour = "#003366", size = 6) +
  geom_text() +
  geom_line(aes(x = actual, y = actual), linetype = "dashed") +
  xlab("Patient CNA diversity") +
  ylab("Patient CNA diversity (predicted)") +
  #xlim(0,0.5) +
  #ylim(0,0.5) +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=24, colour='black'),
        axis.text = element_text(size=24, colour='black'),
        legend.position = "none") +
  geom_smooth(method = 'lm') +
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                        parse = TRUE, label.x = 'right', label.y = 'top', size=8) 

jpeg('tempfig.jpeg', width = 1000, height = 1000)
plot
dev.off()

# BOOTSTRAP on random 27 bins----------------------------------------------------------------------
# Are these 33 bins actually special? try random 33 from candidate and assess model

test33 <- t(car.raw[,-c(1:3)])
test33 <- data.frame(apply(test33, 2, as.character), check.names = FALSE)
test33 <- data.frame(lapply(test33, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
rownames(test33) <- car.info$sampleIDs

test33 <- cbind(ITH=actualIth_train$actual, test33)
names(test33) = paste("bin_", names(test33), sep="")
test33 = test33 %>%
  mutate(across(everything(), as.character))

test33 <- dummy_cols(test33[-1], remove_selected_columns = TRUE, remove_first_dummy = TRUE) 
test33$bin_ITH <- actualIth_train$actual

bins <- unique(as.numeric(gsub("[^\\d]+", "", colnames(test33), perl=TRUE))) #pull bins to select from
bins <- bins[!is.na(bins)]

set.seed(12345)
pR2 <- list()
ll <- list()
tb <- list()
aic2 <- list()
foreach(i = 1:10000) %do% {
  testbins <- sample(bins, 27)
  tb[[i]] <- testbins
  input <- test33[ , which(as.numeric(gsub("[^\\d]+", "", colnames(test33), perl=TRUE)) %in% testbins)]
  input$bin_ITH <- test33$bin_ITH
  try({
    beta <- betareg(bin_ITH ~., data = input)
    pR2[[i]] <- round(beta$pseudo.r.squared,3)
    ll[[i]] <- round(beta$loglik,2)
    aic2[[i]] <- round(AIC(beta),2)
  }, silent = TRUE)
}

length(unlist(aic2))
z <- unlist(aic2)
mean(z)
min(z)

jpeg('tempfig.jpeg', width = 1000, height = 1000)
hist(z, xlim = c(-450,-200))
abline(v = aic, col="blue")
dev.off()


# VALIDATE LOOCV----------------------------------------------------------------
loo <- rbind(t(car.raw[,-c(1:3)]), t(testcnBinned[,-c(1:4)]))
loo <- data.frame(apply(loo, 2, as.character), check.names = FALSE)
loo <- data.frame(lapply(loo, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
rownames(loo) <- c(car.info$sampleIDs, val.info$sampleIDs)
loo <- loo[, which(colnames(loo) %in% hg19predictors$bin)]

loo <- cbind(c(actualIth_train$actual, actualIth_test$actual), loo)
names(loo) = paste("bin_", names(loo), sep="")
loo = loo %>%
  mutate(across(everything(), as.character))

loo <- dummy_cols(loo[-1], remove_selected_columns = TRUE, remove_first_dummy = TRUE) 
loo$bin_ITH <- c(actualIth_train$actual, actualIth_test$actual)
rownames(loo) <- c(car.info$sampleIDs, val.info$sampleIDs)

set.seed(123456)
pR2 <- list()
RMSE <- list()
MAE <- list()
foreach(i = 1:235) %do% {
  input <- sample_n(loo, 234)
  data <- data.frame(predicted=predict(beta_52backrem, input), actual=input$bin_ITH)
  data <- na.omit(data)
  
  r <- lm(actual ~ predicted, data=data)
  pR2[[i]] <- round(summary(r)$r.squared ,3)
  RMSE[[i]] <- round(sqrt(mean((data$actual - data$predicted)^2)), 2)
  MAE[[i]] <- round(mae(data$actual, predict(r)) ,2)
}

mean(unlist(pR2))
mean(unlist(RMSE))
mean(unlist(MAE))

# Model diagnostics-------------------------------------------------------------
plot(beta_back,1) #Residuals vs indices of obs
plot(beta_back,2) #Cook's distance plot
plot(beta_back,3) #Generalized leverage vs predicted values
plot(beta_back,4) #Residuals vs linear predictor
plot(beta_back,5) #Half-normal plot of residuals
plot(beta_back,6) #Predicted vs observed values

# predicted corr with other?----------------------------------------------------
beta <- beta_52backrem
multiReg.in <- multiReg.in_52backrem

actualIth <- data.frame(lapply(data.frame(car.diversity$pic.frac), rep, car.info$sampPerPatient))
predictedIth <- predict(beta, multiReg.in) 
x <- cbind(actualIth, predictedIth)
x$sample <- car.info$sampleIDs
x$patient <- sub('\\..*', '', x$sample)
y <- (car.clonality$patientClo$subclonal/car.clonality$patientClo$CNA)
y <- data.frame(ITH=lapply(data.frame(ITH=y), rep, car.info$sampPerPatient))
x$pcsubclonal <- y$ITH

ggplot(data = x, aes(x = predictedIth, y = pcsubclonal)) +
  geom_line(aes(group = patient), size=0.2, colour = "#003366") +
  geom_point(fill = "#003366", colour = "#003366", size = 6) +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=24, colour='black'),
        axis.text = element_text(size=24, colour='black'),
        legend.position = "none") +
  annotate('text', x = c(0), y = c(0.500,0.450), label = c(paste('pseudoR[adj]^2 ==', psuedoR2), paste('AIC ==', aic)), parse=TRUE, size = 8, hjust = 0, vjust = 1)

# save --------

saveRDS(predictedIth_test, "~/Documents/CNA/Data/predictedIth_test.rds")


