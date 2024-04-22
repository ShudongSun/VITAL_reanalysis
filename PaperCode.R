### import data ----
vital_data <- read.csv(file = '..//data//vital_data_2019.csv')
print(dim(vital_data)) # check dimension of data
print(colSums(is.na(vital_data))) # check missing values
features <- c('Subject_ID','vitdactive','fishoilactive','sex','ageyr','raceth',
              'bmi','currsmk','htnmed','cholmed','diabetes','diabmed','parhxmi',
              'fish1_5wk','Aspirin','statins','N_risk_factors','VitDIntake',
              'VitD20','VitD31')
outcomes <- c('malca','malcayrs','brca','brcayrs','prca','prcayrs','colca',
              'colcayrs','majorcvd','majryears','imporcvd','impryears','totmi',
              'miyears','totst','totchd','chdyears','ischemic','styears',
              'hemorrhg','ptca','ptcayears','cabg','cabgyears','confdeath',
              'cadeath','cvdeath','chdeath','mideath','stdeath','randyrs')


composite_end_point <- c('totmi','totst','cvdeath')
vital_data_features = vital_data[,features]
y_MI = as.factor(vital_data$totmi)
y_comp_vars = vital_data[,composite_end_point]
y_comp = as.factor(ifelse(y_comp_vars[,1]|y_comp_vars[,2]|y_comp_vars[,3],1,0))

y_MajorCVD = as.factor(vital_data$majorcvd)

y_stroke = as.factor(vital_data$totst)

y_cvdeath = as.factor(vital_data$cvdeath)

for(i in 1:length(vital_data_features)){
  if(!(i %in% c(1,5,7))){
    vital_data_features[,i] = as.factor(vital_data_features[,i])
  }
}

### deal with missing ----
library(Hmisc)
library(DescTools)

# FOR MI:
set.seed(2022)

data_MI = vital_data_features[,-1]
data_MI$MI = y_MI

# keep 1,2 in raceth only
data_MI = data_MI[which(data_MI$raceth %in% c(1,2)),]
data_MI$raceth = factor(data_MI$raceth)

data_m3_MI = subset(data_MI, select = -c(VitD20, VitD31) )

# delete the patients with missing rows for all other variables:
data_m3_MI = data_m3_MI[which(rowSums(is.na(data_m3_MI))==0),] 
data = data_m3_MI

### matching data
library("MatchIt")

data = data_m3_MI

data$raceth = ifelse(data$raceth==2,1,0)
names(data)[which(names(data)=="raceth")] = "raceth_black"

##  Propensity Score Matching ----

# 1:1 Nearest Neighbot matching w/o replacement
m.out1 <- matchit(raceth_black ~ vitdactive + fishoilactive + sex + ageyr + bmi + currsmk + 
                    htnmed + cholmed + diabetes + diabmed + parhxmi + fish1_5wk + Aspirin +
                    statins + N_risk_factors + VitDIntake, data = data,
                  method = "nearest", distance = "glm")

summary(m.out1)

# plot(m.out1)
# plot(m.out1, type = "hist", interactive = FALSE)
# plot(summary(m.out1))
# plot(m.out1, type = "qq", interactive = FALSE)

# 1:1 optimal matching w/o replacement
m.out2 <- matchit(raceth_black ~ vitdactive + fishoilactive + sex + ageyr + bmi + currsmk + 
                    htnmed + cholmed + diabetes + diabmed + parhxmi + fish1_5wk + Aspirin +
                    statins + N_risk_factors + VitDIntake, data = data,
                  method = "optimal", distance = "glm")

summary(m.out2)

# plot(m.out2, type = "hist", interactive = F)
# plot(summary(m.out2))


# optimal matching w/o replacement
m.out3 <- matchit(raceth_black ~ vitdactive + fishoilactive + sex + ageyr + bmi + currsmk + 
                    htnmed + cholmed + diabetes + diabmed + parhxmi + fish1_5wk + Aspirin +
                    statins + N_risk_factors + VitDIntake, data = data,
                  method = "full", distance = "glm")

summary(m.out3)

# plot(m.out3, type = "hist", interactive = F)
# plot(summary(m.out3))

# genetic matching
m.out4 <- matchit(raceth_black ~ vitdactive + fishoilactive + sex + ageyr + bmi + currsmk + 
                    htnmed + cholmed + diabetes + diabmed + parhxmi + fish1_5wk + Aspirin +
                    statins + N_risk_factors + VitDIntake, data = data,
                  # distance = "mahalanobis",
                  method = "genetic", pop.size = 1000)
summary(m.out4)

# plot(m.out4, type = "hist", interactive = F)
# plot(summary(m.out4))


# save(m.out1,m.out2,m.out3, file="./m.out123_majorVD.Rdata")
# save(m.out4, file="./m.out4.Rdata")


m.data1 <- match.data(m.out1)
m.data2 <- match.data(m.out2)
m.data3 <- match.data(m.out3)
m.data4 <- match.data(m.out4)




### Choose optimal pair matching ----

data = data_m3_MI

data$raceth = ifelse(data$raceth==2,1,0)
names(data)[which(names(data)=="raceth")] = "raceth_black"

# load(".//m.out123.Rdata")

m.data2 <- match.data(m.out2)
data_om_MI = m.data2[,1:(length(m.data2)-3)]

data_om_MI$raceth_black = factor(data_om_MI$raceth_black)

data_om_MI_black = subset(data_om_MI[which(data_om_MI$raceth_black==1),], select = -c(raceth_black))
data_om_MI_white = subset(data_om_MI[which(data_om_MI$raceth_black==0),], select = -c(raceth_black))


### Decision Tree After matching ----

library(rpart)
library(rpart.plot)
library(rattle)
require(caret)

train_data = data_om_MI_black
# train_data = data_om_MI_white

train_data_x = train_data[1:(length(train_data)-1)]

names(train_data)[length(names(train_data))] = "y"

train_y = train_data$y
table(train_y)
print(table(train_y)[1]/table(train_y)[2])

weight_num = (sum(train_data$y==0))/(sum(train_data$y==1))
# print(weight_num)

weights = ifelse(train_data$y==1,weight_num/(weight_num+1),1/(weight_num+1))
# weights = ifelse(train_data$y==1,weight_num,1) * om_weights  # for optimal full pair matched  

set.seed(2022)

tree <- rpart(y~., data=train_data, method='class',weights = weights,
              control = rpart.control(xval = 5, minsplit = 20, cp = 0, usesurrogate = 2))



# index_selected = min(which( tree$cptable[,"xerror"] < (min(tree$cptable[,"xerror"])+tree$cptable[which.min(tree$cptable[,"xerror"]),"xstd"]) ))
index_selected = as.numeric(which( tree$cptable[,"xerror"] == min(tree$cptable[,"xerror"])))

ptree<- prune(tree, cp= tree$cptable[min(index_selected),"CP"])

fancyRpartPlot(ptree, main="Pruned Classification Tree",sub="")

summary(ptree)


### Lasso model to select variables and refit the logistic regression to reduce bias ----

train_data = data_om_MI


names(train_data)[length(names(train_data))] = "y"

train_data_dummy_raceth = train_data

f <- as.formula(y ~ fishoilactive*.)
train_data_dummy_raceth_fishinteract_x <- data.frame(model.matrix(f, train_data_dummy_raceth)[, -1])

train_y = train_data_dummy_raceth$y

library(glmnet)
require(caret)
library(selectiveInference)

weight_num = (sum(train_data$y==0))/(sum(train_data$y==1))

weights = ifelse(train_data$y==1,weight_num/(weight_num+1),1/(weight_num+1))
# weights = ifelse(train_data$y==1,weight_num,1) * om_weights # for optimal full pair matched

train_x = data.matrix(train_data_dummy_raceth_fishinteract_x)

set.seed(2023)
glmnet_cv <- cv.glmnet(x=train_x, y=train_y, family = "binomial",
                       type.measure = "class",nfolds=5, alpha=1,weights = weights)
plot(glmnet_cv)

lambda=glmnet_cv$lambda.1se

glmnet_fit = glmnet(x=train_x, y=train_y,
                    family="binomial",alpha=1, weights = weights)

cof = coef(glmnet_fit, s= lambda)
coef_pca2_selected = cof@Dimnames[[1]][cof@i+1][-1]
coef_selected = data.frame(selected_var=cof@Dimnames[[1]][cof@i+1], coefficient=cof@x)

train_data_dummy_raceth_fishinteract_x_selected = train_data_dummy_raceth_fishinteract_x[,coef_pca2_selected]

train_data_selected = data.frame(train_data_dummy_raceth_fishinteract_x_selected,train_y)
names(train_data_selected) = c(coef_pca2_selected,"y")

glm_refit = glm(y~., data=train_data_selected, family=binomial, weights = weights)
summary(glm_refit)


### Lasso: Parametric Bootstrap to calculate standard errors and Confidence Interval ----

train_data = data_om_MI

names(train_data)[length(names(train_data))] = "y"

train_data_dummy_raceth = train_data

f <- as.formula(y ~ fishoilactive*.)
train_data_dummy_raceth_fishinteract_x <- data.frame(model.matrix(f, train_data_dummy_raceth)[, -1])

train_y = train_data_dummy_raceth$y

library(glmnet)
require(caret)
library(selectiveInference)
library(HDCI)

weight_num = (sum(train_data$y==0))/(sum(train_data$y==1))
weights = ifelse(train_data$y==1,weight_num/(weight_num+1),1/(weight_num+1))
train_x = data.matrix(train_data_dummy_raceth_fishinteract_x)

set.seed(2023)

x <- as.matrix(train_x)
y <- train_y
n <- dim(train_x)[1]
p <- dim(train_x)[2]

selectset <- rep(0, p)
Beta <- rep(0, p)

cvfit <- cv.glmnet(x, y, family = "binomial", standardize = TRUE,
                   type.measure = "class", nfolds=5, alpha=1,weights = weights)

lambda.opt <- cvfit$lambda.1se

globalfit <- glmnet(x, y,family = "binomial", standardize = TRUE, type.measure = "class",
                    alpha=1, weights = weights)

fitlasso <- coef(globalfit, s = lambda.opt)

coef_selected_vars_global = fitlasso@Dimnames[[1]][fitlasso@i+1][-1]

coef_selected_global = data.frame(selected_var=fitlasso@Dimnames[[1]][fitlasso@i+1], coefficient=fitlasso@x)

x_selected = x[,coef_selected_vars_global]

data_selected = data.frame(x_selected,y)
names(data_selected) = c(coef_selected_vars_global,"y")

glm_refit_global = glm(y~., data=data_selected, family=binomial, weights = weights)
summary(glm_refit_global)
fitglm = coef(glm_refit_global)

beta = fitlasso@x
x_lasso = x[,c(1:p) %in% fitlasso@i]
p_predict = exp(cbind2(1,x_lasso)%*%fitglm)/(1+exp(cbind2(1,x_lasso)%*%fitglm))
p0 = p_predict/(weight_num-(weight_num-1)*p_predict)

B=2000

Beta.boot <- matrix(0, nrow = B, ncol = p+1)
colnames(Beta.boot) = row.names(fitlasso)

out <- list()
i=1
while(i <=B) {
  resam <- sample(1:n, n, replace = TRUE)
  
  rx <- x[resam,]
  
  rx_lasso <- x_lasso[resam,]
  rp_predict = exp(cbind2(1,rx_lasso)%*%fitglm)/(1+exp(cbind2(1,rx_lasso)%*%fitglm))
  rp0 = rp_predict/(weight_num-(weight_num-1)*rp_predict)
  ry = as.factor(rbinom(n,1,rp0))
  
  rweight_num = summary(ry)[1]/summary(ry)[2]
  rweights = ifelse(ry==1,weight_num/(weight_num+1),1/(weight_num+1))
  
  globalfit0 <- glmnet(rx, ry, family = "binomial", standardize = TRUE, type.measure = "class",
                       alpha=1, weights = rweights)
  cof <- coef(globalfit0, s = lambda.opt)
  
  coef_selected_vars = cof@Dimnames[[1]][cof@i+1][-1]
  
  coef_selected = data.frame(selected_var=cof@Dimnames[[1]][cof@i+1], coefficient=cof@x)
  
  rx_selected = rx[,coef_selected_vars]
  
  rdata_selected = data.frame(rx_selected,ry)
  names(rdata_selected) = c(coef_selected_vars,"y")
  
  glm_refit = glm(y~., data=rdata_selected, family=binomial, weights = rweights)
  summary(glm_refit)
  
  print(i)
  if(is.na(coef(glm_refit)["fishoilactive1.raceth_black1"]) |
     (!is.na(coef(glm_refit)["fishoilactive1.raceth_black1"]) & (abs(coef(glm_refit)["fishoilactive1.raceth_black1"]-fitglm["fishoilactive1.raceth_black1"]) <= 10))){
    out[[i]] <- coef(glm_refit)
    i=i+1
  }

}

df_repeat = data.frame(out[[1]])
df_repeat$var = row.names(df_repeat)
colnames(df_repeat)[1] = 1
for (i in 2:B) {
  df_temp = data.frame(out[[i]])
  df_temp$var = row.names(df_temp)
  colnames(df_temp)[1] = i
  
  df_repeat = merge(df_repeat, df_temp, by="var", all=TRUE)

}
df_repeat[is.na(df_repeat)] <- 0

Beta.boot = t(df_repeat[,-1])
colnames(Beta.boot) = df_repeat[,1]

var_names = c("(Intercept)",colnames(x))
Beta.boot = Beta.boot[,setdiff(var_names,setdiff(var_names,colnames(Beta.boot)))]

cov_matrix = matrix(data=0, nrow = dim(Beta.boot)[2], ncol = dim(Beta.boot)[2])

rownames(cov_matrix) <- colnames(Beta.boot)
colnames(cov_matrix) <- colnames(Beta.boot)

standard_errors = numeric(dim(Beta.boot)[2])

for(i in 1:dim(Beta.boot)[2]){
  for(j in 1:dim(Beta.boot)[2]){
    cov_matrix[i,j] = cov(Beta.boot[,i],Beta.boot[,j])
    if(i==j){
      standard_errors[i] = sqrt(cov_matrix[i,j])
    }
  }
}

options(scipen=10)
vars_sds = cbind(var_names, standard_errors)
coef_center = cbind(names(fitglm),fitglm)
colnames(coef_center)[1] = "var_names"
result = merge(coef_center, vars_sds, by="var_names", all=TRUE, sort=FALSE)
result[,2][is.na(result[,2])]=0
row.names(result) = result[,1]
result = result[var_names,c(2,3)]
result = result[result$fitglm!="0",]
keep_digit = 4
result[,1] = round(as.numeric(result[,1]),keep_digit)
result[,2] = round(as.numeric(result[,2]),keep_digit)
result["p-value"] = round(2*pnorm(-abs(as.numeric(result[,1])/as.numeric(result[,2]))),keep_digit)
ORs = round(exp(as.numeric(result[,1])),keep_digit)
CI_lower = round(exp(as.numeric(result[,1])-1.96*as.numeric(result[,2])),keep_digit)
CI_higher = round(exp(as.numeric(result[,1])+1.96*as.numeric(result[,2])),keep_digit)
result["OR (95% CI)"] = paste(ORs, " (",CI_lower,", ", CI_higher,")",sep="")

# write.csv(result, file="parametric_MI.csv")
print(result)


### Lasso: Non-parametric Bootstrap to calculate standard errors and Confidence Interval ----

train_data = data_om_MI

names(train_data)[length(names(train_data))] = "y"

train_data_dummy_raceth = train_data

f <- as.formula(y ~ fishoilactive*.)
train_data_dummy_raceth_fishinteract_x <- data.frame(model.matrix(f, train_data_dummy_raceth)[, -1])

train_y = train_data_dummy_raceth$y


library(glmnet)
require(caret)

weight_num = (sum(train_data$y==0))/(sum(train_data$y==1))
print(weight_num)

weights = ifelse(train_data$y==1,weight_num/(weight_num+1),1/(weight_num+1))

train_x = data.matrix(train_data_dummy_raceth_fishinteract_x)

set.seed(2023)

x <- as.matrix(train_x)
y <- train_y
# y <- train_y
n <- dim(train_x)[1]
p <- dim(train_x)[2]

selectset <- rep(0, p)
Beta <- rep(0, p)

cvfit <- cv.glmnet(x, y, family = "binomial", standardize = TRUE,
                   type.measure = "class", nfolds=5, alpha=1, weights = weights)

lambda.opt <- cvfit$lambda.1se

globalfit <- glmnet(x, y,family = "binomial", standardize = TRUE, type.measure = "class",
                    alpha=1, weights = weights)

fitlasso <- coef(globalfit, s = lambda.opt)

coef_selected_vars_global = fitlasso@Dimnames[[1]][fitlasso@i+1][-1]

coef_selected_global = data.frame(selected_var=fitlasso@Dimnames[[1]][fitlasso@i+1], coefficient=fitlasso@x)

x_selected = x[,coef_selected_vars_global]

data_selected = data.frame(x_selected,y)
names(data_selected) = c(coef_selected_vars_global,"y")

glm_refit_global = glm(y~., data=data_selected, family=binomial, weights = weights)
summary(glm_refit_global)
fitglm = coef(glm_refit_global)

B=2000

Beta.boot <- matrix(0, nrow = B, ncol = p+1)
colnames(Beta.boot) = row.names(fitlasso)

out <- list()
i=1
while(i <= B) {
  # resam <- sample(1:n, n, replace = TRUE)
  resam <- c(sample(which(y==0),sum(y==0),replace = TRUE),sample(which(y==1),sum(y==1),replace = TRUE))
  
  rx <- x[resam,]
  ry <- y[resam]
  rweights <- weights[resam]
  
  globalfit0 <- glmnet(rx, ry, family = "binomial", standardize = TRUE,
                       alpha=1, weights = rweights)
  cof <- coef(globalfit0, s = lambda.opt)
  
  coef_selected_vars = cof@Dimnames[[1]][cof@i+1][-1]
  
  rx_selected = rx[,coef_selected_vars]
  
  rdata_selected = data.frame(rx_selected,ry)
  names(rdata_selected) = c(coef_selected_vars,"y")
  
  glm_refit = glm(y~., data=rdata_selected, family=binomial, weights = rweights)

  summary(glm_refit)
  
  print(i)
  if(is.na(coef(glm_refit)["fishoilactive1.raceth_black1"]) |
     (!is.na(coef(glm_refit)["fishoilactive1.raceth_black1"]) & (abs(coef(glm_refit)["fishoilactive1.raceth_black1"]-fitglm["fishoilactive1.raceth_black1"]) <= 10))){
    out[[i]] <- coef(glm_refit)
    i=i+1
  }
  
}

df_repeat = data.frame(out[[1]])
df_repeat$var = row.names(df_repeat)
colnames(df_repeat)[1] = 1
for (i in 2:B) {
  df_temp = data.frame(out[[i]])
  df_temp$var = row.names(df_temp)
  colnames(df_temp)[1] = i
  
  df_repeat = merge(df_repeat, df_temp, by="var", all=TRUE)
  
}
df_repeat[is.na(df_repeat)] <- 0

Beta.boot = t(df_repeat[,-1])
colnames(Beta.boot) = df_repeat[,1]

var_names = c("(Intercept)",colnames(x))
Beta.boot = Beta.boot[,setdiff(var_names,setdiff(var_names,colnames(Beta.boot)))]

cov_matrix = matrix(data=0, nrow = dim(Beta.boot)[2], ncol = dim(Beta.boot)[2])
rownames(cov_matrix) <- row.names(fitlasso)
colnames(cov_matrix) <- row.names(fitlasso)

standard_errors = numeric(dim(Beta.boot)[2])

for(i in 1:dim(Beta.boot)[2]){
  for(j in 1:dim(Beta.boot)[2]){
    cov_matrix[i,j] = cov(Beta.boot[,i],Beta.boot[,j])
    if(i==j){
      standard_errors[i] = sqrt(cov_matrix[i,j])
    }
  }
}

options(scipen=10)
vars_sds = cbind(var_names, standard_errors)
coef_center = cbind(names(fitglm),fitglm)
colnames(coef_center)[1] = "var_names"
result = merge(coef_center, vars_sds, by="var_names", all=TRUE, sort=FALSE)
result[,2][is.na(result[,2])]=0
row.names(result) = result[,1]
result = result[var_names,c(2,3)]

result = result[result$fitglm!="0",]
keep_digit = 4
result[,1] = round(as.numeric(result[,1]),keep_digit)
result[,2] = round(as.numeric(result[,2]),keep_digit)
result["p-value"] = round(2*pnorm(-abs(as.numeric(result[,1])/as.numeric(result[,2]))),keep_digit)
ORs = round(exp(as.numeric(result[,1])),keep_digit)
CI_lower = round(exp(as.numeric(result[,1])-1.96*as.numeric(result[,2])),keep_digit)
CI_higher = round(exp(as.numeric(result[,1])+1.96*as.numeric(result[,2])),keep_digit)
result["OR (95% CI)"] = paste(ORs, " (",CI_lower,", ", CI_higher,")",sep="")

print(result)
