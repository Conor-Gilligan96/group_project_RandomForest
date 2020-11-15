#Install and load required packages 
install.packages("randomForest")   
install.packages("rsample") 
install.packages("ranger") 
require(tidyverse) 
require(randomForest)
library(rsample)      
library(ranger)       
       
#import data set 
Problem1_expressionMatrix <- read.delim(
  "~/Desktop/PhD_CPG/Satija/Group/teamproject1-team-1-master/Problem1_expressionMatrix.txt") 

prob1.df <- data.frame(Problem1_expressionMatrix)

##rearrange data so each gene is a column

#first, transpose 

tranposed_prob1 <- t(prob1.df) 

tranposed_prob1.df <- data.frame(tranposed_prob1)

reshaped_prob1 <- rownames_to_column(tranposed_prob1.df, var = "celltype")  

#arrange data 

reshaped_prob1.v2 <- reshaped_prob1[order(reshaped_prob1$celltype),] 

#Convertcell type into a factor (6 in total) 
reshaped_prob1.v3 <- reshaped_prob1.v2

reshaped_prob1.v3 %>% mutate(classification = str_extract(celltype, ".*_")) -> reshaped_prob1.v4

reshaped_prob1.v4[50,1241] <- "Hyperdip_50_"

reshaped_prob1.v4$classification <- as.factor(reshaped_prob1.v4$classification) 

#Bring classification column to the front of the data frame 

final.prob1.df<- reshaped_prob1.v4[,c(ncol(reshaped_prob1.v4),1:(ncol(reshaped_prob1.v4)-1))]  

##Running the Random Forest using randomForest package in R 

#create consolidated data frame 

random.forest.df <- final.prob1.df[,-(2)]  

#Split data into test and training- 20 percent and 80 percent respectivley. Set seed for reproducibility 
set.seed(101)
rf.split.data <- initial_split(random.forest.df, prop = .8)
train.df <- training(rf.split.data)
test.df  <- testing(rf.split.data) 

#First, I'll implement random forest on full data set. In this example, I am making 1000 trees, 
#and at each internal node, the square root of the number of variables (genes) is used.
#In this case, 35 genes (1239^0.5) are considered at each node.  
#set seed to make outputs reproducible 

set.seed(1)
full.RF <- randomForest(formula = classification ~ ., data = random.forest.df , ntree = 1000)
plot(full.RF) #It looks like the Error stabilizes after ~750 trees are grown.

print(full.RF) #OOB estimate of error rate: 3.66%

#Now, use the subset training data set 

set.seed(1)

subset.RF <- randomForest(formula = classification ~ ., data = train.df, ntree = 1000)
plot(subset.RF) 

#check out the results 
print(subset.RF) # Our subset RF model has given us an Out of Bag (OOB) Error rate of 3.03% 


##Let's use the subset model to find our best genes. Even though 
#Random Forest does not overfit data, it is reccomended to have test and training data

#We can first work on optimizing ntree as an example 
#Later, we will tune many hyperparameters simultaneously. 

#Subset data
subset.RF.errorrate.df <- subset.RF$err.rate 

head(subset.RF.errorrate.df)

#create plot showing OOB in relation to number of trees produced 

outofbagerror <- data.frame(tree_number=rep(1:nrow(subset.RF.errorrate.df), times=7),
category=rep(c("OOB", "BCR_ABL_", "E2A_PBX1_", "Hyperdip_50_", "MLL_", "T_ALL_", "TEL_AML1_"), each=nrow(subset.RF.errorrate.df)),
Error=c(subset.RF.errorrate.df[,"OOB"], subset.RF.errorrate.df[,"BCR_ABL_"],subset.RF.errorrate.df[,"E2A_PBX1_"],subset.RF.errorrate.df[,"Hyperdip_50_"],
subset.RF.errorrate.df[,"MLL_"], subset.RF.errorrate.df[,"T_ALL_"],subset.RF.errorrate.df[,"TEL_AML1_"])) 

subsetplot<- ggplot(data=outofbagerror, aes(x=tree_number, y=Error)) +
  geom_line(aes(color=category))  

print(subsetplot)

#find tree number w/ lowest oob error 
require(dplyr)
arrange(outofbagerror, Error) -> lowerror  
filter(lowerror, category == "OOB") -> lowOOB 
head(lowOOB, n =1) #69 trees yielded the lowest OOB

#adjust number of trees to optimal number (69)  

set.seed(1)

initial.RF.v2 <- randomForest(formula = classification ~ ., data = train.df, ntree = 69)
plot(initial.RF.v2) 

#check out the results 
print(initial.RF.v2) # we now have essentially the same 3.03% OOB error rate  

#Only adjusting the ntree hyperparemter did not do much, 
#so it will be better to focus on remaining hyperparameters  

##Adjusting mtry (number of predictors (genes) sampled for spliting at each node) 
#Default value is sqrt(p)= 35 genes (1239^0.5). 

#quick method using tuneRF found in randomForest packages

gene_predictors <- setdiff(names(train.df), "classification")

set.seed(1)

tune.random.forest.model <- tuneRF(
  x          = train.df[gene_predictors],
  y          = train.df$classification,
  ntreeTry   = 1000,
  mtryStart  = 35,
  stepFactor = 2,
  improve    = 0.01,
  trace      = TRUE      
)  

#Here, I specified tuneRF to start at the mtry default value of 35. Then, adjusting up and 
#down by a step factor of 5 of 2, tuneRF calculates the OOB error for each RF of 1000 trees, 
#until it stops improving by the improvment threshold I set (0.01) 

#It appears that 35 looks pretty good!

##Slow method using for loop cycling through 
#possible values of (it takes a long time, so I commented it out), but it seems that 34 is a good choice.
# set.seed(1)
# out.of.bag.forloop <- vector(length=36)
# for(i in 1:36) {
#out.of.bag.model <- randomForest(classification ~ ., data=train.df, mtry=i, ntree=1000)
#   out.of.bag.forloop[i] <- out.of.bag.model$err.rate[nrow(out.of.bag.model$err.rate),1]
# }
# out.of.bag.forloop
# #identify optimum value of mtry (in range of 1:36)
# which(out.of.bag.forloop == min(out.of.bag.forloop))  



##Next, we can make this process more efficient by tuning many hyperparameters simultaneously 
#Use expand.grid to cycle through different combinations of hyperparameters.  

hyperparameters <- expand.grid(mtry = seq(20, 66, by = 2),min.node.size  = seq(2, 20, by = 2),
sample.fraction = c(.5, .6, .7, .8, .9, 1),
OOB   = 0
) 

#Specifically, we can choose optimal hyperparameters, which include:  
#mtry = Number of variables (genes) randomly sampled as candidates at each split. 
#node_size = Minimum size of terminal nodes 
#sample_size= Size (or proportion) of sample to draw from the data.

nrow(hyperparameters) #This results in 1440 models (where ntree=1000) we can choose from.

#Ranger is a fast implementation of random forests (Breiman 2001) 

for(i in 1:nrow(hyperparameters)) {
  
  # train model
  ranger.model <- ranger(
    formula         = classification ~ ., 
    data            = train.df, 
    num.trees       = 1000,
    mtry            = hyperparameters$mtry[i],
    min.node.size   = hyperparameters$min.node.size[i],
    sample.fraction = hyperparameters$sample.fraction[i],
    seed            = 1
  )
  
  # add OOB error to grid
  hyperparameters$OOB[i] <- ranger.model$prediction.error
}

hyperparameters %>% arrange(OOB) %>% head(1) -> final_hyperparameters.df

print(final_hyperparameters.df) 

#implementing hyperparameters (ntree (num.trees)= 1000, mtry= 58, min.node.size = 8, sample.fraction = 0.7) 


final.ranger.model <- ranger(
    formula         = classification ~ ., 
    data            = train.df, 
    num.trees       = 1000,
    mtry            = 58,
    min.node.size   = 8,
    sample.fraction = .7,
    importance      = 'impurity')

  

variableimportance.df<- data.frame(final.ranger.model$variable.importance)  
variableimportancev2.df<- rownames_to_column(variableimportance.df, var = "gene") 
colnames(variableimportancev2.df) <- c("gene", "importance")

total_gene_list <- arrange(variableimportancev2.df, importance)  

mostimportant20genes <- tail(total_gene_list , n = 20) 

print(mostimportant20genes) 



#Predictions 

# randomForest
predictiontotal <- predict(final.ranger.model, test.df)
prediction <- predictiontotal$predictions

actual <- test.df$classification 

pred_vec <- prediction == actual  

#accuracy = total correct predictions / total predictions made * 100

accuracy =length(pred_vec[pred_vec== TRUE])/length(pred_vec)

print(accuracy)


