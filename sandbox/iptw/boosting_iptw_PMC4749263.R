F.aac.iter=function(i,data,ps.model,ps.num,rep,criterion) {
  # i: number of iterations (trees)
  # data: dataset containing the treatment and the covariates
  # ps.model: the boosting model to estimate p(T_i|X_i)
  # ps.num: the estimated p(T_i)
  # rep: number of replications in bootstrap
  # criterion: the correlation metric used as the stopping criterion
  GBM.fitted=predict(ps.model,newdata=data,n.trees=floor(i),
                     type="response")
  ps.den=dnorm((data$T-GBM.fitted)/sd(data$T-GBM.fitted),0,1)
  wt=ps.num/ps.den
  aac_iter=rep(NA,rep)
  for (i in 1:rep){
    bo=sample(1:dim(data)[1],replace=TRUE,prob=wt)
    newsample=data[bo,]
    j.drop=match(c("T"),names(data))
    j.drop=j.drop[!is.na(j.drop)]
    x=newsample[,-j.drop]
    if(criterion=="spearman"|criterion=="kendall"){
      ac=apply(x, MARGIN=2, FUN=cor, y=newsample$T,
               method=criterion)
    } else if (criterion=="distance"){
      ac=apply(x, MARGIN=2, FUN=dcor, y=newsample$T)
    } else if (criterion=="pearson"){
      ac=matrix(NA,dim(x)[2],1)
      for (j in 1:dim(x)[2]){
        ac[j]=ifelse (!is.factor(x[,j]), cor(newsample$T, x[,j],
                                             method=criterion),polyserial(newsample$T, x[,j]))
      }
    } else print("The criterion is not correctly specified")
    aac_iter[i]=mean(abs(1/2*log((1+ac)/(1-ac))),na.rm=TRUE)
  }
  aac=mean(aac_iter)
  return(aac)
}

# Create the data frame for the covariates
x = data.frame(BMIZ, factor(DIABETZ1), G1BDESTM, G1WTCON,
               factor(INCOME1), M1AGE1, M1BMI, factor(M1CURLS), factor(M1CURMT),
               M1DEPRS, M1ESTEEM, M1GFATCN, M1GNOW, factor(M1GSATN),
               M1MFATCN, M1MNOW, factor(M1MSAT), factor(M1NOEX), M1OGIBOD,
               M1PCEAFF, M1PCEEFF, M1PCEEXT, M1PCEIMP, M1PCEPER, M1PDSTOT,
               M1RLOAD, factor(M1SMOKE), M1WGTTES, M1WTCON, M1YRED, factor(OBESE1),
               g1discal, g1obcdc, g1ovrcdc, g1pFM, g1wgttes, m1cfqcwt, m1cfqenc,
               m1cfqmon, m1cfqpwt, m1cfqrsp, m1cfqrst, m1cfqwtc, m1dis, m1hung,
               m1lim, m1picky, m1rest, m1zsav, m1zsweet)
# Find the optimal number of trees using Pearson/polyserial correlation
library(gbm)
library(polycor)
mydata=data.frame(T=M2WTCON,X=x)
model.num=lm(T~1,data=mydata)
ps.num=dnorm((mydata$T-model.num$fitted)/(summary(model.num))$sigma,0,1)
model.den=gbm(T~.,data=mydata, shrinkage = 0.0005,
              interaction.depth=4, distribution="gaussian",n.trees=20000)
opt=optimize(F.aac.iter,interval=c(1,20000), data=mydata, ps.model=model.den,
             ps.num=ps.num,rep=50,criterion="pearson")
best.aac.iter=opt$minimum
best.aac=opt$objective
# Calculate the inverse probability weights
model.den$fitted=predict(model.den,newdata=mydata,
                         n.trees=floor(best.aac.iter), type="response")
ps.den=dnorm((mydata$T-model.den$fitted)/sd(mydata$T-model.den$fitted),0,1)
weight.gbm=ps.num/ps.den
# Outcome analysis using survey package
library(survey)
dataset=data.frame(earlydiet,M2WTCON, weight.gbm)
design.b=svydesign(ids= ˜1, weights=˜weight.gbm, data=dataset)
fit=svyglm(earlydiet˜M2WTCON, family=quasibinomial(),design=design.b)
summary(fit)