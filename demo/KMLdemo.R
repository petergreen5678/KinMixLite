library(KinMixLite)
data(test2data)

t1<-Sys.time()

## Fit 2-person mixture - baseline model

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
baseline<-logL(mixD)(pars)
baseline

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
mlD <- mixML(mixD, pars)
print(mlD$mle)
pars<-mlD$mle

baseline<-logL(mixD)(pars)
baseline

## Fit 2-person mixture model in which contributor 1 is parent of a typed individual Cgt

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE)
rpt.IBD(mixD,,list(Cgt)) 
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## Fit 2-person mixture model in which contributor 1 is father of a typed individual Cgt 
## with mother Mgt

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE)
rpt.IBD(mixD,,list(Mgt,Cgt)) 
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## Fit 2-person mixture model in which contributors are two parents of a child with 
## genotype Cgt, and a parent of one of them has genotype Rgt. Note the encoding of allele 
## labels to reduce the complexity of the IBD pattern distribution IBD.

IBD<-list(patt=rbind(c(1,3,2,4,1,2,1,5),c(1,3,2,4,1,2,3,5)))

mixD1<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE)
rpt.IBD(mixD1,IBD,list(Cgt,Rgt),1:2) 
log10LR<-(logL(mixD1)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## the same, with individuals and relationships denoted by character tags

mixD2<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE)
rpt.IBD(mixD2,IBD,list(c=Cgt,gf=Rgt),targets=c('f','m','c','gf'),contribs=c('f','m'),quiet=TRUE) 
log10LR<-(logL(mixD2)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## and with m as U1 and f as U2:

mixD2<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE)
rpt.IBD(mixD2,IBD,list(c=Cgt,gf=Rgt),targets=c('f','m','c','gf'),contribs=c('m','f'),quiet=TRUE) 
log10LR<-(logL(mixD2)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## show 'deconvolution'

summary(map.genotypes(mixD2,pmin=0.001,U=1))

## plot one of the DAGs

if(is(mixD2$dom$D12,'gRv')) plot(mixD2$dom$D12$net,type='dag')

t2<-Sys.time(); print(difftime(t2,t1))
