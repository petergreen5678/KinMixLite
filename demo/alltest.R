require(KinMix)
data(test2data)
data(NGMDyes)

mixmarkers<-NULL

if(!is.null(mixmarkers)) epg<-epg[!is.na(match(epg$marker,mixmarkers)),]
cat('markers used in mixture:',unique(sort(epg$marker)),fill=60)

# set threshold C
C<-0.001

## 2-person mixture - baseline model

mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db,dyes=list(NGMDyes))
startpar<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
findmle<-FALSE
if(findmle)
{
	mlD<-mixML(mixD,startpar,trace=F) 
	pars<-mlD$mle
	cat('\nBaseline model maximised likelihood:',mlD$lik,'\n')
	cat('and MLEs:\n')
	print(mlD$mle)
	pars<-mlD$mle
} else {
	pars<-startpar
	print(pars)
	print(logL(mixD)(pars))
}
baseline<-logL(mixD)(pars)

## (1) Fit mixture model in which two contributors are related - in this case, half-sibs
## (code works in principle for any number of contributors, any number of which can be related in any way)

t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=2,C=rep(list(C),length(list(epg))),database=db,triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.IBD(mixD,c(0.5,0.5,0),,c(1,2))
size(mixD)
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1))

## (2a) Fit 2-person mixture model in which contributor 1 is parent of a typed individual Cgt

t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=2,C=rep(list(C),length(list(epg))),database=db,triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
cgtcaca<<-gt2aca(mixD,Cgt) 
replace.Ui.tables(mixD,cgtcaca,1)
size(mixD)
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1))

## (2b) Fit 2-person mixture model in which contributor 1 is parent of a typed individual Cgt, using rpt.typed.relatives

dadson<-list(pr=1,patt=matrix(c(1,2,3,1),1,4))
t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=2,C=rep(list(C),length(list(epg))),database=db,triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.typed.relatives(mixD,dadson,list(Cgt)) 
size(mixD)
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1))

## (3a) Fit 3-person mixture model in which three contributors are related - in this case, father and two sons 
## (code works in principle for any number of contributors, any number of which can be related in any way)
# need new baseline (and we use fixed parameters, MLE not attempted)

f2sons<-list(pr=rep(0.25,4),patt=matrix(c(1,2,3,1,3,1,1,2,3,1,3,2,1,2,3,1,4,1,1,2,3,1,4,2),4,6,byrow=TRUE))

t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=3,C=rep(list(C),length(list(epg))),database=db)
pars3<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.6,U2=0.3,U3=0.1)))
baseline3<-logL(mixD)(pars3)
size(mixD)

mixD<-DNAmixture(list(epg),k=3,C=rep(list(C),length(list(epg))),database=db,triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.IBD(mixD,IBD=f2sons,,1:3)
size(mixD)

log10LR<-(logL(mixD)(pars3)-baseline3)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1)) 

## (3b) Fit 3-person mixture model in which three contributors are related - in this case, father and two sons 
## - the same as above using loop.rpt.single.pattern.IBD function to loop over markers and IBD patterns - which should 
## enable bigger models
## (code works in principle for any number of contributors, any number of which can be related in any way)

t1<-Sys.time()
listdata<-list(epg)
loop.rpt.IBD(listdata,IBD=f2sons,pars=pars3,k=3,C=rep(list(C),length(listdata)),database=db)
t2<-Sys.time(); print(difftime(t2,t1))

## (3c) - using rpt.typed.relatives

t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=3,C=rep(list(C),length(list(epg))),database=db,triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.typed.relatives(mixD,IBD=f2sons,list(),1:3,jtyped=NULL)
size(mixD)
log10LR<-(logL(mixD)(pars3)-baseline3)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1)) 

## (4a) Fit 2-person mixture model, where half-sib of contributor 1 has been typed, using rpt.typed.relative
## (code works in principle for any number of contributors, and a single typed relative who can be related in any way)

t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db,dyes=list(NGMDyes),triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.typed.relative(mixD,Rgt,c(0.5,0.5,0.0),1)
size(mixD)
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1))

## (4b) .. using rpt.typed.relatives 

hs<-list(pr=c(0.5,0.5),patt=matrix(c(1,1,3,3,1,2,4,4,1,2,4,4),2,6))
t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db,dyes=list(NGMDyes),triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.typed.relatives(mixD,IBD=hs,typed.gts=list(Rgt))
size(mixD)
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1))

## (5a) Fit 2-person mixture model, where both parents of contributor 1 have been typed, using rpt.typed.parents
## (code works in principle for any number of contributors)


t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db,dyes=list(NGMDyes),triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.typed.parents(mixD,Mgt,Fgt,1)
size(mixD)
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1))

## (5b) ... using rpt.typed.relatives

c2pars <-
structure(list(pr = c(0.25, 0.25, 0.25, 0.25), patt = structure(c(1, 
1, 1, 1, 2, 2, 2, 2, 1, 3, 1, 3, 3, 1, 3, 1, 2, 2, 4, 4, 4, 4, 
2, 2), .Dim = c(4L, 6L))), .Names = c("pr", "patt"))
t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db,dyes=list(NGMDyes),triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.typed.relatives(mixD,IBD=c2pars,typed.gts=list(Mgt,Fgt))
size(mixD)
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1))

## (6a) Test for double paternity - fit 4 person mixture with contributor 1 the father of contributors 3 and 4, and 
## then set genotypes of contributors 3 and 4, using rpt.IBD
## (code works in principle for any number of contributors, any number of which can be related in any way, and with 
## any of these typed). Note U3 and U4 excluded from mixture by setting phi=0 for U3 and U4

# new baseline

pars4<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3,U3=0,U4=0)))
sons<-f2sons
sons$patt<-sons$patt[,3:6]

t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=4,C=rep(list(C),length(list(epg))),database=db,triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.IBD(mixD,IBD=sons,,c(3,4))
S1aca<<-gt2aca(mixD,S1gt)# ,1e-8) 
S2aca<<-gt2aca(mixD,S2gt)# ,1e-8) 
extra.findings<-make.findings(list(
	list('Aca',ind=3,aca='S1aca'),
	list('Aca',ind=4,aca='S2aca')
	))
baseline4<-logLX(mixD)(pars4)
t2<-Sys.time(); print(difftime(t2,t1)) 

## (6b) .. using rpt.typed.relatives

t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db,dyes=list(NGMDyes),triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.typed.relatives(mixD,IBD=f2sons,typed.gts=list(S1gt,S2gt))
size(mixD)
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1))

## (7) Fit 2-person mixture model, allowing UAF, using replace.tables.for.UAF
## (code works in principle for any number of contributors)

t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db,dyes=list(NGMDyes),triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
replace.tables.for.UAF(mixD,40)
size(mixD)
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1))

## (8) Fit 2-person mixture model where child of contributor 1 has been typed (paternity test), 
## allowing UAF, using replace.tables.for.UAF
## (code works in principle for any number of contributors, and a typed child or parent of any of them)

t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db,dyes=list(NGMDyes),triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
replace.tables.for.UAF(mixD,40,compile=FALSE)
Cgtaca<<-gt2aca(mixD,Cgt) 
add.relative.likd.node(mixD,Cgtaca,1)
extra.findings<-make.findings(list(
	list('Male',ind=1),
	list('Rlikd',aca='Cgtaca',cgt='Cgt',evid='Revid')
	))
size(mixD)
log10LR<-(logLX(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1))

## (9) Check that you get the same answer including a known individual in DNAmixture, or
## using rpt.typed.relatives to include a typed contributor (overlapping jcontr and jtyped)

S1prof<-make.profile(S1gt,'S1')
mixD<-DNAmixture(list(epg),k=3,K='S1',reference.profile=S1prof,C=list(C),database=db,dyes=list(NGMDyes))
size(mixD)
par1<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.5,U2=0.3,S1=0.2)))
logL(mixD)(par1)

t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=3,C=rep(list(C),length(list(epg))),database=db,triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.typed.relatives(mixD,list(patt=1:6),list(S1gt),inds=1:3,jtyped=3,jcontr=1:3) 
size(mixD)
par2<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.5,U2=0.3,U3=0.2)))
logL(mixD)(par2)
t2<-Sys.time(); print(difftime(t2,t1))

## (10) Check that rpt.typed.relatives works with no one typed - recovering baseline

t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=2,C=rep(list(C),length(list(epg))),database=db,triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.typed.relatives(mixD,list(patt=1:4),list(),inds=1:2,jtyped=NULL,jcontr=1:2) 
size(mixD)
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1))

## (11) Check that rpt.typed.relatives with no one typed gets same answer as rpt.IBD

IBD<-list(pr=1,patt=matrix(c(1,2,1,3,2,4),1,6))

t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=3,C=rep(list(C),length(list(epg))),database=db,triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.IBD(mixD,IBD,,1:3)
size(mixD)
log10LR<-(logL(mixD)(pars3)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1))

t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=3,C=rep(list(C),length(list(epg))),database=db,triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.typed.relatives(mixD,IBD,list(),inds=1:3,jtyped=NULL,jcontr=1:3) 
size(mixD)
log10LR<-(logL(mixD)(pars3)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1))

## (12) Testing that autozygosity works

t1<-Sys.time()
mixD<-DNAmixture(list(epg),k=3,C=rep(list(C),length(list(epg))),database=db,triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.typed.relatives(mixD,list(patt=c(1,1,1,2)),list(),inds=1:2,jtyped=NULL,jcontr=1:2)
size(mixD)
log10LR<-(logL(mixD)(pars3)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
t2<-Sys.time(); print(difftime(t2,t1))
