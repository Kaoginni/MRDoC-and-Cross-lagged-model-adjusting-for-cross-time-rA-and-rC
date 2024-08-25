# Script Name: Causal association between Cognitive ability and cognitive stimulation
# By Olakunle Oginni

# To determine the direction of causation in the relationship between Cognitive ability (Vocab, grammar, parent-administered parca and parent-reported parca)
# and cognitive stimulation (play variables) using several models (see below)
# I will specify two latent factors:
# i.   A Cognitive development/ability factor indicated by sum scores for:
#		a. Overall PARCA scores (at ages 3, 4 years)
# ii.  An Cognitive stimulation factor indicated by sum scores for:
#		a. talkrhym (trhym/ttalk/gtapes)
#		b. playbook (tbook/gbooks)
#		c. playgames (gpuzz/gboard)
# Zygosity variable: 		zyg (1=MZ, 2=DZ)
# ------------------------------------------------------------------------------------

# PART 1: By Zygosity groups: MZ/DZ 
#_______________________________________

# PHENOTYPIC MODELS
# MODEL 1: Constrained Correlation model
# MODEL 2: Bivariate factor model

# BIOMETRIC MODELS
# MODEL 3: UNIVARIATE ACE FACTOR MODEL
# MODEL 4: DRIRECTION OF CAUSATION (DoC) MODEL
# MODEL 5a: MENDELIAN RANDOMISATION-DIRECTION OF CAUSATION (MRDoC) MODEL - Using latent constructs
# MODEL 5b: MENDELIAN RANDOMISATION-DIRECTION OF CAUSATION (MRDoC) MODEL - Using observed variables
# MODEL 6a: CROSS-LAGGED MODEL
# MODEL 6b: CROSS-LAGGED MODEL ADJUSTED FOR CROSS-TIME rA and rC (GAUSSIAN)
# MODEL 6c: CROSS-LAGGED MODEL ADJUSTED FOR CROSS-TIME rA and rC (CHOLESKY DECOMPOSITION)

rm(list=ls())
ls()

#check you have the packages:
search()

library(psych)
library(OpenMx)
library(Hmisc)			
library(dplyr)
library(data.table)
library(stringr)
library(R.utils)
library(tidyverse)
library(foreign)

mxOption(NULL, "Default optimizer", "SLSQP") 

# ******************************************************************************************************************************************
# (1) Read in data file and check descriptive statistics
# ******************************************************************************************************************************************

# Read in data 
# Data <- 
# psych::describe(Data)

# ******************************************************************************************************************************************
# (2) Explore and prepare variables of interest: recode, regress-out age and sex, check distribution and transform if necessary.
# ******************************************************************************************************************************************

# i.   A Cognitive development/ability factor indicated by sum scores for:
#		a. Overall PARCA (at ages 3, 4 years)

# ii.  An Cognitive stimulation factor indicated by sum scores for:
#		a. TalkRhym
#		b. Playbooks
#		c. Playgames
# ---------
# Data preparation
# ---------

# __(1)________________________________________________________________________________________________
# Constrained Correlation Model (Gaussian decomposition)
# Restrictions: means and variances equated across birth-order & zygosity groups; 
# One set of Cross-trait cor, symmetric Cross-trait Cross-twin correlations (MZ and DZ)
# ______________________________________________________________________________________________________

nv		<- 9				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv-1)/2)		# number of free elements in a correlation matrix nv*nv
Groups	<- c("mz", "dz")

Vars		<- c("PRSEA","Parca3","Parca4","TalkRhym3","PlayBook3","PlayGames3","TalkRhym4","PlayBook4","PlayGames4")
selVars	<- c("PRSEA1","Parca31","Parca41","TalkRhym31","PlayBook31","PlayGames31","TalkRhym41","PlayBook41","PlayGames41",
			"PRSEA2","Parca32","Parca42","TalkRhym32","PlayBook32","PlayGames32","TalkRhym42","PlayBook42","PlayGames42")

mzData	<- subset(TWINdata2, zyg1%in%c(1)|zyg2%in%c(1), selVars)
dzData	<- subset(TWINdata2, zyg1%in%c(2)|zyg2%in%c(2), selVars)

psych::describe(mzData)
psych::describe(dzData)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
 
(Stmean		<-colMeans(mzData[,1:nv],na.rm=TRUE))
(Stsd 		<-sapply(mzData[,1:nv],sd, na.rm=TRUE))
(StWithinperson	<-vechs(cor(mzData[,1:nv],use="complete")))


StBetweenMZ <-c(1,.7,.7,.7,.7,.7,.7,.7,.7,			#Cross-twin PRS correlation fixed to 1 for MZ twin pairs and 0.5 for DZ twin pairs
			.7,.7,.7,.7,.7,.7,.7,.7,
			.7,.7,.7,.7,.7,.7,.7,
			.7,.7,.7,.7,.7,.7,
			.7,.7,.7,.7,.7,
			.7,.7,.7,.7,
			.7,.7,.7,
			.7,.7,
			.7)

StBetweenDZ <-c(.5,.4,.4,.4,.4,.4,.4,.4,.4,
			.4,.4,.4,.4,.4,.4,.4,.4,
			.4,.4,.4,.4,.4,.4,.4,
			.4,.4,.4,.4,.4,.4,
			.4,.4,.4,.4,.4,
			.4,.4,.4,.4,
			.4,.4,.4,
			.4,.4,
			.4)

pathsMZ 	<-c(F,T,T,T,T,T,T,T,T,
			T,T,T,T,T,T,T,T,
			T,T,T,T,T,T,T,
			T,T,T,T,T,T,
			T,T,T,T,T,
			T,T,T,T,
			T,T,T,
			T,T,
			T)

# Create Labels for Column and Diagonal Matrices
(mLabs	<- paste("m",1:nv,sep=""))
(sdLabs	<- paste("sd",1:nv,sep=""))

# Create Labels for a Correlation Matrix
(rphLabs	<- paste("r",1:ncor,sep=""))

# Create Labels for Lower Triangular Matrices
(MZbLabs 	<- paste("rmz", do.call(c, sapply(seq(1, nv), function(x){ paste(x:nv, x,sep="") })), sep=""))
(DZbLabs 	<- paste("rdz", do.call(c, sapply(seq(1, nv), function(x){ paste(x:nv, x,sep="") })), sep=""))

# ______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
# ______________________________________________________________________________________________________
 
# Specify Matrices to hold parameter estimates for MZ and DZ data

# elements for one overall Means 
Means		<-mxMatrix("Full", 1, ntv, free=T, values=Stmean, labels=c(mLabs,mLabs), name="expMean") 

# elements for the SD
sdZ		<-mxMatrix("Zero", nv, nv, free=F, name="padding")
sdG		<-mxMatrix("Diag", nv, nv, free = TRUE, values = Stsd, labels=c(sdLabs), name="SD") 
sdT		<-mxAlgebra(rbind(cbind(SD,padding), cbind(padding,SD)), name="SDtwin")

# elements for the correlations
Rph		<-mxMatrix("Stand", nv, nv, free = TRUE, values = .7, labels = rphLabs, name="Rwithin") 
MZb		<-mxMatrix("Symm", nv, nv, free = pathsMZ, values = StBetweenMZ, labels=MZbLabs, name="RbetweenMZ") 
DZb		<-mxMatrix("Symm", nv, nv, free = pathsMZ, values = StBetweenDZ, labels=DZbLabs, name="RbetweenDZ") 
corMZ		<-mxAlgebra(rbind(cbind(Rwithin,RbetweenMZ), cbind(RbetweenMZ,Rwithin)), name="RMZ")
corDZ		<-mxAlgebra(rbind(cbind(Rwithin,RbetweenDZ), cbind(RbetweenDZ,Rwithin)), name="RDZ")

# generate expected covariance matrices
covMZ		<-mxAlgebra(SDtwin %&% RMZ, name="expCovMZ")
covDZ		<-mxAlgebra(SDtwin %&% RDZ, name="expCovDZ")

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="expCovMZ", means="expMean", dimnames=selVars )
objDZ		<- mxExpectationNormal( covariance="expCovDZ", means="expMean", dimnames=selVars )

fitFunction <- mxFitFunctionML()

# Combine Groups
pars		<- list(Means, sdZ, sdG, sdT, Rph)
modelMZ	<- mxModel(pars, MZb, corMZ, covMZ, dataMZ, objMZ, fitFunction,  name="MZ" )
modelDZ	<- mxModel(pars, DZb, corDZ, covDZ, dataDZ, objDZ, fitFunction,  name="DZ" )
minus2ll	<- mxAlgebra(expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<- mxFitFunctionAlgebra( "m2LL" )
ciW		<-mxCI("MZ.Rwithin")
ciBM		<-mxCI("MZ.RbetweenMZ")
ciBD		<-mxCI("DZ.RbetweenDZ")
ConCor4Model	<- mxModel("Con4", pars, modelMZ, modelDZ, minus2ll, obj, ciW, ciBM, ciBD )

# -------------------------------------------------------------------------------
# 1 RUN Constrained Correlation Model by Zygosity
#ConCor4Fit	<- mxTryHard(ConCor4Model, intervals=F, extraTries = 200, showInits=F)
ConCor4Fit	<- mxRun(ConCor4Model, intervals=T)
(ConCor4Sum	<- summary(ConCor4Fit))

# To get the output values 
mxEval(Con.expMean, ConCor3Fit)
mxEval(Con.SD, ConCor3Fit)

mxEval(Con.Rwithin, ConCor3Fit)
mxEval(MZ.RbetweenMZ, ConCor3Fit)
mxEval(DZ.RbetweenDZ, ConCor3Fit)
mxEval(MZ.RMZ, ConCor3Fit)
mxEval(DZ.RDZ, ConCor3Fit)

mxEval(MZ.expCovMZ, ConCor3Fit)
mxEval(DZ.expCovDZ, ConCor3Fit)



#*******************************************************************************************************
# __(2) - Bivariate -TR_______________________________________________________________________________
# Phenotypic Covariance Model across Latent Constructs: Cognitive ability and cognitive stimulation
# Restrictions: means and variances equated across birth-order & zygosity groups;
# One set of factor loadings; one set of correlational paths between the factors; one set of error terms
# We estimate the factor variances, giving them a scale by fixing the loading on the 1st variable to 1
# This model specifies a full var/cov structure between the latent factors for MZ and DZ twins 
#______________________________________________________________________________________________________

nv		<- 4				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 2				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nfcor		<- (nfact*(nfact-1)/2)	# number of free elements in a correlation matrix nfact*nfcat
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv-1)/2)		# number of free elements in a correlation matrix nv*nv
Groups	<-c("mz", "dz")

Vars		<- c("Parca3","Parca4","TalkRhym3","TalkRhym4")
selVars	<- c("Parca31","Parca41","TalkRhym31","TalkRhym41",
			"Parca32","Parca42","TalkRhym32","TalkRhym42")

mzData	<- subset(TWINdata2, zyg1%in%c(1)|zyg2%in%c(1), selVars)
dzData	<- subset(TWINdata2, zyg1%in%c(2)|zyg2%in%c(2), selVars)

psych::describe(mzData)
psych::describe(dzData)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)

(Stmean	<-colMeans(mzData[,1:nv],na.rm=TRUE))
StM 		<-c(Stmean, Stmean)

(LabM		<- paste("m",1:nv,sep=""))
MLabs		<-c(LabM,LabM) 

#(LabEr	<- paste("e",1:nv,sep=""))
(LabEr	<-c("e2","e2","e4","e4"))

MZbErLabs	<-c("remz1","remz2","remz3","remz4")
DZbErLabs	<-c("redz1","redz2","redz3","redz4")

Ste		<- c(.6,.6,.6,.6)

# Create Labels for the Factor parameters
(sdLabs	<- paste("sd",1:nfact,sep=""))	# SD
(rphLabs	<- paste("r",1:nfcor,sep=""))		# Correlation Matrix Factors within person 
(MZbLabs 	<- paste("rmz", do.call(c, sapply(seq(1, nfact), function(x){ paste(x:nfact, x,sep="") })), sep=""))	# Lower Triangular Matrices
(DZbLabs 	<- paste("rdz", do.call(c, sapply(seq(1, nfact), function(x){ paste(x:nfact, x,sep="") })), sep=""))

# Create Labels for the Factor Loadings (1st loadings fixed to 1)
PatFl	<- c(F,T,F,F,
	     F,F,F,T)

StFl	<- c(1,0,0,0,
	     0,0,1,.6)

LabFl	<- c('l1','l2',NA,NA,
	      NA,NA,'l3','l4')

# Free parameters
(Pat  <- c( rep(TRUE, nv)))
(Pate <- c( rep(TRUE, nv)))

# ______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
# ______________________________________________________________________________________________________

Mean		<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=c(Pat, Pat), values=c(StM), labels=c(MLabs), name="expm" ) 

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Ze42		<-mxMatrix("Zero", nv, nfact, free=F, name="Z42")
LoadTw	<-mxAlgebra(rbind(cbind(FactL,Z42), cbind(Z42,FactL)), name="FactLTw")

ErPath	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=Pate, values=Ste, labels=LabEr, name="Erp" )
Er		<-mxAlgebra(Erp %*% t(Erp), name="Error")

# Here I specify crosstwin residual correlations as follows:
MZbEr		<-mxMatrix("Diag", nv, nv, free=c(T,T,T,T), values=c(.4,.4,.4,.4), lbound=-.999, ubound=.999, labels=MZbErLabs, name="RbetweenMZEr") 
DZbEr		<-mxMatrix("Diag", nv, nv, free=c(T,T,T,T), values=c(.2,.2,.2,.2), lbound=-.999, ubound=.999, labels=DZbErLabs, name="RbetweenDZEr")

ErMZ		<-mxAlgebra(Erp %&% RbetweenMZEr, name="BErrorMZ")
ErDZ		<-mxAlgebra(Erp %&% RbetweenDZEr, name="BErrorDZ")

#Ze6		<-mxMatrix("Zero", nv, nv, free=F, name="Z6")
ErTwMZ	<-mxAlgebra(rbind(cbind(Error,BErrorMZ), cbind(BErrorMZ,Error)), name="ErrorTwMZ")
ErTwDZ	<-mxAlgebra(rbind(cbind(Error,BErrorDZ), cbind(BErrorDZ,Error)), name="ErrorTwDZ")

# elements for the SD of Factors
Id2		<-mxMatrix("Iden", 2, 2, free=F, name="I2")
sdF		<-mxMatrix("Diag", nfact, nfact, free=c(T,T), values=c(1,1), labels=sdLabs, name="SDf") 
sdFTw		<-mxAlgebra(I2 %x% SDf, name="SDftwin")

# elements for the correlations of Factors
Rph		<-mxMatrix("Stand", nfact, nfact, free=TRUE, values=c(.5), lbound=-.999, ubound=.999, labels=rphLabs, name="Rwithin") 
MZb		<-mxMatrix("Symm", nfact, nfact, free=c(T,T,T), values=c(.5,.6,.5), lbound=-.999, ubound=.999, labels=MZbLabs, name="RbetweenMZ") 
DZb		<-mxMatrix("Symm", nfact, nfact, free=c(T,T,T), values=c(.3,.3,.2), lbound=-.999, ubound=.999, labels=DZbLabs, name="RbetweenDZ") 
FactCorMZ	<-mxAlgebra(rbind(cbind(Rwithin,RbetweenMZ), cbind(RbetweenMZ,Rwithin)), name="RMZ")
FactCorDZ	<-mxAlgebra(rbind(cbind(Rwithin,RbetweenDZ), cbind(RbetweenDZ,Rwithin)), name="RDZ")

# Generate expected Covariance matrices of Factors
FactCovMZ	<-mxAlgebra(SDftwin %&% RMZ , name="expFactCovMZ")
FactCovDZ	<-mxAlgebra(SDftwin %&% RDZ , name="expFactCovDZ")

## This second step then derives the var/cov matrix of the observed/measured variables in terms of the variance/covariances of the latent factors and the Factor Loadings
covMZ		<-mxAlgebra( expression= FactLTw  %&% expFactCovMZ , name="ExpCovMZ" )
covDZ		<-mxAlgebra( expression= FactLTw  %&% expFactCovDZ , name="ExpCovDZ" )

## Finally, we derive the total expected variance/covariances for the measured variables which go in the models
TOTcovMZ	<-mxAlgebra( expression= ExpCovMZ + ErrorTwMZ , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= ExpCovDZ + ErrorTwDZ , name="TOTexpCovDZ" )

# Standardizing parameters **********************

# Standardized Factor Loadings
StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% expFactCovMZ[1:2,1:2] / TOTexpCovMZ[1:4,1:4])) , name="StandFact" )

# Standardise error paths of the measured variables
StEr		<-mxAlgebra( expression= sqrt(diag2vec( Error/TOTexpCovMZ[1:4,1:4])), name="StandEr" )

# ************************************

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expm", dimnames=selVars)
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expm", dimnames=selVars)

fitFunction <- mxFitFunctionML()

# Combine Groups
pars1		<-list(Mean, Load, Ze42, LoadTw, ErPath, Er, Id2, sdF, sdFTw, Rph)
modelMZ	<-mxModel(pars1, MZb, FactCorMZ, FactCovMZ, covMZ, MZbEr, ErMZ, ErTwMZ, TOTcovMZ, dataMZ, objMZ, fitFunction, StFL, StEr, name="MZ" )
modelDZ	<-mxModel(pars1, DZb, FactCorDZ, FactCovDZ, covDZ, DZbEr, ErDZ, ErTwDZ, TOTcovDZ, dataDZ, objDZ, fitFunction, name="DZ" )
minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
cist1		<-mxCI (c ('MZ.Rwithin[2,1]'))
cist2		<-mxCI (c ('MZ.StandFact'))
cist3		<-mxCI (c ('MZ.StandEr'))
cist4		<-mxCI (c ('MZ.RbetweenMZ'))
cist5		<-mxCI (c ('DZ.RbetweenDZ'))
cist6		<-mxCI (c ('MZ.RbetweenMZEr'))
cist7		<-mxCI (c ('DZ.RbetweenDZEr'))

PhC1Model	<-mxModel("PhC1", modelMZ, modelDZ, minus2ll, obj, cist1, cist2, cist3, cist4, cist5, cist6, cist7) 

# --------------------------------------------------------------------------------------------------------------------------------
# 2 RUN Phenotypic Fact Covariance Model by Zygosity

PhC1Fit	<-mxTryHard(PhC1Model, intervals=T, showInits=F)
#PhC1Fit	<-mxRun(PhC1Model, intervals=T)
(PhC1Summ	<-summary(PhC1Fit, verbose=F))

mgPhC1Fitfun <- mxFitFunctionMultigroup(c("MZ","DZ"))
mgPhC1Model 	<- mxModel("mgPhC1", modelMZ, modelDZ,mgPhCFitfun)
mgPhC1Fit 	<- mxRun(mgPhC1Model)
(mgPhC1Summ	<-summary(mgPhC1Fit))

# This checks identification specifically
(Id		<-mxCheckIdentification(mgPhC1Fit))


mxEval(MZ.Rwithin, PhC1Fit)
mxEval(MZ.RbetweenMZ, PhC1Fit)
mxEval(DZ.RbetweenDZ, PhC1Fit)
mxEval(MZ.SDf, PhC1Fit)

mxEval(MZ.expFactCovMZ, PhC1Fit)
mxEval(DZ.expFactCovDZ, PhC1Fit)

mxEval(MZ.ExpCovMZ, PhC1Fit)
mxEval(DZ.ExpCovDZ, PhC1Fit)

mxEval(MZ.TOTexpCovMZ, PhC1Fit)
mxEval(DZ.TOTexpCovDZ, PhC1Fit)

mxEval(MZ.FactL, PhC1Fit)
mxEval(MZ.StandFact, PhC1Fit)

mxEval(MZ.Error, PhC1Fit)
mxEval(MZ.StandEr, PhC1Fit)
mxEval(MZ.RbetweenMZEr, PhC1Fit)
mxEval(DZ.RbetweenDZEr, PhC1Fit)


#****************************************************************************************************************************
# __(3) _____________________________________________________________________________________________________________________
# ACE Factor MODEL by zygosity
# NO causal paths between Factors; A, C and E latent factors have Cholesky Structure
# + Asp, Csp and Esp in the bottom with constraints to identify the model on top
# Correlation between Phenotypic Factors only due to shared A, C and E influences
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1), 
#_____________________________________________________________________________________________________________________________

nv		<- 2				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 1				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
Groups	<-c("mz", "dz")
Vars		<- c("Parca3","Parca4")
selVars	<- c("Parca31","Parca41","Parca32","Parca42")

mzData	<- subset(TWINdata2, zyg1%in%c(1)|zyg2%in%c(1), selVars)
dzData	<- subset(TWINdata2, zyg1%in%c(2)|zyg2%in%c(2), selVars)

psych::describe(mzData)
psych::describe(dzData)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabs	<- paste("m",1:nv,sep=""))
(Stmean	<- colMeans(mzData[,1:nv],na.rm=TRUE))
(Stsd 	<- sapply(mzData[,1:nv],sd, na.rm=TRUE))
(PatM		<- c(TRUE,TRUE))

# Create Labels for Diagonal Matrices
# To identify this model we need to equate the sp effects of var and 2 and fix the Sp of last variable to 0)
(LabEs	<- c('es1','es1'))
(LabAs	<- c('as1','as1'))
(LabCs	<- c('cs1','cs1'))

PatSp		<- c(TRUE,TRUE)
StSpa		<- c(.5,.5)
StSpc		<- c(.5,.5)
StSpe		<- c(.5,.5)

# all 1st loadings fixed to 1
PatFl		<- c(F,T)

StFl		<- c(1,.5)

LabFl		<- c('l1','l2')

# ______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
# ______________________________________________________________________________________________________

Means		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmean,Stmean), labels=c(mLabs,mLabs), name="expMean") 

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTw	<-mxAlgebra(I2%x%FactL, name="FactLTw")
 
# Define the matrix to hold the A and C effects: Specific 
PathsAs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSp, values=StSpa, labels=LabAs, name="as" )
PathsCs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSp, values=StSpc, labels=LabCs, name="cs" )
PathsEs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSp, values=StSpe, labels=LabEs, name="es" )
covAs		<-mxAlgebra( expression= as %*% t(as), name="As" )
covCs		<-mxAlgebra( expression= cs %*% t(cs), name="Cs" )
covEs		<-mxAlgebra( expression= es %*% t(es), name="Es" )
covPs		<-mxAlgebra( expression= As+Cs+Es, name="Vs" )

# Define the matrices to hold the A and C effects: Common 
PathsAc		<-mxMatrix(type="Lower", nrow=nfact, ncol=nfact, free=TRUE, values=.6, labels=c("a11"), name="a_c" )
PathsCc		<-mxMatrix(type="Lower", nrow=nfact, ncol=nfact, free=TRUE, values=.3, labels=c("c11"), name="c_c" )
PathsEc		<-mxMatrix(type="Lower", nrow=nfact, ncol=nfact, free=TRUE, values=.6, labels=c("e11"), name="e_c" )
covAc		<-mxAlgebra( expression= a_c %*% t(a_c), name="Ac" )
covCc		<-mxAlgebra( expression= c_c %*% t(c_c), name="Cc" )
covEc		<-mxAlgebra( expression= e_c %*% t(e_c), name="Ec" )
covPc		<-mxAlgebra( expression= Ac+Cc+Ec, name="Vc" )

# Var-Cov of measured vars in terms of latent factors and AC, Cc, and Ec
FcovMZ		<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(Vc, Ac+Cc), cbind(Ac+Cc, Vc))) , name="expFCovMZ" )#This traces the path from vars to factors and back to vars
FcovDZ		<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(Vc, .5%x%Ac+Cc), cbind(.5%x%Ac+Cc, Vc))) , name="expFCovDZ" )

SpcovMZ		<-mxAlgebra( expression= rbind (cbind(Vs, As+Cs), cbind(As+Cs, Vs)) , name="expSpCovMZ" )
SpcovDZ		<-mxAlgebra( expression= rbind (cbind(Vs, .5%x%As+Cs), cbind(.5%x%As+Cs, Vs)) , name="expSpCovDZ" )

TOTcovMZ		<-mxAlgebra( expression= expFCovMZ + expSpCovMZ , name="TOTexpCovMZ" )
TOTcovDZ		<-mxAlgebra( expression= expFCovDZ + expSpCovDZ , name="TOTexpCovDZ" )

# *******************************************************************************************************
# Calculator

# Standardize the Total var/covariances matrices of the observed variables
Id4		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I4" )
Id1		<-mxMatrix(type="Iden",	nrow=1, ncol=1, name="I1" )
Rfactmz	<-mxAlgebra( expression= solve(sqrt(I4*TOTexpCovMZ)) %&% TOTexpCovMZ, name="FactcorMZ" )
Rfactdz	<-mxAlgebra( expression= solve(sqrt(I4*TOTexpCovDZ)) %&% TOTexpCovDZ, name="FactcorDZ" )
RfactA	<-mxAlgebra( expression= solve(sqrt(I1*Ac)) %&% Ac, name="rA" )
RfactC	<-mxAlgebra( expression= solve(sqrt(I1*Cc)) %&% Cc, name="rC" )
RfactE	<-mxAlgebra( expression= solve(sqrt(I1*Ec)) %&% Ec, name="rE" )
RfactV	<-mxAlgebra( expression= solve(sqrt(I1*Vc)) %&% Vc, name="rPfact" )

# Standardize the Common Effects
stcovAc	<-mxAlgebra( expression= Ac/Vc, name="stAc" )
stcovCc	<-mxAlgebra( expression= Cc/Vc, name="stCc" )
stcovEc	<-mxAlgebra( expression= Ec/Vc, name="stEc" )

# Standardize the Specific Effects
stcovAs	<-mxAlgebra( expression= sqrt(As/( (FactL %&% Vc) +Vs)), name="stAs" )
stcovCs	<-mxAlgebra( expression= sqrt(Cs/( (FactL %&% Vc) +Vs)), name="stCs" )
stcovEs	<-mxAlgebra( expression= sqrt(Es/( (FactL %&% Vc) +Vs)), name="stEs" )

# Standardized Effects of Individual variables from the factors (Variance components) above
stAvar	<-mxAlgebra( expression= (FactL %&% Ac)/( (FactL %&% Vc) +Vs), name="stAvariables" )
stCvar	<-mxAlgebra( expression= (FactL %&% Cc)/( (FactL %&% Vc) +Vs), name="stCvariables" )
stEvar	<-mxAlgebra( expression= (FactL %&% Ec)/( (FactL %&% Vc) +Vs), name="stEvariables" )

# Standardized Factor Loadings
StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% Vc[1,1] / TOTexpCovMZ[1:2,1:2])) , name="StandFact" )

# *******************************************************************************************************

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expMean", dimnames=selVars)
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expMean", dimnames=selVars)

fitFunction <- mxFitFunctionML()
 
# Combine Groups
pars1		<-list(Means,Load,LoadTw,PathsAs,PathsCs,PathsEs,covAs,covCs,covEs,covPs,Id2,Id4,Id1, stAvar, stCvar, stEvar)
pars2		<-list(PathsAc,PathsCc,PathsEc,covAc,covCc,covEc,covPc,stcovAc,stcovCc,stcovEc, stcovAs, stcovCs, stcovEs)
modelMZ	<-mxModel(pars1, pars2, FcovMZ, SpcovMZ, TOTcovMZ, dataMZ, objMZ, Rfactmz, RfactA, RfactC, RfactE, RfactV, fitFunction, StFL, name="MZ" )
modelDZ	<-mxModel(pars1, pars2, FcovDZ, SpcovDZ, TOTcovDZ, dataDZ, objDZ, Rfactdz, fitFunction, name="DZ" )
minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
cistFL	<-mxCI (c ('MZ.StandFact'))
cistFc	<-mxCI (c ('MZ.stAc','MZ.stCc','MZ.stEc') ) 	# standardized var comp from Common feactors	
cistVs	<-mxCI (c ('MZ.stAs','MZ.stCs','MZ.stEs') ) 	# standardized var comp from specific Factors
cistvars	<-mxCI (c ('MZ.stAvariables','MZ.stCvariables','MZ.stEvariables'))
UniACEfact1Model	<-mxModel("Uniacefact1", modelMZ, modelDZ, minus2ll, obj, cistFc, cistFL, cistVs, cistvars) 

# --------------------------------------------------------------------------------------------------------------------------------
# 4 RUN ACE Factor Model: Cholesky (by Zygosity)

UniACEfact1Fit	<-mxRun(UniACEfact1Model, intervals=T)
(UniACEfact1Summ	<-summary(UniACEfact1Fit))


# __(4) _____________________________________________________________________________________________________________________
#
# DOC model based on Nathan Gillespie's script from the Boulder course
# PARCA and TR
#****************************************************************************************************************************

nv		<- 4				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 2				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
Groups	<-c("mz", "dz")
Vars		<- c("Parca3","Parca4","TalkRhym3","TalkRhym4")
selVars	<- c("Parca31","Parca41","TalkRhym31","TalkRhym41",
			"Parca32","Parca42","TalkRhym32","TalkRhym42")

#Vars		<- c("Parca3","Parca4","TalkRhym34","PlayBook34","PlayGames34")
#selVars	<- c("Parca31","Parca41","TalkRhym341","PlayBook341","PlayGames341",

mzData	<- subset(TWINdata2, zyg1%in%c(1)|zyg2%in%c(1), selVars)
dzData	<- subset(TWINdata2, zyg1%in%c(2)|zyg2%in%c(2), selVars)

psych::describe(mzData)
psych::describe(dzData)

# Print Descriptive Statistics
 round(colMeans(mzData,na.rm=TRUE),4)
 round(colMeans(dzData,na.rm=TRUE),4)
 
# Matrices in the full DCS model

# ? psi 	= A, C & E effects on common pathways (aka innovations in autoregression models)
# ? beta 	= auto-regression 'causal pathways
# I identity
# ? lamba   = factor loadings from latent factor 'true scores' to observed variables
# ? epsilon = meaurement error (specify residual A, C, E effects)

# ? psi - A, C & E effects on common pathways (aka innovations in autoregression models)
 # Sources of A, C & E variance on the latent factors or common pathways
 aFree		<- c(T,T,T,T)
 cFree		<- c(T,T,T,T) 
 eFree		<- c(T,T,T,T)  
 aVals		<- c(0.6,0.3,0.3,0.6)
 cVals		<- c(0.2,0.1,0.1,0.3)
 eVals		<- c(0.6,0.1,0.1,0.6)
 psi_a 		<- mxMatrix(type = "Symm", nrow = nfact,ncol=nfact, labels=c("a11","a21","a21","a22"), free=aFree, values=aVals, name="psi_a")
 psi_c 		<- mxMatrix(type = "Symm", nrow = nfact,ncol=nfact, labels=c("c11","c21","c21","c22"), free=cFree, values=cVals, name="psi_c") 
 psi_e 		<- mxMatrix(type = "Symm", nrow = nfact,ncol=nfact, labels=c("e11","e21","e21","e22"), free=eFree, values=eVals, name="psi_e") 
   


# Beta - causal parameters between latent factors / common pathways
 # NB: direction of causation goes DOWN the column & OUT along the row
 Bfree   	<- matrix(c(F,F,  
 			      F,F),2,2, byrow=T)
 Bvals   	<- matrix(c(0,0,
	                  0,0),2,2, byrow=T)  	
 Blabs	 	<- matrix(c(NA,"CP2toCP1",
	                    "CP1toCP2",NA),2,2, byrow=T)							  		  	   
 beta    	<- mxMatrix(type="Full", nrow=nfact, ncol=nfact, labels=Blabs, free=Bfree, values=Bvals, lbound=-1, ubound=1, name = "beta" ) # 
 I       	<- mxMatrix(type="Iden", nrow=nfact, ncol=nfact, name = "I")
 
 
# ? lamba - factor loadings from the latent true scores or common pathways to the observed variables
 Ffree   	<- matrix(c(  F,F,
 				  T,F,
 				  F,F,
 				  F,T),4,2, byrow=T)						  

 Fvals   	<- matrix(c(  1,0.0,
				  0.85,0.0,
				  0.0,1,
				  0.00,0.85),4,2, byrow=T)
						  					  
 Flabs   	<- matrix(c(  "f11",NA,
 				  "f21",NA,
 				  NA,"f32",
 				  NA,"f42"),4,2, byrow=T)
						  						  
 lambda   	<- mxMatrix(type="Full", nrow=nv, ncol=nfact, free = Ffree, values = Fvals, labels=Flabs, name = "lambda")
 

# ? epsilon - measurement error/residuals or transient variance in each observed/manifest item/variable (not explained by the latent true scores)
 res_a_labs	<- c("res_a1","res_a2","res_a3","res_a4")
 res_c_labs	<- c("res_c1","res_c2","res_c3","res_c4") 
 res_e_labs	<- c("res_e1","res_e2","res_e3","res_e4") 

 epsilon_a 	<- mxMatrix(type = "Diag", nrow=nv, ncol=nv, labels=res_a_labs, free=T, values=0.3, name="epsilon_a") 
 epsilon_c 	<- mxMatrix(type = "Diag", nrow=nv, ncol=nv, labels=res_c_labs, free=T, values=0.2, name="epsilon_c") 
 epsilon_e 	<- mxMatrix(type = "Diag", nrow=nv, ncol=nv, labels=res_e_labs, free=T, values=0.1, name="epsilon_e") 

# Expected covariance for mutiple indicator DOC model
 # ? * (I-?)~ * ? * (I-?)~' * ?' + ?
   # ? = lamba - factor loadings
   # I = Identity matrix
   # ? = causal parameters
   # ? = epsilon - measurement error/residuals

 A 	<- mxAlgebra( lambda %&% (solve(I-beta) %&% psi_a) + epsilon_a, name = "A")
 C 	<- mxAlgebra( lambda %&% (solve(I-beta) %&% psi_c) + epsilon_c, name = "C")
 E 	<- mxAlgebra( lambda %&% (solve(I-beta) %&% psi_e) + epsilon_e, name = "E")
 
 covMZ    	<- mxAlgebra( expression= rbind( cbind(A+C+E, A+C), cbind(A+C, A+C+E)), name="expCovMZ" )										 
 covDZ    	<- mxAlgebra( expression= rbind( cbind(A+C+E, 0.5%x%A+C), cbind(0.5%x%A+C, A+C+E)), name="expCovDZ" )

# Calculators
# Standardized variance Components
 V      	<- mxAlgebra( expression= A+C+E, name="V" )
 rowVC      <- rep('VC', nv)
 colVC      <- rep(c('A','C','E','SA','SC','SE'),each=nv)
 VC      	<- mxAlgebra( expression=cbind(A,C,E,A/V,C/V,E/V), dimnames=list(rowVC,colVC), name="VC")

# Standardise factor-specific ACE influences - psi_a,c and e
 Vc      	<- mxAlgebra( expression= psi_a + psi_c + psi_e, name="Vc" )
 h2c      	<- mxAlgebra( expression= psi_a/Vc, name="h2c" )
 c2c      	<- mxAlgebra( expression= psi_c/Vc, name="c2c" )
 e2c      	<- mxAlgebra( expression= psi_e/Vc, name="e2c" )

# Standardise factor-specific ACE influences - psi_a,c and e
 Acc		<- mxAlgebra( (solve(I-beta) %&% psi_a), name = "Acc")
 Ccc		<- mxAlgebra( (solve(I-beta) %&% psi_c), name = "Ccc")
 Ecc		<- mxAlgebra( (solve(I-beta) %&% psi_e), name = "Ecc")
 Vcc      	<- mxAlgebra( expression= Acc + Ccc + Ecc, name="Vcc" )

Stcp1on2	<-mxAlgebra( expression= (beta[2,1]* sqrt(Vcc[1,1]))/sqrt(Vcc[2,2]), name="Stand_1on2")
Stcp2on1	<-mxAlgebra( expression= (beta[1,2]* sqrt(Vcc[2,2]))/sqrt(Vcc[1,1]), name="Stand_2on1")

# Latent true score means for autoregression component
 # Mean matrix = µ
 # Factor mean = ((I-?)~ * µ)'   Adjusts each mean for contribution from causal pathway
 # Manifest mean = (? * ( (I-?)~ * µ ))'  
 
 Mean 	<- mxMatrix( type="Full", nrow=2, ncol=1, free=c(T,T), labels=c( "m1","m2"), values = c(0.01,0.30), name = "Mean" ) # µ
 FacMean  	<- mxAlgebra( solve(I-beta)  %*% Mean, name="FacMean" ) 
 ManMean  	<- mxAlgebra( t(lambda %*% FacMean ), name="ManMean" ) 	
 expMean  	<- mxAlgebra( expression= cbind(ManMean,ManMean), name="expMean" ) 
  
# Data Objects for Multiple Groups
 dataMZ   	<- mxData( observed=mzData, type="raw" )
 dataDZ   	<- mxData( observed=dzData, type="raw" )

# Objective Objects for Multiple Groups
 # Define how the model expectations are calculated 
 # The mxExpectationNormal function uses the algebra defined by the 'covariance' and 'means' arguments to define 
 # the expected covariance and means under the assumption of multivariate normality
 objMZ    	<- mxExpectationNormal( covariance="expCovMZ", means="expMean", dimnames=selVars )
 objDZ    	<- mxExpectationNormal( covariance="expCovDZ", means="expMean", dimnames=selVars )

# Function to compute -2*(log likelihood) of the data given the current values of the free parameters and the expectation function
 fiML    	<- mxFitFunctionML()
 
# ci 		<- mxCI(c("CP2toCP1","CP1toCP2","h2c","c2c","e2c",Acc,Ccc,Ecc,Vcc,Stand_1on2,Stand_2on1))
 ci 		<- mxCI(c("CP2toCP1","CP1toCP2","h2c","c2c","e2c","Acc","Ccc","Ecc","Vcc","Stand_1on2","Stand_2on1"))

# Final steps
 pars    	<- list( psi_a,psi_c,psi_e, I,beta,lambda, epsilon_a,epsilon_c,epsilon_e,A,C,E,Mean,FacMean,ManMean,expMean,V,VC,Vc,h2c,c2c,e2c,Acc,Ccc,Ecc,Vcc,Stcp1on2,Stcp2on1,ci) # 
 modelMZ 	<- mxModel( pars, dataMZ, covMZ, objMZ, fiML, name="MZ" ) # Function to create MxModel objects 
 modelDZ 	<- mxModel( pars, dataDZ, covDZ, objDZ, fiML, name="DZ" ) # Function to create MxModel objects 
 obj  	<- mxFitFunctionMultigroup(c("MZ","DZ"))				  # Function used to fit a multiple group model
 ACE1_model	<- mxModel( "ACE1", pars, modelMZ, modelDZ, obj )#; omxGetParameters(ACE_model)

 # mxCheckIdentification(ACE1_model, details=TRUE)
 		
 # --------------------------------------------------------------------------------------------------------------------------------
# 6 RUN ACE with causal paths Model

ACE1_Fit	<-mxRun(ACE1_model, intervals=T)
(ACE1_Summ	<-summary(ACE1_Fit))

mxEval(MZ.expCovMZ, ACE1_Fit)
mxEval(DZ.expCovDZ, ACE1_Fit)

mxEval(MZ.psi_a, ACE1_Fit)
mxEval(MZ.psi_c, ACE1_Fit)
mxEval(MZ.psi_e, ACE1_Fit)

mxEval(MZ.beta, ACE1_Fit)
mxEval(MZ.lambda, ACE1_Fit)
mxEval(MZ.epsilon_a, ACE1_Fit)
mxEval(MZ.epsilon_c, ACE1_Fit)
mxEval(MZ.epsilon_e, ACE1_Fit)

mxEval(MZ.A, ACE1_Fit)
mxEval(MZ.C, ACE1_Fit)
mxEval(MZ.E, ACE1_Fit)
mxEval(MZ.V, ACE1_Fit)

mxEval(MZ.VC, ACE1_Fit)
mxEval(MZ.Vc, ACE1_Fit)
mxEval(MZ.h2c, ACE1_Fit)
mxEval(MZ.c2c, ACE1_Fit)
mxEval(MZ.e2c, ACE1_Fit)
mxEval(MZ.Acc, ACE1_Fit)
mxEval(MZ.Ccc, ACE1_Fit)
mxEval(MZ.Ecc, ACE1_Fit)
mxEval(MZ.Stand_1on2, ACE1_Fit)


########## Model comparisons ##########

# 0. ACE0 - No correlations, no causation, drop shared ACE components 
 ACE01 	<- mxModel( ACE1_Fit,name="ACE01")
 ACE01 	<- omxSetParameters( ACE01, label=c("a21","c21","e21"),free=F,values=0)   
 summary( ACE01_Fit	<- mxTryHard( ACE01, extraTries=10, greenOK=TRUE,checkHess=FALSE ))
 summary( ACE01_Fit )
 
# 2. AE - correlated A & E risks on common pathways only, no causation, drop shared environmental residuals 
 AE1 	<- mxModel( ACE1_Fit,name="AE1")
 AE1 	<- omxSetParameters( AE1, label=c("c11","c21","c22"),free=F,values=0)   
 AE1 	<- omxSetParameters( AE1, label=res_c_labs,free=F,values=0)  		
 summary( AE1_Fit	<- mxTryHard( AE1, extraTries=10, greenOK=TRUE,checkHess=FALSE , intervals=))
 summary( AE1_Fit )
 
# 3. CE - correated C & E risks on common pathways only, no causation, drop genetic residuals 
 CE1 	<- mxModel( ACE1_Fit,name="CE1")
 CE1 	<- omxSetParameters( CE1, label=c("a11","a21","a22"),free=F,values=0)   
 CE1 	<- omxSetParameters( CE1, label=res_a_labs,free=F,values=0)  		
 summary( CE1_Fit	<- mxTryHard( CE1, extraTries=10, greenOK=TRUE,checkHess=FALSE ))
 summary( CE1_Fit )

# 3.b. E - correated E risks on common pathways, no causation, drop genetic and shared environmental residuals  
 E1 	<- mxModel( ACE1_Fit,name="E1")
 E1 	<- omxSetParameters( E1, label=c("a11","a21","a22"),free=F,values=0)   
 E1 	<- omxSetParameters( E1, label=c("c11","c21","c22"),free=F,values=0)   
 E1 	<- omxSetParameters( E1, label=res_a_labs,free=F,values=0)  		
 E1 	<- omxSetParameters( E1, label=res_c_labs,free=F,values=0)  		
 summary( E1_Fit	<- mxTryHard( E1, extraTries=10, greenOK=TRUE,checkHess=FALSE ))
 summary( E1_Fit )
  
# 4. Compare ACE, AE, CE and E correlated risks models (no causation) 
 subs 	<- c( AE1_Fit, CE1_Fit, E1_Fit, ACE01_Fit )
 comps	<- mxCompare( ACE1_Fit,subs ); comps
 
# 5. ACE - no correlated A, C and E risks, causation from CP1 to CP2
 CP1_to_CP2_1	<- mxModel( ACE1_Fit,name="ACE_CP1toCP2_1")
 CP1_to_CP2_1 	<- omxSetParameters( CP1_to_CP2_1, label=c("a21","c21","e21"),free=FALSE,values=0)  
 CP1_to_CP2_1 	<- omxSetParameters( CP1_to_CP2_1, label="CP1toCP2",free=TRUE,values=0.3)  
 summary( CP1_to_CP2_1_Fit	<- mxTryHard( CP1_to_CP2_1, extraTries=10, greenOK=TRUE,checkHess=FALSE, intervals=T ))
 summary( CP1_to_CP2_1_Fit )
  
mxCompare(ACE1_Fit,1CP1_to_CP2_1_Fit)

# 6. ACE - no correlated A, C and E risks, causation from CP2 to CP1
 CP2_to_CP1_1 	<- mxModel( ACE1_Fit,name="ACE_CP2toCP1_1")
 CP2_to_CP1_1 	<- omxSetParameters( CP2_to_CP1_1, label=c("a21","c21","e21"),free=FALSE,values=0)  
 CP2_to_CP1_1 	<- omxSetParameters( CP2_to_CP1_1, label="CP2toCP1",free=TRUE,values=0.3)  
 summary( CP2_to_CP1_1_Fit	<- mxTryHard( CP2_to_CP1_1, extraTries=10, greenOK=TRUE,checkHess=FALSE, intervals=T ))
 summary( CP2_to_CP1_1_Fit )
    
# 7. ACE - no correlated A, C and E risks, causation from CP2 to CP1 and causation from CP1 to CP2
 CP12_to_CP12_1 	<- mxModel( ACE1_Fit,name="ACE_CP12toCP12_1")
 CP12_to_CP12_1 	<- omxSetParameters( CP12_to_CP12_1, label=c("a21","c21","e21"),free=FALSE,values=0)  
 CP12_to_CP12_1 	<- omxSetParameters( CP12_to_CP12_1, label=c("CP1toCP2","CP2toCP1"),free=TRUE,values=0.3)  
 summary( CP12_to_CP12_1_Fit	<- mxTryHard( CP12_to_CP12_1, extraTries=10, greenOK=TRUE,checkHess=FALSE, intervals=T ))
 summary( CP12_to_CP12_1_Fit )
  mxCompare(ACE1_Fit,CP12_to_CP12_1_Fit)


# 8. Compare the fit of the causal models against the null
 subs 		<- c( ACE01_Fit, CP1_to_CP2_1_Fit, CP2_to_CP1_1_Fit, CP12_to_CP12_1_Fit)
 comps		<- mxCompare( ACE1_Fit, subs ); comps

mxCompare(CP12_to_CP12_1_Fit, CP1_to_CP2_1_Fit)



#****************************************************************************************************************************
# __(5a)_____________________________________________________________________________________________________________________
# Mendelian Randomisation Direction of Causation (MRDoC) MODEL for latent constructs
# We specify Specific effects on the latent factors(Acsp, Ccsp and Ecsp) and add causal paths:
# Causal paths specified between Phenotypic Factors: F1 (PRS) > F2(PARCA) > F3(TR) & F1>F3;
# Asp, Csp and Esp in the bottom with constraints to Identify the model on top
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1), 
#_____________________________________________________________________________________________________________________________

nv		<- 5			# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 3				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

Groups	<-c("mz", "dz")
Vars		<- c("PRSEA","Parca3","Parca4","TalkRhym3","TalkRhym4")
selVars	<- c("PRSEA1","Parca31","Parca41","TalkRhym31","TalkRhym41",
			"PRSEA2","Parca32","Parca42","TalkRhym32","TalkRhym42")

mzData	<- subset(TWINdata2, zyg1%in%c(1)|zyg2%in%c(1), selVars)
dzData	<- subset(TWINdata2, zyg1%in%c(2)|zyg2%in%c(2), selVars)

psych::describe(mzData)
psych::describe(dzData)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabs	<- paste("m",1:nv,sep=""))
(Stmean	<- colMeans(mzData[,1:nv],na.rm=TRUE))
(Stsd 	<- sapply(mzData[,1:nv],sd, na.rm=TRUE))
(PatM		<- c(TRUE,TRUE,TRUE,TRUE,TRUE))

# Create Labels for Diagonal Matrices
# To identify this model we equate the sp effects of the 2 indicators per factor to be equal)
(LabEs	<- c('es1','es2','es2','es4','es4'))
(LabAs	<- c('as1','as2','as2','as4','as4'))
(LabCs	<- c('cs1','cs2','cs2','cs4','cs4'))

PatSpe	<- c(F,TRUE,TRUE,TRUE,TRUE)
PatSpac	<- c(F,T,TRUE,TRUE,TRUE)
StSpa		<- c(0,.5,.5,.5,.5)
StSpc		<- c(0,.5,.5,.5,.5)
StSpe		<- c(0,.5,.5,.5,.5)

# all 1st loadings fixed to 1
PatFl		<- c(F,F,F,F,F,			
		     F,F,T,F,F,
		     F,F,F,F,T)

StFl		<- c(1,0,0,0,0,
		     0,1,.5,0,0,
		     0,0,0,1,.5)

LabFl		<- c('l1',NA,NA,NA,NA,
	 	     NA,'l2','l3',NA,NA,
	 	     NA,NA,NA,'l4','l5')

PatPhC	<- c(F,TRUE,TRUE,
		     F,F,TRUE,
		     F,F,F)

StPhC		<- c(0,.3,.3,
		     0,0,.3,
		     0,0,0)

LabPhC	<- c(NA,'c1on2','c1on3',
		     NA,NA,'c2on3',	
		     NA,NA,NA)	 

#______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
#______________________________________________________________________________________________________

Means		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmean,Stmean), labels=c(mLabs,mLabs), name="expMean") 

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTw	<-mxAlgebra(I2%x%FactL, name="FactLTw")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 3 latent variables  
PhCaus	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhC, labels=LabPhC, name="PhC" )

# Define the matrix to hold the A and C effects: Specific 
PathsAs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpa, labels=LabAs, name="as" )
PathsCs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpc, labels=LabCs, name="cs" )
PathsEs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpe, labels=LabEs, name="es" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components
covAs		<-mxAlgebra( expression= as %*% t(as), name="As" )
covCs		<-mxAlgebra( expression= cs %*% t(cs), name="Cs" )
covEs		<-mxAlgebra( expression= es %*% t(es), name="Es" )
covPs		<-mxAlgebra( expression= As+Cs+Es, name="Vs" )

# Define the matrices to hold the A and C effects: Common 
PathsAcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.5, labels=c("ac22","ac32","ac33"), name="ac" ) # Component paths for factors 2 and 3
PathsCcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("cc22","cc32","cc33"), name="cc" )
PathsEcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,F,T), values=c(.8,0,.8), labels=c("ec22","ec32","ec33"), name="ec" )
PathsP11	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11", name="pc" ) # SD path for factor 1 (the PRS factor)
Ze11		<-mxMatrix(type="Zero",	nrow=1, ncol=1, free=F, name="Z11" )  

Ze21		<-mxMatrix(type="Zero",	nrow=2, ncol=1, free=F, name="Z21" )
Ze12		<-mxMatrix(type="Zero",	nrow=1, ncol=2, free=F, name="Z12" )
covAcsub	<-mxAlgebra( expression= ac %*% t(ac), name="Acsub" )
covCcsub	<-mxAlgebra( expression= cc %*% t(cc), name="Ccsub" )
covEcsub	<-mxAlgebra( expression= ec %*% t(ec), name="Ecsub" )
covPcsub	<-mxAlgebra( expression= Acsub+Ccsub+Ecsub, name="Vcsub" ) #Matrix for the total variance of factors 2 and 3 (i.e. X and Y)

covPc11	<-mxAlgebra( expression= pc %*% t(pc), name="Pc11" ) # variance for factor 1 (the PRS factor), I specify this separately as I do not want to resolve its variance into ACE components

covPc		<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Vcsub)), name="Vc") #I combine the PRS variance with the var-cov matrix of the other two factors.
covPcMz	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Acsub+Ccsub)), name="Vcmz") #I specify the MZ between-twin covariance - excluding E parameters
covPcDz	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z21), rbind(Z12, .5%x%Acsub+Ccsub)), name="Vcdz") #I specify the DZ between-twin covariance - specifying half of A and excluding E

# Generate Covariance of Latent factor model Including Causal Paths between factors
Id3		<-mxMatrix(type="Iden",	nrow=3, ncol=3, free=F, name="I3" )
covFVc	<-mxAlgebra( expression= solve(I3-PhC) %&% Vc, name ="FVc")
covFcMz	<-mxAlgebra( expression= solve(I3-PhC) %&% Vcmz, name ="Fcmz")
covFcDz	<-mxAlgebra( expression= solve(I3-PhC) %&% Vcdz, name ="Fcdz")

FcovMZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcmz), cbind(Fcmz, FVc))) , name="expFCovMZ" )
FcovDZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcdz), cbind(Fcdz, FVc))) , name="expFCovDZ" )

SpcovMZ	<-mxAlgebra( expression= rbind (cbind(Vs, As+Cs), cbind(As+Cs, Vs)) , name="expSpCovMZ" )
SpcovDZ	<-mxAlgebra( expression= rbind (cbind(Vs, .5%x%As+Cs), cbind(.5%x%As+Cs, Vs)) , name="expSpCovDZ" )

TOTcovMZ	<-mxAlgebra( expression= expFCovMZ + expSpCovMZ , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= expFCovDZ + expSpCovDZ , name="TOTexpCovDZ" )
# *******************************************************************************************************
# Calculator

# Standardize the causal effects
Stcp1on2	<-mxAlgebra( expression= (PhC[2,1]* sqrt(FVc[1,1]))/sqrt(FVc[2,2]) , name="Stand_1on2" )
Stcp1on3	<-mxAlgebra( expression= (PhC[3,1]* sqrt(FVc[1,1]))/sqrt(FVc[3,3]) , name="Stand_1on3" )
Stcp2on3	<-mxAlgebra( expression= (PhC[3,2]* sqrt(FVc[2,2]))/sqrt(FVc[3,3]) , name="Stand_2on3" )

# Standardize the Total var/covariances matrices of the observed variables
Id10		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I10" )
Rfactmz	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovMZ)) %&% TOTexpCovMZ, name="FactcorMZ" )
Rfactdz	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovDZ)) %&% TOTexpCovDZ, name="FactcorDZ" )

# Estimate the Common Effects adjusting for causal paths
covAc		<-mxAlgebra( expression= (2 %x% (Fcmz-Fcdz)), name="Ac" )
covCc		<-mxAlgebra( expression= ((2 %x% Fcdz) - Fcmz), name="Cc" )
covEc		<-mxAlgebra( expression= (FVc-Fcmz), name="Ec" )

# Standardize the Total var/covariances matrices of the latent factors without the causal path
RfactPsub	<-mxAlgebra( expression= solve(sqrt(I2*Vcsub)) %&% Vcsub, name="FactcorPsub" ) 

# Phenotypic, A, C and E correlations	
RfactAc	<-mxAlgebra( expression= solve(sqrt(I2*Acsub)) %&% Acsub, name="Ra" )
RfactCc	<-mxAlgebra( expression= solve(sqrt(I2*Ccsub)) %&% Ccsub, name="Rc" )
RfactEc	<-mxAlgebra( expression= solve(sqrt(I2*Ecsub)) %&% Ecsub, name="Re" )
RfactP	<-mxAlgebra( expression= solve(sqrt(I3*FVc)) %&% FVc, name="Rph" )

# Standardize the Common Effects
covFVc22	<-mxAlgebra( expression= FVc[2:3,2:3], name ="FVc22")
stcovAc	<-mxAlgebra( expression= Acsub/FVc22, name="stAc" )
stcovCc	<-mxAlgebra( expression= Ccsub/FVc22, name="stCc" )
stcovEc	<-mxAlgebra( expression= Ecsub/FVc22, name="stEc" )

# Standardised path estimates
StpathAc	<-mxAlgebra( expression= (sqrt(stAc)), name="stpac" )
StpathCc	<-mxAlgebra( expression= (sqrt(stCc)), name="stpcc" )
StpathEc	<-mxAlgebra( expression= (sqrt(stEc)), name="stpec" )

# Algebra to compute Rph-A, Rph-C and Rph-E
RphA		<-mxAlgebra(expression=sqrt(stAc[1,1])*Ra[2,1]*sqrt(stAc[2,2]), name = 'Rpha')
RphC		<-mxAlgebra(expression=sqrt(stCc[1,1])*Rc[2,1]*sqrt(stCc[2,2]), name = 'Rphc')
RphE		<-mxAlgebra(expression=sqrt(stEc[1,1])*Re[2,1]*sqrt(stEc[2,2]), name = 'Rphe')

# Standardize the Specific Effects
stcovAs	<-mxAlgebra( expression= sqrt(As/( (FactL %&% FVc) +Vs)), name="stAs" )
stcovCs	<-mxAlgebra( expression= sqrt(Cs/( (FactL %&% FVc) +Vs)), name="stCs" )
stcovEs	<-mxAlgebra( expression= sqrt(Es/( (FactL %&% FVc) +Vs)), name="stEs" )

# Standardized Factor Loadings
StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% FVc / TOTexpCovMZ[1:5,1:5])) , name="StandFact" )

# *******************************************************************************************************

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expMean", dimnames=selVars)
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expMean", dimnames=selVars)

fitFunction <- mxFitFunctionML()
 
# Combine Groups
pars1		<-list(Means,Load,LoadTw,PhCaus,PathsAs,PathsCs,PathsEs,covAs,covCs,covEs,covPs,Id10,Id2,Id3)
pars2		<-list(PathsAcsub,PathsCcsub,PathsEcsub,PathsP11,Ze21,Ze12,Ze11,covAcsub,covCcsub,covEcsub,covPcsub,covPc11,covPc,covPcMz,covPcDz,covFVc,covFVc22,covFcMz,covFcDz,covAc,covCc,covEc)
parsst	<-list(stcovAs, stcovCs, stcovEs, stcovAc, stcovCc, stcovEc, RfactAc, RfactCc, RfactEc, RfactP, RfactPsub,RphA,RphC,RphE,StpathAc,StpathCc,StpathEc)
parsmed	<-list(Stcp1on2, Stcp1on3, Stcp2on3)
modelMZ	<-mxModel(pars1, pars2, parsmed, FcovMZ, SpcovMZ, TOTcovMZ, dataMZ, objMZ, Rfactmz, parsst, fitFunction, StFL, name="MZ" )
modelDZ	<-mxModel(pars1, pars2, FcovDZ, SpcovDZ, TOTcovDZ, dataDZ, objDZ, Rfactdz, fitFunction, name="DZ" )
minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
cistFL	<-mxCI (c ('MZ.StandFact','MZ.Stand_1on2','MZ.Stand_1on3','MZ.Stand_2on3','MZ.PhC'))
cistVs	<-mxCI (c ('MZ.stAs[2,2]','MZ.stAs[3,3]','MZ.stAs[4,4]','MZ.stAs[5,5]',
				'MZ.stCs[2,2]','MZ.stCs[3,3]','MZ.stCs[4,4]','MZ.stCs[5,5]',
				'MZ.stEs[1,1]','MZ.stEs[2,2]','MZ.stEs[3,3]','MZ.stEs[4,4]','MZ.stEs[5,5]') ) 	# standardized var comp from specific Factors
cistVc	<-mxCI (c ('MZ.stAc[1,1]','MZ.stAc[2,1]','MZ.stAc[2,2]','MZ.stCc[1,1]','MZ.stCc[2,1]','MZ.stCc[2,2]','MZ.stEc[1,1]','MZ.stEc[2,1]','MZ.stEc[2,2]') ) 	# standardized var comp for ACE on latent Factors
cistRc	<-mxCI (c ('MZ.Rpha','MZ.Rphc','MZ.Rphe','MZ.Ra','MZ.Rc','MZ.Re','MZ.stpac','MZ.stpcc','MZ.stpec') ) 	
ACEMs1Model	<-mxModel("aceMs1", pars1, pars2, modelMZ, modelDZ, minus2ll, obj, cistFL, cistVs, cistVc, cistRc) 

# --------------------------------------------------------------------------------------------------------------------------------
# 7a RUN ACEMs Factor Model with phenotypic causal mediation paths by Zygosity

ACEMs1Fit	<-mxRun(ACEMs1Model, intervals=T)
(ACEMs1Summ	<-summary(ACEMs1Fit, verbose=T))

mxEval(MZ.Stand_1on2[1,1], ACEMs1Fit)
mxEval(MZ.Stand_1on3[1,1], ACEMs1Fit)
mxEval(MZ.Stand_2on3[1,1], ACEMs1Fit)

######################################################################
# CHECK IDENTIFICATION
######################################################################
# 7bi: Full Model - UNIDENTIFIED

mgFitfun1 	<- mxFitFunctionMultigroup(c("MZ","DZ"))
mgModel1 	<- mxModel("MRDoC1", modelMZ, modelDZ, mgFitfun1)
mgFit1 	<- mxRun(mgModel1)
(mgSumm1	<-summary(mgFit1))

# This checks identification specifically
(Id		<-mxCheckIdentification(mgFit1))



#****************************************************************************************************************************
# __(5b) ____________________________________________________________________________________________________________________
# Mendelian Randomisation Direction of Causation (MRDoC) MODEL with observed variables as single indicators 
# We specify Specific effects on the latent factors(Acsp, Ccsp and Ecsp) and add causal paths:
# Causal paths specified between Phenotypic Factors: F1 (PRS) >F2 (PARCA) >F3 (TR) & F1>F3;
# Asp, Csp and Esp in the bottom with constraints to Identify the model on top
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1), 
# this because applying a constraint on the factor variances of 1 is problematic especially when we model the causal paths.
# To identify the model we constrain Asp, Csp and Esp variance components loading on variables 1 and 10 to zero
#_____________________________________________________________________________________________________________________________

nv		<- 3			# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 3				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv

Groups	<-c("mz", "dz")
Vars		<- c("PRSEA","Parca3","TalkRhym3")
selVars	<- c("PRSEA1","Parca31","TalkRhym31",
			"PRSEA2","Parca32","TalkRhym32")

#Vars		<- c("PRSEA","Parca3","Parca4","TalkRhym34","PlayBook34","PlayGames34")
#selVars	<- c("PRSEA1","Parca31","Parca41","TalkRhym341","PlayBook341","PlayGames341",
#			"PRSEA2","Parca32","Parca42","TalkRhym342","PlayBook342","PlayGames342")

mzData	<- subset(TWINdata2, zyg1%in%c(1)|zyg2%in%c(1), selVars)
dzData	<- subset(TWINdata2, zyg1%in%c(2)|zyg2%in%c(2), selVars)

psych::describe(mzData)
psych::describe(dzData)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabs	<- paste("m",1:nv,sep=""))
(Stmean	<- colMeans(mzData[,1:nv],na.rm=TRUE))
(Stsd 	<- sapply(mzData[,1:nv],sd, na.rm=TRUE))
(PatM		<- c(TRUE,TRUE,TRUE))

# Create Labels for Diagonal Matrices
# To identify this model we equate the sp effects of the 2 indicators per factor to be equal)
(LabEs	<- c('es1','es2','es3'))
(LabAs	<- c('as1','as2','as3'))
(LabCs	<- c('cs1','cs2','cs3'))

PatSpe	<- c(F)
PatSpac	<- c(F)
StSpa		<- c(0)
StSpc		<- c(0)
StSpe		<- c(0)

# all 1st loadings fixed to 1
PatFl		<- c(F,F,F,			
		     F,F,F,
		     F,F,F)

StFl		<- c(1,0,0,
		     0,1,0,
		     0,0,1)

LabFl		<- c('l1',NA,NA,
	 	     NA,'l2',NA,
	 	     NA,NA,'l3')

PatPhC	<- c(F,TRUE,TRUE,
		     F,F,TRUE,
		     F,F,F)

StPhC		<- c(0,.3,.3,
		     0,0,.3,
		     0,0,0)

LabPhC	<- c(NA,'c1on2','c1on3',
		     NA,NA,'c2on3',	
		     NA,NA,NA)	 

#______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
#______________________________________________________________________________________________________

Means		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmean,Stmean), labels=c(mLabs,mLabs), name="expMean") 

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTw	<-mxAlgebra(I2%x%FactL, name="FactLTw")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 3 latent variables  
PhCaus	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhC, labels=LabPhC, name="PhC" )

# Define the matrix to hold the A and C effects: Specific 
PathsAs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpa, labels=LabAs, name="as" )
PathsCs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpc, labels=LabCs, name="cs" )
PathsEs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpe, labels=LabEs, name="es" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components
covAs		<-mxAlgebra( expression= as %*% t(as), name="As" )
covCs		<-mxAlgebra( expression= cs %*% t(cs), name="Cs" )
covEs		<-mxAlgebra( expression= es %*% t(es), name="Es" )
covPs		<-mxAlgebra( expression= As+Cs+Es, name="Vs" )

# Define the matrices to hold the A and C effects: Common 
PathsAcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.5, labels=c("ac22","ac32","ac33"), name="ac" ) # Component paths for factors 2 and 3
PathsCcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("cc22","cc32","cc33"), name="cc" )
PathsEcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,F,T), values=c(.8,0,.8), labels=c("ec22","ec32","ec33"), name="ec" )
PathsP11	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11", name="pc" ) # SD path for factor 1 (the PRS factor)
Ze11		<-mxMatrix(type="Zero",	nrow=1, ncol=1, free=F, name="Z11" )  

Ze21		<-mxMatrix(type="Zero",	nrow=2, ncol=1, free=F, name="Z21" )
Ze12		<-mxMatrix(type="Zero",	nrow=1, ncol=2, free=F, name="Z12" )
covAcsub	<-mxAlgebra( expression= ac %*% t(ac), name="Acsub" )
covCcsub	<-mxAlgebra( expression= cc %*% t(cc), name="Ccsub" )
covEcsub	<-mxAlgebra( expression= ec %*% t(ec), name="Ecsub" )
covPcsub	<-mxAlgebra( expression= Acsub+Ccsub+Ecsub, name="Vcsub" ) #Matrix for the total variance of factors 2 and 3 (i.e. X and Y)

covPc11	<-mxAlgebra( expression= pc %*% t(pc), name="Pc11" ) # variance for factor 1 (the PRS factor), I specify this separately as I do not want to resolve its variance into ACE components

covPc		<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Vcsub)), name="Vc") #I combine the PRS variance with the var-cov matrix of the other two factors.
covPcMz	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Acsub+Ccsub)), name="Vcmz") #I specify the MZ between-twin covariance - excluding E parameters
covPcDz	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z21), rbind(Z12, .5%x%Acsub+Ccsub)), name="Vcdz") #I specify the DZ between-twin covariance - specifying half of A and excluding E

# Generate Covariance of Latent factor model Including Causal Paths between factors
Id3		<-mxMatrix(type="Iden",	nrow=3, ncol=3, free=F, name="I3" )
covFVc	<-mxAlgebra( expression= solve(I3-PhC) %&% Vc, name ="FVc")
covFcMz	<-mxAlgebra( expression= solve(I3-PhC) %&% Vcmz, name ="Fcmz")
covFcDz	<-mxAlgebra( expression= solve(I3-PhC) %&% Vcdz, name ="Fcdz")

FcovMZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcmz), cbind(Fcmz, FVc))) , name="expFCovMZ" )
FcovDZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcdz), cbind(Fcdz, FVc))) , name="expFCovDZ" )

SpcovMZ	<-mxAlgebra( expression= rbind (cbind(Vs, As+Cs), cbind(As+Cs, Vs)) , name="expSpCovMZ" )
SpcovDZ	<-mxAlgebra( expression= rbind (cbind(Vs, .5%x%As+Cs), cbind(.5%x%As+Cs, Vs)) , name="expSpCovDZ" )

TOTcovMZ	<-mxAlgebra( expression= expFCovMZ + expSpCovMZ , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= expFCovDZ + expSpCovDZ , name="TOTexpCovDZ" )
# *******************************************************************************************************
# Calculator

# Standardize the causal effects
Stcp1on2	<-mxAlgebra( expression= (PhC[2,1]* sqrt(FVc[1,1]))/sqrt(FVc[2,2]) , name="Stand_1on2" )
Stcp1on3	<-mxAlgebra( expression= (PhC[3,1]* sqrt(FVc[1,1]))/sqrt(FVc[3,3]) , name="Stand_1on3" )
Stcp2on3	<-mxAlgebra( expression= (PhC[3,2]* sqrt(FVc[2,2]))/sqrt(FVc[3,3]) , name="Stand_2on3" )

# Standardize the Total var/covariances matrices of the observed variables
Id10		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I10" )
Rfactmz	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovMZ)) %&% TOTexpCovMZ, name="FactcorMZ" )
Rfactdz	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovDZ)) %&% TOTexpCovDZ, name="FactcorDZ" )

# Estimate the Common Effects adjusting for causal paths
covAc		<-mxAlgebra( expression= (2 %x% (Fcmz-Fcdz)), name="Ac" )
covCc		<-mxAlgebra( expression= ((2 %x% Fcdz) - Fcmz), name="Cc" )
covEc		<-mxAlgebra( expression= (FVc-Fcmz), name="Ec" )

# Standardize the Total var/covariances matrices of the latent factors without the causal path
RfactPsub	<-mxAlgebra( expression= solve(sqrt(I2*Vcsub)) %&% Vcsub, name="FactcorPsub" ) 

# Phenotypic, A, C and E correlations	
RfactAc	<-mxAlgebra( expression= solve(sqrt(I2*Acsub)) %&% Acsub, name="Ra" )
RfactCc	<-mxAlgebra( expression= solve(sqrt(I2*Ccsub)) %&% Ccsub, name="Rc" )
RfactEc	<-mxAlgebra( expression= solve(sqrt(I2*Ecsub)) %&% Ecsub, name="Re" )
RfactP	<-mxAlgebra( expression= solve(sqrt(I3*FVc)) %&% FVc, name="Rph" )

# Standardize the Common Effects
covFVc22	<-mxAlgebra( expression= FVc[2:3,2:3], name ="FVc22")
stcovAc	<-mxAlgebra( expression= Acsub/FVc22, name="stAc" )
stcovCc	<-mxAlgebra( expression= Ccsub/FVc22, name="stCc" )
stcovEc	<-mxAlgebra( expression= Ecsub/FVc22, name="stEc" )

# Standardised path estimates
StpathAc	<-mxAlgebra( expression= (sqrt(stAc)), name="stpac" )
StpathCc	<-mxAlgebra( expression= (sqrt(stCc)), name="stpcc" )
StpathEc	<-mxAlgebra( expression= (sqrt(stEc)), name="stpec" )

# Algebra to compute Rph-A, Rph-C and Rph-E
RphA		<-mxAlgebra(expression=sqrt(stAc[1,1])*Ra[2,1]*sqrt(stAc[2,2]), name = 'Rpha')
RphC		<-mxAlgebra(expression=sqrt(stCc[1,1])*Rc[2,1]*sqrt(stCc[2,2]), name = 'Rphc')
RphE		<-mxAlgebra(expression=sqrt(stEc[1,1])*Re[2,1]*sqrt(stEc[2,2]), name = 'Rphe')

# Standardize the Specific Effects
stcovAs	<-mxAlgebra( expression= sqrt(As/( (FactL %&% FVc) +Vs)), name="stAs" )
stcovCs	<-mxAlgebra( expression= sqrt(Cs/( (FactL %&% FVc) +Vs)), name="stCs" )
stcovEs	<-mxAlgebra( expression= sqrt(Es/( (FactL %&% FVc) +Vs)), name="stEs" )

# Standardized Factor Loadings
StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% FVc / TOTexpCovMZ[1:3,1:3])) , name="StandFact" )

# *******************************************************************************************************

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expMean", dimnames=selVars)
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expMean", dimnames=selVars)

fitFunction <- mxFitFunctionML()
 
# Combine Groups
pars1		<-list(Means,Load,LoadTw,PhCaus,PathsAs,PathsCs,PathsEs,covAs,covCs,covEs,covPs,Id10,Id2,Id3)
pars2		<-list(PathsAcsub,PathsCcsub,PathsEcsub,PathsP11,Ze21,Ze12,Ze11,covAcsub,covCcsub,covEcsub,covPcsub,covPc11,covPc,covPcMz,covPcDz,covFVc,covFVc22,covFcMz,covFcDz,covAc,covCc,covEc)
parsst	<-list(stcovAs, stcovCs, stcovEs, stcovAc, stcovCc, stcovEc, RfactAc, RfactCc, RfactEc, RfactP, RfactPsub,RphA,RphC,RphE,StpathAc,StpathCc,StpathEc)
parsmed	<-list(Stcp1on2, Stcp1on3, Stcp2on3)
modelMZ	<-mxModel(pars1, pars2, parsmed, FcovMZ, SpcovMZ, TOTcovMZ, dataMZ, objMZ, Rfactmz, parsst, fitFunction, StFL, name="MZ" )
modelDZ	<-mxModel(pars1, pars2, FcovDZ, SpcovDZ, TOTcovDZ, dataDZ, objDZ, Rfactdz, fitFunction, name="DZ" )
minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
cistFL	<-mxCI (c ('MZ.StandFact','MZ.Stand_1on2','MZ.Stand_1on3','MZ.Stand_2on3','MZ.PhC'))
cistVc	<-mxCI (c ('MZ.stAc[1,1]','MZ.stAc[2,1]','MZ.stAc[2,2]','MZ.stCc[1,1]','MZ.stCc[2,1]','MZ.stCc[2,2]','MZ.stEc[1,1]','MZ.stEc[2,1]','MZ.stEc[2,2]') ) 	# standardized var comp for ACE on latent Factors
cistRc	<-mxCI (c ('MZ.Rpha','MZ.Rphc','MZ.Rphe','MZ.Ra','MZ.Rc','MZ.Re','MZ.stpac','MZ.stpcc','MZ.stpec') ) 	
ACEMsv1Model	<-mxModel("aceMsv1", pars1, pars2, modelMZ, modelDZ, minus2ll, obj, cistFL, cistVc, cistRc) 

# --------------------------------------------------------------------------------------------------------------------------------
# 7ai RUN ACEMs Factor Model with phenotypic causal mediation paths by Zygosity

ACEMsv1Fit	<-mxRun(ACEMsv1Model, intervals=T)
(ACEMsv1Summ	<-summary(ACEMsv1Fit, verbose=F))

mxEval(MZ.Stand_1on2[1,1], ACEMsv1Fit)
mxEval(MZ.Stand_1on3[1,1], ACEMsv1Fit)
mxEval(MZ.Stand_2on3[1,1], ACEMsv1Fit)

######################################################################
# CHECK IDENTIFICATION
######################################################################
# 7bi: Full Model - UNIDENTIFIED

mgFitfun1 	<- mxFitFunctionMultigroup(c("MZ","DZ"))
mgModel1 	<- mxModel("MRDoC1", modelMZ, modelDZ, mgFitfun1)
mgFit1 	<- mxRun(mgModel1)
(mgSumm1	<-summary(mgFit1))

# This checks identification specifically
(Id		<-mxCheckIdentification(mgFit1))


mxEval(MZ.Tot_eff, ACEMsv1Fit)
mxEval(MZ.Indi_eff, ACEMsv1Fit)

mxEval(MZ.FactcorMZ, ACEMsv1Fit)
mxEval(MZ.FactcorDZ, ACEMsv1Fit)

mxEval(MZ.Acsub, ACEMsv1Fit)
mxEval(MZ.Ccsub, ACEMsv1Fit)
mxEval(MZ.Ecsub, ACEMsv1Fit)
mxEval(MZ.Vcsub, ACEMsv1Fit)
mxEval(MZ.FVc22, ACEMsv1Fit)

mxEval(MZ.Ac, ACEMsv1Fit)
mxEval(MZ.Cc, ACEMsv1Fit)
mxEval(MZ.Ec, ACEMsv1Fit)
mxEval(MZ.FVc, ACEMsv1Fit)

mxEval(MZ.stAc, ACEMsv1Fit)
mxEval(MZ.stCc, ACEMsv1Fit)
mxEval(MZ.stEc, ACEMsv1Fit)

mxEval(MZ.stpac, ACEMsv1Fit)
mxEval(MZ.stpcc, ACEMsv1Fit)
mxEval(MZ.stpec, ACEMsv1Fit)

mxEval(MZ.Ra, ACEMsv1Fit)
mxEval(MZ.Rc, ACEMsv1Fit)
mxEval(MZ.Re, ACEMsv1Fit)
mxEval(MZ.Rph, ACEMsv1Fit)
mxEval(MZ.FactcorPsub, ACEMsv1Fit) 

mxEval(MZ.Rpha, ACEMsv1Fit) 
mxEval(MZ.Rphc, ACEMsv1Fit) 
mxEval(MZ.Rphe, ACEMsv1Fit) 

mxEval(MZ.StandFact, ACEMsv1Fit)
mxEval(MZ.Stand_1on2[1,1], ACEMsv1Fit)
mxEval(MZ.Stand_1on3[1,1], ACEMsv1Fit)



#****************************************************************************************************************************
# __(6a)_____________________________________________________________________________________________________________________
# CROSS-LAG MODEL with var/cov specified with Additive genetic and shared and ind-specific Environ effects;
# And causal paths between observed variables
# For CD and TR
# BURT model
# ___________________________________________________________________________________________________________________________


nv		<- 4				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 2				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nfcor		<- (nfact*(nfact-1)/2)	# number of free elements in a correlation matrix nfact*nfcat
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv-1)/2)		# number of free elements in a correlation matrix nv*nv
Groups	<-c("mz", "dz")

Vars		<- c("Parca3","TalkRhym3","Parca4","TalkRhym4")
selVars	<- c("Parca31","TalkRhym31","Parca41","TalkRhym41",
			"Parca32","TalkRhym32","Parca42","TalkRhym42")

mzData	<- subset(TWINdata2, zyg1%in%c(1)|zyg2%in%c(1), selVars)
dzData	<- subset(TWINdata2, zyg1%in%c(2)|zyg2%in%c(2), selVars)

psych::describe(mzData)
psych::describe(dzData)

# Create labels and starting values

(Stmean	<-colMeans(mzData[,1:nv],na.rm=TRUE))
StM 		<-c(Stmean, Stmean)

(LabM		<- paste("m",1:nv,sep=""))
MLabs		<-c(LabM,LabM) 

# Mean
Mean		<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, labels=MLabs, values=StM, name="ExpMean")

# Matrices to store a, c and e Path Coefficients
pathAw1	<-mxMatrix( type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.7,.7), labels=c("aw111","aw122"), name="a1") 
pathCw1	<-mxMatrix( type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.7,.7), labels=c("cw111","cw122"), name="c1")
pathEw1	<-mxMatrix( type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.7,.7), labels=c("ew111","ew122"), name="e1") 
pathAw2	<-mxMatrix( type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.4,.5), labels=c("ar11","ar22"), name="ar") 
pathCw2	<-mxMatrix( type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.4,.5), labels=c("cr11","cr22"), name="cr")
pathEw2	<-mxMatrix( type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.4,.5), labels=c("er11","er22"), name="er")

# ACE correlations at waves 1 and 2
pathAw1R	<-mxMatrix( type="Stand", nrow=2, ncol=2, free=c(T), values=c(.4), labels=c("raw1"), lbound=-.999, ubound=.999, name="Raw1")
pathCw1R	<-mxMatrix( type="Stand", nrow=2, ncol=2, free=c(T), values=c(.5), labels=c("rcw1"), lbound=-.999, ubound=.999, name="Rcw1") 
pathEw1R	<-mxMatrix( type="Stand", nrow=2, ncol=2, free=c(T), values=c(.5), labels=c("rew1"), lbound=-.999, ubound=.999, name="Rew1") 
pathAw2R	<-mxMatrix( type="Stand", nrow=2, ncol=2, free=c(T), values=c(.4), labels=c("raw2"), lbound=-.999, ubound=.999, name="Rar")
pathCw2R	<-mxMatrix( type="Stand", nrow=2, ncol=2, free=c(T), values=c(.5), labels=c("rcw2"), lbound=-.999, ubound=.999, name="Rcr") 
pathEw2R	<-mxMatrix( type="Stand", nrow=2, ncol=2, free=c(T), values=c(.5), labels=c("rew2"), lbound=-.999, ubound=.999, name="Rer") 

pathME	<-mxMatrix( type="Diag", nrow=nv, ncol=nv, free=c(F,F,F,F), values=c(0,0,0,0), labels=c("me1","me2",NA,NA), name="Me") # Measurement Error for Variable 1 and 2 only -original FR

matZ22	<-mxMatrix( type="Zero", nrow=2, ncol=2, name="Z22" )
matZ42	<-mxMatrix(type="Zero", nrow=nv, ncol=2, name="Z42")
matZ44	<-mxMatrix(type="Zero", nrow=nv, ncol=nv, name="Z44")

# Matrices generated to hold A, C and E computed Variance Components within a person not including causal paths as well as ME variance
CA	<-mxAlgebra( expression=rbind( cbind(a1%&%Raw1,Z22),  cbind(Z22, ar%&%Rar) ) , name="A" )
CC	<-mxAlgebra( expression=rbind( cbind(c1%&%Rcw1,Z22),  cbind(Z22, cr%&%Rcr) ) , name="C" )
CE	<-mxAlgebra( expression=rbind( cbind(e1%&%Rew1,Z22),  cbind(Z22, er%&%Rer) ) , name="E" )
CME	<-mxAlgebra(expression= Me%*%Me, name="ME") 

# Single-headed Causal paths between the 4 variables 
causP	<-mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), labels=c("stabCD","CDonTR","TRonCD","stabTR"), values=c(.6,.1,.1,.6), name="P") 
Caus	<-mxAlgebra(expression= cbind( rbind(Z22, P), Z42), name="SingArInd") 

# To blow it up for a pair
CausS		<-mxAlgebra(expression= rbind( cbind(SingArInd, Z44), cbind(Z44, SingArInd) ), name="SingArPair") 

matI88	<-mxMatrix(type="Iden", nrow=ntv, ncol=ntv, name="I88")
matI22	<-mxMatrix(type="Iden", nrow=2, ncol=2, name="I22")

covA	<-mxAlgebra(expression= (solve(I88-SingArPair)) %&% (rbind( cbind(A,A),cbind(A,A)) ) , 		name="CovA") 
covC	<-mxAlgebra(expression= (solve(I88-SingArPair)) %&% (rbind( cbind(C,C),cbind(C,C)) ) , 		name="CovC") 
covE	<-mxAlgebra(expression= (solve(I88-SingArPair)) %&% (rbind( cbind(E+ME,Z44),cbind(Z44,E+ME)) ) , 	name="CovE") 

covMZ	<-mxAlgebra(expression= (solve(I88-SingArPair)) %&% (rbind( cbind(A+C+E+ME,A+C),cbind(A+C,A+C+E+ME)) ) , 	name="ExpCovMZ") 
covDZ	<-mxAlgebra(expression= (solve(I88-SingArPair)) %&% (rbind( cbind(A+C+E+ME,.5%x%A+C),cbind(.5%x%A+C,A+C+E+ME)) ) , 	name="ExpCovDZ") 

# ----------------------------------------------------------------------------------------------------------------------------------------------
# Calculator in the script to generate standardized estimates (no implications on the model fit)

# Total var/cov of all 4 variables
matI44	<-mxMatrix(type="Iden", nrow=4, ncol=4, name="I44")
covT		<-mxAlgebra(expression= ((solve(I44-SingArInd)) %&% (A+C+E+ME)) , name="V") 		# Total Var/Cov of an indiv (including ME)										
corT		<-mxAlgebra( expression= solve(sqrt(I44*V)) %&% V  , 	name="ReT")				# Implied Cor between the 4 variables
							
# Stand A, C and E (Total)
standA	<-mxAlgebra( expression= ((solve(I44-SingArInd)) %&% A) /V, name="stA" )												
standC	<-mxAlgebra( expression= ((solve(I44-SingArInd)) %&% C) /V, name="stC" )												
standE1	<-mxAlgebra( expression= ((solve(I44-SingArInd)) %&% (E)) /V, name="stEvar12" )
standE2	<-mxAlgebra( expression= ((solve(I44-SingArInd)) %&% (E+ME)) /V, name="stEvar34" )

# Stand A, C and E (Residual, time-2 specific variance in elements 3,3 and 4,4)
standAr	<-mxAlgebra( expression= A/V, name="stAr" )												
standCr	<-mxAlgebra( expression= C/V, name="stCr" )												
standEr	<-mxAlgebra( expression= E/V, name="stEr" )													
standMe	<-mxAlgebra( expression= ME/V, name="stMe" )
standTEr	<-mxAlgebra( expression= E+ME/V, name="stTEr" )

# Stand causal STABILITY paths, square to get variance
StStab1	<-mxAlgebra( expression= ( P[1,1] * sqrt(V[1,1])) /(sqrt(V[3,3]))  , name="stStabCDPath" )
StStab2	<-mxAlgebra( expression= ( P[2,2] * sqrt(V[2,2])) /(sqrt(V[4,4]))  , name="stStabTRPath" )

# Stand causal CROSS-LAGG paths, square to get variance
StCaus1	<-mxAlgebra( expression= ( P[2,1] * sqrt(V[1,1])) /(sqrt(V[4,4])) , name="stCDonTR" )
StCaus2	<-mxAlgebra( expression= ( P[1,2] * sqrt(V[2,2])) /(sqrt(V[3,3])) , name="stTRonCD" )

# ----------------------------------------------------------------------------------------------------------------------------------------------

# Data objects 
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="ExpCovMZ", means="ExpMean", dimnames=selVars)
objDZ		<- mxExpectationNormal( covariance="ExpCovDZ", means="ExpMean", dimnames=selVars)

fitFunction 	<- mxFitFunctionML()
 
# Combine Groups
pars1	   	<- list(pathAw1, pathAw2, pathCw1, pathCw2, pathEw1, pathEw2, pathAw1R, pathAw2R, pathCw1R, pathCw2R, pathEw1R, pathEw2R,
				pathME, CME, CA, CC, CE, causP, Caus, CausS, covA, covC, covE, matZ22, matZ42, matZ44, matI88, matI44, matI22)
pars2	   	<- list(corT, covT, standA, standC,standE1, standE2, StStab1, StStab2, StCaus1, StCaus2, standAr, standCr, standEr, standMe, standTEr)
modelMZ	<- mxModel(pars1, pars2, covMZ, Mean, dataMZ, objMZ, fitFunction, name="MZ" )
modelDZ	<- mxModel(pars1, covDZ, Mean, dataDZ, objDZ, fitFunction, name="DZ" )
minus2ll	<- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<- mxFitFunctionAlgebra("m2LL" )
ciCaus1	<- mxCI (c  ('MZ.stCDonTR[1,1]'))
ciCaus2	<- mxCI (c  ('MZ.stTRonCD[1,1]'))
ciStabRTVb	<- mxCI (c  ('MZ.stStabCDPath[1,1]'))
ciStabAsb	<- mxCI (c  ('MZ.stStabTRPath[1,1]'))
cistA		<- mxCI (c  ('MZ.stA[1,1]', 'MZ.stA[3,3]', 'MZ.stAr[3,3]', 'MZ.stAr[4,4]' ))
cistC		<- mxCI (c  ('MZ.stC[1,1]', 'MZ.stC[3,3]', 'MZ.stCr[3,3]', 'MZ.stCr[4,4]' ))
cistE		<- mxCI (c  ('MZ.stEvar34[1,1]',  'MZ.stEvar34[3,3]', 'MZ.stEvar12[1,1]')) # E total IQ1 and IQ2 and E on IQ1 without ME
cistME	<- mxCI (c  ('MZ.stMe[1,1]','MZ.stMe[2,2]')) # ME on IQ1 
cistE2	<- mxCI (c  ('MZ.stEr[3,3]', 'MZ.stEr[4,4]')) # Eres on IQ2 and AS2
ciRaw1	<- mxCI ('MZ.Raw1[1,2]')
ciRcw1	<- mxCI ('MZ.Rcw1[1,2]')
ciRew1 	<- mxCI ('MZ.Rew1[1,2]')
ciRar 	<- mxCI ('MZ.Rar[2,1]')
ciRcr 	<- mxCI ('MZ.Rcr[2,1]')
ciRer 	<- mxCI ('MZ.Rer[2,1]')

ACEclModel	<- mxModel( "ACEcl", modelMZ, modelDZ, minus2ll, obj, ciCaus1, ciCaus2, ciStabRTVb, ciStabAsb, cistA, cistC, cistE, cistE2, cistME, ciRaw1, ciRcw1, ciRew1, ciRar, ciRcr, ciRer   ) 

#----------------------------------------------------------------------------------------------------------------------------------------------------
# RUN ACEcl model 
ACEclFit		<- mxRun(ACEclModel, intervals=F)
(ACEclSumm		<- summary(ACEclFit,verbose=F))


# OUTPUT specific parameters of interest of sub model

# checks
mxEval(MZ.CovA[1:4,1:4],ACEclFit)
mxEval(MZ.CovC[1:4,1:4],ACEclFit)
mxEval(MZ.CovE[1:4,1:4],ACEclFit)
#mxEval(MZ.CovTEnv[1:4,1:4],ACEclFit)
mxEval(MZ.V,ACEclFit)
mxEval(MZ.ReT,ACEclFit)

mxEval(MZ.Raw1,ACEclFit)
mxEval(MZ.Rcw1,ACEclFit)
mxEval(MZ.Rew1,ACEclFit)
mxEval(MZ.Rar,ACEclFit)
mxEval(MZ.Rcr,ACEclFit)
mxEval(MZ.Rer,ACEclFit)

# Unstandardized x-lag paths
mxEval(MZ.SingArPair,ACEclFit)

# Standardized variance components
# Note that elements 3,3 and 4,4, give the Total F and E variance for the time-two variables (i.e. the transmitted + Residual effects)
mxEval(MZ.stA,ACEclFit)
mxEval(MZ.stC,ACEclFit)
# note that E for var 1 does not include ME
mxEval(MZ.stEvar12,ACEclFit)
mxEval(MZ.stEvar34,ACEclFit)
mxEval(MZ.stMe,ACEclFit)

# Note that elements 3,3 and 4,4, give the Residual A, C and E variance for the time-two variables 
mxEval(MZ.stAr,ACEclFit)
mxEval(MZ.stCr,ACEclFit)
mxEval(MZ.stEr,ACEclFit)

# Standardized M error for 1,1, and 2,2
mxEval(MZ.stMe,ACEclFit)

# Standardized x-lag paths
mxEval(MZ.stStabCDPath,ACEclFit)
mxEval(MZ.stStabTRPath,ACEclFit)
mxEval(MZ.stCDonTR,ACEclFit)
mxEval(MZ.stTRonCD,ACEclFit)


#****************************************************************************************************************************
# __(6b)_____________________________________________________________________________________________________________________
# CROSS-LAG MODEL with var/cov specified with Additive genetic and shared and ind-specific Environ effects
# And causal paths between observed variables
# For CD and TR
# BURT model
# Now including cross-time rA and rC correlations for both the stability and cross-lag paths (Gaussian)
# ___________________________________________________________________________________________________________________________

nv		<- 4				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 2				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nfcor		<- (nfact*(nfact-1)/2)	# number of free elements in a correlation matrix nfact*nfcat
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv-1)/2)		# number of free elements in a correlation matrix nv*nv
Groups	<-c("mz", "dz")

Vars		<- c("Parca3","TalkRhym3","Parca4","TalkRhym4")
selVars	<- c("Parca31","TalkRhym31","Parca41","TalkRhym41",
			"Parca32","TalkRhym32","Parca42","TalkRhym42")

mzData	<- subset(TWINdata2, zyg1%in%c(1)|zyg2%in%c(1), selVars)
dzData	<- subset(TWINdata2, zyg1%in%c(2)|zyg2%in%c(2), selVars)

psych::describe(mzData)
psych::describe(dzData)

# Create labels and starting values

(Stmean	<-colMeans(mzData[,1:nv],na.rm=TRUE))
StM 		<-c(Stmean, Stmean)

(LabM		<- paste("m",1:nv,sep=""))
MLabs		<-c(LabM,LabM) 

# Mean
Mean		<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, labels=MLabs, values=StM, name="ExpMean")

# Matrices to store a, c and e Path Coefficients
pathAw1	<-mxMatrix( type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.7,.7), labels=c("aw111","aw122"), name="a1") 
pathCw1	<-mxMatrix( type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.7,.7), labels=c("cw111","cw122"), name="c1")
pathEw1	<-mxMatrix( type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.7,.7), labels=c("ew111","ew122"), name="e1") 
pathAw2	<-mxMatrix( type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.4,.5), labels=c("ar11","ar22"), name="ar") 
pathCw2	<-mxMatrix( type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.4,.5), labels=c("cr11","cr22"), name="cr")
pathEw2	<-mxMatrix( type="Diag", nrow=2, ncol=2, free=c(T,T), values=c(.4,.5), labels=c("er11","er22"), name="er")

# ACE correlations at waves 1 and 2
pathAw1R	<-mxMatrix( type="Stand", nrow=2, ncol=2, free=c(T), values=c(.4), labels=c("raw1"), lbound=-.999, ubound=.999, name="Raw1")
pathCw1R	<-mxMatrix( type="Stand", nrow=2, ncol=2, free=c(T), values=c(.5), labels=c("rcw1"), lbound=-.999, ubound=.999, name="Rcw1") 
pathEw1R	<-mxMatrix( type="Stand", nrow=2, ncol=2, free=c(T), values=c(.5), labels=c("rew1"), lbound=-.999, ubound=.999, name="Rew1") 
pathAw2R	<-mxMatrix( type="Stand", nrow=2, ncol=2, free=c(T), values=c(.4), labels=c("raw2"), lbound=-.999, ubound=.999, name="Rar")
pathCw2R	<-mxMatrix( type="Stand", nrow=2, ncol=2, free=c(T), values=c(.5), labels=c("rcw2"), lbound=-.999, ubound=.999, name="Rcr") 
pathEw2R	<-mxMatrix( type="Stand", nrow=2, ncol=2, free=c(T), values=c(.5), labels=c("rew2"), lbound=-.999, ubound=.999, name="Rer") 

pathME	<-mxMatrix( type="Diag", nrow=nv, ncol=nv, free=c(F,F,F,F), values=c(0,0,0,0), labels=c(NA,NA,"me1","me2"), name="Me") # Measurement Error for Variable 3 and 4 only -original FR

pathA12R	<-mxMatrix( type="Full", nrow=2, ncol=2, free=T, values=c(.2), labels=c("ra11r","ra12r","ra21r","ra22r"), lbound=-.999, ubound=.999, name="Raxt")
pathC12R	<-mxMatrix( type="Full", nrow=2, ncol=2, free=T, values=c(.2), labels=c("rc11r","rc12r","rc21r","rc22r"), lbound=-.999, ubound=.999, name="Rcxt")

matZ22	<-mxMatrix( type="Zero", nrow=2, ncol=2, name="Z22" )
matZ42	<-mxMatrix(type="Zero", nrow=nv, ncol=2, name="Z42")
matZ44	<-mxMatrix(type="Zero", nrow=nv, ncol=nv, name="Z44")

# Matrices generated to hold A, C and E computed Variance Components within a person not including causal paths as well as ME variance
CA	<-mxAlgebra( expression=rbind( cbind(a1%&%Raw1,a1%*%Raxt%*%t(ar)),  cbind(a1%*%Raxt%*%t(ar), ar%&%Rar) ) , name="A" )
CC	<-mxAlgebra( expression=rbind( cbind(c1%&%Rcw1,c1%*%Rcxt%*%t(cr)),  cbind(c1%*%Rcxt%*%t(cr), cr%&%Rcr) ) , name="C" )
CE	<-mxAlgebra( expression=rbind( cbind(e1%&%Rew1,Z22),  cbind(Z22, er%&%Rer) ) , name="E" )
CME	<-mxAlgebra(expression= Me%*%Me, name="ME") 

# Single-headed Causal paths between the 4 variables 
causP	<-mxMatrix(type="Full", nrow=2, ncol=2, free=c(T,T,T,T), labels=c("stabCD","CDonTR","TRonCD","stabTR"), values=c(.6,.1,.1,.6), name="P") 
Caus	<-mxAlgebra(expression= cbind( rbind(Z22, P), Z42), name="SingArInd") 

# To blow it up for a pair
CausS		<-mxAlgebra(expression= rbind( cbind(SingArInd, Z44), cbind(Z44, SingArInd) ), name="SingArPair") 

matI88	<-mxMatrix(type="Iden", nrow=ntv, ncol=ntv, name="I88")
matI22	<-mxMatrix(type="Iden", nrow=2, ncol=2, name="I22")

# These are the MZ and DZ variance-covariance matrices that we need
covMZ	<-mxAlgebra(expression= (solve(I88-SingArPair)) %&% (rbind( cbind(A+C+E+ME,A+C),cbind(A+C,A+C+E+ME)) ) , 	name="ExpCovMZ") 
covDZ	<-mxAlgebra(expression= (solve(I88-SingArPair)) %&% (rbind( cbind(A+C+E+ME,.5%x%A+C),cbind(.5%x%A+C,A+C+E+ME)) ) , 	name="ExpCovDZ") 

# ----------------------------------------------------------------------------------------------------------------------------------------------
# Calculator in the script to generate standardized estimates (no implications on the model fit)

# Total var/cov of all 4 variables
matI44	<-mxMatrix(type="Iden", nrow=4, ncol=4, name="I44")
covT		<-mxAlgebra(expression= ((solve(I44-SingArInd)) %&% (A+C+E+ME)) , name="V") 		# Total Var/Cov of an indiv (including ME)										
corT		<-mxAlgebra( expression= solve(sqrt(I44*V)) %&% V  , 	name="ReT")				# Implied Cor between the 4 variables
							
# Stand A, C and E (Total)
standA	<-mxAlgebra( expression= ((solve(I44-SingArInd)) %&% A) /V, name="stA" )												
standC	<-mxAlgebra( expression= ((solve(I44-SingArInd)) %&% C) /V, name="stC" )												
standE1	<-mxAlgebra( expression= ((solve(I44-SingArInd)) %&% (E)) /V, name="stEvar12" )
standE2	<-mxAlgebra( expression= ((solve(I44-SingArInd)) %&% (E+ME)) /V, name="stEvar34" )

# Stand A, C and E (Residual, time-2 specific variance in elements 3,3 and 4,4)
standAr	<-mxAlgebra( expression= A/V, name="stAr" )												
standCr	<-mxAlgebra( expression= C/V, name="stCr" )												
standEr	<-mxAlgebra( expression= E/V, name="stEr" )													
standMe	<-mxAlgebra( expression= ME/V, name="stMe" )
standTEr	<-mxAlgebra( expression= E+ME/V, name="stTEr" )

# Stand causal STABILITY paths, square to get variance
StStab1	<-mxAlgebra( expression= ( P[1,1] * sqrt(V[1,1])) /(sqrt(V[3,3]))  , name="stStabCDPath" )
StStab2	<-mxAlgebra( expression= ( P[2,2] * sqrt(V[2,2])) /(sqrt(V[4,4]))  , name="stStabTRPath" )

# Stand causal CROSS-LAGG paths, square to get variance
StCaus1	<-mxAlgebra( expression= ( P[2,1] * sqrt(V[1,1])) /(sqrt(V[4,4])) , name="stCDonTR" )
StCaus2	<-mxAlgebra( expression= ( P[1,2] * sqrt(V[2,2])) /(sqrt(V[3,3])) , name="stTRonCD" )

# Proportion of cross-time correlations due to A, C, E and causal paths
# CD3->TR4
# Through cross-time cross-trait rA and rC
AC3_TR4_xrA		<-mxAlgebra( expression= ( sqrt(stA[1,1])* Raxt[2,1]*sqrt(stAr[4,4]))  , name="CD3_TR4_rxA" )
CC3_TR4_xrC		<-mxAlgebra( expression= ( sqrt(stC[1,1])* Rcxt[2,1]*sqrt(stCr[4,4]))  , name="CD3_TR4_rxC" )
# Through within-time rA, rC and rE (Time 1) and stability paths (TR3->TR4)
RA3_Stab_TR4	<-mxAlgebra( expression= ( sqrt(stA[1,1])* Raw1[2,1]*sqrt(stA[2,2])*stStabTRPath)  , name="CD3_rAw_StabTR4" )
RC3_Stab_TR4	<-mxAlgebra( expression= ( sqrt(stC[1,1])* Rcw1[2,1]*sqrt(stC[2,2])*stStabTRPath)  , name="CD3_rCw_StabTR4" )
RE3_Stab_TR4	<-mxAlgebra( expression= ( sqrt(stEvar12[1,1])* Rew1[2,1]*sqrt(stEvar12[2,2])*stStabTRPath)  , name="CD3_rEw_StabTR4" )

# Proportions of each effect
AC3_TR4_Prop	<-mxAlgebra( expression= ( (CD3_TR4_rxA)/(CD3_TR4_rxA+CD3_TR4_rxC+stCDonTR[1,1]+CD3_rAw_StabTR4+CD3_rCw_StabTR4+CD3_rEw_StabTR4))  , name="CD3_TR4_Aprop" )
CC3_TR4_Prop	<-mxAlgebra( expression= ( (CD3_TR4_rxC)/(CD3_TR4_rxA+CD3_TR4_rxC+stCDonTR[1,1]+CD3_rAw_StabTR4+CD3_rCw_StabTR4+CD3_rEw_StabTR4))  , name="CD3_TR4_Cprop" )
RA3_Stab_TR4_Prop	<-mxAlgebra( expression= ( (CD3_rAw_StabTR4)/(CD3_TR4_rxA+CD3_TR4_rxC+stCDonTR[1,1]+CD3_rAw_StabTR4+CD3_rCw_StabTR4+CD3_rEw_StabTR4))  , name="CD3_Stab_TR4_Aprop" )
RC3_Stab_TR4_Prop	<-mxAlgebra( expression= ( (CD3_rCw_StabTR4)/(CD3_TR4_rxA+CD3_TR4_rxC+stCDonTR[1,1]+CD3_rAw_StabTR4+CD3_rCw_StabTR4+CD3_rEw_StabTR4))  , name="CD3_Stab_TR4_Cprop" )
RE3_Stab_TR4_Prop	<-mxAlgebra( expression= ( (CD3_rEw_StabTR4)/(CD3_TR4_rxA+CD3_TR4_rxC+stCDonTR[1,1]+CD3_rAw_StabTR4+CD3_rCw_StabTR4+CD3_rEw_StabTR4))  , name="CD3_Stab_TR4_Eprop" )
CausC3_TR4_Prop	<-mxAlgebra( expression= ( (stCDonTR[1,1])/(CD3_TR4_rxA+CD3_TR4_rxC+stCDonTR[1,1]+CD3_rAw_StabTR4+CD3_rCw_StabTR4+CD3_rEw_StabTR4))  , name="CD3_TR4_Causprop" )
R_Sum_CDTR		<-mxAlgebra( expression= ( (CD3_TR4_rxA+CD3_TR4_rxC+stCDonTR[1,1]+CD3_rAw_StabTR4+CD3_rCw_StabTR4+CD3_rEw_StabTR4))  , name="r_Sum_CDTR" )

# TR->CD
# Through cross-time cross-trait rA and rC
ATR3_C4_xrA		<-mxAlgebra( expression= ( sqrt(stA[2,2])* Raxt[1,2]*sqrt(stAr[3,3]))  , name="TR3_CD4_rxA" )
CTR3_C4_xrC		<-mxAlgebra( expression= ( sqrt(stC[2,2])* Rcxt[1,2]*sqrt(stCr[3,3]))  , name="TR3_CD4_rxC" )
# Through within-time rA, rC and rE and stability paths (CD3->CD4)
RA3_Stab_CD4	<-mxAlgebra( expression= ( sqrt(stA[2,2])* Raw1[2,1]*sqrt(stA[1,1])*stStabCDPath)  , name="TR3_rAw_StabCD4" )
RC3_Stab_CD4	<-mxAlgebra( expression= ( sqrt(stC[2,2])* Rcw1[2,1]*sqrt(stC[1,1])*stStabCDPath)  , name="TR3_rCw_StabCD4" )
RE3_Stab_CD4	<-mxAlgebra( expression= ( sqrt(stEvar12[2,2])* Rew1[2,1]*sqrt(stEvar12[1,1])*stStabCDPath)  , name="TR3_rEw_StabCD4" )
# Through stability paths (CD3->CD4) and within-time rA, rC and rE (time 2)

# Proportions of each effect
ATR3_C4_Prop	<-mxAlgebra( expression= ( (TR3_CD4_rxA)/(TR3_CD4_rxA+TR3_CD4_rxC+stTRonCD[1,1]+TR3_rAw_StabCD4+TR3_rCw_StabCD4+TR3_rEw_StabCD4))  , name="TR3_CD4_Aprop" )
CTR3_C4_Prop	<-mxAlgebra( expression= ( (TR3_CD4_rxC)/(TR3_CD4_rxA+TR3_CD4_rxC+stTRonCD[1,1]+TR3_rAw_StabCD4+TR3_rCw_StabCD4+TR3_rEw_StabCD4))  , name="TR3_CD4_Cprop" )
RA3_Stab_CD4_Prop	<-mxAlgebra( expression= ( (TR3_rAw_StabCD4)/(TR3_CD4_rxA+TR3_CD4_rxC+stTRonCD[1,1]+TR3_rAw_StabCD4+TR3_rCw_StabCD4+TR3_rEw_StabCD4))  , name="TR3_Stab_CD4_Aprop" )
RC3_Stab_CD4_Prop	<-mxAlgebra( expression= ( (TR3_rCw_StabCD4)/(TR3_CD4_rxA+TR3_CD4_rxC+stTRonCD[1,1]+TR3_rAw_StabCD4+TR3_rCw_StabCD4+TR3_rEw_StabCD4))  , name="TR3_Stab_CD4_Cprop" )
RE3_Stab_CD4_Prop	<-mxAlgebra( expression= ( (TR3_rEw_StabCD4)/(TR3_CD4_rxA+TR3_CD4_rxC+stTRonCD[1,1]+TR3_rAw_StabCD4+TR3_rCw_StabCD4+TR3_rEw_StabCD4))  , name="TR3_Stab_CD4_Eprop" )
CausTR3_C4_Prop	<-mxAlgebra( expression= ( (stTRonCD[1,1])/(TR3_CD4_rxA+TR3_CD4_rxC+stTRonCD[1,1]+TR3_rAw_StabCD4+TR3_rCw_StabCD4+TR3_rEw_StabCD4))  , name="TR3_CD4_Causprop" )
R_Sum_TRCD		<-mxAlgebra( expression= ( (TR3_CD4_rxA+TR3_CD4_rxC+stTRonCD[1,1]+TR3_rAw_StabCD4+TR3_rCw_StabCD4+TR3_rEw_StabCD4))  , name="r_Sum_TRCD" )

# ----------------------------------------------------------------------------------------------------------------------------------------------

# Data objects 
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="ExpCovMZ", means="ExpMean", dimnames=selVars)
objDZ		<- mxExpectationNormal( covariance="ExpCovDZ", means="ExpMean", dimnames=selVars)

fitFunction 	<- mxFitFunctionML()
 
# Combine Groups
pars1	   	<- list(pathAw1, pathAw2, pathCw1, pathCw2, pathEw1, pathEw2, pathAw1R, pathAw2R, pathCw1R, pathCw2R, pathEw1R, pathEw2R,
				pathME, CME, CA, CC, CE, causP, Caus, CausS, pathA12R, pathC12R, matZ22, matZ42, matZ44, matI88, matI44, matI22)
pars2	   	<- list(corT, covT, standA, standC,standE1, standE2, StStab1, StStab2, StCaus1, StCaus2, standAr, standCr, standEr, standMe, standTEr,
				AC3_TR4_xrA,CC3_TR4_xrC,AC3_TR4_Prop,CC3_TR4_Prop,CausC3_TR4_Prop,ATR3_C4_xrA,CTR3_C4_xrC,ATR3_C4_Prop,CTR3_C4_Prop,CausTR3_C4_Prop,
				RA3_Stab_TR4,RC3_Stab_TR4,RE3_Stab_TR4,RA3_Stab_TR4_Prop,RC3_Stab_TR4_Prop,RE3_Stab_TR4_Prop,
				RA3_Stab_CD4,RC3_Stab_CD4,RE3_Stab_CD4,RA3_Stab_CD4_Prop,RC3_Stab_CD4_Prop,RE3_Stab_CD4_Prop, R_Sum_TRCD, R_Sum_CDTR)
modelMZ	<- mxModel(pars1, pars2, covMZ, Mean, dataMZ, objMZ, fitFunction, name="MZ" )
modelDZ	<- mxModel(pars1, covDZ, Mean, dataDZ, objDZ, fitFunction, name="DZ" )
minus2ll	<- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<- mxFitFunctionAlgebra("m2LL" )
ciCaus1	<- mxCI (c  ('MZ.stCDonTR[1,1]'))
ciCaus2	<- mxCI (c  ('MZ.stTRonCD[1,1]'))
ciStabRTVb	<- mxCI (c  ('MZ.stStabCDPath[1,1]'))
ciStabAsb	<- mxCI (c  ('MZ.stStabTRPath[1,1]'))
cistA		<- mxCI (c  ('MZ.stA[1,1]', 'MZ.stA[2,2]', 'MZ.stA[3,1]', 'MZ.stA[4,1]', 'MZ.stA[3,2]', 'MZ.stA[4,2]', 'MZ.stAr[3,3]', 'MZ.stAr[4,4]' ))
cistC		<- mxCI (c  ('MZ.stC[1,1]', 'MZ.stC[2,2]', 'MZ.stC[3,1]', 'MZ.stC[4,1]', 'MZ.stC[3,2]', 'MZ.stC[4,2]', 'MZ.stCr[3,3]', 'MZ.stCr[4,4]' ))
cistE		<- mxCI (c  ('MZ.stEvar34[1,1]',  'MZ.stEvar34[2,2]', 'MZ.stEvar12[1,1]', 'MZ.stEvar12[2,2]')) # E total IQ1 and IQ2 and E on IQ1 without ME
cistE2	<- mxCI (c  ('MZ.stEr[3,3]', 'MZ.stEr[4,4]')) # Eres on IQ2 and AS2
ciRaw1	<- mxCI ('MZ.Raw1[1,2]')
ciRcw1	<- mxCI ('MZ.Rcw1[1,2]')
ciRew1 	<- mxCI ('MZ.Rew1[1,2]')
ciRar 	<- mxCI ('MZ.Rar[2,1]')
ciRcr 	<- mxCI ('MZ.Rcr[2,1]')
ciRer 	<- mxCI ('MZ.Rer[2,1]')
ciRa12 	<- mxCI (c('MZ.Raxt[1,1]','MZ.Raxt[1,2]','MZ.Raxt[2,1]','MZ.Raxt[2,2]'))
ciRc12 	<- mxCI (c('MZ.Rcxt[1,1]','MZ.Rcxt[1,2]','MZ.Rcxt[2,1]','MZ.Rcxt[2,2]'))

ACEcliiModel	<- mxModel( "ACEclii", modelMZ, modelDZ, minus2ll, obj, ciCaus1, ciCaus2, ciStabRTVb, ciStabAsb, cistA, cistC, cistE, cistE2, ciRaw1, ciRcw1, ciRew1, ciRar, ciRcr, ciRer, ciRa12, ciRc12   ) 

#----------------------------------------------------------------------------------------------------------------------------------------------------
# RUN ACEclii model 
ACEcliiFit		<- mxRun(ACEcliiModel, intervals=F)
(ACEcliiSumm	<- summary(ACEcliiFit,verbose=F))

mxCompare(ACEcliiFit, ACEclFit)


# OUTPUT specific parameters of interest of sub model

# checks
mxEval(MZ.A,ACEcliiFit)
mxEval(MZ.C,ACEcliiFit)
mxEval(MZ.E,ACEcliiFit)
mxEval(MZ.V,ACEcliiFit)
mxEval(MZ.ReT,ACEcliiFit)
mxEval(MZ.r_Sum_CDTR,ACEcliiFit)
mxEval(MZ.r_Sum_TRCD,ACEcliiFit)

mxEval(MZ.Raw1,ACEcliiFit)
mxEval(MZ.Rcw1,ACEcliiFit)
mxEval(MZ.Rew1,ACEcliiFit)
mxEval(MZ.Rar,ACEcliiFit)
mxEval(MZ.Rcr,ACEcliiFit)
mxEval(MZ.Rer,ACEcliiFit)

#,,,,,,,,,,

mxEval(MZ.CD3_TR4_rxA,ACEcliiFit)
mxEval(MZ.CD3_TR4_rxC,ACEcliiFit)
mxEval(MZ.CD3_rAw_StabTR4,ACEcliiFit)
mxEval(MZ.CD3_rCw_StabTR4,ACEcliiFit)
mxEval(MZ.CD3_rEw_StabTR4,ACEcliiFit)
mxEval(MZ.stCDonTR,ACEcliiFit)
mxEval(MZ.stStabCDPath,ACEcliiFit)
mxEval(MZ.stStabTRPath,ACEcliiFit)

mxEval(MZ.CD3_TR4_Aprop,ACEcliiFit)
mxEval(MZ.CD3_TR4_Cprop,ACEcliiFit)
mxEval(MZ.CD3_Stab_TR4_Aprop,ACEcliiFit)
mxEval(MZ.CD3_Stab_TR4_Cprop,ACEcliiFit)
mxEval(MZ.CD3_Stab_TR4_Eprop,ACEcliiFit)
mxEval(MZ.CD3_TR4_Causprop,ACEcliiFit)

mxEval(MZ.TR3_CD4_rxA,ACEcliiFit)
mxEval(MZ.TR3_CD4_rxC,ACEcliiFit)
mxEval(MZ.TR3_rAw_StabCD4,ACEcliiFit)
mxEval(MZ.TR3_rCw_StabCD4,ACEcliiFit)
mxEval(MZ.TR3_rEw_StabCD4,ACEcliiFit)
mxEval(MZ.StabTR4_rAr_CD4,ACEcliiFit)
mxEval(MZ.StabTR4_rCr_CD4,ACEcliiFit)
mxEval(MZ.StabTR4_rEr_CD4,ACEcliiFit)
mxEval(MZ.stTRonCD,ACEcliiFit)

mxEval(MZ.TR3_CD4_Aprop,ACEcliiFit)
mxEval(MZ.TR3_CD4_Cprop,ACEcliiFit)
mxEval(MZ.TR3_Stab_CD4_Aprop,ACEcliiFit)
mxEval(MZ.TR3_Stab_CD4_Cprop,ACEcliiFit)
mxEval(MZ.TR3_Stab_CD4_Eprop,ACEcliiFit)
mxEval(MZ.Stab_TR4_CD4_Aprop,ACEcliiFit)
mxEval(MZ.Stab_TR4_CD4_Aprop,ACEcliiFit)
mxEval(MZ.Stab_TR4_CD4_Aprop,ACEcliiFit)
mxEval(MZ.TR3_CD4_Causprop,ACEcliiFit)

# Unstandardized x-lag paths
mxEval(MZ.SingArPair,ACEcliiFit)

# Standardized variance components
# Note that elements 3,3 and 4,4, give the Total F and E variance for the time-two variables (i.e. the transmitted + Residual effects)
mxEval(MZ.stA,ACEcliiFit)
mxEval(MZ.stC,ACEcliiFit)
# note that E for var 1 does not include ME
mxEval(MZ.stEvar12,ACEcliiFit)
mxEval(MZ.stEvar34,ACEcliiFit)
mxEval(MZ.stMe,ACEcliiFit)

# Note that elements 3,3 and 4,4, give the Residual A, C and E variance for the time-two variables 
mxEval(MZ.stAr,ACEcliiFit)
mxEval(MZ.stCr,ACEcliiFit)
mxEval(MZ.stEr,ACEcliiFit)

# Standardized M error for 1,1, and 2,2
mxEval(MZ.stMe,ACEcliiFit)

#------------------------------------------------------------------------
# Submodel 1a: Drop non-significant paths from previous model 
#------------------------------------------------------------------------

ACEcliiaModel	<- mxModel(ACEcliiFit, name="ACEcliia")
ACEcliiaModel	<- omxSetParameters(ACEcliiaModel, labels=c('raw1','rew2','ra21r','TRonCD'), free=FALSE, values=0)
ACEcliiaFit		<- mxRun(ACEcliiaModel, intervals=F)
(ACEcliiaSumm	<- summary(ACEcliiaFit))

mxCompare(ACEcliiFit, ACEcliiaFit)

mxEval(MZ.ReT,ACEcliiaFit)
mxEval(MZ.stStabCDPath,ACEcliiaFit)
mxEval(MZ.stStabTRPath,ACEcliiaFit)

mxEval(MZ.stCDonTR,ACEcliiaFit)
mxEval(MZ.CD3_TR4_rxA,ACEcliiaFit)
mxEval(MZ.CD3_TR4_rxC,ACEcliiaFit)
mxEval(MZ.CD3_rAw_StabTR4,ACEcliiaFit)
mxEval(MZ.CD3_rCw_StabTR4,ACEcliiaFit)
mxEval(MZ.CD3_rEw_StabTR4,ACEcliiaFit)
mxEval(MZ.r_Sum_CDTR,ACEcliiaFit)

mxEval(MZ.CD3_TR4_Causprop,ACEcliiaFit)
mxEval(MZ.CD3_TR4_Aprop,ACEcliiaFit)
mxEval(MZ.CD3_TR4_Cprop,ACEcliiaFit)
mxEval(MZ.CD3_Stab_TR4_Aprop,ACEcliiaFit)
mxEval(MZ.CD3_Stab_TR4_Cprop,ACEcliiaFit)
mxEval(MZ.CD3_Stab_TR4_Eprop,ACEcliiaFit)

mxEval(MZ.stTRonCD,ACEcliiaFit)
mxEval(MZ.TR3_CD4_rxA,ACEcliiaFit)
mxEval(MZ.TR3_CD4_rxC,ACEcliiaFit)
mxEval(MZ.TR3_rAw_StabCD4,ACEcliiaFit)
mxEval(MZ.TR3_rCw_StabCD4,ACEcliiaFit)
mxEval(MZ.TR3_rEw_StabCD4,ACEcliiaFit)
mxEval(MZ.r_Sum_TRCD,ACEcliiaFit)

mxEval(MZ.TR3_CD4_Causprop,ACEcliiaFit)
mxEval(MZ.TR3_CD4_Aprop,ACEcliiaFit)
mxEval(MZ.TR3_CD4_Cprop,ACEcliiaFit)
mxEval(MZ.TR3_Stab_CD4_Aprop,ACEcliiaFit)
mxEval(MZ.TR3_Stab_CD4_Cprop,ACEcliiaFit)
mxEval(MZ.TR3_Stab_CD4_Eprop,ACEcliiaFit)



#****************************************************************************************************************************
# __(6c)_____________________________________________________________________________________________________________________
# CROSS-LAG MODEL with var/cov specified with Additive genetic and shared and ind-specific Environ effects
# And causal paths between observed variables
# For CD and TR
# BURT model
# Including cross-time rA and rC correlations for both the stability and cross-lag paths but using the Cholesky decomposition
# ___________________________________________________________________________________________________________________________


nv		<- 4				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 2				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nfcor		<- (nfact*(nfact-1)/2)	# number of free elements in a correlation matrix nfact*nfcat
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv-1)/2)		# number of free elements in a correlation matrix nv*nv
Groups	<-c("mz", "dz")

Vars		<- c("Parca3","TalkRhym3","Parca4","TalkRhym4")
selVars	<- c("Parca31","TalkRhym31","Parca41","TalkRhym41",
			"Parca32","TalkRhym32","Parca42","TalkRhym42")

mzData	<- subset(TWINdata2, zyg1%in%c(1)|zyg2%in%c(1), selVars)
dzData	<- subset(TWINdata2, zyg1%in%c(2)|zyg2%in%c(2), selVars)

psych::describe(mzData)
psych::describe(dzData)

# Create labels and starting values

(Stmean	<-colMeans(mzData[,1:nv],na.rm=TRUE))
StM 		<-c(Stmean, Stmean)

(LabM		<- paste("m",1:nv,sep=""))
MLabs		<-c(LabM,LabM) 

# Mean
Mean		<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, labels=MLabs, values=StM, name="ExpMean")

# Matrices to store a, c and e Path Coefficients
pathAw1	<-mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=.3, labels=c("a11","a21","a31","a41","a22","a32","a42","a33","a43","a44"), name="a1") 
pathCw1	<-mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=.7, labels=c("c11","c21","c31","c41","c22","c32","c42","c33","c43","c44"), name="c1")
pathEw1	<-mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=.3, labels=c("e11","e21","e31","e41","e22","e32","e42","e33","e43","e44"), name="e1") 
pathME	<-mxMatrix( type="Diag", nrow=nv, ncol=nv, free=c(F,F,F,F), values=c(0,0,0,0), labels=c(NA,NA,"me1","me2"), name="Me") # Measurement Error for Variable 3 and 4 only -original FR

matZ22	<-mxMatrix( type="Zero", nrow=2, ncol=2, name="Z22" )
matZ42	<-mxMatrix(type="Zero", nrow=nv, ncol=2, name="Z42")
matZ44	<-mxMatrix(type="Zero", nrow=nv, ncol=nv, name="Z44")

# Matrices generated to hold A, C and E computed Variance Components within a person not including causal paths as well as ME variance
CA	<-mxAlgebra( expression= a1%*%t(a1) , name="A" )
CC	<-mxAlgebra( expression= c1%*%t(c1) , name="C" )
CE	<-mxAlgebra( expression= e1%*%t(e1) , name="E" )
CME	<-mxAlgebra(expression= Me%*%Me, name="ME") 

# Single-headed Causal paths between the 4 variables 
causP	<-mxMatrix(type="Full", nrow=2, ncol=2, free=F, labels=c("stabCD","CDonTR","TRonCD","stabTR"), values=0, name="P") 
Caus	<-mxAlgebra(expression= cbind( rbind(Z22, P), Z42), name="SingArInd") 

# To blow it up for a pair
CausS		<-mxAlgebra(expression= rbind( cbind(SingArInd, Z44), cbind(Z44, SingArInd) ), name="SingArPair") 

matI88	<-mxMatrix(type="Iden", nrow=ntv, ncol=ntv, name="I88")
matI22	<-mxMatrix(type="Iden", nrow=2, ncol=2, name="I22")

# These are the MZ and DZ variance-covariance matrices that we need
covMZ	<-mxAlgebra(expression= (solve(I88-SingArPair)) %&% (rbind( cbind(A+C+E+ME,A+C),cbind(A+C,A+C+E+ME)) ) , 	name="ExpCovMZ") 
covDZ	<-mxAlgebra(expression= (solve(I88-SingArPair)) %&% (rbind( cbind(A+C+E+ME,.5%x%A+C),cbind(.5%x%A+C,A+C+E+ME)) ) , 	name="ExpCovDZ") 

# ----------------------------------------------------------------------------------------------------------------------------------------------
# Calculator in the script to generate standardized estimates (no implications on the model fit)

# Total var/cov of all 4 variables
matI44	<-mxMatrix(type="Iden", nrow=4, ncol=4, name="I44")
covT		<-mxAlgebra(expression= ((solve(I44-SingArInd)) %&% (A+C+E+ME)) , name="V") 		# Total Var/Cov of an indiv (including ME)										
corT		<-mxAlgebra( expression= solve(sqrt(I44*V)) %&% V  , 	name="ReT")				# Implied Cor between the 4 variables
							
# Stand A, C and E (Total)
standA	<-mxAlgebra( expression= ((solve(I44-SingArInd)) %&% A) /V, name="stA" )												
standC	<-mxAlgebra( expression= ((solve(I44-SingArInd)) %&% C) /V, name="stC" )												
standE1	<-mxAlgebra( expression= ((solve(I44-SingArInd)) %&% (E)) /V, name="stEvar12" )
standE2	<-mxAlgebra( expression= ((solve(I44-SingArInd)) %&% (E+ME)) /V, name="stEvar34" )

# Stand A, C and E (Residual, time-2 specific variance in elements 3,3 and 4,4)
standAr	<-mxAlgebra( expression= A/V, name="stAr" )												
standCr	<-mxAlgebra( expression= C/V, name="stCr" )												
standEr	<-mxAlgebra( expression= E/V, name="stEr" )													
standMe	<-mxAlgebra( expression= ME/V, name="stMe" )
standTEr	<-mxAlgebra( expression= E+ME/V, name="stTEr" )

# Stand causal STABILITY paths, square to get variance
StStab1	<-mxAlgebra( expression= ( P[1,1] * sqrt(V[1,1])) /(sqrt(V[3,3]))  , name="stStabCDPath" )
StStab2	<-mxAlgebra( expression= ( P[2,2] * sqrt(V[2,2])) /(sqrt(V[4,4]))  , name="stStabTRPath" )

# Stand causal CROSS-LAGG paths, square to get variance
StCaus1	<-mxAlgebra( expression= ( P[2,1] * sqrt(V[1,1])) /(sqrt(V[4,4])) , name="stCDonTR" )
StCaus2	<-mxAlgebra( expression= ( P[1,2] * sqrt(V[2,2])) /(sqrt(V[3,3])) , name="stTRonCD" )

# Derive standardised path estimates
SD		<-mxAlgebra( expression= solve(sqrt(I44*V)), name="iSD")
stpthA	<-mxAlgebra( expression= a1%*%iSD, name="sta" )												
stpthC	<-mxAlgebra( expression= c1%*%iSD, name="stc" )												
stpthE	<-mxAlgebra( expression= e1%*%iSD, name="ste" )												

# ----------------------------------------------------------------------------------------------------------------------------------------------

# Data objects 
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="ExpCovMZ", means="ExpMean", dimnames=selVars)
objDZ		<- mxExpectationNormal( covariance="ExpCovDZ", means="ExpMean", dimnames=selVars)

fitFunction 	<- mxFitFunctionML()
 
# Combine Groups
pars1	   	<- list(pathAw1, pathCw1, pathEw1,pathME, CME, CA, CC, CE, causP, Caus, CausS, matZ22, matZ42, matZ44, matI88, matI44, matI22, SD, stpthA, stpthC, stpthE)
pars2	   	<- list(corT, covT, standA, standC,standE1, standE2, StStab1, StStab2, StCaus1, StCaus2, standAr, standCr, standEr, standMe, standTEr)
modelMZ	<- mxModel(pars1, pars2, covMZ, Mean, dataMZ, objMZ, fitFunction, name="MZ" )
modelDZ	<- mxModel(pars1, pars2, covDZ, Mean, dataDZ, objDZ, fitFunction, name="DZ" )
minus2ll	<- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<- mxFitFunctionAlgebra("m2LL" )
ciCaus1	<- mxCI (c  ('MZ.stCDonTR[1,1]'))
ciCaus2	<- mxCI (c  ('MZ.stTRonCD[1,1]'))
ciStabRTVb	<- mxCI (c  ('MZ.stStabCDPath[1,1]'))
ciStabAsb	<- mxCI (c  ('MZ.stStabTRPath[1,1]'))
cistA		<- mxCI (c  ('MZ.stA[1,1]', 'MZ.stA[2,2]', 'MZ.stA[3,1]', 'MZ.stA[4,1]', 'MZ.stA[3,2]', 'MZ.stA[4,2]', 'MZ.stAr[3,3]', 'MZ.stAr[4,4]' ))
cistC		<- mxCI (c  ('MZ.stC[1,1]', 'MZ.stC[2,2]', 'MZ.stC[3,1]', 'MZ.stC[4,1]', 'MZ.stC[3,2]', 'MZ.stC[4,2]', 'MZ.stCr[3,3]', 'MZ.stCr[4,4]' ))
cistE		<- mxCI (c  ('MZ.stEvar34[1,1]',  'MZ.stEvar34[2,2]', 'MZ.stEvar12[1,1]', 'MZ.stEvar12[2,2]')) # E total IQ1 and IQ2 and E on IQ1 without ME
cistE2	<- mxCI (c  ('MZ.stEr[3,3]', 'MZ.stEr[4,4]')) # Eres on IQ2 and AS2

ACEcholliiModel	<- mxModel( "ACEchollii", modelMZ, modelDZ, minus2ll, obj, ciCaus1, ciCaus2, ciStabRTVb, ciStabAsb, cistA, cistC, cistE, cistE2) 

#----------------------------------------------------------------------------------------------------------------------------------------------------
# RUN ACEclii model 
ACEcholliiFit		<- mxRun(ACEcholliiModel, intervals=F)
(ACEcholliiSumm	<- summary(ACEcholliiFit,verbose=F))

mxCompare(ACEcholliiFit, ACEclFit)
mxCompare(ACEcholliiFit, ACEcliiFit)


# OUTPUT specific parameters of interest of sub model

# checks
mxEval(MZ.A,ACEcholliiFit)
mxEval(MZ.C,ACEcholliiFit)
mxEval(MZ.E,ACEcholliiFit)
mxEval(MZ.V,ACEcholliiFit)
mxEval(MZ.ReT,ACEcholliiFit)

mxEval(MZ.stA,ACEcholliiFit)
mxEval(MZ.stC,ACEcholliiFit)
mxEval(MZ.stEvar12,ACEcholliiFit)
mxEval(MZ.stEvar34,ACEcholliiFit)
mxEval(MZ.stAr,ACEcholliiFit)
mxEval(MZ.stCr,ACEcholliiFit)
mxEval(MZ.stEr,ACEcholliiFit)


#------------------------------------------------------------------------
# Submodel 1a: Cross-time E covariances from previous model 
# And estimate causal paths
#------------------------------------------------------------------------

SubACEcholliiModel	<- mxModel(ACEcholliiFit, name="SubACEchollii")
SubACEcholliiModel	<- omxSetParameters(SubACEcholliiModel, labels=c('e31','e32','e41','e42'), free=FALSE, values=0)
SubACEcholliiModel	<- omxSetParameters(SubACEcholliiModel, labels=c('CDonTR','TRonCD','stabCD','stabTR'), free=TRUE, values=.2)
SubACEcholliiModel	<- omxAssignFirstParameters(SubACEcholliiModel)
SubACEcholliiFit		<- mxRun(SubACEcholliiModel, intervals=T)
(SubACEcholliiSumm	<- summary(SubACEcholliiFit))

mxCompare(ACEcholliiFit, SubACEcholliiFit)
mxCompare(SubACEcholliiFit, ACEcliiFit)
mxCompare(SubACEcholliiFit, ACEcliicFit)

mxEval(MZ.ReT,SubACEcholliiFit)
mxEval(MZ.ExpCovMZ,SubACEcholliiFit)
mxEval(DZ.ExpCovDZ,SubACEcholliiFit)
mxEval(MZ.stStabCDPath,SubACEcholliiFit)
mxEval(MZ.stStabTRPath,SubACEcholliiFit)
mxEval(MZ.stCDonTR,SubACEcholliiFit)
mxEval(MZ.stTRonCD,SubACEcholliiFit)

mxEval(MZ.A,SubACEcholliiFit)
mxEval(MZ.C,SubACEcholliiFit)
mxEval(MZ.E,SubACEcholliiFit)
mxEval(MZ.V,SubACEcholliiFit)

mxEval(MZ.stA,SubACEcholliiFit)
mxEval(MZ.stC,SubACEcholliiFit)
mxEval(MZ.stEvar12,SubACEcholliiFit)
mxEval(MZ.stEvar34,SubACEcholliiFit)
mxEval(MZ.stAr,SubACEcholliiFit)
mxEval(MZ.stCr,SubACEcholliiFit)
mxEval(MZ.stEr,SubACEcholliiFit)
mxEval(MZ.sta,SubACEcholliiFit)
mxEval(MZ.stc,SubACEcholliiFit)
mxEval(MZ.ste,SubACEcholliiFit)


