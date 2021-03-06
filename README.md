# UncertaintiesReturnValues

Estimation of return values in the presence of uncertain extreme value model parameters.

Content supports the publication "Uncertainties in return values from extreme value analysis using the generalised Pareto (GP) distribution" by Jonathan, Randell, Wadsworth and Tawn (2020).

# Code files

Code to (a) simulate GP data for samples of different sizes, (b) estimate return values using (i) different estimators and (ii) different estimation schemes, and (c) characterise return value estimates obtained in terms of bias of (i) return value, (ii) exceedance probability and (iii) log exceedance probability.

RunSimulation = main simulation code

gpfitPWM = PWM fit

gpfitMOM = MOM fit

gpfitEB = EB fit (Zhang 2010)

AnnCdf2 = distribution of annual maximum

AnnExcPrb = annual exceedance probability

RtrVal = return value

pLgn = nice legends for plots

pGI = nice png plots

pDflBig = nice font size

pAxsLmt = nice axis limits

(Original versions of EB, MOM and PWM code provided by Ed Mackay, checked by PhJ.)

# Figure files

Figure files provide supporting evidence for the discussion in Section 7 of paper, where different estimation schemes for GP model are considered. Use RunSimulation to generate more!

AllMthPrmSct.png = GP parameter estimate scatter as a function of sample size and estimation method

FrcBiasEB.png = fractional bias in return values for empirical Bayes (EB) approach

FrcBiasMOM.png = fractional bias in return values using method of moments (MOM)

FrcBiasPWM.png = fractional bias in return values using probability weighted moments (PWM)

