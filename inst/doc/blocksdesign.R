## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(fig.width=6, fig.height=4) 
knitr::opts_chunk$set(echo = FALSE)

## ------------------------------------------------------------------------

a1=c("reps",  2,1)
a2=c("reps+sub1", 17, 0.8769)
a3=c("reps+sub1+sub2", 35, 0.7275)
a4=c("reps+sub1+sub2+sub3", 71, 0.4113)
add=data.frame(rbind(a1,a2,a3,a4))
colnames(add)=c("Additive model","DF","Efficiency")
rownames(add)=NULL
knitr::kable(add, caption = "Efficiency factors for nested blocks model")


## ---- echo=FALSE---------------------------------------------------------
Hist1=vector("list", 6)
class(Hist1)="histogram"
Hist1$breaks=c(0,  20,  40,  60,  80, 100, 120, 140, 160)
Hist1$counts=c(148, 313, 283, 168,  64,  18,   4,   2)
Hist1$density=c(0.00740, 0.01565, 0.01415, 0.00840, 0.00320, 0.00090, 0.00020, 0.00010)
Hist1$mids=c(10,  30,  50,  70,  90, 110, 130, 150)
Hist1$xname="ChiSq difference"
Hist1$equidist="TRUE"
plot(Hist1, main="Fig 1. One level of nesting versus complete blocks ")

Hist2=vector("list", 6)
class(Hist2)="histogram"
Hist2$breaks=c(0,  5, 10, 15, 20, 25, 30, 35, 40, 45)
Hist2$counts=c(508, 215, 116,  87,  36,  23,  10,   4,   1)
Hist2$density=c(0.1016, 0.0430, 0.0232, 0.0174, 0.0072, 0.0046, 0.0020, 0.0008, 0.0002)
Hist2$mids=c(2.5,  7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5)
Hist2$xname="ChiSq difference"
Hist2$equidist="TRUE"
plot(Hist2, main="Fig 2. Two levels of nesting versus one level ")

Hist3=vector("list", 6)
class(Hist3)="histogram"
Hist3$breaks=c(0,  5, 10, 15, 20, 25, 30, 35, 40, 45)
Hist3$counts=c( 751, 136,  58,  26,  13,   7,   4,   4,   1)
Hist3$density3=c(0.1502, 0.0272, 0.0116, 0.0052, 0.0026, 0.0014, 0.0008, 0.0008, 0.0002)
Hist3$mids=c(2.5,  7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5)
Hist3$xname="ChiSq difference"
Hist3$equidist="TRUE"
plot(Hist3,  main="Fig 3. Three levels of nesting versus two levels")


## ---- echo=FALSE---------------------------------------------------------

m0=c("(1 | Reps) + (1 | Rows)")
m1=c("(1 | Reps) + (1 | Rows) + (1 | Col1)")
m2=c("(1 | Reps) + (1 | Rows) + (1 | Col1) + (1 | Rows:Col1)")
m3=c("(1 | Reps) + (1 | Rows) + (1 | Col1) + (1 | Rows:Col1) + (1 | Col2)")
m4=c("(1 | Reps) + (1 | Rows) + (1 | Col1) + (1 | Rows:Col1) + (1 | Col2) + (1 | Rows:Col2)")
m5=c("(1 | Reps) + (1 | Rows) + (1 | Col1) + (1 | Rows:Col1) + (1 | Col2) + (1 | Rows:Col2) + (1 | Col3)")
mod=data.frame(rbind(m0,m1,m2,m3,m4,m5))
colnames(mod)=c("Block model terms")
rownames(mod)=NULL
knitr::kable(mod, caption = "Sequentially fitted block model terms")


## ---- echo=FALSE---------------------------------------------------------

a1=c(0, 275, 584.32, 1766.5, -17.161,   34.322, " ", " ", " ")
a2=c(1, 276, 476.77, 1663.3,  37.615,  -75.230, 109.5521,      1,  "< 2.2e-16 ***")
a3=c(2, 277, 429.72, 1620.5,  62.138, -124.276,  49.0459,      1,  "2.500e-12 ***")
a4=c(3, 278, 366.40, 1561.5,  94.801, -189.601,  65.3248,      1,  "6.352e-16 ***")
a5=c(4, 279, 363.16, 1562.6,  97.420, -194.841,   5.2397,      1,  "0.02208 *" )
a6=c(5, 280, 282.81, 1486.5, 138.593, -277.187,  82.3458,      1,  "< 2.2e-16 ***")
add=data.frame(rbind(a1,a2,a3,a4,a5,a6))
colnames(add)=c("Model", "Df",    "AIC",    "BIC",  "logLik", "deviance",    "Chisq", "Chi Df", "Pr(>Chisq)" )
rownames(add)=NULL
knitr::kable(add, caption = "Improvement in model fit due to each added block effect")


## ---- echo=FALSE---------------------------------------------------------
library(blocksdesign)
library(lme4)
data(durban) 
durban=durban[c(3,1,2,4,5)]
# re-orders in plot order
durban = durban[ do.call(order, durban), ]
Reps = factor(rep(1:2,each=272))
Rows = factor(rep(1:16,each=34))
Col1 = factor(rep(rep(1:4,c(9,8,8,9)),16))
Col2 = factor(rep(rep(1:8,c(5,4,4,4,4,4,4,5)),16))
Col3 = factor(rep(1:34,16))
Model1= lmer(yield ~ gen + (1|Reps) + (1|Rows) + (1|Col3), durban)
Model2= lmer(yield ~ gen + (1|Reps) + (1|Rows)  + (1|Col1)  + (1|Rows:Col1) + (1|Col2) + (1|Col3), durban)
ylim=.6
coplot(resid(Model1)~bed|factor(row), data=durban, cex=.5,ylim=c(-ylim, ylim),ylab="Residuals",
       panel=function(x,y,...) panel.smooth(x,y,span=.75,col.smooth = "black",...))
title("Fig 4. Original row-and-column model")
coplot(resid(Model2)~bed|factor(row), data=durban, cex=.5,ylim=c(-ylim, ylim),ylab="Residuals",
       panel=function(x,y,...) panel.smooth(x,y,span=.75,col.smooth = "black",...))
title("Fig 5. Best fitting row-and-column model")

## ---- echo=FALSE---------------------------------------------------------
a1=c("Reps",  1,       1 ,1, 1)
a2=c("Reps+Rows", 15,0.9648,0.9648,0.9648)
a3=c(" Reps+Rows+Col1", 18,  0.9605, 0.9603, 0.9573)
a4=c("Reps+Rows+Col1+Col2", 22,  0.9507,0.9506,0.9476)
a5=c("Reps+Rows+Col1+Col2+Col3", 48,  0.8875,0.8872, 0.887)
add=data.frame(rbind(a1,a2,a3,a4,a5))
colnames(add)=c("Additive model","DF","weight=0","weight=.5","weight=1")
rownames(add)=NULL
m1=c( "Reps" ,   1,       1, 1, 1)
m2=c( "Reps.Rows", 15,0.9648,0.9648,0.9648)
m3=c( "Reps.Rows.Col1",  63,   0.841, 0.8441, 0.8442)
m4=c( "Reps.Rows.Col1.Col2", 127,  0.6749, 0.6784, 0.6785)
m5=c( "Reps.Rows.Col1.Col2.Col3", 543, "NA","NA","NA")
mult=data.frame(rbind(m1,m2,m3,m4,m5))
colnames(mult)=c("Multiplicative model","DF","weight=0","weight=.5","weight=1")
rownames(mult)=NULL
knitr::kable(add, caption = "Efficiency factors for additive effects of crossed blocks model")
knitr::kable(mult, caption = "Efficiency factors for multiplicative effects of crossed blocks model")

