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
knitr::kable(add, caption = "Table 1. Efficiency factors for nested blocks model")


## ---- echo=FALSE---------------------------------------------------------

Hist1=vector("list", 6)
class(Hist1)="histogram"
Hist1$breaks=c( 0,  5, 10, 15, 20, 25, 30, 35)
Hist1$counts=c( 330, 320, 195,  96,  44,  13,   2)
Hist1$density=c(0.0660, 0.0640, 0.0390, 0.0192, 0.0088, 0.0026, 0.0004)
Hist1$mids=c( 2.5,  7.5, 12.5, 17.5, 22.5, 27.5, 32.5)
Hist1$xname="ChiSq difference"
Hist1$equidist="TRUE"
plot(Hist1, main="Fig 1. Blocks of size 8 fitted after blocks of size 4")

Hist2=vector("list", 6)
class(Hist2)="histogram"
Hist2$breaks=c( 0,  5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55)
Hist2$counts=c(483, 240, 140,  66,  36,  20,   8,   5,   0,   1,   1)
Hist2$density=c(0.0966, 0.0480, 0.0280, 0.0132, 0.0072, 0.0040, 0.0016, 0.0010, 0.0000, 0.0002, 0.0002)
Hist2$mids=c( 2.5,  7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5)
Hist2$xname="ChiSq difference"
Hist2$equidist="TRUE"
plot(Hist2, main="Fig 2. Blocks of size 4 fitted after blocks of size 8 ")

Hist3=vector("list", 6)
class(Hist3)="histogram"
Hist3$breaks=c( 0,  5, 10, 15, 20, 25, 30, 35, 40)
Hist3$counts=c( 756, 125,  69,  21,  14,  11,   3,   1)
Hist3$density3=c(0.1512, 0.0250, 0.0138, 0.0042, 0.0028, 0.0022, 0.0006, 0.0002)
Hist3$mids=c(2.5,  7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5)
Hist3$xname="ChiSq difference"
Hist3$equidist="TRUE"
plot(Hist3,  main="Fig 3. Blocks of size 2 fitted after blocks of size 4 and 8")


## ---- echo=FALSE---------------------------------------------------------
a= c("Reps + Rows + Col1 + Rows:Col1 + Col2 + Rows:Col2 + Col3", 280.9, 1484, 139.6, -279.1, 264)
b= c("Reps + Rows + Col1 + Col2 + Rows:Col2 + Col3", 285.8, 1485, 136.1, -272.2, 265)
c= c("Reps + Rows + Col2 + Rows:Col2 + Col3", 284.5, 1479, 135.7, -271.4, 266)
d= c("Reps + Rows + Col1 + Rows:Col1 + Col2 + Col3", 302.6, 1502, 127.7, -255.4, 265) 
e= c("Reps + Rows + Col1 + Rows:Col1 + Col3", 305.9, 1501,  125.0, -250.0,  266)
f = c("Reps + Rows + Col1 + Rows:Col1 + Col2 + Rows:Col2",361.2, 1560, 98.4, -196.8, 265)
add=data.frame(rbind(a,b,c,d,e,f))
colnames(add)=c("Blocks Model", "AIC", "BIC", "logL","dev","df" )
rownames(add)=NULL
knitr::kable(add, caption = "Table 2. Goodness of fit of column block model terms")


## ----message=FALSE, echo=FALSE-------------------------------------------
require(blocksdesign)
require(lme4)
data(durban) 
durban=durban[c(3,1,2,4,5)]
# re-orders in plot order
durban = durban[ do.call(order, durban), ]
Reps = factor(rep(1:2,each=272))
Rows = factor(rep(1:16,each=34))
Col1 = factor(rep(rep(1:4,c(9,8,8,9)),16))
Col2 = factor(rep(rep(1:8,c(5,4,4,4,4,4,4,5)),16))
Col3 = factor(rep(1:34,16))
Model1 = lmer(yield ~ gen + Reps + (1|Rows) + (1|Col3), durban)
Model2 = lmer(yield~  gen + Reps + (1|Rows) + (1|Col2) + (1|Rows:Col2) + (1|Col3) ,durban)
ylim=.6
coplot(resid(Model1)~bed|factor(row), data=durban, cex=.5,ylim=c(-ylim, ylim),ylab="Residuals",
       panel=function(x,y,...) panel.smooth(x,y,span=.75,col.smooth = "black",...))
title("Fig 4. Original row-and-column blocks model")
coplot(resid(Model2)~bed|factor(row), data=durban, cex=.5,ylim=c(-ylim, ylim),ylab="Residuals",
       panel=function(x,y,...) panel.smooth(x,y,span=.75,col.smooth = "black",...))
title("Fig 5. Best fitting multi-level blocks model")


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

