
#### Regression analysis Lake vs. Marine vs. Terrestrial ####
#### Written by C.Raby for J.Mallalieu
#### Written December 2019

#### Packages ####
# If these packages are not installed, use install.packages("...")
# Insert the package name into the ellipses

library(fitdistrplus) # explore distribution
library(lme4) # linear regression package
library(e1071) # skewness
library(bestNormalize) # compare transformation
library(GGally) 
library(robustlmm) # robust lmm (rlmer)
library(ggplot2) # plots
library(reshape2) # melt df
library(effects)
library(ggeffects)
library(ggpubr)
library(cowplot)
library(dplyr)

#### File import ####
# file.choose()   # Run this in order to find the file on a different computer
# Then add the file into the quotations below


IceMargins <- read.csv("E:\\Side projects\\JoeGlacierTerminus2019\\IceMarginData.csv", h=T)
#IceMargins <- read.csv("C:\\Users\\craby\\Documents\\Side projects\\JoeGlacierTerminus2019\\IceMarginData.csv", h=T)
summary(IceMargins) # check data

####--------------------------------------------------------------------------------------------------------####
#### Marine, Terrestrial, Lake comparison ####

#### Subset the data
LMT <- IceMargins[,c(1:2,9:11)]
summary(LMT)

#### Explore the dependent variable
hist(LMT$AnnualChange)    # Plot frequency
#(LMT$AnnualChange, discrete = FALSE, boot=1000) # Raw data VERY skewed

#### Remove the extreme outlier (>2000m)
LMT$AnnualChange[LMT$AnnualChange ==-2112.65] <- NA #new data frame
LMT <- LMT[complete.cases(LMT), ]

#### Check for multicollinearity - LMT ####
LMT_cor <- LMT[,c(2:3)]
LMT_cor1 <- LMT[,c(2:5)]
cor(LMT_cor) # non over |r| = 0.70

# Plot out correlations
boxplot(LMT$Latitude ~ LMT$Type)
plot(LMT$AnnualChange ~ LMT$Latitude)
TypeLat <- lm(Latitude ~ Type, data=LMT)
summary(TypeLat)

#### Scale the variables ####
# or convert into factors

LMT$Year <- scale(as.numeric(as.character(LMT$Year)))
LMT$Latitude <- scale(as.numeric(as.character(LMT$Latitude)))
LMT$Type <- factor(LMT$Type)
LMT$FID <- as.factor(LMT$FID)

#### Transformation options ####

# Cubed transformation
LMT$AC_cubed <- sign(LMT$AnnualChange) * (abs(LMT$AnnualChange)^(1/3))
hist(LMT$AC_cubed)
skewness(LMT$AC_cubed) 

# Hyperbolic transformation
ihs <- function(x) {
  y <- log(x + sqrt(x ^ 2 + 1))
  return(y)
}
LMT$AC_ihs <- ihs(LMT$AnnualChange)
hist(LMT$AC_ihs)
skewness(LMT$AC_ihs)

# Best normalize – to compare the transformation options

AC <- LMT$AnnualChange
#(BNobject <- bestNormalize(AC))

onAC <- orderNorm(AC)
LMT$OrdNormAC <- onAC$x.t


#### Compare linear regressions with different transformations ####

LMT$AnnualChange <- scale(LMT$AnnualChange)
Model_ <- lmer(AnnualChange ~  Year + Type + Latitude + (1|FID),  data = LMT)

Model_cubed <- lmer(AC_cubed ~  Year + Type + Latitude + (1|FID),  data = LMT) # isSingular

Model_ihs <- lmer(AC_ihs ~  Year + Type + Latitude + (1|FID),  data = LMT) # isSingular

Model_ordNorm<- lmer(OrdNormAC ~  Year + Type + Latitude + (1|FID),  data = LMT) # isSingular

AIC(Model_, Model_cubed, Model_ihs, Model_ordNorm)
#               df      AIC
# Model_         7 49026.80
# Model_cubed   7 69257.41
# Model_ihs     7 78796.24
# Model_ordNorm 7 50124.62

plot(Model_)
plot(Model_cubed)
plot(Model_ihs)
plot(Model_ordNorm)
# All plots of fitted vs residuals are similar

hist(resid(Model_, type="pearson"))
hist(resid(Model_cubed, type="pearson"))
hist(resid(Model_ihs, type="pearson"))
hist(resid(Model_ordNorm, type="pearson"))
# residual plots

Model_
Model_cubed
Model_ihs
Model_ordNorm
summary(Model_)
summary(Model_cubed)
summary(Model_ihs)
summary(Model_ordNorm)

drop1(Model_, test="Chi")
anova(Model_)
require(lmerTest)
Model_1 <- lmer(AnnualChange ~  Year + Type + Latitude + (1|FID),  data = LMT)
summary(Model_1)
confint(Model_1)

#### Also explore variables vs. residuals in order norm transformed model

sresid <- resid(Model_ordNorm, type="pearson")
plot(sresid, LMT$Latitude)
plot(sresid, LMT$Type)
plot(sresid, LMT$FID)

#### Main issue is heavy tailed - the non-transformed data has a lower AIC and similar plots

Model_Rob <- rlmer(AnnualChange ~  Year + Type + Latitude + (1|FID),  data = LMT)

# Results and summary
plot(Model_Rob)
summary(Model_Rob)
print(plot(Model_Rob, which = 3)[[1]] + scale_color_continuous(guide = "none"))
compare(Model_Rob, Model_, show.rho.functions = FALSE)


# Confidence intervals

confint.rlmerMod <- function(object,parm,level=0.95) {
  beta <- fixef(object)
  if (missing(parm)) parm <- names(beta)
  se <- sqrt(diag(vcov(object)))
  z <- qnorm((1+level)/2)
  ctab <- cbind(beta-z*se,beta+z*se)
  colnames(ctab) <- stats:::format.perc(c((1-level)/2,(1+level)/2),
                                        digits=3)
  return(ctab[parm,])
}
confint(Model_Rob)

#### Compare with LMM non-transformed model
as.data.frame(effect("Type", Model_))
as.data.frame(effect("Type", Model_Rob))
as.data.frame(effect("Latitude", Model_))
as.data.frame(effect("Latitude", Model_Rob))

compare(Model_Rob, Model_)


####--------------------------------------------------------------------------------------------------------####
#### Graphs  ####

#### Comparison of LMT across time
LMTplot <- IceMargins[,c(1:2,9:11)]
#### Remove the extreme outlier (>2000m)
LMTplot$AnnualChange[LMTplot$AnnualChange ==-2112.65] <- NA #new data frame
LMTplot <- LMTplot[complete.cases(LMTplot), ]
summary(LMTplot)

LMTplot %>%
  group_by(Year, Type) %>%
  summarise(mean = mean(AnnualChange), n = n())

boxplot(AnnualChange ~ Type, data=LMTplot)

LMTplot$Type <- ordered(LMTplot$Type, levels = c("Terrestrial","Lake","Marine"))
LMTplot.long<-melt(LMTplot,id.vars=c("AnnualChange","Type","Year"))


ggplot(LMTplot, aes(as.factor(Year), AnnualChange, fill=as.factor(Type))) +
  geom_boxplot() + 
  coord_cartesian(ylim = c(-300, 120), expand = TRUE) + 
  theme_bw()+
  scale_fill_brewer(palette="Accent") +
  geom_hline(yintercept = 0, colour="red", lty=5, lwd=0.9) + 
  xlab("Year") + ylab("Annual Change") + labs(fill="Margin Type")

LMTplot2 <- LMTplot
LMTplot2$Type[LMTplot2$Type == "Marine"] <- NA #new data frame
LMTplot2 <- LMTplot2[complete.cases(LMTplot2), ]

ggplot(LMTplot2, aes(as.factor(Year), AnnualChange, fill=as.factor(Type))) +
  geom_boxplot() + 
  coord_cartesian(ylim = c(-79, 39), expand = TRUE) + 
  theme_bw()+
  scale_fill_brewer(palette="Accent") +
  geom_hline(yintercept = 0, colour="red", lty=5, lwd=0.9) + 
  xlab("Year") + ylab("Annual Change") + labs(fill="Margin Type")

plot_grid(LMTplot1, LTplot1, align = "h", nrow=2)

ggplot(LMTplot2, aes(as.numeric(Year), AnnualChange, color = Type)) +
  #geom_point() +
  geom_smooth(method = "loess") +
  #coord_cartesian(ylim = c(-15, 5), expand = TRUE) + 
  theme_bw()+
  scale_fill_brewer(palette="Accent") +
  xlab("Year") + ylab("Annual Change") + labs(fill="Margin Type")


####--------------------------------------------------------------------------------------------------------####
#### Variables shaping Lake - glacier change ####

Subplot <- IceMargins[which(IceMargins$Type == "Lake"),]
head(Subplot)
Lakes.df <- Subplot[,c(1:3, 5:6, 9:10, 12)]
Lakes.df <- Lakes.df[complete.cases(Lakes.df),]

#### Scale the variables ####
hist(Subplot$AnnualChange)
Lakes.df$Latitude <- scale(Lakes.df$Latitude)
Lakes.df$Altitude <- scale(Lakes.df$Altitude)
Lakes.df$Area <- scale(Lakes.df$Area)
Lakes.df$Year <- scale(Lakes.df$Year)
Lakes.df$IntersectLength <- scale(Lakes.df$IntersectLength)

#### Check for multicollinearity - Lakes ####
summary(Lakes.df)
Lakes.df$FID <- as.numeric(Lakes.df$FID)
cor(Lakes.df)

#### Multicollinearity present: |r| Area and Intersect Length = 0.7484 ####
#### Therefore will construct two seperate models

#### Transformation options ####

# Cubed transformation
Lakes.df$AC_cubed <- sign(Lakes.df$AnnualChange) * (abs(Lakes.df$AnnualChange)^(1/3))
hist(Lakes.df$AC_cubed)
skewness(Lakes.df$AC_cubed) 

# Hyperbolic transformation
ihs <- function(x) {
  y <- log(x + sqrt(x ^ 2 + 1))
  return(y)
}
Lakes.df$AC_ihs <- ihs(Lakes.df$AnnualChange)
hist(Lakes.df$AC_ihs)
skewness(Lakes.df$AC_ihs)

# Best normalize – to compare the transformation options

Lakes.dfAC <- Lakes.df$AnnualChange
#(BNobject <- bestNormalize(AC))

Lakes.dfonAC <- orderNorm(Lakes.dfAC)
Lakes.df$OrdNormAC <- Lakes.dfonAC$x.t
hist(Lakes.df$OrdNormAC)

#### Compare models: Area ####

Lakes.df$AnnualChange <- scale(Lakes.df$AnnualChange)
Lakes.df$FID <- as.factor(Lakes.df$FID)
Lakes.df$Durability <- factor(Lakes.df$Durability, ordered=TRUE)

LakeModel_ <- lmer(AnnualChange ~ Latitude + Altitude + Area + Year + Durability + (1|FID),  data = Lakes.df) # isSingular

LakeModel_cubed <- lmer(AC_cubed ~ Latitude + Altitude + Area + Year + Durability + (1|FID),  data = Lakes.df) 

LakeModel_ihs <- lmer(AC_ihs ~ Latitude + Altitude + Area + Year + Durability + (1|FID),  data = Lakes.df) 

LakeModel_ordNorm<- lmer(OrdNormAC ~ Latitude + Altitude + Area + Year + Durability + (1|FID),  data = Lakes.df) 

AIC(LakeModel_, LakeModel_cubed, LakeModel_ihs, LakeModel_ordNorm)

LakeModel_Rob <- rlmer(AnnualChange ~ Latitude + Altitude + Area + Year + Durability + (1|FID),  data = Lakes.df)

summary(LakeModel_Rob)
summary(LakeModel_)
summary(LakeModel_ordNorm)

compare(LakeModel_, LakeModel_Rob, LakeModel_ordNorm)

confint(LakeModel_Rob)
confint(LakeModel_)
confint(LakeModel_ordNorm)


#### Compare models: IntersectLength ####
head(Lakes.df)

LakeModel2_ <- lmer(AnnualChange ~ Latitude + Altitude + IntersectLength + Year + Durability + (1|FID),  data = Lakes.df) # isSingular

LakeModel2_cubed <- lmer(AC_cubed ~ Latitude + Altitude + IntersectLength + Year + Durability + (1|FID),  data = Lakes.df) 

LakeModel2_ihs <- lmer(AC_ihs ~ Latitude + Altitude + IntersectLength + Year + Durability + (1|FID),  data = Lakes.df) 

LakeModel2_ordNorm<- lmer(OrdNormAC ~ Latitude + Altitude + IntersectLength + Year + Durability + (1|FID),  data = Lakes.df)

AIC(LakeModel2_, LakeModel2_cubed, LakeModel2_ihs, LakeModel2_ordNorm)

LakeModel2_Rob <- rlmer(AnnualChange ~ Latitude + Altitude + IntersectLength + Year + Durability + (1|FID),  data = Lakes.df)

summary(LakeModel2_Rob)
summary(LakeModel2_)
summary(LakeModel2_ordNorm)

compare(LakeModel2_, LakeModel2_Rob, LakeModel2_ordNorm)

confint(LakeModel2_Rob)
confint(LakeModel2_)
confint(LakeModel2_ordNorm)

####--------------------------------------------------------------------------------------------------------####

#### Regression Line function ####
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

#### GRAPHS ####

# Base plots to explore relationships
head(Lakes.df)
plot(AnnualChange ~ Year, data=Subplot)

ggplot(Subplot, aes(x=Latitude, y=AnnualChange)) + 
  #geom_point()+
  coord_cartesian(ylim = c(-10, 0), expand = TRUE) +
  geom_smooth(method=lm)

# Plot 1: Latitude

Latpred <- ggpredict(LakeModel_, terms = "Latitude")
#Lplot <- ggplot(Latpred, aes(x, predicted)) +
#  geom_line(colour="Blue") +
#  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
#  labs(x="Latitude", y="Annual Change") +
#   theme_bw()
Lattrend <- ggplot(Subplot, aes(x=Latitude, y=AnnualChange)) + 
  geom_smooth(method=lm)+
  theme_bw()

ggplotRegression(lm(AnnualChange ~ Latitude, data = Subplot))

# Plot 2: Altitude

Altpred <- ggpredict(LakeModel_, terms = "Altitude")
#Aplot <- ggplot(Altpred, aes(x, predicted)) +
#  geom_line(colour="Blue") +
#  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
#  labs(x="Altitude", y="Annual Change") +
#  theme_bw()
Alttrend <- ggplot(Subplot, aes(x=Altitude, y=AnnualChange)) + 
  geom_smooth(method=lm)+
  theme_bw()

ggplotRegression(lm(AnnualChange ~ Altitude, data = Subplot))

# Plot 3: Area

Areapred <- ggpredict(LakeModel_, terms = "Area")
#Arplot <- ggplot(Areapred, aes(x, predicted)) +
#  geom_line(colour="Blue") +
#  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
#  labs(x="Area", y="Annual Change") +
#  theme_bw()
Areatrend <- ggplot(Subplot, aes(x=Area, y=AnnualChange)) + 
  geom_smooth(method=lm)+
  theme_bw()

#### Units Raw, convert to km2
Subplot$Area10 <- (Subplot$Area / 1000000)
ggplotRegression(lm(AnnualChange ~ Area10, data = Subplot))
summary(Subplot$Area10)


# Plot 4: Year

Yearpred <- ggpredict(LakeModel_, terms = "Year")
#Yplot <- ggplot(Yearpred, aes(x, predicted)) +
#  geom_line(colour="Blue") +
#  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
#  labs(x="Year", y="Annual Change") +
#  theme_bw()
Yeartrend <- ggplot(Subplot, aes(x=Year, y=AnnualChange)) + 
  geom_smooth(method=lm)+
  theme_bw()
Yeartrend

ggplotRegression(lm(AnnualChange ~ Year, data = Subplot))

# Plot 5: Intersect length

Intpred <- ggpredict(LakeModel2_, terms = "IntersectLength")

#Iplot <- ggplot(Intpred, aes(x, predicted)) +
#  geom_line(colour="Blue") +
# geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
#  labs(x="Intersect Length", y="Annual Change") +
#  theme_bw()
Itrend <- ggplot(Subplot, aes(x=IntersectLength, y=AnnualChange)) + 
  geom_smooth(method=lm) +
  theme_bw()

#### Units raw, convert to km2
Subplot$IL10 <- (Subplot$IntersectLength / 1000)
ggplotRegression(lm(AnnualChange ~ IL10, data = Subplot))

### Combine plots

ggarrange(Lattrend, Alttrend, Areatrend, Yeartrend, Itrend,
          labels = c("A", "B", "C", "D", "E"),
          ncol = 2, nrow = 3)

#### Boxplot raw Annual change and Durability

boxplot(AnnualChange ~ Durability, Lakes.df)
