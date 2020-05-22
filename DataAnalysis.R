#######################
library(ggplot2)
library(car)
library(gridExtra)
library(piecewiseSEM)

############
### FIGURE 1

gg.corr <- ggplot(subset(dat.corr,(effsize==.07|effsize==.25|effsize==.46|effsize==.70)),aes(y=p.500,ymin=p.025,ymax=p.975,x=cquality,group=effsize,color=effsize,fill=effsize))+
	geom_line(size=0.5)+geom_ribbon(alpha=.5,size=0.25)+
	facet_grid(ctype~sample)+theme_bw()+xlab("Quality of Coding")+ylab("Correlation estimates")+ylim(-1,1)+geom_hline(yintercept=0)+ggtitle("95% confidence bands")+scale_x_continuous(breaks=c(0,.25,.50,.75,1.00),limit=c(0,1))+geom_line(aes(y=sig.05,x=cquality),color="red",size=2,alpha=.1)+
	scale_color_manual(values=c("#d53e4f","#fc8d59","#1a9850","#3288bd"))+scale_fill_manual(values=c("#d53e4f","#fc8d59","#1a9850","#3288bd"))+geom_point(aes(y=issig05),size=1.0)+
	theme(panel.grid = element_blank(),axis.text.x=element_text(angle = 90,hjust=0))+
	guides(fill=guide_legend("True\neffect\nsize"),color=guide_legend("True\neffect\nsize"))
			
FIGURE.1 <- gg.corr
			
ggsave(FIGURE.1,file="FIGURE-1.svg",dpi=1200,units="cm",width=22,height=12,scale=1.4)
	

############
### FIGURE 2

gg.corr.spur <- ggplot(subset(dat.corr,ceffsize==0.07))+geom_line(aes(y=p.975,x=cquality))+geom_hline(yintercept=0)+
				geom_line(aes(y=delta.upper,x=cquality),color="#d53e4f")+geom_line(aes(y=sig.05,x=cquality),color="red",size=2,alpha=.5)+
				facet_grid(type~sample)+theme_light()+theme(panel.grid = element_blank(),axis.text.x=element_text(angle = 90,hjust=0))+ylim(-0.1,1)+
				xlab("Coding accuracy")+ylab("Upper 2.5% of correlation estimates, relative to critical value for statistical significance")

FIGURE.2 <- gg.corr.spur

ggsave(FIGURE.2,file="FIGURE-2.svg",dpi=1200,units="cm",width=22,height=14,scale=1.4)

############
### FIGURE 3

gg.corr.s.pred <- ggplot(dat.corr.s,aes(y=pre.500,ymin=pre.025,ymax=pre.975,x=cquality,group=effsize,color=effsize,fill=effsize))+geom_ribbon(alpha=.2)+
	facet_grid(.~sample)+theme_light()+ylim(-1,1)+geom_hline(yintercept=0)+ggtitle("Scale-based guessing, from equation")+ylab("Correlation estimate")+xlab("Coding accuracy")+
	scale_fill_brewer(palette="Spectral")+scale_color_brewer(palette="Spectral")+theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 90,hjust=0))
gg.corr.s.simu <- ggplot(dat.corr.s,aes(y=p.500,ymin=p.025,ymax=p.975,x=cquality,group=effsize,color=effsize,fill=effsize))+geom_ribbon(alpha=.2)+
	facet_grid(.~sample)+theme_light()+ylim(-1,1)+geom_hline(yintercept=0)+ggtitle("Scale-based guessing, from simulation")+ylab("Correlation estimate")+xlab("Coding accuracy")+
	scale_fill_brewer(palette="Spectral")+scale_color_brewer(palette="Spectral")+theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 90,hjust=0))
comp.s <- grid.arrange(gg.corr.s.pred,gg.corr.s.simu,ncol=1)

FIGURE.3 <- comp.s
			
ggsave(FIGURE.3,file="FIGURE-3.svg",dpi=1200,units="cm",width=22,height=12,scale=1.4)
	
############
### FIGURE 4

gg.corr.m.pred <- ggplot(dat.corr.m,aes(y=pre.500,ymin=pre.025,ymax=pre.975,x=cquality,group=effsize,color=effsize,fill=effsize))+geom_ribbon(alpha=.2)+
	facet_grid(.~sample)+theme_light()+ylim(-1,1)+geom_hline(yintercept=0)+ggtitle("Marginal-based guessing, from equation")+ylab("Correlation estimate")+xlab("Coding accuracy")+
	scale_fill_brewer(palette="Spectral")+scale_color_brewer(palette="Spectral")+theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 90,hjust=0))
gg.corr.m.simu <- ggplot(dat.corr.m,aes(y=p.500,ymin=p.025,ymax=p.975,x=cquality,group=effsize,color=effsize,fill=effsize))+geom_ribbon(alpha=.2)+
	facet_grid(.~sample)+theme_light()+ylim(-1,1)+geom_hline(yintercept=0)+ggtitle("Marginal-based guessing, from simulation")+ylab("Correlation estimate")+xlab("Coding accuracy")+
	scale_fill_brewer(palette="Spectral")+scale_color_brewer(palette="Spectral")+theme(panel.grid = element_blank(),axis.text.x = element_text(angle = 90,hjust=0))
comp.m <- grid.arrange(gg.corr.m.pred,gg.corr.m.simu,ncol=1)

FIGURE.4 <- comp.m
			
ggsave(FIGURE.4,file="FIGURE-4.svg",dpi=1200,units="cm",width=22,height=12,scale=1.4)

############
### FIGURE 5

eq2.gg <- ggplot()+
	geom_raster(data=eq2.df,aes(fill=min.sample.size,x=quality,y=effect.size))+
	geom_text(data=eq2.df,size=3,aes(fill=min.sample.size,x=quality,y=effect.size,label=min.sample.size),color="white")+
	geom_line(data=th2.df,aes(x=threshold,y=y,group=author,color=as.character(clr)),size=2.0)+
	geom_point(data=th2.df,aes(x=threshold,y=y,group=author,color=clr),size=3.5)+
	annotate(geom="text",hjust=0,x=c(.10,.10,.10),y=c(-.04,-.015,.01),label=c("Landis & Koch / Altman","Fleiss, Levin & Paik","Krippendorff"))+
	annotate(geom="text",hjust=0,x=c(.05,.25,.45,.65,.85),y=c(-.075),color=c("#cb181d","#ec7014","#dddd33","#6baed6","#1a9850"),label=c("do not use","maybe use","probably use","use","always use"))+
	xlab("Coding Accuracy")+ylab("Effect Size")+
	scale_fill_gradient("Minimum\nsample\nsize",trans="log",breaks=c(10,100,1000,10000,50000),low="#31a354",high="#cb181d",na.value="#dddddd")+
	scale_color_manual(breaks=c("do not use","probably use","be careful","use","always use"),labels=c("use","always use","do not use","be careful","probably use"),values=c("#6baed6","#1a9850","#cb181d","#ec7014","#dddd33"))+
	theme_light()+guides(color=FALSE)+theme(legend.position="bottom")

FIGURE.5 <- eq2.gg

ggsave(FIGURE.5,file="FIGURE-5.svg",units="cm",width=16,height=16,dpi=1200,scale=1.5)

############
### FIGURE 6

eq.gg <- ggplot()+
	geom_raster(data=eq.df,aes(fill=min.sample.size,x=quality,y=effect.size))+
	geom_text(data=eq.df,size=3,aes(fill=min.sample.size,x=quality,y=effect.size,label=min.sample.size),color="white")+
	geom_line(data=th.df,aes(x=threshold,y=y,group=author,color=as.character(clr)),size=2.0)+
	geom_point(data=th.df,aes(x=threshold,y=y,group=author,color=clr),size=3.5)+
	annotate(geom="text",hjust=0,x=c(.10,.10,.10),y=c(-.04,-.015,.01),label=c("Landis & Koch / Altman","Fleiss, Levin & Paik","Krippendorff"))+
	annotate(geom="text",hjust=0,x=c(.05,.25,.45,.65,.85),y=c(-.075),color=c("#cb181d","#ec7014","#dddd33","#6baed6","#1a9850"),label=c("do not use","maybe use","probably use","use","always use"))+
	xlab("Coding Accuracy")+ylab("Effect Size")+
	scale_fill_gradient("Minimum\nsample\nsize",trans="log",breaks=c(10,100,1000,10000,50000),low="#31a354",high="#cb181d",na.value="#dddddd")+
	scale_color_manual(breaks=c("do not use","probably use","be careful","use","always use"),labels=c("use","always use","do not use","be careful","probably use"),values=c("#6baed6","#1a9850","#cb181d","#ec7014","#dddd33"))+
	theme_light()+guides(color=FALSE)+theme(legend.position="bottom")

FIGURE.6 <- eq.gg

ggsave(FIGURE.6,file="FIGURE-6.svg",units="cm",width=16,height=16,dpi=1200,scale=1.5)

### TABLE 1

m1a <- lm(p.025~poly(csampsize,6),data=dat.corr.s)
m1b <- lm(p.025~poly(cquality,3),data=dat.corr.s)
m1c <- lm(p.025~poly(ceffsize,1),data=dat.corr.s)

m2a <- lm(p.025~poly(csampsize,6)+poly(cquality,3),data=dat.corr.s)
m2b <- lm(p.025~poly(ceffsize,1)+poly(csampsize,6),data=dat.corr.s)
m2c <- lm(p.025~poly(cquality,3)+poly(ceffsize,1),data=dat.corr.s)
m3a <- lm(p.025~poly(ceffsize,1)+poly(csampsize,6)+poly(cquality,3),data=dat.corr.s)

m2d <- lm(p.025~poly(csampsize,6)*poly(cquality,3),data=dat.corr.s)
m2e <- lm(p.025~poly(ceffsize,1)*poly(csampsize,6),data=dat.corr.s)
m2f <- lm(p.025~poly(cquality,3)*poly(ceffsize,1),data=dat.corr.s)
m3b <- lm(p.025~poly(ceffsize,1)*poly(csampsize,6)*poly(cquality,3),data=dat.corr.s)

rsquared(m1a)
rsquared(m1b)
rsquared(m1c)

rsquared(m2a)
rsquared(m2b)
rsquared(m2c)
rsquared(m3a)

rsquared(m2d)
rsquared(m2e)
rsquared(m2f)
rsquared(m3b)

##### APPENDIX

############
### FIGURE A1

# Created with MS Excel
# Just for illustration purposes, not based on any data

############
### FIGURE A2

df.reli <- melt(reliability.test,id.vars="accuracy")

df.reli[,c("coef","guess")] <- t(matrix(unlist(str_split(df.reli$variable,pattern="\\.")),nrow=2))

df.reli$Coefficient <- Recode(df.reli$coef,"'gg'='Geiss´ Gamma';'bp'='Kappa Brennan-Prediger';'ka'='Alpha Krippendorff';'fk'='Kappa Fleiss';'ck'='Kappa Cohen'")
df.reli$Guess <- Recode(df.reli$guess,"'s'='Scale based coder guessing';'m'='Marginal based coder guessing'")


gg.coef <- ggplot(df.reli,aes(y=value,x=accuracy))+geom_point()+geom_abline(aes(slope=1,intercept=0))+geom_smooth(color="#d73027",se=FALSE)+ylim(0,1)+xlim(0,1)+
	facet_grid(Coefficient~Guess)+theme_light()+xlab("Coding accuracy")+ylab("Coding agreement coefficient")

FIGURE.A2 <- gg.coef

ggsave(FIGURE.A2,file="FIGURE-A2.svg",units="cm",width=12,height=20,dpi=1200,scale=1.25)

############
### FIGURE A34

eq3.gg <- ggplot()+
	geom_raster(data=eq.df,aes(fill=1+diff.sample.size/min.sample.size,x=quality,y=effect.size))+
	geom_text(data=eq.df,aes(fill=1+diff.sample.size/min.sample.size,x=quality,y=effect.size,label=ifelse((diff.sample.size/min.sample.size==4),"???",1+round(diff.sample.size/min.sample.size,2))),color="white")+
	xlab("Coding Accuracy")+ylab("Effect Size")+
	scale_fill_gradient("Minimum\nsample\nsize",breaks=c(1,2,3,4),low="#31a354",high="#cb181d",na.value="#dddddd")+
	theme_light()

FIGURE.A3 <- eq3.gg

ggsave(FIGURE.A3,file="FIGURE-A3.svg",units="cm",width=16,height=16,dpi=1200,scale=1.5)
