
a <-c("S1","S2","S3","S4","S5","S6","S7","S8","S9") #names
b <-c(0.23,0.26,0.55,0.56,0.36,0.23,0.18,0.06,0.04) #mean t0
c <-c(0.64,0.6,0.81,1.4,0.89,0.55,0.48,0.22,0.09) #mean t1
d <-c(0.20,0.23,0.52,0.53,0.33,0.20,0.15,0.04,0.03) #SD low t0
e <-c(0.26,0.29,0.58,.59,0.39,0.26,0.21,0.08,0.05) #SD high t0
f <-c(0.67,0.63,0.86,1.44,0.93,0.59,0.51,0.25,0.10) #SD high t1
g <-c(0.61,0.57,0.78,1.36,0.85,0.53,0.45,0.19,0.08) #SD low t1
h <-c(0.41,0.34,0.26,0.84,0.53,0.32,0.30,0.16,0.05) #absolute change
data <- data.frame(a,b,c,d,e,f,g,h)

cols <- c("LINE1"="#f04546","LINE2"="#3591d1","BAR"="#62c76b")
ggplot(data=data,aes(x=a)) + 
  geom_bar(stat="identity", aes(y=h, fill = "BAR"),colour="#333333")+ #green
  geom_line(aes(y=b,group=1, colour="LINE1"),size=1.0) +   #red
  geom_point(aes(y=b, colour="LINE1"),size=3) +           #red
  geom_errorbar(aes(ymin=d, ymax=e, colour="LINE1"), width=0.1, size=.8) + 
  geom_line(aes(y=c,group=1,colour="LINE2"),size=1.0) +   #blue 
  geom_point(aes(y=c,colour="LINE2"),size=3) +           #blue
  geom_errorbar(aes(ymin=f, ymax=g,colour="LINE2"), width=0.1, size=.8) + 
  scale_colour_manual(name="Error Bars",values=cols) + scale_fill_manual(name="Bar",values=cols) +
  ylab("Symptom severity") + xlab("PHQ-9 symptoms") +
  ylim(0,1.6) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 15, vjust=0.3))


ggplot(data=data,aes(x=a)) + 
  geom_bar(stat="identity", aes(y=h,fill = "BAR", colour="BAR"))+ #green
  geom_line(aes(y=b,group=1, colour="LINE1"),size=1.0) +   #red
  geom_point(aes(y=b, colour="LINE1", fill="LINE1"),size=3) +           #red
  geom_errorbar(aes(ymin=d, ymax=e, colour="LINE1"), width=0.1, size=.8) + 
  geom_line(aes(y=c,group=1,colour="LINE2"),size=1.0) +   #blue 
  geom_point(aes(y=c,colour="LINE2", fill="LINE2"),size=3) +           #blue
  geom_errorbar(aes(ymin=f, ymax=g,colour="LINE2"), width=0.1, size=.8) + 
  scale_colour_manual(name="Error Bars",values=cols, 
                      guide = guide_legend(override.aes=aes(fill=NA))) + 
  scale_fill_manual(name="Bar",values=cols, guide="none") +
  ylab("Symptom severity") + xlab("PHQ-9 symptoms") +
  ylim(0,1.6) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 15, vjust=0.3))


library(drc)
x=1:2
y=0:1
model = drc::drm(y~x  ,fct = LL.4())

plot(model)
points(x,y, col='red', pch =16,  cex=1.2)

