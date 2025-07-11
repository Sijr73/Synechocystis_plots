library(dplyr)
library(ggplot2)
library(ggtext)
#result of the Extended model
load("Extendedmodel04032024.Rdata")
#result of the Extended model low
load("Extendedmodellow04032024.Rdata")
#result of the Faizi model
load("faizi.Rdata")
#result of the Faizi model low
load("faizi_low.Rdata")
#result of the base model
#result of the base model low
load("BaseModel_results.Rdata")
#Proteomics data
load("synechodata_final.Rdata")

###base model
# Base model growth rate vs. light intensity
# Observational data
xx=c(27.5,55,110,220,440,660,880,1100)
yy=c(0.0254, 0.03882, 0.05869, 0.08112, 0.10436, 0.1038,0.09928,0.09327)
ye=c(0.00178,0.00637,0.00933,0.00961,0.00862,0.01019,0.01252,0.01092)
data=data.frame(xx,yy)
data=data  %>% mutate(percent = (xx[])/max(data$xx)*100)

# Generate plot
g1=ggplot(BaseModel,mapping=aes(x=percent))+
  geom_line(mapping=aes(y=mu),linewidth=1.5,color = "#4477AA",linetype = "solid")+
  geom_point(data=data,mapping=aes(x=percent,y=yy),size=5,shape=18,color="black")+
  geom_errorbar(data=data , aes(x =percent,
                                ymin=yy-ye, ymax=yy+ye),width=2,color="black")+
  geom_path(data=faizi,mapping=aes(x=percent,y=gr),linewidth=1.5,color="#EE6677")+
  geom_path(data=BaseModelLow,mapping=aes(x=percent,y=mu),linewidth=1.5,color="#4477AA",linetype="dashed")+
  geom_path(data=faizi_low,mapping=aes(x=percent,y=gr),linewidth=1.5,color="#EE6677",linetype="dashed")+

  scale_x_continuous(breaks=seq(from=0, to=100,by=10),
                     labels=seq(from=0, to=100,by=10),expand=c(0,0),limits =c(0,102) )+
  scale_y_continuous(limits =c(0,max(BaseModel$mu+0.015)),labels=seq(from=0, to=0.2,by=0.02),
                     breaks=seq(from=0, to=0.2,by=0.02),expand=c(0,0))+
  geom_vline(xintercept = 20, linetype="dotted", color = "#666666", size=1.5)+
  geom_vline(xintercept = 50, linetype="dotted", color = "#666666", size=1.5)+

  annotate("text", x = 10, y = 0.02,label = "(I)",size=7,fontface="bold",family="sans"  )+
  annotate("text", x = 35, y = 0.02,label = "(II)",size=7,fontface="bold",family="sans"  )+
  annotate("text", x = 75, y = 0.02,label = "(III)",size=7,fontface="bold",family="sans" )+

  labs(
    y = bquote(bold("Growth rate") ~ "[h"^-1 * "]"),
    x =  bquote(bold("Light Intensity")~ "[% max]")) +

  theme_bw()+

  theme(

    text = element_text(size=18,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b =
                                                                                            0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x=element_text(size=18,family="sans",face="bold"),
    axis.text.y=element_text(size=18),

    panel.border = element_rect(colour = "black" ,linewidth = 1.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = c(0.7, 0.25),legend.box.background = element_rect(colour = "black",linewidth = 1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank() ,
                                                                                                               panel.background = element_blank(),
                                                                                                               axis.line = element_line(color = '#666666',linewidth = 0.8)
                                                                                                               ,axis.text.x = element_text(colour="#666666"), axis.text.y = element_text(colour="#666666"),
                                                                                                               axis.ticks =element_line(color="#666666",linewidth =1),
                                                                                                               axis.ticks.length=unit(c(-0.3,0.3), "cm"))

g1
ggsave(g1, filename="FigGrowthBase.pdf", path="~/plot/", width=210, height=297/2, units="mm")
###Proteome of base model
# Proteome of base model: Ribosome mass fraction vs. growth rate
xx=c(0.0254, 0.03882, 0.05869, 0.08112, 0.10436, 0.09327)
xe=c(0.00178,0.00637,0.00933,0.00961,0.00862,0.01092)

data1=data.frame(xx,xe)
data1=cbind(data1,syn_Ribosomefinal)
data1$Min=data1$Mean-data1$Min
data1$Max=data1$Max-data1$Mean

g2=ggplot(BaseModel,mapping=aes(x=mu))+ 
  geom_path(mapping=aes(y=phi.Ribosome,color = "High external inorganic carbon",linetype="High external inorganic carbon"),linewidth=1.5)+

  geom_point(data=data1,mapping=aes(x=xx,y=Mean),size=5,shape=18,color="black")+
  geom_path(data=data1, aes(x=xx,y = Mean),color="black",linetype="dashed",linewidth = 1)+
  geom_errorbarh(data=data1, aes(xmin=xx - xe, xmax=xx + xe, y=Mean), color="black",inherit.aes = FALSE) +
  geom_errorbar(data=data1 , aes(x =xx,
                                 ymin=Mean-Min, ymax=Mean+Max),color="black",inherit.aes = FALSE)+  
  geom_path(data=faizi,mapping=aes(x=gr,y=phi_R),linewidth=1.5,color="#EE6677")+

  scale_x_continuous(breaks=seq(from=0, to=0.2,by=0.015),
                     labels=seq(from=0, to=0.2,by=0.015),expand=c(0,0),limits =c(0,0.115) )+
  scale_y_continuous(limits = c(0, max(BaseModel$phi.Ribosome)+0.05),expand=c(0,0),labels=seq(from=0, to=0.6,by=0.1),
                     breaks=seq(from=0, to=0.6,by=0.1))+

  
  labs(
    y = bquote(bold("Proteome Mass Fraction")~ "[%]"),
    x = bquote(bold("Growth rate") ~ "[h"^-1 * "]")) +
  scale_color_manual(name="Ribosome Unit",values = c("Low external inorganic carbon"="#EE6677","High external inorganic carbon"="#4477AA"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  scale_linetype_manual(name="Ribosome Unit",
                        values = c("Low external inorganic carbon"="solid","High external inorganic carbon"="solid"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  theme_bw()+
  
  theme(
    
    text = element_text(size=20,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x=element_text(size=18,family="sans",face="bold"),
    axis.text.y=element_text(size=18),
    axis.title.y.right = element_text(angle = 90,color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 5)),
    
    panel.border = element_rect(colour = "black" ,linewidth = 1.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = c(0.3, 0.8),legend.box.background = element_rect(colour = "black",linewidth = 1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank() ,
                                                                                                              panel.background = element_blank(),
                                                                                                              axis.line = element_line(color = '#666666',linewidth = 0.8)
                                                                                                              ,axis.text.x = element_text(colour="#666666"), axis.text.y = element_text(colour="#666666"),
                                                                                                              axis.ticks =element_line(color="#666666",linewidth =1),
                                                                                                              axis.ticks.length=unit(c(-0.3,0.3), "cm"))


g2

ggsave(g2, filename="FigRibosomebase.pdf", path="~/plot/", width=210, height=297/2, units="mm")


# Proteome of base model: PSU mass fraction vs. growth rate
data1=data.frame(xx,xe)
data1$Meansum=syn_ATpasefinal$Mean+syn_Cytb6ffinal$Mean+syn_PSIfinal$Mean+syn_PSIIfinal$Mean
data1$Minsum=syn_ATpasefinal$Min+syn_Cytb6ffinal$Min+syn_PSIfinal$Min+syn_PSIIfinal$Min
data1$Maxsum=syn_ATpasefinal$Max+syn_Cytb6ffinal$Max+syn_PSIfinal$Max+syn_PSIIfinal$Max
data1$Min=data1$Meansum-data1$Minsum
data1$Max=data1$Maxsum-data1$Meansum
g3 <- ggplot(BaseModel,mapping=aes(x=mu))+ 
  geom_path(mapping=aes(y=phi.PSU,color = "High external inorganic carbon",linetype="High external inorganic carbon"),linewidth=1.5)+
  geom_path(data=faizi,mapping=aes(x=gr,y=phi_P),linewidth=1.5,color="#EE6677")+
  geom_point(data=data1,mapping=aes(x=xx,y=Meansum),size=5,shape=18,color="black")+
  geom_path(data=data1, aes(x=xx,y = Meansum),color="black",linetype="dashed",linewidth = 1)+
  geom_errorbarh(data=data1, aes(xmin=xx - xe, xmax=xx + xe, y=Meansum), color="black",inherit.aes = FALSE) +
  geom_errorbar(data=data1 , aes(x =xx,
                                 ymin=Meansum-Min, ymax=Meansum+Max),color="black",inherit.aes = FALSE)+  

  scale_x_continuous(breaks=seq(from=0, to=0.2,by=0.015),
                     labels=seq(from=0, to=0.2,by=0.015),expand=c(0,0),limits =c(0,0.115) )+
  scale_y_continuous(limits = c(0, max(BaseModel$phi.PSU)+0.05),expand=c(0,0),labels=seq(from=0, to=1,by=0.2),
                     breaks=seq(from=0, to=1,by=0.2))+

  
  labs(
    y = bquote(bold("Proteome Mass Fraction")~ "[%]"),
    x = bquote(bold("Growth rate") ~ "[h"^-1 * "]")) +
  scale_color_manual(name="Photosyntethic Unit",values = c("Low external inorganic carbon"="#CCBB44","High external inorganic carbon"="#4477AA"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  scale_linetype_manual(name="Photosyntethic Unit",
                        values = c("Low external inorganic carbon"="solid","High external inorganic carbon"="solid"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  theme_bw()+
  
  theme(
    
    text = element_text(size=20,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x=element_text(size=18,family="sans",face="bold"),
    axis.text.y=element_text(size=18),
    axis.title.y.right = element_text(angle = 90,color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 5)),
    
    panel.border = element_rect(colour = "black" ,linewidth = 1.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = c(0.3, 0.15),legend.box.background = element_rect(colour = "black",linewidth = 1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank() ,
                                                                                                               panel.background = element_blank(),
                                                                                                               axis.line = element_line(color = '#666666',linewidth = 0.8)
                                                                                                               ,axis.text.x = element_text(colour="#666666"), axis.text.y = element_text(colour="#666666"),
                                                                                                               axis.ticks =element_line(color="#666666",linewidth =1),
                                                                                                               axis.ticks.length=unit(c(-0.3,0.3), "cm"))


g3

ggsave(g3, filename="FigPSUbase.pdf", path="~/plot/", width=210, height=297/2, units="mm")


# Proteome of base model: S mass fraction vs. growth rate
g4=ggplot(BaseModel,mapping=aes(x=mu))+ 
  geom_path(mapping=aes(y=phi.Carbon,color = "High external inorganic carbon",linetype= "High external inorganic carbon"),linewidth=1.5)+
  #geom_path(mapping=aes(y=phi.Carboncarb,color = "Low external inorganic carbon",linetype="Low external inorganic carbon"),linewidth=1.5)+
  geom_path(data=faizi,mapping=aes(x=gr,y=phi_T,color = "Simulation by Faizi et. al, 2018 model",linetype="Simulation by Faizi et. al, 2018 model"),linewidth=1.5)+
  
  #geom_line(aes(y=phi.PSUcarb,color="Increasing external inorganic carbon",linetype="Increasing external inorganic carbon"),linewidth=1.5)+
  #  ylim(0,max(carbonincrease$mu_highlight))+
  scale_x_continuous(breaks=seq(from=0, to=0.2,by=0.015),
                     labels=seq(from=0, to=0.2,by=0.015),expand=c(0,0),limits =c(0,0.112) )+
  scale_y_continuous(limits =c(0,max(BaseModel$phi.Carbon)+0.01),labels=seq(from=0, to=0.5,by=0.01),
                     breaks=seq(from=0, to=0.5,by=0.01),expand = c(0,0))+
  
  
  labs(
    y = bquote(bold("Proteome Mass Fraction")~ "[%]"),
    x = bquote(bold("Growth rate") ~ "[h"^-1 * "]")) +
  
  scale_color_manual(name="Carbon transporter Unit",values = c("Low external inorganic carbon"="#228833","High external inorganic carbon"="#4477AA","Simulation by Faizi et. al, 2018 model"="#EE6677"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  scale_linetype_manual(name="Carbon transporter Unit",
                        values = c("Low external inorganic carbon"="solid","High external inorganic carbon"="solid","Simulation by Faizi et. al, 2018 model"="solid"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  theme_bw()+
  
  theme(
    
    text = element_text(size=20,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x=element_text(size=18,family="sans",face="bold"),
    axis.text.y=element_text(size=18),
    axis.title.y.right = element_text(angle = 90,color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 5)),
    
    panel.border = element_rect(colour = "black" ,linewidth = 1.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = c(0.3, 0.75),legend.box.background = element_rect(colour = "black",linewidth = 1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank() ,
                                                                                                               panel.background = element_blank(),
                                                                                                               axis.line = element_line(color = '#666666',linewidth = 0.8)
                                                                                                               ,axis.text.x = element_text(colour="#666666"), axis.text.y = element_text(colour="#666666"),
                                                                                                               axis.ticks =element_line(color="#666666",linewidth =1),
                                                                                                               axis.ticks.length=unit(c(-0.3,0.3), "cm"))



g4

ggsave(g4, filename="FigCarbonbase.pdf", path="~/plot/", width=210, height=297/2, units="mm")


# Proteome of base model: Metabolism mass fraction vs. growth rate
data1=data.frame(xx,xe)
data1$Meansum=syn_AAfinal$Mean+syn_Calvinfinal$Mean
data1$Minsum=syn_AAfinal$Min+syn_Calvinfinal$Min
data1$Maxsum=syn_AAfinal$Max+syn_Calvinfinal$Max
data1$Min=data1$Meansum-data1$Minsum
data1$Max=data1$Maxsum-data1$Meansum
g5=ggplot(BaseModel,mapping=aes(x=mu))+ 
  geom_path(mapping=aes(y=phi.Enzyme1,color = "High external inorganic carbon",linetype= "High external inorganic carbon"),linewidth=1.5)+
  geom_path(data=faizi,mapping=aes(x=gr,y=phi_M),linewidth=1.5,color="#EE6677")+
  geom_point(data=data1,mapping=aes(x=xx,y=Meansum),size=5,shape=18,color="black")+
  geom_path(data=data1, aes(x=xx,y = Meansum),color="black",linetype="dashed",linewidth = 1)+
  geom_errorbarh(data=data1, aes(xmin=xx - xe, xmax=xx + xe, y=Meansum), color="black",inherit.aes = FALSE) +
  geom_errorbar(data=data1 , aes(x =xx,
                                 ymin=Meansum-Min, ymax=Meansum+Max),color="black",inherit.aes = FALSE)+ 

  scale_x_continuous(breaks=seq(from=0, to=0.2,by=0.015),
                     labels=seq(from=0, to=0.2,by=0.015),expand=c(0,0),limits =c(0,0.115) )+
  scale_y_continuous(limits = c(0, max(BaseModel$phi.Enzyme1)+0.02),expand=c(0,0),labels=seq(from=0, to=0.2,by=0.05),
                     breaks=seq(from=0, to=0.2,by=0.05))+

  
  labs(
    y = bquote(bold("Proteome Mass Fraction")~ "[%]"),
    x = bquote(bold("Growth rate") ~ "[h"^-1 * "]")) +
  scale_color_manual(name="Metabolic enzyme Unit",values = c("Low external inorganic carbon"="#EE6677","High external inorganic carbon"="#4477AA"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  scale_linetype_manual(name="Metabolic enzyme Unit",
                        values = c("Low external inorganic carbon"="solid","High external inorganic carbon"="solid"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  theme_bw()+
  
  theme(
    
    text = element_text(size=20,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x=element_text(size=18,family="sans",face="bold"),
    axis.text.y=element_text(size=18),
    axis.title.y.right = element_text(angle = 90,color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 5)),
    
    panel.border = element_rect(colour = "black" ,linewidth = 1.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = c(0.3, 0.85),legend.box.background = element_rect(colour = "black",linewidth = 1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank() ,
                                                                                                               panel.background = element_blank(),
                                                                                                               axis.line = element_line(color = '#666666',linewidth = 0.8)
                                                                                                               ,axis.text.x = element_text(colour="#666666"), axis.text.y = element_text(colour="#666666"),
                                                                                                               axis.ticks =element_line(color="#666666",linewidth =1),
                                                                                                               axis.ticks.length=unit(c(-0.3,0.3), "cm"))


g5

ggsave(g5, filename="FigMetabolicenzymebase.pdf", path="~/plot/", width=210, height=297/2, units="mm")






####Extended model
# Extended model: growth rate with inhibition vs. light intensity
Extendedmodel=Extendedmodel04032024
Extendedmodellow=Extendedmodellow04032024

Extendedmodel  <-Extendedmodel  %>% mutate(percent = (x_PSU[])/max(Extendedmodel$x_PSU)*100)
Extendedmodellow  <-Extendedmodellow  %>% mutate(percent = (x_PSU[])/max(Extendedmodellow$x_PSU)*100)

#Faizi data points
xx=c(27.5,55,110,220,440,660,880,1100)
yy=c(0.0254, 0.03882, 0.05869, 0.08112, 0.10436, 0.1038,0.09928,0.09327)
ye=c(0.00178,0.00637,0.00933,0.00961,0.00862,0.01019,0.01252,0.01092)
data=data.frame(xx,yy,ye)
data=data  %>% mutate(percent = (xx[])/max(data$xx)*100)
g6=ggplot(Extendedmodel,mapping=aes(x=percent))+ 
  geom_path(mapping=aes(y=mu),linewidth=1.5,color = "#294766",linetype = "solid")+
  geom_point(data=data,mapping=aes(x=percent,y=yy),size=5,shape=18,color="black")+
  geom_errorbar(data=data , aes(x =percent,
                                ymin=yy-ye, ymax=yy+ye),width=2,color="black")+
  
  geom_path(data=Extendedmodellow,mapping=aes(x=percent,y=mu),linewidth=1.5,color="#294766",linetype="dashed")+
  geom_path(data=BaseModelLow,mapping=aes(x=percent,y=mu),linewidth=1.5,color="#4477AA",linetype="dashed")+
  geom_path(data=BaseModel,mapping=aes(x=percent,y=mu),linewidth=1.5,color="#4477AA",linetype="solid")+
  
  scale_x_continuous(breaks=seq(from=0, to=100,by=10),
                     labels=seq(from=0, to=100,by=10),expand=c(0,0),limits =c(0,102) )+
  scale_y_continuous(limits =c(0,max(Extendedmodel$mu)+0.03),labels=seq(from=0, to=0.2,by=0.02),
                     breaks=seq(from=0, to=0.2,by=0.02),expand=c(0,0))+
  geom_vline(xintercept = 20, linetype="dotted", color = "#666666", size=1.5)+
  geom_vline(xintercept = 50, linetype="dotted", color = "#666666", size=1.5)+
  
  annotate("text", x = 10, y = 0.02,label = "(I)",size=7,fontface="bold",family="sans"  )+
  annotate("text", x = 35, y = 0.02,label = "(II)",size=7,fontface="bold",family="sans"  )+
  annotate("text", x = 75, y = 0.02,label = "(III)",size=7,fontface="bold",family="sans" )+
  
  labs(
    y = bquote(bold("Growth rate") ~ "[h"^-1 * "]"),
    x =  bquote(bold("Light Intensity")~ "[% max]")) +

  theme_bw()+
  
  theme(
    
    text = element_text(size=18,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x=element_text(size=18,family="sans",face="bold"),
    axis.text.y=element_text(size=18),
    
    panel.border = element_rect(colour = "black" ,linewidth = 1.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = c(0.3, 0.85),legend.box.background = element_rect(colour = "black",linewidth = 1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank() ,
                                                                                                               panel.background = element_blank(),
                                                                                                               axis.line = element_line(color = '#666666',linewidth = 0.8)
                                                                                                               ,axis.text.x = element_text(colour="#666666"), axis.text.y = element_text(colour="#666666"),
                                                                                                               axis.ticks =element_line(color="#666666",linewidth =1),
                                                                                                               axis.ticks.length=unit(c(-0.3,0.3), "cm"))

g6
ggsave(g6, filename="Figgrowthextended.pdf", path="~/plot/", width=210, height=297/2, units="mm")



# Proteome of Extended model: Ribosome mass fraction vs. growth rate
xx=c(0.0254, 0.03882, 0.05869, 0.08112, 0.10436, 0.09327)
xe=c(0.00178,0.00637,0.00933,0.00961,0.00862,0.01092)
data1=data.frame(xx,xe)
data1=cbind(data1,syn_Ribosomefinal)
data1$Min=data1$Mean-data1$Min
data1$Max=data1$Max-data1$Mean

g7=ggplot(Extendedmodel,mapping=aes(x=mu))+ 
  geom_path(mapping=aes(y=phi.Ribosome,color = "High external inorganic carbon",linetype="High external inorganic carbon"),linewidth=1.5)+
  geom_point(data=data1,mapping=aes(x=xx,y=Mean),size=5,shape=18,color="black")+
  geom_path(data=data1, aes(x=xx,y = Mean),color="black",linetype="dashed",linewidth = 1)+
  geom_errorbarh(data=data1, aes(xmin=xx - xe, xmax=xx + xe, y=Mean), color="black",inherit.aes = FALSE) +
  geom_errorbar(data=data1 , aes(x =xx,
                                 ymin=Mean-Min, ymax=Mean+Max),color="black",inherit.aes = FALSE)+  
  
  geom_path(data=BaseModel, mapping=aes(y=phi.Ribosome,color = "Low external inorganic carbon",linetype="Low external inorganic carbon"),linewidth=1.5)+
  


  scale_x_continuous(breaks=seq(from=0, to=0.2,by=0.015),
                     labels=seq(from=0, to=0.2,by=0.015),expand=c(0,0),limits =c(0,0.115) )+
  
  scale_y_continuous(limits = c(0, max(BaseModel$phi.Ribosome)+0.05),expand=c(0,0),labels=seq(from=0, to=0.5,by=0.1),
                     breaks=seq(from=0, to=0.5,by=0.1))+

  
  labs(
    y = bquote(bold("Proteome Mass Fraction")),
    x = bquote(bold("Growth rate") ~ "[h"^-1 * "]")) +
  scale_color_manual(name="Ribosome Unit",values = c("Low external inorganic carbon"="#4477AA","High external inorganic carbon"="#294766"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  scale_linetype_manual(name="Ribosome Unit",
                        values = c("Low external inorganic carbon"="solid","High external inorganic carbon"="solid"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  theme_bw()+
  
  theme(
    
    text = element_text(size=18,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x=element_text(size=18,family="sans",face="bold"),
    axis.text.y=element_text(size=18),
    axis.title.y.right = element_text(angle = 90,color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 5)),
    panel.border = element_rect(colour = "black" ,linewidth = 1.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = c(0.3, 0.8),legend.box.background = element_rect(colour = "black",linewidth = 1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank() ,
                                                                                                              panel.background = element_blank(),
                                                                                                              axis.line = element_line(color = '#666666',linewidth = 0.8)
                                                                                                              ,axis.text.x = element_text(colour="#666666"), axis.text.y = element_text(colour="#666666"),
                                                                                                              axis.ticks =element_line(color="#666666",linewidth =1),
                                                                                                              axis.ticks.length=unit(c(-0.3,0.3), "cm"))

g7

ggsave(g7, filename="Ribosomeextended.pdf", path="~/plot/", width=210, height=297/2, units="mm")

# Proteome of Extended model: PSI mass fraction vs. growth rate

#PSI unit+PSIcyc+FNR
data1=data.frame(xx,xe)
data1=cbind(data1,syn_PSIfinal)
data1$Min=data1$Mean-data1$Min
data1$Max=data1$Max-data1$Mean
g8=ggplot(Extendedmodel,mapping=aes(x=mu))+ 
  geom_path(mapping=aes(y=phi.PSI,color = "High external inorganic carbon",linetype="High external inorganic carbon"),linewidth=1.5)+
  geom_point(data=data1,mapping=aes(x=xx,y=Mean),size=5,shape=18,color="black")+
  geom_path(data=data1, aes(x=xx,y = Mean),color="black",linetype="dashed",linewidth = 1)+
  geom_errorbarh(data=data1, aes(xmin=xx - xe, xmax=xx + xe, y=Mean), color="black",inherit.aes = FALSE) +
  geom_errorbar(data=data1 , aes(x =xx,
                                 ymin=Mean-Min, ymax=Mean+Max),color="black",inherit.aes = FALSE)+  
  scale_x_continuous(breaks=seq(from=0, to=0.2,by=0.015),
                     labels=seq(from=0, to=0.2,by=0.015),expand=c(0,0),limits =c(0,0.115) )+
  scale_y_continuous(limits = c(0, max(Extendedmodel$phi.PSI)+0.04),expand=c(0,0),labels=seq(from=0, to=0.4,by=0.05),
                     breaks=seq(from=0, to=0.4,by=0.05))+

  
  labs(
    y = bquote(bold("Proteome Mass Fraction")),
    x = bquote(bold("Growth rate") ~ "[h"^-1 * "]")) +
  scale_color_manual(name="Photosystem I",values = c("Low external inorganic carbon"="#757272","High external inorganic carbon"="#294766"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  scale_linetype_manual(name="Photosystem I",
                        values = c("Low external inorganic carbon"="solid","High external inorganic carbon"="solid"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  theme_bw()+
  
  theme(
    
    text = element_text(size=18,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x=element_text(size=18,family="sans",face="bold"),
    axis.text.y=element_text(size=18),
    axis.title.y.right = element_text(angle = 90,color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 5)),
    
    panel.border = element_rect(colour = "black" ,linewidth = 1.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = c(0.25, 0.9),legend.box.background = element_rect(colour = "black",linewidth = 1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank() ,
                                                                                                              panel.background = element_blank(),
                                                                                                              axis.line = element_line(color = '#666666',linewidth = 0.8)
                                                                                                              ,axis.text.x = element_text(colour="#666666"), axis.text.y = element_text(colour="#666666"),
                                                                                                              axis.ticks =element_line(color="#666666",linewidth =1),
                                                                                                              axis.ticks.length=unit(c(-0.3,0.3), "cm"))

g8

ggsave(g8, filename="PSIExtended.pdf", path="~/plot/", width=210, height=297/2, units="mm")


# Proteome of Extended model: PSII mass fraction vs. growth rate

data1=data.frame(xx,xe)
data1=cbind(data1,syn_PSIIfinal)
data1$Min=data1$Mean-data1$Min
data1$Max=data1$Max-data1$Mean

g9=ggplot(Extendedmodel,mapping=aes(x=mu))+ 
  geom_path(mapping=aes(y=phi.PSII,color = "High external inorganic carbon",linetype="High external inorganic carbon"),linewidth=1.5)+
  geom_point(data=data1,mapping=aes(x=xx,y=Mean),size=5,shape=18,color="black")+
  geom_path(data=data1, aes(x=xx,y = Mean),color="black",linetype="dashed",linewidth = 1)+
  geom_errorbarh(data=data1, aes(xmin=xx - xe, xmax=xx + xe, y=Mean), color="black",inherit.aes = FALSE) +
  geom_errorbar(data=data1 , aes(x =xx,
                                 ymin=Mean-Min, ymax=Mean+Max),color="black",inherit.aes = FALSE)+ 
  scale_x_continuous(breaks=seq(from=0, to=0.2,by=0.015),
                     labels=seq(from=0, to=0.2,by=0.015),expand=c(0,0),limits =c(0,0.115) )+
  scale_y_continuous(limits = c(0, max(Extendedmodel$phi.PSII) +0.18),expand=c(0,0),labels=seq(from=0, to=1,by=0.2),
                     breaks=seq(from=0, to=1,by=0.2))+

  
  labs(
    y = bquote(bold("Proteome Mass Fraction")~ "[%]"),
    x = bquote(bold("Growth rate") ~ "[h"^-1 * "]")) +
  scale_color_manual(name="Photosystem II",values = c("Low external inorganic carbon"="#757272","High external inorganic carbon"="#294766"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  scale_linetype_manual(name="Photosystem II",
                        values = c("Low external inorganic carbon"="solid","High external inorganic carbon"="solid"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  theme_bw()+
  
  theme(
    
    text = element_text(size=18,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x=element_text(size=18,family="sans",face="bold"),
    axis.title.y.right = element_text(angle = 90,color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 5)),
    axis.text.y=element_text(size=18),
    
    panel.border = element_rect(colour = "black" ,linewidth = 1.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = c(0.25, 0.9),legend.box.background = element_rect(colour = "black",linewidth = 1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank() ,
                                                                                                                panel.background = element_blank(),
                                                                                                                axis.line = element_line(color = '#666666',linewidth = 0.8)
                                                                                                                ,axis.text.x = element_text(colour="#666666"), axis.text.y = element_text(colour="#666666"),
                                                                                                                axis.ticks =element_line(color="#666666",linewidth =1),
                                                                                                                axis.ticks.length=unit(c(-0.3,0.3), "cm"))

g9

ggsave(g9, filename="PSIIExtended.pdf", path="~/plot/", width=210, height=297/2, units="mm")


# Proteome of Extended model: ATPsyn mass fraction vs. growth rate
data1=data.frame(xx,xe)
data1=cbind(data1,syn_ATpasefinal)
data1$Min=data1$Mean-data1$Min
data1$Max=data1$Max-data1$Mean

g10=ggplot(Extendedmodel,mapping=aes(x=mu))+ 
  geom_path(mapping=aes(y=phi.ATPase,color = "High external inorganic carbon",linetype="High external inorganic carbon"),linewidth=1.5)+
  geom_point(data=data1,mapping=aes(x=xx,y=Mean),size=5,shape=18,color="black")+
  geom_path(data=data1, aes(x=xx,y = Mean),color="black",linetype="dashed",linewidth = 1)+
  geom_errorbarh(data=data1, aes(xmin=xx - xe, xmax=xx + xe, y=Mean), color="black",inherit.aes = FALSE) +
  geom_errorbar(data=data1 , aes(x =xx,
                                 ymin=Mean-Min, ymax=Mean+Max),color="black",inherit.aes = FALSE)+
  scale_x_continuous(breaks=seq(from=0, to=0.2,by=0.015),
                     labels=seq(from=0, to=0.2,by=0.015),expand=c(0,0),limits =c(0,0.115) )+
  scale_y_continuous(limits = c(0, max(Extendedmodel$phi.ATPase)+0.01),expand=c(0,0),labels=seq(from=0, to=0.05,by=0.01),
                     breaks=seq(from=0, to=0.05,by=0.01))+

  
  labs(
    y = bquote(bold("Proteome Mass Fraction")),
    x = bquote(bold("Growth rate") ~ "[h"^-1 * "]")) +
  scale_color_manual(name="ATPsynthase unit",values = c("Low external inorganic carbon"="#757272","High external inorganic carbon"="#294766"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  scale_linetype_manual(name="ATPsynthase unit",
                        values = c("Low external inorganic carbon"="solid","High external inorganic carbon"="solid"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  theme_bw()+
  
  theme(
    
    text = element_text(size=18,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x=element_text(size=18,family="sans",face="bold"),
    axis.text.y=element_text(size=18),
    axis.title.y.right = element_text(angle = 90,color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 5)),
    
    panel.border = element_rect(colour = "black" ,linewidth = 1.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = c(0.3, 0.9),legend.box.background = element_rect(colour = "black",linewidth = 1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank() ,
                                                                                                              panel.background = element_blank(),
                                                                                                              axis.line = element_line(color = '#666666',linewidth = 0.8)
                                                                                                              ,axis.text.x = element_text(colour="#666666"), axis.text.y = element_text(colour="#666666"),
                                                                                                              axis.ticks =element_line(color="#666666",linewidth =1),
                                                                                                              axis.ticks.length=unit(c(-0.3,0.3), "cm"))

g10

ggsave(g10, filename="ATPsynextended.pdf", path="~/plot/", width=210, height=297/2, units="mm")

# Proteome of Extended model: Cytb6f mass fraction vs. growth rate

#Cytb6f+NDH
data1=data.frame(xx,xe)
data1=cbind(data1,syn_Cytb6ffinal)
data1$Min=data1$Mean-data1$Min
data1$Max=data1$Max-data1$Mean

g11=ggplot(Extendedmodel,mapping=aes(x=mu))+ 
  geom_path(mapping=aes(y=phi.Cyt6bf,color = "High external inorganic carbon",linetype="High external inorganic carbon"),linewidth=1.5)+
  geom_point(data=data1,mapping=aes(x=xx,y=Mean),size=5,shape=18,color="black")+
  geom_path(data=data1, aes(x=xx,y = Mean),color="black",linetype="dashed",linewidth = 1)+
  geom_errorbarh(data=data1, aes(xmin=xx - xe, xmax=xx + xe, y=Mean), color="black",inherit.aes = FALSE) +
  geom_errorbar(data=data1 , aes(x =xx,
                                 ymin=Mean-Min, ymax=Mean+Max),color="black",inherit.aes = FALSE)+
  scale_x_continuous(breaks=seq(from=0, to=0.2,by=0.015),
                     labels=seq(from=0, to=0.2,by=0.015),expand=c(0,0),limits =c(0,0.115) )+

  scale_y_continuous(limits =c(0,max(Extendedmodel$phi.Cyt6bf)+0.004),labels=seq(from=0, to=0.03,by=0.002),
                     breaks=seq(from=0, to=0.03,by=0.002),expand = c(0,0))+

  
  labs(
    y = bquote(bold("Proteome Mass Fraction")),
    x = bquote(bold("Growth rate") ~ "[h"^-1 * "]")) +
  scale_color_manual(name="Cytochrome b6 unit",values = c("Low external inorganic carbon"="#757272","High external inorganic carbon"="#294766"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  scale_linetype_manual(name="Cytochrome b6 unit",
                        values = c("Low external inorganic carbon"="solid","High external inorganic carbon"="solid"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  theme_bw()+
  
  theme(
    
    text = element_text(size=18,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x=element_text(size=18,family="sans",face="bold"),
    axis.text.y=element_text(size=18),
    axis.title.y.right = element_text(angle = 90,color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 5)),
    
    panel.border = element_rect(colour = "black" ,linewidth = 1.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = c(0.32, 0.9),legend.box.background = element_rect(colour = "black",linewidth = 1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank() ,
                                                                                                               panel.background = element_blank(),
                                                                                                               axis.line = element_line(color = '#666666',linewidth = 0.8)
                                                                                                               ,axis.text.x = element_text(colour="#666666"), axis.text.y = element_text(colour="#666666"),
                                                                                                               axis.ticks =element_line(color="#666666",linewidth =1),
                                                                                                               axis.ticks.length=unit(c(-0.3,0.3), "cm"))

g11

ggsave(g11, filename="Cytb6extended.pdf", path="~/plot/", width=210, height=297/2, units="mm")


# Proteome of Extended model: Carbon_T mass fraction vs. growth rate

g12=ggplot(Extendedmodel,mapping=aes(x=mu))+ 
  geom_path(mapping=aes(y=phi.Carbon,color = "High external inorganic carbon",linetype="High external inorganic carbon"),linewidth=1.5)+
  geom_path(data=BaseModel, mapping=aes(y=phi.Carbon,color = "Low external inorganic carbon",linetype="Low external inorganic carbon"),linewidth=1.5)+
  
  
  scale_x_continuous(breaks=seq(from=0, to=0.2,by=0.015),
                     labels=seq(from=0, to=0.2,by=0.015),expand=c(0,0),limits =c(0,0.112) )+
  scale_y_continuous(limits =c(0,max(BaseModel$phi.Carbon)+0.002),labels=seq(from=0, to=0.04,by=0.005),
                     breaks=seq(from=0, to=0.04,by=0.005),expand = c(0,0))+
  
  
  labs(
    y = bquote(bold("Proteome Mass Fraction")),
    x = bquote(bold("Growth rate") ~ "[h"^-1 * "]")) +
  scale_color_manual(name="Carbon transporter Unit",values = c("Low external inorganic carbon"="#4477AA","High external inorganic carbon"="#294766"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  scale_linetype_manual(name="Carbon transporter Unit",
                        values = c("Low external inorganic carbon"="solid","High external inorganic carbon"="solid"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  theme_bw()+
  
  theme(
    
    text = element_text(size=18,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x=element_text(size=18,family="sans",face="bold"),
    axis.text.y=element_text(size=18),
    
    panel.border = element_rect(colour = "black" ,linewidth = 1.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = c(0.3, 0.8),legend.box.background = element_rect(colour = "black",linewidth = 1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank() ,
                                                                                                              panel.background = element_blank(),
                                                                                                              axis.line = element_line(color = '#666666',linewidth = 0.8)
                                                                                                              ,axis.text.x = element_text(colour="#666666"), axis.text.y = element_text(colour="#666666"),
                                                                                                              axis.ticks =element_line(color="#666666",linewidth =1),
                                                                                                              axis.ticks.length=unit(c(-0.3,0.3), "cm"))

g12

ggsave(g12, filename="Carbonextended.pdf", path="~/plot/", width=210, height=297/2, units="mm")


# Proteome of Extended model: C_fix mass fraction vs. growth rate

data1=data.frame(xx,xe)
data1=cbind(data1,syn_Calvinfinal)
data1$Min=data1$Mean-data1$Min
data1$Max=data1$Max-data1$Mean

g13=ggplot(Extendedmodel,mapping=aes(x=mu))+ 
  geom_path(mapping=aes(y=phi.Rubisco,color = "High external inorganic carbon",linetype="High external inorganic carbon"),linewidth=1.5)+
  geom_point(data=data1,mapping=aes(x=xx,y=Mean),size=5,shape=18,color="black")+
  geom_path(data=data1, aes(x=xx,y = Mean),color="black",linetype="dashed",linewidth = 1)+
  geom_errorbarh(data=data1, aes(xmin=xx - xe, xmax=xx + xe, y=Mean), color="black",inherit.aes = FALSE) +
  geom_errorbar(data=data1 , aes(x =xx,
                                 ymin=Mean-Min, ymax=Mean+Max),color="black",inherit.aes = FALSE)+
  scale_x_continuous(breaks=seq(from=0, to=0.2,by=0.015),
                     labels=seq(from=0, to=0.2,by=0.015),expand=c(0,0),limits =c(0,0.115) )+
  scale_y_continuous(limits = c(0, max(Extendedmodel$phi.Rubisco)+0.01),expand=c(0,0),labels=seq(from=0, to=0.1,by=0.01),
                     breaks=seq(from=0, to=0.1,by=0.01))+

  
  labs(
    y = bquote(bold("Proteome Mass Fraction")),
    x = bquote(bold("Growth rate") ~ "[h"^-1 * "]")) +
  scale_color_manual(name="Carbon fixation Unit",values = c("Low external inorganic carbon"="#757272","High external inorganic carbon"="#294766"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  scale_linetype_manual(name="Carbon fixation Unit",
                        values = c("Low external inorganic carbon"="solid","High external inorganic carbon"="solid"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  theme_bw()+
  
  theme(
    
    text = element_text(size=18,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x=element_text(size=18,family="sans",face="bold"),
    axis.text.y=element_text(size=18),
    axis.title.y.right = element_text(angle = 90,color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 5)),
    
    panel.border = element_rect(colour = "black" ,linewidth = 1.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = c(0.3, 0.9),legend.box.background = element_rect(colour = "black",linewidth = 1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank() ,
                                                                                                              panel.background = element_blank(),
                                                                                                              axis.line = element_line(color = '#666666',linewidth = 0.8)
                                                                                                              ,axis.text.x = element_text(colour="#666666"), axis.text.y = element_text(colour="#666666"),
                                                                                                              axis.ticks =element_line(color="#666666",linewidth =1),
                                                                                                              axis.ticks.length=unit(c(-0.3,0.3), "cm"))

g13

ggsave(g13, filename="CarbonfixationExtended.pdf", path="~/plot/", width=210, height=297/2, units="mm")


# Proteome of Extended model: Metabolism mass fraction vs. growth rate

data1=data.frame(xx,xe)
data1=cbind(data1,syn_AAfinal)
data1$Min=data1$Mean-data1$Min
data1$Max=data1$Max-data1$Mean

g14=ggplot(Extendedmodel,mapping=aes(x=mu))+ 
  geom_path(mapping=aes(y=phi.AA_s,color = "High external inorganic carbon",linetype="High external inorganic carbon"),linewidth=1.5)+
  geom_point(data=data1,mapping=aes(x=xx,y=Mean),size=5,shape=18,color="black")+
  geom_path(data=data1, aes(x=xx,y = Mean),color="black",linetype="dashed",linewidth = 1)+
  geom_errorbarh(data=data1, aes(xmin=xx - xe, xmax=xx + xe, y=Mean), color="black",inherit.aes = FALSE) +
  geom_errorbar(data=data1 , aes(x =xx,
                                 ymin=Mean-Min, ymax=Mean+Max),color="black",inherit.aes = FALSE)+  
  
  scale_x_continuous(breaks=seq(from=0, to=0.2,by=0.015),
                     labels=seq(from=0, to=0.2,by=0.015),expand=c(0,0),limits =c(0,0.115) )+
  scale_y_continuous(limits = c(0, max(Extendedmodel$phi.AA_s)+0.02),expand=c(0,0),labels=seq(from=0, to=0.5,by=0.02),
                     breaks=seq(from=0, to=0.5,by=0.02))+
  
  
  labs(
    y = bquote(bold("Proteome Mass Fraction")),
    x = bquote(bold("Growth rate") ~ "[h"^-1 * "]")) +
  scale_color_manual(name="Metabolism",values = c("Low external inorganic carbon"="#757272","High external inorganic carbon"="#294766"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  scale_linetype_manual(name="Metabolism",
                        values = c("Low external inorganic carbon"="solid","High external inorganic carbon"="solid"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  theme_bw()+
  
  theme(
    
    text = element_text(size=18,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x=element_text(size=18,family="sans",face="bold"),
    axis.text.y=element_text(size=18),
    
    axis.title.y.right = element_text(angle = 90,color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 5)),
    
    panel.border = element_rect(colour = "black" ,linewidth = 1.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = c(0.3, 0.9),legend.box.background = element_rect(colour = "black",linewidth = 1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank() ,
                                                                                                              panel.background = element_blank(),
                                                                                                              axis.line = element_line(color = '#666666',linewidth = 0.8)
                                                                                                              ,axis.text.x = element_text(colour="#666666"), axis.text.y = element_text(colour="#666666"),
                                                                                                              axis.ticks =element_line(color="#666666",linewidth =1),
                                                                                                              axis.ticks.length=unit(c(-0.3,0.3), "cm"))

g14

ggsave(g14, filename="MetabolismExtended.pdf", path="~/plot/", width=210, height=297/2, units="mm")


##Total Protein concentration, Extended model
xx=c(0.0254, 0.03882, 0.05869, 0.08112, 0.10436, 0.09327)
yy=c(155,130.52,126.07,105.3,99.37,87.53)
data=data.frame(xx,yy)

faizi=faizi  %>% mutate(totalprotein = (total_prot[]*1e+16)/(2.24*6.022*1e+23))
mycolor <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")
mycolor <- c("#294766")

g15=ggplot(Extendedmodel,mapping=aes(x=mu))+ 
  geom_path(mapping=aes(y=PROTEIN,color =x_PSU,linetype="High external inorganic carbon"),linewidth=1.5)+
  
  geom_point(data=data, aes(x=xx,y = yy),size=5,shape=18,color="black")+
  geom_path(data=data, aes(x=xx,y = yy),color="black",linetype="dashed",linewidth = 1)+
  
  scale_x_continuous(breaks=seq(from=0, to=0.2,by=0.015),
                     labels=seq(from=0, to=0.2,by=0.015),expand=c(0,0),limits =c(0,0.115) )+
  scale_y_continuous(limits = c(0, max(Extendedmodel$PROTEIN)+45),expand=c(0,0),labels=seq(from=0, to=200,by=30),
                     breaks=seq(from=0, to=200,by=30))+
  
  
  labs(
    y = bquote(bold("Total Proteins")~ "[g L"^-1 * "]"),
    x = bquote(bold("Growth rate") ~ "[h"^-1 * "]")) +
  scale_color_gradientn(name="x_PSU", colors=mycolor) +
  scale_linetype_manual(name="Metabolism",
                        values = c("Low external inorganic carbon"="solid","High external inorganic carbon"="solid"),guide = guide_legend(keywidth = 0, keyheight = 0),label=NULL)+
  theme_bw()+
  
  theme(
    
    text = element_text(size=18,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text.x=element_text(size=18,family="sans",face="bold"),
    axis.text.y=element_text(size=18),
    
    axis.title.y.right = element_text(angle = 90,color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 5)),
    
    panel.border = element_rect(colour = "black" ,linewidth = 1.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 16, face = "bold",margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = "none",legend.box.background = element_rect(colour = "black",linewidth = 1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank() ,
                                                                                                         panel.background = element_blank(),
                                                                                                         axis.line = element_line(color = '#666666',linewidth = 0.8)
                                                                                                         ,axis.text.x = element_text(colour="#666666"), axis.text.y = element_text(colour="#666666"),
                                                                                                         axis.ticks =element_line(color="#666666",linewidth =1),
                                                                                                         axis.ticks.length=unit(c(-0.3,0.3), "cm"))

g15

ggsave(g15, filename="Totalproteinextended.pdf", path="~/plot/", width=210, height=297/2, units="mm")
