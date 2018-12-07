setwd("~/doc/cours/UPMC/Master/S2/stage/")

# acquisition des données
flore<-read.csv("flore_14.csv",h=T)
flore$plot=as.factor(flore$plot)
levels(flore$plot) <- c("ctrl", "K", "P", "N","PK","NK","NP","NPK","clo_ctrl","clo_NPK")
flore[,87]=apply(flore[,41:60],1,sum)
names(flore)[87]="litiere"
flore=flore[,c(2,3,5:8,9,10:40,61:80,82:87)]

bm<-read.csv("bm_mai_14_0.csv",h=T)
bm$plot=as.factor(bm$plot)
levels(bm$plot) <- c("ctrl", "K", "P", "N","PK","NK","NP","NPK","clo_ctrl","clo_NPK")
bm$bm=bm$bm*4/29 # passe de BM a productivité (*4 pour 1 metre)(/32 = par jours)

# Analyse des communautés
#########################

## ACP ##
library(ade4)
#avril
acpavril<-dudi.pca(flore[1:30,8:64],scale=F,scannf=F,nf=5)
pvpavril=100*acpavril$eig/sum(acpavril$eig);pvpavril;
cumsum(pvpavril)
x11();s.class(acpavril$li,as.factor(flore[1:30,]$plot),cell = 1.5, axesell = F, csta = 1,col=as.numeric(flore[1:30,]$plot),grid=F)#;s.arrow(acpavril$co,boxes=F,add.plot=T)
#dev.print(pdf,"acpavril.pdf") # graph enregistrer
#jpeg(filename = "acpa.jpeg",width=600,height=600,bg = "white")
#x11();s.class(acpavril$li,as.factor(flore[1:30,]$plot),xax=1,yax=3,col=as.numeric(flore[1:30]$plot));s.arrow(acpavril$co,xax=1,yax=3,boxes=F,add.plot=T)
#x11();s.class(acpavril$li,as.factor(flore[1:30,]$plot),xax=2,yax=3,col=as.numeric(flore[1:30]$plot));s.arrow(acpavril$co,xax=2,yax=3,boxes=F,add.plot=T)
# la cloturation deplace vers le haut
# la fertilisation vers le bas et la droite
# effet de fertilisation prime sur effet cloture

#x11();s.class(acpavril$li,as.factor(flore[1:30,]$plot),cell = 0, axesell = F, csta = 1,col=as.numeric(flore[1:30,]$plot),grid=F)
#x11();s.chull(acpavril$li,as.factor(flore[1:30,]$plot), cpoi = 1,col=as.numeric(flore[1:30,]$plot),grid=F)

# Rajouter le % de chq axes
# representation des individus et non des parametre
# un pts = un releve

#mai
acpmai<-dudi.pca(flore[31:60,8:64],scale=F,scannf=F,nf=3)
pvpmai=100*acpmai$eig/sum(acpmai$eig);pvpmai
cumsum(pvpmai)
x11();s.class(acpmai$li,as.factor(flore[31:60,]$plot),cell = 1.5, axesell = F, csta = 1,col=as.numeric(flore[31:60,]$plot),grid=F)#;s.arrow(acpavril$co,boxes=F,add.plot=T)
#dev.print(pdf,"acpmai.pdf")
#x11();s.class(acpmai$li,as.factor(flore[31:60,]$plot),xax=1,yax=3,col=as.numeric(flore[31:60,]$plot));s.arrow(acpmai$co,xax=1,yax=3,boxes=F,add.plot=T)
#x11();s.class(acpmai$li,as.factor(flore[31:60,]$plot),xax=2,yax=3,col=as.numeric(flore[31:60,]$plot));s.arrow(acpmai$co,xax=2,yax=3,boxes=F,add.plot=T)
# effet de fertilisation deplace vers haut et gauche
# cloture vers le bas
# fertilisation remonte


# diversité
#############

library(multcomp)
library(nlme)
library(vegan)

div=data.frame(block=flore[,1],plot=flore[,2],periode=flore[,7],shannon=diversity(flore[,9:63]),nbsp=specnumber(flore[,9:63]))
div$pielou = div$shannon/log(div$nbsp)

# SHANNON
lme1 = lme(shannon ~ plot,rand=~ 1|block/periode,data=div) # lme = modele lineaire mixte. car donnees possiblement non inedependante entre block et periode. on considere ici un effet aleatoire de la periode dans le block
shapiro.test(residuals(lme1)) # normal
bartlett.test(residuals(lme1), div$plot) # homoscedasticite
anova(lme1) # effet du traitement sur shannon; effet ni di au block ni à la periode ( par contre possible effet du block et de la periode mais non tester ici)
summary(glht(lme1,linfct=mcp(plot="Tukey")))
grp1 = cld(glht(lme1,linfct=mcp(plot="Tukey")));grp1 #  NP different de ctrl
grpbardiv=grp1$mcletters$Letters
m=aggregate(diversity(flore[,9:63]),by=list(flore$plot),mean)
x11();bardiv=barplot(m$x,ylim=c(0,2),names.arg=m$Group.1,xlab="Traitements",ylab="Indice de Shannon moyen")
sd=aggregate(diversity(flore[,9:63]),by=list(flore$plot),sd)$x/sqrt(6) # /sqrt(6) car on veut passer de l.ecart type (standard deviation) a l.erreur type (standard error). Pour cela, il faut diviser par le nombre d'observation par traitement soit 6 (avril et mai cf flore$plot)
segments(bardiv,m$x,bardiv,m$x+sd);segments(bardiv,m$x,bardiv,m$x-sd)
text(bardiv,m$x+sd+0.1,grpbardiv)
#title("Evolution de l'indice de Shannon moyen en fonction du traitement")
#dev.print(pdf,"shannon.pdf")
# Mais :
# Richesse specifique
lme2 = lme(nbsp ~ plot,rand=~ 1|block/periode,data=div)
anova(lme2)
shapiro.test(residuals(lme2)) # normal
bartlett.test(residuals(lme2), div$plot) #homo
summary(glht(lme2,linfct=mcp(plot="Tukey")))
m=aggregate(div$nbsp,by=list(div$plot),mean)
x11();bardiv=barplot(m$x,ylim=c(0,15),names.arg=m$Group.1,xlab="Traitements",ylab="Richesse spécifique moyenne")
sd=aggregate(div$nbsp,by=list(div$plot),sd)$x/sqrt(6)
segments(bardiv,m$x,bardiv,m$x+sd);segments(bardiv,m$x,bardiv,m$x-sd)
#title("Evolution de la richesse spécifique moyenne en fonction du traitement")
#dev.print(pdf,"richsp.pdf")
# pas d'effet du traitement sur la richesse sp
# PIELOU
lme3 = lme(pielou ~ plot,rand=~ 1|block/periode,data=div) # lme = modele lineaire mixte ? car donnees possiblement non inedependante entre block et periode
shapiro.test(residuals(lme3)) # normal
bartlett.test(residuals(lme3), div$plot) # HETERO ..............
anova(lme3) # effet du traitement sur Pielou
summary(glht(lme3,linfct=mcp(plot="Tukey")))
grp3 = cld(glht(lme3,linfct=mcp(plot="Tukey")));grp3 #  NP different de ctrl
grpbardiv=grp3$mcletters$Letters
m=aggregate(div$pielou,by=list(div$plot),mean)
x11();bardiv=barplot(m$x,ylim=c(0,1),names.arg=m$Group.1,xlab="Traitements",ylab="Indice de Pielou moyen")
sd=aggregate(div$pielou,by=list(div$plot),sd)$x/sqrt(6)
segments(bardiv,m$x,bardiv,m$x+sd);segments(bardiv,m$x,bardiv,m$x-sd)
text(bardiv,m$x+sd+0.1,grpbardiv)
#title("Evolution de l'indice de Pielou moyen en fonction du traitement")
#dev.print(pdf,"pielou.pdf")
# NP < Ctrl  equitabilité diminue avec fertilisation ->  sp dominante encore plus dominante
# Pas d'effet de clo_NPK ?? Tilman 92
# effet du fertilisant sur lherbivore ?


# PRODUCTIVITE
##############

library(nlme)
library(multcomp)

anova(lm(bm$bm ~ bm$plot)) # effet du traitement sur la BM
anova(lme(bm ~ plot,rand=~ 1|block,data=bm)) # effet du traitement, pas d'effet block (ordonné à l'origine par traitement non differente entre block) difference toujours significative
shapiro.test(residuals(lme(bm ~ plot,rand=~ 1|block,data=bm))) # normal
bartlett.test(residuals(lme(bm ~ plot,rand=~ 1|block,data=bm)), bm$plot) # homo
prod = lme(bm ~ plot,rand=~ 1|block,data=bm) # lme = modele lineaire mixte ?
summary(glht(prod,linfct=mcp(plot="Tukey")))
#prodaov = aov(bm ~ plot,data=bm)
#TukeyHSD(prodaov)
plot(glht(prod,linfct=mcp(plot="Tukey")))
grp = cld(glht(prod,linfct=mcp(plot="Tukey"))) # lettre de difference par barre
grpbar=grp$mcletters$Letters # lettre de difference par barre
mbm = aggregate(bm$bm,by=list( bm$plot),mean)  # moyenne des bm par traitement
sdbm= aggregate(bm$bm,by=list( bm$plot),sd)$x/sqrt(3)
x11();p=barplot(mbm$x,ylim = c(0,6),ylab="Productivité moyenne (g/m²/j)",xlab="Traitements",names.arg=mbm$Group.1)
segments(p,mbm$x,p,mbm$x+sdbm)
segments(p,mbm$x,p,mbm$x-sdbm)
text(p,mbm$x+sdbm+0.2,grpbar)
#dev.print(pdf,"proddiv.pdf")
#title("Evolution de la productivité aérienne entre Avril et Mai en fonction des traitements")
# productivite augmente avec NP, NPK ( N et P sont co limitant)
# clo_ctrl: pas d'effet des herbivores sur la prod
# Mais clo_NPK devrait alors etre tres forte
# dynamique chaotique tilman 92 ?

# pas de  difference entre cloture & cloture+fertilisant  -> pas d'effet de la fertilisation quand les herbivores sont exclu
# difference entre controle & controle+NPK -> effet de la fertilisation en présence des herbivores


div$prod = bm$bm
anova(lm(div$prod~ div$plot+div$nbsp))
anova(lm(div$prod~ div$plot+div$shannon))
plot(div$prod~div$shannon)
anova(lm(div$prod~ div$plot+div$pielou))
plot(div$prod~div$pielou)



# Litiere
##########

summary(lm(formula = flore$litiere ~ div$nbsp)) # pas d'effet
summary(lm(flore$litiere~ div$shannon)) # pas d'effet
summary(lm(flore$litiere~ div$pielou)) # pas d'effet
anova(lm(flore$litiere~ flore$plot)) # effet
summary(lm(flore$litiere~ flore$plot)) # effet de NK, NP, NPK   pas de N ... synergie seulement ?
anova(lm(flore[flore$periode=="avr",]$litiere~flore[flore$periode=="avr",]$plot )) # en Avril effet du traitement sur la litiere
anova(lm(flore[flore$periode=="mai",]$litiere~flore[flore$periode=="mai",]$plot )) # en Mai pas d'effet

lit = lme(litiere ~ plot,rand=~ 1|block,data=flore);anova(lit) # lme = modele lineaire mixte ? 
shapiro.test(residuals(lit)) # normal
bartlett.test(residuals(lit), flore$plot) # HETERO ........;
summary(glht(lit,linfct=mcp(plot="Tukey")))
grplit = cld(glht(lit,linfct=mcp(plot="Tukey"))) # lettre de difference par barre
grpbar=grplit$mcletters$Letters # lettre de difference par barre
mbmlit = aggregate(flore$litiere,by=list( flore$plot),mean)  # moyenne des bm par traitement
sdbmlit= aggregate(flore$litiere,by=list( flore$plot),sd)$x/sqrt(6)
x11();p=barplot(mbmlit$x,ylim = c(0,70),ylab="Litière moyenne (couvert absolu en %)",xlab="Traitements",names.arg=mbmlit$Group.1)
segments(p,mbmlit$x,p,mbmlit$x+sdbmlit)
segments(p,mbmlit$x,p,mbmlit$x-sdbmlit)
text(p,mbmlit$x+sdbmlit+2,grpbar)
#dev.print(pdf,"littraitavrilmai.pdf")

# difference entre cloture & cloture+fertilisant  -> effet de la fertilisation quand les herbivores sont exclu
# pas de difference entre controle & controle+NPK -> pas d'effet de la fertilisation en présence des herbivores

################ AVRIL
floreav=flore[1:30,]
lit = lme(litiere ~ plot,rand=~ 1|block,data=floreav);anova(lit)
shapiro.test(residuals(lit)) # normal
bartlett.test(residuals(lit), floreav$plot) # homo
summary(glht(lit,linfct=mcp(plot="Tukey")))
grplit = cld(glht(lit,linfct=mcp(plot="Tukey"))) # lettre de difference par barre
grpbar=grplit$mcletters$Letters # lettre de difference par barre
mbmlit = aggregate(floreav$litiere,by=list( floreav$plot),mean)
sdbmlit= aggregate(floreav$litiere,by=list( floreav$plot),sd)$x/sqrt(3)
x11();p=barplot(mbmlit$x,ylim = c(0,100),ylab="Litière moyenne en Avril (couvert absolu en %)",xlab="Traitements",names.arg=mbmlit$Group.1)
segments(p,mbmlit$x,p,mbmlit$x+sdbmlit)
segments(p,mbmlit$x,p,mbmlit$x-sdbmlit)
text(p,mbmlit$x+sdbmlit+2,grpbar);
#dev.print(pdf,"littraitavril.pdf")
################


summary(lm((bm$bm)~ flore[flore$periode=="avr",]$litiere))  # effet de la litiere sur la productivite
shapiro.test(residuals(lm((bm$bm)~ flore[flore$periode=="avr",]$litiere))) # normal
bartlett.test(residuals(lm((bm$bm)~ flore[flore$periode=="avr",]$litiere)), bm$plot) # homo

# plus il y a de litiere et moins la productivite est importante
# plus de litiere = moins de surface disponible au sol

# litiere avril vs produtivité
x11();plot(bm$bm~ flore[flore$periode=="avr",]$litiere ,log="",ylab="Productivité moyenne (g/m²/j)",xlab="Litière moyenne en Avril (couvert absolu en %)") # plus il y a de litiere, moins la productivite est forte
abline(lm((bm$bm)~ flore[flore$periode=="avr",]$litiere))
#dev.print(pdf,"prodlit.pdf")


summary(lm(flore[flore$periode=="avr",]$litiere~flore[flore$periode=="mai",]$litiere ))
# pas de lien entre la litiere de Mai et d'avril ... 

# indice d'abondance tres subjectif:
# comparaison entre periode difficil car :
# si bcp de gd plante on vois moins la litiere = sous estimation
# + necessite d'etre effectué par un meme observateur



####
####

aa=data.frame(bm=bm$bm,nbsp=div$nbsp[1:30],shannon=div$shannon[1:30],pielou=div$pielou[1:30])
lm(aa$bm~ aa$nbsp)
anova(lm(aa$bm~ aa$nbsp))

anova(lm(aa$bm~ aa$shannon))
x11();plot(aa$bm~ aa$shannon);abline(lm(aa$bm~ aa$shannon))

anova(lm(aa$bm~ aa$pielou))
x11();plot(aa$bm~ aa$pielou);abline(lm(aa$bm~ aa$pielou))

# effet prod sur diversité
cc=data.frame(bm=bm$bm[c(1:6,9:16,19:26,29:30)],nbsp=div$nbsp[c(31:36,39:46,49:56,59:60)],shannon=div$shannon[c(31:36,39:46,49:56,59:60)],pielou=div$pielou[c(31:36,39:46,49:56,59:60)],lit=floreav$litiere[c(1:6,9:16,19:26,29:30)])
anova(lm(cc$bm~ cc$nbsp))
#shapiro.test(residuals(lm(cc$bm~ cc$nbsp)))
anova(lm(cc$bm~ cc$shannon))
#shapiro.test(residuals(lm(cc$bm~ cc$shannon)))
anova(lm(cc$bm~ cc$pielou))
#shapiro.test(residuals(lm(cc$bm~ cc$pielou)))
plot(cc$bm~ cc$shannon)
plot(cc$bm~ cc$pielou)
plot(cc$bm~ cc$nbsp)
# pas d'effet si on retire les traitements de fertilisation avec effet
# pas d'effet de la biodiversité sur la productivité

summary(lm((div$nbsp[31:60])~ flore[flore$periode=="avr",]$litiere))


anova(lm(cc$bm~cc$lit))
x11();plot(cc$bm~cc$lit)
abline(lm(cc$bm~cc$lit))
