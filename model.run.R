# Code for the paper
# "Estimating weekly excess mortality at sub-national level in Italy during the COVID-19 pandemic"
# by Marta Blangiardo, Michela Cameletti, Monica Pirani, Gianni Corsetti, Marco Battaglini, Gianluca Baio

# Last version: 13/08/2020

library(dplyr)
library(INLA)

# Load extra functions
source("make.functions.R")

############################################################################
# 1. PREPARE DATA 
############################################################################
## Load the macro areas data
load("./Data/MacroRegions.Rdata")

# Now prepare the data
Sex = "Females" #other possible choice: Males
area = "NordOvest" #other possible choices: NordEst, Sud, Centro, Lombardia
data = make.data(macro.regions,area,Sex=Sex)
graph = paste0("./Graphs/",tolower(area),".graph")

############################################################################


############################################################################
# 2. RUN INLA & SAVE THE OUTPUT 
############################################################################

## Formula with temperature
formula = morti ~ 1 + 
  f(ID1,model="bym",graph=graph,scale.model=T,
    hyper=list(theta1=list(prior="loggamma",param=c(1,0.1)),theta2=list(prior="loggamma",param=c(1,0.1)))) +
  f(week,model="rw1",replicate=ID_prov,scale.model=TRUE,hyper=list(prec=list(prior="loggamma",param=c(1,0.1)))) +
  f(Anno,model="iid") +
  f(IDtemp,model="rw2",scale.model=TRUE,hyper=list(theta=list(prior="loggamma",param=c(1,0.1))))


# INLA SET UP
# Under Poisson uses default set up
control.family=inla.set.control.family.default()
# Defines the correct variable to offset the rates in the log-linear predictor
offset = data$E

m = inla(formula,
         data=data,
         E=offset,
         family="Poisson",
         control.family=control.family,
         verbose = TRUE,
         num.threads = round(parallel::detectCores()*.8),
         control.compute=list(config = TRUE)
)

file=paste0("Output/",Sex,"/output",area,".Rdata")
save(m,file=file)

############################################################################
# 3. GIVEN THE INLA OUTPUTS, 
# SIMULATE FROM THE POSTERIOR DISTRIBUTION & SAVE THE OUTPUT 
############################################################################
# 
make.posteriors(area,Sex)
  
############################################################################
# 4. GIVEN THE INLA OUTPUTS AND SIMULATIONS,
# COMPUTE THE PREDICTIONS FOR 2020 & SAVE THE OUTPUT 
############################################################################
make.predictions(area,Sex)
  
############################################################################
# 5. COMBINE THE OUTCOMES 
############################################################################
# Analysis of all output
# Loads all the results
require(dplyr)

pred.italy=list()
for (i in names(macro.regions)){
  file=paste0("Output/",Sex,"/predictions",i,".Rdata")
  load(file)
  pred.italy$SMR=bind_rows(pred.italy$SMR,res$SMR)
  pred.italy$prediction=bind_rows(pred.italy$prediction,res$prediction)
  rm(res)
}

####
# Removes the records that don't have values in the simulations (comuni that don't match the 2016-2019 data)
# also makes 'week' a numeric variable (not a factor)
pred.italy$SMR=pred.italy$SMR %>% filter(!is.na(Obs)) %>% mutate(week=as.numeric(as.character(week)))
pred.italy$prediction=pred.italy$prediction %>% filter(!is.na(Obs)) %>% mutate(week=as.numeric(as.character(week)))

# Aggregates by province
pred.italy$SMR_prov=(pred.italy$SMR %>% group_by(DEN_UTS,week) %>% arrange(COD_PROV,week) %>% slice(1) %>% 
                       select(COD_PROV,DEN_UTS,SIGLA,DEN_REG,COD_REG,week) %>% ungroup()) %>% 
  bind_cols(pred.italy$SMR %>% group_by(DEN_UTS,week) %>% arrange(week) %>% 
              summarise_at(vars("Obs",starts_with("sim")),funs(mean)) %>% ungroup()) %>% select(-c(week...8,DEN_UTS...7)) %>%
  rename(week=week...6,DEN_UTS=DEN_UTS...2)
pred.italy$SMR_prov=pred.italy$SMR_prov %>% 
  mutate(mean=apply(as.matrix(pred.italy$SMR_prov %>% select(contains("sim"))),1,mean,na.rm=TRUE),
         median=apply(as.matrix(pred.italy$SMR_prov %>% select(contains("sim"))),1,median,na.rm=TRUE),
         sd=apply(as.matrix(pred.italy$SMR_prov %>% select(contains("sim"))),1,sd,na.rm=TRUE),
         low=apply(as.matrix(pred.italy$SMR_prov %>% select(contains("sim"))),1,quantile,.025,na.rm=TRUE),
         upp=apply(as.matrix(pred.italy$SMR_prov %>% select(contains("sim"))),1,quantile,.975,na.rm=TRUE)) %>% 
  select(COD_PROV,DEN_UTS,SIGLA,DEN_REG,COD_REG,week,Obs,mean,sd,low,median,upp,everything())

pred.italy$prediction_prov=(pred.italy$prediction %>% group_by(DEN_UTS,week) %>% arrange(COD_PROV,week) %>% slice(1) %>% 
                              select(COD_PROV,DEN_UTS,SIGLA,DEN_REG,COD_REG,week) %>% ungroup()) %>% 
  bind_cols(pred.italy$prediction %>% group_by(DEN_UTS,week) %>% arrange(COD_PROV,week) %>% 
              summarise_at(vars("Obs",starts_with("sim")),funs(sum)) %>% ungroup()) %>% 
  select(-c(week...8,DEN_UTS...7)) %>% rename(week=week...6,DEN_UTS=DEN_UTS...2)
pred.italy$prediction_prov=pred.italy$prediction_prov %>% 
  mutate(mean=apply(as.matrix(pred.italy$prediction_prov %>% select(contains("sim"))),1,mean,na.rm=TRUE),
         median=apply(as.matrix(pred.italy$prediction_prov %>% select(contains("sim"))),1,median,na.rm=TRUE),
         sd=apply(as.matrix(pred.italy$prediction_prov %>% select(contains("sim"))),1,sd,na.rm=TRUE),
         low=apply(as.matrix(pred.italy$prediction_prov %>% select(contains("sim"))),1,quantile,.025,na.rm=TRUE),
         upp=apply(as.matrix(pred.italy$prediction_prov %>% select(contains("sim"))),1,quantile,.975,na.rm=TRUE)) %>% 
  select(COD_PROV,DEN_UTS,SIGLA,DEN_REG,COD_REG,week,Obs,mean,sd,low,median,upp,everything())

# Aggregates by region
pred.italy$prediction_reg=(pred.italy$prediction %>% group_by(DEN_REG,week) %>% arrange(COD_PROV,week) %>% slice(1) %>% 
                             select(COD_REG,DEN_REG,week) %>% ungroup()) %>% 
  bind_cols(pred.italy$prediction %>% group_by(DEN_REG,week) %>% arrange(COD_PROV,week) %>% 
              summarise_at(vars("Obs",starts_with("sim")),funs(sum)) %>% ungroup()) %>% 
  select(-c(week...5,DEN_REG...4)) %>% rename(week=week...3,DEN_REG=DEN_REG...2)
pred.italy$prediction_reg=pred.italy$prediction_reg %>% 
  mutate(mean=apply(as.matrix(pred.italy$prediction_reg %>% select(contains("sim"))),1,mean,na.rm=TRUE),
         median=apply(as.matrix(pred.italy$prediction_reg %>% select(contains("sim"))),1,median,na.rm=TRUE),
         sd=apply(as.matrix(pred.italy$prediction_reg %>% select(contains("sim"))),1,sd,na.rm=TRUE),
         low=apply(as.matrix(pred.italy$prediction_reg %>% select(contains("sim"))),1,quantile,.025,na.rm=TRUE),
         upp=apply(as.matrix(pred.italy$prediction_reg %>% select(contains("sim"))),1,quantile,.975,na.rm=TRUE)) %>% 
  select(DEN_REG,COD_REG,week,Obs,mean,sd,low,median,upp,everything())
####

# Creates excess deaths [= (est-obs)/obs]
xs=-1*sweep(as.matrix(pred.italy$prediction_prov %>% select(contains("sim"))),
            1,pred.italy$prediction_prov$Obs,FUN ="-")
xs=sweep(as.matrix(xs),1,pred.italy$prediction_prov$Obs,FUN ="/")
colnames(xs)=stringr::str_replace(colnames(xs),"sim","xs")
pred.italy$prediction_prov=pred.italy$prediction_prov %>% bind_cols(as_tibble(xs)) 
pred.italy$prediction_prov=pred.italy$prediction_prov %>% mutate(
  mean.excess=apply(as.matrix(pred.italy$prediction_prov %>% select(contains("xs"))),1,mean,na.rm=T),
  sd.excess=apply(as.matrix(pred.italy$prediction_prov %>% select(contains("xs"))),1,sd,na.rm=T),
  low.excess=apply(as.matrix(pred.italy$prediction_prov %>% select(contains("xs"))),1,quantile,.025,na.rm=T),
  upp.excess=apply(as.matrix(pred.italy$prediction_prov %>% select(contains("xs"))),1,quantile,.975,na.rm=T)
) %>% select(COD_PROV,DEN_UTS,SIGLA,DEN_REG,COD_REG,week,Obs,mean.excess,sd.excess,low.excess,upp.excess,everything())


# Creates excess deaths [= (est-obs)/obs]
xs=-1*sweep(as.matrix(pred.italy$prediction %>% select(contains("sim"))),
            1,pred.italy$prediction$Obs,FUN ="-")
xs=sweep(as.matrix(xs),1,pred.italy$prediction$Obs,FUN ="/")
colnames(xs)=stringr::str_replace(colnames(xs),"sim","xs")
pred.italy$prediction=pred.italy$prediction %>% bind_cols(as_tibble(xs)) 
pred.italy$prediction=pred.italy$prediction %>% 
  mutate(
    mean.excess=apply(as.matrix(pred.italy$prediction %>% select(contains("xs"))),1,mean,na.rm=T),
    sd.excess=apply(as.matrix(pred.italy$prediction %>% select(contains("xs"))),1,sd,na.rm=T),
    low.excess=apply(as.matrix(pred.italy$prediction %>% select(contains("xs"))),1,quantile,.025,na.rm=T),
    upp.excess=apply(as.matrix(pred.italy$prediction %>% select(contains("xs"))),1,quantile,.975,na.rm=T)
  ) %>% select(COD_PROV,DEN_UTS,SIGLA,DEN_REG,COD_REG,week,Obs,mean.excess,sd.excess,low.excess,upp.excess,everything())


# Add population from official data and extrapolation at 2020
pop.tot=as_tibble(read.table("p20_att.txt",header = TRUE))
# Aggregate by age groups
pop.tot=pop.tot %>% rename(COD_PROVCOM=PRO_COM_ATT) %>% group_by(COD_PROVCOM) %>% mutate(Total.male=sum(n_m),Total.female=sum(n_f)) %>% 
  slice(1) %>% ungroup() %>% select(COD_PROVCOM,contains("Total"))
# Add population to predictions
pred.italy$prediction=pred.italy$prediction %>% left_join(pop.tot) %>% select(COD_PROVCOM,COMUNE,Total.male,Total.female,everything()) %>% 
  mutate(mean.excess=case_when(Obs==0~NA_real_,
                               TRUE~mean.excess),
         sd.excess=case_when(Obs==0~NA_real_,
                             TRUE~sd.excess),
         low.excess=case_when(Obs==0~NA_real_,
                              TRUE~low.excess),
         upp.excess=case_when(Obs==0~NA_real_,
                              TRUE~upp.excess)
  )
save(pred.italy,file=paste0("Output/",Sex,"/predItaly.Rdata"))



############################################################################
# 6. VISUALISE THE RESULTS
############################################################################
library(gridExtra)

file=paste0("Output/",Sex,"/predItaly.Rdata")
load(file)

# Plots the predicted deaths by comune & week

# Selects comuni with the highest percentage excess deaths
comuni2plot=as.character((pred.italy$prediction %>% filter(Total.male+Total.female>=100000) %>% arrange(desc(mean.excess)) %>% 
                            select(mean.excess,COMUNE,DEN_UTS,DEN_REG) %>% group_by(COMUNE) %>% slice(1) %>% ungroup() %>% 
                            arrange(desc(mean.excess)))$COMUNE)[1:9]
comuni2plot=c("Bergamo","Piacenza","Brescia","Parma","Milano","Cremona","Pesaro","Verona","Bologna")



plot.list=lapply(comuni2plot,function(x) vis.xs.rate(pred.italy,x, sex="Females"))
pdf("excess_rate_females.pdf",width=30,height=24)
do.call(grid.arrange,c(plot.list,ncol=3))
dev.off()

plot.list2=lapply(comuni2plot,function(x) vis.rate(pred.italy,x,sex="Females", multi=TRUE))
pdf("excess_mortality_females.pdf",width=30,height=24)
do.call(grid.arrange,c(plot.list2,ncol=3))
dev.off()


pdf("WeeklyTrendProv_Females.pdf")
lemon::grid_arrange_shared_legend(
  make.map(1),make.map(2),make.map(3),make.map(4),make.map(5),make.map(6),
  make.map(7),make.map(8),make.map(9),
  make.map(10),make.map(11),make.map(12),
  make.map(13),nrow=4,ncol=4,position="bottom" 
)
dev.off()


# Overall trends
# BUT MAKE TRENDS x10,000

vis.trends(pred.italy,"Lombardia")
regions=as.character((pred.italy$prediction_reg %>% group_by(DEN_REG) %>% slice(1) %>% ungroup())$DEN_REG)
plot.list=lapply(regions,function(x) vis.trends(pred.italy,x))
pdf("total_trends.pdf",width=38,height=30)
do.call(grid.arrange,c(plot.list,ncol=4))
dev.off()

relevant=list()
relevant[[1]]=c("Piemonte","Valle d'Aosta", "Liguria")
relevant[[2]]=c("Lombardia")
relevant[[3]]=c("Emilia-Romagna", "Trentino-Alto Adige","Veneto","Friuli Venezia Giulia")
relevant[[4]]=c("Abruzzo","Lazio","Marche","Molise","Toscana","Umbria")
relevant[[5]]=c("Basilicata","Calabria","Campania","Puglia","Sardegna","Sicilia")
names(relevant)=c("NordOvest","Lombardia","NordEst","Centro","Sud")

groups=list()
groups[[1]]=c("Piemonte","Valle d'Aosta", "Liguria")
groups[[2]]="Lombardia"
groups[[3]]=c("Emilia-Romagna", "Trentino-Alto Adige","Veneto","Friuli Venezia Giulia")
groups[[4]]=c("Lazio","Marche","Toscana","Umbria")
groups[[5]]=c("Abruzzo","Basilicata","Calabria","Campania","Molise","Puglia","Sardegna","Sicilia")
names(groups)=c("North-western Italy","Lombardia","North-eastern Italy","Central Italy","Southern Italy + Islands")

plot.list=lapply(groups, function(x) vis.trends(pred.italy,x,ylim=c(500,4000)))
pdf(paste0("total_trends_macroarea_",Sex,".pdf"),width=18,height=16)
do.call(grid.arrange,c(plot.list,ncol=2))
dev.off()

grid.arrange(arrangeGrob(plot.list,ncol=2,heights=c(4,1),widths=c(2,1)))


# This will tell me the posterior probability that Lombardia is outside the observed range in week 9 to complement the text
datatemp=males$prediction_reg %>% filter(DEN_REG %in% region) %>% arrange(week) 
  datatemp=datatemp %>% 
    mutate(low=apply(as.matrix(datatemp %>% select(contains("sim"))),1,quantile,q[1]),
           upp=apply(as.matrix(datatemp %>% select(contains("sim"))),1,quantile,q[2]))
  datatemp=datatemp %>% left_join(pred$prediction %>% filter(week==1) %>% group_by(DEN_REG) %>% summarise(pop=sum(Total.male)) %>% ungroup()) %>% 
    select(DEN_REG,pop,everything())
datatemp %>% mutate(prob=apply(as.matrix(datatemp %>% select(contains("sim"))),1,function(x) sum(x/pop*scale > Obs/pop*scale)/1000)) %>% select(prob,everything())



