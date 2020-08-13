# Code for the paper
# "Estimating weekly excess mortality at sub-national level in Italy during the COVID-19 pandemic"
# by Marta Blangiardo, Michela Cameletti, Monica Pirani, Gianni Corsetti, Marco Battaglini, Gianluca Baio

# Last version: 13/08/2020

# This file contains some functions used in the Model_Run.R file
################################################################################
# Functions for preparing the data

make.data=function(macro.regions=macro.regions,area,Sex) {
  # Makes the data by subsetting the relevant area and removing the data for 2020 and creating relevant indices
  data=macro.regions[[area]] %>% filter(Anno<max(Anno),sex==Sex) %>% mutate(SIGLA=droplevels(SIGLA))
  data$ID1=(data %>% mutate(IDarea=group_indices(.,ID_Ita)))$IDarea
  data=data %>% mutate(ID_prov=group_indices(.,SIGLA)) %>% select(-c(nord,centro,sud)) %>% select(sex,everything())
  data=data %>% mutate(IDtemp=as.numeric(as.character(IDtemp)))
  return(data)
}

################################################################################
# Function for simulating from the posterior distributions

make.posteriors = function(area=c("NordOvest","Lombardia","NordEst","Centro","Sud"),
                           Sex,nsim=1000){
  # area = macro region to analyse
  # sex = Females or Males
  # nsim = number of simulations from the posterior distribution
  
  require(INLA)
  require(dplyr)
  
  # Loads the INLA object with the output
  file=paste0("Output/",Sex,"/output",area,".Rdata")
  load(file)
  
  # Data file
  data=as_tibble(m$.args$data)
  distr=m$.args$family
  
  # Now samples from the approximate joint posterior distribution of all the parameters
  tic=proc.time()
  post.samp=inla.posterior.sample(n=nsim,m,selection=list(
    "(Intercept)"=1,
    "week"=1:nrow(m$summary.random$week),
    "ID1"=1:(nrow(m$summary.random$ID1)/2), #NB Only needs the first half as there are 2 components to BYM
    "IDtemp"=1:nrow(m$summary.random$IDtemp),
    "Anno"=4
  ),
  num.threads = round(parallel::detectCores()*.8),
  verbose=FALSE
  )
  toc=proc.time()-tic
  # Formats simulations in a matrix (like BUGS would do)
  sim=matrix(unlist(lapply(post.samp,function(x) x$latent[1:nrow(post.samp[[1]]$latent),])),ncol=nrow(post.samp[[1]]$latent),byrow=T)
  colnames(sim)=rownames(post.samp[[1]]$latent)
  
  # Now simulates from the posterior of the hyperparameters (to get sd for the 'Anno' effect)
  sd.anno=sqrt(1/inla.hyperpar.sample(nsim,m,improve.marginals=TRUE)[,"Precision for Anno"])
  # Then simulates from the predictive distribution for Anno=5
  anno.2020=rnorm(nsim,0,sd.anno)
  # Then substitutes the relevant column in 'sim'
  sim[,grep("Anno",colnames(sim))]=anno.2020
  file=paste0("Output/",Sex,"/posteriors",area,".Rdata")
  save(sim,data,distr,file=file)
}

################################################################################
# Function for computing the predictions for year 2020

make.predictions=function(area=c("NordOvest","Lombardia","NordEst","Centro","Sud"),Sex){
  # area = macro region to analyse
  # sex = Females or Males
  # nsim = number of simulations from the posterior distribution
  
  require(dplyr)
  
  # Loads mortality & temperature data + data used to run INLA
  load("Data/mortality_temperature_data.Rdata")
  
  # Loads posteriors simulations
  file=paste0("Output/",Sex,"/posteriors",area,".Rdata")
  load(file)
  nsim=nrow(sim)
  
  # matrix of covariates (profiles) to compute the SMRs and predicted
  data.pred=data %>% filter(Anno==1) %>% select(-temperature,-temp_grp,-IDtemp,-morti,-E) 
  data.pred=data.pred %>% left_join(tempdata) %>% 
    mutate(ID1=as.factor(ID1),week=as.factor(week),ID_prov=as.factor(ID_prov))
  
  formula.covs=~ID1+week:ID_prov+Anno+IDtemp
  X=as_tibble(model.matrix(formula.covs,
                           data=data.pred,contrasts.arg=lapply(data.pred[,c("ID1","week","IDtemp")],contrasts,contrasts=FALSE))) %>% 
    select(contains("ID1"),contains("week"),contains("Anno"),contains("IDtemp"),`(Intercept)`)
  # NB: in Sud, the temperature is *always* higher than the first bin, so need to remove on column 
  if(area=="Sud") {
    X=X %>% select(-IDtemp1)
  }
  
  # Now computes logSMR as the linear predictor (on the log scale) = X\beta
  log.smr=as.matrix(X)%*%t(sim)
  # Renames the columns to be sim1,sim2,...,simnsim
  colnames(log.smr)=paste0("sim",1:nrow(sim))
  if(distr=="poisson") {
    SMR=bind_cols(as_tibble(exp(log.smr)),data.pred) %>% 
      select(COD_PROVCOM,COMUNE,COD_PROV,DEN_UTS,SIGLA,DEN_REG,COD_REG,week,Anno,everything()) %>% 
      mutate(Anno=5) %>% 
      left_join(morti_ordered %>% filter(sex==Sex,Anno==max(Anno)) %>% mutate(week=as.factor(week))) %>% 
      select(COD_PROVCOM,COMUNE,COD_PROV,DEN_UTS,SIGLA,DEN_REG,COD_REG,week,morti,E,everything()) %>% rename(Obs=morti)
    rm(log.smr)
    SMR=SMR %>% mutate(mean=apply(as.matrix(SMR %>% select(contains("sim"))),1,mean),
                       sd=apply(as.matrix(SMR %>% select(contains("sim"))),1,sd),
                       low=apply(as.matrix(SMR %>% select(contains("sim"))),1,quantile,.025),
                       upp=apply(as.matrix(SMR %>% select(contains("sim"))),1,quantile,.975)) 
    SMR=SMR %>% select(COD_PROVCOM,COMUNE,COD_PROV,DEN_UTS,SIGLA,DEN_REG,COD_REG,week,E,Obs,mean,sd,low,upp,contains("sim"))
    
    # Now makes matrix of simulations from the predictive distribution of outcome. First creates the linear predictor on the natural scale
    linpred=SMR$E*as.matrix(SMR %>% select(contains("sim")),nrow=nrow(SMR),byrow=TRUE)
    # Different simulation scheme, depending on whether model is ZIP or Poisson
    # Simulates from a simple Poisson with parameter linpred
    pred=t(apply(linpred,1,function(x) rpois(nsim,x)))
    colnames(pred)=paste0("sim",1:nsim)
    prediction=bind_cols(as_tibble(pred),SMR %>% select(-contains("sim")))
    prediction=prediction %>% mutate(mean=apply(as.matrix(prediction %>% select(contains("sim"))),1,mean,na.rm=T),
                                     median=apply(as.matrix(prediction %>% select(contains("sim"))),1,median,na.rm=T),
                                     sd=apply(as.matrix(prediction %>% select(contains("sim"))),1,sd,na.rm=T),
                                     low=apply(as.matrix(prediction %>% select(contains("sim"))),1,quantile,.025,na.rm=T),
                                     upp=apply(as.matrix(prediction %>% select(contains("sim"))),1,quantile,.975,na.rm=T)) %>% 
      select(COD_PROVCOM,COMUNE,COD_PROV,DEN_UTS,SIGLA,DEN_REG,COD_REG,week,Obs,E,mean,sd,low,median,upp,contains("sim"))
  }
  
  # Selects relevant regions for each macro-area (to be used to filter out relevant comuni)
  relevant=list()
  relevant[[1]]=c("Piemonte","Valle d'Aosta", "Liguria")
  relevant[[2]]=c("Lombardia")
  relevant[[3]]=c("Emilia-Romagna", "Trentino-Alto Adige","Veneto","Friuli Venezia Giulia")
  relevant[[4]]=c("Abruzzo","Lazio","Marche","Molise","Toscana","Umbria")
  relevant[[5]]=c("Basilicata","Calabria","Campania","Puglia","Sardegna","Sicilia")
  names(relevant)=c("NordOvest","Lombardia","NordEst","Centro","Sud")
  
  # Selects output of the function --- removes comuni outside the relevant regions
  res=list(SMR=SMR %>% filter(DEN_REG%in%relevant[[area]]),
           prediction=prediction %>% filter(DEN_REG%in%relevant[[area]]))
  file=paste0("Output/",Sex,"/predictions",area,".Rdata")
  save(res,file=file)
}


####################################################################################
#Functions for visualizing the results

vis.comuni=function(pred,comune,...) {
  require(ggplot2)
  rg=range((pred$prediction %>% filter(COMUNE==comune) %>% arrange(week) %>% select(mean,Obs)))
  dates=as.Date("2020-01-01")
  for (i in 2:max(pred$prediction$week)) {
    dates[i]=dates[(i-1)]+7
  }
  dates=format(dates,format="%d%b")
  theme_set(theme_bw())
  ggplot(pred$prediction %>% filter(COMUNE==comune) %>% arrange(week),aes(week,mean)) + 
    geom_line(aes(y=mean),color="blue") + 
    geom_ribbon(aes(ymin=low,ymax=upp),alpha=.2) +
    scale_x_continuous("Week", breaks=1:length(dates), labels=dates) +
    geom_point(data=pred$prediction %>% filter(COMUNE==comune) %>% arrange(week),aes(week,Obs),color="red") +
    geom_line(data=pred$prediction %>% filter(COMUNE==comune) %>% arrange(week),aes(week,Obs),color="red") +
    theme(axis.text.x = element_text(color="black",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(color="black",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(color="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.title.y = element_text(color="black",size=14,angle=90,hjust=.5,vjust=.5,face="plain")) +
    labs(y="All causes deaths",title=paste0("Municipality of ",comune," (",(pred$prediction %>% filter(COMUNE==comune))$DEN_REG,")")) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size=18)) +
    geom_point(shape=20,colour="red",size=3,aes(x=.92,y=rg[2]*1.035)) +
    annotate(geom="segment",x=.8,y=rg[2]*1.052,xend=1.04,yend=rg[2]*1.052,color="grey",size=3) +
    annotate(geom="segment",x=.8,y=rg[2]*1.048,xend=1.04,yend=rg[2]*1.048,color="grey",size=3) +
    annotate(geom="segment",x=.8,y=rg[2]*1.05,xend=1.04,yend=rg[2]*1.05,color="blue",size=.8) +
    annotate(geom="text",x=1.15,y=rg[2]*1.05,label="Predicted deaths (mean and 95% interval)",size=5,vjust=.5,hjust=0) +
    annotate(geom="text",x=1.15,y=rg[2]*1.035,label="Observed deaths",size=5,vjust=.5,hjust=0)
}

vis.rate=function(pred,comune,sex="Male",...) {
  # pred = a list with the model predictions
  # region = a string or a vector of strings with the names of the region(s) to plot
  exArgs=list(...)
  if(exists("title",exArgs)){title=exArgs$title} else {title=paste0("Municipality of ",comune," (",(pred$prediction %>% filter(COMUNE==comune))$DEN_REG,")")}
  if(exists("scale",exArgs)){scale=exArgs$scale} else {scale=100000}
  
  require(ggplot2)
  dates=as.Date("2020-01-01")
  for (i in 2:max(pred$prediction$week)) {
    dates[i]=dates[(i-1)]+7
  }
  dates=format(dates,format="%e\n%b")
  if(sex=="Male"){
    pred$prediction=pred$prediction %>% mutate(Total=Total.male)
  } else {
    pred$prediction=pred$prediction %>% mutate(Total=Total.female)
  }
  
  datatemp=pred$prediction %>% filter(COMUNE==comune) %>% arrange(week) %>% 
    mutate(rate=mean/Total*scale,rate.low=low/Total*scale,rate.upp=upp/Total*scale,rateobs=Obs/Total*scale)
  
  # Defines plot range
  rg=range(datatemp %>% select(rate.low,rate.upp,rateobs))
  
  # Make the plot
  theme_set(theme_bw())
  pl=ggplot(data=datatemp,aes(week,rate)) + 
    geom_line(aes(y=rate),color="blue") + 
    geom_ribbon(aes(ymin=rate.low,ymax=rate.upp),alpha=.2) +
    scale_x_continuous("Week", breaks=1:length(dates), labels=dates) 
  if(exists("ylim",exArgs)){
    rg=exArgs$ylim; rg[2]=rg[2]*1.055
    pl+scale_y_continuous(limits=rg)
  }
  pl+geom_point(data=datatemp,aes(week,rateobs),color="red") +
    geom_line(data=datatemp,aes(week,rateobs),color="red") +
    theme(axis.text.x = element_text(color="black",size=10,angle=0,hjust=.5,vjust=.5,face="bold", family="Times"),
          axis.text.y = element_text(color="black",size=10,angle=0,hjust=.5,vjust=.5,face="bold", family="Times"),  
          axis.title.x = element_text(color="black",size=11,angle=0,hjust=.5,vjust=.5,face="bold", family="Times"),
          axis.title.y = element_text(color="black",size=11,angle=90,hjust=.5,vjust=.5,face="bold", family="Times")) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size=10, face="bold", family="Times")) + 
    labs(y=paste0("Mortality rate for all causes (x",format(scale,big.mark=",",scientific=999),")"),title=title) +
    geom_point(shape=20,colour="red",size=3,aes(x=.92,y=rg[2]*1.01)) +
    annotate(geom="segment",x=.8,y=rg[2]*1.052,xend=1.04,yend=rg[2]*1.052,color="grey",size=3) +
    annotate(geom="segment",x=.8,y=rg[2]*1.048,xend=1.04,yend=rg[2]*1.048,color="grey",size=3) +
    annotate(geom="segment",x=.8,y=rg[2]*1.05,xend=1.04,yend=rg[2]*1.05,color="blue",size=.8) +
    annotate(geom="text",x=1.15,y=rg[2]*1.05,family="Times",
             label=paste0("Predicted mortality rate x",format(scale,big.mark=",",scientific=999)," (mean and 95% interval)"),size=3.8,vjust=.5,hjust=0) +
    annotate(geom="text",x=1.15,y=rg[2]*1.01,family="Times",
             label=paste0("Observed mortality rate x" ,format(scale,big.mark=",",scientific=999)),size=3.8,vjust=.5,hjust=0)+
    coord_cartesian(ylim = c(rg[1],rg[2]*1.05)) 
}


vis.xs.rate=function(pred,comune,...) {
  require(ggplot2)
  dates=as.Date("2020-01-01")
  for (i in 2:max(pred$prediction$week)) {
    dates[i]=dates[(i-1)]+7
  }
  dates=format(dates,format="%e\n%b")
  
  # Defines plot range
  rg=100*range((pred$prediction %>% filter(COMUNE==comune) %>% arrange(week) %>% select(low.excess,mean.excess)),na.rm=TRUE)
  
  # Make the plot
  theme_set(theme_bw())
  ggplot(pred$prediction %>% filter(COMUNE==comune) %>% arrange(week),aes(week,100*mean.excess)) + 
    geom_line(aes(y=100*mean.excess),color="blue") + 
    geom_ribbon(aes(ymin=100*low.excess,ymax=100*upp.excess),alpha=.2) +
    scale_y_continuous("Percentage excess mortality",breaks=seq(-200,100,by=50),labels=seq(-200,100,by=50),limits=c(-290,105)) +
    scale_x_continuous("Week", breaks=1:length(dates), labels=dates) +
    geom_hline(yintercept=0,size=.8, linetype="dashed") +
    geom_hline(yintercept=50,size=.8, linetype="dashed",color="red") +
    geom_hline(yintercept=100,size=.8, linetype="dashed",color="red") +
    theme(axis.text.x = element_text(color="black",size=10,angle=0,hjust=.5,vjust=.5,face="bold", family="Times"),
          axis.text.y = element_text(color="black",size=10,angle=0,hjust=.5,vjust=.5,face="bold", family="Times"),  
          axis.title.x = element_text(color="black",size=11,angle=0,hjust=.5,vjust=.5,face="bold", family="Times"),
          axis.title.y = element_text(color="black",size=11,angle=90,hjust=.5,vjust=.5,face="bold", family="Times")) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size=12, face="bold", family="Times")) +
    labs(y="Percentage excess mortality",title=paste0("Municipality of ",comune," (",(pred$prediction %>% filter(COMUNE==comune))$DEN_REG,")")) +
    annotate(geom="text",x=1,y=50,label="1.5-fold increase",size=3.8,vjust=-1,hjust=0,family="Times") +
    annotate(geom="text",x=1,y=100,label="2-fold increase",size=3.8,vjust=-1,hjust=0,family="Times") +
    coord_cartesian(ylim = c(-100,105)) 
}

plot.smr=function(pred,comune) {
  require(ggplot2)
  dates=as.Date("2020-01-01")
  for (i in 2:max(pred$prediction$week)) {
    dates[i]=dates[(i-1)]+7
  }
  dates=format(dates,format="%d%b")
  rg=range(pred$SMR %>% filter(COMUNE==comune) %>% arrange(week) %>% mutate(rate=Obs/E) %>% select(mean,rate))
  theme_set(theme_bw())
  ggplot(pred$SMR %>% filter(COMUNE==comune) %>% arrange(week),aes(week,mean)) + 
    geom_line(aes(y=mean),color="blue") + 
    geom_ribbon(aes(ymin=low,ymax=upp),alpha=.2) +
    scale_x_continuous("Week", breaks=1:length(dates), labels=dates) +
    geom_point(data=pred$SMR %>% filter(COMUNE==comune) %>% arrange(week) %>% mutate(rate=Obs/E),aes(week,rate),color="red") +
    geom_line(data=pred$SMR %>% filter(COMUNE==comune) %>% arrange(week) %>% mutate(rate=Obs/E),aes(week,rate),color="red") +
    theme(axis.text.x = element_text(color="black",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(color="black",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(color="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.title.y = element_text(color="black",size=14,angle=90,hjust=.5,vjust=.5,face="plain")) +
    labs(y="Standardised Mortality Ratio",title=paste0("Municipality of ",comune," (",(pred$prediction %>% filter(COMUNE==comune))$DEN_REG,")")) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title= element_text(size=18)) +
    geom_point(shape=20,colour="red",size=3,aes(x=1.0,y=rg[2]*1.035)) +
    annotate(geom="segment",x=.8,y=rg[2]*1.052,xend=1.04,yend=rg[2]*1.052,color="grey",size=3) +
    annotate(geom="segment",x=.8,y=rg[2]*1.048,xend=1.04,yend=rg[2]*1.048,color="grey",size=3) +
    annotate(geom="segment",x=.8,y=rg[2]*1.05,xend=1.04,yend=rg[2]*1.05,color="blue",size=.8) +
    annotate(geom="text",x=1.15,y=rg[2]*1.05,label="Predicted SMR (mean and 95% interval)",size=5,vjust=.5,hjust=0) +
    annotate(geom="text",x=1.15,y=rg[2]*1.035,label="Observed SMR",size=5,vjust=.5,hjust=0)
}



make.map=function(w) {
  # Make map using ggplot
  require(ggplot2)
  toplot=prov.shp
  toplot@data=as_tibble(prov.shp@data %>% mutate(COD_REG=as.numeric(as.character(COD_REG)),
                                                 COD_PROV=as.numeric(as.character(COD_PROV)))) %>% 
    left_join(pred.italy$prediction_prov %>% filter(week==w) %>% 
                # NB: Need to recode the COD_PROV for the 4 new provinces otherwise won't merge
                mutate(COD_PROV=case_when(
                  COD_PROV==104~108,
                  COD_PROV==105~109,
                  COD_PROV==106~110,
                  COD_PROV==107~111,
                  TRUE~as.numeric(COD_PROV)))
    )
  
  theme_set(theme_bw())
  dates=as.Date("2020-01-01")
  for (i in 2:max(pred.italy$prediction_prov$week)) {
    dates[i]=dates[(i-1)]+7
  }
  dates=format(dates,format="%d%b")
  shades=colorRampPalette(c("white","blue","red"))#c("grey100", "grey30"))
  ###shades=RColorBrewer::brewer.pal(4,"Greens")
  ggplot(data=sf::st_as_sf(toplot,coords=c(SHAPE_AREA,SHAPE_LEN))) +
    geom_sf()+
    geom_sf(data = sf::st_as_sf(toplot,coords=c(SHAPE_AREA,SHAPE_LEN)), aes(fill = mean.excess),
            ## This controls the thickness of the borders
            size=.4) +
    ### Check this: https://stackoverflow.com/questions/49497182/ggplot-create-a-border-overlay-on-top-of-map
    ### to also plot the regions borders on top of the provinces
    geom_sf(data = sf::st_as_sf(toplot,coords=c(SHAPE_AREA,SHAPE_LEN)) %>% group_by(COD_REG) %>% summarise(),            
            fill="transparent",color="gray20",size=1.1) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_fill_gradientn(
      colours=shades(20),
      name="Posterior mean of excess mortality",
      na.value = "grey100", 
      limits=c(-2.5,1.5)
      #trans = "log",
      #breaks = my_breaks, labels = my_breaks
    ) +  
    ggtitle(paste0("Week of ",dates[w])) 
  #theme(legend.position = pos)
}


vis.trends=function(pred,region,...) {
  # pred = a list with the model predictions
  # region = a string or a vector of strings with the names of the region(s) to plot
  exArgs=list(...)
  if(exists("quants",exArgs)){q=exArgs$quants} else {q=c(.025,.975)}
  if(exists("title",exArgs)){title=exArgs$title} else {title=paste0(region,collapse=", ")}
  
  require(ggplot2)
  # creates temporary dataset with only the relevant regions
  datatemp=pred$prediction_reg %>% filter(DEN_REG %in% region) %>% arrange(week) 
  datatemp=datatemp %>% 
    mutate(low=apply(as.matrix(datatemp %>% select(contains("sim"))),1,quantile,q[1]),
           upp=apply(as.matrix(datatemp %>% select(contains("sim"))),1,quantile,q[2]))
  # if more than one region, then sum up the numbers
  if(length(region)>1) {
    datatemp=datatemp %>% group_by(week) %>% 
      summarise_at(vars("Obs",starts_with("sim")),funs(sum)) %>% ungroup()
    datatemp=datatemp %>% 
      mutate(mean=apply(as.matrix(datatemp %>% select(contains("sim"))),1,mean),
             median=apply(as.matrix(datatemp %>% select(contains("sim"))),1,median),
             sd=apply(as.matrix(datatemp %>% select(contains("sim"))),1,sd),
             low=apply(as.matrix(datatemp %>% select(contains("sim"))),1,quantile,q[1]),
             upp=apply(as.matrix(datatemp %>% select(contains("sim"))),1,quantile,q[2])) %>% 
      select(week,Obs,mean,sd,low,median,upp,everything())
  }
  
  rg=range((datatemp %>% select(low,upp,Obs)))
  
  dates=as.Date("2020-01-01")
  for (i in 2:max(pred$prediction$week)) {
    dates[i]=dates[(i-1)]+7
  }
  dates=format(dates,format="%d-%b")
  
  theme_set(theme_bw())
  pl=ggplot(datatemp,aes(week,mean)) + 
    geom_line(aes(y=mean),color="blue") + 
    geom_ribbon(aes(ymin=low,ymax=upp),alpha=.2) +
    scale_x_continuous("Week", breaks=1:length(dates), labels=dates) 
  if(exists("ylim",exArgs)){
    rg=exArgs$ylim; rg[2]=rg[2]*1.055
    pl+scale_y_continuous(limits=rg)
  }
  pl + geom_point(data=datatemp,aes(week,Obs),color="red") +
    geom_line(data=datatemp,aes(week,Obs),color="red") +
    theme(axis.text.x = element_text(color="black",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.text.y = element_text(color="black",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),  
          axis.title.x = element_text(color="black",size=14,angle=0,hjust=.5,vjust=.5,face="plain"),
          axis.title.y = element_text(color="black",size=14,angle=90,hjust=.5,vjust=.5,face="plain")) +
    labs(y="All causes deaths",title=title) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size=18)) +
    geom_point(shape=20,colour="red",size=3,aes(x=.92,y=rg[2]*1.015)) +
    annotate(geom="segment",x=.8,y=rg[2]*1.052,xend=1.04,yend=rg[2]*1.052,color="grey",size=3) +
    annotate(geom="segment",x=.8,y=rg[2]*1.048,xend=1.04,yend=rg[2]*1.048,color="grey",size=3) +
    annotate(geom="segment",x=.8,y=rg[2]*1.05,xend=1.04,yend=rg[2]*1.05,color="blue",size=.8) +
    annotate(geom="text",x=1.15,y=rg[2]*1.05,label=paste0("Predicted deaths (mean and ",100*(q[2]-q[1]),"% interval)"),
             size=5,vjust=.5,hjust=0) +
    annotate(geom="text",x=1.15,y=rg[2]*1.015,label="Observed deaths",size=5,vjust=.5,hjust=0) + 
    geom_vline(xintercept=(10+5/7),linetype="dashed",color="black",size=.6) + 
    annotate(geom="text",x=(10+4.5/7),y=rg[2]*.85,label="Italy goes in \nlockdown on 9 March",size=3.8,vjust=0,hjust=1, family="Times", face="bold") + 
    # Add this to ensure multiple graphs in grid.arrange have the same scale
    # source: https://stackoverflow.com/questions/41768273/create-ggplots-with-the-same-scale-in-r
    coord_cartesian(ylim = c(rg[1],rg[2]*1.05)) 
}
