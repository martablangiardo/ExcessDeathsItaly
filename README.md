# R Code for the paper "Estimating weekly excess mortality at sub-national level in Italy during the COVID-19 pandemic"
by Marta Blangiardo, Michela Cameletti, Monica Pirani, Gianni Corsetti, Marco Battaglini, Gianluca Baio
https://www.medrxiv.org/content/10.1101/2020.06.08.20125211v2

Last version: 13/08/2020

----------

The following data are required for running the model:
- ./Data/MacroRegions.Rdata
- ./Data/mortality_temperature_data.Rdata
- ./Data/p20_att.txt
- ./Graphs/centro.graph
- ./Graphs/lombardia.graph
- ./Graphs/nordest.graph
- ./Graphs/nordovest.graph
- ./Graphs/sud.graph


For saving the outputs, create the following folders:
- ./Output/Females
- ./Output/Males

Requested R-packages: INLA, ggplot2, dyplyr, gridExtra.
The file make.functions.R contains additional functions written specifically for the case study.


The code contained in **model.run.R** has to be run separately for each gender ("Females" or "Males") and area ("NordOvest","Lombardia","NordEst","Centro","Sud"). 
In particular, it: 
1) prepares the data;
2) estimates the model by using r-inla and saves the output (creating files as ./Output/Sex/outputarea.Rdata, where Sex is "Females" or "Males" and area is one among "NordOvest","Lombardia","NordEst","Centro","Sud");
3) simulates from the posterior distributions by using the **make.posteriors** function (creating files as ./Output/Sex/posteriorsarea.Rdata, where Sex is "Females" or "Males" and area is one among "NordOvest","Lombardia","NordEst","Centro","Sud");
4) computes the predictions for 2020 by using the **make.predictions** function (creating files as ./Output/Sex/predictionsarea.Rdata, where Sex is "Females" or "Males" and area is one among "NordOvest","Lombardia","NordEst","Centro","Sud");
5) combines all the outcomes to produce the results for Italy. This creates the file ./Output/Sex/predItaly.Rdata" where Sex is "Females" or "Males";
6) produces plots saved as pdf files.

