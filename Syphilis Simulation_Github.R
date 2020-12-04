#################
# Syphilis and doxycycline simulation study
# Citation: Tran NK, Goldstein ND, Welles SL. Countering the Rise of Syphilis in the Age of HIV PrEP: A Role for Doxycycline Prophylaxis? Manuscript in preparation.
# Note: Simulation datasets may be downloaded from: https://drive.google.com/file/d/19Qd6CDay0qHLa-A0FutmlRRPvVTZczwt/view?usp=sharing
# 12/4/20 -- Nguyen Tran 
#################


### FUNCTIONS ###

prevention_paradigm = function(rand_vector, prevention_var, dataset, intervention_level=NA)
{
  if (prevention_var=="TAPVL") {
    #TAPVL (viral load under TAP, 1=suppressed, 2=not suppressed)
    starting_vals = ifelse(dataset$HIV_STATUS==3 & rand_vector<=0.42, 1, ifelse(dataset$HIV_STATUS==3 & rand_vector>0.42, 2, ifelse(dataset$HIV_STATUS==2 & dataset$HIV==1, 2, -1)))
  } else if (prevention_var=="PROBCOND") {
    #PROBCOND (condom probability)
    starting_vals = ifelse(dataset$RACE==1 & rand_vector<=0.23, 1, ifelse(dataset$RACE==1 & rand_vector<=0.51, 0, ifelse(dataset$RACE==1 & rand_vector<=0.71, 0.75, ifelse(dataset$RACE==1 & rand_vector<=0.86, 0.50, ifelse(dataset$RACE==1, 0.25, ifelse(dataset$RACE==2 & rand_vector<=0.27, 1, ifelse(dataset$RACE==2 & rand_vector<=0.47, 0, ifelse(dataset$RACE==2 & rand_vector<=0.69, 0.75, ifelse(dataset$RACE==2 & rand_vector<=0.86, 0.50, ifelse(dataset$RACE==2, 0.25, NA))))))))))
  } else if (prevention_var=="SEROTYPE") {
    #SEROTYPE (1=seropositioning/insertive, 2=seropositioning/receptive, 3=serosorting, 4=no seroadaptive behavior)
    starting_vals = ifelse(dataset$HIV_STATUS==3 & rand_vector<=0.133, 2, ifelse(dataset$HIV_STATUS==3 & rand_vector<=0.347, 3, ifelse(dataset$HIV_STATUS==3 & rand_vector>0.347, 4, ifelse(dataset$HIV_STATUS==2 & rand_vector<=0.053, 1, ifelse(dataset$HIV_STATUS==2 & rand_vector<=0.374, 3, ifelse(dataset$HIV_STATUS==2 & rand_vector>0.374, 4, ifelse(dataset$HIV_STATUS==1 & rand_vector<=0.095, 1, ifelse(dataset$HIV_STATUS==1 & rand_vector<=0.402, 3, ifelse(dataset$HIV_STATUS==1 & rand_vector>0.402, 4, NA)))))))))
  } else if (prevention_var=="ONPREP") {
    #ONPREP (0=not on prep, 1=on prep)
    starting_vals = ifelse(dataset$HIV_STATUS==1 & rand_vector<=(intervention_level/((1-0.19)-((1-0.19)*0.369))), 1, 0)
  } else if (prevention_var=="HPVDOSES") {
    #HPVDOSES (number of doses of SY vaccine)
    starting_vals = ifelse(rand_vector<=(1-intervention_level), 0, ifelse(rand_vector<=(intervention_level*(4/13)+(1-intervention_level)), 1, ifelse(rand_vector<=(intervention_level*((4+2)/13)+(1-intervention_level)), 2, 3)))
  } else if (prevention_var=="ONDOXYPREP"){
    #ONDOXYPREP (0=not taking to doxy prep, 1=20% adherent, 2=40% adherent, 3=60% adherent, 4=72% adherent, 5=80% adherent ,6=100% adherent)
    # starting_vals = ifelse(rand_vector<=(1-intervention_level), 0, ifelse(rand_vector<=(intervention_level*0.2+(1-intervention_level)),1, ifelse(rand_vector<=(intervention_level*0.4+(1-intervention_level)), 2, ifelse(rand_vector<=(intervention_level*0.6+(1-intervention_level)), 3, ifelse(rand_vector<=(intervention_level*0.72+(1-intervention_level)), 4, ifelse(rand_vector<=(intervention_level*0.8+(1-intervention_level)), 5, 6))))))
    #ONDOXYPREP (0=not taking doxy prep, 1=93% adherent)
    starting_vals = ifelse(rand_vector<=(1-intervention_level), 0, ifelse(rand_vector<=(intervention_level*0.72+(1-intervention_level)),1, 0))
  } else if (prevention_var=="ONDOXYPEP"){
    #ONDOXYPEP (0=not taking to doxy pep, 1=20% adherent, 2=40% adherent, 3=60% adherent, 4=80% adherent, 5=93%; 6=100% adherent)
    # starting_vals = ifelse(rand_vector<=(1-intervention_level), 0, ifelse(rand_vector<=(intervention_level*0.2+(1-intervention_level)),1, ifelse(rand_vector<=(intervention_level*0.4+(1-intervention_level)), 2, ifelse(rand_vector<=(intervention_level*0.6+(1-intervention_level)), 3, ifelse(rand_vector<=(intervention_level*0.8+(1-intervention_level)), 4, ifelse(rand_vector<=(intervention_level*0.93+(1-intervention_level)), 5, 6))))))
    #ONDOXYPEP (0=not taking doxy pep, 1=72% adherent)
    starting_vals = ifelse(rand_vector<=(1-intervention_level), 0, ifelse(rand_vector<=(intervention_level*0.93+(1-intervention_level)),1, 0))
  }
  
  return(starting_vals)
}

#network modeling
library("EpiModel") 


### STEP0: GLOBAL and NETWORK INITIALIZATION ###

#benchmarking
#time0 = Sys.time()

#ensure that all runs start equivalently
set.seed(777)

#number of simulations
nsims=25

#number of prevention paradigm scenarios evaluated per simulation
nparadigms = 13

#ensure that each simulation is unique
seed_val_population = sample(1:(nsims*10),nsims,replace=F)

#number of agents in model/initial YMSM population size in Philadelphia
nagents_start = 10230

#number of agents in model total including migration (value set later)
nagents_end = NA

#length of simulation (in days); change year multiplier to change length
length_sim = round(365 * 20, 0)

#DEBUG code
#nsims = nparadigms = 1
#length_sim = 365*1
#sims = day = current_data = 1

#initialize container for results of all simulations
abm_sims_50 = list(NA) #abm_sims is simulation data
long_sims_50 = list(NA) #long_sims is longitudinal summary data
network_sims_50 = list(NA) #network_sims is sexual network data

#initialize the network with number of agents
nw = network.initialize(n=nagents_start, directed=F)

#add assortative mixing characteristics
#age distribution: 1= <=24, 2= >24
#nw = set.vertex.attribute(nw, attrname="age", value=ifelse(runif(nagents_start,0,1)<=0.35010,1,2))
#race distribution: 1=White, 2=Non-white
nw = set.vertex.attribute(nw, attrname="race", value=ifelse(runif(nagents_start,0,1)<=0.35912, 1, 2))

#define network parameters and estimate the network
est = suppressWarnings(netest(nw, formation=~edges+nodefactor("race")+nodematch("race"), target.stats=c(round(nagents_start*0.76/2), round(nagents_start*0.4/2), round(nagents_start*0.6/2/2)), coef.diss=dissolution_coefs(dissolution = ~offset(edges), duration=100)))
ns = netsim(est, param=param.net(inf.prob=0, act.rate=0, a.rate=(0.0466/365), ds.rate=(0.0466/365), di.rate=(0.0466/365)), init=init.net(i.num=0), control=control.net(type="SI", nsteps=(length_sim-2), nsims=nsims, ncores=6, verbose=T))

#diagnostics
#dx = netdx(est, nsims=nsims, nsteps=(length_sim-2), ncores=20, keep.tedgelist=T, nwstats.formula=~edges+nodefactor("race")+nodematch("race"))
#summary(dx)
#plot(dx)

#cleanup
rm(nw,est)
gc()

#benchmarking
#time1 = Sys.time()


### SAVE INITIALIZATION: NETWORK ###

save.image("simulation initialization network_20y_25sim.RData")
load("simulation initialization network_20y_25sim.RData")

#write individual networkDynamic objects for optimization
for (sims in 1:nsims) {
  simnw = ns$network[[sims]]
  save(simnw, file=paste("simulation initialization network_25sims/sim", sims, ".RData", sep=""))
}

rm(simnw,sims,ns)
save.image("simulation initialization network_25sims/global.RData")


### STEP1-3: POPULATION INITIALIZATION ###

library("EpiModel") 
load("simulation initialization network_25sims/global.RData")

for (sims in 1:nsims)
{
  cat("\n\n************** ","Simulation: ",sims," **************\n",sep="")
  
  
  load(file=paste("simulation initialization network_25sims/sim", sims, ".RData", sep=""))
  
  #extract edge activity (i.e. sexual network)
  nw = get.edge.activity(simnw, as.spellList=T)
  nw$onset[nw$onset==-Inf] = 0
  nw$terminus[nw$terminus==Inf] = (length_sim-1)
  
  #extract vertex activity (i.e. migration)
  migration = get.vertex.activity(simnw, as.spellList=T)
  migration$onset[1:nagents_start] = 0 #all original agents start at day 0 to align to edgelist
  nagents_end = c(nagents_end, nrow(migration))
  
  #build sexual network dataset
  sex_network = data.frame(matrix(data=NA, nrow=nagents_end[sims+1], ncol=(6+length_sim)), stringsAsFactors=F)
  colnames(sex_network) = c("ID","ENTRY","EXIT","N_PARTNERS","RELATIONSHIP_DAYS","N_SEX_ACTS", paste("DAY",1:length_sim,sep=""))
  sex_network$ID = migration$vertex.id
  
  #add migration data
  sex_network$ENTRY = migration$onset + 1
  sex_network$EXIT = migration$terminus + 1
  
  for (i in 1:nrow(sex_network))
  {
    #check for partnerships
    partners = nw[nw$tail==i | nw$head==i, ]
    
    if (nrow(partners)>0) {
      #add each partner to daily network
      #note: concurrent relationships will be dissolved so that most recent partner overrides concurrency
      for (j in 1:nrow(partners))
      {
        ids = c(partners$tail[j], partners$head[j])
        sex_network[i, (7+partners$onset[j]):(7+partners$terminus[j])] = ids[!i==ids]
      }
    }
    
    #count partners, relationship days, sex acts (computed based on number of months in a relationship)
    sex_network$N_PARTNERS[i] = length(unique(na.omit(as.numeric(sex_network[i,7:ncol(sex_network)]))))
    sex_network$RELATIONSHIP_DAYS[i] = sum(!is.na(sex_network[i,7:ncol(sex_network)]))
    sex_network$N_SEX_ACTS[i] = round((rpois(1,80.6)/12 * sex_network$RELATIONSHIP_DAYS[i]/30))
  }
  rm(i,j,ids,partners,nw,migration)
  
  
  ### STEP2: INITIALIZING BASE POPULATION ###
  
  #counterfactual within each simulation: ensures that characteristics of the population are the same within each simulation
  set.seed(seed_val_population[sims])
  
  #generate an initial population of MSM, with same characteristics for all prevention scenarios
  baseline_pop = data.frame("ID"=sex_network$ID, "ACTIVE"=NA, "AGE"=NA, "RACE"=NA, "HIV"=NA, "HIV_VL"=NA, "HIV_STATUS"=NA, "SY_ANAL"=NA, "SY_ORAL"=NA, "SY_ANAL_TYPE"=NA, "SY_ORAL_TYPE"=NA, "CIRC"=NA, "SEX_ACT"=NA, "SEX_POSITION_ANAL_PREF"=NA, "SEX_POSITION_ORAL_PREF"=NA, "SEX_POSITION_ANAL"=NA, "SEX_POSITION_ORAL"=NA, "HIV_TEST_DAY"=NA, "SY_TEST_DAY"=NA, "MAX_SEX"=NA, "SEX_COUNT"=0, "PARTNERS"=NA, "ORIGINAL_STATUS_HIV"=NA, "ORIGINAL_STATUS_SY_ANAL"=NA, "ORIGINAL_STATUS_SY_ORAL"=NA, "PREP_PREVENT"=0, "TAP_PREVENT"=0, "SERO_PREVENT"=0, "COND_PREVENT_ANAL"=0, "COND_PREVENT_ORAL"=0, "DOXYPEP_PREVENT_ANAL"=0, "DOXYPEP_PREVENT_ORAL"=0, "DOXYPREP_PREVENT_ANAL"=0, "DOXYPREP_PREVENT_ORAL"=0, "OVERALL_PREVENT_ANAL"=0,"OVERALL_PREVENT_ORAL"=0, "DISCORDANT"=0, "CAUSE_INFECT"=0, "DAYS_KNOWN_HIV_POSITIVE"=-1, "INCIDENCE_DAYS_HIV"=-1, "INCIDENCE_DAYS_SY_ANAL"=-1, "INCIDENCE_DAYS_SY_ORAL"=-1, stringsAsFactors=F)
  
  #age distribution: 1=<=24, 2=>24
  baseline_pop$AGE = ifelse(runif(nagents_end[sims+1],0,1)<=0.3501, 1, 2)
  
  #race distribution: 1=White, 2=Non-white
  #note: inherit from epimodel based on assortative mixing characteristic
  baseline_pop$RACE = get.vertex.attribute(simnw, "race")
  
  #infect initial population with HIV
  baseline_pop$HIV = ifelse(baseline_pop$AGE==1 & baseline_pop$RACE==1 & runif(nagents_end[sims+1],0,1)<=0.023, 1,ifelse(baseline_pop$AGE==2 & baseline_pop$RACE==1 & runif(nagents_end[sims+1],0,1)<=0.136,1, ifelse(baseline_pop$AGE==1 & baseline_pop$RACE==2 & runif(nagents_end[sims+1],0,1)<=0.193, 1, ifelse(baseline_pop$AGE==2 & baseline_pop$RACE==2 & runif(nagents_end[sims+1],0,1)<=0.324, 1,0))))
  
  #HIV viral load category: 1=chronic, 2=acute
  baseline_pop$HIV_VL = ifelse(baseline_pop$HIV==1 & runif(nagents_end[sims+1],0,1)<=0.025, 2, ifelse(baseline_pop$HIV==1, 1, -1))
  
  #knowledge of HIV status: 1=HIV negative and aware (63%), 2=unknown, 3=HIV positive and aware (66%)
  baseline_pop$HIV_STATUS = ifelse(baseline_pop$HIV==0 & runif(nagents_end[sims+1],0,1)<=0.369, 2, ifelse(baseline_pop$HIV==0, 1, NA))
  baseline_pop$HIV_STATUS = ifelse(is.na(baseline_pop$HIV_STATUS) & baseline_pop$HIV==1 & runif(nagents_end[sims+1],0,1)<=0.440, 2, ifelse(baseline_pop$HIV==1, 3, baseline_pop$HIV_STATUS))
  
  #syphilis prevalence: 0=not infectious, 1=primary/secondary, 2=early latent
  sy_type = ifelse(runif(nagents_end[sims+1],0,1)<=0.56422,1,2)
  baseline_pop$SY_ANAL = ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==1 & baseline_pop$RACE==1 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.0409, 1,ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==2 & baseline_pop$RACE==1 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.0758, 1, ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==1 & baseline_pop$RACE==2 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.0729, 1, ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==2 & baseline_pop$RACE==2 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.1354, 1, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==1 & baseline_pop$RACE==1 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.0265,1, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==2 & baseline_pop$RACE==1 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.0492, 1, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==1 & baseline_pop$RACE==2 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.0473,1, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==2 & baseline_pop$RACE==2 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.0877,1, ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==1 & baseline_pop$RACE==1 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.0205, 2,ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==2 & baseline_pop$RACE==1 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.0380, 2, ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==1 & baseline_pop$RACE==2 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.0365, 2, ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==2 & baseline_pop$RACE==2 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.0678, 2, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==1 & baseline_pop$RACE==1 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.0316,2, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==2 & baseline_pop$RACE==1 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.0586, 2, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==1 & baseline_pop$RACE==2 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.0563,2, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==2 & baseline_pop$RACE==2 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.1045,2,0))))))))))))))))
  baseline_pop$SY_ORAL = ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==1 & baseline_pop$RACE==1 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.0348, 1,ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==2 & baseline_pop$RACE==1 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.0646, 1, ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==1 & baseline_pop$RACE==2 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.0621, 1, ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==2 & baseline_pop$RACE==2 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.1153, 1, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==1 & baseline_pop$RACE==1 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.0226,1, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==2 & baseline_pop$RACE==1 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.0419, 1, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==1 & baseline_pop$RACE==2 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.0403,1, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==2 & baseline_pop$RACE==2 & sy_type==1 & runif(nagents_end[sims+1],0,1)<=0.0747,1, ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==1 & baseline_pop$RACE==1 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.0174, 2,ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==2 & baseline_pop$RACE==1 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.0323, 2, ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==1 & baseline_pop$RACE==2 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.0311, 2, ifelse(baseline_pop$HIV==1 & baseline_pop$AGE==2 & baseline_pop$RACE==2 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.0577, 2, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==1 & baseline_pop$RACE==1 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.0269,2, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==2 & baseline_pop$RACE==1 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.0499, 2, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==1 & baseline_pop$RACE==2 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.0480,2, ifelse(baseline_pop$HIV==0 & baseline_pop$AGE==2 & baseline_pop$RACE==2 & sy_type==2 & runif(nagents_end[sims+1],0,1)<=0.0891,2,0))))))))))))))))
  rm(sy_type)
  
  #syphilis transmission dynamics category: 1=primary/secondary (56%), 2=early latent (44%)
  baseline_pop$SY_ANAL_TYPE = ifelse(baseline_pop$SY_ANAL==1, 1, ifelse(baseline_pop$SY_ANAL==2, 2, -1))
  baseline_pop$SY_ORAL_TYPE = ifelse(baseline_pop$SY_ORAL==1, 1, ifelse(baseline_pop$SY_ORAL==2, 2, -1))
  
  #set proportion of prevalent cases as incident
  baseline_pop$INCIDENCE_DAYS_HIV = ifelse(baseline_pop$HIV_VL==2, round(runif(sum(baseline_pop$HIV_VL==2),0,63)), baseline_pop$INCIDENCE_DAYS_HIV)
  baseline_pop$INCIDENCE_DAYS_SY_ANAL = ifelse(baseline_pop$SY_ANAL==1, round(runif(sum(baseline_pop$SY_ANAL),0,365)), baseline_pop$INCIDENCE_DAYS_SY_ANAL)
  baseline_pop$INCIDENCE_DAYS_SY_ORAL = ifelse(baseline_pop$SY_ORAL==1, round(runif(sum(baseline_pop$SY_ORAL),0,365)), baseline_pop$INCIDENCE_DAYS_SY_ORAL)
  
  #circumcision status
  #baseline_pop$CIRC = ifelse(runif(nagents_end[sims+1],0,1)<=0.762, 1, 0)
  baseline_pop$CIRC = 0 #set everyone as uncircumcised as inconclusive evidence to benefit in US population
  
  #sex act preference, 1=anal only, 2=oral only, 3=both
  randact = runif(nagents_end[sims+1],0,1)
  baseline_pop$SEX_ACT = ifelse(randact<=0.22, 1, ifelse(randact<=0.66, 2, 3))
  rm(randact)
  
  #sex_position preference, 1=insertive, 2=receptive, 3=both
  analpref = runif(nagents_end[sims+1],0,1)
  oralpref = runif(nagents_end[sims+1],0,1)
  baseline_pop$SEX_POSITION_ANAL_PREF = ifelse(analpref<=0.24, 1, ifelse(analpref<=0.56, 2, 3))
  baseline_pop$SEX_POSITION_ORAL_PREF = ifelse(oralpref<=0.16, 1, ifelse(oralpref<=0.30, 2, 3))
  rm(analpref,oralpref)
  
  #sex position (scale of dominance), if preference cannot be honored
  baseline_pop$SEX_POSITION_ANAL = runif(nagents_end[sims+1],0,1)
  baseline_pop$SEX_POSITION_ORAL = runif(nagents_end[sims+1],0,1)
  
  #day of year will be tested for HIV, based on a future test probability of 1=tested 1/3 of year, 2=tested 2/3 of year, 3=tested 3/3 of year, 4=not tested
  testprob = runif(nagents_end[sims+1],0,1)
  futuretest = ifelse(baseline_pop$HIV_STATUS==3, 4, ifelse(testprob<=0.19, 1, ifelse(testprob<=0.30, 2, ifelse(testprob<=0.70, 3, 4))))
  probdaytest = runif(nagents_end[sims+1],0,1)
  baseline_pop$HIV_TEST_DAY = ifelse(futuretest==1, round(1+(110-1)*probdaytest), ifelse(futuretest==2, round(111+(220-111)*probdaytest), ifelse(futuretest==3, round(221+(330-221)*probdaytest), NA)))
  rm(testprob,probdaytest,futuretest)
  
  #day of year will be tested for syphilis among HIV+ MSM, based on a future test probability of 1=tested 1/4 of year, 2=tested 1/2 of year, 3=tested in the year, 4=not tested
  testprob = runif(nagents_end[sims+1],0,1)
  futuretest = ifelse(baseline_pop$HIV_STATUS==3 & testprob<=0.22, 1, ifelse(baseline_pop$HIV_STATUS==3 & testprob<=0.43, 2, ifelse(baseline_pop$HIV_STATUS==3 & testprob<=0.71, 3, 4)))
  neg = ifelse(baseline_pop$HIV_STATUS==1 | baseline_pop$HIV_STATUS==2,1,0)
  neg1 = ifelse(neg==1 & futuretest==4,1,0)
  futuretest2 = ifelse(neg1==1 & testprob<=0.31, 1, 2)
  probdaytest = runif(nagents_end[sims+1],0,1)
  baseline_pop$SY_TEST_DAY = ifelse(futuretest==1, round(1+(91-1)*probdaytest), ifelse(futuretest==2, round(92+(182-92)*probdaytest), ifelse(futuretest==3, round(183+(365-183)*probdaytest), NA)))
  baseline_pop$SY_TEST_DAY = ifelse(is.na(baseline_pop$SY_TEST_DAY) & futuretest2==1, round(1+(365-1)*probdaytest), NA)
  rm(testprob,probdaytest,futuretest,futuretest2,neg,neg1)
  
  #original status of infections
  baseline_pop$ORIGINAL_STATUS_HIV = baseline_pop$HIV
  baseline_pop$ORIGINAL_STATUS_SY_ANAL = baseline_pop$SY_ANAL
  baseline_pop$ORIGINAL_STATUS_SY_ORAL = baseline_pop$SY_ORAL
  
  #sex and partner distribution
  baseline_pop$MAX_SEX = sex_network$N_SEX_ACTS
  baseline_pop$PARTNERS = sex_network$N_PARTNERS
  
  
  ### STEP3: ASSIGNMENT OF PREVENTION PARADIGMS ###
  
  #baseline population
  randprep = runif(nagents_end[sims+1],0,1)
  randtap = runif(nagents_end[sims+1],0,1)
  randcond = runif(nagents_end[sims+1],0,1)
  randsero = runif(nagents_end[sims+1],0,1)
  randvax =  runif(nagents_end[sims+1],0,1)
  randdoxy = runif(nagents_end[sims+1],0,1)
  
  #create prevention paradigm datasets: abcdefg where a=prep, b=treatment as prevention, c=condom use, d=seroadaption, e=vaccination, f=doxyprep, g=doxypep, h=syphilis detection/treatment
  paradigm00000000 = baseline_pop
  paradigm11111101 = baseline_pop
  paradigm11111011 = baseline_pop
  
  #PrEP, TasP, condoms, seroadaption, vaccination, doxy PrEP
  paradigm00000000$PREP = 0
  paradigm00000000$TAP = 0
  paradigm00000000$COND = 0
  paradigm00000000$SERO = 0
  paradigm00000000$VAX = 1
  paradigm00000000$DOXYPREP = 0
  paradigm00000000$DOXYPEP = 0
  paradigm00000000$SY_TREAT = 0
  paradigm00000000$TAPVL = prevention_paradigm(randtap,"TAPVL",paradigm00000000)
  paradigm00000000$PROBCOND = prevention_paradigm(randcond,"PROBCOND",paradigm00000000)
  paradigm00000000$SEROTYPE = prevention_paradigm(randsero,"SEROTYPE",paradigm00000000)
  
  #PrEP, TasP, condoms, seroadaption, vaccination, doxy PrEP
  paradigm11111101$PREP = 1 #turn off for calibration to historic data
  paradigm11111101$TAP = 1
  paradigm11111101$COND = 1
  paradigm11111101$SERO = 1
  paradigm11111101$VAX = 1
  paradigm11111101$DOXYPREP = 1
  paradigm11111101$DOXYPEP = 0
  paradigm11111101$SY_TREAT = 1
  paradigm11111101$ONPREP = prevention_paradigm(randprep,"ONPREP",paradigm11111101,0.12)
  paradigm11111101$TAPVL = prevention_paradigm(randtap,"TAPVL",paradigm11111101)
  paradigm11111101$PROBCOND = prevention_paradigm(randcond,"PROBCOND",paradigm11111101)
  paradigm11111101$SEROTYPE = prevention_paradigm(randsero,"SEROTYPE",paradigm11111101)
  
  #PrEP, TasP, condoms, seroadaption, vaccination, doxy PEP
  paradigm11111011$PREP = 1 #turn off for calibration to historic data
  paradigm11111011$TAP = 1
  paradigm11111011$COND = 1
  paradigm11111011$SERO = 1
  paradigm11111011$VAX = 1
  paradigm11111011$DOXYPREP = 0
  paradigm11111011$DOXYPEP = 1
  paradigm11111011$SY_TREAT = 1 
  paradigm11111011$ONPREP = prevention_paradigm(randprep,"ONPREP",paradigm11111011,0.12)
  paradigm11111011$TAPVL = prevention_paradigm(randtap,"TAPVL",paradigm11111011)
  paradigm11111011$PROBCOND = prevention_paradigm(randcond,"PROBCOND",paradigm11111011)
  paradigm11111011$SEROTYPE = prevention_paradigm(randsero,"SEROTYPE",paradigm11111011)
  
  #assign HIV PrEP and doxy into lists
  abm_sims = list("paradigm00000000doxy0"=cbind(paradigm00000000,"ONDOXYPREP"=prevention_paradigm(randdoxy,"ONDOXYPREP",paradigm00000000,0)),
                  "paradigm11111101doxy0"=cbind(paradigm11111101,"ONDOXYPREP"=prevention_paradigm(randdoxy,"ONDOXYPREP",paradigm11111101,0)),
                  "paradigm11111101doxy10"=cbind(paradigm11111101,"ONDOXYPREP"=prevention_paradigm(randdoxy,"ONDOXYPREP",paradigm11111101,0.10)),
                  "paradigm11111101doxy20"=cbind(paradigm11111101,"ONDOXYPREP"=prevention_paradigm(randdoxy,"ONDOXYPREP",paradigm11111101,0.20)),
                  "paradigm11111101doxy30"=cbind(paradigm11111101,"ONDOXYPREP"=prevention_paradigm(randdoxy,"ONDOXYPREP",paradigm11111101,0.30)),
                  "paradigm11111101doxy40"=cbind(paradigm11111101,"ONDOXYPREP"=prevention_paradigm(randdoxy,"ONDOXYPREP",paradigm11111101,0.40)),
                  "paradigm11111101doxy50"=cbind(paradigm11111101,"ONDOXYPREP"=prevention_paradigm(randdoxy,"ONDOXYPREP",paradigm11111101,0.50)),
                  "paradigm11111011doxy0"=cbind(paradigm11111011,"ONDOXYPEP"=prevention_paradigm(randdoxy,"ONDOXYPEP",paradigm11111011,0)),
                  "paradigm11111011doxy10"=cbind(paradigm11111011,"ONDOXYPEP"=prevention_paradigm(randdoxy,"ONDOXYPEP",paradigm11111011,0.10)),
                  "paradigm11111011doxy20"=cbind(paradigm11111011,"ONDOXYPEP"=prevention_paradigm(randdoxy,"ONDOXYPEP",paradigm11111011,0.20)),
                  "paradigm11111011doxy30"=cbind(paradigm11111011,"ONDOXYPEP"=prevention_paradigm(randdoxy,"ONDOXYPEP",paradigm11111011,0.30)),
                  "paradigm11111011doxy40"=cbind(paradigm11111011,"ONDOXYPEP"=prevention_paradigm(randdoxy,"ONDOXYPEP",paradigm11111011,0.40)),
                  "paradigm11111011doxy50"=cbind(paradigm11111011,"ONDOXYPEP"=prevention_paradigm(randdoxy,"ONDOXYPEP",paradigm11111011,0.50)))
  
  
  #initialize a list for longitudinal tracking of outcomes; must match length of abm_sims
  long_sims = list("paradigm00000000doxy0"=data.frame("DAY"=1:length_sim,"HIV"=NA, "SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm11111101doxy0"=data.frame("DAY"=1:length_sim,"HIV"=NA, "SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm11111101doxy10"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm11111101doxy20"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm11111101doxy30"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm11111101doxy40"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm11111101doxy50"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm11111011doxy0"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm11111011doxy10"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm11111011doxy20"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm11111011doxy30"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm11111011doxy40"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm11111011doxy50"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F))
  
  #cleanup
  rm(randtap,randcond,randsero,randprep,randvax,randdoxy,baseline_pop,paradigm11111101,paradigm11111011,paradigm00000000,simnw)
  gc()
  
  
  #add this initialization to the lists
  names(abm_sims) = paste("simulation",sims,names(abm_sims),sep="_")
  abm_sims_50 = c(abm_sims_50,abm_sims)
  names(long_sims) = paste("simulation",sims,names(long_sims),sep="_")
  long_sims_50 = c(long_sims_50,long_sims)
  sex_network = list(sex_network)
  names(sex_network) = paste("simulation",sims,sep="_")
  network_sims_50 = c(network_sims_50,sex_network)
}
rm(abm_sims,long_sims,sex_network,seed_val_population,prevention_paradigm,sims)
abm_sims_50[[1]] = NULL
long_sims_50[[1]] = NULL
network_sims_50[[1]] = NULL
nagents_end = nagents_end[-1]
gc()

#benchmarking
#time2 = Sys.time()

#migration: per time step, % of population expected to enter/exit
#sum(network_sims_50[[1]]$ENTRY==1 & network_sims_50[[1]]$EXIT==length_sim)
#sum(network_sims_50[[1]]$ENTRY!=1 | network_sims_50[[1]]$EXIT!=length_sim)
#sum(network_sims_50[[1]]$ENTRY!=1 | network_sims_50[[1]]$EXIT!=length_sim)/length_sim=
#sum(network_sims_50[[1]]$ENTRY!=1 | network_sims_50[[1]]$EXIT!=length_sim)/length_sim/nagents_start


# Assessment of Baseline Syphilis 
# t = table(abm_sims_50$simulation_1_paradigm1111110doxy0$SY_ANAL); t
# prop.table(t)
# t1 = table(abm_sims_50$simulation_1_paradigm1111110doxy0$SY_ANAL, abm_sims_50$simulation_1_paradigm1111110doxy0$AGE); t1
# prop.table(t1)
# t2 = table(abm_sims_50$simulation_1_paradigm1111110doxy0$SY_ANAL, abm_sims_50$simulation_1_paradigm1111110doxy0$RACE); t2
# prop.table(t2)
# t3 = table(abm_sims_50$simulation_1_paradigm1111110doxy0$SY_ORAL); t3
# prop.table(t3)
# t4 = table(abm_sims_50$simulation_1_paradigm1111110doxy0$SY_ORAL, abm_sims_50$simulation_1_paradigm1111110doxy0$AGE); t4
# prop.table(t4)
# t5 = table(abm_sims_50$simulation_1_paradigm1111110doxy0$SY_ORAL, abm_sims_50$simulation_1_paradigm1111110doxy0$RACE); t5
# prop.table(t5)


### SAVE INITIALIZATION: POPULATION ###

save.image("simulation initialization population_20y_25sim.RData")


### STEP4: SIMULATION ###

load("simulation initialization population_20y_25sim.RData")

l#ensure that all runs start equivalently
set.seed(777)
# nparadigms=22

for (sims in 1:nsims)
{
  #load previously initialized objects
  abm_sims = abm_sims_50[((sims-1)*nparadigms+1):(((sims-1)*nparadigms)+nparadigms)]
  long_sims = long_sims_50[((sims-1)*nparadigms+1):(((sims-1)*nparadigms)+nparadigms)]
  sex_network = network_sims_50[[sims]]
  
  #ensure that each day is unique
  seed_val_day = sample(1:(length_sim*10),length_sim,replace=F)
  
  for (day in 1:length_sim)
  {
    cat("\n\n************** ","Simulation: ",sims," *** Day: ",day," **************\n",sep="")
    
    #active agents (i.e. those in Philadelphia today)
    active = sex_network$ID[sex_network$ENTRY<=day & sex_network$EXIT>=day]
    
    #go through each prevention paradigm scenario
    for (current_data in 1:nparadigms)
    {
      #cat("\n\n************** ","Scenario: ",current_data," **************\n",sep="")
      
      ## SET DAILY VARIABLES ##
      
      #ensure each scenario is identical for a given day
      set.seed(seed_val_day[day])
      
      #set active agents
      abm_sims[[current_data]]$ACTIVE = 0
      abm_sims[[current_data]]$ACTIVE[active] = 1
      
      #check for annual STI testing today
      abm_sims[[current_data]]$HIV_STATUS = ifelse((!is.na(abm_sims[[current_data]]$HIV_TEST_DAY) & ((day-abm_sims[[current_data]]$HIV_TEST_DAY) %% 365)==0) & (abm_sims[[current_data]]$HIV==1), 3, abm_sims[[current_data]]$HIV_STATUS)
      
      #check if testing identify as positive, used for TAP, after 30 days you may be virally suppressed
      abm_sims[[current_data]]$DAYS_KNOWN_HIV_POSITIVE = ifelse((!is.na(abm_sims[[current_data]]$HIV_TEST_DAY)& (day==abm_sims[[current_data]]$HIV_TEST_DAY) & (abm_sims[[current_data]]$HIV==1)), 1, abm_sims[[current_data]]$DAYS_KNOWN_HIV_POSITIVE)
      
      #increment days positive, if positive
      abm_sims[[current_data]]$DAYS_KNOWN_HIV_POSITIVE = ifelse(abm_sims[[current_data]]$DAYS_KNOWN_HIV_POSITIVE>=0, abm_sims[[current_data]]$DAYS_KNOWN_HIV_POSITIVE+1, abm_sims[[current_data]]$DAYS_KNOWN_HIV_POSITIVE)
      
      #in a TAP scenario, if put on treatment, after 30 days viral load goes to undetectable
      abm_sims[[current_data]]$TAPVL = ifelse(abm_sims[[current_data]]$TAP==1 & abm_sims[[current_data]]$DAYS_KNOWN_HIV_POSITIVE==30 & abm_sims[[current_data]]$HIV==1 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.42, 1, abm_sims[[current_data]]$TAPVL)
      
      #increment HIV incident days
      abm_sims[[current_data]]$INCIDENCE_DAYS_HIV = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HIV>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HIV+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HIV)
      
      #resolve acute phase of HIV infection
      abm_sims[[current_data]]$HIV_VL = ifelse(abm_sims[[current_data]]$HIV_VL==2 & abm_sims[[current_data]]$INCIDENCE_DAYS_HIV==63, 1, abm_sims[[current_data]]$HIV_VL)
      
      #if HIV positive and aware and set as receptive
      abm_sims[[current_data]]$SEROTYPE = ifelse(abm_sims[[current_data]]$HIV_STATUS==3 & abm_sims[[current_data]]$SEROTYPE==1, 2, abm_sims[[current_data]]$SEROTYPE) 
      
      #increment SYPhilis INCIDENCE days
      abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ANAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ANAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ANAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ANAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ORAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ORAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ORAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ORAL)
      
      #transition from primary/secondary to early latent syphilis
      abm_sims[[current_data]]$SY_ANAL_TYPE = ifelse(abm_sims[[current_data]]$SY_ANAL_TYPE==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ANAL==366, 2, abm_sims[[current_data]]$SY_ANAL_TYPE)
      abm_sims[[current_data]]$SY_ORAL_TYPE = ifelse(abm_sims[[current_data]]$SY_ORAL_TYPE==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ORAL==366, 2, abm_sims[[current_data]]$SY_ORAL_TYPE)
      
      #treatment for syphilis 
      if (unique(abm_sims[[current_data]]$SY_TREAT)==1) {
        abm_sims[[current_data]]$SY_ANAL = ifelse((!is.na(abm_sims[[current_data]]$SY_TEST_DAY) & ((day-abm_sims[[current_data]]$SY_TEST_DAY) %% 365)==0) & (abm_sims[[current_data]]$SY_ANAL!=0), 0, abm_sims[[current_data]]$SY_ANAL)
        abm_sims[[current_data]]$SY_ORAL = ifelse((!is.na(abm_sims[[current_data]]$SY_TEST_DAY) & ((day-abm_sims[[current_data]]$SY_TEST_DAY) %% 365)==0) & (abm_sims[[current_data]]$SY_ORAL!=0), 0, abm_sims[[current_data]]$SY_ORAL)
      }
      
      ## ENGAGE IN SEX ##
      
      #active partnerships today
      sex_selection = data.frame("Ego"=sex_network$ID[which(!is.na(sex_network[,6+day]))], "Partner"=as.numeric(na.omit(sex_network[,6+day])), stringsAsFactors=F)
      
      #sex today, multiplier is to make the algorithm more conservative from calibration (i.e. less sex)
      sex_selection$Ego_sex = (sex_network$N_SEX_ACTS[sex_selection$Ego] / sex_network$RELATIONSHIP_DAYS[sex_selection$Ego]) <= (runif(nrow(sex_selection),0,1)*1)
      sex_selection$Partner_sex = (sex_network$N_SEX_ACTS[sex_selection$Partner] / sex_network$RELATIONSHIP_DAYS[sex_selection$Partner]) <= (runif(nrow(sex_selection),0,1)*1)
      sex_selection = sex_selection[-1*which(sex_selection$Ego_sex==F | sex_selection$Partner_sex==F), ]
      
      #remove the duplicates (if ego has sex with partner, partner has sex with ego)
      dedup=T; i=1
      while(dedup) {
        dup = which(sex_selection$Ego[i]==sex_selection$Partner)
        if (length(dup)>0) { sex_selection = sex_selection[-dup, ] }
        i=i+1
        if (i>nrow(sex_selection)) { dedup=F }
      }
      rm(dedup,i,dup)
      
      #ensuring at least one partnership
      if (nrow(sex_selection)>0)
      {
        #create the exposure matrix for each pathogen
        sex_exposed = data.frame("Ego"=sex_selection$Ego, "Partner"=sex_selection$Partner, "Exposed"=0, "HIV"=0, "SY_ANAL"=0, "SY_ORAL"=0, stringsAsFactors=F)
        sex_exposed$HIV[which((abm_sims[[current_data]]$HIV[sex_selection$Ego] + abm_sims[[current_data]]$HIV[sex_selection$Partner])==1)] = 1
        sex_exposed$SY_ANAL[which((abm_sims[[current_data]]$SY_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$SY_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$SY_ORAL[which((abm_sims[[current_data]]$SY_ORAL[sex_selection$Ego] + abm_sims[[current_data]]$SY_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$SY_ORAL[which((abm_sims[[current_data]]$SY_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$SY_ORAL[sex_selection$Partner])==1)] = 1
        
        #HIV/SY anal coinfections (used for seroadaption scenario)
        sex_exposed$HIV_SY_ANAL = sex_exposed$HIV * sex_exposed$SY_ANAL
        
        #check for unexposed
        sex_exposed$Exposed = rowSums(sex_exposed[,c("HIV","SY_ANAL","SY_ORAL")])
        sex_unexposed = c(sex_exposed$Ego[sex_exposed$Exposed==0], sex_exposed$Partner[sex_exposed$Exposed==0])
        sex_exposed = sex_exposed[sex_exposed$Exposed!=0, ]
        rm(sex_selection)
        
        #assign for anal intercourse (AI) sex position, 1=ego receptive/partner insertive, 2=ego insertive/partner receptive, 3=flip fuck 
        #based on preference first, then scale of dominance if preference cannot be honored
        sex_exposed$simAI = ifelse(abm_sims[[current_data]]$SEX_POSITION_ANAL_PREF[sex_exposed$Ego]==1 & abm_sims[[current_data]]$SEX_POSITION_ANAL_PREF[sex_exposed$Partner]==2, 2, ifelse(abm_sims[[current_data]]$SEX_POSITION_ANAL_PREF[sex_exposed$Ego]==2 & abm_sims[[current_data]]$SEX_POSITION_ANAL_PREF[sex_exposed$Partner]==1, 1, ifelse(abm_sims[[current_data]]$SEX_POSITION_ANAL_PREF[sex_exposed$Ego]==3 & abm_sims[[current_data]]$SEX_POSITION_ANAL_PREF[sex_exposed$Partner]==3, 3, NA)))
        sex_exposed$simAI = ifelse(is.na(sex_exposed$simAI), ifelse(abs(abm_sims[[current_data]]$SEX_POSITION_ANAL[sex_exposed$Ego]-abm_sims[[current_data]]$SEX_POSITION_ANAL[sex_exposed$Partner])<=0.15, 3, ifelse(abm_sims[[current_data]]$SEX_POSITION_ANAL[sex_exposed$Ego]<abm_sims[[current_data]]$SEX_POSITION_ANAL[sex_exposed$Partner], 1, 2)), sex_exposed$simAI)
        
        #assign for oral intercourse (OI) sex position, 1=ego receptive/partner insertive, 2=ego insertive/partner receptive, 3=flip fuck 
        #based on preference first, then scale of dominance if preference cannot be honored
        sex_exposed$simOI = ifelse(abm_sims[[current_data]]$SEX_POSITION_ORAL_PREF[sex_exposed$Ego]==1 & abm_sims[[current_data]]$SEX_POSITION_ORAL_PREF[sex_exposed$Partner]==2, 2, ifelse(abm_sims[[current_data]]$SEX_POSITION_ORAL_PREF[sex_exposed$Ego]==2 & abm_sims[[current_data]]$SEX_POSITION_ORAL_PREF[sex_exposed$Partner]==1, 1, ifelse(abm_sims[[current_data]]$SEX_POSITION_ORAL_PREF[sex_exposed$Ego]==3 & abm_sims[[current_data]]$SEX_POSITION_ORAL_PREF[sex_exposed$Partner]==3, 3, NA)))
        sex_exposed$simOI = ifelse(is.na(sex_exposed$simOI), ifelse(abs(abm_sims[[current_data]]$SEX_POSITION_ORAL[sex_exposed$Ego]-abm_sims[[current_data]]$SEX_POSITION_ORAL[sex_exposed$Partner])<=0.15, 3, ifelse(abm_sims[[current_data]]$SEX_POSITION_ORAL[sex_exposed$Ego]<abm_sims[[current_data]]$SEX_POSITION_ORAL[sex_exposed$Partner], 1, 2)), sex_exposed$simOI)
        
        #determine sex type, 1=anal (both consent to anal only), 2=oral (one does not consent to anal), 3=both (both consent to anal but one, or both, would also like oral)
        sex_exposed$sex = ifelse(abm_sims[[current_data]]$SEX_ACT[sex_exposed$Ego]==1 & abm_sims[[current_data]]$SEX_ACT[sex_exposed$Partner]==1, 1, ifelse(abm_sims[[current_data]]$SEX_ACT[sex_exposed$Ego]==2 & abm_sims[[current_data]]$SEX_ACT[sex_exposed$Partner]==2, 2, ifelse(abm_sims[[current_data]]$SEX_ACT[sex_exposed$Ego]==3 & abm_sims[[current_data]]$SEX_ACT[sex_exposed$Partner]==3, 3, ifelse(abm_sims[[current_data]]$SEX_ACT[sex_exposed$Ego]==3 & abm_sims[[current_data]]$SEX_ACT[sex_exposed$Partner]==2, 3, ifelse(abm_sims[[current_data]]$SEX_ACT[sex_exposed$Ego]==2 & abm_sims[[current_data]]$SEX_ACT[sex_exposed$Partner]==3, 3, 2)))))
        
        if (sum(sex_exposed$Exposed)>0)
        {
          #pathogen discordance in at least one partnership
          
          ## DETERMINE INFECTION STATUS ##
          #create exposure lists of ego/partner IDs
          exposed_hiv_ego = sex_exposed$Ego[sex_exposed$HIV==1]
          exposed_hiv_partner = sex_exposed$Partner[sex_exposed$HIV==1]
          exposed_sy_anal_ego = sex_exposed$Ego[sex_exposed$SY_ANAL==1]
          exposed_sy_anal_partner = sex_exposed$Partner[sex_exposed$SY_ANAL==1]
          exposed_sy_oral_ego = sex_exposed$Ego[sex_exposed$SY_ORAL==1]
          exposed_sy_oral_partner = sex_exposed$Partner[sex_exposed$SY_ORAL==1]
          
          #check for who is exposed, 1=ego, 2=partner, 3=both (oral)
          sex_exposed$exposure_hiv[sex_exposed$HIV==1] = ifelse(abm_sims[[current_data]]$HIV[exposed_hiv_ego]<abm_sims[[current_data]]$HIV[exposed_hiv_partner], 1, 2)
          sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1] = ifelse(abm_sims[[current_data]]$SY_ANAL[exposed_sy_anal_ego]<abm_sims[[current_data]]$SY_ANAL[exposed_sy_anal_partner], 1, 2)
          sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1] = ifelse((abm_sims[[current_data]]$SY_ORAL[exposed_sy_oral_ego]<abm_sims[[current_data]]$SY_ANAL[exposed_sy_oral_partner]) & (abm_sims[[current_data]]$SY_ANAL[exposed_sy_oral_ego]>abm_sims[[current_data]]$SY_ORAL[exposed_sy_oral_partner]), 3, ifelse((abm_sims[[current_data]]$SY_ORAL[exposed_sy_oral_ego]>abm_sims[[current_data]]$SY_ANAL[exposed_sy_oral_partner]) & (abm_sims[[current_data]]$SY_ANAL[exposed_sy_oral_ego]<abm_sims[[current_data]]$SY_ORAL[exposed_sy_oral_partner]), 3, ifelse((abm_sims[[current_data]]$SY_ORAL[exposed_sy_oral_ego]<abm_sims[[current_data]]$SY_ANAL[exposed_sy_oral_partner]) | (abm_sims[[current_data]]$SY_ANAL[exposed_sy_oral_ego]<abm_sims[[current_data]]$SY_ORAL[exposed_sy_oral_partner]), 1, ifelse((abm_sims[[current_data]]$SY_ORAL[exposed_sy_oral_ego]>abm_sims[[current_data]]$SY_ANAL[exposed_sy_oral_partner]) | (abm_sims[[current_data]]$SY_ANAL[exposed_sy_oral_ego]>abm_sims[[current_data]]$SY_ORAL[exposed_sy_oral_partner]), 2, NA))))
          
          #indicators for various metrics
          sex_exposed$condct_hiv = 0            #HIV infection blocked from condom use, ego
          sex_exposed$condct_sy_anal = 0        #SY anal infection blocked from condom use, ego
          sex_exposed$condct_sy_oral = 0        #SY oral infection blocked from condom use, ego
          sex_exposed$prepct = 0                #infection blocked from PrEP, ego
          sex_exposed$tapct = 0                 #infection blocked from TAP, ego
          sex_exposed$seroct_hiv = 0            #HIV infection blocked from seroadapation, ego
          sex_exposed$seroct_sy = 0             #SY anal infection blocked from seroadapation, ego
          sex_exposed$doxypepct_anal = 0        #SY anal infection blocked from doxycycline pep, ego
          sex_exposed$doxyprepct_anal = 0       #SY anal infection blocked from doxycycline prep, ego
          sex_exposed$doxypepct_oral = 0        #SY oral infection blocked from doxycycline pep, ego
          sex_exposed$doxyprepct_oral = 0       #SY oral infection blocked from doxycycline prep, ego
          sex_exposed$prevct_sy_anal = 0        #SY anal infection blocked from any prevention strategy, ego
          sex_exposed$prevct_sy_oral = 0        #SY oral infection blocked from any prevention strategy, ego
          sex_exposed$infectct = 0              #ego caused a new SY anal infection to partner
          sex_exposed$newHIV = 0                #new HIV infection, ego
          sex_exposed$newVL = -1                #set HIV_VL to acute if new infection, ego
          sex_exposed$newSY_ANAL = 0            #any new SY anal infection, ego
          sex_exposed$newSY_ORAL = 0            #any new SY oral infection, ego
          sex_exposed$newSY_ANAL_PS = 0         #new primary/secondary SY anal infection, ego
          sex_exposed$newSY_ORAL_PS = 0         #new primary/secondary SY oral infection, ego
          sex_exposed$newSY_ANAL_EL = 0         #new early latent SY anal infection, ego
          sex_exposed$newSY_ORAL_EL = 0         #new early latent SY oral infection, ego
          sex_exposed$newSY_ANAL_TYPE = -1      #set SY_ANAL_TYPE to primary/secondary if new infection, ego
          sex_exposed$newSY_ORAL_TYPE = -1      #set SY_ORAL_TYPE to primary/secondary if new infection, ego
          sex_exposed$partcondct_hiv = 0        #HIV infection blocked from condom use, partner
          sex_exposed$partcondct_sy_anal = 0    #SY anal infection blocked from condom use, partner
          sex_exposed$partcondct_sy_oral = 0    #SY oral infection blocked from condom use, partner
          sex_exposed$partprepct = 0            #infection blocked from PrEP, partner
          sex_exposed$parttapct = 0             #infection blocked from TAP, partner
          sex_exposed$partseroct_hiv = 0        #HIV infection blocked from seroadapation, partner
          sex_exposed$partseroct_sy = 0         #SY anal infection blocked from seroadapation, partner
          sex_exposed$partdoxypepct_anal = 0    #SY anal infection blocked from doxycycline pep, partner
          sex_exposed$partdoxyprepct_anal = 0   #SY anal infection blocked from doxycycline prep, partner
          sex_exposed$partdoxypepct_oral = 0    #SY oral infection blocked from doxycycline pep, partner
          sex_exposed$partdoxyprepct_oral = 0   #SY oral infection blocked from doxycycline prep, partner
          sex_exposed$partprevct_sy_anal = 0    #SY anal infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_sy_oral = 0    #SY oral infection blocked from any prevention strategy, partner
          sex_exposed$partinfectct = 0          #partner caused a new SY anal infection to partner
          sex_exposed$partnewHIV = 0            #new HIV infection, partner
          sex_exposed$partnewVL = -1            #set HIV_VL to acute if new infection, partner
          sex_exposed$partnewSY_ANAL_TYPE = -1  #set SY_ANAL_TYPE to primary/secondary if new infection, partner
          sex_exposed$partnewSY_ORAL_TYPE = -1  #set SY_ORAL_TYPE to primary/secondary if new infection, partner 
          sex_exposed$partnewSY_ANAL = 0        #new SY anal infection, partner
          sex_exposed$partnewSY_ORAL = 0        #new SY oral infection, partner
          sex_exposed$partnewSY_ANAL_PS = 0     #new primary/secondary SY anal infection, partner
          sex_exposed$partnewSY_ORAL_PS = 0     #new primary/secondary SY oral infection, partner
          sex_exposed$partnewSY_ANAL_EL = 0     #new early latent SY anal infection, partner
          sex_exposed$partnewSY_ORAL_EL = 0     #new early latent SY oral infection, partner
          
          #calculate probabilities of infection
          probinfect1_hiv = runif(nrow(sex_exposed),0,1)
          probinfect2_hiv = runif(nrow(sex_exposed),0,1)
          probinfect1_sy_anal = runif(nrow(sex_exposed),0,1)
          probinfect2_sy_anal = runif(nrow(sex_exposed),0,1)
          probinfect1_sy_oral = runif(nrow(sex_exposed),0,1)
          probinfect2_sy_oral = runif(nrow(sex_exposed),0,1)
          
          #determine infection under no prevention for each pathogen
          
          #hiv anal, ego then partner for chronic then acute
          sex_exposed$infection_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0134, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0010, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==0 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0059, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0134, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0010, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==0 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0059, 1, 0))))))
          sex_exposed$infection_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.1284, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0101, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==0 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0569, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.1284, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0101, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==0 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0569, 1, sex_exposed$infection_hiv[sex_exposed$HIV==1]))))))
          sex_exposed$partinfection_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0134, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0010, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==0 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0059, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0134, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0010, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==0 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0059, 1, 0))))))
          sex_exposed$partinfection_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.1284, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0101, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==0 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0569, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.1284, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0101, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==0 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0569, 1, sex_exposed$partinfection_hiv[sex_exposed$HIV==1]))))))
          
          # Simulation 9 ####
          #syphilis anal, HIV negative then positive, ego then partner for primary/secondary then latent syphilis
          sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==1 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==0 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==1 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==0 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, 0))))))
          sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==1 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==2 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.005, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.005, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==0 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.005, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==2 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.005, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==1 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.005, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==0 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.005, 1, 0))))))
          sex_exposed$partinfection_sy_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==1 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==0 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==1 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==0 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, 0))))))
          sex_exposed$partinfection_sy_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==1 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==2 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.005, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.005, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==0 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.005, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==2 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.005, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==1 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.005, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==0 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.005, 1, 0))))))
          
          #Include mulitplicative increase for HIV coinfection 
          sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==1 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==1 & abm_sims[[current_data]]$HIV[exposed_sy_anal_partner]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.020, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==1 & abm_sims[[current_data]]$HIV[exposed_sy_anal_partner]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.020, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==0 & abm_sims[[current_data]]$HIV[exposed_sy_anal_partner]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.020, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==1 & abm_sims[[current_data]]$HIV[exposed_sy_anal_partner]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.020, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==1 & abm_sims[[current_data]]$HIV[exposed_sy_anal_partner]==1 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.020, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==0 & abm_sims[[current_data]]$HIV[exposed_sy_anal_partner]==1 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.020, 1, sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1]))))))
          sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==1 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==2 & abm_sims[[current_data]]$HIV[exposed_sy_anal_partner]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==1 & abm_sims[[current_data]]$HIV[exposed_sy_anal_partner]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==0 & abm_sims[[current_data]]$HIV[exposed_sy_anal_partner]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==2 & abm_sims[[current_data]]$HIV[exposed_sy_anal_partner]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==1 & abm_sims[[current_data]]$HIV[exposed_sy_anal_partner]==1 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_ego]==0 & abm_sims[[current_data]]$HIV[exposed_sy_anal_partner]==1 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1]))))))
          sex_exposed$partinfection_sy_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==1 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==1 & abm_sims[[current_data]]$HIV[exposed_sy_anal_ego]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.020, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==1 & abm_sims[[current_data]]$HIV[exposed_sy_anal_ego]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.020, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==0 & abm_sims[[current_data]]$HIV[exposed_sy_anal_ego]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.020, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==1 & abm_sims[[current_data]]$HIV[exposed_sy_anal_ego]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.020, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==1 & abm_sims[[current_data]]$HIV[exposed_sy_anal_ego]==1 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.020, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==0 & abm_sims[[current_data]]$HIV[exposed_sy_anal_ego]==1 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.020, 1, sex_exposed$partinfection_sy_anal[sex_exposed$SY_ANAL==1]))))))
          sex_exposed$partinfection_sy_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==1 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==2 & abm_sims[[current_data]]$HIV[exposed_sy_anal_ego]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==1 & abm_sims[[current_data]]$HIV[exposed_sy_anal_ego]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==2 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==0 & abm_sims[[current_data]]$HIV[exposed_sy_anal_ego]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==2 & abm_sims[[current_data]]$HIV[exposed_sy_anal_ego]==1 & probinfect1_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==1 & abm_sims[[current_data]]$HIV[exposed_sy_anal_ego]==1 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ANAL==1]!=2 & sex_exposed$exposure_sy_anal[sex_exposed$SY_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$SY_ANAL==1]==3 & abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_anal_partner]==0 & abm_sims[[current_data]]$HIV[exposed_sy_anal_ego]==1 & probinfect2_sy_anal[sex_exposed$SY_ANAL==1] <= 0.010, 1, sex_exposed$partinfection_sy_anal[sex_exposed$SY_ANAL==1]))))))
          
          #syphilis type specific oral, ego then partner, 1=new oral SYPHILIS infection, 2=new genital SYPHILIS infection, 3=new oral & genital SYPHILIS infection
          sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==2 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==2 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==0 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, ifelse((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004) & ((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004) | (sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==0 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004)), 3, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==0 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, 0)))))))
          sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==2 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==2 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==0 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002, 2, ifelse((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002) & ((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002) | (sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==0 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002)), 3, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==0 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002, 2, 0)))))))
          sex_exposed$partinfection_sy_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==2 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==0 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, ifelse((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004) & ((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004) | (sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==0 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004)), 3, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==0 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, 0)))))))
          sex_exposed$partinfection_sy_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==2 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==0 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002, 2, ifelse((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002) & ((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002) | (sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==0 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002)), 3, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==0 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.002, 2, 0)))))))
          
          #Include multiplicative increase for HIV coninfection
          sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==2 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==2 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==0 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008, 2, ifelse((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008) & ((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==0 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008)), 3, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==0 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008, 2, sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1])))))))
          sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==2 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==2 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==0 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, ifelse((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004) & ((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004) | (sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==0 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004)), 3, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_ego]==0 & abm_sims[[current_data]]$HIV[exposed_sy_oral_partner]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1])))))))
          sex_exposed$partinfection_sy_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==2 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==0 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008, 2, ifelse((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008) & ((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==0 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008)), 3, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==0 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.008, 2, sex_exposed$partinfection_sy_oral[sex_exposed$SY_ORAL==1])))))))
          sex_exposed$partinfection_sy_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==2 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==0 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, ifelse((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004) & ((sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004) | (sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==0 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004)), 3, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect1_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 1, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==1 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, ifelse(sex_exposed$sex[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$exposure_sy_oral[sex_exposed$SY_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$SY_ORAL==1]==3 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner]==2 & abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_sy_oral_partner]==0 & abm_sims[[current_data]]$HIV[exposed_sy_oral_ego]==1 & probinfect2_sy_oral[sex_exposed$SY_ORAL==1] <= 0.004, 2, sex_exposed$partinfection_sy_oral[sex_exposed$SY_ORAL==1])))))))
          
          #determine infection under seroadaption scenarios
          if (unique(abm_sims[[current_data]]$SERO)==1)
          {
            #check if avoiding sex: both insertive seroadapters, or serosorter (ego/partner) and partner is aware HIV+
            avoid = ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==3 & abm_sims[[current_data]]$HIV_STATUS[exposed_hiv_partner]==3, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==3 & abm_sims[[current_data]]$HIV_STATUS[exposed_hiv_partner]==3, 1, 0)))
            partavoid = ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==3 & abm_sims[[current_data]]$HIV_STATUS[exposed_hiv_ego]==3, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==3 & abm_sims[[current_data]]$HIV_STATUS[exposed_hiv_ego]==3, 1, 0)))
            
            #seroadaption strategies for HIV
            seroblock_hiv = ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0010, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0059, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0101, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0569, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0010, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0059, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0101, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0569, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0134, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.1284, 1, 0))))))))))
            partseroblock_hiv = ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0010, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0059, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0101, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0569, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0010, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0059, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0101, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0569, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0134, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.1284, 1, 0))))))))))
            
            #check if infection blocked for HIV
            sex_exposed$seroct_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==1, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & seroblock_hiv==1, 1, 0))
            sex_exposed$partseroct_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==1, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & partseroblock_hiv==1, 1, 0))
            
            #seroadaption strategies for HIV, SYPHILIS implications; NAs (and possibly 0s) get returned for HIV infections without correspoding SY infection
            seroblock_sy = ifelse(sex_exposed$infection_sy_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_sy_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_sy_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_sy_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_sy_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            partseroblock_sy = ifelse(sex_exposed$partinfection_sy_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_sy_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_sy_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_sy_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_sy_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            
            #check if infection blocked for SYPHILIS anal, removing NAs/0s from earlier that are not coinfections
            if (sum(sex_exposed$HIV_SY_ANAL)>0) { sex_exposed$seroct_sy[sex_exposed$HIV_SY_ANAL==1] = (ifelse(sex_exposed$infection_sy_anal[sex_exposed$HIV==1]==1 & avoid==1, 1, ifelse(sex_exposed$infection_sy_anal[sex_exposed$HIV==1]==1 & avoid==0 & seroblock_sy==1, 1, 0)))[which(!is.na(sex_exposed$infection_sy_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_SY_ANAL)>0) { sex_exposed$partseroct_sy[sex_exposed$HIV_SY_ANAL==1] = (ifelse(sex_exposed$partinfection_sy_anal[sex_exposed$HIV==1]==1 & partavoid==1, 1, ifelse(sex_exposed$partinfection_sy_anal[sex_exposed$HIV==1]==1 & partavoid==0 & partseroblock_sy==1, 1, 0)))[which(!is.na(sex_exposed$partinfection_sy_anal[sex_exposed$HIV==1]))] }
            
            #decrement sex counter for those who have avoided sex (because it later gets incremented irrespective of seroadaption outcome)
            abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Ego[sex_exposed$HIV==1]] = ifelse(avoid==1, abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Ego[sex_exposed$HIV==1]] - 1, abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Ego[sex_exposed$HIV==1]])
            abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Partner[sex_exposed$HIV==1]] = ifelse(partavoid==1, abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Partner[sex_exposed$HIV==1]] - 1, abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Partner[sex_exposed$HIV==1]])
            
            rm(avoid,partavoid,seroblock_hiv,partseroblock_hiv,seroblock_sy,partseroblock_sy)
          }
          
          #determine infection under condom scenarios
          if (unique(abm_sims[[current_data]]$COND)==1)
          {
            #check if a condom was worn by either partner
            randcond = runif(nrow(sex_exposed),0,1)
            condom_worn_anal = (randcond<=abm_sims[[current_data]]$PROBCOND[sex_exposed$Ego]) | (randcond<=abm_sims[[current_data]]$PROBCOND[sex_exposed$Partner])
            condom_worn_oral = (randcond<=0.04)
            
            #check if the condom failed
            randfail = runif(nrow(sex_exposed),0,1)
            condom_success_anal = condom_worn_anal & (randfail<=0.705)
            condom_success_oral = condom_worn_oral & (randfail<=0.705)
            
            #check if infection blocked for HIV
            sex_exposed$condct_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & condom_success_anal[sex_exposed$HIV==1], 1, 0)
            sex_exposed$partcondct_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & condom_success_anal[sex_exposed$HIV==1], 1, 0)
            
            #check if infection blocked for SY anal
            sex_exposed$condct_sy_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1]==1 & condom_success_anal[sex_exposed$SY_ANAL==1], 1, 0)
            sex_exposed$partcondct_sy_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$partinfection_sy_anal[sex_exposed$SY_ANAL==1]==1 & condom_success_anal[sex_exposed$SY_ANAL==1], 1, 0)
            
            #check if infection blocked for SY oral
            sex_exposed$condct_sy_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1]>=1 & condom_success_oral[sex_exposed$SY_ORAL==1], 1, 0)
            sex_exposed$partcondct_sy_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$partinfection_sy_oral[sex_exposed$SY_ORAL==1]>=1 & condom_success_oral[sex_exposed$SY_ORAL==1], 1, 0)
            
            rm(randcond,randfail,condom_worn_anal,condom_success_anal,condom_worn_oral,condom_success_oral)
          }
          
          #determine infection under TAP scenarios
          if (unique(abm_sims[[current_data]]$TAP)==1)
          {
            #check if infection blocked
            sex_exposed$tapct[sex_exposed$HIV==1] = ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$TAPVL[exposed_hiv_partner]==1, 1, 0)
            sex_exposed$parttapct[sex_exposed$HIV==1] = ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$TAPVL[exposed_hiv_ego]==1, 1, 0)
          }
          
          #determine infection under PrEP scenarios
          if (unique(abm_sims[[current_data]]$PREP)==1)
          {
            #check if infection blocked
            sex_exposed$prepct[sex_exposed$HIV==1] = ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$ONPREP[exposed_hiv_ego]==1, 1, 0)
            sex_exposed$partprepct[sex_exposed$HIV==1] = ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$ONPREP[exposed_hiv_partner]==1, 1, 0)
          }
          
          #determine infection under doxy PEP scenarios 
          if (unique(abm_sims[[current_data]]$DOXYPEP)==1)
          {
            #calculate probabilities of doxycycline success
            randdoxy = runif(nrow(sex_exposed),0,1)

            #check if infection blocked
            sex_exposed$doxypepct_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1]==1 & abm_sims[[current_data]]$ONDOXYPEP[exposed_sy_anal_ego]==1 & randdoxy[sex_exposed$SY_ANAL==1]<=0.73, 1, 0)
            sex_exposed$partdoxypepct_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1]==1 & abm_sims[[current_data]]$ONDOXYPEP[exposed_sy_anal_partner]==1 & randdoxy[sex_exposed$SY_ANAL==1]<=0.73, 1, 0)
            sex_exposed$doxypepct_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$ONDOXYPEP[exposed_sy_oral_ego]==1 & randdoxy[sex_exposed$SY_ORAL==1]<=0.73, 1, 0)
            sex_exposed$partdoxypepct_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$ONDOXYPEP[exposed_sy_oral_partner]==1 & randdoxy[sex_exposed$SY_ORAL==1]<=0.73, 1, 0)
            
            rm(randdoxy)
          }
          
          #determine infection under doxy PREP scenarios 
          if (unique(abm_sims[[current_data]]$DOXYPREP)==1)
          {
            #calculate probabilities of doxycycline success
            randdoxy = runif(nrow(sex_exposed),0,1)
            
            #check if infection blocked
            sex_exposed$doxyprepct_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1]==1 & abm_sims[[current_data]]$ONDOXYPREP[exposed_sy_anal_ego]==1 & randdoxy[sex_exposed$SY_ANAL==1]<=0.73, 1, 0)
            sex_exposed$partdoxyprepct_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1]==1 & abm_sims[[current_data]]$ONDOXYPREP[exposed_sy_anal_partner]==1 & randdoxy[sex_exposed$SY_ANAL==1]<=0.73, 1, 0)
            sex_exposed$doxyprepct_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$ONDOXYPREP[exposed_sy_oral_ego]==1 & randdoxy[sex_exposed$SY_ORAL==1]<=0.73, 1, 0)
            sex_exposed$partdoxyprepct_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$ONDOXYPREP[exposed_sy_oral_partner]==1 & randdoxy[sex_exposed$SY_ORAL==1]<=0.73, 1, 0)
            
            rm(randdoxy)
          }
          
          
          #resolve infection statistics, HIV
          sex_exposed$newHIV[sex_exposed$HIV==1] = ifelse((sex_exposed$condct_hiv[sex_exposed$HIV==1] + sex_exposed$prepct[sex_exposed$HIV==1] + sex_exposed$tapct[sex_exposed$HIV==1] + sex_exposed$seroct_hiv[sex_exposed$HIV==1])==0 & sex_exposed$infection_hiv[sex_exposed$HIV==1]==1, 1, sex_exposed$newHIV[sex_exposed$HIV==1])
          sex_exposed$newVL[sex_exposed$HIV==1] = ifelse((sex_exposed$condct_hiv[sex_exposed$HIV==1] + sex_exposed$prepct[sex_exposed$HIV==1] + sex_exposed$tapct[sex_exposed$HIV==1] + sex_exposed$seroct_hiv[sex_exposed$HIV==1])==0 & sex_exposed$infection_hiv[sex_exposed$HIV==1]==1, 2, sex_exposed$newVL[sex_exposed$HIV==1])
          sex_exposed$partnewHIV[sex_exposed$HIV==1] = ifelse((sex_exposed$partcondct_hiv[sex_exposed$HIV==1] + sex_exposed$partprepct[sex_exposed$HIV==1] + sex_exposed$parttapct[sex_exposed$HIV==1] + sex_exposed$partseroct_hiv[sex_exposed$HIV==1])==0 & sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1, 1, sex_exposed$partnewHIV[sex_exposed$HIV==1])
          sex_exposed$partnewVL[sex_exposed$HIV==1] = ifelse((sex_exposed$partcondct_hiv[sex_exposed$HIV==1] + sex_exposed$partprepct[sex_exposed$HIV==1] + sex_exposed$parttapct[sex_exposed$HIV==1] + sex_exposed$partseroct_hiv[sex_exposed$HIV==1])==0 & sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1, 2, sex_exposed$partnewVL[sex_exposed$HIV==1])
          
          #resolve infection statistics, SYPHILIS anal
          sex_exposed$newSY_ANAL[sex_exposed$SY_ANAL==1] = ifelse((sex_exposed$condct_sy_anal[sex_exposed$SY_ANAL==1] + sex_exposed$seroct_sy[sex_exposed$SY_ANAL==1] + sex_exposed$doxypepct_anal[sex_exposed$SY_ANAL==1] + sex_exposed$doxyprepct_anal[sex_exposed$SY_ANAL==1])==0 & sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1]==1, 1, sex_exposed$newSY_ANAL[sex_exposed$SY_ANAL==1])
          sex_exposed$newSY_ANAL_TYPE[sex_exposed$SY_ANAL==1] = ifelse((sex_exposed$condct_sy_anal[sex_exposed$SY_ANAL==1] + sex_exposed$seroct_sy[sex_exposed$SY_ANAL==1] + sex_exposed$doxypepct_anal[sex_exposed$SY_ANAL==1] + sex_exposed$doxyprepct_anal[sex_exposed$SY_ANAL==1])==0 & sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1]==1, 1, sex_exposed$newSY_ANAL_TYPE[sex_exposed$SY_ANAL==1])
          sex_exposed$partnewSY_ANAL[sex_exposed$SY_ANAL==1] = ifelse((sex_exposed$partcondct_sy_anal[sex_exposed$SY_ANAL==1] + sex_exposed$partseroct_sy[sex_exposed$SY_ANAL==1] + sex_exposed$partdoxypepct_anal[sex_exposed$SY_ANAL==1] + sex_exposed$partdoxyprepct_anal[sex_exposed$SY_ANAL==1])==0 & sex_exposed$partinfection_sy_anal[sex_exposed$SY_ANAL==1]==1, 1, sex_exposed$partnewSY_ANAL[sex_exposed$SY_ANAL==1])
          sex_exposed$partnewSY_ANAL_TYPE[sex_exposed$SY_ANAL==1] = ifelse((sex_exposed$partcondct_sy_anal[sex_exposed$SY_ANAL==1] + sex_exposed$partseroct_sy[sex_exposed$SY_ANAL==1] + sex_exposed$partdoxypepct_anal[sex_exposed$SY_ANAL==1] + sex_exposed$partdoxyprepct_anal[sex_exposed$SY_ANAL==1])==0 & sex_exposed$partinfection_sy_anal[sex_exposed$SY_ANAL==1]==1, 1, sex_exposed$partnewSY_ANAL_TYPE[sex_exposed$SY_ANAL==1])
          sex_exposed$prevct_sy_anal[sex_exposed$SY_ANAL==1] = ifelse((sex_exposed$condct_sy_anal[sex_exposed$SY_ANAL==1] + sex_exposed$seroct_sy[sex_exposed$SY_ANAL==1] + sex_exposed$doxypepct_anal[sex_exposed$SY_ANAL==1] + sex_exposed$doxyprepct_anal[sex_exposed$SY_ANAL==1])>0 & sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1]==1, 1, sex_exposed$prevct_sy_anal[sex_exposed$SY_ANAL==1])
          sex_exposed$partprevct_sy_anal[sex_exposed$SY_ANAL==1] = ifelse((sex_exposed$partcondct_sy_anal[sex_exposed$SY_ANAL==1] + sex_exposed$partseroct_sy[sex_exposed$SY_ANAL==1] + sex_exposed$partdoxypepct_anal[sex_exposed$SY_ANAL==1] + sex_exposed$partdoxyprepct_anal[sex_exposed$SY_ANAL==1])>0 & sex_exposed$partinfection_sy_anal[sex_exposed$SY_ANAL==1]==1, 1, sex_exposed$partprevct_sy_anal[sex_exposed$SY_ANAL==1])
          
          #resolve infection statistics, SYPHILIS oral
          sex_exposed$newSY_ORAL[sex_exposed$SY_ORAL==1] = ifelse((sex_exposed$condct_sy_oral[sex_exposed$SY_ORAL==1] + sex_exposed$doxypepct_oral[sex_exposed$SY_ORAL==1] + sex_exposed$doxyprepct_oral[sex_exposed$SY_ORAL==1])==0 & sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1]>=1, sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1], sex_exposed$newSY_ORAL[sex_exposed$SY_ORAL==1])
          sex_exposed$newSY_ORAL_TYPE[sex_exposed$SY_ORAL==1] = ifelse((sex_exposed$condct_sy_oral[sex_exposed$SY_ORAL==1] + sex_exposed$doxypepct_oral[sex_exposed$SY_ORAL==1] + sex_exposed$doxyprepct_oral[sex_exposed$SY_ORAL==1])==0 & sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1]==1, 1, sex_exposed$newSY_ORAL_TYPE[sex_exposed$SY_ORAL==1])
          sex_exposed$partnewSY_ORAL[sex_exposed$SY_ORAL==1] = ifelse((sex_exposed$partcondct_sy_oral[sex_exposed$SY_ORAL==1] + sex_exposed$partdoxypepct_oral[sex_exposed$SY_ORAL==1] + sex_exposed$partdoxyprepct_oral[sex_exposed$SY_ORAL==1])==0 & sex_exposed$partinfection_sy_oral[sex_exposed$SY_ORAL==1]>=1, sex_exposed$partinfection_sy_oral[sex_exposed$SY_ORAL==1], sex_exposed$partnewSY_ORAL[sex_exposed$SY_ORAL==1])
          sex_exposed$partnewSY_ORAL_TYPE[sex_exposed$SY_ORAL==1] = ifelse((sex_exposed$partcondct_sy_oral[sex_exposed$SY_ORAL==1] + sex_exposed$partdoxypepct_oral[sex_exposed$SY_ORAL==1] + sex_exposed$partdoxyprepct_oral[sex_exposed$SY_ORAL==1])==0 & sex_exposed$partinfection_sy_oral[sex_exposed$SY_ORAL==1]==1, 1, sex_exposed$partnewSY_ORAL_TYPE[sex_exposed$SY_ORAL==1])
          sex_exposed$prevct_sy_oral[sex_exposed$SY_ORAL==1] = ifelse((sex_exposed$condct_sy_oral[sex_exposed$SY_ORAL==1] + sex_exposed$doxypepct_oral[sex_exposed$SY_ORAL==1] + sex_exposed$doxyprepct_oral[sex_exposed$SY_ORAL==1])>0 & sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1]>=1, 1, sex_exposed$prevct_sy_oral[sex_exposed$SY_ORAL==1])
          sex_exposed$partprevct_sy_oral[sex_exposed$SY_ORAL==1] = ifelse((sex_exposed$partcondct_sy_oral[sex_exposed$SY_ORAL==1] + sex_exposed$partdoxypepct_oral[sex_exposed$SY_ORAL==1] + sex_exposed$partdoxyprepct_oral[sex_exposed$SY_ORAL==1])>0 & sex_exposed$partinfection_sy_oral[sex_exposed$SY_ORAL==1]>=1, 1, sex_exposed$partprevct_sy_oral[sex_exposed$SY_ORAL==1])
          
          
          ## UPDATE STATS IN MAIN DATASET ##
          
          #HIV
          #abm_sims[[current_data]]$PREP_PREVENT[exposed_hiv_ego] = abm_sims[[current_data]]$PREP_PREVENT[exposed_hiv_ego] + sex_exposed$prepct[sex_exposed$HIV==1]
          #abm_sims[[current_data]]$PREP_PREVENT[exposed_hiv_partner] = abm_sims[[current_data]]$PREP_PREVENT[exposed_hiv_partner] + sex_exposed$partprepct[sex_exposed$HIV==1]
          #abm_sims[[current_data]]$TAP_PREVENT[exposed_sy_anal_ego] = abm_sims[[current_data]]$TAP_PREVENT[exposed_sy_anal_ego] + sex_exposed$tapct[sex_exposed$HIV==1]
          #abm_sims[[current_data]]$TAP_PREVENT[exposed_sy_anal_partner] = abm_sims[[current_data]]$TAP_PREVENT[exposed_sy_anal_partner] + sex_exposed$parttapct[sex_exposed$HIV==1]
          #abm_sims[[current_data]]$CAUSE_INFECT[exposed_hiv_ego] = ifelse(sex_exposed$infectct[sex_exposed$HIV==1]==1, abm_sims[[current_data]]$CAUSE_INFECT[exposed_hiv_ego] + sex_exposed$infectct[sex_exposed$HIV==1], abm_sims[[current_data]]$CAUSE_INFECT[exposed_hiv_ego])
          #abm_sims[[current_data]]$CAUSE_INFECT[exposed_hiv_partner] = ifelse(sex_exposed$partinfectct[sex_exposed$HIV==1]==1, abm_sims[[current_data]]$CAUSE_INFECT[exposed_hiv_partner] + sex_exposed$partinfectct[sex_exposed$HIV==1], abm_sims[[current_data]]$CAUSE_INFECT[exposed_hiv_partner])
          abm_sims[[current_data]]$HIV[exposed_hiv_ego] = ifelse(sex_exposed$newHIV[sex_exposed$HIV==1]==1, 1, abm_sims[[current_data]]$HIV[exposed_hiv_ego])
          abm_sims[[current_data]]$HIV[exposed_hiv_partner] = ifelse(sex_exposed$partnewHIV[sex_exposed$HIV==1]==1, 1, abm_sims[[current_data]]$HIV[exposed_hiv_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HIV[exposed_hiv_ego] = ifelse(sex_exposed$newHIV[sex_exposed$HIV==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HIV[exposed_hiv_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HIV[exposed_hiv_partner] = ifelse(sex_exposed$partnewHIV[sex_exposed$HIV==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HIV[exposed_hiv_partner])
          abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego] = ifelse(sex_exposed$newVL[sex_exposed$HIV==1]==2, 2, abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego])
          abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner] = ifelse(sex_exposed$partnewVL[sex_exposed$HIV==1]==2, 2, abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner])
          #abm_sims[[current_data]]$DISCORDANT[exposed_hiv_ego] = abm_sims[[current_data]]$DISCORDANT[exposed_hiv_ego] + 1
          #abm_sims[[current_data]]$DISCORDANT[exposed_hiv_partner] = abm_sims[[current_data]]$DISCORDANT[exposed_hiv_partner] + 1
          
          #SYPHILIS anal
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_sy_anal_ego] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_sy_anal_ego] + sex_exposed$condct_sy_anal[sex_exposed$SY_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_sy_anal_partner] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_sy_anal_partner] + sex_exposed$partcondct_sy_anal[sex_exposed$SY_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_sy_anal_ego] = abm_sims[[current_data]]$SERO_PREVENT[exposed_sy_anal_ego] + sex_exposed$seroct_sy[sex_exposed$SY_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_sy_anal_partner] = abm_sims[[current_data]]$SERO_PREVENT[exposed_sy_anal_partner] + sex_exposed$partseroct_sy[sex_exposed$SY_ANAL==1]
          # mistake in earlier part of code: change NA to 0 
          abm_sims[[current_data]]$DOXYPREP_PREVENT_ANAL[exposed_sy_anal_ego] = ifelse(is.na(abm_sims[[current_data]]$DOXYPREP_PREVENT_ANAL[exposed_sy_anal_ego]), 0 , abm_sims[[current_data]]$DOXYPREP_PREVENT_ANAL[exposed_sy_anal_ego])
          abm_sims[[current_data]]$DOXYPREP_PREVENT_ANAL[exposed_sy_anal_partner] = ifelse(is.na(abm_sims[[current_data]]$DOXYPREP_PREVENT_ANAL[exposed_sy_anal_partner]), 0, abm_sims[[current_data]]$DOXYPREP_PREVENT_ANAL[exposed_sy_anal_partner])
          abm_sims[[current_data]]$DOXYPEP_PREVENT_ANAL[exposed_sy_anal_ego] = ifelse(is.na(abm_sims[[current_data]]$DOXYPEP_PREVENT_ANAL[exposed_sy_anal_ego]),0,abm_sims[[current_data]]$DOXYPEP_PREVENT_ANAL[exposed_sy_anal_ego])
          abm_sims[[current_data]]$DOXYPEP_PREVENT_ANAL[exposed_sy_anal_partner] = ifelse(is.na(abm_sims[[current_data]]$DOXYPEP_PREVENT_ANAL[exposed_sy_anal_partner]),0,abm_sims[[current_data]]$DOXYPEP_PREVENT_ANAL[exposed_sy_anal_partner])
          # update new prevented infections
          abm_sims[[current_data]]$DOXYPEP_PREVENT_ANAL[exposed_sy_anal_ego] = abm_sims[[current_data]]$DOXYPEP_PREVENT_ANAL[exposed_sy_anal_ego] + sex_exposed$doxypepct_anal[sex_exposed$SY_ANAL==1]
          abm_sims[[current_data]]$DOXYPEP_PREVENT_ANAL[exposed_sy_anal_partner] = abm_sims[[current_data]]$DOXYPEP_PREVENT_ANAL[exposed_sy_anal_partner] + sex_exposed$partdoxypepct_anal[sex_exposed$SY_ANAL==1]
          abm_sims[[current_data]]$DOXYPREP_PREVENT_ANAL[exposed_sy_anal_ego] = abm_sims[[current_data]]$DOXYPREP_PREVENT_ANAL[exposed_sy_anal_ego] + sex_exposed$doxyprepct_anal[sex_exposed$SY_ANAL==1]  
          abm_sims[[current_data]]$DOXYPREP_PREVENT_ANAL[exposed_sy_anal_partner] = abm_sims[[current_data]]$DOXYPREP_PREVENT_ANAL[exposed_sy_anal_partner] + sex_exposed$partdoxyprepct_anal[sex_exposed$SY_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_sy_anal_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_sy_anal_ego] + sex_exposed$prevct_sy_anal[sex_exposed$SY_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_sy_anal_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_sy_anal_partner] + sex_exposed$partprevct_sy_anal[sex_exposed$SY_ANAL==1]
          abm_sims[[current_data]]$SY_ANAL[exposed_sy_anal_ego] = ifelse(sex_exposed$newSY_ANAL[sex_exposed$SY_ANAL==1]==1, 1, abm_sims[[current_data]]$SY_ANAL[exposed_sy_anal_ego])
          abm_sims[[current_data]]$SY_ANAL[exposed_sy_anal_partner] = ifelse(sex_exposed$partnewSY_ANAL[sex_exposed$SY_ANAL==1]==1, 1, abm_sims[[current_data]]$SY_ANAL[exposed_sy_anal_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ANAL[exposed_sy_anal_ego] = ifelse(sex_exposed$newSY_ANAL[sex_exposed$SY_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ANAL[exposed_sy_anal_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ANAL[exposed_sy_anal_partner] = ifelse(sex_exposed$partnewSY_ANAL[sex_exposed$SY_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ANAL[exposed_sy_anal_partner])
          abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego] = ifelse(sex_exposed$newSY_ANAL_TYPE[sex_exposed$SY_ANAL==1]==1, 1, abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_ego])
          abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner] = ifelse(sex_exposed$partnewSY_ANAL_TYPE[sex_exposed$SY_ANAL==1]==1, 1, abm_sims[[current_data]]$SY_ANAL_TYPE[exposed_sy_anal_partner])
          
          #SYPHILIS oral
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_sy_oral_ego] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_sy_oral_ego] + sex_exposed$condct_sy_oral[sex_exposed$SY_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_sy_oral_partner] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_sy_oral_partner] + sex_exposed$partcondct_sy_oral[sex_exposed$SY_ORAL==1]
          # mistake in earlier part of code: change NA to 0 
          abm_sims[[current_data]]$DOXYPREP_PREVENT_ORAL[exposed_sy_oral_ego] = ifelse(is.na(abm_sims[[current_data]]$DOXYPREP_PREVENT_ORAL[exposed_sy_oral_ego]), 0 , abm_sims[[current_data]]$DOXYPREP_PREVENT_ORAL[exposed_sy_oral_ego])
          abm_sims[[current_data]]$DOXYPREP_PREVENT_ORAL[exposed_sy_oral_partner] = ifelse(is.na(abm_sims[[current_data]]$DOXYPREP_PREVENT_ORAL[exposed_sy_oral_partner]), 0, abm_sims[[current_data]]$DOXYPREP_PREVENT_ORAL[exposed_sy_oral_partner])
          abm_sims[[current_data]]$DOXYPEP_PREVENT_ORAL[exposed_sy_oral_ego] = ifelse(is.na(abm_sims[[current_data]]$DOXYPEP_PREVENT_ORAL[exposed_sy_oral_ego]),0,abm_sims[[current_data]]$DOXYPEP_PREVENT_ORAL[exposed_sy_oral_ego])
          abm_sims[[current_data]]$DOXYPEP_PREVENT_ORAL[exposed_sy_oral_partner] = ifelse(is.na(abm_sims[[current_data]]$DOXYPEP_PREVENT_ORAL[exposed_sy_oral_partner]),0,abm_sims[[current_data]]$DOXYPEP_PREVENT_ORAL[exposed_sy_oral_partner])
          # update new prevented infections
          abm_sims[[current_data]]$DOXYPEP_PREVENT_ORAL[exposed_sy_oral_ego] = abm_sims[[current_data]]$DOXYPEP_PREVENT_ORAL[exposed_sy_oral_ego] + sex_exposed$doxypepct_oral[sex_exposed$SY_ORAL==1]
          abm_sims[[current_data]]$DOXYPEP_PREVENT_ORAL[exposed_sy_oral_partner] = abm_sims[[current_data]]$DOXYPEP_PREVENT_ORAL[exposed_sy_oral_partner] + sex_exposed$partdoxypepct_oral[sex_exposed$SY_ORAL==1]
          abm_sims[[current_data]]$DOXYPREP_PREVENT_ORAL[exposed_sy_oral_ego] = abm_sims[[current_data]]$DOXYPREP_PREVENT_ORAL[exposed_sy_oral_ego] + sex_exposed$doxyprepct_oral[sex_exposed$SY_ORAL==1]
          abm_sims[[current_data]]$DOXYPREP_PREVENT_ORAL[exposed_sy_oral_partner] = abm_sims[[current_data]]$DOXYPREP_PREVENT_ORAL[exposed_sy_oral_partner] + sex_exposed$partdoxyprepct_oral[sex_exposed$SY_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_sy_oral_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_sy_oral_ego] + sex_exposed$prevct_sy_oral[sex_exposed$SY_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_sy_oral_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_sy_oral_partner] + sex_exposed$partprevct_sy_oral[sex_exposed$SY_ORAL==1]
          abm_sims[[current_data]]$SY_ORAL[exposed_sy_oral_ego] = ifelse(sex_exposed$newSY_ORAL[sex_exposed$SY_ORAL==1]==1 | sex_exposed$newSY_ORAL[sex_exposed$SY_ORAL==1]==3, 1, abm_sims[[current_data]]$SY_ORAL[exposed_sy_oral_ego])
          abm_sims[[current_data]]$SY_ORAL[exposed_sy_oral_partner] = ifelse(sex_exposed$partnewSY_ORAL[sex_exposed$SY_ORAL==1]==1 | sex_exposed$partnewSY_ORAL[sex_exposed$SY_ORAL==1]==3, 1, abm_sims[[current_data]]$SY_ORAL[exposed_sy_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ORAL[exposed_sy_oral_ego] = ifelse(sex_exposed$newSY_ORAL[sex_exposed$SY_ORAL==1]==1 | sex_exposed$newSY_ORAL[sex_exposed$SY_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ORAL[exposed_sy_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ORAL[exposed_sy_oral_partner] = ifelse(sex_exposed$partnewSY_ORAL[sex_exposed$SY_ORAL==1]==1 | sex_exposed$partnewSY_ORAL[sex_exposed$SY_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ORAL[exposed_sy_oral_partner])
          abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego] = ifelse(sex_exposed$newSY_ORAL_TYPE[sex_exposed$SY_ORAL==1]==1, 1, abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_ego])
          abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner] = ifelse(sex_exposed$partnewSY_ORAL_TYPE[sex_exposed$SY_ORAL==1]==1, 1, abm_sims[[current_data]]$SY_ORAL_TYPE[exposed_sy_oral_partner])
          
          # SYPHILIS oral to anal infection??
          # abm_sims[[current_data]]$SY_ANAL[exposed_sy_oral_ego] = ifelse(sex_exposed$newSY_ORAL[sex_exposed$SY_ORAL==1]==2 | sex_exposed$newSY_ORAL[sex_exposed$SY_ORAL==1]==3, 1, abm_sims[[current_data]]$SY_ANAL[exposed_sy_oral_ego])
          # abm_sims[[current_data]]$SY_ANAL[exposed_sy_oral_partner] = ifelse(sex_exposed$partnewSY_ORAL[sex_exposed$SY_ORAL==1]==2 | sex_exposed$partnewSY_ORAL[sex_exposed$SY_ORAL==1]==3, 1, abm_sims[[current_data]]$SY_ANAL[exposed_sy_oral_partner])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ANAL[exposed_sy_oral_ego] = ifelse(sex_exposed$newSY_ORAL[sex_exposed$SY_ORAL==1]==2 | sex_exposed$newSY_ORAL[sex_exposed$SY_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ANAL[exposed_sy_oral_ego])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ANAL[exposed_sy_oral_partner] = ifelse(sex_exposed$partnewSY_ORAL[sex_exposed$SY_ORAL==1]==2 | sex_exposed$partnewSY_ORAL[sex_exposed$SY_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ANAL[exposed_sy_oral_partner])
          
          #update number of sex acts
          abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Ego] = abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Ego] + 1
          abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Partner] = abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Partner] + 1
          
          if (length(sex_unexposed)>0)
            abm_sims[[current_data]]$SEX_COUNT[sex_unexposed] = abm_sims[[current_data]]$SEX_COUNT[sex_unexposed] + 1
          
          ## CLEAN UP ##
          
          rm(sex_exposed,sex_unexposed,exposed_hiv_ego,exposed_hiv_partner,probinfect1_hiv,probinfect2_hiv,exposed_sy_anal_ego,exposed_sy_anal_partner,probinfect1_sy_anal,probinfect2_sy_anal,exposed_sy_oral_ego,exposed_sy_oral_partner,probinfect1_sy_oral,probinfect2_sy_oral)
          #gc()
          
        } else {
          
          #pathogen concordance in every partnership
          
          ## UPDATE STATS IN MAIN DATASET ##
          
          abm_sims[[current_data]]$SEX_COUNT[sex_unexposed] = abm_sims[[current_data]]$SEX_COUNT[sex_unexposed] + 1
          
          ## CLEAN UP ##
          
          rm(sex_exposed,sex_unexposed)
          #gc()
          
        }
        
      } else {
        
        #no partnerships this day
        rm(sex_selection)   
        
      }
      
      #track longitudinal data using incidence days as indicator of new infection
      long_sims[[current_data]]$HIV[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HIV==0)
      long_sims[[current_data]]$SY_ANAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ANAL==0)
      long_sims[[current_data]]$SY_ORAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_SY_ORAL==0)
      
    }
    
  }
  rm(day,current_data,active)
  
  #return this simulation to the list
  abm_sims_50[((sims-1)*nparadigms+1):(((sims-1)*nparadigms)+nparadigms)] = abm_sims
  long_sims_50[((sims-1)*nparadigms+1):(((sims-1)*nparadigms)+nparadigms)] = long_sims
  
}
rm(sims,abm_sims,sex_network,network_sims_50,long_sims,seed_val_day)
gc()

#benchmarking
#time3 = Sys.time()


### SAVE SIMULATION RESULTS ###

save.image("simulation results_25sims.RData")


### STEP5.1: CALIBRATION - HPV Simulation ###

# load("simulation results calibration_10y.RData")          

# #create an annual HIV incidence dataframe
# annual_hiv = data.frame(matrix(data=NA, nrow=0, ncol=(length_sim/365)), stringsAsFactors=F)
# for (i in seq(1,nsims*nparadigms,by=nparadigms))
# {
#   long_sims_50[[i]]$YEAR = (floor((long_sims_50[[i]]$DAY-1)/365) + 1)
#   annual_hiv = rbind(annual_hiv, matrix(data=as.numeric(by(long_sims_50[[i]]$HIV, long_sims_50[[i]]$YEAR, FUN=sum)), nrow=1, ncol=length(unique(long_sims_50[[i]]$YEAR))))
#   
# }
# rm(i)
# 
# colMeans(annual_hiv)
# 
# #plot predicted vs surveillance data
# 
# years = 2013:2017
# surveillance = c(149,138,136,124,126)
# #output to tif for publication 
# #tiff("Figure1.tif",height=6,width=10,units='in',res=1200) 
# plot(x=years, y=surveillance, ylim=c(0,250), pch=18, cex=2, xlab="Surveillance Year", ylab="Reported Cases")
# if (ncol(annual_hiv)>5)
# {
#   annual_hiv = annual_hiv[,(ncol(annual_hiv)-4):ncol(annual_hiv)]
# }
# for (i in 1:length(surveillance)) 
# {
#   arrows(x0=years[i], y0=min(annual_hiv[,i]), x1=years[i], y1=max(annual_hiv[,i]), angle=90, length=0.05, code=3, lwd=2)
# }
# legend("topright", c("surveillance","predicted"), pch=c(18,NA), lty=c(NA,1))
# #close file 
# #dev.off() 
# 
# #check amount of sex
# summary(abm_sims_50[[1]]$MAX_SEX)
# summary(abm_sims_50[[1]]$SEX_COUNT)
# hist(abm_sims_50[[1]]$MAX_SEX,breaks="fd")
# hist(abm_sims_50[[1]]$SEX_COUNT,breaks="fd")
# sum(abm_sims_50[[1]]$SEX_COUNT <= abm_sims_50[[1]]$MAX_SEX)
# sum(abm_sims_50[[1]]$SEX_COUNT > abm_sims_50[[1]]$MAX_SEX)
# boxplot(abm_sims_50[[1]]$SEX_COUNT/abm_sims_50[[1]]$MAX_SEX)
# hist(abm_sims_50[[1]]$PARTNERS,breaks="fd")
# 
# 
# ### STEP6: TALLY RESULTS for PAPER ###
# 
# load("simulation results 045.RData")
# #load("simulation results 045 sens vax95.RData")
# #load("simulation results 045 sens condom50.RData")
# 
# #scenario labels
# hpv_type = c(6, 11, 16, 18, 31, 33, 45, 52, 58)
# vax_level = c(0, 13, 25, 50, 80, 100)
# 
# #create results data frame
# results_infections_anal = data.frame("HPV_type"=NA,"Vax_level"=NA,"Prev_mean"=NA,"Prev_SE"=NA,"Incidence_mean"=NA,"Incidence_SE"=NA,"New_mean"=NA,"New_SE"=NA,"Total_mean"=NA,"Total_SE"=NA,"Prevent_mean"=NA,"Prevent_SE"=NA,"Vax_prevent_mean"=NA,"Vax_prevent_SE"=NA,"Sero_prevent_mean"=NA,"Sero_prevent_SE"=NA,"Cond_prevent_mean"=NA,"Cond_prevent_SE"=NA,stringsAsFactors=F)
# results_infections_oral = data.frame("HPV_type"=NA,"Vax_level"=NA,"Prev_mean"=NA,"Prev_SE"=NA,"Incidence_mean"=NA,"Incidence_SE"=NA,"New_mean"=NA,"New_SE"=NA,"Total_mean"=NA,"Total_SE"=NA,"Prevent_mean"=NA,"Prevent_SE"=NA,"Vax_prevent_mean"=NA,"Vax_prevent_SE"=NA,"Cond_prevent_mean"=NA,"Cond_prevent_SE"=NA,stringsAsFactors=F)
# results_longitudinal = data.frame("HPV_type"=NA,"Vax_level"=NA,"Yr1_anal_new_mean"=NA,"Yr1_anal_new_SE"=NA,"Yr5_anal_new_mean"=NA,"Yr5_anal_new_SE"=NA,"Yr10_anal_new_mean"=NA,"Yr10_anal_new_SE"=NA,"Yr1_oral_new_mean"=NA,"Yr1_oral_new_SE"=NA,"Yr5_oral_new_mean"=NA,"Yr5_oral_new_SE"=NA,"Yr10_oral_new_mean"=NA,"Yr10_oral_new_SE"=NA,stringsAsFactors=F)
# for (i in 1:length(hpv_type))
# {
#   for (j in 1:length(vax_level))
#   {
#     stat_prev_anal = NA         #point prevalance at post simulation
#     stat_incidence_anal = NA    #cumulative incidence (new + spontaneous clearance)
#     stat_new_anal = NA          #new + spontaneous clearance
#     stat_total_anal = NA        #infections at post simulation
#     stat_prevent_anal = NA      #total preventions (aggregate across all serotypes)
#     stat_prevent_vax_anal = NA  #total vax preventions (aggregate across all serotypes)
#     stat_prevent_sero_anal = NA #total sero preventions (aggregate across all serotypes)
#     stat_prevent_cond_anal = NA #total condom preventions (aggregate across all serotypes)
#     stat_yr1_new_anal = NA      #1 yr new infections
#     stat_yr5_new_anal = NA      #5 yr new infections
#     stat_yr10_new_anal = NA      #10 yr new infections
#     stat_prev_oral = NA         #point prevalance at post simulation
#     stat_incidence_oral = NA    #cumulative incidence (new + spontaneous clearance)
#     stat_new_oral = NA          #new + spontaneous clearance
#     stat_total_oral = NA        #infections at post simulation
#     stat_prevent_oral = NA      #total preventions (aggregate across all serotypes)
#     stat_prevent_vax_oral = NA  #total vax preventions (aggregate across all serotypes)
#     stat_prevent_cond_oral = NA #total condom preventions (aggregate across all serotypes)
#     stat_yr1_new_oral = NA      #1 yr new infections
#     stat_yr5_new_oral = NA      #5 yr new infections
#     stat_yr10_new_oral = NA      #10 yr new infections
#     
#     for (k in 0:(nsims-1))
#     {
#       current = k*length(vax_level)
#       
#       if (hpv_type[i]==6) {
#         stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV6_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV6_ANAL) + sum(abm_sims_50[[j+current]]$HPV6_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV6_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV6_ANAL) + sum(abm_sims_50[[j+current]]$HPV6_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV6_ANAL)))
#         stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV6_ANAL)))
#         stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
#         stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
#         stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
#         stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
#         stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV6_ANAL[1:(365*1)])))
#         stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV6_ANAL[1:(365*5)])))
#         stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV6_ANAL[1:(365*10)])))
#         stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV6_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV6_ORAL) + sum(abm_sims_50[[j+current]]$HPV6_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV6_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV6_ORAL) + sum(abm_sims_50[[j+current]]$HPV6_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV6_ORAL)))
#         stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV6_ORAL)))
#         stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
#         stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
#         stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
#         stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV6_ORAL[1:(365*1)])))
#         stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV6_ORAL[1:(365*5)])))
#         stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV6_ORAL[1:(365*10)])))
#       } else if (hpv_type[i]==11) {
#         stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV11_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV11_ANAL) + sum(abm_sims_50[[j+current]]$HPV11_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV11_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV11_ANAL) + sum(abm_sims_50[[j+current]]$HPV11_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV11_ANAL)))
#         stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV11_ANAL)))
#         stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
#         stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
#         stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
#         stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
#         stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV11_ANAL[1:(365*1)])))
#         stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV11_ANAL[1:(365*5)])))
#         stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV11_ANAL[1:(365*10)])))
#         stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV11_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV11_ORAL) + sum(abm_sims_50[[j+current]]$HPV11_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV11_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV11_ORAL) + sum(abm_sims_50[[j+current]]$HPV11_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV11_ORAL)))
#         stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV11_ORAL)))
#         stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
#         stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
#         stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
#         stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV11_ORAL[1:(365*1)])))
#         stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV11_ORAL[1:(365*5)])))
#         stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV11_ORAL[1:(365*10)])))
#       } else if (hpv_type[i]==16) {
#         stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV16_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV16_ANAL) + sum(abm_sims_50[[j+current]]$HPV16_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV16_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV16_ANAL) + sum(abm_sims_50[[j+current]]$HPV16_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV16_ANAL)))
#         stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV16_ANAL)))
#         stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
#         stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
#         stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
#         stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
#         stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV16_ANAL[1:(365*1)])))
#         stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV16_ANAL[1:(365*5)])))
#         stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV16_ANAL[1:(365*10)])))
#         stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV16_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV16_ORAL) + sum(abm_sims_50[[j+current]]$HPV16_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV16_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV16_ORAL) + sum(abm_sims_50[[j+current]]$HPV16_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV16_ORAL)))
#         stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV16_ORAL)))
#         stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
#         stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
#         stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
#         stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV16_ORAL[1:(365*1)])))
#         stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV16_ORAL[1:(365*5)])))
#         stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV16_ORAL[1:(365*10)])))
#       } else if (hpv_type[i]==18) {
#         stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV18_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV18_ANAL) + sum(abm_sims_50[[j+current]]$HPV18_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV18_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV18_ANAL) + sum(abm_sims_50[[j+current]]$HPV18_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV18_ANAL)))
#         stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV18_ANAL)))
#         stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
#         stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
#         stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
#         stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
#         stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV18_ANAL[1:(365*1)])))
#         stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV18_ANAL[1:(365*5)])))
#         stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV18_ANAL[1:(365*10)])))
#         stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV18_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV18_ORAL) + sum(abm_sims_50[[j+current]]$HPV18_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV18_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV18_ORAL) + sum(abm_sims_50[[j+current]]$HPV18_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV18_ORAL)))
#         stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV18_ORAL)))
#         stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
#         stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
#         stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
#         stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV18_ORAL[1:(365*1)])))
#         stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV18_ORAL[1:(365*5)])))
#         stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV18_ORAL[1:(365*10)])))
#       } else if (hpv_type[i]==31) {
#         stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV31_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV31_ANAL) + sum(abm_sims_50[[j+current]]$HPV31_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV31_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV31_ANAL) + sum(abm_sims_50[[j+current]]$HPV31_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV31_ANAL)))
#         stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV31_ANAL)))
#         stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
#         stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
#         stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
#         stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
#         stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV31_ANAL[1:(365*1)])))
#         stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV31_ANAL[1:(365*5)])))
#         stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV31_ANAL[1:(365*10)])))
#         stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV31_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV31_ORAL) + sum(abm_sims_50[[j+current]]$HPV31_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV31_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV31_ORAL) + sum(abm_sims_50[[j+current]]$HPV31_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV31_ORAL)))
#         stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV31_ORAL)))
#         stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
#         stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
#         stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
#         stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV31_ORAL[1:(365*1)])))
#         stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV31_ORAL[1:(365*5)])))
#         stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV31_ORAL[1:(365*10)])))
#       } else if (hpv_type[i]==33) {
#         stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV33_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV33_ANAL) + sum(abm_sims_50[[j+current]]$HPV33_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV33_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV33_ANAL) + sum(abm_sims_50[[j+current]]$HPV33_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV33_ANAL)))
#         stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV33_ANAL)))
#         stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
#         stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
#         stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
#         stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
#         stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV33_ANAL[1:(365*1)])))
#         stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV33_ANAL[1:(365*5)])))
#         stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV33_ANAL[1:(365*10)])))
#         stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV33_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV33_ORAL) + sum(abm_sims_50[[j+current]]$HPV33_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV33_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV33_ORAL) + sum(abm_sims_50[[j+current]]$HPV33_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV33_ORAL)))
#         stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV33_ORAL)))
#         stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
#         stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
#         stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
#         stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV33_ORAL[1:(365*1)])))
#         stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV33_ORAL[1:(365*5)])))
#         stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV33_ORAL[1:(365*10)])))
#       } else if (hpv_type[i]==45) {
#         stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV45_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV45_ANAL) + sum(abm_sims_50[[j+current]]$HPV45_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV45_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV45_ANAL) + sum(abm_sims_50[[j+current]]$HPV45_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV45_ANAL)))
#         stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV45_ANAL)))
#         stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
#         stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
#         stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
#         stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
#         stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV45_ANAL[1:(365*1)])))
#         stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV45_ANAL[1:(365*5)])))
#         stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV45_ANAL[1:(365*10)])))
#         stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV45_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV45_ORAL) + sum(abm_sims_50[[j+current]]$HPV45_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV45_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV45_ORAL) + sum(abm_sims_50[[j+current]]$HPV45_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV45_ORAL)))
#         stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV45_ORAL)))
#         stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
#         stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
#         stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
#         stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV45_ORAL[1:(365*1)])))
#         stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV45_ORAL[1:(365*5)])))
#         stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV45_ORAL[1:(365*10)])))
#       } else if (hpv_type[i]==52) {
#         stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV52_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV52_ANAL) + sum(abm_sims_50[[j+current]]$HPV52_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV52_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV52_ANAL) + sum(abm_sims_50[[j+current]]$HPV52_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV52_ANAL)))
#         stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV52_ANAL)))
#         stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
#         stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
#         stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
#         stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
#         stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV52_ANAL[1:(365*1)])))
#         stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV52_ANAL[1:(365*5)])))
#         stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV52_ANAL[1:(365*10)])))
#         stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV52_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV52_ORAL) + sum(abm_sims_50[[j+current]]$HPV52_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV52_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV52_ORAL) + sum(abm_sims_50[[j+current]]$HPV52_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV52_ORAL)))
#         stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV52_ORAL)))
#         stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
#         stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
#         stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
#         stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV52_ORAL[1:(365*1)])))
#         stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV52_ORAL[1:(365*5)])))
#         stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV52_ORAL[1:(365*10)])))
#       } else if (hpv_type[i]==58) {
#         stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV58_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV58_ANAL) + sum(abm_sims_50[[j+current]]$HPV58_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV58_ANAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV58_ANAL) + sum(abm_sims_50[[j+current]]$HPV58_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV58_ANAL)))
#         stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV58_ANAL)))
#         stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
#         stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
#         stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
#         stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
#         stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV58_ANAL[1:(365*1)])))
#         stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV58_ANAL[1:(365*5)])))
#         stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV58_ANAL[1:(365*10)])))
#         stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV58_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV58_ORAL) + sum(abm_sims_50[[j+current]]$HPV58_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV58_ORAL))/nrow(abm_sims_50[[j+current]]))
#         stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV58_ORAL) + sum(abm_sims_50[[j+current]]$HPV58_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV58_ORAL)))
#         stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV58_ORAL)))
#         stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
#         stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
#         stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
#         stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV58_ORAL[1:(365*1)])))
#         stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV58_ORAL[1:(365*5)])))
#         stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV58_ORAL[1:(365*10)])))
#       }
#     }
#     
#     stat_prev_anal = stat_prev_anal[-1]
#     stat_incidence_anal = stat_incidence_anal[-1]
#     stat_new_anal = stat_new_anal[-1]
#     stat_total_anal = stat_total_anal[-1]
#     stat_prevent_anal = stat_prevent_anal[-1]
#     stat_prevent_vax_anal = stat_prevent_vax_anal[-1]
#     stat_prevent_sero_anal = stat_prevent_sero_anal[-1]
#     stat_prevent_cond_anal = stat_prevent_cond_anal[-1]
#     stat_yr1_new_anal = stat_yr1_new_anal[-1]
#     stat_yr5_new_anal = stat_yr5_new_anal[-1]
#     stat_yr10_new_anal = stat_yr10_new_anal[-1]
#     stat_prev_oral = stat_prev_oral[-1]
#     stat_incidence_oral = stat_incidence_oral[-1]
#     stat_new_oral = stat_new_oral[-1]
#     stat_total_oral = stat_total_oral[-1]
#     stat_prevent_oral = stat_prevent_oral[-1]
#     stat_prevent_vax_oral = stat_prevent_vax_oral[-1]
#     stat_prevent_cond_oral = stat_prevent_cond_oral[-1]
#     stat_yr1_new_oral = stat_yr1_new_oral[-1]
#     stat_yr5_new_oral = stat_yr5_new_oral[-1]
#     stat_yr10_new_oral = stat_yr10_new_oral[-1]
#     
#     results_infections_anal = rbind(results_infections_anal, data.frame("HPV_type"=hpv_type[i],"Vax_level"=vax_level[j],"Prev_mean"=mean(stat_prev_anal),"Prev_SE"=(sd(stat_prev_anal)/sqrt(length(stat_prev_anal))),"Incidence_mean"=mean(stat_incidence_anal),"Incidence_SE"=(sd(stat_incidence_anal)/sqrt(length(stat_incidence_anal))),"New_mean"=mean(stat_new_anal),"New_SE"=(sd(stat_new_anal)/sqrt(length(stat_new_anal))),"Total_mean"=mean(stat_total_anal),"Total_SE"=(sd(stat_total_anal)/sqrt(length(stat_total_anal))),"Prevent_mean"=mean(stat_prevent_anal),"Prevent_SE"=(sd(stat_prevent_anal)/sqrt(length(stat_prevent_anal))),"Vax_prevent_mean"=mean(stat_prevent_vax_anal),"Vax_prevent_SE"=(sd(stat_prevent_vax_anal)/sqrt(length(stat_prevent_vax_anal))),"Sero_prevent_mean"=mean(stat_prevent_sero_anal),"Sero_prevent_SE"=(sd(stat_prevent_sero_anal)/sqrt(length(stat_prevent_sero_anal))),"Cond_prevent_mean"=mean(stat_prevent_cond_anal),"Cond_prevent_SE"=(sd(stat_prevent_cond_anal)/sqrt(length(stat_prevent_cond_anal))),stringsAsFactors=F))
#     results_infections_oral = rbind(results_infections_oral, data.frame("HPV_type"=hpv_type[i],"Vax_level"=vax_level[j],"Prev_mean"=mean(stat_prev_oral),"Prev_SE"=(sd(stat_prev_oral)/sqrt(length(stat_prev_oral))),"Incidence_mean"=mean(stat_incidence_oral),"Incidence_SE"=(sd(stat_incidence_oral)/sqrt(length(stat_incidence_oral))),"New_mean"=mean(stat_new_oral),"New_SE"=(sd(stat_new_oral)/sqrt(length(stat_new_oral))),"Total_mean"=mean(stat_total_oral),"Total_SE"=(sd(stat_total_oral)/sqrt(length(stat_total_oral))),"Prevent_mean"=mean(stat_prevent_oral),"Prevent_SE"=(sd(stat_prevent_oral)/sqrt(length(stat_prevent_oral))),"Vax_prevent_mean"=mean(stat_prevent_vax_oral),"Vax_prevent_SE"=(sd(stat_prevent_vax_oral)/sqrt(length(stat_prevent_vax_oral))),"Cond_prevent_mean"=mean(stat_prevent_cond_oral),"Cond_prevent_SE"=(sd(stat_prevent_cond_oral)/sqrt(length(stat_prevent_cond_oral))),stringsAsFactors=F))
#     results_longitudinal = rbind(results_longitudinal, data.frame("HPV_type"=hpv_type[i],"Vax_level"=vax_level[j],"Yr1_anal_new_mean"=mean(stat_yr1_new_anal),"Yr1_anal_new_SE"=(sd(stat_yr1_new_anal)/sqrt(length(stat_yr1_new_anal))),"Yr5_anal_new_mean"=mean(stat_yr5_new_anal),"Yr5_anal_new_SE"=(sd(stat_yr5_new_anal)/sqrt(length(stat_yr5_new_anal))),"Yr10_anal_new_mean"=mean(stat_yr10_new_anal),"Yr10_anal_new_SE"=(sd(stat_yr10_new_anal)/sqrt(length(stat_yr10_new_anal))),"Yr1_oral_new_mean"=mean(stat_yr1_new_oral),"Yr1_oral_new_SE"=(sd(stat_yr1_new_oral)/sqrt(length(stat_yr1_new_oral))),"Yr5_oral_new_mean"=mean(stat_yr5_new_oral),"Yr5_oral_new_SE"=(sd(stat_yr5_new_oral)/sqrt(length(stat_yr5_new_oral))),"Yr10_oral_new_mean"=mean(stat_yr10_new_oral),"Yr10_oral_new_SE"=(sd(stat_yr10_new_oral)/sqrt(length(stat_yr10_new_oral))),stringsAsFactors=F))
#   }
# }
# rm(i,j,k,current,stat_prev_anal,stat_incidence_anal,stat_new_anal,stat_total_anal,stat_prevent_anal,stat_prevent_vax_anal,stat_prevent_sero_anal,stat_prevent_cond_anal,stat_yr1_new_anal,stat_yr5_new_anal,stat_yr10_new_anal,stat_prev_oral,stat_incidence_oral,stat_new_oral,stat_total_oral,stat_prevent_oral,stat_prevent_vax_oral,stat_prevent_cond_oral,stat_yr1_new_oral,stat_yr5_new_oral,stat_yr10_new_oral)
# results_infections_anal = results_infections_anal[-1, ]
# results_infections_oral = results_infections_oral[-1, ]
# results_longitudinal = results_longitudinal[-1, ]
# 
# #overall starting prevalence (all simulations start out the same)
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV_ANAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV_ORAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HIV[1:nagents_start])/nagents_start
# 
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV6_ANAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV11_ANAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV16_ANAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV18_ANAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV31_ANAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV33_ANAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV45_ANAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV52_ANAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV58_ANAL[1:nagents_start])/nagents_start
# 
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV6_ORAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV11_ORAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV16_ORAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV18_ORAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV31_ORAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV33_ORAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV45_ORAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV52_ORAL[1:nagents_start])/nagents_start
# sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV58_ORAL[1:nagents_start])/nagents_start
# 
# #population at 10yr simulation end
# mean(nagents_end); sd(nagents_end)
# 
# #quantify sex
# sex = NA
# partners = NA
# for (i in 1:nsims) {
#   sex = c(sex, abm_sims_50[[i*nparadigms]]$SEX_COUNT)
#   partners = c(partners, abm_sims_50[[i*nparadigms]]$PARTNERS)
# }
# summary(partners, na.rm=T)
# summary(sex, na.rm=T)
# rm(i,sex,partners)
# 
# #% reduction in 10yr prevalence overall at each vaccination level compared to present day vaccination levels (13%)
# mean((results_infections_anal$Total_mean[results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$Vax_level==25]) / results_infections_anal$Total_mean[results_infections_anal$Vax_level==13]) * 100
# mean((results_infections_anal$Total_mean[results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$Vax_level==50]) / results_infections_anal$Total_mean[results_infections_anal$Vax_level==13]) * 100
# mean((results_infections_anal$Total_mean[results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$Vax_level==13]) * 100
# mean((results_infections_anal$Total_mean[results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$Vax_level==100]) / results_infections_anal$Total_mean[results_infections_anal$Vax_level==13]) * 100
# 
# mean((results_infections_oral$Total_mean[results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$Vax_level==25]) / results_infections_oral$Total_mean[results_infections_oral$Vax_level==13]) * 100
# mean((results_infections_oral$Total_mean[results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$Vax_level==50]) / results_infections_oral$Total_mean[results_infections_oral$Vax_level==13]) * 100
# mean((results_infections_oral$Total_mean[results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$Vax_level==13]) * 100
# mean((results_infections_oral$Total_mean[results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$Vax_level==100]) / results_infections_oral$Total_mean[results_infections_oral$Vax_level==13]) * 100
# 
# #% reduction in 10yr prevalence by serotype at 80% vaccination (healthy people goal) compared to present day vaccination levels (13%)
# mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==6 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==6 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==6 & results_infections_anal$Vax_level==13]) * 100
# mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==11 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==11 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==11 & results_infections_anal$Vax_level==13]) * 100
# mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==16 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==16 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==16 & results_infections_anal$Vax_level==13]) * 100
# mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==18 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==18 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==18 & results_infections_anal$Vax_level==13]) * 100
# mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==31 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==31 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==31 & results_infections_anal$Vax_level==13]) * 100
# mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==33 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==33 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==33 & results_infections_anal$Vax_level==13]) * 100
# mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==45 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==45 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==45 & results_infections_anal$Vax_level==13]) * 100
# mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==52 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==52 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==52 & results_infections_anal$Vax_level==13]) * 100
# mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==58 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==58 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==58 & results_infections_anal$Vax_level==13]) * 100
# 
# mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==6 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==6 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==6 & results_infections_oral$Vax_level==13]) * 100
# mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==11 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==11 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==11 & results_infections_oral$Vax_level==13]) * 100
# mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==16 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==16 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==16 & results_infections_oral$Vax_level==13]) * 100
# mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==18 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==18 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==18 & results_infections_oral$Vax_level==13]) * 100
# mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==31 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==31 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==31 & results_infections_oral$Vax_level==13]) * 100
# mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==33 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==33 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==33 & results_infections_oral$Vax_level==13]) * 100
# mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==45 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==45 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==45 & results_infections_oral$Vax_level==13]) * 100
# mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==52 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==52 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==52 & results_infections_oral$Vax_level==13]) * 100
# mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==58 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==58 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==58 & results_infections_oral$Vax_level==13]) * 100
# 
# #percent reduction in new infections at 1yr, 5yr, 10yr time points at each vaccination level compared to present day vaccination levels (13%)
# summary(((results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==25])) / (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==50])) / (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==80])) / (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==100])) / (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
# 
# summary(((results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==25])) / (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==50])) / (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==80])) / (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==100])) / (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
# 
# summary(((results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==25])) / (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==50])) / (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==80])) / (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==100])) / (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
# 
# summary(((results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==25])) / (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==50])) / (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==80])) / (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==100])) / (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
# 
# summary(((results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==25])) / (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==50])) / (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==80])) / (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==100])) / (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
# 
# summary(((results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==25])) / (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==50])) / (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==80])) / (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
# summary(((results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==100])) / (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
# 
# #% prevention attributed to each strategy, anal
# mean(results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==0]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==0] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==0] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==0]) * 100
# mean(results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==13]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==13] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==13] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==13]) * 100
# mean(results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==25]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==25] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==25] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==25]) * 100
# mean(results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==50]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==50] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==50] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==50]) * 100
# mean(results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==80]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==80] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==80] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==80]) * 100
# mean(results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==100]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==100] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==100] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==100]) * 100
# 
# mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==0]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==0] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==0] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==0]) * 100
# mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==13]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==13] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==13] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==13]) * 100
# mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==25]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==25] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==25] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==25]) * 100
# mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==50]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==50] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==50] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==50]) * 100
# mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==80]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==80] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==80] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==80]) * 100
# mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==100]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==100] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==100] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==100]) * 100
# 
# mean(results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==0]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==0] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==0] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==0]) * 100
# mean(results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==13]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==13] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==13] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==13]) * 100
# mean(results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==25]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==25] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==25] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==25]) * 100
# mean(results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==50]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==50] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==50] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==50]) * 100
# mean(results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==80]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==80] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==80] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==80]) * 100
# mean(results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==100]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==100] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==100] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==100]) * 100
# 
# #% prevention attributed to each strategy, oral
# mean(results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==0]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==0] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==0]) * 100
# mean(results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==13]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==13] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==13]) * 100
# mean(results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==25]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==25] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==25]) * 100
# mean(results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==50]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==50] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==50]) * 100
# mean(results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==80]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==80] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==80]) * 100
# mean(results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==100]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==100] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==100]) * 100
# 
# mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==0]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==0] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==0]) * 100
# mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==13]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==13] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==13]) * 100
# mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==25]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==25] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==25]) * 100
# mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==50]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==50] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==50]) * 100
# mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==80]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==80] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==80]) * 100
# mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==100]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==100] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==100]) * 100


### STEP5.2: CALIBRATION - Syphilis Simulation ###
library(ggplot2); library(ggpubr)
load("simulation results_25sims.RData")

#create an annual HIV incidence dataframe
annual_hiv = data.frame(matrix(data=NA, nrow=0, ncol=(length_sim/365)), stringsAsFactors=F)
for (i in seq(2,nsims*nparadigms,by=nparadigms))
{
  long_sims_50[[i]]$YEAR = (floor((long_sims_50[[i]]$DAY-1)/365) + 1)
  annual_hiv = rbind(annual_hiv, matrix(data=as.numeric(by(long_sims_50[[i]]$HIV, long_sims_50[[i]]$YEAR, FUN=sum)), nrow=1, ncol=length(unique(long_sims_50[[i]]$YEAR))))
  
}
rm(i)

colMeans(annual_hiv)

# hiv calibration
surveillance = c(302,285,303,323,288,310,268,272,207)
annual_hiv1 = annual_hiv[,c(2:10)]
simulation = colMeans(annual_hiv[,c(2:10)])
ll = numeric(9)
for (i in 1:length(ll)) {
  ll[i] = quantile(annual_hiv1[,i], 0.025)
}
ul = numeric(9)
for (i in 1:length(ul)) {
  ul[i] = quantile(annual_hiv1[,i], 0.975)
}

hivdat = data.frame(years=rep(2010:2018), group=c(rep("Simulated",9),rep("Surveillance",9)), estimates=c(simulation,surveillance))
hivdat$ll = c(ll,rep(NA,9))
hivdat$ul = c(ul,rep(NA,9))

# hiv calibration plot
c1 = ggplot(hivdat, aes(years, estimates, linetype = group, shape = group)) +
  geom_line() +
  scale_linetype_manual(values=c("longdash", "solid")) +
  geom_point(size = 2) + 
  scale_shape_manual(values = c(1,16)) +
  geom_pointrange(data = hivdat, aes(ymin=ll, ymax=ul), width = 0.1, linetype = "solid") +
  ylab("Reported Cases") +
  xlab("Surveillance Years") +
  coord_cartesian(ylim = c(0, 600)) +
  ggtitle("HIV Calibration Results") +
  theme_bw() +
  theme(text = element_text(color = "#22211d", size = 10),legend.text=element_text(size=13), plot.title = element_text(size = 15, face = "bold"), 
        axis.title = element_text(size = 13), axis.text = element_text(size = 10), legend.title = element_blank()) 


#create an annual anal syphilis incidence dataframe
annual_sy_anal = data.frame(matrix(data=NA, nrow=0, ncol=(length_sim/365)), stringsAsFactors=F)
for (i in seq(2,nsims*nparadigms,by=nparadigms))
{
  long_sims_50[[i]]$YEAR = (floor((long_sims_50[[i]]$DAY-1)/365) + 1)
  annual_sy_anal = rbind(annual_sy_anal, matrix(data=as.numeric(by(long_sims_50[[i]]$SY_ANAL, long_sims_50[[i]]$YEAR, FUN=sum)), nrow=1, ncol=length(unique(long_sims_50[[i]]$YEAR))))
  
}
rm(i)

colMeans(annual_sy_anal)

#create an annual anal syphilis incidence dataframe
annual_sy_oral = data.frame(matrix(data=NA, nrow=0, ncol=(length_sim/365)), stringsAsFactors=F)
for (i in seq(2,nsims*nparadigms,by=nparadigms))
{
  long_sims_50[[i]]$YEAR = (floor((long_sims_50[[i]]$DAY-1)/365) + 1)
  annual_sy_oral = rbind(annual_sy_oral, matrix(data=as.numeric(by(long_sims_50[[i]]$SY_ORAL, long_sims_50[[i]]$YEAR, FUN=sum)), nrow=1, ncol=length(unique(long_sims_50[[i]]$YEAR))))
  
}
rm(i)

colMeans(annual_sy_oral)

#sum annual anal and oral syphilis 
annual_sy = annual_sy_anal + annual_sy_oral
colMeans(annual_sy)

# syphilis calibration 
surveillance1 = c(218,217,268,248,252,416,568,659,626)
annual_sy1 = annual_sy[,c(2:10)]
simulation1 = colMeans(annual_sy[,c(2:10)])
ll1 = numeric(9)
for (i in 1:length(ll1)) {
  ll1[i] = quantile(annual_sy1[,i], 0.025)
}
ul1 = numeric(9)
for (i in 1:length(ul1)) {
  ul1[i] = quantile(annual_sy1[,i], 0.975)
}

sydat = data.frame(years=rep(2010:2018), group=c(rep("Simulated",9),rep("Surveillance",9)), estimates=c(simulation1,surveillance1))
sydat$ll1 = c(ll1,rep(NA,9))
sydat$ul1 = c(ul1,rep(NA,9))

# syphilis calibration plot
c2 = ggplot(sydat, aes(years, estimates, linetype = group, shape = group)) +
  geom_line() +
  scale_linetype_manual(values=c("longdash", "solid")) +
  geom_point(size = 2) + 
  scale_shape_manual(values = c(1,16)) +
  geom_pointrange(data = sydat, aes(ymin=ll1, ymax=ul1), width = 0.1) +
  ylab("Reported Cases") +
  xlab("Surveillance Years") +
  coord_cartesian(ylim = c(0, 1100)) +
  ggtitle("Syphilis Calibration Results") +
  theme_bw() +
  theme(text = element_text(color = "#22211d", size = 10),legend.text=element_text(size=13), plot.title = element_text(size = 15, face = "bold"), 
        axis.title = element_text(size = 13), axis.text = element_text(size = 10), legend.title = element_blank()) 

# FIG 1. Calibration plots
c3 = ggarrange(c1, c2, ncol=2, nrow=1, labels = c("A", "B"), common.legend = TRUE, legend="top") 
ggsave("calibration_plot.jpeg", plot=c3, width = 10, height = 4, units = "in", dpi = 1200)

# clean up environment 
rm(c1,c2,c3, annual_hiv, annual_hiv1, annual_sy, annual_sy_anal, annual_sy_oral, annual_sy1, hivdat, sydat, i, ll, ll1, ul, ul1, simulation, simulation1, surveillance, surveillance1)


### STEP6. Cumulative incidence Estimation ###
# Separate abm_sims and long_sims into respective doxy PrEP and PEP scenarios 
abm_sims_prep = c(abm_sims_50[c(1:6,14:19,27:32,40:45,53:58,66:71,79:84,92:97,105:110,118:123,
                                131:136,144:149,157:162,170:175,183:188,196:201,209:214,222:227,
                                235:240,248:253,261:266,274:279,287:292,300:305,313:318)+1])

abm_sims_pep = c(abm_sims_50[c(8:13,21:26,34:39,47:52,60:65,73:78,86:91,99:104,112:117,125:130,
                               138:143,151:156,164:169,177:182,190:195,203:208,216:221,229:234,
                               242:247,255:260,268:273,281:286,294:299,307:312,320:325)])

long_sims_prep = c(long_sims_50[c(1:6,14:19,27:32,40:45,53:58,66:71,79:84,92:97,105:110,118:123,
                                  131:136,144:149,157:162,170:175,183:188,196:201,209:214,222:227,
                                  235:240,248:253,261:266,274:279,287:292,300:305,313:318)+1])

long_sims_pep = c(long_sims_50[c(8:13,21:26,34:39,47:52,60:65,73:78,86:91,99:104,112:117,125:130,
                                 138:143,151:156,164:169,177:182,190:195,203:208,216:221,229:234,
                                 242:247,255:260,268:273,281:286,294:299,307:312,320:325)])


### Cumulative incidence 

paradigm = 1:2
uptake_level = c(0,10,20,30,40,50)
day = 365.24

analwide = data.frame("Paradigm"=NA,"Uptake_level"=NA,"yr1"=NA,"yr2"=NA,"yr3"=NA,"yr4"=NA,"yr5"=NA,"yr6"=NA,"yr7"=NA,"yr8"=NA,"yr9"=NA,"yr10"=NA,stringsAsFactors=F)
oralwide = data.frame("Paradigm"=NA,"Uptake_level"=NA,"yr1"=NA,"yr2"=NA,"yr3"=NA,"yr4"=NA,"yr5"=NA,"yr6"=NA,"yr7"=NA,"yr8"=NA,"yr9"=NA,"yr10"=NA,stringsAsFactors=F)


for (i in 1:length(paradigm))
{
  for (j in 1:length(uptake_level))
  {
    # anal
    stat_yr1_new_anal = NA           #1 yr new infections
    stat_yr2_new_anal = NA           #2 yr new infections
    stat_yr3_new_anal = NA           #3 yr new infections
    stat_yr4_new_anal = NA           #4 yr new infections
    stat_yr5_new_anal = NA           #5 yr new infections
    stat_yr6_new_anal = NA           #6 yr new infections
    stat_yr7_new_anal = NA           #7 yr new infections
    stat_yr8_new_anal = NA           #8 yr new infections
    stat_yr9_new_anal = NA           #9 yr new infections
    stat_yr10_new_anal = NA          #10 yr new infections
    # oral
    stat_yr1_new_oral = NA           #1 yr new infections
    stat_yr2_new_oral = NA           #2 yr new infections
    stat_yr3_new_oral = NA           #3 yr new infections
    stat_yr4_new_oral = NA           #4 yr new infections
    stat_yr5_new_oral = NA           #5 yr new infections
    stat_yr6_new_oral = NA           #6 yr new infections
    stat_yr7_new_oral = NA           #7 yr new infections
    stat_yr8_new_oral = NA           #8 yr new infections
    stat_yr9_new_oral = NA           #9 yr new infections
    stat_yr10_new_oral = NA          #10 yr new infections
    
    for (k in 0:(nsims-1))
    {
      current = k*length(uptake_level)
      
      if (paradigm[i]==1) {
        # DOXY PREP
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_prep[[j+current]]$SY_ANAL[(365*9+1):(365*10)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new_anal = c(stat_yr2_new_anal, (sum(long_sims_prep[[j+current]]$SY_ANAL[(365*9+1):(365*11)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new_anal = c(stat_yr3_new_anal, (sum(long_sims_prep[[j+current]]$SY_ANAL[(365*9+1):(365*12)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new_anal = c(stat_yr4_new_anal, (sum(long_sims_prep[[j+current]]$SY_ANAL[(365*9+1):(365*13)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_prep[[j+current]]$SY_ANAL[(365*9+1):(365*14)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new_anal = c(stat_yr6_new_anal, (sum(long_sims_prep[[j+current]]$SY_ANAL[(365*9+1):(365*15)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new_anal = c(stat_yr7_new_anal, (sum(long_sims_prep[[j+current]]$SY_ANAL[(365*9+1):(365*16)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new_anal = c(stat_yr8_new_anal, (sum(long_sims_prep[[j+current]]$SY_ANAL[(365*9+1):(365*17)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new_anal = c(stat_yr9_new_anal, (sum(long_sims_prep[[j+current]]$SY_ANAL[(365*9+1):(365*18)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_prep[[j+current]]$SY_ANAL[(365*9+1):(365*19)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*20))))
        
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_prep[[j+current]]$SY_ORAL[(365*9+1):(365*10)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new_oral = c(stat_yr2_new_oral, (sum(long_sims_prep[[j+current]]$SY_ORAL[(365*9+1):(365*11)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new_oral = c(stat_yr3_new_oral, (sum(long_sims_prep[[j+current]]$SY_ORAL[(365*9+1):(365*12)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new_oral = c(stat_yr4_new_oral, (sum(long_sims_prep[[j+current]]$SY_ORAL[(365*9+1):(365*13)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_prep[[j+current]]$SY_ORAL[(365*9+1):(365*14)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new_oral = c(stat_yr6_new_oral, (sum(long_sims_prep[[j+current]]$SY_ORAL[(365*9+1):(365*15)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new_oral = c(stat_yr7_new_oral, (sum(long_sims_prep[[j+current]]$SY_ORAL[(365*9+1):(365*16)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new_oral = c(stat_yr8_new_oral, (sum(long_sims_prep[[j+current]]$SY_ORAL[(365*9+1):(365*17)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new_oral = c(stat_yr9_new_oral, (sum(long_sims_prep[[j+current]]$SY_ORAL[(365*9+1):(365*18)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_prep[[j+current]]$SY_ORAL[(365*9+1):(365*19)]) / (nagents_start + (((nrow(abm_sims_prep[[j+current]]) - nagents_start)/(length_sim/day))*20))))
      } else if (paradigm[i]==2) {
        # DOXY PEP 
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_pep[[j+current]]$SY_ANAL[(365*9+1):(365*10)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new_anal = c(stat_yr2_new_anal, (sum(long_sims_pep[[j+current]]$SY_ANAL[(365*9+1):(365*11)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new_anal = c(stat_yr3_new_anal, (sum(long_sims_pep[[j+current]]$SY_ANAL[(365*9+1):(365*12)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new_anal = c(stat_yr4_new_anal, (sum(long_sims_pep[[j+current]]$SY_ANAL[(365*9+1):(365*13)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_pep[[j+current]]$SY_ANAL[(365*9+1):(365*14)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new_anal = c(stat_yr6_new_anal, (sum(long_sims_pep[[j+current]]$SY_ANAL[(365*9+1):(365*15)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new_anal = c(stat_yr7_new_anal, (sum(long_sims_pep[[j+current]]$SY_ANAL[(365*9+1):(365*16)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new_anal = c(stat_yr8_new_anal, (sum(long_sims_pep[[j+current]]$SY_ANAL[(365*9+1):(365*17)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new_anal = c(stat_yr9_new_anal, (sum(long_sims_pep[[j+current]]$SY_ANAL[(365*9+1):(365*18)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_pep[[j+current]]$SY_ANAL[(365*9+1):(365*19)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*20))))
        
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_pep[[j+current]]$SY_ORAL[(365*9+1):(365*10)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new_oral = c(stat_yr2_new_oral, (sum(long_sims_pep[[j+current]]$SY_ORAL[(365*9+1):(365*11)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new_oral = c(stat_yr3_new_oral, (sum(long_sims_pep[[j+current]]$SY_ORAL[(365*9+1):(365*12)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new_oral = c(stat_yr4_new_oral, (sum(long_sims_pep[[j+current]]$SY_ORAL[(365*9+1):(365*13)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_pep[[j+current]]$SY_ORAL[(365*9+1):(365*14)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new_oral = c(stat_yr6_new_oral, (sum(long_sims_pep[[j+current]]$SY_ORAL[(365*9+1):(365*15)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new_oral = c(stat_yr7_new_oral, (sum(long_sims_pep[[j+current]]$SY_ORAL[(365*9+1):(365*16)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new_oral = c(stat_yr8_new_oral, (sum(long_sims_pep[[j+current]]$SY_ORAL[(365*9+1):(365*17)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new_oral = c(stat_yr9_new_oral, (sum(long_sims_pep[[j+current]]$SY_ORAL[(365*9+1):(365*18)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_pep[[j+current]]$SY_ORAL[(365*9+1):(365*19)]) / (nagents_start + (((nrow(abm_sims_pep[[j+current]]) - nagents_start)/(length_sim/day))*20))))
      }
    }
    
    # anal
    stat_yr1_new_anal = stat_yr1_new_anal[-1]
    stat_yr2_new_anal = stat_yr2_new_anal[-1]
    stat_yr3_new_anal = stat_yr3_new_anal[-1]
    stat_yr4_new_anal = stat_yr4_new_anal[-1]
    stat_yr5_new_anal = stat_yr5_new_anal[-1]
    stat_yr6_new_anal = stat_yr6_new_anal[-1]
    stat_yr7_new_anal = stat_yr7_new_anal[-1]
    stat_yr8_new_anal = stat_yr8_new_anal[-1]
    stat_yr9_new_anal = stat_yr9_new_anal[-1]
    stat_yr10_new_anal = stat_yr10_new_anal[-1]
    # oral
    stat_yr1_new_oral = stat_yr1_new_oral[-1]
    stat_yr2_new_oral = stat_yr2_new_oral[-1]
    stat_yr3_new_oral = stat_yr3_new_oral[-1]
    stat_yr4_new_oral = stat_yr4_new_oral[-1]
    stat_yr5_new_oral = stat_yr5_new_oral[-1]
    stat_yr6_new_oral = stat_yr6_new_oral[-1]
    stat_yr7_new_oral = stat_yr7_new_oral[-1]
    stat_yr8_new_oral = stat_yr8_new_oral[-1]
    stat_yr9_new_oral = stat_yr9_new_oral[-1]
    stat_yr10_new_oral = stat_yr10_new_oral[-1]
    
    analwide = rbind(analwide, data.frame("Paradigm"=paradigm[i],"Uptake_level"=uptake_level[j],"yr1"=stat_yr1_new_anal,"yr2"=stat_yr2_new_anal,"yr3"=stat_yr3_new_anal,"yr4"=stat_yr4_new_anal,"yr5"=stat_yr5_new_anal,"yr6"=stat_yr6_new_anal,"yr7"=stat_yr7_new_anal,"yr8"=stat_yr8_new_anal,"yr9"=stat_yr9_new_anal,"yr10"=stat_yr10_new_anal,stringsAsFactors=F))
    oralwide = rbind(oralwide, data.frame("Paradigm"=paradigm[i],"Uptake_level"=uptake_level[j],"yr1"=stat_yr1_new_oral,"yr2"=stat_yr2_new_oral,"yr3"=stat_yr3_new_oral,"yr4"=stat_yr4_new_oral,"yr5"=stat_yr5_new_oral,"yr6"=stat_yr6_new_oral,"yr7"=stat_yr7_new_oral,"yr8"=stat_yr8_new_oral,"yr9"=stat_yr9_new_oral,"yr10"=stat_yr10_new_oral,stringsAsFactors=F))
    
  }
}

analwide = analwide[-1,]
oralwide = oralwide[-1,]
rm(paradigm, uptake_level, day, current, i, j, k, stat_yr1_new_anal, stat_yr2_new_anal, stat_yr3_new_anal, stat_yr4_new_anal, stat_yr5_new_anal, stat_yr6_new_anal, stat_yr7_new_anal, stat_yr8_new_anal, stat_yr9_new_anal, stat_yr10_new_anal, stat_yr1_new_oral, stat_yr2_new_oral, stat_yr3_new_oral, stat_yr4_new_oral, stat_yr5_new_oral, stat_yr6_new_oral, stat_yr7_new_oral, stat_yr8_new_oral, stat_yr9_new_oral, stat_yr10_new_oral)

# convert to long form
analLong = reshape(analwide, timevar = "year", times = 2019:2028, v.names = "incidence", varying = list(3:12), direction = "long")
oralLong = reshape(oralwide, timevar = "year", times = 2019:2028, v.names = "incidence", varying = list(3:12), direction = "long")

analLong$year1 = analLong$year - 2018
oralLong$year1 = oralLong$year - 2018

# incidence per 1000
analLong$ir = analLong$incidence * 1000
oralLong$ir = oralLong$incidence * 1000


rm(oralwide, analwide)

# FIG 2. Incidence plots
p1 = ggplot(analLong[analLong$Paradigm==1,], aes(x=as.factor(year1), y=ir, col = as.factor(Uptake_level), group = as.factor(Uptake_level))) +
  geom_pointrange(mapping = aes(x = as.factor(year1), y = ir, group = as.factor(Uptake_level)),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.025)},
                  fun.max = function(z) {quantile(z,0.975)},
                  fun = mean) +
  stat_summary(fun=mean, geom="line", linetype = 1) +
  labs(x="Years", y="Cumulative Incidence per 1,000") + ggtitle("Doxy PrEP for Anal Infections") + 
  scale_color_grey(name = "Doxy Uptake",labels = c("None","10%","20%","30%","40%","50%"), start = 0.0, end = 0.7) +
  theme_bw() +
  theme(text = element_text(color = "#22211d", size = 10),legend.text=element_text(size=13), plot.title = element_text(size = 15, face = "bold"), 
        axis.title = element_text(size = 13), axis.text = element_text(size = 10), legend.title = element_text(size=13), legend.position = "top") +
  guides(col = guide_legend(nrow = 1))


p2 = ggplot(analLong[analLong$Paradigm==2,], aes(x=as.factor(year1), y=ir, col = as.factor(Uptake_level), group = as.factor(Uptake_level))) +
  geom_pointrange(mapping = aes(x = as.factor(year1), y = ir, group = as.factor(Uptake_level)),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.025)},
                  fun.max = function(z) {quantile(z,0.975)},
                  fun = mean) +
  stat_summary(fun=mean, geom="line", linetype = 1) +
  labs(x="Years", y="Cumulative Incidence per 1,000") + ggtitle("Doxy PEP for Anal Infections") +
  scale_color_grey(name = "Doxy Uptake",labels = c("None","10%","20%","30%","40%","50%"), start = 0.0, end = 0.7) +
  theme_bw() +
  theme(text = element_text(color = "#22211d", size = 10),legend.text=element_text(size=13), plot.title = element_text(size = 15, face = "bold"), 
        axis.title = element_text(size = 13), axis.text = element_text(size = 10), legend.title = element_text(size=13), legend.position = "top") +
  guides(col = guide_legend(nrow = 1))

p3 = ggplot(oralLong[oralLong$Paradigm==1,], aes(x=as.factor(year1), y=ir, col = as.factor(Uptake_level), group = as.factor(Uptake_level))) +
  geom_pointrange(mapping = aes(x = as.factor(year1), y = ir, group = as.factor(Uptake_level)),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.025)},
                  fun.max = function(z) {quantile(z,0.975)},
                  fun = mean) +
  stat_summary(fun=mean, geom="line", linetype = 1) +
  labs(x="Years", y="Cumulative Incidence per 1,000") + ggtitle("Doxy PrEP for Oral Infections") + 
  scale_color_grey(name = "Doxy Uptake",labels = c("None","10%","20%","30%","40%","50%"), start = 0.0, end = 0.7) +
  theme_bw() +
  theme(text = element_text(color = "#22211d", size = 10),legend.text=element_text(size=13), plot.title = element_text(size = 15, face = "bold"), 
        axis.title = element_text(size = 13), axis.text = element_text(size = 10), legend.title = element_text(size=13), legend.position = "top") +  
  guides(col = guide_legend(nrow = 1))

p4 = ggplot(oralLong[oralLong$Paradigm==2,], aes(x=as.factor(year1), y=ir, col = as.factor(Uptake_level), group = as.factor(Uptake_level))) +
  geom_pointrange(mapping = aes(x = as.factor(year1), y = ir, group = as.factor(Uptake_level)),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.025)},
                  fun.max = function(z) {quantile(z,0.975)},
                  fun = mean) +
  stat_summary(fun=mean, geom="line", linetype = 1) +
  labs(x="Years", y="Cumulative Incidence per 1,000") + ggtitle("Doxy PEP for Oral Infections") + 
  scale_color_grey(name = "Doxy Uptake",labels = c("None","10%","20%","30%","40%","50%"), start = 0.0, end = 0.7) +
  theme_bw() +
  theme(text = element_text(color = "#22211d", size = 10),legend.text=element_text(size=13), plot.title = element_text(size = 15, face = "bold"), 
        axis.title = element_text(size = 13), axis.text = element_text(size = 10), legend.title = element_text(size=13), legend.position = "top") +  
  guides(col = guide_legend(nrow = 1))


p5 = ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend="top")
ggsave("CI_plot.jpeg", plot=p5, width = 10, height = 10, units = "in", dpi = 1200)

# risk difference in anal infections for each level of doxy prep uptake at year 1, 5, and 10
m1 = lm(ir ~ as.factor(Uptake_level), data = analLong[analLong$Paradigm==1 & analLong$year1==1,])
m2 = lm(ir ~ as.factor(Uptake_level), data = analLong[analLong$Paradigm==1 & analLong$year1==5,])
m3 = lm(ir ~ as.factor(Uptake_level), data = analLong[analLong$Paradigm==1 & analLong$year1==10,])

# calculate % difference
abs(coef(m1)[6]/coef(m1)[1])
abs(coef(m2)[6]/coef(m2)[1])
abs(coef(m3)[6]/coef(m3)[1])

summary(m3)

# risk difference in anal infections for each level of doxy pep uptake at year 1, 5, and 10
m1 = lm(ir ~ as.factor(Uptake_level), data = analLong[analLong$Paradigm==2 & analLong$year1==1,])
m2 = lm(ir ~ as.factor(Uptake_level), data = analLong[analLong$Paradigm==2 & analLong$year1==5,])
m3 = lm(ir ~ as.factor(Uptake_level), data = analLong[analLong$Paradigm==2 & analLong$year1==10,])

# calculate % difference
abs(coef(m1)[6]/coef(m1)[1])
abs(coef(m2)[6]/coef(m2)[1])
abs(coef(m3)[6]/coef(m3)[1])

summary(m3)

# risk difference in oral infections for each level of doxy prep uptake at year 1, 5, and 10
m1 = lm(ir ~ as.factor(Uptake_level), data = oralLong[oralLong$Paradigm==1 & oralLong$year1==1,])
m2 = lm(ir ~ as.factor(Uptake_level), data = oralLong[oralLong$Paradigm==1 & oralLong$year1==5,])
m3 = lm(ir ~ as.factor(Uptake_level), data = oralLong[oralLong$Paradigm==1 & oralLong$year1==10,])

# calculate % difference
abs(coef(m1)[6]/coef(m1)[1])
abs(coef(m2)[6]/coef(m2)[1])
abs(coef(m3)[6]/coef(m3)[1])

summary(m3)

# risk difference in oral infections for each level of doxy pep uptake at year 1, 5, and 10
m1 = lm(ir ~ as.factor(Uptake_level), data = oralLong[oralLong$Paradigm==2 & oralLong$year1==1,])
m2 = lm(ir ~ as.factor(Uptake_level), data = oralLong[oralLong$Paradigm==2 & oralLong$year1==5,])
m3 = lm(ir ~ as.factor(Uptake_level), data = oralLong[oralLong$Paradigm==2 & oralLong$year1==10,])

# calculate % difference
abs(coef(m1)[6]/coef(m1)[1])
abs(coef(m2)[6]/coef(m2)[1])
abs(coef(m3)[6]/coef(m3)[1])

summary(m3)

rm(m1,m2,m3,oralLong,analLong)

### STEP7: Estimate PIA for Condoms and Doxy ###

paradigm = 1:2
uptake_level = c(0,10,20,30,40,50)

results_infections_anal = data.frame("Paradigm"=NA,"Uptake_level"=NA,"Total_mean"=NA,"Total_SE"=NA,"Prev_mean"=NA,"Prev_SE"=NA,"New_mean"=NA,"New_SE"=NA,"Incidence_mean"=NA,"Incidence_SE"=NA,"Prevent_mean"=NA,"Prevent_SE"=NA,"Doxy_prevent_mean"=NA,"Doxy_prevent_SE"=NA,"Sero_prevent_mean"=NA,"Sero_prevent_SE"=NA,"Cond_prevent_mean"=NA,"Cond_prevent_SE"=NA,"Total_prevent"=NA,"Doxy_prevent_pct"=NA,"Doxy_prevent_pct_SE"=NA,"Doxy_prevent_pct_ll"=NA,"Doxy_prevent_pct_ul"=NA,"Cond_prevent_pct"=NA,"Cond_prevent_pct_SE"=NA,"Cond_prevent_pct_ll"=NA,"Cond_prevent_pct_ul"=NA,stringsAsFactors=F)
results_infections_oral = data.frame("Paradigm"=NA,"Uptake_level"=NA,"Total_mean"=NA,"Total_SE"=NA,"Prev_mean"=NA,"Prev_SE"=NA,"New_mean"=NA,"New_SE"=NA,"Incidence_mean"=NA,"Incidence_SE"=NA,"Prevent_mean"=NA,"Prevent_SE"=NA,"Doxy_prevent_mean"=NA,"Doxy_prevent_SE"=NA,"Cond_prevent_mean"=NA,"Cond_prevent_SE"=NA,"Total_prevent"=NA,"Doxy_prevent_pct"=NA,"Doxy_prevent_pct_SE"=NA,"Doxy_prevent_pct_ll"=NA,"Doxy_prevent_pct_ul"=NA,"Cond_prevent_pct"=NA,"Cond_prevent_pct_SE"=NA,"Cond_prevent_pct_ll"=NA,"Cond_prevent_pct_ul"=NA,stringsAsFactors=F)

for (i in 1:length(paradigm))
{
  for (j in 1:length(uptake_level))
  {
    stat_prev_anal = NA               #point prevalance at post simulation
    stat_incidence_anal = NA          #cumulative incidence (primary/secondary + early latent)
    stat_new_anal = NA                #primary/secondary + early latent
    stat_total_anal = NA              #infections at post simulation
    stat_prevent_anal = NA            #total preventions 
    stat_prevent_doxy_anal = NA       #total doxy prep preventions 
    stat_prevent_sero_anal = NA       #total sero preventions 
    stat_prevent_cond_anal = NA       #total condom preventions
    stat_prevent_doxy_anal_pct = NA   #percentage of anal infections prevented due to doxy
    stat_prevent_cond_anal_pct = NA   #percentage of anal infections prevented due to condoms 
    stat_prev_oral = NA               #point prevalance at post simulation
    stat_incidence_oral = NA          #cumulative incidence (primary/secondary + early latent)
    stat_new_oral = NA                #primary/secondary + early latent
    stat_total_oral = NA              #infections at post simulation
    stat_prevent_oral = NA            #total preventions 
    stat_prevent_doxy_oral = NA       #total doxy prep preventions 
    stat_prevent_cond_oral = NA       #total condom preventions 
    stat_prevent_doxy_oral_pct = NA   #percentage of oral infections prevented due to doxy
    stat_prevent_cond_oral_pct = NA   #percentage of oral infections prevented due to condoms 
    
    for (k in 0:(nsims-1))
    {
      current = k*length(uptake_level)
      if (paradigm[i]==1) {
        # DOXY PREP
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_prep[[j+current]]$SY_ANAL)))
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_prep[[j+current]]$SY_ANAL))/nrow(abm_sims_prep[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_prep[[j+current]]$SY_ANAL) - sum(abm_sims_prep[[j+current]]$ORIGINAL_STATUS_SY_ANAL)))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_prep[[j+current]]$SY_ANAL) - sum(abm_sims_prep[[j+current]]$ORIGINAL_STATUS_SY_ANAL))/nrow(abm_sims_prep[[j+current]]))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_prep[[j+current]]$DOXYPREP_PREVENT_ANAL, na.rm = T) + sum(abm_sims_prep[[j+current]]$SERO_PREVENT) + sum(abm_sims_prep[[j+current]]$COND_PREVENT_ANAL)))
        stat_prevent_doxy_anal = c(stat_prevent_doxy_anal, (sum(abm_sims_prep[[j+current]]$DOXYPREP_PREVENT_ANAL, na.rm = T)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_prep[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_prep[[j+current]]$COND_PREVENT_ANAL)))
        stat_prevent_doxy_anal_pct = c(stat_prevent_doxy_anal_pct, ((sum(abm_sims_prep[[j+current]]$DOXYPREP_PREVENT_ANAL, na.rm = T))/(sum(abm_sims_prep[[j+current]]$DOXYPREP_PREVENT_ANAL, na.rm = T) + sum(abm_sims_prep[[j+current]]$SERO_PREVENT) + sum(abm_sims_prep[[j+current]]$COND_PREVENT_ANAL))))
        stat_prevent_cond_anal_pct = c(stat_prevent_cond_anal_pct, ((sum(abm_sims_prep[[j+current]]$COND_PREVENT_ANAL))/(sum(abm_sims_prep[[j+current]]$DOXYPREP_PREVENT_ANAL, na.rm = T) + sum(abm_sims_prep[[j+current]]$SERO_PREVENT) + sum(abm_sims_prep[[j+current]]$COND_PREVENT_ANAL))))
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_prep[[j+current]]$SY_ORAL)))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_prep[[j+current]]$SY_ORAL))/nrow(abm_sims_prep[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_prep[[j+current]]$SY_ORAL) - sum(abm_sims_prep[[j+current]]$ORIGINAL_STATUS_SY_ORAL)))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_prep[[j+current]]$SY_ORAL) - sum(abm_sims_prep[[j+current]]$ORIGINAL_STATUS_SY_ORAL))/nrow(abm_sims_prep[[j+current]]))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_prep[[j+current]]$DOXYPREP_PREVENT_ORAL, na.rm = T) + sum(abm_sims_prep[[j+current]]$COND_PREVENT_ORAL)))
        stat_prevent_doxy_oral = c(stat_prevent_doxy_oral, (sum(abm_sims_prep[[j+current]]$DOXYPREP_PREVENT_ORAL, na.rm = T)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_prep[[j+current]]$COND_PREVENT_ORAL)))
        stat_prevent_doxy_oral_pct = c(stat_prevent_doxy_oral_pct, ((sum(abm_sims_prep[[j+current]]$DOXYPREP_PREVENT_ORAL, na.rm = T))/(sum(abm_sims_prep[[j+current]]$DOXYPREP_PREVENT_ORAL, na.rm = T) + sum(abm_sims_prep[[j+current]]$COND_PREVENT_ORAL))))
        stat_prevent_cond_oral_pct = c(stat_prevent_cond_oral_pct, ((sum(abm_sims_prep[[j+current]]$COND_PREVENT_ORAL))/(sum(abm_sims_prep[[j+current]]$DOXYPREP_PREVENT_ORAL, na.rm = T) + sum(abm_sims_prep[[j+current]]$COND_PREVENT_ORAL))))
      } else if (paradigm[i]==2) {
        # DOXY PEP
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_pep[[j+current]]$SY_ANAL)))
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_pep[[j+current]]$SY_ANAL))/nrow(abm_sims_pep[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_pep[[j+current]]$SY_ANAL) - sum(abm_sims_pep[[j+current]]$ORIGINAL_STATUS_SY_ANAL)))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_pep[[j+current]]$SY_ANAL) - sum(abm_sims_pep[[j+current]]$ORIGINAL_STATUS_SY_ANAL))/nrow(abm_sims_pep[[j+current]]))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_pep[[j+current]]$DOXYPEP_PREVENT_ANAL, na.rm = T) + sum(abm_sims_pep[[j+current]]$SERO_PREVENT) + sum(abm_sims_pep[[j+current]]$COND_PREVENT_ANAL)))
        stat_prevent_doxy_anal = c(stat_prevent_doxy_anal, (sum(abm_sims_pep[[j+current]]$DOXYPEP_PREVENT_ANAL, na.rm = T)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_pep[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_pep[[j+current]]$COND_PREVENT_ANAL)))
        stat_prevent_doxy_anal_pct = c(stat_prevent_doxy_anal_pct, ((sum(abm_sims_pep[[j+current]]$DOXYPEP_PREVENT_ANAL, na.rm = T))/(sum(abm_sims_pep[[j+current]]$DOXYPEP_PREVENT_ANAL, na.rm = T) + sum(abm_sims_pep[[j+current]]$SERO_PREVENT) + sum(abm_sims_pep[[j+current]]$COND_PREVENT_ANAL))))
        stat_prevent_cond_anal_pct = c(stat_prevent_cond_anal_pct, ((sum(abm_sims_pep[[j+current]]$COND_PREVENT_ANAL))/(sum(abm_sims_pep[[j+current]]$DOXYPEP_PREVENT_ANAL, na.rm = T) + sum(abm_sims_pep[[j+current]]$SERO_PREVENT) + sum(abm_sims_pep[[j+current]]$COND_PREVENT_ANAL))))
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_pep[[j+current]]$SY_ORAL)))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_pep[[j+current]]$SY_ORAL))/nrow(abm_sims_pep[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_pep[[j+current]]$SY_ORAL) - sum(abm_sims_pep[[j+current]]$ORIGINAL_STATUS_SY_ORAL)))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_pep[[j+current]]$SY_ORAL) - sum(abm_sims_pep[[j+current]]$ORIGINAL_STATUS_SY_ORAL))/nrow(abm_sims_pep[[j+current]]))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_pep[[j+current]]$DOXYPEP_PREVENT_ORAL, na.rm = T) + sum(abm_sims_pep[[j+current]]$COND_PREVENT_ORAL)))
        stat_prevent_doxy_oral = c(stat_prevent_doxy_oral, (sum(abm_sims_pep[[j+current]]$DOXYPEP_PREVENT_ORAL, na.rm = T)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_pep[[j+current]]$COND_PREVENT_ORAL)))
        stat_prevent_doxy_oral_pct = c(stat_prevent_doxy_oral_pct, ((sum(abm_sims_pep[[j+current]]$DOXYPEP_PREVENT_ORAL, na.rm = T))/(sum(abm_sims_pep[[j+current]]$DOXYPEP_PREVENT_ORAL, na.rm = T) + sum(abm_sims_pep[[j+current]]$COND_PREVENT_ORAL))))
        stat_prevent_cond_oral_pct = c(stat_prevent_cond_oral_pct, ((sum(abm_sims_pep[[j+current]]$COND_PREVENT_ORAL))/(sum(abm_sims_pep[[j+current]]$DOXYPEP_PREVENT_ORAL, na.rm = T) + sum(abm_sims_pep[[j+current]]$COND_PREVENT_ORAL))))
      }
    }
    stat_total_anal = stat_total_anal[-1]
    stat_prev_anal = stat_prev_anal[-1]
    stat_new_anal = stat_new_anal[-1]
    stat_incidence_anal = stat_incidence_anal[-1]
    stat_prevent_anal = stat_prevent_anal[-1]
    stat_prevent_doxy_anal = stat_prevent_doxy_anal[-1]
    stat_prevent_sero_anal = stat_prevent_sero_anal[-1]
    stat_prevent_cond_anal = stat_prevent_cond_anal[-1]
    stat_prevent_doxy_anal_pct = stat_prevent_doxy_anal_pct[-1]
    stat_prevent_cond_anal_pct = stat_prevent_cond_anal_pct[-1]
    stat_total_oral = stat_total_oral[-1]
    stat_prev_oral = stat_prev_oral[-1]
    stat_new_oral = stat_new_oral[-1]
    stat_incidence_oral = stat_incidence_oral[-1]
    stat_prevent_oral = stat_prevent_oral[-1]
    stat_prevent_doxy_oral = stat_prevent_doxy_oral[-1]
    stat_prevent_cond_oral = stat_prevent_cond_oral[-1]
    stat_prevent_doxy_oral_pct = stat_prevent_doxy_oral_pct[-1]
    stat_prevent_cond_oral_pct = stat_prevent_cond_oral_pct[-1]
    
    results_infections_anal = rbind(results_infections_anal, data.frame("Paradigm"=paradigm[i],"Uptake_level"=uptake_level[j],"Total_mean"=mean(stat_total_anal),"Total_SE"=(sd(stat_total_anal)/sqrt(length(stat_total_anal))),"Prev_mean"=mean(stat_prev_anal),"Prev_SE"=(sd(stat_prev_anal)/sqrt(length(stat_prev_anal))),"New_mean"=mean(stat_new_anal),"New_SE"=(sd(stat_new_anal)/sqrt(length(stat_new_anal))),"Incidence_mean"=mean(stat_incidence_anal),"Incidence_SE"=(sd(stat_incidence_anal)/sqrt(length(stat_incidence_anal))),"Prevent_mean"=mean(stat_prevent_anal),"Prevent_SE"=(sd(stat_prevent_anal)/sqrt(length(stat_prevent_anal))),"Doxy_prevent_mean"=mean(stat_prevent_doxy_anal),"Doxy_prevent_SE"=(sd(stat_prevent_doxy_anal)/sqrt(length(stat_prevent_doxy_anal))),"Sero_prevent_mean"=mean(stat_prevent_sero_anal),"Sero_prevent_SE"=(sd(stat_prevent_sero_anal)/sqrt(length(stat_prevent_sero_anal))),"Cond_prevent_mean"=mean(stat_prevent_cond_anal),"Cond_prevent_SE"=(sd(stat_prevent_cond_anal)/sqrt(length(stat_prevent_cond_anal))),"Total_prevent"=mean(stat_prevent_all_anal),"Doxy_prevent_pct"=mean(stat_prevent_doxy_anal_pct),"Doxy_prevent_pct_SE"=(sd(stat_prevent_doxy_anal_pct)/sqrt(length(stat_prevent_doxy_anal_pct))),"Doxy_prevent_pct_ll"=quantile(stat_prevent_doxy_anal_pct, 0.025),"Doxy_prevent_pct_ul"=quantile(stat_prevent_doxy_anal_pct, 0.975),"Cond_prevent_pct"=mean(stat_prevent_cond_anal_pct),"Cond_prevent_pct_SE"=(sd(stat_prevent_cond_anal_pct)/sqrt(length(stat_prevent_cond_anal_pct))),"Cond_prevent_pct_ll"=quantile(stat_prevent_cond_anal_pct,0.025),"Cond_prevent_pct_ul"=quantile(stat_prevent_cond_anal_pct,0.975),stringsAsFactors=F))
    results_infections_oral = rbind(results_infections_oral, data.frame("Paradigm"=paradigm[i],"Uptake_level"=uptake_level[j],"Total_mean"=mean(stat_total_oral),"Total_SE"=(sd(stat_total_oral)/sqrt(length(stat_total_oral))),"Prev_mean"=mean(stat_prev_oral),"Prev_SE"=(sd(stat_prev_oral)/sqrt(length(stat_prev_oral))),"New_mean"=mean(stat_new_oral),"New_SE"=(sd(stat_new_oral)/sqrt(length(stat_new_oral))),"Incidence_mean"=mean(stat_incidence_oral),"Incidence_SE"=(sd(stat_incidence_oral)/sqrt(length(stat_incidence_oral))),"Prevent_mean"=mean(stat_prevent_oral),"Prevent_SE"=(sd(stat_prevent_oral)/sqrt(length(stat_prevent_oral))),"Doxy_prevent_mean"=mean(stat_prevent_doxy_oral),"Doxy_prevent_SE"=(sd(stat_prevent_doxy_oral)/sqrt(length(stat_prevent_doxy_oral))),"Cond_prevent_mean"=mean(stat_prevent_cond_oral),"Cond_prevent_SE"=(sd(stat_prevent_cond_oral)/sqrt(length(stat_prevent_cond_oral))),"Total_prevent"=mean(stat_prevent_all_oral),"Doxy_prevent_pct"=mean(stat_prevent_doxy_oral_pct),"Doxy_prevent_pct_SE"=(sd(stat_prevent_doxy_oral_pct)/sqrt(length(stat_prevent_doxy_oral_pct))),"Doxy_prevent_pct_ll"=quantile(stat_prevent_doxy_oral_pct,0.025),"Doxy_prevent_pct_ul"=quantile(stat_prevent_doxy_oral_pct, 0.975),"Cond_prevent_pct"=mean(stat_prevent_cond_oral_pct),"Cond_prevent_pct_SE"=(sd(stat_prevent_cond_oral_pct)/sqrt(length(stat_prevent_cond_oral_pct))),"Cond_prevent_pct_ll"=quantile(stat_prevent_cond_oral_pct, 0.025),"Cond_prevent_pct_ul"=quantile(stat_prevent_cond_oral_pct,0.975),stringsAsFactors=F))
  }
}

rm(paradigm, uptake_level, current, i, j, k, stat_incidence_anal, stat_new_anal, stat_prev_anal, stat_total_anal, stat_prevent_anal, stat_prevent_cond_anal, stat_prevent_doxy_anal, stat_prevent_sero_anal, stat_prevent_cond_anal_pct, stat_prevent_doxy_anal_pct, stat_incidence_oral, stat_new_oral, stat_prev_oral, stat_total_oral, stat_prevent_oral, stat_prevent_cond_oral, stat_prevent_doxy_oral, stat_prevent_cond_oral_pct, stat_prevent_doxy_oral_pct)
results_infections_anal = results_infections_anal[-1, ]
results_infections_oral = results_infections_oral[-1, ]

### Agents characteristics post-simulation 

#population at 10yr simulation end
delta = (nagents_end - nagents_start)/20
mean(nagents_start + delta*10); sd(nagents_start + delta*10)

#quantify sex
sex = NA
partners = NA
for (i in 1:nsims) {
  sex = c(sex, abm_sims_50[[i*nparadigms]]$SEX_COUNT)
  partners = c(partners, abm_sims_50[[i*nparadigms]]$PARTNERS)
}
summary(partners, na.rm=T)
summary(sex, na.rm=T)
rm(i,sex,partners)
