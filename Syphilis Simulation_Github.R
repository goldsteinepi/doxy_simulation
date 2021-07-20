#################
# Syphilis and doxycycline simulation study
# Citation: Tran NK, Goldstein ND, Welles SL. Countering the Rise of Syphilis: A Role for Post-Exposure Doxycycline Prophylaxis? Manuscript in preparation.
# Note: Simulation datasets may be downloaded from: https://drive.google.com/file/d/14ccrB6zC6JWIOYozXjo05cVvrjeX8HlT/view?usp=sharing
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
    starting_vals = ifelse(dataset$HIV_STATUS==1 & rand_vector<=intervention_level, 1, 0)
  } else if (prevention_var=="HPVDOSES") {
    #HPVDOSES (number of doses of SY vaccine)
    starting_vals = ifelse(rand_vector<=(1-intervention_level), 0, ifelse(rand_vector<=(intervention_level*(4/13)+(1-intervention_level)), 1, ifelse(rand_vector<=(intervention_level*((4+2)/13)+(1-intervention_level)), 2, 3)))
  } else if (prevention_var=="ONDOXY"){
    #ONDOXY (0=not taking doxy pep, 1=taking doxy)
    starting_vals = ifelse(rand_vector<=(1-intervention_level), 0, ifelse(rand_vector<=(intervention_level+(1-intervention_level)),1, 0))
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
nparadigms = 26

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
  baseline_pop = data.frame("ID"=sex_network$ID, "ACTIVE"=NA, "AGE"=NA, "RACE"=NA, "HIV"=NA, "HIV_VL"=NA, "HIV_STATUS"=NA, "SY_ANAL"=NA, "SY_ORAL"=NA, "SY_ANAL_TYPE"=NA, "SY_ORAL_TYPE"=NA, "CIRC"=NA, "SEX_ACT"=NA, "SEX_POSITION_ANAL_PREF"=NA, "SEX_POSITION_ORAL_PREF"=NA, "SEX_POSITION_ANAL"=NA, "SEX_POSITION_ORAL"=NA, "HIV_TEST_DAY"=NA, "SY_TEST_DAY"=NA, "MAX_SEX"=NA, "SEX_COUNT"=0, "PARTNERS"=NA, "ORIGINAL_STATUS_HIV"=NA, "ORIGINAL_STATUS_SY_ANAL"=NA, "ORIGINAL_STATUS_SY_ORAL"=NA, "PREP_PREVENT"=0, "TAP_PREVENT"=0, "SERO_PREVENT"=0, "COND_PREVENT_ANAL"=0, "COND_PREVENT_ORAL"=0, "DOXY_PREVENT_ANAL"=0, "DOXY_PREVENT_ORAL"=0, "OVERALL_PREVENT_ANAL"=0,"OVERALL_PREVENT_ORAL"=0, "DISCORDANT"=0, "CAUSE_INFECT"=0, "DAYS_KNOWN_HIV_POSITIVE"=-1, "INCIDENCE_DAYS_HIV"=-1, "INCIDENCE_DAYS_SY_ANAL"=-1, "INCIDENCE_DAYS_SY_ORAL"=-1, stringsAsFactors=F)
  
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
  
  #day of syphilis testing will be dependent based on HIV status
  #HIV+ MSM: based on  future test probability 1=tested quarterly, 2=tested bi-annually, 3=tested annually, 4=not tested
  testprob = runif(nagents_end[sims+1],0,1)
  futuretest = ifelse(baseline_pop$HIV_STATUS==3 & testprob<=0.22, 1, ifelse(baseline_pop$HIV_STATUS==3 & testprob<=0.43, 2, ifelse(baseline_pop$HIV_STATUS==3 & testprob<=0.71, 3, 4)))
  # HIV- or HIV unknown: only get tested annually 
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
  
  #create prevention paradigm datasets: abcdef where a=prep, b=treatment as prevention, c=condom use, d=seroadaption, e=vaccination, f=doxycycline
  paradigm111111 = baseline_pop
  
  #PrEP, TasP, condoms, seroadaption, vaccination, doxy PEP
  paradigm111111$PREP = 1 
  paradigm111111$TAP = 1
  paradigm111111$COND = 1
  paradigm111111$SERO = 1
  paradigm111111$VAX = 1
  paradigm111111$DOXY = 1
  paradigm111111$SYTREAT = 1 
  paradigm111111$ONPREP = prevention_paradigm(randprep,"ONPREP",paradigm111111,0.184)
  paradigm111111$TAPVL = prevention_paradigm(randtap,"TAPVL",paradigm111111)
  paradigm111111$PROBCOND = prevention_paradigm(randcond,"PROBCOND",paradigm111111)
  paradigm111111$SEROTYPE = prevention_paradigm(randsero,"SEROTYPE",paradigm111111)
  
  
  #assign HIV PrEP and doxy into lists
  abm_sims = list("paradigm111111doxy00"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0)),
                  "paradigm111111doxy11"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.2)),
                  "paradigm111111doxy12"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.2)),
                  "paradigm111111doxy13"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.2)),
                  "paradigm111111doxy14"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.2)),
                  "paradigm111111doxy15"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.2)),
                  "paradigm111111doxy21"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.4)),
                  "paradigm111111doxy22"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.4)),
                  "paradigm111111doxy23"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.4)),
                  "paradigm111111doxy24"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.4)),
                  "paradigm111111doxy25"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.4)),
                  "paradigm111111doxy31"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.6)),
                  "paradigm111111doxy32"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.6)),
                  "paradigm111111doxy33"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.6)),
                  "paradigm111111doxy34"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.6)),
                  "paradigm111111doxy35"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.6)),
                  "paradigm111111doxy41"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.8)),
                  "paradigm111111doxy42"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.8)),
                  "paradigm111111doxy43"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.8)),
                  "paradigm111111doxy44"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.8)),
                  "paradigm111111doxy45"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,0.8)),
                  "paradigm111111doxy51"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,1)),
                  "paradigm111111doxy52"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,1)),
                  "paradigm111111doxy53"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,1)),
                  "paradigm111111doxy54"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,1)),
                  "paradigm111111doxy55"=cbind(paradigm111111,"ONDOXY"=prevention_paradigm(randdoxy,"ONDOXY",paradigm111111,1)))
  
  # determine doxycycline adherence for each level of uptake
  abm_sims$paradigm111111doxy00$ADHERE = 0
  abm_sims$paradigm111111doxy11$ADHERE = ifelse(abm_sims$paradigm111111doxy11$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.2, 1, 0)
  abm_sims$paradigm111111doxy12$ADHERE = ifelse(abm_sims$paradigm111111doxy12$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.4, 1, 0)
  abm_sims$paradigm111111doxy13$ADHERE = ifelse(abm_sims$paradigm111111doxy13$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.6, 1, 0)
  abm_sims$paradigm111111doxy14$ADHERE = ifelse(abm_sims$paradigm111111doxy14$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.8, 1, 0)
  abm_sims$paradigm111111doxy15$ADHERE = ifelse(abm_sims$paradigm111111doxy15$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 1.0, 1, 0)
  abm_sims$paradigm111111doxy21$ADHERE = ifelse(abm_sims$paradigm111111doxy21$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.2, 1, 0)
  abm_sims$paradigm111111doxy22$ADHERE = ifelse(abm_sims$paradigm111111doxy22$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.4, 1, 0)
  abm_sims$paradigm111111doxy23$ADHERE = ifelse(abm_sims$paradigm111111doxy23$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.6, 1, 0)
  abm_sims$paradigm111111doxy24$ADHERE = ifelse(abm_sims$paradigm111111doxy24$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.8, 1, 0)
  abm_sims$paradigm111111doxy25$ADHERE = ifelse(abm_sims$paradigm111111doxy25$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 1.0, 1, 0)
  abm_sims$paradigm111111doxy31$ADHERE = ifelse(abm_sims$paradigm111111doxy31$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.2, 1, 0)
  abm_sims$paradigm111111doxy32$ADHERE = ifelse(abm_sims$paradigm111111doxy32$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.4, 1, 0)
  abm_sims$paradigm111111doxy33$ADHERE = ifelse(abm_sims$paradigm111111doxy33$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.6, 1, 0)
  abm_sims$paradigm111111doxy34$ADHERE = ifelse(abm_sims$paradigm111111doxy34$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.8, 1, 0)
  abm_sims$paradigm111111doxy35$ADHERE = ifelse(abm_sims$paradigm111111doxy35$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 1.0, 1, 0)
  abm_sims$paradigm111111doxy41$ADHERE = ifelse(abm_sims$paradigm111111doxy41$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.2, 1, 0)
  abm_sims$paradigm111111doxy42$ADHERE = ifelse(abm_sims$paradigm111111doxy42$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.4, 1, 0)
  abm_sims$paradigm111111doxy43$ADHERE = ifelse(abm_sims$paradigm111111doxy43$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.6, 1, 0)
  abm_sims$paradigm111111doxy44$ADHERE = ifelse(abm_sims$paradigm111111doxy44$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.8, 1, 0)
  abm_sims$paradigm111111doxy45$ADHERE = ifelse(abm_sims$paradigm111111doxy45$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 1.0, 1, 0)
  abm_sims$paradigm111111doxy51$ADHERE = ifelse(abm_sims$paradigm111111doxy51$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.2, 1, 0)
  abm_sims$paradigm111111doxy52$ADHERE = ifelse(abm_sims$paradigm111111doxy52$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.4, 1, 0)
  abm_sims$paradigm111111doxy53$ADHERE = ifelse(abm_sims$paradigm111111doxy53$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.6, 1, 0)
  abm_sims$paradigm111111doxy54$ADHERE = ifelse(abm_sims$paradigm111111doxy54$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 0.8, 1, 0)
  abm_sims$paradigm111111doxy55$ADHERE = ifelse(abm_sims$paradigm111111doxy55$ONDOXY==1 & runif(nagents_end[sims+1],0,1) <= 1.0, 1, 0)
  
  #initialize a list for longitudinal tracking of outcomes; must match length of abm_sims
  long_sims = list("paradigm111111doxy00"=data.frame("DAY"=1:length_sim,"HIV"=NA, "SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy11"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy12"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy13"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy14"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy15"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy21"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy22"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy23"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy24"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy25"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy31"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy32"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy33"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy34"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy35"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy41"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy42"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy43"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy44"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy45"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy51"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy52"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy53"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy54"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F),
                   "paradigm111111doxy55"=data.frame("DAY"=1:length_sim,"HIV"=NA,"SY_ANAL"=NA, "SY_ORAL"=NA, stringsAsFactors=F))
  
  #cleanup
  rm(randtap,randcond,randsero,randprep,randvax,randdoxy,baseline_pop,paradigm111111,simnw)
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
# t = table(abm_sims_50$simulation_1_paradigm11111101doxy0$SY_ANAL); t
# prop.table(t)
# t1 = table(abm_sims_50$simulation_1_paradigm11111101doxy0$SY_ANAL, abm_sims_50$simulation_1_paradigm11111101doxy0$AGE); t1
# prop.table(t1)
# t2 = table(abm_sims_50$simulation_1_paradigm11111101doxy0$SY_ANAL, abm_sims_50$simulation_1_paradigm11111101doxy0$RACE); t2
# prop.table(t2)
# t3 = table(abm_sims_50$simulation_1_paradigm11111101doxy0$SY_ORAL); t3
# prop.table(t3)
# t4 = table(abm_sims_50$simulation_1_paradigm11111101doxy0$SY_ORAL, abm_sims_50$simulation_1_paradigm11111101doxy0$AGE); t4
# prop.table(t4)
# t5 = table(abm_sims_50$simulation_1_paradigm11111101doxy0$SY_ORAL, abm_sims_50$simulation_1_paradigm11111101doxy0$RACE); t5
# prop.table(t5)


### SAVE INITIALIZATION: POPULATION ###

save.image("simulation initialization population_20y_25sim.RData")


### STEP4: SIMULATION ###

load("simulation initialization population_20y_25sim.RData")

#ensure that all runs start equivalently
set.seed(777)

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
      if (unique(abm_sims[[current_data]]$SYTREAT)==1) {
        abm_sims[[current_data]]$SY_ANAL = ifelse((!is.na(abm_sims[[current_data]]$SY_TEST_DAY) & ((day-abm_sims[[current_data]]$SY_TEST_DAY) %% 365)==0) & (abm_sims[[current_data]]$SY_ANAL!=0) & runif(nrow(abm_sims[[current_data]]),0,1)<=0.95, 0, abm_sims[[current_data]]$SY_ANAL)
        abm_sims[[current_data]]$SY_ORAL = ifelse((!is.na(abm_sims[[current_data]]$SY_TEST_DAY) & ((day-abm_sims[[current_data]]$SY_TEST_DAY) %% 365)==0) & (abm_sims[[current_data]]$SY_ORAL!=0) & runif(nrow(abm_sims[[current_data]]),0,1)<=0.95, 0, abm_sims[[current_data]]$SY_ORAL)
      }
      
      ## ENGAGE IN SEX ##
      
      #active partnerships today
      sex_selection = data.frame("Ego"=sex_network$ID[which(!is.na(sex_network[,6+day]))], "Partner"=as.numeric(na.omit(sex_network[,6+day])), stringsAsFactors=F)
      
      #sex today, multiplier is to make the algorithm more conservative from calibration (i.e. less sex)
      sex_selection$Ego_sex = (sex_network$N_SEX_ACTS[sex_selection$Ego] / sex_network$RELATIONSHIP_DAYS[sex_selection$Ego]) <= (runif(nrow(sex_selection),0,1)*1.0)
      sex_selection$Partner_sex = (sex_network$N_SEX_ACTS[sex_selection$Partner] / sex_network$RELATIONSHIP_DAYS[sex_selection$Partner]) <= (runif(nrow(sex_selection),0,1)*1.0)
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
          sex_exposed$doxyct_anal = 0           #SY anal infection blocked from doxycycline, ego
          sex_exposed$doxyct_oral = 0           #SY oral infection blocked from doxycycline, ego
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
          sex_exposed$partdoxyct_anal = 0       #SY anal infection blocked from doxycycline, partner
          sex_exposed$partdoxyct_oral = 0       #SY oral infection blocked from doxycycline, partner
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
            seroblock_hiv = ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0010, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==0 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0059, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0101, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==0 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0569, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0010, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==0 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0059, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0101, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==0 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0569, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0134, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.1284, 1, 0))))))))))
            partseroblock_hiv = ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0010, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==0 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0059, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0101, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==0 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0569, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0010, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==0 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0059, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0101, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==0 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0569, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0134, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.1284, 1, 0))))))))))
            
            #check if infection blocked for HIV
            sex_exposed$seroct_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==1, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & seroblock_hiv==1, 1, 0))
            sex_exposed$partseroct_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==1, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & partseroblock_hiv==1, 1, 0))
            
            #seroadaption strategies for HIV, SYPHILIS implications; NAs (and possibly 0s) get returned for HIV infections without correspoding SY infection
            seroblock_sy = ifelse(sex_exposed$infection_sy_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_sy_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==0 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_sy_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_sy_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==0 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_sy_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            partseroblock_sy = ifelse(sex_exposed$partinfection_sy_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_sy_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==0 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_sy_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_sy_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==0 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_sy_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & probinfect1_sy_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            
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
          
          #determine infection under doxy adherence scenarios  
          if (unique(abm_sims[[current_data]]$DOXY)==1)
          {
            #calculate probabilities of doxycycline success
            randdoxy = runif(nrow(sex_exposed),0,1)
            
            #calculate probabilities of doxycycline treatment efficacy
            betadoxy = rbeta(nrow(sex_exposed), shape1 = 6.621, shape2 = 2.405)
            
            #check if infection blocked
            sex_exposed$doxyct_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1]==1 & abm_sims[[current_data]]$ADHERE[exposed_sy_anal_ego]==1 & randdoxy[sex_exposed$SY_ANAL==1]<=betadoxy[sex_exposed$SY_ANAL==1], 1, 0)
            sex_exposed$partdoxyct_anal[sex_exposed$SY_ANAL==1] = ifelse(sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1]==1 & abm_sims[[current_data]]$ADHERE[exposed_sy_anal_partner]==1 & randdoxy[sex_exposed$SY_ANAL==1]<=betadoxy[sex_exposed$SY_ANAL==1], 1, 0)
            sex_exposed$doxyct_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$ADHERE[exposed_sy_oral_ego]==1 & randdoxy[sex_exposed$SY_ORAL==1]<=betadoxy[sex_exposed$SY_ORAL==1], 1, 0)
            sex_exposed$partdoxyct_oral[sex_exposed$SY_ORAL==1] = ifelse(sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1]==1 & abm_sims[[current_data]]$ADHERE[exposed_sy_oral_partner]==1 & randdoxy[sex_exposed$SY_ORAL==1]<=betadoxy[sex_exposed$SY_ORAL==1], 1, 0)
            
            rm(randdoxy,betadoxy)
          }
          
          
          #resolve infection statistics, HIV
          sex_exposed$newHIV[sex_exposed$HIV==1] = ifelse((sex_exposed$condct_hiv[sex_exposed$HIV==1] + sex_exposed$prepct[sex_exposed$HIV==1] + sex_exposed$tapct[sex_exposed$HIV==1] + sex_exposed$seroct_hiv[sex_exposed$HIV==1])==0 & sex_exposed$infection_hiv[sex_exposed$HIV==1]==1, 1, sex_exposed$newHIV[sex_exposed$HIV==1])
          sex_exposed$newVL[sex_exposed$HIV==1] = ifelse((sex_exposed$condct_hiv[sex_exposed$HIV==1] + sex_exposed$prepct[sex_exposed$HIV==1] + sex_exposed$tapct[sex_exposed$HIV==1] + sex_exposed$seroct_hiv[sex_exposed$HIV==1])==0 & sex_exposed$infection_hiv[sex_exposed$HIV==1]==1, 2, sex_exposed$newVL[sex_exposed$HIV==1])
          sex_exposed$partnewHIV[sex_exposed$HIV==1] = ifelse((sex_exposed$partcondct_hiv[sex_exposed$HIV==1] + sex_exposed$partprepct[sex_exposed$HIV==1] + sex_exposed$parttapct[sex_exposed$HIV==1] + sex_exposed$partseroct_hiv[sex_exposed$HIV==1])==0 & sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1, 1, sex_exposed$partnewHIV[sex_exposed$HIV==1])
          sex_exposed$partnewVL[sex_exposed$HIV==1] = ifelse((sex_exposed$partcondct_hiv[sex_exposed$HIV==1] + sex_exposed$partprepct[sex_exposed$HIV==1] + sex_exposed$parttapct[sex_exposed$HIV==1] + sex_exposed$partseroct_hiv[sex_exposed$HIV==1])==0 & sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1, 2, sex_exposed$partnewVL[sex_exposed$HIV==1])
          
          #resolve infection statistics, SYPHILIS anal
          sex_exposed$newSY_ANAL[sex_exposed$SY_ANAL==1] = ifelse((sex_exposed$condct_sy_anal[sex_exposed$SY_ANAL==1] + sex_exposed$seroct_sy[sex_exposed$SY_ANAL==1] + sex_exposed$doxyct_anal[sex_exposed$SY_ANAL==1])==0 & sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1]==1, 1, sex_exposed$newSY_ANAL[sex_exposed$SY_ANAL==1])
          sex_exposed$newSY_ANAL_TYPE[sex_exposed$SY_ANAL==1] = ifelse((sex_exposed$condct_sy_anal[sex_exposed$SY_ANAL==1] + sex_exposed$seroct_sy[sex_exposed$SY_ANAL==1] + sex_exposed$doxyct_anal[sex_exposed$SY_ANAL==1])==0 & sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1]==1, 1, sex_exposed$newSY_ANAL_TYPE[sex_exposed$SY_ANAL==1])
          sex_exposed$prevct_sy_anal[sex_exposed$SY_ANAL==1] = ifelse((sex_exposed$condct_sy_anal[sex_exposed$SY_ANAL==1] + sex_exposed$seroct_sy[sex_exposed$SY_ANAL==1] + sex_exposed$doxyct_anal[sex_exposed$SY_ANAL==1])>0 & sex_exposed$infection_sy_anal[sex_exposed$SY_ANAL==1]==1, 1, sex_exposed$prevct_sy_anal[sex_exposed$SY_ANAL==1])
          sex_exposed$partnewSY_ANAL[sex_exposed$SY_ANAL==1] = ifelse((sex_exposed$partcondct_sy_anal[sex_exposed$SY_ANAL==1] + sex_exposed$partseroct_sy[sex_exposed$SY_ANAL==1] + sex_exposed$partdoxyct_anal[sex_exposed$SY_ANAL==1])==0 & sex_exposed$partinfection_sy_anal[sex_exposed$SY_ANAL==1]==1, 1, sex_exposed$partnewSY_ANAL[sex_exposed$SY_ANAL==1])
          sex_exposed$partnewSY_ANAL_TYPE[sex_exposed$SY_ANAL==1] = ifelse((sex_exposed$partcondct_sy_anal[sex_exposed$SY_ANAL==1] + sex_exposed$partseroct_sy[sex_exposed$SY_ANAL==1] + sex_exposed$partdoxyct_anal[sex_exposed$SY_ANAL==1])==0 & sex_exposed$partinfection_sy_anal[sex_exposed$SY_ANAL==1]==1, 1, sex_exposed$partnewSY_ANAL_TYPE[sex_exposed$SY_ANAL==1])
          sex_exposed$partprevct_sy_anal[sex_exposed$SY_ANAL==1] = ifelse((sex_exposed$partcondct_sy_anal[sex_exposed$SY_ANAL==1] + sex_exposed$partseroct_sy[sex_exposed$SY_ANAL==1] + sex_exposed$partdoxyct_anal[sex_exposed$SY_ANAL==1])>0 & sex_exposed$partinfection_sy_anal[sex_exposed$SY_ANAL==1]==1, 1, sex_exposed$partprevct_sy_anal[sex_exposed$SY_ANAL==1])
          
          #resolve infection statistics, SYPHILIS oral
          sex_exposed$newSY_ORAL[sex_exposed$SY_ORAL==1] = ifelse((sex_exposed$condct_sy_oral[sex_exposed$SY_ORAL==1] + sex_exposed$doxyct_oral[sex_exposed$SY_ORAL==1])==0 & sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1]>=1, sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1], sex_exposed$newSY_ORAL[sex_exposed$SY_ORAL==1])
          sex_exposed$newSY_ORAL_TYPE[sex_exposed$SY_ORAL==1] = ifelse((sex_exposed$condct_sy_oral[sex_exposed$SY_ORAL==1] + sex_exposed$doxyct_oral[sex_exposed$SY_ORAL==1])==0 & sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1]==1, 1, sex_exposed$newSY_ORAL_TYPE[sex_exposed$SY_ORAL==1])
          sex_exposed$prevct_sy_oral[sex_exposed$SY_ORAL==1] = ifelse((sex_exposed$condct_sy_oral[sex_exposed$SY_ORAL==1] + sex_exposed$doxyct_oral[sex_exposed$SY_ORAL==1])>0 & sex_exposed$infection_sy_oral[sex_exposed$SY_ORAL==1]>=1, 1, sex_exposed$prevct_sy_oral[sex_exposed$SY_ORAL==1])
          sex_exposed$partnewSY_ORAL[sex_exposed$SY_ORAL==1] = ifelse((sex_exposed$partcondct_sy_oral[sex_exposed$SY_ORAL==1] + sex_exposed$partdoxyct_oral[sex_exposed$SY_ORAL==1])==0 & sex_exposed$partinfection_sy_oral[sex_exposed$SY_ORAL==1]>=1, sex_exposed$partinfection_sy_oral[sex_exposed$SY_ORAL==1], sex_exposed$partnewSY_ORAL[sex_exposed$SY_ORAL==1])
          sex_exposed$partnewSY_ORAL_TYPE[sex_exposed$SY_ORAL==1] = ifelse((sex_exposed$partcondct_sy_oral[sex_exposed$SY_ORAL==1] + sex_exposed$partdoxyct_oral[sex_exposed$SY_ORAL==1])==0 & sex_exposed$partinfection_sy_oral[sex_exposed$SY_ORAL==1]==1, 1, sex_exposed$partnewSY_ORAL_TYPE[sex_exposed$SY_ORAL==1])
          sex_exposed$partprevct_sy_oral[sex_exposed$SY_ORAL==1] = ifelse((sex_exposed$partcondct_sy_oral[sex_exposed$SY_ORAL==1] + sex_exposed$partdoxyct_oral[sex_exposed$SY_ORAL==1])>0 & sex_exposed$partinfection_sy_oral[sex_exposed$SY_ORAL==1]>=1, 1, sex_exposed$partprevct_sy_oral[sex_exposed$SY_ORAL==1])
          
          
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
          abm_sims[[current_data]]$DOXY_PREVENT_ANAL[exposed_sy_anal_ego] = ifelse(is.na(abm_sims[[current_data]]$DOXY_PREVENT_ANAL[exposed_sy_anal_ego]), 0 , abm_sims[[current_data]]$DOXY_PREVENT_ANAL[exposed_sy_anal_ego])
          abm_sims[[current_data]]$DOXY_PREVENT_ANAL[exposed_sy_anal_partner] = ifelse(is.na(abm_sims[[current_data]]$DOXY_PREVENT_ANAL[exposed_sy_anal_partner]), 0, abm_sims[[current_data]]$DOXY_PREVENT_ANAL[exposed_sy_anal_partner])
          # update new prevented infections
          abm_sims[[current_data]]$DOXY_PREVENT_ANAL[exposed_sy_anal_ego] = abm_sims[[current_data]]$DOXY_PREVENT_ANAL[exposed_sy_anal_ego] + sex_exposed$doxyct_anal[sex_exposed$SY_ANAL==1]
          abm_sims[[current_data]]$DOXY_PREVENT_ANAL[exposed_sy_anal_partner] = abm_sims[[current_data]]$DOXY_PREVENT_ANAL[exposed_sy_anal_partner] + sex_exposed$partdoxyct_anal[sex_exposed$SY_ANAL==1]
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
          abm_sims[[current_data]]$DOXY_PREVENT_ORAL[exposed_sy_oral_ego] = ifelse(is.na(abm_sims[[current_data]]$DOXY_PREVENT_ORAL[exposed_sy_oral_ego]), 0 , abm_sims[[current_data]]$DOXY_PREVENT_ORAL[exposed_sy_oral_ego])
          abm_sims[[current_data]]$DOXY_PREVENT_ORAL[exposed_sy_oral_partner] = ifelse(is.na(abm_sims[[current_data]]$DOXY_PREVENT_ORAL[exposed_sy_oral_partner]), 0, abm_sims[[current_data]]$DOXY_PREVENT_ORAL[exposed_sy_oral_partner])
          # update new prevented infections
          abm_sims[[current_data]]$DOXY_PREVENT_ORAL[exposed_sy_oral_ego] = abm_sims[[current_data]]$DOXY_PREVENT_ORAL[exposed_sy_oral_ego] + sex_exposed$doxyct_oral[sex_exposed$SY_ORAL==1]
          abm_sims[[current_data]]$DOXY_PREVENT_ORAL[exposed_sy_oral_partner] = abm_sims[[current_data]]$DOXY_PREVENT_ORAL[exposed_sy_oral_partner] + sex_exposed$partdoxyct_oral[sex_exposed$SY_ORAL==1]
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
library(ggplot2); library(ggpubr); library(directlabels)
load("simulation results_25sims.RData")

#create an annual HIV incidence dataframe
annual_hiv = data.frame(matrix(data=NA, nrow=0, ncol=(length_sim/365)), stringsAsFactors=F)
for (i in seq(1,nsims*nparadigms,by=nparadigms))
{
  long_sims_50[[i]]$YEAR = (floor((long_sims_50[[i]]$DAY-1)/365) + 1)
  annual_hiv = rbind(annual_hiv, matrix(data=as.numeric(by(long_sims_50[[i]]$HIV, long_sims_50[[i]]$YEAR, FUN=sum)), nrow=1, ncol=length(unique(long_sims_50[[i]]$YEAR))))
  
}
rm(i)

colMeans(annual_hiv)

#create an annual syphilis incidence dataframe
annual_sy = data.frame(matrix(data=NA, nrow=0, ncol=(length_sim/365)), stringsAsFactors=F)
for (i in seq(1,nsims*nparadigms,by=nparadigms))
{
  long_sims_50[[i]]$YEAR = (floor((long_sims_50[[i]]$DAY-1)/365) + 1)
  annual_sy = rbind(annual_sy, matrix(data=as.numeric(by(rowSums(cbind(long_sims_50[[i]]$SY_ORAL, long_sims_50[[i]]$SY_ANAL)), long_sims_50[[i]]$YEAR, FUN=sum)), nrow=1, ncol=length(unique(long_sims_50[[i]]$YEAR))))
}
rm(i)

colMeans(annual_sy)

# subset calibration data for plotting
hiv_calibrate = annual_hiv[,c(2:10)]
sy_calibrate = annual_sy[,c(2:10)]

# initialize necessary estimates for plotting
hiv_sur = c(305,286,304,328,288,313,271,278,215)
hiv_sim = colMeans(hiv_calibrate)
hiv_ll = numeric(9)
hiv_ul = numeric(9)
sy_sur = c(218,217,268,248,252,416,568,659,626)
sy_sim = colMeans(sy_calibrate)
sy_ll = numeric(9)
sy_ul = numeric(9)
for (i in 1:length(hiv_ll)) {
  hiv_ll[i] = quantile(hiv_calibrate[,i], 0.025)
  hiv_ul[i] = quantile(hiv_calibrate[,i], 0.975)
  sy_ll[i] = quantile(sy_calibrate[,i], 0.025)
  sy_ul[i] = quantile(sy_calibrate[,i], 0.975)
}

# calibration dataframe 
dta_calibrate = data.frame(years=rep(2010:2018), group=rep(c("Simulated","Surveillance"), each=9, length=36), sti=rep(c("HIV","Syphilis"), each=18), estimates=c(hiv_sim, hiv_sur, sy_sim, sy_sur), ll = c(hiv_ll, rep(NA, 9), sy_ll, rep(NA, 9)), ul = c(hiv_ul, rep(NA, 9), sy_ul, rep(NA, 9)))

# calibration plot 
c_plot = ggplot(dta_calibrate, aes(years, estimates, linetype = group, shape = group)) +
  geom_line() +
  scale_linetype_manual(values=c("longdash", "solid")) +
  geom_point(size = 2) + 
  scale_shape_manual(values = c(1,16)) +
  geom_pointrange(data = dta_calibrate, aes(ymin=ll, ymax=ul), width = 0.1) +
  ylab("Reported Cases") +
  xlab("Surveillance Years") +
  coord_cartesian(ylim = c(0, 875)) +
  scale_x_continuous(breaks = 0:2100) +
  facet_wrap(~sti, scales = "free") +
  theme_bw() +
  theme(text = element_text(color = "#22211d", size = 12),legend.text=element_text(size=12), plot.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 12), legend.title = element_blank(),legend.position = "top") 

ggsave("calibration_plot.jpeg", plot=c_plot, width = 9.5, height = 5.5, units = "in", dpi = 1200)

# clean up environment 
rm(c_plot, annual_hiv, annual_sy, hiv_calibrate, sy_calibrate, dta_calibrate, i, hiv_ll, hiv_ul, sy_ll, sy_ul, hiv_sim, hiv_sur, sy_sim, sy_sur)


### STEP6. Cumulative incidence Estimation ###
# Separate abm_sims and long_sims into respective doxy uptake scenarios 
abm_sims_20 = c(abm_sims_50[rep(c(TRUE, FALSE), c(6, 20))])
abm_sims_40 = c(abm_sims_50[rep(c(TRUE, FALSE, TRUE, FALSE), c(1, 5, 5, 15))])
abm_sims_60 = c(abm_sims_50[rep(c(TRUE, FALSE, TRUE, FALSE), c(1, 10, 5, 10))])
abm_sims_80 = c(abm_sims_50[rep(c(TRUE, FALSE, TRUE, FALSE), c(1, 15, 5, 5))])
abm_sims_100 = c(abm_sims_50[rep(c(TRUE, FALSE, TRUE), c(1, 20, 5))])

long_sims_20 = c(long_sims_50[rep(c(TRUE, FALSE), c(6, 20))])
long_sims_40 = c(long_sims_50[rep(c(TRUE, FALSE, TRUE, FALSE), c(1, 5, 5, 15))])
long_sims_60 = c(long_sims_50[rep(c(TRUE, FALSE, TRUE, FALSE), c(1, 10, 5, 10))])
long_sims_80 = c(long_sims_50[rep(c(TRUE, FALSE, TRUE, FALSE), c(1, 15, 5, 5))])
long_sims_100 = c(long_sims_50[rep(c(TRUE, FALSE, TRUE), c(1, 20, 5))])


### Cumulative incidence 

uptake = c(20,40,60,80,100)
adhere = c(0,20,40,60,80,100)
day = 365.24

analwide = data.frame("Uptake"=NA,"Adhere"=NA,"yr1"=NA,"yr2"=NA,"yr3"=NA,"yr4"=NA,"yr5"=NA,"yr6"=NA,"yr7"=NA,"yr8"=NA,"yr9"=NA,"yr10"=NA,stringsAsFactors=F)
oralwide = data.frame("Uptake"=NA,"Adhere"=NA,"yr1"=NA,"yr2"=NA,"yr3"=NA,"yr4"=NA,"yr5"=NA,"yr6"=NA,"yr7"=NA,"yr8"=NA,"yr9"=NA,"yr10"=NA,stringsAsFactors=F)
dtawide = data.frame("Uptake"=NA,"Adhere"=NA,"yr1"=NA,"yr2"=NA,"yr3"=NA,"yr4"=NA,"yr5"=NA,"yr6"=NA,"yr7"=NA,"yr8"=NA,"yr9"=NA,"yr10"=NA,stringsAsFactors=F)

for (i in 1:length(uptake))
{
  for (j in 1:length(adhere))
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
    # anal + oral infections 
    stat_yr1_new = NA           #1 yr new infections
    stat_yr2_new = NA           #2 yr new infections
    stat_yr3_new = NA           #3 yr new infections
    stat_yr4_new = NA           #4 yr new infections
    stat_yr5_new = NA           #5 yr new infections
    stat_yr6_new = NA           #6 yr new infections
    stat_yr7_new = NA           #7 yr new infections
    stat_yr8_new = NA           #8 yr new infections
    stat_yr9_new = NA           #9 yr new infections
    stat_yr10_new = NA          #10 yr new infections
    
    for (k in 0:(nsims-1))
    {
      current = k*length(adhere)
      
      if (uptake[i]==20) {
        # 20% Doxy Uptake
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*10)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new_anal = c(stat_yr2_new_anal, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*11)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new_anal = c(stat_yr3_new_anal, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*12)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new_anal = c(stat_yr4_new_anal, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*13)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*14)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new_anal = c(stat_yr6_new_anal, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*15)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new_anal = c(stat_yr7_new_anal, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*16)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new_anal = c(stat_yr8_new_anal, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*17)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new_anal = c(stat_yr9_new_anal, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*18)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*19)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*20))))
        
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*10)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new_oral = c(stat_yr2_new_oral, (sum(long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*11)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new_oral = c(stat_yr3_new_oral, (sum(long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*12)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new_oral = c(stat_yr4_new_oral, (sum(long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*13)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*14)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new_oral = c(stat_yr6_new_oral, (sum(long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*15)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new_oral = c(stat_yr7_new_oral, (sum(long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*16)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new_oral = c(stat_yr8_new_oral, (sum(long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*17)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new_oral = c(stat_yr9_new_oral, (sum(long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*18)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*19)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*20))))
        
        stat_yr1_new = c(stat_yr1_new, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*10)], long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*10)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new = c(stat_yr2_new, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*11)], long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*11)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new = c(stat_yr3_new, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*12)], long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*12)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new = c(stat_yr4_new, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*13)], long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*13)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new = c(stat_yr5_new, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*14)], long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*14)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new = c(stat_yr6_new, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*15)], long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*15)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new = c(stat_yr7_new, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*16)], long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*16)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new = c(stat_yr8_new, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*17)], long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*17)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new = c(stat_yr9_new, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*18)], long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*18)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new = c(stat_yr10_new, (sum(long_sims_20[[j+current]]$SY_ANAL[(day*9+1):(day*19)], long_sims_20[[j+current]]$SY_ORAL[(day*9+1):(day*19)]) / (nagents_start + (((nrow(abm_sims_20[[j+current]]) - nagents_start)/(length_sim/day))*20))))
      } else if (uptake[i]==40) {
        # 40% Doxy Uptake 
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*10)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new_anal = c(stat_yr2_new_anal, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*11)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new_anal = c(stat_yr3_new_anal, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*12)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new_anal = c(stat_yr4_new_anal, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*13)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*14)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new_anal = c(stat_yr6_new_anal, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*15)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new_anal = c(stat_yr7_new_anal, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*16)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new_anal = c(stat_yr8_new_anal, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*17)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new_anal = c(stat_yr9_new_anal, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*18)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*19)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*20))))
        
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*10)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new_oral = c(stat_yr2_new_oral, (sum(long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*11)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new_oral = c(stat_yr3_new_oral, (sum(long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*12)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new_oral = c(stat_yr4_new_oral, (sum(long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*13)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*14)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new_oral = c(stat_yr6_new_oral, (sum(long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*15)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new_oral = c(stat_yr7_new_oral, (sum(long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*16)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new_oral = c(stat_yr8_new_oral, (sum(long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*17)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new_oral = c(stat_yr9_new_oral, (sum(long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*18)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*19)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*20))))
        
        stat_yr1_new = c(stat_yr1_new, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*10)], long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*10)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new = c(stat_yr2_new, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*11)], long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*11)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new = c(stat_yr3_new, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*12)], long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*12)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new = c(stat_yr4_new, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*13)], long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*13)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new = c(stat_yr5_new, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*14)], long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*14)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new = c(stat_yr6_new, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*15)], long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*15)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new = c(stat_yr7_new, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*16)], long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*16)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new = c(stat_yr8_new, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*17)], long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*17)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new = c(stat_yr9_new, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*18)], long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*18)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new = c(stat_yr10_new, (sum(long_sims_40[[j+current]]$SY_ANAL[(day*9+1):(day*19)], long_sims_40[[j+current]]$SY_ORAL[(day*9+1):(day*19)]) / (nagents_start + (((nrow(abm_sims_40[[j+current]]) - nagents_start)/(length_sim/day))*20))))
      } else if (uptake[i]==60) {
        # 60% Doxy Uptake
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*10)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new_anal = c(stat_yr2_new_anal, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*11)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new_anal = c(stat_yr3_new_anal, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*12)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new_anal = c(stat_yr4_new_anal, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*13)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*14)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new_anal = c(stat_yr6_new_anal, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*15)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new_anal = c(stat_yr7_new_anal, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*16)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new_anal = c(stat_yr8_new_anal, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*17)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new_anal = c(stat_yr9_new_anal, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*18)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*19)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*20))))
        
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*10)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new_oral = c(stat_yr2_new_oral, (sum(long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*11)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new_oral = c(stat_yr3_new_oral, (sum(long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*12)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new_oral = c(stat_yr4_new_oral, (sum(long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*13)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*14)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new_oral = c(stat_yr6_new_oral, (sum(long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*15)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new_oral = c(stat_yr7_new_oral, (sum(long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*16)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new_oral = c(stat_yr8_new_oral, (sum(long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*17)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new_oral = c(stat_yr9_new_oral, (sum(long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*18)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*19)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*20))))   
        
        stat_yr1_new = c(stat_yr1_new, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*10)], long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*10)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new = c(stat_yr2_new, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*11)], long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*11)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new = c(stat_yr3_new, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*12)], long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*12)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new = c(stat_yr4_new, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*13)], long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*13)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new = c(stat_yr5_new, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*14)], long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*14)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new = c(stat_yr6_new, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*15)], long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*15)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new = c(stat_yr7_new, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*16)], long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*16)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new = c(stat_yr8_new, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*17)], long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*17)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new = c(stat_yr9_new, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*18)], long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*18)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new = c(stat_yr10_new, (sum(long_sims_60[[j+current]]$SY_ANAL[(day*9+1):(day*19)], long_sims_60[[j+current]]$SY_ORAL[(day*9+1):(day*19)]) / (nagents_start + (((nrow(abm_sims_60[[j+current]]) - nagents_start)/(length_sim/day))*20))))
      } else if (uptake[i]==80) {
        # 80% Doxy Uptake 
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*10)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new_anal = c(stat_yr2_new_anal, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*11)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new_anal = c(stat_yr3_new_anal, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*12)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new_anal = c(stat_yr4_new_anal, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*13)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*14)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new_anal = c(stat_yr6_new_anal, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*15)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new_anal = c(stat_yr7_new_anal, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*16)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new_anal = c(stat_yr8_new_anal, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*17)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new_anal = c(stat_yr9_new_anal, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*18)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*19)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*20))))
        
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*10)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new_oral = c(stat_yr2_new_oral, (sum(long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*11)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new_oral = c(stat_yr3_new_oral, (sum(long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*12)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new_oral = c(stat_yr4_new_oral, (sum(long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*13)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*14)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new_oral = c(stat_yr6_new_oral, (sum(long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*15)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new_oral = c(stat_yr7_new_oral, (sum(long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*16)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new_oral = c(stat_yr8_new_oral, (sum(long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*17)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new_oral = c(stat_yr9_new_oral, (sum(long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*18)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*19)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*20))))
        
        stat_yr1_new = c(stat_yr1_new, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*10)], long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*10)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new = c(stat_yr2_new, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*11)], long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*11)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new = c(stat_yr3_new, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*12)], long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*12)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new = c(stat_yr4_new, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*13)], long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*13)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new = c(stat_yr5_new, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*14)], long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*14)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new = c(stat_yr6_new, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*15)], long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*15)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new = c(stat_yr7_new, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*16)], long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*16)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new = c(stat_yr8_new, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*17)], long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*17)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new = c(stat_yr9_new, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*18)], long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*18)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new = c(stat_yr10_new, (sum(long_sims_80[[j+current]]$SY_ANAL[(day*9+1):(day*19)], long_sims_80[[j+current]]$SY_ORAL[(day*9+1):(day*19)]) / (nagents_start + (((nrow(abm_sims_80[[j+current]]) - nagents_start)/(length_sim/day))*20))))
      } else if (uptake[i]==100) {
        # 100% Doxy Uptake
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*10)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new_anal = c(stat_yr2_new_anal, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*11)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new_anal = c(stat_yr3_new_anal, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*12)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new_anal = c(stat_yr4_new_anal, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*13)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*14)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new_anal = c(stat_yr6_new_anal, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*15)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new_anal = c(stat_yr7_new_anal, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*16)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new_anal = c(stat_yr8_new_anal, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*17)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new_anal = c(stat_yr9_new_anal, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*18)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*19)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*20))))
        
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*10)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new_oral = c(stat_yr2_new_oral, (sum(long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*11)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new_oral = c(stat_yr3_new_oral, (sum(long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*12)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new_oral = c(stat_yr4_new_oral, (sum(long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*13)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*14)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new_oral = c(stat_yr6_new_oral, (sum(long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*15)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new_oral = c(stat_yr7_new_oral, (sum(long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*16)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new_oral = c(stat_yr8_new_oral, (sum(long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*17)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new_oral = c(stat_yr9_new_oral, (sum(long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*18)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*19)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*20))))
        
        stat_yr1_new = c(stat_yr1_new, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*10)], long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*10)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*11))))
        stat_yr2_new = c(stat_yr2_new, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*11)], long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*11)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*12))))
        stat_yr3_new = c(stat_yr3_new, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*12)], long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*12)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*13))))
        stat_yr4_new = c(stat_yr4_new, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*13)], long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*13)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*14))))
        stat_yr5_new = c(stat_yr5_new, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*14)], long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*14)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*15))))
        stat_yr6_new = c(stat_yr6_new, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*15)], long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*15)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*16))))
        stat_yr7_new = c(stat_yr7_new, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*16)], long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*16)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*17))))
        stat_yr8_new = c(stat_yr8_new, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*17)], long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*17)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*18))))
        stat_yr9_new = c(stat_yr9_new, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*18)], long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*18)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*19))))
        stat_yr10_new = c(stat_yr10_new, (sum(long_sims_100[[j+current]]$SY_ANAL[(day*9+1):(day*19)], long_sims_100[[j+current]]$SY_ORAL[(day*9+1):(day*19)]) / (nagents_start + (((nrow(abm_sims_100[[j+current]]) - nagents_start)/(length_sim/day))*20))))
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
    # anal + oral 
    stat_yr1_new = stat_yr1_new[-1]
    stat_yr2_new = stat_yr2_new[-1]
    stat_yr3_new = stat_yr3_new[-1]
    stat_yr4_new = stat_yr4_new[-1]
    stat_yr5_new = stat_yr5_new[-1]
    stat_yr6_new = stat_yr6_new[-1]
    stat_yr7_new = stat_yr7_new[-1]
    stat_yr8_new = stat_yr8_new[-1]
    stat_yr9_new = stat_yr9_new[-1]
    stat_yr10_new = stat_yr10_new[-1]
    
    analwide = rbind(analwide, data.frame("Uptake"=uptake[i],"Adhere"=adhere[j],"yr1"=stat_yr1_new_anal,"yr2"=stat_yr2_new_anal,"yr3"=stat_yr3_new_anal,"yr4"=stat_yr4_new_anal,"yr5"=stat_yr5_new_anal,"yr6"=stat_yr6_new_anal,"yr7"=stat_yr7_new_anal,"yr8"=stat_yr8_new_anal,"yr9"=stat_yr9_new_anal,"yr10"=stat_yr10_new_anal,stringsAsFactors=F))
    oralwide = rbind(oralwide, data.frame("Uptake"=uptake[i],"Adhere"=adhere[j],"yr1"=stat_yr1_new_oral,"yr2"=stat_yr2_new_oral,"yr3"=stat_yr3_new_oral,"yr4"=stat_yr4_new_oral,"yr5"=stat_yr5_new_oral,"yr6"=stat_yr6_new_oral,"yr7"=stat_yr7_new_oral,"yr8"=stat_yr8_new_oral,"yr9"=stat_yr9_new_oral,"yr10"=stat_yr10_new_oral,stringsAsFactors=F))
    dtawide = rbind(dtawide, data.frame("Uptake"=uptake[i],"Adhere"=adhere[j],"yr1"=stat_yr1_new,"yr2"=stat_yr2_new,"yr3"=stat_yr3_new,"yr4"=stat_yr4_new,"yr5"=stat_yr5_new,"yr6"=stat_yr6_new,"yr7"=stat_yr7_new,"yr8"=stat_yr8_new,"yr9"=stat_yr9_new,"yr10"=stat_yr10_new,stringsAsFactors=F))
  }
}

analwide = analwide[-1,]
oralwide = oralwide[-1,]
dtawide = dtawide[-1,]
rm(uptake, adhere, day, current, i, j, k, stat_yr1_new_anal, stat_yr2_new_anal, stat_yr3_new_anal, stat_yr4_new_anal, stat_yr5_new_anal, stat_yr6_new_anal, stat_yr7_new_anal, stat_yr8_new_anal, stat_yr9_new_anal, stat_yr10_new_anal, stat_yr1_new_oral, stat_yr2_new_oral, stat_yr3_new_oral, stat_yr4_new_oral, stat_yr5_new_oral, stat_yr6_new_oral, stat_yr7_new_oral, stat_yr8_new_oral, stat_yr9_new_oral, stat_yr10_new_oral, stat_yr1_new, stat_yr2_new, stat_yr3_new, stat_yr4_new, stat_yr5_new, stat_yr6_new, stat_yr7_new, stat_yr8_new, stat_yr9_new, stat_yr10_new, long_sims_20, long_sims_40, long_sims_60, long_sims_80, long_sims_100)

# convert to long form
analLong = reshape(analwide, timevar = "year", times = 2019:2028, v.names = "incidence", varying = list(3:12), direction = "long")
oralLong = reshape(oralwide, timevar = "year", times = 2019:2028, v.names = "incidence", varying = list(3:12), direction = "long")
dtaLong = reshape(dtawide, timevar = "year", times = 2019:2028, v.names = "incidence", varying = list(3:12), direction = "long")

analLong$year1 = analLong$year - 2018
oralLong$year1 = oralLong$year - 2018
dtaLong$year1 = dtaLong$year - 2018

# incidence per 1000
analLong$ir = analLong$incidence * 1000
oralLong$ir = oralLong$incidence * 1000
dtaLong$ir = dtaLong$incidence * 1000

# rename uptake categories for plots
analLong$Uptake = ifelse(analLong$Uptake==20, "20% Uptake", ifelse(analLong$Uptake==40, "40% Uptake", ifelse(analLong$Uptake==60, "60% Uptake", ifelse(analLong$Uptake==80, "80% Uptake", "100% Uptake"))))
oralLong$Uptake = ifelse(oralLong$Uptake==20, "20% Uptake", ifelse(oralLong$Uptake==40, "40% Uptake", ifelse(oralLong$Uptake==60, "60% Uptake", ifelse(oralLong$Uptake==80, "80% Uptake", "100% Uptake"))))
dtaLong$Uptake = ifelse(dtaLong$Uptake==20, "20% Uptake", ifelse(dtaLong$Uptake==40, "40% Uptake", ifelse(dtaLong$Uptake==60, "60% Uptake", ifelse(dtaLong$Uptake==80, "80% Uptake", "100% Uptake"))))


analLong$Uptake = factor(analLong$Uptake, levels=c("20% Uptake","40% Uptake","60% Uptake","80% Uptake","100% Uptake"))
oralLong$Uptake = factor(oralLong$Uptake, levels=c("20% Uptake","40% Uptake","60% Uptake","80% Uptake","100% Uptake"))
dtaLong$Uptake = factor(dtaLong$Uptake, levels=c("20% Uptake","40% Uptake","60% Uptake","80% Uptake","100% Uptake"))

rm(oralwide, analwide, dtawide)

# FIG 2. Incidence plots
p1 = ggplot(analLong, aes(x=as.factor(year1), y=ir, col = as.factor(Adhere), group = as.factor(Adhere))) +
  geom_pointrange(mapping = aes(x = as.factor(year1), y = ir, group = as.factor(Adhere)),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.025)},
                  fun.max = function(z) {quantile(z,0.975)},
                  fun = mean) +
  stat_summary(fun=mean, geom="line", linetype = 1) +
  labs(x="Years", y="Cumulative Incidence of Anal Infections per 1,000") + 
  scale_color_grey(name = "Adherence",labels = c("None","20%","40%","60%","80%","100%"), start = 0.0, end = 0.7) +
  facet_wrap(~Uptake) +
  theme_bw() +
  theme(text = element_text(color = "#22211d", size = 12),legend.text=element_text(size=12), plot.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 12), legend.title = element_text(size=12)) 

p2 = ggplot(oralLong, aes(x=as.factor(year1), y=ir, col = as.factor(Adhere), group = as.factor(Adhere))) +
  geom_pointrange(mapping = aes(x = as.factor(year1), y = ir, group = as.factor(Adhere)),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.025)},
                  fun.max = function(z) {quantile(z,0.975)},
                  fun = mean) +
  stat_summary(fun=mean, geom="line", linetype = 1) +
  labs(x="Years", y="Cumulative Incidence of Oral Infections per 1,000") + 
  scale_color_grey(name = "Adherence",labels = c("None","20%","40%","60%","80%","100%"), start = 0.0, end = 0.7) +
  facet_wrap(~Uptake) +
  theme_bw() +
  theme(text = element_text(color = "#22211d", size = 12),legend.text=element_text(size=12), plot.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 10), legend.title = element_text(size=12)) 

p3 = ggplot(dtaLong, aes(x=as.factor(year1), y=ir, col = as.factor(Adhere), group = as.factor(Adhere))) +
  geom_pointrange(mapping = aes(x = as.factor(year1), y = ir, group = as.factor(Adhere)),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.025)},
                  fun.max = function(z) {quantile(z,0.975)},
                  fun = mean) +
  stat_summary(fun=mean, geom="line", linetype = 1) +
  labs(x="Years", y="Cumulative Incidence of Oral Infections per 1,000") + 
  scale_color_grey(name = "Adherence",labels = c("None","20%","40%","60%","80%","100%"), start = 0.0, end = 0.7) +
  facet_wrap(~Uptake) +
  theme_bw() +
  theme(text = element_text(color = "#22211d", size = 12),legend.text=element_text(size=12), plot.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 12), legend.title = element_text(size=12)) 

ggsave("ci_anal.jpeg", plot=p1, width = 10, height = 6.7, units = "in", dpi = 1200)
ggsave("ci_oral.jpeg", plot=p2, width = 10, height = 6.7, units = "in", dpi = 1200)
ggsave("ci_overall.jpeg", plot=p3, width = 10, height = 6.7, units = "in", dpi = 1200)

# regression 

dtaLong$Uptake = as.factor(dtaLong$Uptake)
dtaLong$Adhere = as.factor(dtaLong$Adhere)

# risk difference in anal infections for each level of doxy pep adherence at year 1, 5, and 10 under the 20% uptake scenario
m1 = lm(ir ~ as.factor(Adhere), data = dtaLong[dtaLong$Uptake=="20% Uptake" & dtaLong$year1==1,])
m2 = lm(ir ~ as.factor(Adhere), data = dtaLong[dtaLong$Uptake=="20% Uptake" & dtaLong$year1==5,])
m3 = lm(ir ~ as.factor(Adhere), data = dtaLong[dtaLong$Uptake=="20% Uptake" & dtaLong$year1==10,])

# calculate % difference
abs(coef(m3)[2]/coef(m3)[1])
abs(coef(m3)[3]/coef(m3)[1])
abs(coef(m3)[4]/coef(m3)[1])
abs(coef(m3)[5]/coef(m3)[1])
abs(coef(m3)[6]/coef(m3)[1])

summary(m3)

# risk difference in anal infections for each level of doxy pep uptake at year 1, 5, and 10 under the 100% uptake scenario 
m1 = lm(ir ~ as.factor(Adhere), data = dtaLong[dtaLong$Uptake=="100% Uptake" & dtaLong$year1==1,])
m2 = lm(ir ~ as.factor(Adhere), data = dtaLong[dtaLong$Uptake=="100% Uptake" & dtaLong$year1==5,])
m3 = lm(ir ~ as.factor(Adhere), data = dtaLong[dtaLong$Uptake=="100% Uptake" & dtaLong$year1==10,])

# calculate % difference
abs(coef(m3)[2]/coef(m3)[1])
abs(coef(m3)[3]/coef(m3)[1])
abs(coef(m3)[4]/coef(m3)[1])
abs(coef(m3)[5]/coef(m3)[1])
abs(coef(m3)[6]/coef(m3)[1])

summary(m3)

rm(m1,m2,m3,oralLong,analLong,dtaLong,p1,p2,p3)

### STEP7: Estimate PIA for Condoms and Doxy ###

uptake = c(20,40,60,80,100)
adhere = c(0,20,40,60,80,100)

results_infections_anal = data.frame("Uptake"=NA,"Adhere"=NA,"Total_mean"=NA,"Total_SE"=NA,"Prev_mean"=NA,"Prev_SE"=NA,"New_mean"=NA,"New_SE"=NA,"Incidence_mean"=NA,"Incidence_SE"=NA,"Prevent_mean"=NA,"Prevent_SE"=NA,"Doxy_prevent_mean"=NA,"Doxy_prevent_SE"=NA,"Sero_prevent_mean"=NA,"Sero_prevent_SE"=NA,"Cond_prevent_mean"=NA,"Cond_prevent_SE"=NA,"Total_prevent"=NA,"Doxy_prevent_pct"=NA,"Doxy_prevent_pct_ll"=NA,"Doxy_prevent_pct_ul"=NA,"Sero_prevent_pct"=NA,"Sero_prevent_pct_ll"=NA,"Sero_prevent_pct_ul"=NA,"Cond_prevent_pct"=NA,"Cond_prevent_pct_ll"=NA,"Cond_prevent_pct_ul"=NA,stringsAsFactors=F)
results_infections_oral = data.frame("Uptake"=NA,"Adhere"=NA,"Total_mean"=NA,"Total_SE"=NA,"Prev_mean"=NA,"Prev_SE"=NA,"New_mean"=NA,"New_SE"=NA,"Incidence_mean"=NA,"Incidence_SE"=NA,"Prevent_mean"=NA,"Prevent_SE"=NA,"Doxy_prevent_mean"=NA,"Doxy_prevent_SE"=NA,"Cond_prevent_mean"=NA,"Cond_prevent_SE"=NA,"Total_prevent"=NA,"Doxy_prevent_pct"=NA,"Doxy_prevent_pct_ll"=NA,"Doxy_prevent_pct_ul"=NA,"Cond_prevent_pct"=NA,"Cond_prevent_pct_ll"=NA,"Cond_prevent_pct_ul"=NA,stringsAsFactors=F)
results_infections = data.frame("Uptake"=NA,"Adhere"=NA,"Total_mean"=NA,"Total_SE"=NA,"Prev_mean"=NA,"Prev_SE"=NA,"New_mean"=NA,"New_SE"=NA,"Incidence_mean"=NA,"Incidence_SE"=NA,"Prevent_mean"=NA,"Prevent_SE"=NA,"Doxy_prevent_mean"=NA,"Doxy_prevent_SE"=NA,"Sero_prevent_mean"=NA,"Sero_prevent_SE"=NA,"Cond_prevent_mean"=NA,"Cond_prevent_SE"=NA,"Total_prevent"=NA,"Doxy_prevent_pct"=NA,"Doxy_prevent_pct_ll"=NA,"Doxy_prevent_pct_ul"=NA,"Sero_prevent_pct"=NA,"Sero_prevent_pct_ll"=NA,"Sero_prevent_pct_ul"=NA,"Cond_prevent_pct"=NA,"Cond_prevent_pct_ll"=NA,"Cond_prevent_pct_ul"=NA,stringsAsFactors=F)

for (i in 1:length(uptake))
{
  for (j in 1:length(adhere))
  {
    stat_prev_anal = NA               #point prevalance at post simulation
    stat_incidence_anal = NA          #cumulative incidence (primary/secondary + early latent)
    stat_new_anal = NA                #primary/secondary + early latent
    stat_total_anal = NA              #infections at post simulation
    stat_prevent_anal = NA            #total preventions 
    stat_prevent_doxy_anal = NA       #total doxy pep preventions 
    stat_prevent_sero_anal = NA       #total sero preventions 
    stat_prevent_cond_anal = NA       #total condom preventions
    stat_prevent_doxy_anal_pct = NA   #percentage of anal infections prevented due to doxy
    stat_prevent_sero_anal_pct = NA   #percentage of anal infections prevented due to doxy
    stat_prevent_cond_anal_pct = NA   #percentage of anal infections prevented due to condoms 
    stat_prev_oral = NA               #point prevalance at post simulation
    stat_incidence_oral = NA          #cumulative incidence (primary/secondary + early latent)
    stat_new_oral = NA                #primary/secondary + early latent
    stat_total_oral = NA              #infections at post simulation
    stat_prevent_oral = NA            #total preventions 
    stat_prevent_doxy_oral = NA       #total doxy pep preventions 
    stat_prevent_cond_oral = NA       #total condom preventions 
    stat_prevent_doxy_oral_pct = NA   #percentage of oral infections prevented due to doxy
    stat_prevent_cond_oral_pct = NA   #percentage of oral infections prevented due to condoms 
    stat_prev = NA                    #point prevalance at post simulation
    stat_incidence = NA               #cumulative incidence (primary/secondary + early latent)
    stat_new = NA                     #primary/secondary + early latent
    stat_total = NA                   #infections at post simulation
    stat_prevent = NA                 #total preventions
    stat_prevent_doxy = NA            #total doxy pep preventions
    stat_prevent_sero = NA            #total sero preventions 
    stat_prevent_cond = NA            #total condom preventions
    stat_prevent_doxy_pct = NA        #percentage of anal + oral infections prevented due to doxy
    stat_prevent_sero_pct = NA        #percentage of anal infections prevented due to seroadaption (no oral) 
    stat_prevent_cond_pct = NA        #percentage of oral infections prevented due to condoms
    
    for (k in 0:(nsims-1))
    {
      current = k*length(adhere)
      if (uptake[i]==20) {
        # 20% uptake
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_20[[j+current]]$SY_ANAL)))
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_20[[j+current]]$SY_ANAL))/nrow(abm_sims_20[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_20[[j+current]]$SY_ANAL) - sum(abm_sims_20[[j+current]]$ORIGINAL_STATUS_SY_ANAL)))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_20[[j+current]]$SY_ANAL) - sum(abm_sims_20[[j+current]]$ORIGINAL_STATUS_SY_ANAL))/nrow(abm_sims_20[[j+current]]))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_20[[j+current]]$SERO_PREVENT, abm_sims_20[[j+current]]$COND_PREVENT_ANAL, na.rm = T)))
        stat_prevent_doxy_anal = c(stat_prevent_doxy_anal, (sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ANAL, na.rm = T)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_20[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_20[[j+current]]$COND_PREVENT_ANAL)))
        stat_prevent_doxy_anal_pct = c(stat_prevent_doxy_anal_pct, ((sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ANAL, na.rm = T))/(sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_20[[j+current]]$SERO_PREVENT, abm_sims_20[[j+current]]$COND_PREVENT_ANAL, na.rm = T))))
        stat_prevent_sero_anal_pct = c(stat_prevent_sero_anal_pct, ((sum(abm_sims_20[[j+current]]$SERO_PREVENT, na.rm = T))/(sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_20[[j+current]]$SERO_PREVENT, abm_sims_20[[j+current]]$COND_PREVENT_ANAL, na.rm = T))))
        stat_prevent_cond_anal_pct = c(stat_prevent_cond_anal_pct, ((sum(abm_sims_20[[j+current]]$COND_PREVENT_ANAL))/(sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_20[[j+current]]$SERO_PREVENT, abm_sims_20[[j+current]]$COND_PREVENT_ANAL, na.rm = T))))
        
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_20[[j+current]]$SY_ORAL)))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_20[[j+current]]$SY_ORAL))/nrow(abm_sims_20[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_20[[j+current]]$SY_ORAL) - sum(abm_sims_20[[j+current]]$ORIGINAL_STATUS_SY_ORAL)))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_20[[j+current]]$SY_ORAL) - sum(abm_sims_20[[j+current]]$ORIGINAL_STATUS_SY_ORAL))/nrow(abm_sims_20[[j+current]]))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_20[[j+current]]$COND_PREVENT_ORAL, na.rm = T)))
        stat_prevent_doxy_oral = c(stat_prevent_doxy_oral, (sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_20[[j+current]]$COND_PREVENT_ORAL)))
        stat_prevent_doxy_oral_pct = c(stat_prevent_doxy_oral_pct, ((sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T))/(sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_20[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        stat_prevent_cond_oral_pct = c(stat_prevent_cond_oral_pct, ((sum(abm_sims_20[[j+current]]$COND_PREVENT_ORAL))/(sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_20[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        
        stat_total = c(stat_total, (sum(abm_sims_20[[j+current]]$SY_ANAL, abm_sims_20[[j+current]]$SY_ORAL)))
        stat_prev = c(stat_prev, (sum(abm_sims_20[[j+current]]$SY_ANAL, abm_sims_20[[j+current]]$SY_ORAL))/nrow(abm_sims_20[[j+current]]))
        stat_new = c(stat_new, (sum(abm_sims_20[[j+current]]$SY_ANAL, abm_sims_20[[j+current]]$SY_ORAL) - sum(abm_sims_20[[j+current]]$ORIGINAL_STATUS_SY_ANAL, abm_sims_20[[j+current]]$ORIGINAL_STATUS_SY_ORAL)))
        stat_incidence = c(stat_incidence, (sum(abm_sims_20[[j+current]]$SY_ANAL, abm_sims_20[[j+current]]$SY_ORAL) - sum(abm_sims_20[[j+current]]$ORIGINAL_STATUS_SY_ANAL, abm_sims_20[[j+current]]$ORIGINAL_STATUS_SY_ORAL))/nrow(abm_sims_20[[j+current]]))
        stat_prevent = c(stat_prevent, (sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_20[[j+current]]$SERO_PREVENT, abm_sims_20[[j+current]]$COND_PREVENT_ANAL, abm_sims_20[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_20[[j+current]]$COND_PREVENT_ORAL, na.rm = T)))
        stat_prevent_doxy = c(stat_prevent_doxy, (sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_20[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T)))
        stat_prevent_sero = c(stat_prevent_sero, (sum(abm_sims_20[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond = c(stat_prevent_cond, (sum(abm_sims_20[[j+current]]$COND_PREVENT_ANAL, abm_sims_20[[j+current]]$COND_PREVENT_ORAL)))
        stat_prevent_doxy_pct = c(stat_prevent_doxy_pct, ((sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_20[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T))/(sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_20[[j+current]]$SERO_PREVENT, abm_sims_20[[j+current]]$COND_PREVENT_ANAL, abm_sims_20[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_20[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        stat_prevent_sero_pct = c(stat_prevent_sero_pct, ((sum(abm_sims_20[[j+current]]$SERO_PREVENT))/(sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_20[[j+current]]$SERO_PREVENT, abm_sims_20[[j+current]]$COND_PREVENT_ANAL, abm_sims_20[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_20[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        stat_prevent_cond_pct = c(stat_prevent_cond_pct, ((sum(abm_sims_20[[j+current]]$COND_PREVENT_ANAL, abm_sims_20[[j+current]]$COND_PREVENT_ORAL))/(sum(abm_sims_20[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_20[[j+current]]$SERO_PREVENT, abm_sims_20[[j+current]]$COND_PREVENT_ANAL, abm_sims_20[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_20[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
      } else if (uptake[i]==40) {
        # 40% Uptake
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_40[[j+current]]$SY_ANAL)))
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_40[[j+current]]$SY_ANAL))/nrow(abm_sims_40[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_40[[j+current]]$SY_ANAL) - sum(abm_sims_40[[j+current]]$ORIGINAL_STATUS_SY_ANAL)))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_40[[j+current]]$SY_ANAL) - sum(abm_sims_40[[j+current]]$ORIGINAL_STATUS_SY_ANAL))/nrow(abm_sims_40[[j+current]]))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_40[[j+current]]$SERO_PREVENT, abm_sims_40[[j+current]]$COND_PREVENT_ANAL, na.rm = T)))
        stat_prevent_doxy_anal = c(stat_prevent_doxy_anal, (sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ANAL, na.rm = T)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_40[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_40[[j+current]]$COND_PREVENT_ANAL)))
        stat_prevent_doxy_anal_pct = c(stat_prevent_doxy_anal_pct, ((sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ANAL, na.rm = T))/(sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_40[[j+current]]$SERO_PREVENT, abm_sims_40[[j+current]]$COND_PREVENT_ANAL, na.rm = T))))
        stat_prevent_sero_anal_pct = c(stat_prevent_sero_anal_pct, ((sum(abm_sims_40[[j+current]]$SERO_PREVENT))/(sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_40[[j+current]]$SERO_PREVENT, abm_sims_40[[j+current]]$COND_PREVENT_ANAL, na.rm = T))))
        stat_prevent_cond_anal_pct = c(stat_prevent_cond_anal_pct, ((sum(abm_sims_40[[j+current]]$COND_PREVENT_ANAL))/(sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_40[[j+current]]$SERO_PREVENT, abm_sims_40[[j+current]]$COND_PREVENT_ANAL, na.rm = T))))
        
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_40[[j+current]]$SY_ORAL)))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_40[[j+current]]$SY_ORAL))/nrow(abm_sims_40[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_40[[j+current]]$SY_ORAL) - sum(abm_sims_40[[j+current]]$ORIGINAL_STATUS_SY_ORAL)))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_40[[j+current]]$SY_ORAL) - sum(abm_sims_40[[j+current]]$ORIGINAL_STATUS_SY_ORAL))/nrow(abm_sims_40[[j+current]]))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_40[[j+current]]$COND_PREVENT_ORAL, na.rm = T)))
        stat_prevent_doxy_oral = c(stat_prevent_doxy_oral, (sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_40[[j+current]]$COND_PREVENT_ORAL)))
        stat_prevent_doxy_oral_pct = c(stat_prevent_doxy_oral_pct, ((sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T))/(sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_40[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        stat_prevent_cond_oral_pct = c(stat_prevent_cond_oral_pct, ((sum(abm_sims_40[[j+current]]$COND_PREVENT_ORAL))/(sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_40[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        
        stat_total = c(stat_total, (sum(abm_sims_40[[j+current]]$SY_ANAL, abm_sims_40[[j+current]]$SY_ORAL)))
        stat_prev = c(stat_prev, (sum(abm_sims_40[[j+current]]$SY_ANAL, abm_sims_40[[j+current]]$SY_ORAL))/nrow(abm_sims_40[[j+current]]))
        stat_new = c(stat_new, (sum(abm_sims_40[[j+current]]$SY_ANAL, abm_sims_40[[j+current]]$SY_ORAL) - sum(abm_sims_40[[j+current]]$ORIGINAL_STATUS_SY_ANAL, abm_sims_40[[j+current]]$ORIGINAL_STATUS_SY_ORAL)))
        stat_incidence = c(stat_incidence, (sum(abm_sims_40[[j+current]]$SY_ANAL, abm_sims_40[[j+current]]$SY_ORAL) - sum(abm_sims_40[[j+current]]$ORIGINAL_STATUS_SY_ANAL, abm_sims_40[[j+current]]$ORIGINAL_STATUS_SY_ORAL))/nrow(abm_sims_40[[j+current]]))
        stat_prevent = c(stat_prevent, (sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_40[[j+current]]$SERO_PREVENT, abm_sims_40[[j+current]]$COND_PREVENT_ANAL, abm_sims_40[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_40[[j+current]]$COND_PREVENT_ORAL, na.rm = T)))
        stat_prevent_doxy = c(stat_prevent_doxy, (sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_40[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T)))
        stat_prevent_sero = c(stat_prevent_sero, (sum(abm_sims_40[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond = c(stat_prevent_cond, (sum(abm_sims_40[[j+current]]$COND_PREVENT_ANAL, abm_sims_40[[j+current]]$COND_PREVENT_ORAL)))
        stat_prevent_doxy_pct = c(stat_prevent_doxy_pct, ((sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_40[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T))/(sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_40[[j+current]]$SERO_PREVENT, abm_sims_40[[j+current]]$COND_PREVENT_ANAL, abm_sims_40[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_40[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        stat_prevent_sero_pct = c(stat_prevent_sero_pct, ((sum(abm_sims_40[[j+current]]$SERO_PREVENT))/(sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_40[[j+current]]$SERO_PREVENT, abm_sims_40[[j+current]]$COND_PREVENT_ANAL, abm_sims_40[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_40[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        stat_prevent_cond_pct = c(stat_prevent_cond_pct, ((sum(abm_sims_40[[j+current]]$COND_PREVENT_ANAL, abm_sims_40[[j+current]]$COND_PREVENT_ORAL))/(sum(abm_sims_40[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_40[[j+current]]$SERO_PREVENT, abm_sims_40[[j+current]]$COND_PREVENT_ANAL, abm_sims_40[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_40[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
      } else if (uptake[i]==60) {
        # 60% Uptake
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_60[[j+current]]$SY_ANAL)))
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_60[[j+current]]$SY_ANAL))/nrow(abm_sims_60[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_60[[j+current]]$SY_ANAL) - sum(abm_sims_60[[j+current]]$ORIGINAL_STATUS_SY_ANAL)))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_60[[j+current]]$SY_ANAL) - sum(abm_sims_60[[j+current]]$ORIGINAL_STATUS_SY_ANAL))/nrow(abm_sims_60[[j+current]]))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_60[[j+current]]$SERO_PREVENT, abm_sims_60[[j+current]]$COND_PREVENT_ANAL, na.rm = T)))
        stat_prevent_doxy_anal = c(stat_prevent_doxy_anal, (sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ANAL, na.rm = T)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_60[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_60[[j+current]]$COND_PREVENT_ANAL)))
        stat_prevent_doxy_anal_pct = c(stat_prevent_doxy_anal_pct, ((sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ANAL, na.rm = T))/(sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_60[[j+current]]$SERO_PREVENT, abm_sims_60[[j+current]]$COND_PREVENT_ANAL, na.rm = T))))
        stat_prevent_sero_anal_pct = c(stat_prevent_sero_anal_pct, ((sum(abm_sims_60[[j+current]]$SERO_PREVENT))/(sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_60[[j+current]]$SERO_PREVENT, abm_sims_60[[j+current]]$COND_PREVENT_ANAL, na.rm = T))))
        stat_prevent_cond_anal_pct = c(stat_prevent_cond_anal_pct, ((sum(abm_sims_60[[j+current]]$COND_PREVENT_ANAL))/(sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_60[[j+current]]$SERO_PREVENT, abm_sims_60[[j+current]]$COND_PREVENT_ANAL, na.rm = T))))
        
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_60[[j+current]]$SY_ORAL)))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_60[[j+current]]$SY_ORAL))/nrow(abm_sims_60[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_60[[j+current]]$SY_ORAL) - sum(abm_sims_60[[j+current]]$ORIGINAL_STATUS_SY_ORAL)))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_60[[j+current]]$SY_ORAL) - sum(abm_sims_60[[j+current]]$ORIGINAL_STATUS_SY_ORAL))/nrow(abm_sims_60[[j+current]]))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_60[[j+current]]$COND_PREVENT_ORAL, na.rm = T)))
        stat_prevent_doxy_oral = c(stat_prevent_doxy_oral, (sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_60[[j+current]]$COND_PREVENT_ORAL)))
        stat_prevent_doxy_oral_pct = c(stat_prevent_doxy_oral_pct, ((sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T))/(sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_60[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        stat_prevent_cond_oral_pct = c(stat_prevent_cond_oral_pct, ((sum(abm_sims_60[[j+current]]$COND_PREVENT_ORAL))/(sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_60[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        
        stat_total = c(stat_total, (sum(abm_sims_60[[j+current]]$SY_ANAL, abm_sims_60[[j+current]]$SY_ORAL)))
        stat_prev = c(stat_prev, (sum(abm_sims_60[[j+current]]$SY_ANAL, abm_sims_60[[j+current]]$SY_ORAL))/nrow(abm_sims_60[[j+current]]))
        stat_new = c(stat_new, (sum(abm_sims_60[[j+current]]$SY_ANAL, abm_sims_60[[j+current]]$SY_ORAL) - sum(abm_sims_60[[j+current]]$ORIGINAL_STATUS_SY_ANAL, abm_sims_60[[j+current]]$ORIGINAL_STATUS_SY_ORAL)))
        stat_incidence = c(stat_incidence, (sum(abm_sims_60[[j+current]]$SY_ANAL, abm_sims_60[[j+current]]$SY_ORAL) - sum(abm_sims_60[[j+current]]$ORIGINAL_STATUS_SY_ANAL, abm_sims_60[[j+current]]$ORIGINAL_STATUS_SY_ORAL))/nrow(abm_sims_60[[j+current]]))
        stat_prevent = c(stat_prevent, (sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_60[[j+current]]$SERO_PREVENT, abm_sims_60[[j+current]]$COND_PREVENT_ANAL, abm_sims_60[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_60[[j+current]]$COND_PREVENT_ORAL, na.rm = T)))
        stat_prevent_doxy = c(stat_prevent_doxy, (sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_60[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T)))
        stat_prevent_sero = c(stat_prevent_sero, (sum(abm_sims_60[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond = c(stat_prevent_cond, (sum(abm_sims_60[[j+current]]$COND_PREVENT_ANAL, abm_sims_60[[j+current]]$COND_PREVENT_ORAL)))
        stat_prevent_doxy_pct = c(stat_prevent_doxy_pct, ((sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_60[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T))/(sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_60[[j+current]]$SERO_PREVENT, abm_sims_60[[j+current]]$COND_PREVENT_ANAL, abm_sims_60[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_60[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        stat_prevent_sero_pct = c(stat_prevent_sero_pct, ((sum(abm_sims_60[[j+current]]$SERO_PREVENT))/(sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_60[[j+current]]$SERO_PREVENT, abm_sims_60[[j+current]]$COND_PREVENT_ANAL, abm_sims_60[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_60[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        stat_prevent_cond_pct = c(stat_prevent_cond_pct, ((sum(abm_sims_60[[j+current]]$COND_PREVENT_ANAL, abm_sims_60[[j+current]]$COND_PREVENT_ORAL))/(sum(abm_sims_60[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_60[[j+current]]$SERO_PREVENT, abm_sims_60[[j+current]]$COND_PREVENT_ANAL, abm_sims_60[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_60[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
      } else if (uptake[i]==80) {
        # 8% Uptake
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_80[[j+current]]$SY_ANAL)))
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_80[[j+current]]$SY_ANAL))/nrow(abm_sims_80[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_80[[j+current]]$SY_ANAL) - sum(abm_sims_80[[j+current]]$ORIGINAL_STATUS_SY_ANAL)))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_80[[j+current]]$SY_ANAL) - sum(abm_sims_80[[j+current]]$ORIGINAL_STATUS_SY_ANAL))/nrow(abm_sims_80[[j+current]]))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_80[[j+current]]$SERO_PREVENT, abm_sims_80[[j+current]]$COND_PREVENT_ANAL, na.rm = T)))
        stat_prevent_doxy_anal = c(stat_prevent_doxy_anal, (sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ANAL, na.rm = T)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_80[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_80[[j+current]]$COND_PREVENT_ANAL)))
        stat_prevent_doxy_anal_pct = c(stat_prevent_doxy_anal_pct, ((sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ANAL, na.rm = T))/(sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_80[[j+current]]$SERO_PREVENT, abm_sims_80[[j+current]]$COND_PREVENT_ANAL, na.rm = T))))
        stat_prevent_sero_anal_pct = c(stat_prevent_sero_anal_pct, ((sum(abm_sims_80[[j+current]]$SERO_PREVENT))/(sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_80[[j+current]]$SERO_PREVENT, abm_sims_80[[j+current]]$COND_PREVENT_ANAL, na.rm = T))))
        stat_prevent_cond_anal_pct = c(stat_prevent_cond_anal_pct, ((sum(abm_sims_80[[j+current]]$COND_PREVENT_ANAL))/(sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_80[[j+current]]$SERO_PREVENT, abm_sims_80[[j+current]]$COND_PREVENT_ANAL, na.rm = T))))
        
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_80[[j+current]]$SY_ORAL)))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_80[[j+current]]$SY_ORAL))/nrow(abm_sims_80[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_80[[j+current]]$SY_ORAL) - sum(abm_sims_80[[j+current]]$ORIGINAL_STATUS_SY_ORAL)))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_80[[j+current]]$SY_ORAL) - sum(abm_sims_80[[j+current]]$ORIGINAL_STATUS_SY_ORAL))/nrow(abm_sims_80[[j+current]]))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_80[[j+current]]$COND_PREVENT_ORAL, na.rm = T)))
        stat_prevent_doxy_oral = c(stat_prevent_doxy_oral, (sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_80[[j+current]]$COND_PREVENT_ORAL)))
        stat_prevent_doxy_oral_pct = c(stat_prevent_doxy_oral_pct, ((sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T))/(sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_80[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        stat_prevent_cond_oral_pct = c(stat_prevent_cond_oral_pct, ((sum(abm_sims_80[[j+current]]$COND_PREVENT_ORAL))/(sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_80[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        
        stat_total = c(stat_total, (sum(abm_sims_80[[j+current]]$SY_ANAL, abm_sims_80[[j+current]]$SY_ORAL)))
        stat_prev = c(stat_prev, (sum(abm_sims_80[[j+current]]$SY_ANAL, abm_sims_80[[j+current]]$SY_ORAL))/nrow(abm_sims_80[[j+current]]))
        stat_new = c(stat_new, (sum(abm_sims_80[[j+current]]$SY_ANAL, abm_sims_80[[j+current]]$SY_ORAL) - sum(abm_sims_80[[j+current]]$ORIGINAL_STATUS_SY_ANAL, abm_sims_80[[j+current]]$ORIGINAL_STATUS_SY_ORAL)))
        stat_incidence = c(stat_incidence, (sum(abm_sims_80[[j+current]]$SY_ANAL, abm_sims_80[[j+current]]$SY_ORAL) - sum(abm_sims_80[[j+current]]$ORIGINAL_STATUS_SY_ANAL, abm_sims_80[[j+current]]$ORIGINAL_STATUS_SY_ORAL))/nrow(abm_sims_80[[j+current]]))
        stat_prevent = c(stat_prevent, (sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_80[[j+current]]$SERO_PREVENT, abm_sims_80[[j+current]]$COND_PREVENT_ANAL, abm_sims_80[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_80[[j+current]]$COND_PREVENT_ORAL, na.rm = T)))
        stat_prevent_doxy = c(stat_prevent_doxy, (sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_80[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T)))
        stat_prevent_sero = c(stat_prevent_sero, (sum(abm_sims_80[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond = c(stat_prevent_cond, (sum(abm_sims_80[[j+current]]$COND_PREVENT_ANAL, abm_sims_80[[j+current]]$COND_PREVENT_ORAL)))
        stat_prevent_doxy_pct = c(stat_prevent_doxy_pct, ((sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_80[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T))/(sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_80[[j+current]]$SERO_PREVENT, abm_sims_80[[j+current]]$COND_PREVENT_ANAL, abm_sims_80[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_80[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        stat_prevent_sero_pct = c(stat_prevent_sero_pct, ((sum(abm_sims_80[[j+current]]$SERO_PREVENT))/(sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_80[[j+current]]$SERO_PREVENT, abm_sims_80[[j+current]]$COND_PREVENT_ANAL, abm_sims_80[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_80[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        stat_prevent_cond_pct = c(stat_prevent_cond_pct, ((sum(abm_sims_80[[j+current]]$COND_PREVENT_ANAL, abm_sims_80[[j+current]]$COND_PREVENT_ORAL))/(sum(abm_sims_80[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_80[[j+current]]$SERO_PREVENT, abm_sims_80[[j+current]]$COND_PREVENT_ANAL, abm_sims_80[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_80[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
      } else if (uptake[i]==100) {
        # 100% Uptake
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_100[[j+current]]$SY_ANAL)))
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_100[[j+current]]$SY_ANAL))/nrow(abm_sims_100[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_100[[j+current]]$SY_ANAL) - sum(abm_sims_100[[j+current]]$ORIGINAL_STATUS_SY_ANAL)))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_100[[j+current]]$SY_ANAL) - sum(abm_sims_100[[j+current]]$ORIGINAL_STATUS_SY_ANAL))/nrow(abm_sims_100[[j+current]]))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_100[[j+current]]$SERO_PREVENT, abm_sims_100[[j+current]]$COND_PREVENT_ANAL, na.rm = T)))
        stat_prevent_doxy_anal = c(stat_prevent_doxy_anal, (sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ANAL, na.rm = T)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_100[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_100[[j+current]]$COND_PREVENT_ANAL)))
        stat_prevent_doxy_anal_pct = c(stat_prevent_doxy_anal_pct, ((sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ANAL, na.rm = T))/(sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_100[[j+current]]$SERO_PREVENT, abm_sims_100[[j+current]]$COND_PREVENT_ANAL, na.rm = T))))
        stat_prevent_sero_anal_pct = c(stat_prevent_sero_anal_pct, ((sum(abm_sims_100[[j+current]]$SERO_PREVENT))/(sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_100[[j+current]]$SERO_PREVENT, abm_sims_100[[j+current]]$COND_PREVENT_ANAL, na.rm = T))))
        stat_prevent_cond_anal_pct = c(stat_prevent_cond_anal_pct, ((sum(abm_sims_100[[j+current]]$COND_PREVENT_ANAL))/(sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_100[[j+current]]$SERO_PREVENT, abm_sims_100[[j+current]]$COND_PREVENT_ANAL, na.rm = T))))
        
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_100[[j+current]]$SY_ORAL)))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_100[[j+current]]$SY_ORAL))/nrow(abm_sims_100[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_100[[j+current]]$SY_ORAL) - sum(abm_sims_100[[j+current]]$ORIGINAL_STATUS_SY_ORAL)))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_100[[j+current]]$SY_ORAL) - sum(abm_sims_100[[j+current]]$ORIGINAL_STATUS_SY_ORAL))/nrow(abm_sims_100[[j+current]]))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_100[[j+current]]$COND_PREVENT_ORAL, na.rm = T)))
        stat_prevent_doxy_oral = c(stat_prevent_doxy_oral, (sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_100[[j+current]]$COND_PREVENT_ORAL)))
        stat_prevent_doxy_oral_pct = c(stat_prevent_doxy_oral_pct, ((sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T))/(sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_100[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        stat_prevent_cond_oral_pct = c(stat_prevent_cond_oral_pct, ((sum(abm_sims_100[[j+current]]$COND_PREVENT_ORAL))/(sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_100[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        
        stat_total = c(stat_total, (sum(abm_sims_100[[j+current]]$SY_ANAL, abm_sims_100[[j+current]]$SY_ORAL)))
        stat_prev = c(stat_prev, (sum(abm_sims_100[[j+current]]$SY_ANAL, abm_sims_100[[j+current]]$SY_ORAL))/nrow(abm_sims_100[[j+current]]))
        stat_new = c(stat_new, (sum(abm_sims_100[[j+current]]$SY_ANAL, abm_sims_100[[j+current]]$SY_ORAL) - sum(abm_sims_100[[j+current]]$ORIGINAL_STATUS_SY_ANAL, abm_sims_100[[j+current]]$ORIGINAL_STATUS_SY_ORAL)))
        stat_incidence = c(stat_incidence, (sum(abm_sims_100[[j+current]]$SY_ANAL, abm_sims_100[[j+current]]$SY_ORAL) - sum(abm_sims_100[[j+current]]$ORIGINAL_STATUS_SY_ANAL, abm_sims_100[[j+current]]$ORIGINAL_STATUS_SY_ORAL))/nrow(abm_sims_100[[j+current]]))
        stat_prevent = c(stat_prevent, (sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_100[[j+current]]$SERO_PREVENT, abm_sims_100[[j+current]]$COND_PREVENT_ANAL, abm_sims_100[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_100[[j+current]]$COND_PREVENT_ORAL, na.rm = T)))
        stat_prevent_doxy = c(stat_prevent_doxy, (sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_100[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T)))
        stat_prevent_sero = c(stat_prevent_sero, (sum(abm_sims_100[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond = c(stat_prevent_cond, (sum(abm_sims_100[[j+current]]$COND_PREVENT_ANAL, abm_sims_100[[j+current]]$COND_PREVENT_ORAL)))
        stat_prevent_doxy_pct = c(stat_prevent_doxy_pct, ((sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_100[[j+current]]$DOXY_PREVENT_ORAL, na.rm = T))/(sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_100[[j+current]]$SERO_PREVENT, abm_sims_100[[j+current]]$COND_PREVENT_ANAL, abm_sims_100[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_100[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        stat_prevent_sero_pct = c(stat_prevent_sero_pct, ((sum(abm_sims_100[[j+current]]$SERO_PREVENT))/(sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_100[[j+current]]$SERO_PREVENT, abm_sims_100[[j+current]]$COND_PREVENT_ANAL, abm_sims_100[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_100[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
        stat_prevent_cond_pct = c(stat_prevent_cond_pct, ((sum(abm_sims_100[[j+current]]$COND_PREVENT_ANAL, abm_sims_100[[j+current]]$COND_PREVENT_ORAL))/(sum(abm_sims_100[[j+current]]$DOXY_PREVENT_ANAL, abm_sims_100[[j+current]]$SERO_PREVENT, abm_sims_100[[j+current]]$COND_PREVENT_ANAL, abm_sims_100[[j+current]]$DOXY_PREVENT_ORAL, abm_sims_100[[j+current]]$COND_PREVENT_ORAL, na.rm = T))))
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
    stat_prevent_sero_anal_pct = stat_prevent_sero_anal_pct[-1]
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
    stat_total = stat_total[-1]
    stat_prev = stat_prev[-1]
    stat_new = stat_new[-1]
    stat_incidence = stat_incidence[-1]
    stat_prevent = stat_prevent[-1]
    stat_prevent_doxy = stat_prevent_doxy[-1]
    stat_prevent_sero = stat_prevent_sero[-1]
    stat_prevent_cond = stat_prevent_cond[-1]
    stat_prevent_doxy_pct = stat_prevent_doxy_pct[-1]
    stat_prevent_sero_pct = stat_prevent_sero_pct[-1]
    stat_prevent_cond_pct = stat_prevent_cond_pct[-1]
    
    results_infections_anal = rbind(results_infections_anal, data.frame("Uptake"=uptake[i],"Adhere"=adhere[j],"Total_mean"=mean(stat_total_anal),"Total_SE"=(sd(stat_total_anal)/sqrt(length(stat_total_anal))),"Prev_mean"=mean(stat_prev_anal),"Prev_SE"=(sd(stat_prev_anal)/sqrt(length(stat_prev_anal))),"New_mean"=mean(stat_new_anal),"New_SE"=(sd(stat_new_anal)/sqrt(length(stat_new_anal))),"Incidence_mean"=mean(stat_incidence_anal),"Incidence_SE"=(sd(stat_incidence_anal)/sqrt(length(stat_incidence_anal))),"Prevent_mean"=mean(stat_prevent_anal),"Prevent_SE"=(sd(stat_prevent_anal)/sqrt(length(stat_prevent_anal))),"Doxy_prevent_mean"=mean(stat_prevent_doxy_anal),"Doxy_prevent_SE"=(sd(stat_prevent_doxy_anal)/sqrt(length(stat_prevent_doxy_anal))),"Sero_prevent_mean"=mean(stat_prevent_sero_anal),"Sero_prevent_SE"=(sd(stat_prevent_sero_anal)/sqrt(length(stat_prevent_sero_anal))),"Cond_prevent_mean"=mean(stat_prevent_cond_anal),"Cond_prevent_SE"=(sd(stat_prevent_cond_anal)/sqrt(length(stat_prevent_cond_anal))),"Total_prevent"=mean(stat_prevent_anal),"Doxy_prevent_pct"=mean(stat_prevent_doxy_anal_pct),"Doxy_prevent_pct_ll"=quantile(stat_prevent_doxy_anal_pct, 0.025),"Doxy_prevent_pct_ul"=quantile(stat_prevent_doxy_anal_pct, 0.975),"Sero_prevent_pct"=mean(stat_prevent_sero_anal_pct),"Sero_prevent_pct_ll"=quantile(stat_prevent_sero_anal_pct, 0.025),"Sero_prevent_pct_ul"=quantile(stat_prevent_sero_anal_pct, 0.975),"Cond_prevent_pct"=mean(stat_prevent_cond_anal_pct),"Cond_prevent_pct_ll"=quantile(stat_prevent_cond_anal_pct, 0.025),"Cond_prevent_pct_ul"=quantile(stat_prevent_cond_anal_pct, 0.975),stringsAsFactors=F))
    results_infections_oral = rbind(results_infections_oral, data.frame("Uptake"=uptake[i],"Adhere"=adhere[j],"Total_mean"=mean(stat_total_oral),"Total_SE"=(sd(stat_total_oral)/sqrt(length(stat_total_oral))),"Prev_mean"=mean(stat_prev_oral),"Prev_SE"=(sd(stat_prev_oral)/sqrt(length(stat_prev_oral))),"New_mean"=mean(stat_new_oral),"New_SE"=(sd(stat_new_oral)/sqrt(length(stat_new_oral))),"Incidence_mean"=mean(stat_incidence_oral),"Incidence_SE"=(sd(stat_incidence_oral)/sqrt(length(stat_incidence_oral))),"Prevent_mean"=mean(stat_prevent_oral),"Prevent_SE"=(sd(stat_prevent_oral)/sqrt(length(stat_prevent_oral))),"Doxy_prevent_mean"=mean(stat_prevent_doxy_oral),"Doxy_prevent_SE"=(sd(stat_prevent_doxy_oral)/sqrt(length(stat_prevent_doxy_oral))),"Cond_prevent_mean"=mean(stat_prevent_cond_oral),"Cond_prevent_SE"=(sd(stat_prevent_cond_oral)/sqrt(length(stat_prevent_cond_oral))),"Total_prevent"=mean(stat_prevent_oral),"Doxy_prevent_pct"=mean(stat_prevent_doxy_oral_pct),"Doxy_prevent_pct_ll"=quantile(stat_prevent_doxy_oral_pct,0.025),"Doxy_prevent_pct_ul"=quantile(stat_prevent_doxy_oral_pct, 0.975),"Cond_prevent_pct"=mean(stat_prevent_cond_oral_pct),"Cond_prevent_pct_ll"=quantile(stat_prevent_cond_oral_pct,0.025),"Cond_prevent_pct_ul"=quantile(stat_prevent_cond_oral_pct, 0.975),stringsAsFactors=F))
    results_infections = rbind(results_infections, data.frame("Uptake"=uptake[i],"Adhere"=adhere[j],"Total_mean"=mean(stat_total),"Total_SE"=(sd(stat_total)/sqrt(length(stat_total))),"Prev_mean"=mean(stat_prev),"Prev_SE"=(sd(stat_prev)/sqrt(length(stat_prev))),"New_mean"=mean(stat_new),"New_SE"=(sd(stat_new)/sqrt(length(stat_new))),"Incidence_mean"=mean(stat_incidence),"Incidence_SE"=(sd(stat_incidence)/sqrt(length(stat_incidence))),"Prevent_mean"=mean(stat_prevent),"Prevent_SE"=(sd(stat_prevent)/sqrt(length(stat_prevent))),"Doxy_prevent_mean"=mean(stat_prevent_doxy),"Doxy_prevent_SE"=(sd(stat_prevent_doxy)/sqrt(length(stat_prevent_doxy))),"Sero_prevent_mean"=mean(stat_prevent_sero),"Sero_prevent_SE"=(sd(stat_prevent_sero)/sqrt(length(stat_prevent_sero))),"Cond_prevent_mean"=mean(stat_prevent_cond),"Cond_prevent_SE"=(sd(stat_prevent_cond)/sqrt(length(stat_prevent_cond))),"Total_prevent"=mean(stat_prevent),"Doxy_prevent_pct"=mean(stat_prevent_doxy_pct),"Doxy_prevent_pct_ll"=quantile(stat_prevent_doxy_pct, 0.025),"Doxy_prevent_pct_ul"=quantile(stat_prevent_doxy_pct, 0.975),"Sero_prevent_pct"=mean(stat_prevent_sero_pct),"Sero_prevent_pct_ll"=quantile(stat_prevent_sero_pct, 0.025),"Sero_prevent_pct_ul"=quantile(stat_prevent_sero_pct, 0.975),"Cond_prevent_pct"=mean(stat_prevent_cond_pct),"Cond_prevent_pct_ll"=quantile(stat_prevent_cond_pct, 0.025),"Cond_prevent_pct_ul"=quantile(stat_prevent_cond_pct, 0.975),stringsAsFactors=F))
    
  }
}

rm(uptake, adhere, current, i, j, k, abm_sims_20, abm_sims_40, abm_sims_60, abm_sims_80, abm_sims_100, stat_incidence_anal, stat_new_anal, stat_prev_anal, stat_total_anal, stat_prevent_anal, stat_prevent_cond_anal, stat_prevent_doxy_anal, stat_prevent_sero_anal, stat_prevent_sero_anal_pct, stat_prevent_cond_anal_pct, stat_prevent_doxy_anal_pct, stat_incidence_oral, stat_new_oral, stat_prev_oral, stat_total_oral, stat_prevent_oral, stat_prevent_cond_oral, stat_prevent_doxy_oral, stat_prevent_cond_oral_pct, stat_prevent_doxy_oral_pct, stat_incidence, stat_new, stat_prev, stat_total, stat_prevent, stat_prevent_cond, stat_prevent_doxy, stat_prevent_sero, stat_prevent_sero_pct, stat_prevent_cond_pct, stat_prevent_doxy_pct)
results_infections_anal = results_infections_anal[-1, ]
results_infections_oral = results_infections_oral[-1, ]
results_infections = results_infections[-1, ]

#contour plots of infections prevented
p = ggplot(results_infections, aes(Uptake, Adhere, z= Doxy_prevent_pct)) +
  geom_contour(aes(colour = ..level..)) +
  stat_contour(color = "black") +
  labs(x="Uptake Coverage (%)", y="Adherence Level (%)") + 
  theme_bw() +
  theme(text = element_text(color = "#22211d", size = 12),legend.text=element_text(size=12), plot.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 12), legend.title = element_text(size=12)) 
p = direct.label(p, list("far.from.others.borders", colour='black', hjust = 1, vjust = 1))

p1 = ggplot(results_infections_anal, aes(Uptake, Adhere, z= Doxy_prevent_pct)) +
  geom_contour(aes(colour = ..level..)) +
  stat_contour(color = "black") +
  labs(x="Uptake Coverage (%)", y="Adherence Level (%)") + 
  theme_bw() +
  theme(text = element_text(color = "#22211d", size = 12),legend.text=element_text(size=12), plot.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 12), legend.title = element_text(size=12)) 
p1 = direct.label(p1, list("far.from.others.borders", colour='black', hjust = 1, vjust = 1))

p2 = ggplot(results_infections_oral, aes(Uptake, Adhere, z= Doxy_prevent_pct)) +
  geom_contour(aes(colour = ..level..)) +
  stat_contour(color = "black") +
  labs(x="Uptake Coverage (%)", y="Adherence Level (%)") + 
  theme_bw() +
  theme(text = element_text(color = "#22211d", size = 12),legend.text=element_text(size=12), plot.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 12), legend.title = element_text(size=12)) 
p2 = direct.label(p2, list("far.from.others.borders", colour='black', hjust = 1, vjust = 1))

ggsave("contour_anal.jpeg", plot=p1, width = 8, height = 5.4, units = "in", dpi = 1200)
ggsave("contour_oral.jpeg", plot=p2, width = 8, height = 5.4, units = "in", dpi = 1200)
ggsave("contour_overall.jpeg", plot=p, width = 8, height = 5.4, units = "in", dpi = 1200)

rm(p, p1, p2, results_infections, results_infections_anal, results_infections_oral)

# Agents characteristics post-simulation 

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

#% with 10 or more partners
partners10 = ifelse(partners[-1]>=28.23, 1, 0)
mean(partners10)
rm(i,sex,partners,delta,partners10)

#sensitivity analysis: restrict sample to quartiles 
abm_sims_20 = c(abm_sims_50[rep(c(TRUE, FALSE), c(6, 20))])
abm_sims_40 = c(abm_sims_50[rep(c(TRUE, FALSE, TRUE, FALSE), c(1, 5, 5, 15))])
abm_sims_60 = c(abm_sims_50[rep(c(TRUE, FALSE, TRUE, FALSE), c(1, 10, 5, 10))])
abm_sims_80 = c(abm_sims_50[rep(c(TRUE, FALSE, TRUE, FALSE), c(1, 15, 5, 5))])
abm_sims_100 = c(abm_sims_50[rep(c(TRUE, FALSE, TRUE), c(1, 20, 5))])

# quartile 1: 0 - 7 
# quartile 2: 8 - 17
# quartile 3: 18 - 31
# quartile 4: 32 - 156

abm_sims_20 = lapply(abm_sims_20, subset, c(PARTNERS>31))
abm_sims_40 = lapply(abm_sims_40, subset,  c(PARTNERS>31))
abm_sims_60 = lapply(abm_sims_60, subset,  c(PARTNERS>31))
abm_sims_80 = lapply(abm_sims_80, subset,  c(PARTNERS>31))
abm_sims_100 = lapply(abm_sims_100, subset, c(PARTNERS>31))

# re-run STEP7: Estimate PIA for Condoms and Doxy for each quartile 

# save results_infection for each quartile into a separate datafile
# dta_qrt1 = results_infections[,c("Uptake","Adhere","Doxy_prevent_pct")]
# dta_qrt2 = results_infections[,c("Uptake","Adhere","Doxy_prevent_pct")]
# dta_qrt3 = results_infections[,c("Uptake","Adhere","Doxy_prevent_pct")]
# dta_qrt4 = results_infections[,c("Uptake","Adhere","Doxy_prevent_pct")]

# merge each quartile dataframe together
dta = rbind(dta_qrt1, dta_qrt2, dta_qrt3, dta_qrt4)
dta$part_qrt = rep(1:4, each = 30)

# plot contour plot for PIA by quartile of sexual partners 
p = ggplot(dta, aes(Uptake, Adhere, z= Doxy_prevent_pct)) +
  geom_contour(aes(colour = ..level..)) +
  stat_contour(color = "black") +
  labs(x="Uptake Coverage (%)", y="Adherence Level (%)") + 
  facet_wrap(~Partner_Q) +
  theme_bw() +
  theme(text = element_text(color = "#22211d", size = 12),legend.text=element_text(size=12), plot.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 12), legend.title = element_text(size=12)) 
p = direct.label(p, list("far.from.others.borders", colour='black', hjust = 1, vjust = 1))
ggsave("supple_fig_contour_plot.jpeg", plot=p, width = 8.5, height = 8.5, units = "in", dpi = 1200)
