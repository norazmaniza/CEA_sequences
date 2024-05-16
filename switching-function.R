rm(list=ls())
library(openxlsx)
library(diagram)
library(dampack)
library(darthtools)

#Import excel which contain input parameter------
excel_source<- file.choose()

#Import model parameter from excel
import<-read.xlsx(
  excel_source,
  "MALAYSIA",
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
  detectDates = FALSE,
  skipEmptyRows = TRUE,
  skipEmptyCols = TRUE,
  rows = NULL,
  cols = seq(2,15,1),
  check.names = FALSE,
  sep.names = "_",
  namedRegion = NULL,
  na.strings = "NA",
  fillMergedCells = FALSE
)
#mortality tables
mortality<-read.xlsx(
  excel_source,
  "Mortality_Malaysia",
  startRow = 2,
  colNames = TRUE,
  rowNames = FALSE,
  detectDates = FALSE,
  skipEmptyRows = TRUE,
  skipEmptyCols = TRUE,
  rows = NULL,
  cols = NULL,
  check.names = FALSE,
  sep.names = "_",
  namedRegion = NULL,
  na.strings = "NA",
  fillMergedCells = FALSE
)
#Simulation list (permutations)
per<-read.xlsx(
  excel_source,
  "Permutation",
  startRow = 1,
  colNames = FALSE,
  rowNames = FALSE,
  detectDates = FALSE,
  skipEmptyRows = TRUE,
  skipEmptyCols = TRUE,
  rows = NULL,
  cols = NULL,
  check.names = FALSE,
  sep.names = "_",
  namedRegion = NULL,
  na.strings = "NA",
  fillMergedCells = FALSE
)
sensitivity<-read.xlsx(
  excel_source,
  "Sensitivity",
  startRow = 1,
  colNames = TRUE,
  rowNames = FALSE,
  detectDates = FALSE,
  skipEmptyRows = TRUE,
  skipEmptyCols = TRUE,
  rows = NULL,
  cols = NULL,
  check.names = FALSE,
  sep.names = ".",
  namedRegion = NULL,
  na.strings = "NA",
  fillMergedCells = FALSE
)

## Store the all permutation (simulation)
#number of simulation
n_sim<-length(per[,1])
for (i in 1:n_sim) {
  curr_data<-per[i,]
  curr_data<-curr_data[!is.na(curr_data)]
  if(exists("v_names_Str")){
    v_names_Str <- append(v_names_Str,paste(curr_data,collapse='-'))
  } else {
    v_names_Str<-paste(curr_data,collapse='-')
  }
}




######## sensitivity parameter-----------

#for one way DSA
my_params_range<-sensitivity[,c("pars","min","max")]
my_params_basecase<-as.list(sensitivity[,"base"])
names(my_params_basecase)<-sensitivity[,"pars"]

#For PSA
my_dists<-sensitivity[,"Distribution"]
my_parameterization_types<-sensitivity[,"Params"]
my_dists_params<-type.convert(strsplit(sensitivity[,11],","), as.is=TRUE)
my_params<-as.vector(sensitivity[,"pars"])

#Main function for switching simulation----
switching_sim<-function(
    my_params_basecase
  )
{
  with(as.list(my_params_basecase), {
## General setup

### Cycle
cycle_length<-import["cycle_length",1] # cycle length 
n_cycles<-ceiling(import["n_cycles",1]) # time horizon

### Strategies 
v_names_sim_str<-colnames(import) #names of strategies available
n_str<-length(v_names_sim_str) #number of strategies available

#search for row names for transition probability
rownum<-which(rownames(import)=="p_RR"):which(rownames(import)=="p_NN")
t_prob_str<-import[rownum,] #matrix of all probabilities for each strategy

p_R_sen<-c(p_R_sec,p_R_ust,p_R_ada,p_R_ada_bs,p_R_sys,p_R_photo)### sort manually ## #assign value from sensitivity run
rr_discon_sen<-c(rr_discon_sec,rr_discon_ust,rr_discon_ada,rr_discon_ada_bs,rr_discon_sys,rr_discon_photo) #assign value from sensitivity run
r_PASI75<- -log(1-p_R_sen)/(16/52)
r_16wk_discon<-r_PASI75-(r_PASI75*rr_discon_sen)

t_prob_str[1,]<-1-exp(-r_16wk_discon*(16/52)) #
t_prob_str[2,]<-1-t_prob_str[1,]

# number of health state per strategy
v_names_hs<-c("R","N")
n_hs<-length(v_names_hs)

# Create empty vectors to store total utilities and costs 
v_tot_qaly<-v_tot_cost<-vector(mode ="numeric", length = length(per[,1]))
names(v_tot_qaly)<-names(v_tot_cost)<-v_names_Str

#create empty list to store cohort trace, vector of cost and utility
l_m_M<-l_u<-l_c<-sapply(v_names_Str, function(x) NULL)  

###Run simulation

for (run in 1:n_sim) {
  ##  switching simulation to run 
  #if(run==2){break} #limit the number of simulation
  sim=per[run,]
  sim<-sim[!is.na(sim)] #remove NA values
  n_str_sim=length(sim) #number of switching strategy for the current run
  
  #set default strategy matrix for this simulation, names
  v_names_states<-paste(rep("str",n_str_sim),1:n_str_sim,sep = "") #create default strategy names
  v_names_states<-paste(sort(rep(v_names_states,2)),rep(v_names_hs,n_str_sim),sep="_") #repeat and combine with health states
  v_names_states<-c(v_names_states,"Bsc","Death") #add BSC and death 
  n_names_sim<-length(v_names_states) #total number of health states
  
  #create template of matrix probability for the current simulation
  m_p<-matrix(data=0,nrow = n_names_sim, ncol = n_names_sim, dimnames = list(v_names_states,v_names_states)) #create probability matrix
  
  #enter individual simulated strategy probability matrix
  
  loop<-n_str_sim*n_hs+2 
  nloop=matrix(data=1:loop,nrow = n_hs) #matrix for finding position in m_p
  
  #loop for each strategy,store probability matrix
  for (i in 1:n_str_sim) {
    curr_str<-sim[i]
    curr_data<-t_prob_str[,curr_str]
    curr_t_prob<-matrix(data=curr_data,nrow = n_hs,ncol = n_hs, byrow=TRUE) 
    j<-nloop[1,i]
    k<-nloop[2,i] 
    m_p[j:k,j:k]<-curr_t_prob
  }
  
  #include BSC transition probability and death as default
  m_p["Bsc","Bsc"]<-m_p["Death","Death"]<-1

  #account for strategy switching
  
  #remove transition from NR strN to next
  for (i in 1:n_str_sim)
  {
    m_p[nloop[n_hs,i],nloop[n_hs,i]]<-0 # 
  }
  
  #get induction probability and enter at the NR row of each strategy
  for (i in 1:n_str_sim)
  {
    j<-nloop[2,i]
    k<-nloop[1,i+1]
    
    if (i==n_str_sim) { #if this is the last strategy, all cohort moves to BSC
      curr_data1<-1
      m_p[j,k]<-curr_data1
      break
    }
    curr_str<-as.matrix(sim[i+1]) #select the next strategy name
    curr_data1<-p_R_sen[match(curr_str,names(import))] #select induction probability of the strategy
    curr_data2<-1-curr_data1 #calculate non-response value
    
    m_p[j,k]<-curr_data1 #input to matrix
    m_p[j,k+1]<-curr_data2
  }

  #expand matrix to include death for each cycle (third dimension)
  a_m_p <- array(m_p, dim = c(n_names_sim, n_names_sim, n_cycles), 
                 dimnames = list(v_names_states, v_names_states, 0:(n_cycles - 1)))
  
  #store death probability 
  start_age<-import["n_age_init", 1]
  end_age<-import["n_age_max", 1]
  v_dead<-mortality[1:n_cycles,"p_w16"]
  
  #assign death probability to all 
  for (i in 1:n_cycles){

    #attribute death to all response
    rstate<-rowSums(a_m_p[,,i]!=0) #get number of state present
    a_m_p[,,i]<-replace(a_m_p[,,i], a_m_p[,,i]==0, NA) #change 0 to NA
    for (n in 1:n_names_sim) { #for each row
      #minus probability with death(divided by num of probability present)
      a_m_p[n,,i]<-a_m_p[n,,i]-(v_dead[i]/rstate[n])
    }
    
    #replace NA with 0
    a_m_p[,,i]<-replace(a_m_p[,,i], is.na(a_m_p[,,i]), 0)
    
    #insert death probability at cycle i
    a_m_p[,'Death',i]<-v_dead[i] 
    
    #consideration for last strategy transition to BSC
    a_m_p[n_str_sim*2,'Bsc',i]<-1-v_dead[i]
    
    #consideration for BSC and Death probability
    a_m_p['Bsc','Bsc',i]<-1-v_dead[i]
    a_m_p['Death','Death',i]<-1
    
  }
  
  #matrix probability done
  
  #set initial cohort trace vector 
  ## Initial state vector
  v_m_init<-rep(0,n_names_sim) #set a blank vector
  v_m_init[1]<-import["p_R",as.character(sim[1])]*import["n_cohort",as.character(sim[1])] # initial state vector based on proportion of patients achieved PASI75 for chosen simulation
  v_m_init[2]<-import["n_cohort",as.character(sim[1])]-v_m_init[1]
  
  ## Initialize cohort trace for strategy
  m_M<-matrix(NA,nrow = (n_cycles+1), ncol = n_names_sim, dimnames = list(0:n_cycles, v_names_states))                       
  m_M[1,]<-v_m_init #Store the initial state vector in the first row of the cohort trace

  ### Discount weight for costs and effects 
  v_dwe     <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
  v_dwc     <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
  
  
  # Run Markov model
  # Iterative solution of time-independent cSTM
  for (t in 1:n_cycles){  
    
    m_M [t + 1, ] <- round(m_M [t, ] %*% a_m_p[,,t],7) #also consider third dimension (death with respect to cycle)
  }
  
  ###utilities
  
  #import utility 
  
  #create matrix for utility per cycle
  v_u<-matrix(0, nrow = n_cycles+1, ncol = n_names_sim, dimnames = list(0:n_cycles, v_names_states))
  
  #search for row names for utility data
  rownum<-which(rownames(import)=="u_R"):which(rownames(import)=="u_N")
  t_util_str<-import[rownum,] #matrix of all utility for each strategy
  t_util_str[1,]<-c(rep(u_R,4),rep(u_sys_R,2)) #assign value from sensitivity run
  t_util_str[2,]<-c(rep(u_N,4),rep(u_sys_N,2)) #assign value from sensitivity run
  for (i in 1:n_str_sim){
    j<-nloop[1,i]
    v_u[,j]<-t_util_str["u_R",as.character(sim[i])]
    v_u[,j+1]<-t_util_str["u_N",as.character(sim[i])]
  }
  #utility for BSC
  v_u[,n_names_sim-1]<-u_bsc
  
  #scale by the cycle length
  v_u<-v_u*cycle_length
  
  ###costs
  
  #cost is different for each state and transition to next drug state
  #create template of matrix cost for the current scenario
  #m_c<-matrix(data=0,nrow = n_names_sim, ncol = n_names_sim, dimnames = list(v_names_states,v_names_states))
  c_unit_drug<-c(c_unit_drug_sec,c_unit_drug_ust,c_unit_drug_ada,c_unit_drug_ada_bs,c_unit_drug_sys,c_unit_drug_photo) #assign value from sensitivity run

  c_drug_Id<-c_unit_drug*import["f_drug_ind",] #assign value from sensitivity run

  c_admin_Id<-outp_cost*import["f_outp_ind",] #assign value from sensitivity run

  c_test_unit<-c(rep(c_test_unit_bio,4),c_test_unit_sys,c_test_unit_pho) #assign value from sensitivity run
  c_tests_Id<-outp_cost+ (c_test_unit*import["f_tests_ind",]) 
  
  c_soc_Id<-import["c_soc_Id",]
  
  c_tot_Id<-c_drug_Id+c_admin_Id+c_tests_Id+c_soc_Id

  c_drugs_Mt16<-c_unit_drug*import["f_drug_mt",] #assign value from sensitivity run

  c_admin_Mt16<-outp_cost*import["f_outp_mt",] #assign value from sensitivity run

  c_tests_Mt16<-c_tests_Id
  
  c_soc_Mt16<-import["c_soc_Mt16",]
  
  c_tot_Mt16<-c_drugs_Mt16+c_admin_Mt16+c_tests_Mt16+c_soc_Mt16
  
  c_bsc<-import["c_bsc",]
  c_death<-import["c_death",]
  
  #set matrix cost to be same structure as matrix probability array
  a_m_c<-a_m_p
  
  #loop for each strategy,store cost matrix
  for (i in 1:n_str_sim)
  {
    curr_str<-as.matrix(sim[i])
    j<-nloop[1,i]
    a_m_c[j,j,]<-a_m_c[j,j+1,]<-c_tot_Mt16[,curr_str] #assign maintenance R and NR cost
    if(i==1) {next} #if first strategy, no induction cost in the matrix
    a_m_c[j-1,j,]<-a_m_c[j-1,j+1,]<-c_tot_Id[,curr_str] #assign induction R and NR cost
  }

  #assign BSC cost and death cost, assume same similar cost for bsc and death cost in each state
  a_m_c["Bsc","Bsc",]<-a_m_c[n_str_sim*n_hs,"Bsc",]<-c_bsc[,1] #for last strategy transition to bsc and bsc maintenance 
  a_m_c[,"Death",]<-c_death[,1] #cost of death in all states
  
  #create vector for cost 
  v_c<-matrix(0, nrow = n_cycles+1, ncol = n_names_sim, dimnames = list(0:n_cycles, v_names_states))
  
  v_c[1,]<-v_m_init%*% a_m_c[,,1]
  #scale cost matrix with cohort distribution
  #calculate cost for all cycles
  for (t in 1:n_cycles){  
    v_c [t+1, ] <- m_M [t+1, ] %*% a_m_c[,,t] #
  }
  
  
  ### Expected QALYs and costs per cycle 
  ## Vector of QALYs and Costs
  # Apply state rewards 
  v_qaly_str<-m_M*v_u # sum the utilities of all states for each cycle
  v_qaly_str<-rowSums(v_qaly_str) 
  
  # sum the costs of all states for each cycle
  v_cost_str<-rowSums(v_c)
   
  #within cycle correction Simpson’s 1/3rd rule 
  
  # First, we define two functions to identify if a number is even or odd
  is_even <- function(x) x %% 2 == 0
  is_odd <- function(x) x %% 2 != 0
  ## Vector with cycles
  v_cycles <- seq(1, n_cycles + 1)
  ## Generate 2/3 and 4/3 multipliers for even and odd entries, respectively
  v_wcc <- is_even(v_cycles)*(2/3) + is_odd(v_cycles)*(4/3)
  ## Substitute 1/3 in first and last entries
  v_wcc[1] <- v_wcc[n_cycles + 1] <- 1/3
  
  ### Discounted total expected QALYs and Costs per strategy
  # QALYs
  v_tot_qaly[run]<-t(v_qaly_str) %*%(v_dwe * v_wcc)
  
  # Costs
  v_tot_cost[run]<-t(v_cost_str) %*%(v_dwc * v_wcc)
  
  #name of current strategy
  v_names_sim<-paste(sim,collapse='-')
  
}

##########ICER

## Incremental cost-effectiveness ratios (ICERs) 
df_cea<-dampack::calculate_icers(cost       = v_tot_cost, 
                                 effect     = v_tot_qaly,
                                 strategies = v_names_Str)
#print(" ICER table ")

#create  result table
df_cea_out=data.frame(Strategy= v_names_Str,Cost = v_tot_cost, Effect= v_tot_qaly,Inc_Cost=rep(0,n_sim),Inc_Effect=rep(0,n_sim), ICER = rep(0,n_sim), Status = rep(NA,n_sim),row.names = NULL)
df_cea_out$ICER<-df_cea$ICER[order(match(df_cea$Strategy,df_cea_out$Strategy))]
df_cea_out<-replace(df_cea_out,is.na(df_cea_out),0)

#ICER table for common comparator 
df_cea_ref<-df_cea_out

for (i in 2:length(v_names_Str)) {
  df_cea_ref[i,"Inc_Cost"]<-df_cea_ref[i,"Cost"]-df_cea_ref[1,"Cost"]
  df_cea_ref[i,"Inc_Effect"]<-df_cea_ref[i,"Effect"]-df_cea_ref[1,"Effect"]
  df_cea_ref[i,"ICER"]<-df_cea_ref[i,"Inc_Cost"] / df_cea_ref[i,"Inc_Effect"] 

}
#df_cea_ref<-df_cea_ref[order(df_cea_ref$Cost),] 

return(df_cea_ref)

  })
}

#show basecase
base_result<-switching_sim(my_params_basecase)
print("Base result analysis")
print(base_result)
Str<-make.names(v_names_Str)
write.xlsx(base_result, file = "Output_BaseResultICERTable.xlsx")

Sys.sleep(5)
stop("switching sim run with basecase is done ")





#ONE WAY DSA-------------------

dsa_oneway<-run_owsa_det(
  params_range = my_params_range,
  params_basecase = my_params_basecase,
  nsamp = 50,
  FUN = switching_sim,
  outcomes = c("Cost", "Effect", "ICER"),
  strategies = Str,
  progress = TRUE
)



#seperate dataframe
dsa_icer<-dsa_oneway[["owsa_ICER"]]
dsa_cost<-dsa_oneway[["owsa_Cost"]]
dsa_effect<-dsa_oneway[["owsa_Effect"]]

#start dsa oneway

wb <- createWorkbook()
for (i in 1:length(names(dsa_oneway))){
  addWorksheet(wb, names(dsa_oneway)[i])
  writeData(wb,i,dsa_oneway[[i]])
  #plot graph oneway
  dsa_plot<-plot(
    dsa_oneway[[i]],
    txtsize = 10,
    col = "full",
    facet_scales = "free",
    facet_nrow = NULL,
    facet_ncol = NULL,
    size = 1,
    n_x_ticks = 3,
    n_y_ticks = 3,
  )
  print(dsa_plot+ylab(unlist(strsplit(names(dsa_oneway)[i],"_"))[2])) #modify axis label
  print(paste("Saving DSA plot ... ",names(dsa_oneway)[i]))
  Sys.sleep(5)
  insertPlot(wb,i,startCol = 6, width = 29, height = 21 ,units = "cm")
}
saveWorkbook(wb, "DSA_oneway.xlsx", overwrite = TRUE)

#end of dsa oneway

#start dsa tornado

wb <- createWorkbook()
dsa_icer<-dsa_oneway[["owsa_ICER"]]
#plot tornado
for (i in 1:n_sim) {
sel_str<-Str[i]
addWorksheet(wb, sel_str)
tor<-owsa_tornado(
  dsa_icer[dsa_icer$strategy==sel_str,],
  return = "plot",
  txtsize = 16,
  min_rel_diff = 0.01,
  col = "full",
  n_y_ticks = 5,
  ylim = NULL,
  ybreaks = NULL
)
#plot name
plotname<-paste(paste("Tornado", chartr(".","-",sel_str),sep="_"),".png",sep="")
print(tor+xlab("Parameter")+ylab("ICER")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))  
print(paste("Saving Tornado plot ... ",sel_str))
Sys.sleep(2)
insertPlot(wb,sel_str,startCol = 9, width = 20, height = 20 ,units = "cm")

#export data
tor_data<-owsa_tornado(
  dsa_icer[dsa_icer$strategy==sel_str,],
  return = "data",
  txtsize = 16,
  min_rel_diff = 0.01, #filter out parameter which has small effect
  col = "full",
  n_y_ticks = 5,
  ylim = NULL,
  ybreaks = NULL
)
writeData(wb,sel_str,tor_data)
}
addWorksheet(wb, "Baseline_ICER")
writeData(wb,"Baseline_ICER",base_result)
saveWorkbook(wb, "DSA_tornado.xlsx", overwrite = TRUE)

#end dsa tornado

#stop()

#PSA ---------------------------

#generate sample
psa_samp<-gen_psa_samp(
  params = my_params,
  dists = my_dists,
  parameterization_types = my_parameterization_types,
  dists_params = my_dists_params,
  nsamp = 10000
)
wb <- createWorkbook()
addWorksheet(wb, "PSA_sample")
writeData(wb,"PSA_sample",psa_samp)

#run PSA
psa_output <- run_psa(psa_samp = psa_samp,
                      params_basecase = my_params_basecase,
                      FUN = switching_sim,
                      outcomes = c("Cost","Effect","ICER"),
                      strategies = v_names_Str,
                      progress = TRUE)

addWorksheet(wb, "PSA_output_ICER")
writeData(wb,"PSA_output_ICER",psa_output$ICER$other_outcome)

#make PSA object
cea_psa <- make_psa_obj(cost = psa_output$Cost$other_outcome, 
                        effect = psa_output$Effect$other_outcome, 
                        parameters = psa_output$Cost$parameters,
                        strategies = psa_output$Cost$strategies,
                        currency = "USD")

addWorksheet(wb, "PSA_cea")
writeData(wb,"PSA_cea",cea_psa$cost,startCol = 1)
writeData(wb,"PSA_cea",cea_psa$effectiveness,startCol = cea_psa$n_strategies+3)
writeData(wb,"PSA_cea",cea_psa$parameters,startCol = (cea_psa$n_strategies+3)*2)

#PSA scatterplot
scatterplot<-plot(
    cea_psa,
     center = TRUE,
     ellipse = TRUE,
     alpha = 0.2,
     txtsize = 18,
     col = "full"
     ) 
print(scatterplot+xlab("Effectiveness (QALY)") )
addWorksheet(wb, "Scatterplot")
print(paste("Saving Scatterplot ... "))
Sys.sleep(3)
insertPlot(wb,"Scatterplot",startCol = 1, width = 25, height = 25 ,units = "cm")

#PSA summary
psa_sum <- summary(cea_psa, 
                   calc_sds = TRUE)

addWorksheet(wb, "PSA_summary")
writeData(wb,"PSA_summary",psa_sum)

#PSA Icers
psa_icers <- calculate_icers(cost = psa_sum$meanCost, 
                         effect = psa_sum$meanEffect, 
                         strategies = psa_sum$Strategy)

#recalculate for comparator
psa_cea_out<-data.frame(Strategy= psa_sum$Strategy,Cost = psa_sum$meanCost, Effect= psa_sum$meanEffect,Inc_Cost=rep(0,n_sim),Inc_Effect=rep(0,n_sim), ICER = rep(0,n_sim), Status = rep(NA,n_sim),row.names = NULL)
psa_cea_out$ICER<-psa_icers$ICER[order(match(psa_icers$Strategy,psa_cea_out$Strategy))]
psa_cea_out$Status<-psa_icers$Status[order(match(psa_icers$Strategy,psa_cea_out$Strategy))]
psa_cea_out<-replace(psa_cea_out,is.na(psa_cea_out),0)

for (i in 2:length(v_names_Str)) {
  psa_cea_out[i,"Inc_Cost"]<-psa_cea_out[i,"Cost"]-psa_cea_out[1,"Cost"]
  psa_cea_out[i,"Inc_Effect"]<-psa_cea_out[i,"Effect"]-psa_cea_out[1,"Effect"]
  psa_cea_out[i,"ICER"]<-psa_cea_out[i,"Inc_Cost"] / psa_cea_out[i,"Inc_Effect"] 

}
addWorksheet(wb, "PSA_ICERS")
writeData(wb,"PSA_ICERS",psa_cea_out)

#set WTP vector
my_wtp<-seq(0, 1e6, by = 1e2) #vector of 0 WTP until 1 million

#Cost-effectiveness Acceptability Curve
ceac_obj <- ceac(wtp = my_wtp, # WTP 
                 psa = cea_psa)
summary(ceac_obj)

addWorksheet(wb, "CEAC_data")
writeData(wb,"CEAC_data",ceac_obj)

#Plot CEAC
ceac_plot<-plot(ceac_obj, 
     frontier = FALSE,
     points=FALSE,
     currency = "USD",
     txtsize = 22,
     min_prob = 0.1,
     n_x_ticks=10,
     col = "full") #+ scale_x_continuous(expand = c(0, 0), limits = c(0, 1e3)) #rescale x axis according to max WTP
ceac_plot+ylab("Probability Cost-Effective") +geom_path(linewidth = 1)
Sys.sleep(5)
addWorksheet(wb, "CEAC_plot")
insertPlot(wb,"CEAC_plot",startCol = 1, width = 25, height = 25 ,units = "cm")

#expected loss
el <- calc_exp_loss(wtp = my_wtp, #enter WTP value
                    psa = cea_psa)

addWorksheet(wb, "EL_data")
writeData(wb,"EL_data",el)

el_plot<-plot(el, 
     log_y = FALSE,
     frontier = FALSE,
     points = FALSE,
     lsize = 1,
     txtsize = 18,
     currency = "USD",
     effect_units = "QALY",
     n_y_ticks = 8,
     n_x_ticks = 20,
     xbreaks = NULL,
     ybreaks = NULL,
     xlim = c(0, NA),
     ylim = NULL,
     col = "full")
el_plot+ scale_x_continuous(expand = c(0, 0), limits = c(0, 1e3))
Sys.sleep(5)
addWorksheet(wb, "EL_plot")
insertPlot(wb,"EL_plot",startCol = 1, width = 25, height = 25 ,units = "cm")

saveWorkbook(wb, "PSA.xlsx", overwrite = TRUE)
print("End of PSA analysis")










stop()




##TWO WAY NOT USED

#TWO WAY DSA-------------------

#set parameter
my_twsa_params_range <- data.frame(pars = c("c_unit_drug_sec", "c_unit_drug_ust"), #set two parameter to vary
                                   min = c(4493.76, 13118.78), # input min value for each
                                   max = c(6740.64, 19678.16 )) #input max value for each

#run two way
l_twsa_det <- run_twsa_det(params_range = my_twsa_params_range,
                           params_basecase = my_params_basecase,
                           nsamp = 50,
                           FUN = switching_sim,
                           outcomes = c("Cost", "Effect", "ICER"),
                           strategies = v_names_Str,
                           progress = TRUE)
#select dataset to plot
my_twsa_ICER <- l_twsa_det$twsa_ICER

# plot optimal strategy as a function of the two parameters varied in the two-way DSA
plot(my_twsa_ICER)



#not used
twoway<-create_dsa_twoway(
  my_twsa_params_range,
  effectiveness = NULL,
  v_names_Str,
  cost = NULL,
  currency = "$",
  other_outcome = NULL
)
#not used
twsa(
  twoway,
  param1 = "c_unit_drug_sec",
  param2 = "p_R_sec",
  ranges = NULL,
  nsamp = 100,
  outcome = "nmb",
  wtp = 104337, #CET in USD, estimated to be 3 times GDP
  strategies = NULL,
  poly.order = 2
)


