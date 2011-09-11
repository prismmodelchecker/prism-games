// Model implementing demand-side energy management algorithm of:
// H. Hildmann and F. Saffre. Influence of Variable Supply and Load Flexibility on Demand-Side Management
// IEEE 8th International Conference on the European Energy Market (EEM). 2011. 
//
// In contrast with original (DTMC) model, some agents are allowed to 
// make a decision whether to execute their jobs modeled by non-determinism.
// 
// Model has to be built using PRISM preprocessor (http://www.prismmodelchecker.org/prismpp/)
// using the following command: prismpp dsm_mdp-dtmc.pp <N> <D> <d> <L> <PS> > dsm_mdp-dtmc.pm, where 
//
// <N> - number of households
// <D> - number of days
// <d> - number of deterministic households
// <L> - maximum job duration
// <PS> - 0.Pstart
//
//
// Aistis Simaitis 25/08/11 

#const N#
#const D#
#const d#
#const L#
#const PS#

smg

// number of households
const int N = #N#;

// number of days
const int D = #D#;

// number of time intervals in the day
const int K = 16;

// expected number of jobs per household per day
const int Exp_J = 9;

// cost limits for households
const double price_limit = 1.5;

// initiation probabilities for jobs (uuniform distribution)
#for i=1:L#
const double P_J#i# = 1/#L#;
#end# 

// probability of starting a task independently of the cost
const double P_start = 0.#PS#;

// distribution of the expected demand across intervals
const double D_K1 = 0.0614;
const double D_K2 = 0.0392;
const double D_K3 = 0.0304;
const double D_K4 = 0.0304;
const double D_K5 = 0.0355;
const double D_K6 = 0.0518;
const double D_K7 = 0.0651;
const double D_K8 = 0.0643;
const double D_K9 = 0.0625;
const double D_K10 = 0.0618;
const double D_K11 = 0.0614;
const double D_K12 = 0.0695;
const double D_K13 = 0.0887;
const double D_K14 = 0.1013;
const double D_K15 = 0.1005;
const double D_K16 = 0.0762;

// time limit
const int max_time = K*D+1;

// ---------------------------------------------------

// time counter
global time : [1..max_time];

// jobs of households
#for i=1:N#
global job#i# : [0..#L#];
#end#

// scheduling variable
global sched : [0..N];

// definition of scheduling module
module player0

	[] sched = 0 -> 1/N : (sched'=1)
#for i=2:N#
		      + 1/N : (sched'=#i#)
#end#;

endmodule

// definitions of deterministic households
#for i=1:d#
player p#i#
	player#i#, [start#i#], [backoff#i#]
endplayer

module player#i#
	
	job_arrived#i# : [0..#L#];

	[] sched=#i# & active = 0 & job#i# > 0 & time < max_time -> (job1'=new_j1)#for j=2:N# & (job#j#'=new_j#j#)#end# & (time'=time+1) & (sched'=0);

	// initiate the job with probability P_init
	[] sched=#i# & active = 0 & job#i# = 0 & time < max_time -> P_init*P_J1 : (job_arrived#i#'=1)
				       			    #for j=2:L#
							    + P_init*P_J#j# : (job_arrived#i#'=#j#)
							    #end#
				       			    + (1-P_init) : (job1'=new_j1)#for j=2:N# & (job#j#'=new_j#j#)#end# & (time'=time+1) & (sched'=0);

	// start job if cost below the limit
	[start#i#] sched=#i# & job_arrived#i# > 0 & price <= price_limit & time < max_time -> (job#i#'=job_arrived#i#) #for j=1:i-1# & (job#j#'=new_j#j#)#end##for j=i+1:N# & (job#j#'=new_j#j#)#end# & (job_arrived#i#'=0) & (time'=time+1) & (sched'=0);

	// back-off with probability 1-P_start
	[backoff#i#] sched=#i# & job_arrived#i# > 0 & price > price_limit & time < max_time ->   P_start : (job#i#'=job_arrived#i#) #for j=1:i-1# & (job#j#'=new_j#j#)#end##for j=i+1:N# & (job#j#'=new_j#j#)#end# & (job_arrived#i#'=0) & (time'=time+1) & (sched'=0)
						   + (1-P_start) : (job1'=new_j1)#for j=2:N# & (job#j#'=new_j#j#)#end# & (job_arrived#i#'=0) & (time'=time+1) & (sched'=0);
	// finished
	[] sched=#i# & time=max_time -> (job1'=new_j1)#for j=2:N# & (job#j#'=new_j#j#)#end# & (sched'=0);

endmodule

#end#

// definitions of non-deterministic households
#for i=d+1:N#
player p#i#
	player#i#, [start#i#], [backoff#i#], [nbackoff#i#]
endplayer

module player#i#
	
	job_arrived#i# : [0..#L#];

	[] sched=#i# & active = 0 & job#i# > 0 & time < max_time -> (job1'=new_j1)#for j=2:N# & (job#j#'=new_j#j#)#end# & (time'=time+1) & (sched'=0);

	// initiate the job with probability P_init
	[] sched=#i# & active = 0 & job#i# = 0 & time < max_time -> P_init*P_J1 : (job_arrived#i#'=1)
				       			    #for j=2:L#
							    + P_init*P_J#j# : (job_arrived#i#'=#j#)
							    #end#
				       			    + (1-P_init) : (job1'=new_j1)#for j=2:N# & (job#j#'=new_j#j#)#end# & (time'=time+1) & (sched'=0);

	// start job if cost below the limit
	[start#i#] sched=#i# & job_arrived#i# > 0 & price <= price_limit & time < max_time-> (job#i#'=job_arrived#i#) #for j=1:i-1# & (job#j#'=new_j#j#)#end##for j=i+1:N# & (job#j#'=new_j#j#)#end# & (job_arrived#i#'=0) & (time'=time+1) & (sched'=0);

	// back-off with probability 1-P_start
	[backoff#i#] sched=#i# & job_arrived#i# > 0 & price > price_limit & time < max_time->   P_start : (job#i#'=job_arrived#i#) #for j=1:i-1# & (job#j#'=new_j#j#)#end##for j=i+1:N# & (job#j#'=new_j#j#)#end# & (job_arrived#i#'=0) & (time'=time+1) & (sched'=0)
						   + (1-P_start) : (job1'=new_j1)#for j=2:N# & (job#j#'=new_j#j#)#end# & (job_arrived#i#'=0) & (time'=time+1) & (sched'=0);

	// don't back-off 
	[nbackoff#i#] sched=#i# & job_arrived#i# > 0 & price > price_limit & time < max_time -> (job#i#'=job_arrived#i#) #for j=1:i-1# & (job#j#'=new_j#j#)#end##for j=i+1:N# & (job#j#'=new_j#j#)#end# & (job_arrived#i#'=0) & (time'=time+1) & (sched'=0);

	// finished
	[] sched=#i# & time=max_time -> (job1'=new_j1)#for j=2:N# & (job#j#'=new_j#j#)#end# & (sched'=0);

endmodule

#end#



// probability to initiate the load
formula P_init = Exp_J *
		 (mod(time,K) = 1 ? D_K1 : 
		 (mod(time,K) = 2 ? D_K2 : 
		 (mod(time,K) = 3 ? D_K3 : 
		 (mod(time,K) = 4 ? D_K4 : 
		 (mod(time,K) = 5 ? D_K5 :
		 (mod(time,K) = 6 ? D_K6 :
		 (mod(time,K) = 7 ? D_K7 :
		 (mod(time,K) = 8 ? D_K8 :
		 (mod(time,K) = 9 ? D_K9 :
		 (mod(time,K) = 10 ? D_K10 :
		 (mod(time,K) = 11 ? D_K11 : 
		 (mod(time,K) = 12 ? D_K12 :
		 (mod(time,K) = 13 ? D_K13 :
		 (mod(time,K) = 14 ? D_K14 :
		 (mod(time,K) = 15 ? D_K15 : 
		 D_K16)))))))))))))));

// formula to compute current cost
formula jobs_running = (job1>0?1:0) #for i=2:N# + (job#i#>0?1:0)#end#;

// formula to identify say that only one agent is active
formula active = job_arrived1#for i=2:N# + job_arrived#i# #end#;

// formula to update job status
#for i=1:N#
formula new_j#i# = job#i#=0?0:job#i#-1;
#end#

formula price = jobs_running+1;


rewards "cost"
	sched!=0 : jobs_running*jobs_running;
endrewards

rewards "tasks_started"
#for i=1:N#
	sched!=0 & job#i#=1 : 1;
#end#
endrewards



#for i=1:N#
rewards "value#for j=1:i##j##end#"
	#for j=1:i#
	sched!=0 & job#j#>0 : 1/jobs_running;
	#end#
endrewards
#end#

rewards "common_value"
	sched!=0 : jobs_running=0?0:1/jobs_running;
endrewards

rewards "upfront_cost1"

	[start1] true : 1/(job_arrived1*price);
	[backoff1] true : P_start/(job_arrived1*price);
	[nbackoff1] true : 1/(job_arrived1*price);
	
endrewards

rewards "upfront_tcost"

#for i=1:d#
	[start#i#] true : 1/(job_arrived#i#*price);
	[backoff#i#] true : P_start/(job_arrived#i#*price);
#end#

#for i=d+1:N#
	[start#i#] true : 1/(job_arrived#i#*price);
	[backoff#i#] true : P_start/(job_arrived#i#*price);
	[nbackoff#i#] true : 1/(job_arrived#i#*price);

#end#
endrewards














