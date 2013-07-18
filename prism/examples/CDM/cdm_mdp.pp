// MDP model implementing collective decision making algorithm of:
// F. Saffre and A. Simaitis. Host Selection through Collective Decision
// ACM Transactions on Autonomous and Adaptive Systems (TAAS). 2011. 
//
// In contrast with original (DTMC) model, all agents are allowed to 
// make a decision whether to explore or communicate modeled by non-determinism.
// 
// Model has to be built using PRISM preprocessor (http://www.prismmodelchecker.org/prismpp/)
// using the following command: prismpp cdm_dtmc.pp <N> <K> <L> > cmd_dtmc.pm, where 
// <N> - number of agents,
// <K> - number of hosting sites,
// <L> - number of confidence levels.
//
// Hosting site qualities and other constant model parameters should
// be adjusted directly in the model file.
//
// Expected model size is: K^2N * L^N
//
// Aistis Simaitis 23/06/11 

#const N#
#const K#
#const L#
mdp

// number of agents
const int N = #N#;

// number of sites
const int K = #K#;

// number of confidence levels
const int L = #L#;

// model parameters
const double Pexp = 0.5;
const double eta = 2.0;
const double gamma = 2.0;
const double lambda = 2.0;

// quality of the sites
#for i=1:K#
const double Q#i# = 1.0;
#end#

// confidence levels of agents
#for i=1:N#
global confidence#i# : [1..L];
#end#

// site preferences of agents
#for i=1:N#
global preference#i# : [1..K];
#end#

// agent definitions
#for i=1:N#
module player#i#

	// location of the player
	location#i# : [1..K];

	// exploring sites
	[] true ->	  0 : true
			#for j=1:K#			
			// -- evaluating site and changing preference with probability Pswitchxy
			  + 1/K * Pswitch#i#_#j# : (location#i#'=#j#) & (preference#i#'=#j#) & (confidence#i#'=1)
			  + 1/K * (1-Pswitch#i#_#j#) : (location#i#'=#j#)

			#end#;
	// communicating with other agents
	[] true ->	  0 : true
			#for j=1:i-1#
			// -- trying to communicate with agent
			  + Pmeet_p#i# * (location#i#=location#j#?1:0) * (preference#i#=preference#j#?1:0) : (confidence#i#'=inc_conf#i#) & (confidence#j#'=inc_conf#j#) // same site
			  + Pmeet_p#i# * (location#i#=location#j#?1:0) * (preference#i#=preference#j#?0:1) * Pwin#i#_#j# : (confidence#i#'=inc_conf#i#) & (preference#j#'=preference#i#) & (confidence#j#'=1) // win
			  + Pmeet_p#i# * (location#i#=location#j#?1:0) * (preference#i#=preference#j#?0:1) * (1-Pwin#i#_#j#) : (confidence#i#'=1) & (preference#i#'=preference#j#) & (confidence#j#'=inc_conf#j#) // lose
			
			#end#
			#for j=i+1:N#
			// -- trying to communicate with agent
			  + Pmeet_p#i# * (location#i#=location#j#?1:0) * (preference#i#=preference#j#?1:0) : (confidence#i#'=inc_conf#i#) & (confidence#j#'=inc_conf#j#) // same site
			  + Pmeet_p#i# * (location#i#=location#j#?1:0) * (preference#i#=preference#j#?0:1) * Pwin#i#_#j# : (confidence#i#'=inc_conf#i#) & (preference#j#'=preference#i#) & (confidence#j#'=1) // win
			  + Pmeet_p#i# * (location#i#=location#j#?1:0) * (preference#i#=preference#j#?0:1) * (1-Pwin#i#_#j#) : (confidence#i#'=1) & (preference#i#'=preference#j#) & (confidence#j#'=inc_conf#j#) // lose
			
			#end#
			// -- no other agents at the site
			  + (count_p#i#=1?1:0) : true;

endmodule

#end#





// formulae to increase agents' confidence levels
	
	#for i=1:N#
	formula inc_conf#i# = confidence#i#=L ? L : (confidence#i#+1);
	#end#

// formulae to count agents at sites

	// probability for agent to meet another agent at its location
	#for i=1:N#
	formula Pmeet_p#i# = (count_p#i#=1?0:1/(count_p#i#-1));
	#end#

	// the number of agents at agent's current location
	#for i=1:N#
	formula count_p#i# = #for j=1:K-1# location#i#=#j# ? count#j# : (#end#count#K##for j=1:K-1#) #end#;
	#end#

	// number of agents at each location
	#for i=1:K# 
	formula count#i# = #+ j=1:N# (location#j#=#i#?1:0) #end#;
	#end#

// formulae to get qualities of agents' preferred sites
	#for i=1:N#
	formula Q_p#i# = #for j=1:K-1# preference#i#=1 ? Q#j# : (#end#Q#K##for j=1:K-1#) #end#;
	#end#

// formulae for evaluating the sites (Pswitchij = prob of to switch from size i to site j).
	
	#for i=1:N# 
	#for j=1:K#
	formula Pswitch#i#_#j# = preference#i#=#j# ? 0 : pow(Q#j#, eta) / (pow(Q#j#, eta) + pow(Q_p#i#, eta));	
	#end# 

	#end#

// formulae for conducting tournaments

	#for i=1:N#
	#for j=1:i-1#
	formula Pwin#i#_#j# = 1-Pwin#j#_#i#;
	#end#
	#for j=i+1:N#
	formula Pwin#i#_#j# = (pow(Q_p#i#, lambda) * pow(confidence#i#, gamma)) / 
		((pow(Q_p#i#, lambda) * pow(confidence#i#, gamma))+(pow(Q_p#j#, lambda) * pow(confidence#j#, gamma)));
	#end#

	#end#

// labeling states
	
// -- formulae to generate labels
	
	// agreement on site
	#for i=1:K#
	formula all_prefer_#i# = #& j=1:N# preference#j#=#i# #end#;
	#end#

	// compute total confidence
	formula total_confidence = #+ j=1:N# confidence#j# #end#;
	
	// confidence measures
	formula all_max_conf = total_confidence/N = L;
	formula half_max_conf = ((#+ j=1:N# confidence#j#=L?1:0 #end#)/N) >= 0.5;
		
// -- labels

	// agreement on particular sites
	#for i=1:K#
	label "all_prefer_#i#" = all_prefer_#i#;
	#end#

	// all agents have max confidence
	label "all_max_conf" = all_max_conf;

	label "half_max_conf" = half_max_conf;	

	// agreement on a site
	label "decision_made" = #| j=1:K# all_prefer_#j# #end#;

	





