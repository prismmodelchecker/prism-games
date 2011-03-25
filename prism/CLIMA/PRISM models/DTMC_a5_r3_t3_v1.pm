dtmc

// parameters
const int n_resources = 3;
const int n_tasks = 3;
const int n_sensors = 5;


// sensor resources
const int resource1 = 1;
const int resource2 = 2;
const int resource3 = 1;
const int resource4 = 2;
const int resource5 = 3;


// network configuration
const int e12=1;
const int e13=1;
const int e14=1;
const int e15=1;

const int e21=e12;
const int e23=1;
const int e24=1;
const int e25=1;

const int e31=e13;
const int e32=e23;
const int e34=1;
const int e35=1;

const int e41=e14;
const int e42=e24;
const int e43=e34;
const int e45=1;

const int e51=e15;
const int e52=e25;
const int e53=e35;
const int e54=e45;




formula s1_sched = (turn1=1 | turn2=1 | turn3=1 | turn4=1 | turn5=1);
formula s2_sched = (turn1=2 | turn2=2 | turn3=2 | turn4=2 | turn5=2);
formula s3_sched = (turn1=3 | turn2=3 | turn3=3 | turn4=3 | turn5=3);
formula s4_sched = (turn1=4 | turn2=4 | turn3=4 | turn4=4 | turn5=4);
formula s5_sched = (turn1=5 | turn2=5 | turn3=5 | turn4=5 | turn5=5);
formula all_not_sched = !(s1_sched | s2_sched | s3_sched | s4_sched | s5_sched); 
formula all_sched = (s1_sched & s2_sched & s3_sched & s4_sched & s5_sched);

rewards "scheduling_steps"
    [] true : 1;
endrewards 

module scheduler // schedules the algorithm

	// 0-not initialised; 1-tasks initialised; 2-player scheduled; 3-player executing
	status : [0..5];

	// schedule
	turn1 : [0..5];
	turn2 : [0..5];
	turn3 : [0..5];
	turn4 : [0..5];
	turn5 : [0..5];

	
	// scheduling first turn
	[] status=0 & all_not_sched -> 1/5 : (turn1'=1) & (status'=1) 
				     + 1/5 : (turn1'=2) & (status'=1)
				     + 1/5 : (turn1'=3) & (status'=1) 
				     + 1/5 : (turn1'=4) & (status'=1) 
				     + 1/5 : (turn1'=5) & (status'=1);
	
	// scheduling second turn
	[] status=1 & s1_sched -> 1/4 : (turn2'=2) & (status'=2) + 1/4 : (turn2'=3) & (status'=2) + 1/4 : (turn2'=4) & (status'=2) + 1/4 : (turn2'=5) & (status'=2); 
	[] status=1 & s2_sched -> 1/4 : (turn2'=1) & (status'=2) + 1/4 : (turn2'=3) & (status'=2) + 1/4 : (turn2'=4) & (status'=2) + 1/4 : (turn2'=5) & (status'=2) ; 
	[] status=1 & s3_sched -> 1/4 : (turn2'=2) & (status'=2) + 1/4 : (turn2'=1) & (status'=2) + 1/4 : (turn2'=4) & (status'=2) + 1/4 : (turn2'=5) & (status'=2) ; 
	[] status=1 & s4_sched -> 1/4 : (turn2'=2) & (status'=2) + 1/4 : (turn2'=3) & (status'=2) + 1/4 : (turn2'=1) & (status'=2) + 1/4 : (turn2'=5) & (status'=2) ; 
	[] status=1 & s5_sched -> 1/4 : (turn2'=2) & (status'=2) + 1/4 : (turn2'=3) & (status'=2) + 1/4 : (turn2'=4) & (status'=2) + 1/4 : (turn2'=1) & (status'=2) ; 

	// scheduling third turn
	[] status=2 & s1_sched & s2_sched -> 1/3 : (turn3'=3) & (status'=3) + 1/3 : (turn3'=4) & (status'=3) + 1/3 : (turn3'=5) & (status'=3); 
	[] status=2 & s1_sched & s3_sched -> 1/3 : (turn3'=2) & (status'=3) + 1/3 : (turn3'=4) & (status'=3) + 1/3 : (turn3'=5) & (status'=3);
	[] status=2 & s1_sched & s4_sched -> 1/3 : (turn3'=3) & (status'=3) + 1/3 : (turn3'=2) & (status'=3) + 1/3 : (turn3'=5) & (status'=3);
	[] status=2 & s1_sched & s5_sched -> 1/3 : (turn3'=3) & (status'=3) + 1/3 : (turn3'=4) & (status'=3) + 1/3 : (turn3'=2) & (status'=3); 
	[] status=2 & s2_sched & s3_sched -> 1/3 : (turn3'=1) & (status'=3) + 1/3 : (turn3'=4) & (status'=3) + 1/3 : (turn3'=5) & (status'=3); 
	[] status=2 & s2_sched & s4_sched -> 1/3 : (turn3'=1) & (status'=3) + 1/3 : (turn3'=3) & (status'=3) + 1/3 : (turn3'=5) & (status'=3);
	[] status=2 & s2_sched & s5_sched -> 1/3 : (turn3'=1) & (status'=3) + 1/3 : (turn3'=3) & (status'=3) + 1/3 : (turn3'=4) & (status'=3);
	[] status=2 & s3_sched & s4_sched -> 1/3 : (turn3'=1) & (status'=3) + 1/3 : (turn3'=2) & (status'=3) + 1/3 : (turn3'=5) & (status'=3);
	[] status=2 & s3_sched & s5_sched -> 1/3 : (turn3'=1) & (status'=3) + 1/3 : (turn3'=2) & (status'=3) + 1/3 : (turn3'=4) & (status'=3);
	[] status=2 & s4_sched & s5_sched -> 1/3 : (turn3'=1) & (status'=3) + 1/3 : (turn3'=2) & (status'=3) + 1/3 : (turn3'=3) & (status'=3);

	
	// scheduling fourth turn
	[] status=3 & s1_sched & s2_sched & s3_sched -> 1/2 : (turn4'=4) & (status'=4) + 1/2 : (turn4'=5) & (status'=4);
	[] status=3 & s1_sched & s2_sched & s4_sched -> 1/2 : (turn4'=3) & (status'=4) + 1/2 : (turn4'=5) & (status'=4);
	[] status=3 & s1_sched & s2_sched & s5_sched -> 1/2 : (turn4'=3) & (status'=4) + 1/2 : (turn4'=4) & (status'=4);
	[] status=3 & s1_sched & s3_sched & s4_sched -> 1/2 : (turn4'=2) & (status'=4) + 1/2 : (turn4'=5) & (status'=4);
	[] status=3 & s1_sched & s3_sched & s5_sched -> 1/2 : (turn4'=2) & (status'=4) + 1/2 : (turn4'=4) & (status'=4);
	[] status=3 & s1_sched & s4_sched & s5_sched -> 1/2 : (turn4'=2) & (status'=4) + 1/2 : (turn4'=3) & (status'=4);
	[] status=3 & s2_sched & s3_sched & s4_sched -> 1/2 : (turn4'=1) & (status'=4) + 1/2 : (turn4'=5) & (status'=4);
	[] status=3 & s2_sched & s3_sched & s5_sched -> 1/2 : (turn4'=1) & (status'=4) + 1/2 : (turn4'=4) & (status'=4);
	[] status=3 & s2_sched & s4_sched & s5_sched -> 1/2 : (turn4'=1) & (status'=4) + 1/2 : (turn4'=3) & (status'=4);
	[] status=3 & s3_sched & s4_sched & s5_sched -> 1/2 : (turn4'=1) & (status'=4) + 1/2 : (turn4'=2) & (status'=4);

	// scheduling fifth turn
	[] status=4 & !s1_sched -> (turn5'=1) & (status'=5);
	[] status=4 & !s2_sched -> (turn5'=2) & (status'=5);
	[] status=4 & !s3_sched -> (turn5'=3) & (status'=5);
	[] status=4 & !s4_sched -> (turn5'=4) & (status'=5);
	[] status=4 & !s5_sched -> (turn5'=5) & (status'=5);


	[] status=5 -> (status'=5);

endmodule
















