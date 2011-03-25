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


module scheduler // schedules the algorithm

	// 0-not initialised; 1-tasks initialised; 2-player scheduled; 3-player executing
	status : [0..3];

	// sensor's turn 0-no 1-yes.
	turn1 : [0..1];
	turn2 : [0..1];
	turn3 : [0..1];
	turn4 : [0..1];
	turn5 : [0..1];


	// initialising tasks
	[initialise] status=0 -> (status'=1);

	// scheduling sensors randomly
	[scheduling] status=1 -> 0.2 : (turn1'=1) & (turn2'=0) & (turn3'=0) & (turn4'=0) & (turn5'=0) & (status'=2) 
			       + 0.2 : (turn1'=0) & (turn2'=1) & (turn3'=0) & (turn4'=0) & (turn5'=0) & (status'=2)
			       + 0.2 : (turn1'=0) & (turn2'=0) & (turn3'=1) & (turn4'=0) & (turn5'=0) & (status'=2)
			       + 0.2 : (turn1'=0) & (turn2'=0) & (turn3'=0) & (turn4'=1) & (turn5'=0) & (status'=2)
 			       + 0.2 : (turn1'=0) & (turn2'=0) & (turn3'=0) & (turn4'=0) & (turn5'=1) & (status'=2);

	// executing the schedule
	[str1] status=2 & turn1=1 -> (status'=3);
	[fin1] status=3 & turn1=1 -> (status'=1) & (turn1'=0) & (turn2'=0) & (turn3'=0) & (turn4'=0) & (turn5'=0);

	[str2] status=2 & turn2=1 -> (status'=3);
	[fin2] status=3 & turn2=1 -> (status'=1) & (turn1'=0) & (turn2'=0) & (turn3'=0) & (turn4'=0) & (turn5'=0);

	[str3] status=2 & turn3=1 -> (status'=3);
	[fin3] status=3 & turn3=1 -> (status'=1) & (turn1'=0) & (turn2'=0) & (turn3'=0) & (turn4'=0) & (turn5'=0);

	[str4] status=2 & turn4=1 -> (status'=3);
	[fin4] status=3 & turn4=1 -> (status'=1) & (turn1'=0) & (turn2'=0) & (turn3'=0) & (turn4'=0) & (turn5'=0);

	[str5] status=2 & turn5=1 -> (status'=3);
	[fin5] status=3 & turn5=1 -> (status'=1) & (turn1'=0) & (turn2'=0) & (turn3'=0) & (turn4'=0) & (turn5'=0);

endmodule

module task_generator // generates random set of tasks

	initial : [0..2];
	
	// resources required for sub-tasks for task 1
	T1_1 : [1..n_resources];
	T1_2 : [1..n_resources];

	// resources required for sub-tasks for task 2
	T2_1 : [1..n_resources];
	T2_2 : [1..n_resources];


	// initialising resources needed for task 1 at random
	[] initial=0 -> 1/9 : (T1_1'=1) & (T1_2'=1) & (initial'=1) 
		      + 1/9 : (T1_1'=1) & (T1_2'=2) & (initial'=1) 
		      + 1/9 : (T1_1'=1) & (T1_2'=3) & (initial'=1) 
		      + 1/9 : (T1_1'=2) & (T1_2'=1) & (initial'=1)
		      +	1/9 : (T1_1'=2) & (T1_2'=2) & (initial'=1) 
		      + 1/9 : (T1_1'=2) & (T1_2'=3) & (initial'=1) 
		      + 1/9 : (T1_1'=3) & (T1_2'=1) & (initial'=1) 
	 	      + 1/9 : (T1_1'=3) & (T1_2'=2) & (initial'=1)
		      + 1/9 : (T1_1'=3) & (T1_2'=3) & (initial'=1);

	// initialising resources needed for task 2 at random
	[initialise] initial=1 -> 1/9 : (T2_1'=1) & (T2_2'=1) 
		      		+ 1/9 : (T2_1'=1) & (T2_2'=2) 
		      		+ 1/9 : (T2_1'=1) & (T2_2'=3) 
		      		+ 1/9 : (T2_1'=2) & (T2_2'=1)
		      		+ 1/9 : (T2_1'=2) & (T2_2'=2) 
		      		+ 1/9 : (T2_1'=2) & (T2_2'=3) 
		      		+ 1/9 : (T2_1'=3) & (T2_2'=1) 
	 	      		+ 1/9 : (T2_1'=3) & (T2_2'=2)
		      		+ 1/9 : (T2_1'=3) & (T2_2'=3);
	

endmodule

module sensor1

	// sensor's status 0-UNCOMMITTED; 1-COMMITTED
	status1 : [0..1];

	// sensor's resource
	//resource1 : [1..n_resources];

	// sensor's teams, 0-non member; 1-member.
	M1_1 : [0..1];
	M1_2 : [0..1];


	// 1 - started turn; 2-selected task; 3 - finished turn
	state1 : [0..3];

	// number of tasks tested so far
	tasks_tested1 : [0..n_tasks];
	task_selected1 : [1..n_tasks];

	// initialising sensors resource at random
	//[initialise] true -> 1/n_resources : (resource1'=1) + 1/n_resources : (resource1'=2);

	// starting my turn
	[str1] true -> (state1'=1);

	// executing actions

	// if UNCOMMITED schedule the next task to test
	[] status1=0 & state1=1 & tasks_tested1<n_tasks -> 1/n_tasks : (task_selected1'=1) & (state1'=2) & (tasks_tested1'=tasks_tested1+1)
							 + 1/n_tasks : (task_selected1'=2) & (state1'=2) & (tasks_tested1'=tasks_tested1+1);
	

//// TASK 1 CHOSEN ///////
	// if sensor does not have resource required or the resource is already filled for task 1 choosing another task
	[] status1=0 & state1=2 & task_selected1=1 & (!has_resource_T1 | resource_filled_T1) -> (state1'=1);
	
	// if sensor has resource required for task 1 then with prob IP1 comitting to the team
	[] status1=0 & state1=2 & task_selected1=1 & team_size_T1=0 & has_resource_T1 -> (1-IP) : (state1'=1) + IP : (M1_1'=1) & (status1'=1) & (state1'=1);

	// if sensor has resource required for task1 and it is not yet filled committing to the team
	[] status1=0 & state1=2 & task_selected1=1 & team_size_T1!=0 & has_resource_T1 & !resource_filled_T1-> (M1_1'=1) & (status1'=1) & (state1'=1);	
//////////////////////////

	
//// TASK 2 CHOSEN ///////
	// if sensor does not have resource required or the resource is already filled for task 2 choosing another task
	[] status1=0 & state1=2 & task_selected1=2 & (!has_resource_T2 | resource_filled_T2) -> (state1'=1);
	
	// if sensor has resource required for task 2 then with prob IP1 comitting to the team
	[] status1=0 & state1=2 & task_selected1=2 & team_size_T2=0 & has_resource_T2 -> (1-IP) : (state1'=1) + IP : (M1_2'=1) & (status1'=1) & (state1'=1);

	// if sensor has resource required for task1 and it is not yet filled committing to the team
	[] status1=0 & state1=2 & task_selected1=2 & team_size_T2!=0 & has_resource_T2 & !resource_filled_T2-> (M1_2'=1) & (status1'=1) & (state1'=1);	
//////////////////////////


	// if tested all tasks - finish the turn
	[] status1=0 & state1=1 & tasks_tested1=n_tasks -> (state1'=3);
  

	// if COMMITTED finish the turn
	[] status1=1 & (state1=2 | state1=1) -> (state1'=3);

	// finishing my turn
	[fin1] state1=3 -> (state1'=0); 

endmodule

module sensor2 = sensor1 
[ 
	status1=status2, 
	M1_1=M2_1, 
	M1_2=M2_2, 
	state1=state2, 
	tasks_tested1=tasks_tested2, 
	task_selected1=task_selected2, 

	str1=str2,
	fin1=fin2,

	e12=e21,
	e13=e23,
	e14=e24,
	e15=e25,
	status2=status1,
	resource1=resource2,
	resource2=resource1,
	M2_1=M1_1,
	M2_2=M1_2
] 
endmodule

module sensor3 = sensor1 
[ 
	status1=status3, 
	M1_1=M3_1, 
	M1_2=M3_2, 	
	state1=state3, 
	tasks_tested1=tasks_tested3, 
	task_selected1=task_selected3, 

	str1=str3,
	fin1=fin3,

	e12=e31,
	e13=e32,
	e14=e34,
	e15=e35,
	status2=status1,
	status3=status2,
	resource1=resource3,
	M2_1=M1_1,
	M3_1=M2_1,
	M2_2=M1_2,
	M3_2=M2_2,
	resource2=resource1,
	resource3=resource2
] 
endmodule

module sensor4 = sensor1 
[ 
	status1=status4, 
	M1_1=M4_1, 
	M1_2=M4_2,
	state1=state4, 
	tasks_tested1=tasks_tested4, 
	task_selected1=task_selected4, 

	str1=str4,
	fin1=fin4,

	e12=e41,
	e13=e42,
	e14=e43,
	e15=e45,
	status2=status1,
	status3=status2,
	status4=status3,
	resource1=resource4,

	M2_1=M1_1,
	M3_1=M2_1,
	M4_1=M3_1,
	M2_2=M1_2,
	M3_2=M2_2,
	M4_2=M3_2,
	resource2=resource1,
	resource3=resource2,
	resource4=resource3,
	resource5=resource4
] 
endmodule


module sensor5 = sensor1 
[ 
	status1=status5, 
	M1_1=M5_1, 
	M1_2=M5_2, 
	state1=state5, 
	tasks_tested1=tasks_tested5, 
	task_selected1=task_selected5, 

	str1=str5,
	fin1=fin5,

	e12=e51,
	e13=e52,
	e14=e53,
	e15=e54,
	status2=status1,
	status3=status2,
	status4=status3,
	status5=status4,
	resource1=resource5,

	M2_1=M1_1,
	M3_1=M2_1,
	M4_1=M3_1,
	M5_1=M4_1,
	M2_2=M1_2,
	M3_2=M2_2,
	M4_2=M3_2,
	M5_2=M4_2,
	resource2=resource1,
	resource3=resource2,
	resource4=resource3,
	resource5=resource4
] 
endmodule

formula team_size_T1 = M1_1+M2_1+M3_1+M4_1+M5_1;
formula team_size_T2 = M1_2+M2_2+M3_2+M4_2+M5_2;

formula IP = (e12*(1-status2)+e13*(1-status3)+e14*(1-status4)+e15*(1-status5)) / (e12+e13+e14+e15);

formula has_resource_T1 = (T1_1=resource1 | T1_2=resource1);
formula has_resource_T2 = (T2_1=resource1 | T2_2=resource1);

formula exists_team_T1 = (e12*M2_1+e13*M3_1+e14*M4_1+e15*M5_1) > 0;
formula exists_team_T2 = (e12*M2_2+e13*M3_2+e14*M4_2+e15*M5_2) > 0;

formula resource_filled_T1 = (M2_1=1 & resource1=resource2) | (M3_1=1 & resource1=resource3) | (M4_1=1 & resource1=resource4) | (M5_1=1 & resource1=resource5);
formula resource_filled_T2 = (M2_2=1 & resource1=resource2) | (M3_2=1 & resource1=resource3) | (M4_2=1 & resource1=resource4) | (M5_2=1 & resource1=resource5);

// labels
label "T1_covered" = (M1_1*((resource1=T1_1)?1:0)+M2_1*((resource2=T1_1)?1:0)+M3_1*((resource3=T1_1)?1:0))+M4_1*((resource4=T1_1)?1:0)+M5_1*((resource5=T1_1)?1:0) > 0
		   & (M1_1*((resource1=T1_2)?1:0)+M2_1*((resource2=T1_2)?1:0)+M3_1*((resource3=T1_2)?1:0))+M4_1*((resource4=T1_2)?1:0)+M5_1*((resource5=T1_2)?1:0) > 0;

label "T2_covered" = (M1_2*((resource1=T2_1)?1:0)+M2_2*((resource2=T2_1)?1:0)+M3_2*((resource3=T2_1)?1:0))+M4_2*((resource4=T2_1)?1:0)+M5_2*((resource5=T2_1)?1:0) > 0
		   & (M1_2*((resource1=T2_2)?1:0)+M2_2*((resource2=T2_2)?1:0)+M3_2*((resource3=T2_2)?1:0))+M4_2*((resource4=T2_2)?1:0)+M5_2*((resource5=T2_2)?1:0) > 0;


label "s1_commited" = (status1=1); 
label "team_s1_s2" = (M1_1=1 & M2_1=1) | (M1_2=1 & M2_2=1);
label "team_s1_s2_s3" = (M1_1=1 & M2_1=1 & M3_1=1) | (M1_2=1 & M2_2=1 & M3_2=1);
label "team_s1_s2_s3_s4" = (M1_1=1 & M2_1=1 & M3_1=1 & M4_1=1) | (M1_2=1 & M2_2=1 & M3_2=1 & M4_2=1);
label "team_s1_s2_s3_s4_s5" = (M1_1=1 & M2_1=1 & M3_1=1 & M4_1=1 & M5_1=1) | (M1_2=1 & M2_2=1 & M3_2=1 & M4_2=1 & M5_2=1);

label "All_sensors_committed" = status1=1 & status2=1 & status3=1 & status4=1 & status5=1;

label "All_sensors_uncommitted" = status1=0 & status2=0 & status3=0 & status4=0 & status5=0;















