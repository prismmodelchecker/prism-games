// model of 2 sensors

mdp

// parameters
const int T_size = 2;
const int N_skills = 2;
const int N_tasks = 1;
const int N_sensors = 2;


// sensor skills
const int skill1 = 1;
const int skill2 = 2;


// network configuration (constant)
const int e12=1;

const int e21=e12;

// scheduling mode 0-random sensor, 1 random sequence.
const int sch_mode = 0;

module scheduler // schedules the algorithm

	// 0-not initialised; 1-tasks initialised; 2-player scheduled; 3-player executing
	status : [0..3];

	// sensor's turn 0-no 1-yes.
	turn1 : [0..1];
	turn2 : [0..1];

	// initialising tasks
	[initialise] status=0 -> (status'=1);

	// scheduling sensors mode 0
	[scheduling] sch_mode=0 & status=1 -> 0.5 : (turn1'=1) & (turn2'=0) & (status'=2)
				  	    + 0.5 : (turn1'=0) & (turn2'=1) & (status'=2);

	// first sensor's turn
	[str1] status=2 & turn1=1 -> (status'=3);
	[fin1] status=3 & turn1=1 -> (status'=1) & (turn1'=0) & (turn2'=0);


	// second sensor's turn
	[str2] status=2 & turn2=1 -> (status'=3);
	[fin2] status=3 & turn2=1 -> (status'=1) & (turn1'=0) & (turn2'=0) ;


endmodule

module task_generator // generates random set of tasks

	initial : [0..1];
	
	// skills required for sub-tasks for task 1
	T1_1 : [1..N_skills];
	T1_2 : [1..N_skills];

	// initialising skills needed for task 2 at random
//	[initialise] initial=0 -> 1/4 : (T1_1'=1) & (T1_2'=1) 
//		      		+ 1/4 : (T1_1'=1) & (T1_2'=2) 
//		      		+ 1/4 : (T1_1'=2) & (T1_2'=2) 
//		      		+ 1/4 : (T1_1'=2) & (T1_2'=1);
	

	[initialise] initial=0 ->  (T1_1'=1) & (T1_2'=2); 
		      		

endmodule

module sensor1

	// sensor's status 0-UNCOMMITTED; 1-COMMITTED
	status1 : [0..1];

	// sensor's teams, 0-non member; 1-member.
	M1_1 : [0..1];

	// 1 - started turn; 3 - finished turn
	state1 : [0..3];

	// starting turn
	[str1] true -> (state1'=1);
	

	// if sensor has skill required for a task then it can initiate the new team
	[] status1=0 & state1=1 & team_size_T1=0 & has_skill_T1 -> (M1_1'=1) & (status1'=1) & (state1'=3);

	// if sensor has skill required for a task and it is not yet filled committing it can join the team
	[] status1=0 & state1=1 & team_size_T1!=0 & has_skill_T1 & !skill_filled_T1 -> (M1_1'=1) & (status1'=1) & (state1'=3);
  
	[] state1=1 -> (state1'=3);

	// finishing turn
	[fin1] state1=3 -> (state1'=0); 

endmodule

module sensor2 = sensor1 
[ 
	status1=status2, 
	M1_1=M2_1, 
	state1=state2, 

	str1=str2,
	fin1=fin2,

	e12=e21,
	status2=status1,
	skill1=skill2,
	skill2=skill1,
	M2_1=M1_1
] 
endmodule


formula team_size_T1 = M1_1+M2_1;

formula IP = (e12*(1-status2)) / (e12);

formula has_skill_T1 = (T1_1=skill1 | T1_2=skill1);

formula exists_team_T1 = (e12*M2_1) > 0;

formula skill_filled_T1 = (M2_1=1 & skill1=skill2);
















