dtmc

// parameters
const int n_resources = 3;
const int n_tasks = 3;
const int n_sensors = 5;


// sensor resources
const int resource1 = 2;
const int resource2 = 1;
const int resource3 = 1;
const int resource4 = 1;
const int resource5 = 1;

module controller // schedules the algorithm

	// algorithm status
	status : [0..100];

	// task resource indicator variables
	t1_r1 : [0..1];
	t1_r2 : [0..1];
	t1_r3 : [0..1];
	
	t2_r1 : [0..1];
	t2_r2 : [0..1];
	t2_r3 : [0..1];

	// schedule placeholders
	turn1 : [0..5];
	turn2 : [0..5];
	turn3 : [0..5];
	turn4 : [0..5];
	turn5 : [0..5];

	// selecting schedule uniformly at random
	[] status=0 -> 1/120 : (turn1'=1) & (turn2'=2) & (turn3'=3) & (turn4'=4) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=2) & (turn3'=3) & (turn4'=5) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=2) & (turn3'=4) & (turn4'=3) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=2) & (turn3'=4) & (turn4'=5) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=2) & (turn3'=5) & (turn4'=3) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=2) & (turn3'=5) & (turn4'=4) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=3) & (turn3'=2) & (turn4'=4) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=3) & (turn3'=2) & (turn4'=5) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=3) & (turn3'=4) & (turn4'=2) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=3) & (turn3'=4) & (turn4'=5) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=3) & (turn3'=5) & (turn4'=2) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=3) & (turn3'=5) & (turn4'=4) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=4) & (turn3'=2) & (turn4'=3) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=4) & (turn3'=2) & (turn4'=5) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=4) & (turn3'=3) & (turn4'=2) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=4) & (turn3'=3) & (turn4'=5) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=4) & (turn3'=5) & (turn4'=2) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=4) & (turn3'=5) & (turn4'=3) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=5) & (turn3'=2) & (turn4'=3) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=5) & (turn3'=2) & (turn4'=4) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=5) & (turn3'=3) & (turn4'=2) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=5) & (turn3'=3) & (turn4'=4) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=5) & (turn3'=4) & (turn4'=2) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=1) & (turn2'=5) & (turn3'=4) & (turn4'=3) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=1) & (turn3'=3) & (turn4'=4) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=1) & (turn3'=3) & (turn4'=5) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=1) & (turn3'=4) & (turn4'=3) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=1) & (turn3'=4) & (turn4'=5) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=1) & (turn3'=5) & (turn4'=3) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=1) & (turn3'=5) & (turn4'=4) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=3) & (turn3'=1) & (turn4'=4) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=3) & (turn3'=1) & (turn4'=5) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=3) & (turn3'=4) & (turn4'=1) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=3) & (turn3'=4) & (turn4'=5) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=3) & (turn3'=5) & (turn4'=1) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=3) & (turn3'=5) & (turn4'=4) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=4) & (turn3'=1) & (turn4'=3) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=4) & (turn3'=1) & (turn4'=5) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=4) & (turn3'=3) & (turn4'=1) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=4) & (turn3'=3) & (turn4'=5) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=4) & (turn3'=5) & (turn4'=1) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=4) & (turn3'=5) & (turn4'=3) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=5) & (turn3'=1) & (turn4'=3) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=5) & (turn3'=1) & (turn4'=4) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=5) & (turn3'=3) & (turn4'=1) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=5) & (turn3'=3) & (turn4'=4) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=5) & (turn3'=4) & (turn4'=1) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=2) & (turn2'=5) & (turn3'=4) & (turn4'=3) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=1) & (turn3'=2) & (turn4'=4) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=1) & (turn3'=2) & (turn4'=5) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=1) & (turn3'=4) & (turn4'=2) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=1) & (turn3'=4) & (turn4'=5) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=1) & (turn3'=5) & (turn4'=2) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=1) & (turn3'=5) & (turn4'=4) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=2) & (turn3'=1) & (turn4'=4) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=2) & (turn3'=1) & (turn4'=5) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=2) & (turn3'=4) & (turn4'=1) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=2) & (turn3'=4) & (turn4'=5) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=2) & (turn3'=5) & (turn4'=1) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=2) & (turn3'=5) & (turn4'=4) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=4) & (turn3'=1) & (turn4'=2) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=4) & (turn3'=1) & (turn4'=5) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=4) & (turn3'=2) & (turn4'=1) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=4) & (turn3'=2) & (turn4'=5) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=4) & (turn3'=5) & (turn4'=1) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=4) & (turn3'=5) & (turn4'=2) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=5) & (turn3'=1) & (turn4'=2) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=5) & (turn3'=1) & (turn4'=4) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=5) & (turn3'=2) & (turn4'=1) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=5) & (turn3'=2) & (turn4'=4) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=5) & (turn3'=4) & (turn4'=1) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=3) & (turn2'=5) & (turn3'=4) & (turn4'=2) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=1) & (turn3'=2) & (turn4'=3) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=1) & (turn3'=2) & (turn4'=5) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=1) & (turn3'=3) & (turn4'=2) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=1) & (turn3'=3) & (turn4'=5) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=1) & (turn3'=5) & (turn4'=2) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=1) & (turn3'=5) & (turn4'=3) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=2) & (turn3'=1) & (turn4'=3) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=2) & (turn3'=1) & (turn4'=5) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=2) & (turn3'=3) & (turn4'=1) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=2) & (turn3'=3) & (turn4'=5) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=2) & (turn3'=5) & (turn4'=1) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=2) & (turn3'=5) & (turn4'=3) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=3) & (turn3'=1) & (turn4'=2) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=3) & (turn3'=1) & (turn4'=5) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=3) & (turn3'=2) & (turn4'=1) & (turn5'=5) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=3) & (turn3'=2) & (turn4'=5) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=3) & (turn3'=5) & (turn4'=1) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=3) & (turn3'=5) & (turn4'=2) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=5) & (turn3'=1) & (turn4'=2) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=5) & (turn3'=1) & (turn4'=3) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=5) & (turn3'=2) & (turn4'=1) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=5) & (turn3'=2) & (turn4'=3) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=5) & (turn3'=3) & (turn4'=1) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=4) & (turn2'=5) & (turn3'=3) & (turn4'=2) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=1) & (turn3'=2) & (turn4'=3) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=1) & (turn3'=2) & (turn4'=4) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=1) & (turn3'=3) & (turn4'=2) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=1) & (turn3'=3) & (turn4'=4) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=1) & (turn3'=4) & (turn4'=2) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=1) & (turn3'=4) & (turn4'=3) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=2) & (turn3'=1) & (turn4'=3) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=2) & (turn3'=1) & (turn4'=4) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=2) & (turn3'=3) & (turn4'=1) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=2) & (turn3'=3) & (turn4'=4) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=2) & (turn3'=4) & (turn4'=1) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=2) & (turn3'=4) & (turn4'=3) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=3) & (turn3'=1) & (turn4'=2) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=3) & (turn3'=1) & (turn4'=4) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=3) & (turn3'=2) & (turn4'=1) & (turn5'=4) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=3) & (turn3'=2) & (turn4'=4) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=3) & (turn3'=4) & (turn4'=1) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=3) & (turn3'=4) & (turn4'=2) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=4) & (turn3'=1) & (turn4'=2) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=4) & (turn3'=1) & (turn4'=3) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=4) & (turn3'=2) & (turn4'=1) & (turn5'=3) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=4) & (turn3'=2) & (turn4'=3) & (turn5'=1) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=4) & (turn3'=3) & (turn4'=1) & (turn5'=2) & (status'=1)
		 + 1/120 : (turn1'=5) & (turn2'=4) & (turn3'=3) & (turn4'=2) & (turn5'=1) & (status'=1);


	// initialising non-empty tasks uniformly at random
	[] status=1 -> 1/49 : (t1_r1'=0) & (t1_r2'=0) & (t1_r3'=1) & (t2_r1'=0) & (t2_r2'=0) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=0) & (t1_r3'=1) & (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=0) & (t1_r3'=1) & (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=0) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=0) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=0) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=0) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=0) & (t2_r1'=0) & (t2_r2'=0) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=0) & (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=0) & (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=0) & (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=0) & (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=0) & (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=0) & (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=1) & (t2_r1'=0) & (t2_r2'=0) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=1) & (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=1) & (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=0) & (t2_r1'=0) & (t2_r2'=0) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=0) & (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=0) & (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=0) & (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=0) & (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=0) & (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=0) & (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=1) & (t2_r1'=0) & (t2_r2'=0) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=1) & (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=1) & (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=0) & (t2_r1'=0) & (t2_r2'=0) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=0) & (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=0) & (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=0) & (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=0) & (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=0) & (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=0) & (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=1) & (t2_r1'=0) & (t2_r2'=0) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=1) & (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=1) & (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=1) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=0) & (status'=2)
		 + 1/49 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=1) & (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=1) & (status'=2);

	// executing the schedule

	// 1st round
	[str1] status=2 & turn1=1 -> (status'=2);
	[fin1] status=2 & turn1=1 -> (status'=3);
	[str2] status=2 & turn1=2 -> (status'=2);
	[fin2] status=2 & turn1=2 -> (status'=3);
	[str3] status=2 & turn1=3 -> (status'=2);
	[fin3] status=2 & turn1=3 -> (status'=3);
	[str4] status=2 & turn1=4 -> (status'=2);
	[fin4] status=2 & turn1=4 -> (status'=3);
	[str5] status=2 & turn1=5 -> (status'=2);
	[fin5] status=2 & turn1=5 -> (status'=3);

	// 2nd round
	[str1] status=3 & turn2=1 -> (status'=3);
	[fin1] status=3 & turn2=1 -> (status'=4);
	[str2] status=3 & turn2=2 -> (status'=3);
	[fin2] status=3 & turn2=2 -> (status'=4);
	[str3] status=3 & turn2=3 -> (status'=3);
	[fin3] status=3 & turn2=3 -> (status'=4);
	[str4] status=3 & turn2=4 -> (status'=3);
	[fin4] status=3 & turn2=4 -> (status'=4);
	[str5] status=3 & turn2=5 -> (status'=3);
	[fin5] status=3 & turn2=5 -> (status'=4);

	// 3rd round
	[str1] status=4 & turn3=1 -> (status'=4);
	[fin1] status=4 & turn3=1 -> (status'=5);
	[str2] status=4 & turn3=2 -> (status'=4);
	[fin2] status=4 & turn3=2 -> (status'=5);
	[str3] status=4 & turn3=3 -> (status'=4);
	[fin3] status=4 & turn3=3 -> (status'=5);
	[str4] status=4 & turn3=4 -> (status'=4);
	[fin4] status=4 & turn3=4 -> (status'=5);
	[str5] status=4 & turn3=5 -> (status'=4);
	[fin5] status=4 & turn3=5 -> (status'=5);

	// 4th round
	[str1] status=5 & turn4=1 -> (status'=5);
	[fin1] status=5 & turn4=1 -> (status'=6);
	[str2] status=5 & turn4=2 -> (status'=5);
	[fin2] status=5 & turn4=2 -> (status'=6);
	[str3] status=5 & turn4=3 -> (status'=5);
	[fin3] status=5 & turn4=3 -> (status'=6);
	[str4] status=5 & turn4=4 -> (status'=5);
	[fin4] status=5 & turn4=4 -> (status'=6);
	[str5] status=5 & turn4=5 -> (status'=5);
	[fin5] status=5 & turn4=5 -> (status'=6);

	// 5th round
	[str1] status=6 & turn5=1 -> (status'=6);
	[fin1] status=6 & turn5=1 -> (status'=7);
	[str2] status=6 & turn5=2 -> (status'=6);
	[fin2] status=6 & turn5=2 -> (status'=7);
	[str3] status=6 & turn5=3 -> (status'=6);
	[fin3] status=6 & turn5=3 -> (status'=7);
	[str4] status=6 & turn5=4 -> (status'=6);
	[fin4] status=6 & turn5=4 -> (status'=7);
	[str5] status=6 & turn5=5 -> (status'=6);
	[fin5] status=6 & turn5=5 -> (status'=7);

	[] status=7 -> (status'=7);


endmodule

module sensor1

	state1 : [0..4];

	// team membership indicators
	m1_t1 : [0..1];
	m1_t2 : [0..1];

	// task scheduling
	turn1_1 : [0..2];
	turn2_1 : [0..2];



	// starting turn, selecting order of tasks
	[str1] state1=0 -> 1/2 : (turn1_1'=1) & (turn2_1'=2) & (state1'=1)
			 + 1/2 : (turn1_1'=2) & (turn2_1'=1) & (state1'=1);

	// if there is not team and has required skill - initiating the team
	[] state1=1 & turn1_1=1 & !committed & team_size_t1=0 & has_resource_t1 -> (m1_t1'=1) & (state1'=2);
	[] state1=1 & turn1_1=1 & !committed & team_size_t1=0 & !has_resource_t1 -> (state1'=2);
	[] state1=1 & turn1_1=2 & !committed & team_size_t2=0 & has_resource_t2 -> (m1_t2'=1) & (state1'=2);
	[] state1=1 & turn1_1=2 & !committed & team_size_t2=0 & !has_resource_t2 -> (state1'=2);

	[] state1=2 & turn2_1=1 & !committed & team_size_t1=0 & has_resource_t1 -> (m1_t1'=1) & (state1'=3);
	[] state1=2 & turn2_1=1 & !committed & team_size_t1=0 & !has_resource_t1 -> (state1'=3);
	[] state1=2 & turn2_1=2 & !committed & team_size_t2=0 & has_resource_t2 -> (m1_t2'=1) & (state1'=3);
	[] state1=2 & turn2_1=2 & !committed & team_size_t2=0 & !has_resource_t2 -> (state1'=3);

	// if team already exists - joining the team 
	[] state1=1 & turn1_1=1 & !committed & team_size_t1>0 & has_resource_t1 & !resource_filled_t1 -> (m1_t1'=1) & (state1'=2);
	[] state1=1 & turn1_1=1 & !committed & team_size_t1>0 & (!has_resource_t1 | resource_filled_t1) -> (state1'=2);
	[] state1=1 & turn1_1=2 & !committed & team_size_t2>0 & has_resource_t2 & !resource_filled_t2 -> (m1_t2'=1) & (state1'=2);
	[] state1=1 & turn1_1=2 & !committed & team_size_t2>0 & (!has_resource_t2 | resource_filled_t2) -> (state1'=2);

	[] state1=2 & turn2_1=1 & !committed & team_size_t1>0 & has_resource_t1 & !resource_filled_t1 -> (m1_t1'=1) & (state1'=3);
	[] state1=2 & turn2_1=1 & !committed & team_size_t1>0 & (!has_resource_t1 | resource_filled_t1) -> (state1'=3);
	[] state1=2 & turn2_1=2 & !committed & team_size_t2>0 & has_resource_t2 & !resource_filled_t2 -> (m1_t2'=1) & (state1'=3);
	[] state1=2 & turn2_1=2 & !committed & team_size_t2>0 & (!has_resource_t2 | resource_filled_t2) -> (state1'=3);


	[fin1] state1=3 | committed  -> true; 

endmodule

module sensor2 = sensor1 
[ 
	state1=state2, 

	str1=str2,
	fin1=fin2,

	m1_t1=m2_t1,
	m1_t2=m2_t2,

	m2_t1=m1_t1,
	m2_t2=m1_t2,

	turn1_1=turn1_2,
	turn2_1=turn2_2,

	resource1=resource2,	
	resource2=resource1
] 
endmodule

module sensor3 = sensor1 
[ 
	state1=state3, 

	str1=str3,
	fin1=fin3,

	m1_t1=m3_t1,
	m1_t2=m3_t2,
	m3_t1=m1_t1,
	m3_t2=m1_t2,
	turn1_1=turn1_3,
	turn2_1=turn2_3,

	resource1=resource3,	
	resource3=resource1
] 
endmodule

module sensor4 = sensor1 
[ 
	state1=state4, 

	str1=str4,
	fin1=fin4,

	m1_t1=m4_t1,
	m1_t2=m4_t2,

	m4_t1=m1_t1,
	m4_t2=m1_t2,

	turn1_1=turn1_4,
	turn2_1=turn2_4,

	resource1=resource4,	
	resource4=resource1
] 
endmodule

module sensor5 = sensor1 
[ 
	state1=state5, 

	str1=str5,
	fin1=fin5,

	m1_t1=m5_t1,
	m1_t2=m5_t2,

	m5_t1=m1_t1,
	m5_t2=m1_t2,

	turn1_1=turn1_5,
	turn2_1=turn2_5,

	resource1=resource5,	
	resource5=resource1
] 
endmodule

// formulae for scheduling
formula s1_sched = (turn1=1 | turn2=1 | turn3=1 | turn4=1 | turn5=1);
formula s2_sched = (turn1=2 | turn2=2 | turn3=2 | turn4=2 | turn5=2);
formula s3_sched = (turn1=3 | turn2=3 | turn3=3 | turn4=3 | turn5=3);
formula s4_sched = (turn1=4 | turn2=4 | turn3=4 | turn4=4 | turn5=4);
formula s5_sched = (turn1=5 | turn2=5 | turn3=5 | turn4=5 | turn5=5);
formula all_not_sched = !(s1_sched | s2_sched | s3_sched | s4_sched | s5_sched); 
formula all_sched = (s1_sched & s2_sched & s3_sched & s4_sched & s5_sched);


// agent is committed to some team
formula committed = (m1_t1+m1_t2) > 0;

// formulae to compute team sizes
formula team_size_t1 = m1_t1+m2_t1+m3_t1+m4_t1+m5_t1;
formula team_size_t2 = m1_t2+m2_t2+m3_t2+m4_t2+m5_t2;

// formulae to check whether agent has the resource required by the task
formula has_resource_t1 = ( (t1_r1=1&resource1=1) | (t1_r2=1&resource1=2) | (t1_r3=1&resource1=3) );
formula has_resource_t2 = ( (t2_r1=1&resource1=1) | (t2_r2=1&resource1=2) | (t2_r3=1&resource1=3) );

// formulae to check whether the resource of an agent has been already filled in the team
formula resource_filled_t1 = (m2_t1=1 & resource1=resource2) | (m3_t1=1 & resource1=resource3) | (m4_t1=1 & resource1=resource4) | (m5_t1=1 & resource1=resource5);
formula resource_filled_t2 = (m2_t2=1 & resource1=resource2) | (m3_t2=1 & resource1=resource3) | (m4_t2=1 & resource1=resource4) | (m5_t2=1 & resource1=resource5);



// rewards
rewards "scheduling_steps"
    [] true : 1;
endrewards 

// labels
label "sensor1_joins_successful_team" = 
					(m1_t1=1 & (t1_r1=1 => 
							( resource1=1 | (m2_t1=1 & resource2=1) | (m3_t1=1 & resource3=1) | (m4_t1=1 & resource4=1) | (m5_t1=1 & resource5=1) )
						   )
						 & (t1_r2=1 =>
							( resource1=2 | (m2_t1=1 & resource2=2) | (m3_t1=1 & resource3=2) | (m4_t1=1 & resource4=2) | (m5_t1=1 & resource5=2) )
						   )
						 & (t1_r3=1 =>
							( resource1=3 | (m2_t1=1 & resource2=3) | (m3_t1=1 & resource3=3) | (m4_t1=1 & resource4=3) | (m5_t1=1 & resource5=3) )
						   )
  					) 
					|
					(m1_t2=1 & (t2_r1=1 => 
							( resource1=1 | (m2_t2=1 & resource2=1) | (m3_t2=1 & resource3=1) | (m4_t2=1 & resource4=1) | (m5_t2=1 & resource5=1) )
						   )
						 & (t2_r2=1 =>
							( resource1=2 | (m2_t2=1 & resource2=2) | (m3_t2=1 & resource3=2) | (m4_t2=1 & resource4=2) | (m5_t2=1 & resource5=2) )
						   )
						 & (t2_r3=1 =>
							( resource1=3 | (m2_t2=1 & resource2=3) | (m3_t2=1 & resource3=3) | (m4_t2=1 & resource4=3) | (m5_t2=1 & resource5=3) )
						   )
  					);
