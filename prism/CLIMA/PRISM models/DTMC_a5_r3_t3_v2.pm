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

	t3_r1 : [0..1];
	t3_r2 : [0..1];
	t3_r3 : [0..1];

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


	// selecting task 1 uniformly at random
	[] status=1 -> 	  1/7 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=1) & (status'=2)
			+ 1/7 : (t1_r1'=1) & (t1_r2'=1) & (t1_r3'=0) & (status'=2)
			+ 1/7 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=1) & (status'=2)
			+ 1/7 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=1) & (status'=2)
			+ 1/7 : (t1_r1'=1) & (t1_r2'=0) & (t1_r3'=0) & (status'=2)
			+ 1/7 : (t1_r1'=0) & (t1_r2'=1) & (t1_r3'=0) & (status'=2)
			+ 1/7 : (t1_r1'=0) & (t1_r2'=0) & (t1_r3'=1) & (status'=2);

	[] status=2 -> 	  1/7 : (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=1) & (status'=3)
			+ 1/7 : (t2_r1'=1) & (t2_r2'=1) & (t2_r3'=0) & (status'=3)
			+ 1/7 : (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=1) & (status'=3)
			+ 1/7 : (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=1) & (status'=3)
			+ 1/7 : (t2_r1'=1) & (t2_r2'=0) & (t2_r3'=0) & (status'=3)
			+ 1/7 : (t2_r1'=0) & (t2_r2'=1) & (t2_r3'=0) & (status'=3)
			+ 1/7 : (t2_r1'=0) & (t2_r2'=0) & (t2_r3'=1) & (status'=3);

	[] status=3 -> 	  1/7 : (t3_r1'=1) & (t3_r2'=1) & (t3_r3'=1) & (status'=4)
			+ 1/7 : (t3_r1'=1) & (t3_r2'=1) & (t3_r3'=0) & (status'=4)
			+ 1/7 : (t3_r1'=1) & (t3_r2'=0) & (t3_r3'=1) & (status'=4)
			+ 1/7 : (t3_r1'=0) & (t3_r2'=1) & (t3_r3'=1) & (status'=4)
			+ 1/7 : (t3_r1'=1) & (t3_r2'=0) & (t3_r3'=0) & (status'=4)
			+ 1/7 : (t3_r1'=0) & (t3_r2'=1) & (t3_r3'=0) & (status'=4)
			+ 1/7 : (t3_r1'=0) & (t3_r2'=0) & (t3_r3'=1) & (status'=4);

	[] status=4 -> (status'=4);

endmodule


// formulae
formula s1_sched = (turn1=1 | turn2=1 | turn3=1 | turn4=1 | turn5=1);
formula s2_sched = (turn1=2 | turn2=2 | turn3=2 | turn4=2 | turn5=2);
formula s3_sched = (turn1=3 | turn2=3 | turn3=3 | turn4=3 | turn5=3);
formula s4_sched = (turn1=4 | turn2=4 | turn3=4 | turn4=4 | turn5=4);
formula s5_sched = (turn1=5 | turn2=5 | turn3=5 | turn4=5 | turn5=5);
formula all_not_sched = !(s1_sched | s2_sched | s3_sched | s4_sched | s5_sched); 
formula all_sched = (s1_sched & s2_sched & s3_sched & s4_sched & s5_sched);

// rewards
rewards "scheduling_steps"
    [] true : 1;
endrewards 














