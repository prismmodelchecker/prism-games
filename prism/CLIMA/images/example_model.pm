mdp

module scheduler
	turn : [1..3] init 1;
	n_tasks : [-1..2] init -1;
	[] n_tasks=-1 -> 0.5 : (n_tasks'=1) + 0.5 : (n_tasks'=2);
	[go1] n_tasks>0 & turn=1 -> (turn'=2);
	[go2] n_tasks>0 & turn=2 -> (turn'=3);
	[do] n_tasks>0 & turn=3 -> (turn'=1) & (n_tasks'=n_tasks-1);
endmodule

module agent1
	team1 : [1..2] init 1;
	[go1] true -> (team1'=1);
	[go1] true -> (team1'=2);
endmodule

module agent2 = agent1 [go1=go2,team1=team2,team2=team1]
endmodule

rewards "total"
	turn=3 & team1!=team2 : 0.3;	
	turn=3 & team1=team2 : 1;
endrewards
