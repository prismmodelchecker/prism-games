#!/usr/bin/python -O

import sys
N = 4; # default number of modules
if(len(sys.argv)>=2):
    N = int(sys.argv[1]);

smg = open("car.gen.smg", "w");
prop = open("car.gen.props", "w");

smg.write("// PRISM model for case study of car driving\n")
smg.write("// - The game consists of two players and will have three goals.\n")
smg.write("//\n")
smg.write("// This is a machine-generated file\n")
smg.write("// Author: Clemens Wiltsche\n")
smg.write("// Last modified: 25/01/2013\n")

smg.write("\n\n")

smg.write("smg\n\n")

# initialize player actions
smg.write("player p1\n")
smg.write("\t[move_slow_1_init], [move_fast_1_init], [move_slow_2_init], [move_fast_2_init],\n");
for i in xrange(1,N):
    if(i!=1):
        smg.write(",\n")
    smg.write("\t[move_slow_1_%i], [move_fast_1_%i], [move_slow_2_%i], [move_fast_2_%i], [honk_%i], [brake_%i], [ignore_%i]" % (i,i,i,i,i,i,i))
smg.write("\nendplayer\n\n")

smg.write("player p2\n")
for i in xrange(1,N):
    if(i!=1):
        smg.write(",\n")
    smg.write("\t[adult_%i], [kid_%i], [none_%i], [stay_%i], [retreat_%i], [dash_%i]" % (i,i,i,i,i,i))
smg.write(",\n\t[term]")
smg.write("\nendplayer\n\n")



# pedestrian positions
smg.write("const int SIDEWALK = 0;\n")
smg.write("const int STREET = 1;\n\n")

# reaction of car to pedestrian
smg.write("const int HONKED = 1;\n")
smg.write("const int BRAKED = 2;\n")
smg.write("const int STOPPED = 3;\n")
smg.write("const int IGNORED = 4;\n")
smg.write("const int NONE = 0;\n\n")

# positions in topology
for i in xrange(0,N+1):
    smg.write("const int POS_%i = %i;\n" % (i, i))


# the topology
# linear for now
for i in xrange(0,N):
    smg.write("\nconst int DEST_1_%i = POS_%i;\n" % (i, i+1))
    smg.write("const int DEST_2_%i = POS_%i;\n" % (i, i+1))


smg.write("\nglobal p : [1..2] init 1;\n") # players
smg.write("global car_position : [0..%i] init 0;\n" % (N))
smg.write("global car_slow : bool init false;\n")
smg.write("global accident : bool init false;\n")
smg.write("global stall : bool init false;\n")
smg.write("global pedestrian_position : [0..1] init 0;\n")
smg.write("global reaction : [0..4] init 0;\n\n")

# the initial field
smg.write("module field_init\n")
smg.write("\t[move_slow_1_init] p=1 & car_position=POS_0 & accident=false & stall=false -> (car_position'=DEST_1_0) & (car_slow'=true) & (p'=2);\n")
smg.write("\t[move_fast_1_init] p=1 & car_position=POS_0 & accident=false & stall=false -> (car_position'=DEST_1_0) & (car_slow'=false) & (p'=2);\n")
smg.write("\t[move_slow_2_init] p=1 & car_position=POS_0 & accident=false & stall=false -> (car_position'=DEST_2_0) & (car_slow'=true) & (p'=2);\n")
smg.write("\t[move_fast_2_init] p=1 & car_position=POS_0 & accident=false & stall=false -> (car_position'=DEST_2_0) & (car_slow'=false) & (p'=2);\n")
smg.write("endmodule\n\n")


# the base case
smg.write("module field_1\n\n")

smg.write("\t[adult_1] p=2 & car_position=POS_1 & pedestrian_position=SIDEWALK & accident=false & stall=false -> 0.05 : (pedestrian_position'=STREET) & (p'=1) + 0.95 : (p'=1);\n")
smg.write("\t[kid_1] p=2 & car_position=POS_1 & pedestrian_position=SIDEWALK & accident=false & stall=false -> 0.2 : (pedestrian_position'=STREET) & (p'=1) + 0.8 : (p'=1);\n")
smg.write("\t[none_1] p=2 & car_position=POS_1 & pedestrian_position=SIDEWALK & accident=false & stall=false -> (p'=1);\n\n")

smg.write("\t[honk_1] p=1 & car_position=POS_1 & pedestrian_position=STREET & accident=false & stall=false -> (reaction'=HONKED) & (p'=2);\n")
smg.write("\t[brake_1] p=1 & car_position=POS_1 & pedestrian_position=STREET & car_slow=false & accident=false & stall=false -> (reaction'=BRAKED) & (car_slow'=true) & (p'=2);\n")
smg.write("\t[brake_1] p=1 & car_position=POS_1 & pedestrian_position=STREET & car_slow=true & accident=false & stall=false -> (reaction'=STOPPED) & (p'=2);\n")
smg.write("\t[ignore_1] p=1 & car_position=POS_1 & pedestrian_position=STREET & accident=false & stall=false -> (reaction'=IGNORED) & (p'=2);\n\n")

#	//[stay] p=2 & car_position=1 & pedestrian_position=STREET & reaction=HONKED & accident=false & stall=false -> (accident'=true);
smg.write("\t[retreat_1] p=2 & car_position=POS_1 & pedestrian_position=STREET & reaction=HONKED & accident=false & stall=false -> 0.3 : (accident'=true) + 0.7 : (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n")
smg.write("\t[dash_1] p=2 & car_position=POS_1 & pedestrian_position=STREET & reaction=HONKED & accident=false & stall=false -> 0.2 : (accident'=true) + 0.8 : (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n\n")

smg.write("\t[stay_1] p=2 & car_position=POS_1 & pedestrian_position=STREET & reaction=BRAKED & accident=false & stall=false -> (accident'=true);\n")
smg.write("\t[retreat_1] p=2 & car_position=POS_1 & pedestrian_position=STREET & reaction=BRAKED & accident=false & stall=false -> 0.1 : (accident'=true) + 0.9 : (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n")
smg.write("\t[dash_1] p=2 & car_position=POS_1 & pedestrian_position=STREET & reaction=BRAKED & accident=false & stall=false -> 0.05 : (accident'=true) + 0.95 : (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n\n")

smg.write("\t[stay_1] p=2 & car_position=POS_1 & pedestrian_position=STREET & reaction=STOPPED & accident=false & stall=false -> (stall'=true);\n")
smg.write("\t[retreat_1] p=2 & car_position=POS_1 & pedestrian_position=STREET & reaction=STOPPED & accident=false & stall=false -> (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n")
smg.write("\t[dash_1] p=2 & car_position=POS_1 & pedestrian_position=STREET & reaction=STOPPED & accident=false & stall=false -> (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n\n")

smg.write("\t[stay_1] p=2 & car_position=POS_1 & pedestrian_position=STREET & reaction=IGNORED & accident=false & stall=false -> (accident'=true);\n")
smg.write("\t[retreat_1] p=2 & car_position=POS_1 & pedestrian_position=STREET & reaction=IGNORED & accident=false & stall=false -> 0.95 : (accident'=true) + 0.05 : (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n")
smg.write("\t[dash_1] p=2 & car_position=POS_1 & pedestrian_position=STREET & reaction=IGNORED & accident=false & stall=false -> 0.9 : (accident'=true) + 0.1 : (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n\n")

smg.write("\t[move_slow_1_1] p=1 & car_position=POS_1 & pedestrian_position=SIDEWALK & accident=false & stall=false -> (car_position'=DEST_1_1) & (car_slow'=true) & (p'=2);\n")
smg.write("\t[move_fast_1_1] p=1 & car_position=POS_1 & pedestrian_position=SIDEWALK & accident=false & stall=false -> (car_position'=DEST_1_1) & (car_slow'=false) & (p'=2);\n")
smg.write("\t[move_slow_2_1] p=1 & car_position=POS_1 & pedestrian_position=SIDEWALK & accident=false & stall=false -> (car_position'=DEST_2_1) & (car_slow'=true) & (p'=2);\n")
smg.write("\t[move_fast_2_1] p=1 & car_position=POS_1 & pedestrian_position=SIDEWALK & accident=false & stall=false -> (car_position'=DEST_2_1) & (car_slow'=false) & (p'=2);\n\n")

smg.write("endmodule\n\n")



# iterate through topology
for i in xrange(2,N):
    smg.write("module field_%i = field_1 [POS_1=POS_%i, DEST_1_1=DEST_1_%i, DEST_2_1=DEST_2_%i, adult_1=adult_%i, kid_1=kid_%i, none_1=none_%i, honk_1=honk_%i, brake_1=brake_%i, ignore_1=ignore_%i, retreat_1=retreat_%i, dash_1=dash_%i, stay_1=stay_%i, move_slow_1_1=move_slow_1_%i, move_fast_1_1=move_fast_1_%i, move_slow_2_1=move_slow_2_%i, move_fast_2_1=move_fast_2_%i] endmodule\n" % (i, i, i, i, i, i, i, i, i, i, i, i, i, i, i, i, i))

# terminal states
smg.write("\nmodule terminals\n")
smg.write("\t[term] accident=true -> (accident'=true);\n")
smg.write("\t[term] stall=true -> (stall'=true);\n")
smg.write("\t[term] car_position=POS_%i -> (car_position'=POS_%i);\n" % (N, N))
smg.write("endmodule\n\n")

# properties
smg.write("formula goal1 = (car_position=POS_%i);\n" % (N))
smg.write("formula goal2 = (accident=false);\n")
smg.write("formula goal3 = (stall=false);\n\n")

# rewards
smg.write("rewards\n")
smg.write("\taccident & !car_slow : 1000;\n")
smg.write("\taccident & car_slow : 50;\n")
smg.write("\tstall : 2;\n")
smg.write("endrewards\n\n")

# write the property - a single three-objective goal
prop.write("<<1>> P>=0.5 [ goal1 U P>=0.5 [ goal2 U goal3 ] ]");

# close down files - flush
smg.close()
prop.close()
