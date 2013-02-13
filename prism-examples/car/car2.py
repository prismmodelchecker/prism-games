#!/usr/bin/python -O

import sys

import osm2graph

filename = "map(4).osm"
if(len(sys.argv)>=2):
    filename = sys.argv[1];

G = osm2graph.read_osm(filename)

# calculate distances and assign to links
#edges_to_add = []
for e in G.edges(data=True):
    dist = osm2graph.calculateDistance(G.node[e[0]]["data"], G.node[e[1]]["data"])
    # fill in relevant tags
    name = ''
    if "name" in e[2]["data"].tags:
        name = e[2]["data"].tags["name"]
    oneway = False
    if "oneway" in e[2]["data"].tags and e[2]["data"].tags["oneway"]=="yes":
        oneway = True

    print name
    print oneway

    e[2]["data"] = {"dist": dist, "oneway": oneway}
    if not name=='':
        e[2]["data"]["name"] = name

# collapse graph to not contain several edges where no intersection is
to_remove = []
to_add = []
changed = True
while changed:
    changed = False
    for e in G.edges(data=True):
        edges = G.edges(e[1],data=True)
        if len(edges)==1 and not edges[0]==e:
            name = ''
            if "name" in e[2]["data"]:
                name = e[2]["data"]["name"]
            elif "name" in edges[0][2]["data"]:
                name = edges[0][2]["data"]["name"]
            oneway = e[2]["data"]["oneway"] or edges[0][2]["data"]["oneway"]
            G.add_edge(e[0],edges[0][1],attr_dict={'data': {'dist': e[2]["data"]["dist"]+edges[0][2]["data"]["dist"], 'name': name, 'oneway': oneway}})
            G.remove_edge(e[0],e[1])
            G.remove_edge(edges[0][0],edges[0][1])
            changed = True
            break

# for each edge also introduce the reverse edge if it does not already exist
for e in G.edges(data=True):
    print e
    if not e[2]["data"]["oneway"] and not G.has_edge(e[1],e[0]):
        G.add_edge(e[1],e[0],attr_dict=e[2])

for e in G.edges(data=True):
    print e

# one field per edge
N = len(G.edges())

# TODO: properly identify goal
goal = N-1

# build a dictionary of edges and the ids i use
edgeid = {}
i = 0
for e in G.edges():
    edgeid[e] = i
    i = i + 1


smg = open("car2.gen.smg", "w");
prop = open("car2.gen.props", "w");

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
smg.write("\t[move_slow_init], [move_fast_init],\n");
first = True
for e in G.edges():
    i = edgeid[e]
    if(not first):
        smg.write(",\n")
    first = False
    smg.write("\t")
    for j in xrange(0, len(G.edges([e[1]]))):
        smg.write("[move_slow_%i_%i], [move_fast_%i_%i], " % (j, i, j, i))
    if len(G.edges([e[1]]))==0:
        smg.write("[term_%i], " % (i))
    smg.write("[honk_%i], [brake_%i], [ignore_%i]" % (i,i,i))
smg.write("\nendplayer\n\n")

smg.write("player p2\n")
for i in xrange(0,N):
    if(i!=0):
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
smg.write("const int POS_init = %i;\n" % (N))
smg.write("const int POS_goal = %i;\n" % (N+1))
smg.write("const int POS_term = %i;\n" % (N+2))
for e in G.edges(data=True):
    i = edgeid[(e[0],e[1])]
    name = ''
    if "name" in e[2]["data"]:
        name = e[2]["data"]["name"]
    if not "name" == '':
        smg.write("const int POS_%i = %i;\t// name: %s\n" % (i, i, name))
    else:
        smg.write("const int POS_%i = %i;\n" % (i, i))

smg.write("\n")

# the topology
for e in G.edges(data=True):
    ie = edgeid[(e[0],e[1])]
    if len(G.edges([e[1]]))==0:
        if ie==goal:
            smg.write("const int DEST_%i = POS_goal;\n" % (ie))
        else:
            smg.write("const int DEST_%i = POS_term;\n" % (ie))
    else:
        j = 0
        for f in G.edges([e[1]], data=True):
            i_f = edgeid[(f[0],f[1])]
            if ie==goal:
                smg.write("const int DEST_%i_%i = POS_goal;\n" % (j, ie))
            else:
                name_to = ''
                if "name" in f[2]["data"]:
                    name_to = f[2]["data"]["name"]
                name_from = ''
                if "name" in e[2]["data"]:
                    name_from = e[2]["data"]["name"]
                if not name_to == '' and not name_from == '':
                    smg.write("const int DEST_%i_%i = POS_%i;\t// from %s to %s\n" % (j, ie, i_f, name_from, name_to))
                elif not name_to == '':
                    smg.write("const int DEST_%i_%i = POS_%i;\t// to %s\n" % (j, ie, i_f, name_to))
                elif not name_from == '':
                    smg.write("const int DEST_%i_%i = POS_%i;\t// from %s\n" % (j, ie, i_f, name_from))
                else:
                    smg.write("const int DEST_%i_%i = POS_%i;\n" % (j, ie, i_f))
            j = j + 1


# find maximum number of successors of any edge
maxSucc = 0
for e in G.edges():
    if len(G.edges([e[1]])) > maxSucc:
        maxSucc = len(G.edges([e[1]]))

# identify which edge has how many successors
successors = {}
for e in G.edges():
    succ = len(G.edges([e[1]]))
    if succ not in successors:
        successors[succ] = [e]
    else:
        new_succs = successors[succ]
        new_succs.append(e)
        successors[succ] = new_succs

smg.write("\nglobal p : [1..2] init 1;\n") # players
smg.write("global car_position : [0..%i] init POS_init;\n" % (N+2))
smg.write("global car_slow : bool init false;\n")
smg.write("global accident : bool init false;\n")
smg.write("global stall : bool init false;\n")
smg.write("global pedestrian_position : [0..1] init 0;\n")
smg.write("global reaction : [0..4] init 0;\n\n")

# the initial field
smg.write("module field_init\n")
smg.write("\t[move_slow_init] p=1 & car_position=POS_init & accident=false & stall=false -> (car_position'=POS_0) & (car_slow'=true) & (p'=2);\n")
smg.write("\t[move_fast_init] p=1 & car_position=POS_init & accident=false & stall=false -> (car_position'=POS_0) & (car_slow'=false) & (p'=2);\n")
smg.write("endmodule\n\n")


# the base cases - one for each number of successors
for i, succs in successors.iteritems():
    # k is the first edge with i successors
    k = edgeid[succs[0]]
    smg.write("module field_%i\n\n" % (k))

    smg.write("\t[adult_%i] p=2 & car_position=POS_%i & pedestrian_position=SIDEWALK & accident=false & stall=false -> 0.05 : (pedestrian_position'=STREET) & (p'=1) + 0.95 : (p'=1);\n" % (k, k))
    smg.write("\t[kid_%i] p=2 & car_position=POS_%i & pedestrian_position=SIDEWALK & accident=false & stall=false -> 0.2 : (pedestrian_position'=STREET) & (p'=1) + 0.8 : (p'=1);\n" % (k, k))
    smg.write("\t[none_%i] p=2 & car_position=POS_%i & pedestrian_position=SIDEWALK & accident=false & stall=false -> (p'=1);\n\n" % (k, k))

    smg.write("\t[honk_%i] p=1 & car_position=POS_%i & pedestrian_position=STREET & accident=false & stall=false -> (reaction'=HONKED) & (p'=2);\n" % (k, k))
    smg.write("\t[brake_%i] p=1 & car_position=POS_%i & pedestrian_position=STREET & car_slow=false & accident=false & stall=false -> (reaction'=BRAKED) & (car_slow'=true) & (p'=2);\n" % (k, k))
    smg.write("\t[brake_%i] p=1 & car_position=POS_%i & pedestrian_position=STREET & car_slow=true & accident=false & stall=false -> (reaction'=STOPPED) & (p'=2);\n" % (k, k))
    smg.write("\t[ignore_%i] p=1 & car_position=POS_%i & pedestrian_position=STREET & accident=false & stall=false -> (reaction'=IGNORED) & (p'=2);\n\n" % (k, k))

    #	//[stay] p=2 & car_position=1 & pedestrian_position=STREET & reaction=HONKED & accident=false & stall=false -> (accident'=true);
    smg.write("\t[retreat_%i] p=2 & car_position=POS_%i & pedestrian_position=STREET & reaction=HONKED & accident=false & stall=false -> 0.3 : (accident'=true) + 0.7 : (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n" % (k, k))
    smg.write("\t[dash_%i] p=2 & car_position=POS_%i & pedestrian_position=STREET & reaction=HONKED & accident=false & stall=false -> 0.2 : (accident'=true) + 0.8 : (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n\n" % (k, k))

    smg.write("\t[stay_%i] p=2 & car_position=POS_%i & pedestrian_position=STREET & reaction=BRAKED & accident=false & stall=false -> (accident'=true);\n" % (k, k))
    smg.write("\t[retreat_%i] p=2 & car_position=POS_%i & pedestrian_position=STREET & reaction=BRAKED & accident=false & stall=false -> 0.1 : (accident'=true) + 0.9 : (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n" % (k, k))
    smg.write("\t[dash_%i] p=2 & car_position=POS_%i & pedestrian_position=STREET & reaction=BRAKED & accident=false & stall=false -> 0.05 : (accident'=true) + 0.95 : (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n\n" % (k, k))

    smg.write("\t[stay_%i] p=2 & car_position=POS_%i & pedestrian_position=STREET & reaction=STOPPED & accident=false & stall=false -> (stall'=true);\n" %(k, k))
    smg.write("\t[retreat_%i] p=2 & car_position=POS_%i & pedestrian_position=STREET & reaction=STOPPED & accident=false & stall=false -> (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n" % (k, k))
    smg.write("\t[dash_%i] p=2 & car_position=POS_%i & pedestrian_position=STREET & reaction=STOPPED & accident=false & stall=false -> (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n\n" % (k, k))

    smg.write("\t[stay_%i] p=2 & car_position=POS_%i & pedestrian_position=STREET & reaction=IGNORED & accident=false & stall=false -> (accident'=true);\n" % (k, k))
    smg.write("\t[retreat_%i] p=2 & car_position=POS_%i & pedestrian_position=STREET & reaction=IGNORED & accident=false & stall=false -> 0.95 : (accident'=true) + 0.05 : (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n" % (k, k))
    smg.write("\t[dash_%i] p=2 & car_position=POS_%i & pedestrian_position=STREET & reaction=IGNORED & accident=false & stall=false -> 0.9 : (accident'=true) + 0.1 : (pedestrian_position'=SIDEWALK) & (reaction'=NONE) & (p'=1);\n\n" % (k, k))

    for j in xrange(0,i):
        smg.write("\t[move_slow_%i_%i] p=1 & car_position=POS_%i & pedestrian_position=SIDEWALK & accident=false & stall=false -> (car_position'=DEST_%i_%i) & (car_slow'=true) & (p'=2);\n" % (j, k, k, j, k))
        smg.write("\t[move_fast_%i_%i] p=1 & car_position=POS_%i & pedestrian_position=SIDEWALK & accident=false & stall=false -> (car_position'=DEST_%i_%i) & (car_slow'=false) & (p'=2);\n" % (j, k, k, j, k))

    if i == 0:
        smg.write("\t[term_%i] p=1 & car_position=POS_%i & pedestrian_position=SIDEWALK & accident=false & stall=false -> (car_position'=DEST_%i) & (p'=2);\n" % (k, k, k))
    
    smg.write("endmodule\n\n")



# iterate through topology
for i, succs in successors.iteritems():
    # if just one edge with i successors, then have it covered already
    if len(succs) > 1:
        l = edgeid[succs[0]]
        for e in succs[1:]: # start iterating at seccond edge
            k = edgeid[e]
            smg.write("module field_%i = field_%i [POS_%i=POS_%i, adult_%i=adult_%i, kid_%i=kid_%i, none_%i=none_%i, honk_%i=honk_%i, brake_%i=brake_%i, ignore_%i=ignore_%i, retreat_%i=retreat_%i, dash_%i=dash_%i, stay_%i=stay_%i" % (k, l, l, k, l, k, l, k, l, k, l, k, l, k, l, k, l, k, l, k, l, k))
            for j in xrange(0,i):
                smg.write(", move_slow_%i_%i=move_slow_%i_%i, move_fast_%i_%i=move_fast_%i_%i, DEST_%i_%i = DEST_%i_%i" % (j, l, j, k, j, l, j, k, j, l, j, k))
            if i == 0:
                smg.write(", term_%i=term_%i, DEST_%i=DEST_%i" % (l, k, l, k))
            smg.write("] endmodule\n")

# terminal states
smg.write("\nmodule terminals\n")
smg.write("\t[term] accident=true -> (accident'=true);\n")
smg.write("\t[term] stall=true -> (stall'=true);\n")
smg.write("\t[term] car_position=POS_goal -> (car_position'=POS_goal);\n")
smg.write("\t[term] car_position=POS_term -> (car_position'=POS_term);\n")
smg.write("endmodule\n\n")

# properties
smg.write("formula goal1 = (car_position=POS_goal);\n")
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
