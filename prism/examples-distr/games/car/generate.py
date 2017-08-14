#!/usr/bin/python -O
from __future__ import division
import datetime, math, sys
import copy
import osm2graph

# default map:
filename = "islip"
goal = 46 # Middle Street (right) ... islip
init = 124 # Kidlington Road (top left) ... islip

if(len(sys.argv)>=2):
    if(int(sys.argv[1])==1):
        filename = "islip"
        goal = 46 # Middle Street (right) ... islip
        init = 124 # Kidlington Road (top left) ... islip
    elif(int(sys.argv[1])==2):
        filename = "charlton"
        goal = 33 # road right top ... charlton
        init = 43 # road left bottom ... charlton
    else:
        print "Unknown map ID specified ... usage:"
        print "Call generate.py n, where n stands for one of the following:"
        print "\t1 ... islip"
        print "\t2 ... charlton"
        print "If n is not specified, islip is taken as default."
        sys.exit(1);

smg = open("%s.prism" % (filename), "w");
prop = open("%s.props" % (filename), "w");
bash = open("%s.sh" % (filename), "w");

mapfile = "%s.osm" % (filename)

#########################################################################
# Prepare Graph                                                         #
#########################################################################

G = osm2graph.read_osm(mapfile)
usedvalues = {}
# calculate distances and assign to links
for e in G.edges(data=True):
    db = osm2graph.calculateDistanceAndBearing(G.node[e[0]]["data"], G.node[e[1]]["data"])
    # fill in relevant tags
    name = ''
    if "name" in e[2]["data"].tags:
        name = e[2]["data"].tags["name"]
    oneway = False
    if "oneway" in e[2]["data"].tags and e[2]["data"].tags["oneway"]=="yes":
        oneway = True
    lanes = 1
    if "lanes" in e[2]["data"].tags:
        lanes = e[2]["data"].tags["lanes"]
    value = 0; # give the road a value depending on its type
    values = {'motorway': 20,
              'motorway_link': 19,
              'trunk': 15,
              'trunk_link': 14,
              'primary': 10,
              'primary_link': 9,
              'secondary': 8,
              'secondary_link': 7,
              'tertiary': 6,
              'tertiary_link': 5,
              'living_street': 0.5,
              'pedestrian': 0, # never use
              'residential': 3,
              'unclassified': 1,
              'service': 1,
              'track': 0.2,
              'bus_guideway': 0, # never use
              'raceway': 0, # never use
              'road': 0.5}
    if "highway" in e[2]["data"].tags:
        if e[2]["data"].tags["highway"] in values:
            value = values[e[2]["data"].tags["highway"]]
            usedvalues[e[2]["data"].tags["highway"]] = True
    
    e[2]["data"] = {"dist": db[0], "init_bearing": db[1], "final_bearing": db[1], "oneway": oneway, "value": value, "name": name}

# first, remove zero-value edges (they are never used, e.g. train tracks)
to_remove = []
for e in G.edges(data=True):
    if e[2]["data"]["value"] == 0:
        to_remove.append(e)
for e in to_remove:
    G.remove_edge(e[0],e[1])

# collapse graph to not contain several edges where no intersection is
changed = True
while changed:
    changed = False
    for e in G.edges(data=True):
        edges = G.edges(e[1],data=True)
        if len(edges)==1 and not edges[0]==e:
            # need to prevent roads coming in to be ignored
            collapse = True
            for f in G.edges(data=True):
                if (f[0] != e[0]) and f[1]==e[1]:
                    collapse = False
            if(collapse):
                name = ''
                if e[2]["data"]["name"] is not "":
                    name = e[2]["data"]["name"]
                elif edges[0][2]["data"]["name"] is not "":
                    name = edges[0][2]["data"]["name"]
                oneway = e[2]["data"]["oneway"] or edges[0][2]["data"]["oneway"]
                # connect e to edges[0]
                init_bearing = e[2]["data"]["init_bearing"]
                final_bearing = edges[0][2]["data"]["final_bearing"]
                value = min(e[2]["data"]["value"], edges[0][2]["data"]["value"])
                G.add_edge(e[0],edges[0][1],attr_dict={'data': {'dist': e[2]["data"]["dist"]+edges[0][2]["data"]["dist"], 'name': name, 'oneway': oneway, 'init_bearing': init_bearing, 'final_bearing': final_bearing, 'value': value}})
                G.remove_edge(e[0],e[1])
                G.remove_edge(edges[0][0],edges[0][1])
                changed = True
                break

# for each edge also introduce the reverse edge
# moreover, assign which edge to go to for a u-turn
# note that lanes do not get assigned separate edges
for e in G.edges(data=True):
    if not e[2]["data"]["oneway"] and not G.has_edge(e[1],e[0]):
        new_data = copy.deepcopy(e[2])
        new_data["data"]["init_bearing"] = ((180 - e[2]["data"]["final_bearing"]) % 360)
        new_data["data"]["final_bearing"] = ((180 - e[2]["data"]["init_bearing"]) % 360)
        G.add_edge(e[1],e[0],attr_dict=e[2])

# identify which edge has how many successors
successors = {}
for e in G.edges(data=True):
    succ = 0;
    for f in G.edges([e[1]]):
        if (f[1],f[0]) == (e[0],e[1]): # no uturn in intersection
            continue
        succ = succ + 1
    if succ not in successors:
        successors[succ] = [e]
    else:
        new_succs = successors[succ]
        new_succs.append(e)
        successors[succ] = new_succs

# here identify which transitions are left turns, which ones right turns, and which ones go straight.
for e in G.edges(data=True):
    for f in G.edges([e[1]], data=True):
        # going from e to f
        angle = (e[2]["data"]["final_bearing"] - f[2]["data"]["init_bearing"] + 180) % 360
        if angle < 135: # right turn
            pass
        elif angle >= 135 and angle < 225: # straight
            pass
        else: # left
            pass

#########################################################################
# Hazards and Reactions                                                 #
#########################################################################

hazards = ['pedestrian', 'obstacle', 'jam']
# control ocurrence probability, one for each hazard, the larger, the more likely
alphas = {'roadblock': 0.002,
          'pedestrian': 0.05,
          'obstacle': 0.02,
          'jam': 0.1} 
reactions = {'roadblock': ['uturn'],
             'pedestrian' : ['brake', 'honk', 'changelane'],
             'obstacle': ['changelane', 'uturn'],
             'jam': ['honk', 'uturn']}
accident_prob = {('pedestrian', 'brake')      : 0.01,
                 ('pedestrian', 'honk')       : 0.04,
                 ('pedestrian', 'changelane') : 0.03,
                 ('obstacle', 'changelane')   : 0.02,
                 ('obstacle', 'uturn')        : 0.02,
                 ('jam', 'honk')              : 0.01,
                 ('jam', 'uturn')             : 0.02 }

# all possible reactions
reacts = set([])
for r in reactions.itervalues():
    reacts = reacts | set(r)
reacts = list(reacts)

def powerset(S):
    if len(S) <= 1:
        yield S
        yield []
    else:
        for s in powerset(S[1:]):
            yield [S[0]]+s
            yield s
# all possible combinations of hazards
haz = []
for i1 in xrange(0,len(hazards)):
    haz.append([hazards[i1]])
    for i2 in xrange(i1+1, len(hazards)):
        haz.append([hazards[i1],hazards[i2]])
haz.append([])

#########################################################################
# Output PRISM Model                                                    #
#########################################################################

# one field per edge
N = len(G.edges())
M = 3*(len(haz)-1) - 2*len(hazards)

# build a dictionary of edges and the associated ids
# do the same in reverse as well
edgeid = {}
idedge = {}
# also record lanes of each road
lanes = {}
i = 0
for e in G.edges(data=True):
    edgeid[(e[0],e[1])] = i
    idedge[i] = (e[0],e[1])
    if "lanes" in e[2]["data"]:
        lanes[i] = e[2]["data"]["lanes"]
    else:
        lanes[i] = 1
    i = i + 1

smg.write("// PRISM model for autonomous car case study\n")
smg.write("// Two player stochastic game with three objectives.\n")
smg.write("//\n")
smg.write("// This is a machine-generated file\n")
smg.write("// Author: Clemens Wiltsche\n")
today = datetime.date.today()
smg.write("// Generated: %s\n" % (today.strftime('%d/%m/%Y')))
smg.write("\n\n")
smg.write("smg\n\n")

# initialize player actions
smg.write("player p1\n")
first = True
for e in G.edges():
    i = edgeid[e]
    if(not first):
        smg.write(",\n")
    first = False
    smg.write("\t")
    # reactions
    for r in reacts:
        smg.write("[%s_%i], " % (r, i))
    # steering
    j = 0
    first_ = True
    onlyuturn = True
    for f in G.edges([e[1]]):
        if (f[1],f[0]) == (e[0],e[1]): # no uturn in intersection
            continue
        onlyuturn = False
        if(not first_):
            smg.write(", ")
        first_ = False
        smg.write("[move_%i_%i]" % (j, i))
        j = j + 1
    # alternatively, termination
    if onlyuturn or len(G.edges([e[1]]))==0:
        smg.write("[term_%i]" % (i))
smg.write("\nendplayer\n\n")

smg.write("player p2\n")
first = True
for e in G.edges():
    i = edgeid[e]
    if(not first):
        smg.write(",\n")
    first = False
    smg.write("\t")
    # hazard
    for h in hazards:
        smg.write("[%s_%i], " % (h, i))
    smg.write("[hazard_%i]" % (i))
smg.write(",\n\t[term]")
smg.write("\nendplayer\n\n")


# positions in topology
smg.write("const int POS_init = %i;\n" % (init))
smg.write("const int POS_goal = %i;\n" % (N))
smg.write("const int POS_term = %i;\n" % (N+1))
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
        onlyone = False
        if len(G.edges([e[1]]))==1:
            onlyone = True
        for f in G.edges([e[1]], data=True):
            if ie != goal and (f[1],f[0]) == (e[0],e[1]): # no uturn in intersection
                if onlyone:
                    smg.write("const int DEST_%i = POS_term;\n" % (ie))
                continue
            i_f = edgeid[(f[0],f[1])]
            if ie==goal and onlyone:
                smg.write("const int DEST_%i = POS_goal;\n" % (ie))
            elif ie==goal:
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
smg.write("\n")

# returns a distribution of hazard probabilities given the length of a road
def distr(length):
    # how many hazard combinations are there?
    dist = {} # nothing in distribution yet
    cumulative = 0.0
    empty_i = 0
    # fill distribution - one entry for each hazard combination
    for i in xrange(0,len(haz)):
        hs = haz[i]
        if hs==[]:
            empty_i = i
        else:
            Alpha = 1
            for h in hs:
                Alpha = Alpha*alphas[h]
            dist[i] = math.tanh(Alpha*length)/len(haz)
            cumulative = cumulative + dist[i]
    # probability of nothing happening
    dist[empty_i] = 1 - cumulative
    # return distribution
    return dist

# write hazard probabilities for each edge
for e in G.edges(data=True):
    ie = edgeid[(e[0],e[1])]
    d = distr(e[2]["data"]["dist"])
    for j in xrange(0, len(haz)):
        smg.write("const double DISTR_%i_%i = %f;\n" % (j, ie, d[j]))
smg.write("\n")

for e in G.edges(data=True):
    ie = edgeid[(e[0],e[1])]
    if (e[1],e[0]) in edgeid or lanes[ie]>1: # changelane possible if not single-lane oneway
        smg.write("const int CHANGELANE_0_%i = -1;\n" % (ie)) # accident first option
        smg.write("const int CHANGELANE_1_%i = -2;\n" % (ie)) # changing is second option
    else: # changelane illegal - go to violation state in any case
        smg.write("const int CHANGELANE_0_%i = -3;\n" % (ie))
        smg.write("const int CHANGELANE_1_%i = -3;\n" % (ie))

# write uturn reactions for each edge
for e in G.edges(data=True):
    ie = edgeid[(e[0],e[1])]
    if (e[1],e[0]) in edgeid: # uturn possible if not one-way (this also means highways)
        revedgeid = edgeid[(e[1],e[0])]
        smg.write("const int UTURN_0_%i = -1;\n" % (ie)) # accident first option
        smg.write("const int UTURN_1_%i = -2;\n" % (ie)) # turning is second option
        smg.write("const int UTURNPLAYER_%i = 1;\n" % (ie)) # if in -2, need player 1
        smg.write("const int REVEDGE_%i = %i;\n" % (ie, revedgeid)) # reverse edge exists
    else: # uturn illegal - go to violation state in any case
        smg.write("const int UTURN_0_%i = -3;\n" % (ie))
        smg.write("const int UTURN_1_%i = -3;\n" % (ie))
        smg.write("const int UTURNPLAYER_%i = 2;\n" % (ie)) # if in -3, need player 2
        smg.write("const int REVEDGE_%i = %i;\n" % (ie, ie)) # reverse edge doesn't exist but is irrelevant

smg.write("\nglobal p : [1..2] init 2;\n") # players
smg.write("global car_position : [0..%i] init POS_init;\n" % (N+2))
# special states:
# -1: accident (terminal)
# -2: sink - goto next edge
# -3: traffic violation (terminal)
# -4: braking
# -5: uturning
# -11: game terminated
smg.write("global s : [-11..%i] init 0;\n" % (M))

def subhazard(from_s, to_s, subhaz, k):
    if len(subhaz)==1:
        carreact(from_s, subhaz[0], k)
    else:
        for i in xrange(0, len(subhaz)):
            smg.write("\t[%s_%i] p=2 & s=%i & car_position=POS_%i-> (p'=1) & (s'=%i);\n" % (subhaz[i], k, from_s, k, to_s[i]))

def carreact(from_s, hazz, k):
    for i in xrange(0, len(reactions[hazz])):
        if (reactions[hazz][i]=='brake'): # brake action is special (need to insert state to let time pass)
            smg.write("\t[brake_%i] p=1 & s=%i & car_position=POS_%i -> %f : (p'=2) & (s'=-1) + %f : (p'=1) & (s'=-2);\n" % (k, from_s, k, accident_prob[(hazz,'brake')], 1.0-accident_prob[(hazz,'brake')])) 
        elif(reactions[hazz][i]=='uturn'): # uturn action is special
            smg.write("\t[uturn_%i] p=1 & s=%i & car_position=POS_%i -> %f : (p'=2) & (s'=UTURN_0_%i) + %f : (p'=UTURNPLAYER_%i) & (s'=UTURN_1_%i) & (car_position'=REVEDGE_%i);\n" % (k, from_s, k, accident_prob[(hazz,'uturn')], k, 1.0-accident_prob[(hazz,'uturn')], k, k, k))
        elif(reactions[hazz][i]=='changelane'): # changelane action is special
            smg.write("\t[changelane_%i] p=1 & s=%i & car_position=POS_%i -> %f : (p'=2) & (s'=CHANGELANE_0_%i) + %f : (p'=1) & (s'=CHANGELANE_1_%i);\n" % (k, from_s, k, accident_prob[(hazz,'changelane')], k, 1.0-accident_prob[(hazz,'changelane')], k)) 
        else: # standard actions
            smg.write("\t[%s_%i] p=1 & s=%i & car_position=POS_%i -> %f : (p'=2) & (s'=-1) + %f : (p'=1) & (s'=-2);\n" % (reactions[hazz][i], k, from_s, k, accident_prob[(hazz,reactions[hazz][i])], 1.0-accident_prob[(hazz,reactions[hazz][i])])) 

# the base cases - one for each number of successors (taking no uturn at intersection into account)
for i, succs in successors.iteritems():
    # k is the first edge with i successors
    k = edgeid[(succs[0][0],succs[0][1])]
    # the corresponding distribution
    smg.write("\n\nmodule field_%i\n\n" % (k))
    smg.write("\t[hazard_%i] p=2 & s=0 & car_position=POS_%i -> DISTR_0_%i+DISTR_1_%i+DISTR_2_%i+DISTR_5_%i : (s'=-6) + DISTR_3_%i+DISTR_4_%i+DISTR_6_%i : (s'=-7);\n" % (k, k, k, k, k, k, k, k, k))

    # jam or pedestrian
    smg.write("\t[hazard_%i] p=2 & s=-6 & car_position=POS_%i -> (DISTR_2_%i+DISTR_5_%i)/(DISTR_0_%i+DISTR_1_%i+DISTR_2_%i+DISTR_5_%i) : (s'=-8) + (DISTR_0_%i+DISTR_1_%i)/(DISTR_0_%i+DISTR_1_%i+DISTR_2_%i+DISTR_5_%i) : (s'=-9);\n" % (k, k, k, k, k, k, k, k, k, k, k, k, k, k))
    # obstacle or nothing
    smg.write("\t[hazard_%i] p=2 & s=-7 & car_position=POS_%i -> (DISTR_3_%i+DISTR_4_%i)/(DISTR_3_%i+DISTR_4_%i+DISTR_6_%i) : (s'=-10) + (DISTR_6_%i)/(DISTR_3_%i+DISTR_4_%i+DISTR_6_%i) : (s'=-2) & (p'=1);\n" % (k, k, k, k, k, k, k, k, k, k, k))

    # jam&pedestrian or jam only
    pj = 0
    if ['pedestrian','jam'] in haz:
        pj = haz.index(['pedestrian','jam'])
    else:
        pj = haz.index(['jam','pedestrian'])
    j = haz.index(['jam'])
    smg.write("\t[hazard_%i] p=2 & s=-8 & car_position=POS_%i -> DISTR_2_%i/(DISTR_2_%i+DISTR_5_%i) : (s'=%i) + DISTR_5_%i/(DISTR_2_%i+DISTR_5_%i) : (s'=%i) & (p'=1);\n" % (k, k, k, k, k, pj+1, k, k, k, j+1))

    # pedestrian&obstacle or pedestrian only
    po = 0
    if ['pedestrian','obstacle'] in haz:
        po = haz.index(['pedestrian','obstacle'])
    else:
        po = haz.index(['obstacle','pedestrian'])
    p = haz.index(['pedestrian'])
    smg.write("\t[hazard_%i] p=2 & s=-9 & car_position=POS_%i -> DISTR_1_%i/(DISTR_0_%i+DISTR_1_%i) : (s'=%i) + DISTR_0_%i/(DISTR_0_%i+DISTR_1_%i) : (s'=%i) & (p'=1);\n" % (k, k, k, k, k, po+1, k, k, k, p+1))

    # obstacle&jam or obstacle only
    oj = 0
    if ['jam','obstacle'] in haz:
        oj = haz.index(['jam','obstacle'])
    else:
        oj = haz.index(['obstacle','jam'])
    o = haz.index(['obstacle'])
    smg.write("\t[hazard_%i] p=2 & s=-10 & car_position=POS_%i -> DISTR_4_%i/(DISTR_3_%i+DISTR_4_%i) : (s'=%i) + DISTR_3_%i/(DISTR_3_%i+DISTR_4_%i) : (s'=%i) & (p'=1);\n" % (k, k, k, k, k, oj+1, k, k, k, o+1))

    # instantiate hazard and pick reaction - bad guy followed by good guy
    for j in xrange(0, len(haz)-1):
        if len(haz[j]) == 1:
            subhazard(j+1, [], haz[j], k)
        if len(haz[j]) == 2:
            subhazard(j+1, [haz.index([haz[j][0]])+1, haz.index([haz[j][1]])+1], haz[j], k)

    # pick destination - good guy
    for j in xrange(0,i):
        smg.write("\t[move_%i_%i] p=1 & s=%i & car_position=POS_%i-> (car_position'=DEST_%i_%i) & (s'=0) & (p'=2);\n" % (j, k, -2, k, j, k))

    if i == 0:
        smg.write("\t[term_%i] p=1 & s=%i & car_position=POS_%i -> (car_position'=DEST_%i) & (s'=0) & (p'=2);\n" % (k, -2, k, k))
    
    smg.write("\nendmodule\n\n\n")

# iterate through topology
for i, succs in successors.iteritems():
    # if just one edge with i successors, then have it covered already
    if len(succs) > 1:
        l = edgeid[(succs[0][0],succs[0][1])]
        for e in succs[1:]: # start iterating at seccond edge
            k = edgeid[(e[0],e[1])]
            # module and position
            smg.write("module field_%i = field_%i [POS_%i=POS_%i" % (k, l, l, k))
            # hazard
            for h in hazards:
                smg.write(", %s_%i=%s_%i" % (h, l, h, k))
            smg.write(", hazard_%i=hazard_%i" % (l, k))
            # reactions
            for r in reacts:
                smg.write(", %s_%i=%s_%i" % (r, l, r, k))
            # steering
            for j in xrange(0,i):
                smg.write(", move_%i_%i=move_%i_%i, DEST_%i_%i = DEST_%i_%i" % (j, l, j, k, j, l, j, k))
            # hazard-set distribution
            for j in xrange(0,len(haz)):
                smg.write(", DISTR_%i_%i=DISTR_%i_%i" % (j, l, j, k))
            # uturn reactions
            smg.write(", UTURN_0_%i=UTURN_0_%i, UTURN_1_%i=UTURN_1_%i, REVEDGE_%i=REVEDGE_%i, UTURNPLAYER_%i=UTURNPLAYER_%i" % (l, k, l, k, l, k, l, k))
            # changelane reactions
            smg.write(", CHANGELANE_0_%i=CHANGELANE_0_%i, CHANGELANE_1_%i=CHANGELANE_1_%i" % (l, k, l, k))
            # alternatively, termination
            if i == 0:
                smg.write(", term_%i=term_%i, DEST_%i=DEST_%i" % (l, k, l, k))
            smg.write("] endmodule\n")

# terminal states
smg.write("\nmodule terminals\n")
smg.write("\t[term] s=-1 -> (car_position'=POS_term) & (p'=2) & (s'=-11);\n") # accident
smg.write("\t[term] s=-3 -> (car_position'=POS_term) & (p'=2) & (s'=-11);\n") # illegal action
smg.write("\t[term] car_position=POS_goal -> (car_position'=POS_term) & (p'=2) & (s'=-11);\n")
smg.write("\t[term] car_position=POS_term -> (car_position'=POS_term) & (p'=2) & (s'=-11);\n")
smg.write("endmodule\n\n")

# REWARDS

# reaching goal
smg.write("\nrewards \"reach_goal\"\n")
smg.write("\tcar_position=POS_goal : 1;\n")
smg.write("endrewards\n\n")

# avoiding accidents
smg.write("\nrewards \"accidents\"\n")
smg.write("\ts=-1 : 1;\n")
smg.write("endrewards\n\n")

# road value
smg.write("\nrewards \"road_quality\"\n")
for e in G.edges(data=True):
    ie = edgeid[(e[0],e[1])]
    smg.write("\tcar_position=%i & s=-2: %f;\n" % (ie, e[2]["data"]["value"]*e[2]["data"]["dist"]/1000))
smg.write("endrewards\n\n")

#########################################################################
# Properties File                                                       #
#########################################################################

# write the property - a three-objective conjunctive goal
if(len(sys.argv)>=2):
    if(int(sys.argv[1])==1): # islip
        prop.write("<<1>> (R{\"reach_goal\"}>=0.7 [ C ] &  R{\"accidents\"}<=0.3 [ C ] &  R{\"road_quality\"}>=6.0 [ C ])");
    elif(int(sys.argv[1])==2): # islip
        prop.write("<<1>> (R{\"reach_goal\"}>=0.7 [ C ] & R{\"accidents\"}<=0.3 [ C ] & R{\"road_quality\"}>=6.0 [ C ])");
else: # islip is default
    prop.write("<<1>> (R{\"reach_goal\"}>=0.7 [ C ] & R{\"accidents\"}<=0.3 [ C ] & R{\"road_quality\"}>=6.0 [ C ])");

#########################################################################
# Bash File                                                             #
#########################################################################

# write the command to execute
bash.write("#!/bin/bash\n\n")
bash.write("../../bin/prism %s{.prism,.props} -prop 1 -multirounding -multimaxciter 500 -baselineaccuracy 200 -increasefactor 1.01 -paretoepsilon 0.001 -logcpareto -gs -exportstrat %s.strat 2>&1 1> %s.log &\n" % (filename, filename, filename))

#########################################################################
# Wrap up                                                               #
#########################################################################

# close down files - flush
smg.close()
prop.close()
bash.close()
