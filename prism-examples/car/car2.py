#!/usr/bin/python -O

from __future__ import division
import datetime, math, sys
import osm2graph

smg = open("car2.gen.smg", "w");
prop = open("car2.gen.props", "w");

mapfile = "map(4).osm"
if(len(sys.argv)>=2):
    mapfile = sys.argv[1];

#########################################################################
# Prepare Graph                                                         #
#########################################################################

G = osm2graph.read_osm(mapfile)

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
    if not e[2]["data"]["oneway"] and not G.has_edge(e[1],e[0]):
        G.add_edge(e[1],e[0],attr_dict=e[2])


#########################################################################
# Hazards and Reactions                                                 #
#########################################################################

#hazards = ['tree', 'person', 'debris', 'pothole', 'aliens']
hazards = ['tree', 'person']
alphas = {'tree': 0.002, 'person': 0.05} # control ocurrence probability, one for each hazard, the larger, the more likely
maxVelocity = 2;
#reactions = {'tree' : ['dodge', 'brake'], 'person' : ['brake', 'honk', 'dodge'], 'debris': ['dodge', 'brake'], 'pothole': ['dodge', 'brake'], 'aliens': ['honk', 'lazergun']}
reactions = {'tree' : ['dodge', 'brake'], 'person' : ['brake', 'honk', 'dodge']}
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
haz = list(powerset(hazards))

#########################################################################
# Output PRISM Model                                                    #
#########################################################################

# one field per edge
N = len(G.edges())
n = len(hazards)
M = 1 + 2**n + n*2**(n-1)

# TODO: properly identify goal
goal = 104 # hilltop gardens for now

# build a dictionary of edges and the associated ids
edgeid = {}
i = 0
for e in G.edges():
    edgeid[e] = i
    i = i + 1

smg.write("// PRISM model for case study of car driving\n")
smg.write("// - The game consists of two players and will have three goals.\n")
smg.write("//\n")
smg.write("// This is a machine-generated file\n")
smg.write("// Author: Clemens Wiltsche\n")
today = datetime.date.today()
smg.write("// Last modified: %s\n" % (today.strftime('%d/%m/%Y')))

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
    # velocity
    for v in xrange(1,maxVelocity+1):
        smg.write("[velocity_%i_%i], " % (v, i))
    # reactions
    for r in reacts:
        smg.write("[%s_%i], " % (r, i))
    # steering
    for j in xrange(0, len(G.edges([e[1]]))):
        if not j==0:
            smg.write(", ")
        smg.write("[move_%i_%i]" % (j, i))
    # alternatively, termination
    if len(G.edges([e[1]]))==0:
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
    smg.write("[none_%i], [hazard_%i]" % (i, i))
smg.write(",\n\t[term]")
smg.write("\nendplayer\n\n")


# positions in topology
smg.write("const int POS_init = %i;\n" % (0))
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
smg.write("\n")

# returns a distribution of hazard probabilities given the length of a road
def distr(length):
    # how many hazard combinations are there?
    numhaz = n*2**(n-1) # adding all cardinalities of elements in powerset of hazards
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
            dist[i] = math.tanh(Alpha*length)/numhaz
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





# identify which edge has how many successors
successors = {}
for e in G.edges(data=True):
    succ = len(G.edges([e[1]]))
    if succ not in successors:
        successors[succ] = [e]
    else:
        new_succs = successors[succ]
        new_succs.append(e)
        successors[succ] = new_succs

smg.write("\nglobal p : [1..2] init 1;\n") # players
smg.write("global car_position : [0..%i] init POS_init;\n" % (N+2))
smg.write("global s : [-2..%i] init 0;\n\n" % (M))
smg.write("global car_velocity : [1..%i] init 1;\n" % (maxVelocity))


def subhazard(from_s, to_s, subhaz, k):
    if len(subhaz)==0:
        smg.write("\t[none_%i] p=2 & s=%i & car_position=POS_%i -> (p'=1) & (s'=%i);\n" % (k, from_s, k, -2))
    new_s = to_s
    for i in xrange(0, len(subhaz)):
        smg.write("\t[%s_%i] p=2 & s=%i & car_position=POS_%i-> (p'=1) & (s'=%i);\n" % (subhaz[i], k, from_s, k, new_s))
        carreact(new_s, subhaz[i], k)
        new_s = new_s + 1
    return new_s

def carreact(from_s, hazz, k):
    # TODO: actually insert probability here - possibly stored in reaction dict
    acc_prob_num = 1
    acc_prob_den = 3
    for i in xrange(0, len(reactions[hazz])):
        smg.write("\t[%s_%i] p=1 & s=%i -> %i/%i : (p'=2) & (s'=%i) + %i/%i : (p'=1) & (s'=%i);\n" % (reactions[hazz][i], k, from_s, acc_prob_num, acc_prob_den, -1, acc_prob_den-acc_prob_num, acc_prob_den, -2)) 


# the base cases - one for each number of successors
for i, succs in successors.iteritems():

    # k is the first edge with i successors
    k = edgeid[(succs[0][0],succs[0][1])]
    # the corresponding distribution

    smg.write("module field_%i\n\n" % (k))

    # pick speed - good guy
    for v in xrange(1,maxVelocity+1):
        smg.write("\t[velocity_%i_%i] p=1 & s=0 & car_position=POS_%i-> (car_velocity'=%i) & (s'=1) & (p'=2);\n" % (v, k, k, v))

    # pick available hazards - stochastic guy
    smg.write("\t[hazard_%i] p=2 & s=1 & car_position=POS_%i-> " % (k, k))
    init = True
    for j in xrange(0,len(haz)):
        if not init:
            smg.write(" + ")
        init = False
        smg.write("DISTR_%i_%i : (s'=%i) & (p'=2)" % (j, k, 1+j+1))
    smg.write(";\n")

    # instantiate hazard - bad guy
    new_s = 1+len(haz)+1
    for j in xrange(0, len(haz)):
        new_s = subhazard(1+j+1, new_s, haz[j], k)

    # pick destination - good guy
    for j in xrange(0,i):
        smg.write("\t[move_%i_%i] p=1 & s=%i & car_position=POS_%i-> (car_position'=DEST_%i_%i) & (s'=0) & (p'=1);\n" % (j, k, -2, k, j, k))

    if i == 0:
        smg.write("\t[term_%i] p=1 & s=%i & car_position=POS_%i -> (car_position'=DEST_%i) & (s'=0) & (p'=2);\n" % (k, -2, k, k))
    
    smg.write("endmodule\n\n")

# iterate through topology
for i, succs in successors.iteritems():
    # if just one edge with i successors, then have it covered already
    if len(succs) > 1:
        l = edgeid[(succs[0][0],succs[0][1])]
        for e in succs[1:]: # start iterating at seccond edge
            k = edgeid[(e[0],e[1])]
            # module and position
            smg.write("module field_%i = field_%i [POS_%i=POS_%i" % (k, l, l, k))
            # velocity
            for v in xrange(1,maxVelocity+1):
                smg.write(", velocity_%i_%i=velocity_%i_%i" % (v, l, v, k))
            # hazard
            for h in hazards:
                smg.write(", %s_%i=%s_%i" % (h, l, h, k))
            smg.write(", none_%i=none_%i, hazard_%i=hazard_%i" % (l, k, l, k))
            # reactions
            for r in reacts:
                smg.write(", %s_%i=%s_%i" % (r, l, r, k))
            # steering
            for j in xrange(0,i):
                smg.write(", move_%i_%i=move_%i_%i, DEST_%i_%i = DEST_%i_%i" % (j, l, j, k, j, l, j, k))
            # hazard-set distribution
            for j in xrange(0,len(haz)):
                smg.write(", DISTR_%i_%i=DISTR_%i_%i" % (j, l, j, k))
            # alternatively, termination
            if i == 0:
                smg.write(", term_%i=term_%i, DEST_%i=DEST_%i" % (l, k, l, k))
            smg.write("] endmodule\n")

# terminal states
smg.write("\nmodule terminals\n")
smg.write("\t[term] s=-1 -> (s'=-1);\n") # accident
smg.write("\t[term] car_position=POS_goal -> (car_position'=POS_goal);\n")
smg.write("\t[term] car_position=POS_term -> (car_position'=POS_term);\n")
smg.write("endmodule\n\n")

# properties
smg.write("formula goal1 = (car_position=POS_goal);\n")
smg.write("formula goal2 = (s==-1);\n")


'''
# rewards
smg.write("rewards\n")
smg.write("\taccident & !car_slow : 1000;\n")
smg.write("\taccident & car_slow : 50;\n")
smg.write("\tstall : 2;\n")
smg.write("endrewards\n\n")
'''


#########################################################################
# Properties                                                            #
#########################################################################


# write the property - a single two-objective goal
prop.write("<<1>> P>=0.5 [ goal1 U goal2 ]");

# close down files - flush
smg.close()
prop.close()
