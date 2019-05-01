This case study is based on a distributed team formation protocol
proposed in [GdJ05] and modelled/verified in [CKPS11].
It was also used as a stochastic game benchmark in [CFK+13b].

The folder contains several variants of the models and properties.
The models are in files team-form-X-Y.prism where X is the algorithm type
(X = online,offline) and Y is the topology of the agents (Y=fc,ia,s,r), where:
fc=fully connected, ia=isolated agent, s=star, r=ring.
The properties for the models are in corresponding file team-form-Y.props.

=====================================================================================

[CFK+13b]
Taolue Chen, Vojtěch Forejt, Marta Kwiatkowska, David Parker and Aistis Simaitis.
Automatic Verification of Competitive Stochastic Systems.
Formal Methods in System Design, 43(1), pages 61-92, Springer. August 2013. [pdf] [bib]

[CKPS11]
Taolue Chen, Marta Kwiatkowska, David Parker and Aistis Simaitis.
Verifying Team Formation Protocols with Probabilistic Model Checking.
In Proc. 12th International Workshop on Computational Logic in Multi-Agent Systems (CLIMA XII 2011), volume 6814 of LNCS, pages 190-297, Springer, 2011.

[GdJ05]
M. Gaston and M. desJardins.
Agent-organized networks for dynamic team formation.
In Proc. 4th International Joint Conference on Autonomous Agents and Multi-Agent Systems (AAMAS’05), pages 230-237, ACM, 2005.

