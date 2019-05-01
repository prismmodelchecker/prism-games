Microgrid Demand-Side Management

This case study is based on the distributed energy management algorithm for smart grids
by Hildmann & Saffre [HS11], modelled and verified in [CFK+13b]. The files here
correspond to a model of the protocol where households are allowed to "cheat"
(deviate from the protocol). In [CFK+13b], a fix is proposed where households
are punished for deviating. These models are the "no-punishment" variant.

File name explanation:
* mdsm2304.{prism,props} - 2 households, 3 day time horizon, 0 deterministic households, 4 jobs

For the full set of files, see: http://www.prismmodelchecker.org/files/fmsd-smg/

For more information, see: https://www.prismmodelchecker.org/casestudies/mdsm.php

=====================================================================================

[CFK+13b]
Taolue Chen, Vojtěch Forejt, Marta Kwiatkowska, David Parker and Aistis Simaitis.
Automatic Verification of Competitive Stochastic Systems.
Formal Methods in System Design, 43(1), pages 61-92, Springer, 2013.

[HS11]
H. Hildmann and F. Saffre.
Influence of variable supply and load flexibility on Demand-Side Management.
In Proc. 8th International Conference on the European Energy Market (EEM'11), pages 63-68, 2011.
