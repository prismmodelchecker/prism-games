This case study is concerned with controlling the temperature in multiple rooms in a building. The specific problem formulation is inspired by [FI04] and [SGA15]. The model is formulated as stochastic differential equations, and we use formal discretisation methods similar to those in [SGA15] to obtain a stochastic game model in PRISM-games.

The model is in temperature.prism, and the properties are in temperature.props. The discretisation can be re-run using the MATLAB script temperature.m, where also the model parameters can be changed.


For more information, see: http://www.prismmodelchecker.org/games/examples.php

=====================================================================================

[FI04]
Fehnker, A. and Ivancic, F.
Benchmarks for hybrid systems verification.
In Proc. HSCC'04, volume 2993 of LNCS, pages 326-341, Springer. 2004.

[SGA15]
Soudjani, S.E.Z. and Gevaerts, C. and Abate, A.
FAUST2: Formal Abstractions of Uncountable-STate STochastic Processes
In Proc. TACAS'15, volume 9035 of LNCS, pages 272-286, Springer. 2015.
