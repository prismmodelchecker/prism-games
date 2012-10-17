# Futures market investor example
# Based on: 
#   A. McIver and C. Morgan
#   Results on the quantitative mu-calculus qMu
#   ACM Transactions on Computational Logic, 8(1), 2007


# One investor (two player-game)
prism-games investor.prism investor.props -const vmax=10,vinit=5

# Two investors (three-player game)
prism-games two_investors.prism two_investors.props

# Hard-coded "seat of the pants" strategy (MDP)
prism-games investor_sotp.prism investor_sotp.props -const vmax=10,vinit=5

# One investor, with infinite loop removed, allowing us to test F rather than Fc
prism-games investor_mod.prism investor_mod.props -const vmax=10,vinit=5

# One investor, with conversion to min, allowing us to test F rather than Fc (but breaks because of ECs?)
prism-games investor_min.nm -pctl '<<1>> Rmin=? [ F i=2 ]' -const vmax=10,vinit=5


# Experiments: change: -const vmax=10,vinit=5
# to: -const vmax=10,vinit=0:10 -exportresults stdout
# to get results for a range of initial values