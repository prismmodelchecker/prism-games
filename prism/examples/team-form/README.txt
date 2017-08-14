The folder contains the SMG models and properties from []
The file team-form-X-Y.prism conains a model for X\in{online, offline}
algorithms and agents arranged in  Y\in{fc,ia,s,r}  topology,
where fc=fully connected, ia=isolated agent, s=star, r=ring.
the properties for the models are in corresponding team-form-Y.props
file.
Script run-experiments.sh runs model checking for all models and properties 
and outputs the results.