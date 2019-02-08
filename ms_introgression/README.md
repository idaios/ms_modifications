# ms_modifications
A version of ms with the following modifications

NOTE: This version of ms works for cases of only 1 individual per population and no recombination (if the option -x is used). 

1. add an argument -x <POPULATION> <TIME>
This argument will force sampling of population <POPULATION> at time <TIME>.  In fact, what it will do, is to avoid mutations if they are more recent than TIME in population POPULATION. Since, however, coalescent and recombination events may also happen (and this will be wrong in our case), I restrict its usage for cases when the sample size per population is 1 and there is no recombination.


2. add an argument -ev <TIME> popi popj
At time TIME a line will migrate (backwards in time) from popi to popj. This argument is different than usual migration (when a rate is given) because it forces a migration event at a certain time point.


