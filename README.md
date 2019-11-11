# Human Online Adaptation to Changes in Prior Probability

This repository accompanies the article by Norton et al. (2019) [[1](https://github.com/lacerbi/ChangeProb/blob/master/README.md#reference)].
It includes human subjects' data and the code used for fitting and comparing the models reported in the paper.

**Data:** The data are stored in the `./data` subfolder, in separate `.mat` files for each of the `N = 11` subjects used in the paper. Each file contains a `data` struct, whose fields contain information about the sessions. Each struct contains information for both the covert-criterion task and for the overt-criterion task (800 trials each).

**Models:** All model files are stored in the `./matlab` subfolder. The main file to run a model fit is `changeprob_runfit.m`. See the [function description](https://github.com/lacerbi/ChangeProb/blob/master/matlab/changeprob_runfit.m) for more information about the inputs and outputs. An example run on one dataset and model (identified by model and data id `1`) can be obtained with:
```
> changepoint_runfit(1,'marglike');
```
Example scripts to run the model fits on a computer cluster (using the [Slurm](https://slurm.schedmd.com/) workload manager) are in the `./scripts` subfolder.


## Reference

1. Norton EH, Acerbi L, Ma WJ, Landy MS (2019) Human online adaptation to changes in prior probability. *PLoS Comput Biol* 15(7): e1006681. https://doi.org/10.1371/journal.pcbi.1006681 ([link](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006681))

### License

This code is released under the terms of the [MIT License](https://github.com/lacerbi/ChangeProb/blob/master/LICENSE.txt).
