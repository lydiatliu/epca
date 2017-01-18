# *e*PCA
Software implementation for *e*PCA and denoising as described in the [paper](http://arxiv.org/abs/1611.05550):
* Lydia T. Liu, Edgar Dobriban and Amit Singer. ***e*****PCA: High Dimensional Exponential Family PCA.**  *arXiv preprint arXiv:1611.05550*, 2016. 


## Contents
* ```software/``` : software for applying ePCA on one's own datasets. In particular, ```exp_fam_pca.m``` implements ePCA and ```wiener_filter.m``` implements the generalized wiener filter (or EBLP) for denoising, as introduced in the aforementioned paper.
* ```experiments/``` : scripts for reproducing experimental results in the aforementioned paper. Large datasets are excluded.

## Requirements
MATLAB. No other downloads are required.

## Acknowledgements
This implementation uses ```standard_spiked_forward.m``` from [EigenEdge](https://github.com/dobriban/EigenEdge).

