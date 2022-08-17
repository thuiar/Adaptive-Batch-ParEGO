# Adaptive Batch-ParEGO
A paralled version of ParEGO algorithm for expensive multi-objective problems (MOPs).

## Introduction
 This repository contains Matlab implementation of the algorithm framework for adaptive batch-ParEGO in the research paper [An Adaptive Batch Bayesian Optimization Approach for Expensive Multi-Objective Optimization](https://www.sciencedirect.com/science/article/pii/S0020025522009057) (**Accepted by [Information Sciences 2022](https://www.sciencedirect.com/journal/information-sciences)**).

## Algorithm Preparation
We use [PlatEMO-V2.9.0](https://github.com/BIMK/PlatEMO/releases/tag/PlatEMO_v2.9.0), an evolutionary multi-objective optimization platform, to implement all the related experiments. Details on how to use PlatEMO can be found in [manual.pdf](https://github.com/BIMK/PlatEMO/blob/master/PlatEMO/manual.pdf). Before starting our methods, we recommend to  carefully study how to use the PlatEMO platform.

## Get Started
We recommend runing all the related experiments with GUI of PlatEMO. To invoke the interface, use the function:
```
main()
```
More details can be found in [manual.pdf](https://github.com/BIMK/PlatEMO/blob/master/PlatEMO/manual.pdf). 

## Baselines
The baseline methods in Adaptive Batch-ParEGO include ParEGO [[1]](#parego), MOEA/D-EGO [[2]](#moeadego), ReMO [[3]](#remo), Multi-LineBO, SparseEA [[4]](#sparseea) and MOEA/PSL [[5]](#moeapsl). ReMO is an optimization architecture that can be equiped with any well known derivative-free MO algorithm. We equip ReMO with ParEGO in this paper to make comparisons with EGO-based methods. Multi-LineBO is a version of single-objective LineBO [[6]](#linebo).

|    Algorithm Name      | Characteristics|    Published     |
| :---------: | :-----------------------------------------------------: | :------------------------------------------------------------------------------------: | 
[ParEGO](https://www.cs.bham.ac.uk/~jdk/parego/) | Multi-objective, Sequential |        IEEE Transactions on Evolutionary Computation 2006         |  
[MOEA/D-EGO](https://ieeexplore.ieee.org/document/5353656) | Multi-objective, Batch |        IEEE Transactions on Evolutionary Computation 2010         | 
[ReMO](https://ojs.aaai.org/index.php/AAAI/article/view/10664) | Multi-objective, Sequential |        AAAI 2017         |   
[Multi-LineBO](http://proceedings.mlr.press/v97/kirschner19a/kirschner19a.pdf) | Multi-objective, Sequential |        ICML 2019         |
[SparseEA](https://ieeexplore.ieee.org/document/8720021) | Multi-objective, Large-scale |IEEE Transactions on Evolutionary Computation 2019 |  
[MOEA/PSL](https://ieeexplore.ieee.org/document/9047876) | Multi-objective, Large-scale |IEEE Transactions on Cybernetics 2020 | 

## Benchmark Problems
Benchmark problems contain six three-objective benchmark problems taken from the DTLZ test suite [[7]](#dtlz), seven two-objective benchmark problems from the UF test suite [[8]](#uf), nine three-objective benchmark problems from WFG test suite [[9]](#wfg) and a real-world hyper-parameter tuning of neural network task [[10]](#nn). The source code of hyper-parameter tuning task of neural networks can be found [here](https://github.com/rasmusbergpalm/DeepLearnToolbox).
|    Problems      |                   M| D                           | Characteristics
| :---------: | :-----------------------------------------------------: | :------------------------------------------------------------------------------------: | :---------------: |
[DTLZ](https://www.cs.bham.ac.uk/~jdk/parego/) | 3 |  10  |  DTLZ11, DTLZ2, DTLZ3,DTLZ5, DTLZ6, DTLZ7|
[WFG](https://ieeexplore.ieee.org/document/5353656) | 3 | 10  | WFG1-9
[UF](https://ojs.aaai.org/index.php/AAAI/article/view/10664) | 2 | 10| UF1-7       |   
[Hyper-parameter Tuning](http://www2.imm.dtu.dk/pubdb/edoc/imm6284.pdf) | 2 |  5 |Objectives include error and prediction time. Hyper-parameters include the number of hidden layers, number of neurons per hidden layer, learning rate, dropout rate and L2 regularization weight penalties

## Results
All the baseline results recorded in our paper are reported in [Google Cloud Drive](https://drive.google.com/drive/folders/1ANE701izoLUNoADnfkngrapyTqHCHyGS).


## Citation
Please cite our paper if you find our work useful for your research:
```
@article{WANG2022,
title = {An Adaptive Batch Bayesian Optimization Approach for Expensive Multi-Objective Optimization},
author = {Hongyan Wang and Hua Xu and Yuan Yuan and Zeqiu Zhang},
journal = {Information Sciences},
year = {2022},
issn = {0020-0255},
}
```
If there's any question, please feel free to contact why17@mails.tsinghua.edu.cn and zzhang77@gwmail.gwu.edu.

## References

<a name="1">
</a>

[1] J. Knowles, [ParEGO: a hybrid algorithm with on-line landscape approximation for expensive multi-objective optimization problems](https://ieeexplore.ieee.org/document/1583627), IEEE Transactions on Evolutionary Computation 10 (1) (2006) 50–66.

<a name="2">
</a>

[2] Q. Zhang, W. Liu, E. Tsang, B. Virginas, [Expensive multi-objective optimization by MOEA/D with Gaussian process model](https://ieeexplore.ieee.org/document/5353656), IEEE Transactions on Evolutionary Computation 14 (3) (2010) 456–474.

<a name="3">
</a>

[3] H. Qian, Y. Yu, [Solving high-dimensional multi-objective optimization problems with low effective di- mensions](https://ojs.aaai.org/index.php/AAAI/article/view/10664), in: Proceedings of the 31st AAAI Conference on Artificial Intelligence, AAAI’17, AAAI Press, 2017, p. 875–881.

<a name="4">
</a>

[4] J. Kirschner, M. Mutny, N. Hiller, R. Ischebeck, A. Krause, [Adaptive and safe bayesian optimization in high dimensions via one-dimensional subspaces](http://proceedings.mlr.press/v97/kirschner19a/kirschner19a.pdf), in: Proceedings of the 36th International Conference on Machine Learning, Vol. 97 of ICML’19, PMLR, Long Beach, California, USA, 2019, pp. 3429–3438.

<a name="5">
</a>

[5] Y. Tian, X. Zhang, C. Wang, Y. Jin, [An Evolutionary Algorithm for Large-Scale Sparse Multiobjective Optimization Problems](https://ieeexplore.ieee.org/document/8720021), IEEE Transactions on Evolutionary Computation 24 (2) (2020) 380–393.


<a name="6">
</a>

[6] Y. Tian, C. Lu, X. Zhang, K. C. Tan, Y. Jin, [Solving large-scale multiobjective optimization problems with sparse optimal solutions via unsupervised neural networks](https://ieeexplore.ieee.org/document/9047876), IEEE Transactions on Cybernetics 51 (6) (2021) 3115–3128.


<a name="7">
</a>

[7] K. Deb, L. Thiele, M. Laumanns, E. Zitzler, [Scalable Test Problems for Evolutionary Multiobjective Optimization](https://link.springer.com/chapter/10.1007/1-84628-137-7_6), Springer London, London, 2005, pp. 105–145.

<a name="8">
</a>

[8] Q. Zhang, A. Zhou, S. Zhao, P. Suganthan, W. Liu, S. Tiwari, [Multiobjective optimization test instances for the cec 2009 special session and competition](https://www.researchgate.net/publication/265432807_Multiobjective_optimization_Test_Instances_for_the_CEC_2009_Special_Session_and_Competition), Mechanical Engineering (2008) 1–30.

<a name="9">
</a>

[9] Huband, S., Barone, L., While, L., Hingston, P. [A Scalable Multi-objective Test Problem Toolkit](https://link.springer.com/chapter/10.1007/978-3-540-31880-4_20). In: Coello Coello, C.A., Hernández Aguirre, A., Zitzler, E. (eds) Evolutionary Multi-Criterion Optimization. EMO 2005.

<a name="10">
</a>

[10] Palm, Rasmus Berg. [Prediction as a candidate for learning deep hierarchical models of data](http://www2.imm.dtu.dk/pubdb/edoc/imm6284.pdf). Technical University of Denmark 5 2012.
