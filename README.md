### Improved Approximation Algorithms for k-Submodular Maximization under a Knapsack Constraint
## The source code is organized as follows:
- algs.cpp: contains experimental algorithms, including Alg3, Alg4, SmkDetAcc, SmkRanAcc, SmkStream, Fantom, and SampleGreedy.
- obj.cpp: stores functions for calculating function f.
- binheap.cpp: implements a heap structure.
- mygraph.cpp: includes the data structure.
- gen_er.cpp: used for generating ER datasets.
- Source code of our 2 streaming algorithms, greedy and SGr.
- Facebook dataset (in "data" folder) for testing the algorithms. Due to file size restriction in Github, please find the Sensor dataset in Intel Lab Data (http://db.csail.mit.edu/labdata/labdata.html).

Since estimating F in influence maximization is very time consuming, our code uses OpenMP for parallelization (https://en.wikipedia.org/wiki/OpenMP).

## To build our code, run:

```
	g++ -std=c++11 *.cpp -o ksub -DIL_STD -fopenmp -g
```

## After building, to run our code, run:

```
	./ksub -f <data filename> -c <cost filename>
		-V <size of V>
		-t <type of experiment, 0: influence maximization, 1: sensor placement>
		-k <value of k>
		-B <value of B>
		-b <value of beta>
		-r <value of rho>
		-e <value of epsilon>
		-n <value of eta - denoise step for RStream>
		-g <value of gamma>
		-a <algorithm, 0: Greedy, 1: DStream, 2: RStream, 3: SGr, 4: SampleRstream. Please use SSA source code for testing IM algorithm>
		-p <number of threads (OpenMP) to running algorithms>
```

### Dependencies
- GNU `g++`, `make`
