#pragma once
#include "Framework.h"

class Greedy : public Framework
{
public:
    Greedy(Network *g);
    ~Greedy();
    double get_solution(kseeds &seedsf,bool is_ds,estimate &e_max1,kpoint &fe_max,double &fe);
    double get_solution2(kseeds &seedsf, bool is_ds);
    double get_solution3(kseeds &seedsf, bool is_ds);
    double test_mpi(kseeds &seedsf, bool is_ds);
    double get_solutionT1(kseeds &seedsf);
    double get_solutionT2(kseeds &seedsf);
    int get_no_queries();

private:
    int no_queries;
};
