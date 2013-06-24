
//
//  global_problem.h
//  InstanceGenerator
//
//  Created by Dragos Ciocan on 2/14/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#ifndef DistributedLPSolver_global_problem_h
#define DistributedLPSolver_global_problem_h

#include <cmath>
#include <iostream>
#include <vector>

#include "subproblem.h"

using namespace std;

namespace distributed_solver {
    class GlobalProblem {
    public:
        int num_partitions_;
        long double budget_;
        int num_iterations_;
        bool binary_;
        
        vector<pair<int, long double> > budget_allocation_;
        vector<Subproblem> subproblems_;
        vector<long double> slacks_;
        
        vector<__gnu_cxx::hash_map<int, pair<long double, long double> > >* solution_;
        
        GlobalProblem(int num_partitions, long double max_bid, long double advertiser_indegree, long double numerical_accuracy_tolerance,
                      vector<__gnu_cxx::hash_map<int, pair<long double, long double> > >* solution, bool binary);
        void InitializeInstance();
        void InitializeBudgetAllocation();
        void ConstructPrimal(vector<__gnu_cxx::hash_map<int, long double> >* primal_sol, int iteration);
        
    private:
        void FindOptimalBudgetAllocation();
        void FindOptimalBudgetAllocationBinSearch(long double lower, long double upper,  int num_ratios);
        vector<long double> CalculateAllocationDeltaBin(vector<long double>* critical_ratios, vector<long double>* budget_usage);
        void AllocateCurrentRatio(long double critical_ratio);
        void AllocateCurrentRatios(long double lower_bound, long double upper_bound);
        void ConstructSubproblemPrimal(vector<__gnu_cxx::hash_map<int, long double> >* primal_sol,
                                       int subproblem_index, long double budget_allocation, int opt_region);
        long double numerical_accuracy_tolerance_;
        long double primal_assignment_test_;
    };
    
    class Slope {
    public:
        long double slope_;
        int subproblem_index_;
        int region_index_;
        Slope(long double slope, int subproblem_index, int region_index);
    };
    
    struct compare_Slope
    {
        bool operator() (const Slope & lhs, const Slope & rhs) {
            return lhs.slope_ > rhs.slope_;
        }
    };

}

#endif
