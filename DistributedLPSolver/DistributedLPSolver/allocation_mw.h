//
//  allocation_mw.h
//  InstanceGenerator
//
//  Created by Dragos Ciocan on 6/11/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#ifndef __DistributedLPSolver__allocation_mw__
#define __DistributedLPSolver__allocation_mw__

#include <ext/hash_map>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "instance.h"
#include "global_problem.h"

using namespace std;

namespace distributed_solver {
    class AllocationMW {
        int num_advertisers_;
        int num_impressions_;
        int num_slots_; // Number of ad slots per impression.
        long double bid_sparsity_; // Percentage of impressions advertiser bids for.
        int num_shards_;
        long double width_;
        long double epsilon_;
        long double max_bid_;
        long double numerical_accuracy_tolerance_;

        vector<__gnu_cxx::hash_map<int, long double> >* bids_matrix_;
        vector<__gnu_cxx::hash_map<int, long double> >* transpose_bids_matrix_;
        vector<long double>* budgets_;
        vector<long double> slacks_;
        vector<long double> avg_slacks_;
        vector<__gnu_cxx::hash_map<int, pair<long double, long double> > >* solution_;

        // Multiplicative weights related vars.
        int iteration_count_;
        vector<long double> weights_;
        GlobalProblem global_problem_;

    public:
        AllocationMW(int num_advertisers, int num_impressions, int num_slots, long double bid_sparsity,
                     long double max_bid, long double epsilon, long double numerical_accuracy_tolerance,
                     vector<__gnu_cxx::hash_map<int, long double> >* bids_matrix,
                     vector<__gnu_cxx::hash_map<int, long double> >* transpose_bids_matrix,
                     vector<long double>* budgets,
                     vector<__gnu_cxx::hash_map<int, pair<long double, long double> > >* solution,
                     bool binary);
        AllocationMW(int num_advertisers, int num_impressions, int num_slots, long double bid_sparsity,
                     long double max_bid, long double epsilon, long double numerical_accuracy_tolerance,
                     vector<__gnu_cxx::hash_map<int, long double> >* bids_matrix,
                     vector<__gnu_cxx::hash_map<int, long double> >* transpose_bids_matrix,
                     vector<long double>* budgets,
                     vector<__gnu_cxx::hash_map<int, pair<long double, long double> > >* solution,
                     bool binary, long double scale, int intervals);

        // Generation and output functions.
        void GenerateInstance();
        void WriteInstanceToCSV(std::string file_name_handle);
        void GenerateAndWriteInstance(std::string file_name_handle);
        void SetBudgets(long double factor);

        // Creates current global problem.
        void RunMultiplicativeWeights(int num_iterations, long double numerical_accuracy_tolerance);
        void CreateGlobalProblem();
        void ComputeWeightedBudget();
        void UpdateGlobalProblem();
        void UpdateWeights();
        void CalculateAllocationMWWidth();
        void CalculateSlacks();
        void UpdateAvgSlacks(int t);
        void ReportWorstInfeasibility(int t);
        void ReportWeightStats();
        void ComputeCPLEXRevenue();
        void RunAllocationMW(int num_iterations);

    private:
        void VerifySolution();
        long double CalculateGlobalMWProblemOpt();
    };
}

#endif /* defined(__DistributedLPSolver__allocation_mw__) */
