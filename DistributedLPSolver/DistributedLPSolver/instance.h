//
//  instance.h
//  DistributedLPSolver
//
//  Created by Dragos Ciocan on 6/24/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#ifndef __DistributedLPSolver__instance__
#define __DistributedLPSolver__instance__
#include <ext/hash_map>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "allocation_mw.h"

using namespace std;

namespace distributed_solver {

    class Instance {
        int num_advertisers_;
        int num_impressions_;
        int num_slots_; // Number of ad slots per impression.
        long double bid_sparsity_; // Percentage of impressions advertiser bids for.
        int num_shards_;
        long double width_;
        long double epsilon_;
        long double scaling_factor_;
        vector<__gnu_cxx::hash_map<int, long double> > bids_matrix_; // Matrix of bids of advertisers for impressions.
        vector<__gnu_cxx::hash_map<int, long double> > transpose_bids_matrix_;
        vector<long double> budgets_;
        vector<__gnu_cxx::hash_map<int, pair<long double, long double> > >* solution_;

        // Multiplicative weights related vars.
        int iteration_count_;
        long double numerical_accuracy_tolerance_;

    public:
        Instance(int num_advertisers, int num_impressions, int num_slots, long double bid_sparsity, long double epsilon,
                 long double scaling_factor, long double numerical_accuracy_tolerance);
        long double max_bid_;

        // Generation and output functions.
        void GenerateInstance();
        void WriteInstanceToCSV(std::string file_name_handle);
        void GenerateAndWriteInstance(std::string file_name_handle);
        void SetBudgets();

        // Creates current global problem.
        void RunMultiplicativeWeights(long double num_iterations, long double numerical_accuracy_tolerance, bool binary);
        void RunMultiplicativeWeights(long double num_iterations, long double numerical_accuracy_tolerance, bool binary, long double scale, int intervals);

        static void UpdateAvgPrimal(int t, vector<__gnu_cxx::hash_map<int, pair<long double, long double> > >* solution);
        void BuildPrimals();
        static void ResetCurrentPrimal(vector<__gnu_cxx::hash_map<int, pair<long double, long double> > >* sol);

    private:
        void ReportGraphTopology();
    };
}

#endif /* defined(__DistributedLPSolver__instance__) */
