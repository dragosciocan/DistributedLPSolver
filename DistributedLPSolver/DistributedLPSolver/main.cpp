//
//  main.cpp
//  DistributedLPSolver
//
//  Created by Dragos Ciocan on 6/24/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#include <iostream>
#include "instance.h"
#include "time.h"

int main(int argc, const char * argv[])
{
	clock_t t1 = clock();
    using namespace distributed_solver;

    // Add scenarios and path to store instances here...
    int advertiser_dimensions [] = {1000};
    int impression_dimensions [] = {1000};
    long double bid_sparsity_scenario [] = {0.1};

    //int advertiser_dimensions [] = {2};
    //int impression_dimensions [] = {10};
    //long double bid_sparsity_scenario [] = {0.5};

    //int advertiser_dimensions [] = {100000};
    //int impression_dimensions [] = {1000000};
    //long double bid_sparsity_scenario [] = {0.0001};

    std::string file_name_path = "/Users/ciocan/Documents/Google/data/experiment_";
    int num_iterations = 300;
    long double epsilon = 0.01;
    long double numerical_accuracy_tolerance = 0.000000000000000001;

    bool use_binary_search = true;	// Setting this to false runs the sort method
    int num_bin_intervals = 3;	// Sets the number of critical ratios to calculate
    long double cr_transition_scale = 1 - epsilon * 0.001;


    for (int i = 0; i < (sizeof(advertiser_dimensions) / sizeof(int)); ++i) {
        for (int j = 0; j < (sizeof(impression_dimensions) / sizeof(int)); ++j) {
            for (int k = 0; k < (sizeof(bid_sparsity_scenario) / sizeof(long double)); ++k) {
                Instance inst = Instance(advertiser_dimensions[i],
                                         impression_dimensions[j],
                                         1,
                                         bid_sparsity_scenario[k],
                                         epsilon,
                                         0.25,
                                         numerical_accuracy_tolerance);
                inst.GenerateInstance();
                // std::cout << file_name_path + "\n";
                // inst.WriteInstanceToCSV(file_name_path);
                // inst.GenerateAndWriteInstance(file_name_path);

                // std::cout << "created global problem \n";

                if(!use_binary_search){
                	inst.RunMultiplicativeWeights(num_iterations, numerical_accuracy_tolerance, use_binary_search);
                }
                else{
                	inst.RunMultiplicativeWeights(num_iterations, numerical_accuracy_tolerance, use_binary_search,
                		cr_transition_scale, num_bin_intervals);
                }
                std::cout << "finished \n";
            }
        }
    }
    clock_t t2 = clock();

    cout << "Total running time: " << (t2-t1) << endl;

    return 0;
}
