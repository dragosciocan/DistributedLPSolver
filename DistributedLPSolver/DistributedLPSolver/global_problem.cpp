//
//  global_problem.cpp
//  InstanceGenerator
//
//  Created by Dragos Ciocan on 2/14/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#include <algorithm>
#include <cmath>

#include "global_problem.h"
#include "instance.h"
#include "time.h"
#include "float.h"

namespace distributed_solver {
    GlobalProblem::GlobalProblem(int num_partitions, long double max_bid, long double advertiser_indegree,
                                 long double numerical_accuracy_tolerance,
                                 vector<__gnu_cxx::hash_map<int, pair<long double, long double> > >* solution, bool binary) {
        num_partitions_ = num_partitions;
        budget_ = 0;
        solution_ = solution;
        num_iterations_ = 0;
        numerical_accuracy_tolerance_ = numerical_accuracy_tolerance;
        binary_ = binary;

    }

    void GlobalProblem::FindOptimalBudgetAllocationBinSearch(long double lower, long double upper, int num_ratios) {
        vector<long double>* budget_usage = new vector<long double>(num_partitions_);
        long double lower_bound = lower;
        long double upper_bound = upper;
           int number_ratios = num_ratios;
           vector<long double>* critical_ratios = new vector<long double>(number_ratios);
           long double critical_ratio_tolerance = 0.0000000000000001;	//terminating interval

           long double difference = (upper_bound - lower_bound) / (double) number_ratios;
           long double ratio = lower_bound - difference;
           critical_ratios->clear();
           for(int i = 0; i < number_ratios; i++){
        	   ratio += difference;
        	   (*critical_ratios).push_back(ratio);
           }

           //figure out (budget_usage - budget_) for critical ratios
           vector<long double> delta = CalculateAllocationDeltaBin(critical_ratios, budget_usage);


           for (int critical_index = 0; critical_index < delta.size(); critical_index++) {
        	   //critical ratio will start out being too loose (delta > 0), then at index k switches to
        	   //too low, which is when we want to capture the CR bounds.
				if (delta[critical_index] > 0.0) {
					// Allocation > budget_. Need higher CR
					lower_bound = critical_ratios->at(critical_index);
				} else if (delta[critical_index] < 0.0) {
					// Allocation < budget_. Need lower CR
					upper_bound = critical_ratios->at(critical_index);
					break;
				} else {
					//budget perfectly assigned, allocate as is.
					AllocateCurrentRatio(critical_ratios->at(critical_index));
					return;
				}
           }


           if (upper_bound - lower_bound < critical_ratio_tolerance) {
        	   cout << "CR: (" << lower_bound << ", " << upper_bound << ")\n";
        	   //terminate if critical ratio is closely bounded
        	   AllocateCurrentRatios(lower_bound, upper_bound);
        	   return;
           }
           else{
        	   //recurse
        	   FindOptimalBudgetAllocationBinSearch(lower_bound, upper_bound, number_ratios);
           }
       }

       //This method takes in the vector of critical ratios and returns budget allocation difference for all those ratios.
       vector<long double> GlobalProblem::CalculateAllocationDeltaBin(vector<long double>* critical_ratios,
                                                    vector<long double>* budget_usage) {

    	   vector<long double>* total_budget_usage = new vector<long double>(critical_ratios->size());

           for(int critical_ratio_index = 0; critical_ratio_index < critical_ratios->size(); critical_ratio_index++){
        	   (*total_budget_usage)[critical_ratio_index] = 0.0;
        	   for (int i = 0; i < num_partitions_; ++i) {
				   (*budget_usage)[i] = 0;
				   for (int j = 0; j < subproblems_[i].envelope_points_.size() - 1; ++j) {
					   if ((long double)subproblems_[i].envelope_points_[j].first >= (long double)critical_ratios->at(critical_ratio_index)) {
						   (*budget_usage)[i] += subproblems_[i].budget_cutoffs_[j+1] - subproblems_[i].budget_cutoffs_[j];
					   }
				   }
				   (*total_budget_usage)[critical_ratio_index] += (*budget_usage)[i];
			   }
			   //subtract budget_ to get the delta
			   (*total_budget_usage)[critical_ratio_index] -= budget_;
           }

           return (*total_budget_usage);
       }

	void GlobalProblem::AllocateCurrentRatio(long double critical_ratio) {
		// This method updates optimal budget_allocation_ based on critical ratio method.
		for (int i = 0; i < num_partitions_; ++i) {
			budget_allocation_[i].second = 0;
			for (int j = 0; j < subproblems_[i].envelope_points_.size() - 1; ++j) {
				if (subproblems_[i].envelope_points_[j].first >= critical_ratio) {
					budget_allocation_[i] = make_pair(j,
							budget_allocation_[i].second
									+ subproblems_[i].budget_cutoffs_[j + 1]
									- subproblems_[i].budget_cutoffs_[j]);
				}
			}
		}
	}

	void GlobalProblem::AllocateCurrentRatios(long double lower, long double upper){
		//CR is in range (lower, upper)
		long double remaining_budget = budget_;

		//first assign for slope > upper
		for (int i = 0; i < num_partitions_; ++i) {
			budget_allocation_[i].second = 0;
					for (int j = 0; j < subproblems_[i].envelope_points_.size() - 1; ++j) {
						if (subproblems_[i].envelope_points_[j].first > upper) {
							long double allocation_increase =
									subproblems_[i].budget_cutoffs_[j+1]
									- subproblems_[i].budget_cutoffs_[j];
							budget_allocation_[i] = make_pair(j, budget_allocation_[i].second + allocation_increase);
							remaining_budget -= allocation_increase;
						}
					}
		}
		//while remaining_budget positive, assign for slope /in (lower, upper)
		int totalRatiosInRange = 0;
		for (int i = 0; i < num_partitions_; ++i){
			for (int j = 0; j < subproblems_[i].envelope_points_.size() - 1; ++j) {
				if((subproblems_[i].envelope_points_[j].first > lower) && (subproblems_[i].envelope_points_[j].first <= upper)){
					long double allocation_increase = min(remaining_budget,
														subproblems_[i].budget_cutoffs_[j + 1]
														- subproblems_[i].budget_cutoffs_[j]);
					budget_allocation_[i] = make_pair(j, budget_allocation_[i].second + allocation_increase);
					remaining_budget -= allocation_increase;
					totalRatiosInRange++;
					if(remaining_budget == 0.0)
					{
						cout << "total ratios in range " << totalRatiosInRange << endl;
						return;
					}
				}
			}
		}
		cout << "total ratios in range " << totalRatiosInRange << endl;

	}

    void GlobalProblem::FindOptimalBudgetAllocation() {
        long double remaining_budget = budget_;
        
        // Create list of slopes.
        vector<Slope> slopes;
        for (int i = 0; i < num_partitions_; ++i) {
            // Reset budget allocation.
            budget_allocation_[i].second = 0;
            
            if (subproblems_[i].envelope_points_.size() > 1) {
                // -1 because at last envelope return on budget is 0.
                for (int j = 0; j < subproblems_[i].envelope_points_.size() - 1; ++j) {
                    slopes.push_back(Slope(subproblems_[i].envelope_points_[j].first, i, j));
                }
            }
        }
        // Sort list of slopes.
        sort(slopes.begin(), slopes.end(), compare_Slope());
        
        // Allocate budgets in decreasing order of slopes.
        int sp_index;
        for (int k = 0; k < slopes.size(); ++k) {
            sp_index = slopes[k].subproblem_index_;
            long double allocation_increase = min(remaining_budget,
                                                       subproblems_[sp_index].budget_cutoffs_[slopes[k].region_index_ + 1] -
                                                       subproblems_[sp_index].budget_cutoffs_[slopes[k].region_index_]);
            if ((allocation_increase > 0) && (remaining_budget > 0)) {
                budget_allocation_[sp_index] = make_pair(slopes[k].region_index_, budget_allocation_[sp_index].second + allocation_increase);
            }
            remaining_budget = remaining_budget - allocation_increase;
        }
    }
    
    void GlobalProblem::ConstructPrimal(vector<__gnu_cxx::hash_map<int, long double> >* primal_sol, int iteration) {
        num_iterations_++;
        
        // Reset primal solution.
        //Instance::ResetPrimal(primal_sol);
        Instance::ResetCurrentPrimal(solution_);
        
        clock_t t1, t2;
        float diff;
        t1 = clock();
        // Find optimal budget allocations to problems.
        cout << "Solving subproblems \n";
        for (int i = 0; i < num_partitions_; ++i) {
            //subproblems_[i].SolveSubproblem(iteration, i);
            //subproblems_[i].SolveSubproblemConvexHull(iteration, i);
            subproblems_[i].SolveSubproblemConvexHullOptimized(iteration, i);
        }
        t2 = clock();
        diff = ((float)t2-(float)t1);
        cout << "subproblems took  " << diff << "\n";
        
        
        cout << "Find optimal budget allocation \n";
        t1 = clock();

        if(binary_){
        	FindOptimalBudgetAllocation();
        }

        /**
        vector<pair<int, long double> > check_budget = budget_allocation_;
        long double dual_val1 = 0.0;
        for (int i = 0; i < num_partitions_; ++i) {
                    if (budget_allocation_[i].second > 0) {
                        long double u = subproblems_[i].envelope_points_[budget_allocation_[i].first].first;
                        long double v = subproblems_[i].envelope_points_[budget_allocation_[i].first].second;
                        dual_val1 += u * budget_allocation_[i].second + v;
                    }
                }
**/
        else{
        FindOptimalBudgetAllocationBinSearch(0 , 10000, 3);
        }
        /**
        long double dual_val2 = 0.0;
                for (int i = 0; i < num_partitions_; ++i) {
                            if (budget_allocation_[i].second > 0) {
                                long double u = subproblems_[i].envelope_points_[budget_allocation_[i].first].first;
                                long double v = subproblems_[i].envelope_points_[budget_allocation_[i].first].second;
                                dual_val2 += u * budget_allocation_[i].second + v;
                            }
                        }

        double sum_difference = 0.0;

        for(int z = 0; z < check_budget.size(); z++){
        	if(budget_allocation_[z].second - check_budget[z].second!= 0){
        		sum_difference += abs(budget_allocation_[z].second - check_budget[z].second);
        		//cout << "Difference at location " << z << " is " << budget_allocation_[z].second - check_budget[z].second << "\n";
        	}
        }

        cout << "total absolute value of difference in B_i is: " <<sum_difference <<endl;
        cout << "total difference in dual values: " << dual_val1 - dual_val2 << endl;
	**/


        t2 = clock();
        diff = ((float)t2-(float)t1);
        cout << "budget allocation took  " << diff << "\n";
        
        // Calculate primal solution for each subproblem.
        cout << "Constructing primal \n";
        t1 = clock();
        long double dual_val = 0;
        primal_assignment_test_ = 0;
        for (int i = 0; i < num_partitions_; ++i) {
            if (budget_allocation_[i].second > 0) {
                long double u = subproblems_[i].envelope_points_[budget_allocation_[i].first].first;
                long double v = subproblems_[i].envelope_points_[budget_allocation_[i].first].second;
                dual_val += u * budget_allocation_[i].second + v;
                ConstructSubproblemPrimal(primal_sol, i, budget_allocation_[i].second, budget_allocation_[i].first);
            }
        }
        t2 = clock();
        diff = ((float)t2-(float)t1);
        cout << "constructing primal took  " << diff << "\n";
        
        cout << "Dual Value = ";
        cout << dual_val;
        cout << "\n";
    }
    
    void GlobalProblem::ConstructSubproblemPrimal(vector<__gnu_cxx::hash_map<int, long double> >* primal_sol,
                                                  int subproblem_index, long double budget_allocation, int opt_region) {
        // Figure out opt u, v.
        long double u = subproblems_[subproblem_index].envelope_points_[opt_region].first;
        long double v = subproblems_[subproblem_index].envelope_points_[opt_region].second;
        
        long double allocation_value = 0;
         
        // If optimum is u = 0, optimal allocation is greedy wrt price.
        if (u == 0) {
            int max_price_index;
            long double max_price = 0;
            for (int i = 0; i < subproblems_[subproblem_index].num_vars_; ++i) {
                if ((max_price < subproblems_[subproblem_index].constraints_[i].price_) and (subproblems_[subproblem_index].constraints_[i].is_active_)) {
                    max_price_index = i;
                    max_price = subproblems_[subproblem_index].constraints_[i].price_;
                }
            }
            
            (*solution_)[subproblems_[subproblem_index].advertiser_index_->at(max_price_index)][subproblem_index].first = 1;
            allocation_value = (*solution_)[subproblems_[subproblem_index].advertiser_index_->at(max_price_index)][subproblem_index].first * subproblems_[subproblem_index].constraints_[max_price_index].price_;
            
            primal_assignment_test_ += allocation_value;
        }
        
        // If optimum is v = 0, optimal allocation is greedy wrt price/weight ratio.
        if (v == 0) {
            int max_ratio_index;
            long double max_ratio = 0;
            for (int i = 0; i < subproblems_[subproblem_index].num_vars_; ++i) {
                if (max_ratio < (subproblems_[subproblem_index].constraints_[i].price_ /
                                 subproblems_[subproblem_index].constraints_[i].coefficient_) and
                    (subproblems_[subproblem_index].constraints_[i].is_active_)) {
                    max_ratio_index = i;
                    max_ratio = (subproblems_[subproblem_index].constraints_[i].price_ /
                                 subproblems_[subproblem_index].constraints_[i].coefficient_);
                }
            }
             
            (*solution_)[subproblems_[subproblem_index].advertiser_index_->at(max_ratio_index)][subproblem_index].first = budget_allocation / subproblems_[subproblem_index].constraints_[max_ratio_index].coefficient_;
            allocation_value = (*solution_)[subproblems_[subproblem_index].advertiser_index_->at(max_ratio_index)][subproblem_index].first * subproblems_[subproblem_index].constraints_[max_ratio_index].price_;
            primal_assignment_test_ += allocation_value;
        }
        
        // If opt is u, v > 0, the optimum whp only has two positive allocations, which are the
        // solutions to a system of 2 equations.
        if ((u > 0) and (v > 0)) {
            vector<int> tight_constraint_indices;
            for (int i = 0; i < subproblems_[subproblem_index].num_vars_; ++i) {
                if (subproblems_[subproblem_index].constraints_[i].is_active_) {
                    long double slack = subproblems_[subproblem_index].constraints_[i].price_ -
                    (u * subproblems_[subproblem_index].constraints_[i].coefficient_ + v);
                    if (slack < 0) { slack = (-1) * slack;}
                    if (slack < numerical_accuracy_tolerance_) {
                        tight_constraint_indices.push_back(i);
                    }
                }
            }
            if (tight_constraint_indices.size() > 2) {
                cout << "ERROR, PERTURB PRICES \n";
            }
            if (tight_constraint_indices.size() == 1) {
                (*solution_)[subproblems_[subproblem_index].advertiser_index_->
                              at(tight_constraint_indices[0])][subproblem_index].first = fmin(budget_allocation / subproblems_[subproblem_index].constraints_[tight_constraint_indices[0]].coefficient_, 1);
                allocation_value = (*solution_)[subproblems_[subproblem_index].advertiser_index_->
                                                at(tight_constraint_indices[0])][subproblem_index].first *
                                   subproblems_[subproblem_index].constraints_[tight_constraint_indices[0]].coefficient_;
                primal_assignment_test_ += allocation_value;
            }
            if (tight_constraint_indices.size() == 2) {
                int first_index = tight_constraint_indices[0];
                int second_index = tight_constraint_indices[1];
                
                long double x_1 = (budget_allocation -
                              subproblems_[subproblem_index].constraints_[second_index].coefficient_) / (subproblems_[subproblem_index].constraints_[first_index].coefficient_ - subproblems_[subproblem_index].constraints_[second_index].coefficient_);
                (*solution_)[subproblems_[subproblem_index].advertiser_index_->at(first_index)][subproblem_index].first = x_1;
                (*solution_)[subproblems_[subproblem_index].advertiser_index_->at(second_index)][subproblem_index].first = 1 - x_1;
                
                allocation_value = x_1 * subproblems_[subproblem_index].constraints_[first_index].price_ + (1-x_1)*subproblems_[subproblem_index].constraints_[second_index].price_;
                primal_assignment_test_ += allocation_value;
                
            }
        }
        if (((allocation_value - (u * budget_allocation + v)) > numerical_accuracy_tolerance_) ||
            ((allocation_value - (u * budget_allocation + v)) < -numerical_accuracy_tolerance_)) {
            cout << "**************Error at problem************ " << subproblem_index << " equal to " <<
                         allocation_value << " - " << (u * budget_allocation + v) << "\n";
        }
    }
    
    void GlobalProblem::InitializeBudgetAllocation() {
        for (int i = 0; i < num_partitions_; ++i) {
            budget_allocation_.push_back(make_pair(0, 0.0));
        }
    }
    
    Slope::Slope(long double slope, int subproblem_index, int region_index) {
        slope_ = slope;
        subproblem_index_ = subproblem_index;
        region_index_ = region_index;
    }
}
