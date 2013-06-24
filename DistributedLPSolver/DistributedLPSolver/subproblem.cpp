  //
//  subproblem.cpp
//  InstanceGenerator
//
//  Created by Dragos Ciocan on 2/14/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#include <cmath>
#include <list>
#include <float.h>
#include "convex_hull.h"
#include "upper_envelope.h"
#include "subproblem.h"
#include "time.h"

#include <fstream>

using namespace std;

namespace distributed_solver {
    Constraint::Constraint(long double price, long double coefficient, long double weight) {
        price_ = price;
        coefficient_ = coefficient;
        weight_ = weight;
    }
    
    Constraint::Constraint() {}
    
    void Constraint::set_active(bool value) {
        is_active_ = value;
    }
    
    Subproblem::Subproblem(int num_vars,
                           std::vector<std::pair<long double, long double> >* coefficients,
                           std::vector<int>* advertiser_index) {
        num_vars_ = num_vars;
        for (int i = 0; i < num_vars_; ++i) {
            constraints_.push_back(Constraint((*coefficients)[i].first,
                                              (*coefficients)[i].first * (*coefficients)[i].second,
                                              (*coefficients)[i].second));
            constraints_[i].set_active(true);
        }
        advertiser_index_ = advertiser_index;
    }

    void Subproblem::SolveSubproblem(int iteration, int index) {
        // Eliminate all superfluous constraints.
        std::vector<Constraint> active_constraints;
        for (int i = 0; i < num_vars_; ++i) {
            for (int j = i + 1; j < num_vars_; ++j) {
                if ((constraints_[j].is_active_ == true) &&
                    (constraints_[j].price_ >= constraints_[i].price_) &&
                    (constraints_[j].weight_ <= constraints_[i].weight_)) {
                    constraints_[i].set_active(false);
                }
                if ((constraints_[i].is_active_ == true) &&
                    (constraints_[i].price_ >= constraints_[j].price_) &&
                    (constraints_[i].weight_ <= constraints_[j].weight_)) {
                    constraints_[j].set_active(false);
                }
            }
        }
        for (int i = 0; i < num_vars_; ++i) {
            if (constraints_[i].is_active_ == true) {
                active_constraints.push_back(Constraint(constraints_[i].price_, constraints_[i].coefficient_, constraints_[i].weight_));
            }
        }
                
        if (active_constraints.size() < 1) {
            return;
        }
       
        // Sort list of active constraints. Note that the current implementation creates a clone of the
        // constraints vector. This isn't necessary, should be fixed.
        std::sort(active_constraints.begin(), active_constraints.end(), compare_Constraint_by_weight());
               
        // Go through active constraints, create envelope points.
        envelope_points_.clear();
        // First create intersection with x axis
        envelope_points_.push_back(std::make_pair((active_constraints[0].price_ / active_constraints[0].coefficient_), 0.0));
        for (int k = 0; (k < (active_constraints.size() - 1)); ++k) {
            long double u_value = ((active_constraints[k].price_ - active_constraints[k + 1].price_) /
                              (active_constraints[k].coefficient_ - active_constraints[k + 1].coefficient_));
            envelope_points_.push_back(std::make_pair(u_value,
                                                      (active_constraints[k].price_ -
                                                       active_constraints[k].coefficient_ * u_value)));
        }
        // Last, create intersection with x axis
        envelope_points_.push_back(std::make_pair(0.0, active_constraints[active_constraints.size() - 1].price_));
        // Find the budget ranges where there is a change of basis.
        budget_cutoffs_.clear();
        // Start by addding 0 for convenience.
        budget_cutoffs_.push_back(0.0);
        for (int k = 0; k < envelope_points_.size() - 1; ++k) {
            if ((envelope_points_[k].first - envelope_points_[k + 1].first) > 0.00000000000001) {
                budget_cutoffs_.push_back((envelope_points_[k + 1].second - envelope_points_[k].second) /
                                          (envelope_points_[k].first - envelope_points_[k + 1].first));
            } else {
                budget_cutoffs_.push_back(budget_cutoffs_[k]);
            }
            
            
        }
        // Add +\infty for convenience.
        budget_cutoffs_.push_back(DBL_MAX);
        active_constraints.clear();
    };

    void Subproblem::SolveSubproblemConvexHull(int iteration, int index) {
        if (constraints_.size() < 1) {
            return;
        }
        
        //if (iteration == 6 && index == 715) {
        if (iteration == 2 && index == 715) {
            int y = 0;
            y++;
            cout << y;
        }
        
        // Eliminate all superfluous constraints.
        std::vector<Constraint> active_constraints;
       
        vector<Point> points;
        for (int p = 0;  p < constraints_.size(); ++p) {
            points.push_back(Point());
            points[p].x = (-1.0) * constraints_[p].coefficient_;
            points[p].y = (-1.0) * constraints_[p].price_;
            points[p].weight = constraints_[p].weight_;
            points[p].constraint_id = p;
            constraints_[p].is_active_ = false;
        }
        points.push_back(Point());
        points[constraints_.size()].x = 0.0;
        points[constraints_.size()].y = 0.0;
        points[constraints_.size()].weight = -1;
        points[constraints_.size()].constraint_id = -1;
        
        vector<Point> convex_hull_points = convex_hull(points);
        
        for (int p = 0;  p < convex_hull_points.size(); ++p) {
            if ((convex_hull_points[p].x != 0.0) || (convex_hull_points[p].y != 0.0)) {
                active_constraints.push_back(Constraint(convex_hull_points[p].y * (-1.0),
                                                        convex_hull_points[p].x * (-1.0),
                                                        convex_hull_points[p].x / convex_hull_points[p].y));
                constraints_[convex_hull_points[p].constraint_id].is_active_ = true;
            }
        }
        
        // Sort list of active constraints. Note that the current implementation creates a clone of the
        // constraints vector. This isn't necessary, should be fixed.
        std::sort(active_constraints.begin(), active_constraints.end(), compare_Constraint_by_weight());
        
        // Go through active constraints, create envelope points.
        envelope_points_.clear();
        // First create intersection with x axis
        envelope_points_.push_back(std::make_pair((active_constraints[0].price_ / active_constraints[0].coefficient_),
                                                  0.0));
        for (int k = 0; (k < (active_constraints.size() - 1)); ++k) {
            long double u_value = ((active_constraints[k].price_ - active_constraints[k + 1].price_) /
                                   (active_constraints[k].coefficient_ - active_constraints[k + 1].coefficient_));
            envelope_points_.push_back(std::make_pair(u_value,
                                                      (active_constraints[k].price_ -
                                                       active_constraints[k].coefficient_ * u_value)));
        }
        // Last, create intersection with y axis
        envelope_points_.push_back(std::make_pair(0.0, active_constraints[active_constraints.size() - 1].price_));
        // Find the budget ranges where there is a change of basis.
        budget_cutoffs_.clear();
        // Start by addding 0 for convenience.
        budget_cutoffs_.push_back(0.0);
        for (int k = 0; k < envelope_points_.size() - 1; ++k) {
            if ((envelope_points_[k].first - envelope_points_[k + 1].first) > 0.00000000000001) {
                budget_cutoffs_.push_back((envelope_points_[k + 1].second - envelope_points_[k].second) /
                                          (envelope_points_[k].first - envelope_points_[k + 1].first));
            } else {
                budget_cutoffs_.push_back(budget_cutoffs_[k]);
            }
        }
        // Add +\infty for convenience.
        budget_cutoffs_.push_back(DBL_MAX);
        active_constraints.clear();
    }
    
    void Subproblem::SolveSubproblemConvexHullOptimized(int iteration, int index) {
        if (constraints_.size() < 1) {
            return;
        }
        
        //if (iteration == 6 && index == 715) {
        if (iteration == 2 && index == 715) {
            int y = 0;
            y++;
            cout << y;
        }
        clock_t t1, t2, t3;
        float diff;
        t1 = clock();
        std::vector<Constraint> active_constraints = upper_envelope(constraints_, 0.00000000000001);
        t2 = clock();
        diff = ((float)t2-(float)t1);
        t3 = clock();
        diff = ((float)t3-(float)t2);
        
        // Go through active constraints, create envelope points.
        envelope_points_.clear();
        envelope_points_.push_back(std::make_pair((active_constraints[active_constraints.size() - 2].price_ /
                                                   active_constraints[active_constraints.size() - 2].coefficient_), 0.0));
        for (int k = (active_constraints.size() - 2); k > 0; --k) {
            long double u_value = ((active_constraints[k - 1].price_ - active_constraints[k].price_) /
                                   (active_constraints[k - 1].coefficient_ - active_constraints[k].coefficient_));
            envelope_points_.push_back(std::make_pair(u_value, (active_constraints[k - 1].price_ -
                                                                active_constraints[k - 1].coefficient_ * u_value)));
        }
        envelope_points_.push_back(std::make_pair(0.0, active_constraints[0].price_));
        
        // Find the budget ranges where there is a change of basis.
        budget_cutoffs_.clear();
        // Start by addding 0 for convenience.
        budget_cutoffs_.push_back(0.0);
        for (int k = 0; k < envelope_points_.size() - 1; ++k) {
            if ((envelope_points_[k].first - envelope_points_[k + 1].first) > 0.00000000000001) {
                budget_cutoffs_.push_back((envelope_points_[k + 1].second - envelope_points_[k].second) /
                                          (envelope_points_[k].first - envelope_points_[k + 1].first));
            } else {
                budget_cutoffs_.push_back(budget_cutoffs_[k]);
            }
        }
        // Add +\infty for convenience.
        budget_cutoffs_.push_back(DBL_MAX);
        active_constraints.clear();
    };
}
