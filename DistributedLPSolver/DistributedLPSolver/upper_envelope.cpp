//
//  upper_envelope.cpp
//  InstanceGenerator
//
//  Created by Dragos Ciocan on 6/13/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#include "upper_envelope.h"

namespace distributed_solver {
    // 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    // Returns a positive value, if OAB makes a counter-clockwise turn,
    // negative for clockwise turn, and zero if the points are collinear.
    coord2_t cross(const Constraint &O, const Constraint &A, const Constraint &B)
    {
        return (A.coefficient_ - O.coefficient_) * (B.price_ - O.price_) - (A.price_ - O.price_) * (B.coefficient_ - O.coefficient_);
    }
    
    // Returns a list of points on the convex hull in counter-clockwise order.
    vector<Constraint> upper_envelope(vector<Constraint> P, long double numerical_tolerance)
    {
        // Add (0, 0)
        P.push_back(Constraint(0.0, 0.0, -1.0));
        int n = P.size(), k = 0;
        vector<Constraint> H(n);
        // Sort points lexicographically
        sort(P.begin(), P.end(), compare_Constraint_lexicographically());
        // Build lower hull
        for (int i = 0; i < n; i++) {
            while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= numerical_tolerance) k--;
            H[k] = P[i];
            k++;
        }
        H.resize(k);
        //std::sort(H.begin(), H.end(), compare_Constraint_by_weight());
        return H;
    }
}