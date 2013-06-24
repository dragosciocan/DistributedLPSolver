//
//  convex_hull.cpp
//  InstanceGenerator
//
//  Created by Dragos Ciocan on 5/21/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#include "convex_hull.h"

// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
coord2_t cross(const Point &O, const Point &A, const Point &B)
{
    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

// Returns a list of points on the convex hull in counter-clockwise order.
vector<Point> convex_hull(vector<Point> P)
{
    int n = P.size(), k = 0;
    vector<Point> H(2*n);
    
    // Sort points lexicographically
    sort(P.begin(), P.end());
    
    // Build lower hull
    for (int i = 0; i < n; i++) {
        while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0.000000000001) k--;
        H[k] = P[i];
        k++;
    }
    
    H.resize(k);
    return H;
}