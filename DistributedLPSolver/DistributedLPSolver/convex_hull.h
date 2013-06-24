//
//  convex_hull.h
//  InstanceGenerator
//
//  Created by Dragos Ciocan on 5/21/13.
//  Copyright (c) 2013 Dragos Ciocan. All rights reserved.
//

#ifndef __DistributedLPSolver__convex_hull__
#define __DistributedLPSolver__convex_hull__

// Implementation of Andrew's monotone chain 2D convex hull algorithm.
// Asymptotic complexity: O(n log n).
// Practical performance: 0.5-1.0 seconds for n=1000000 on a 1GHz machine.
#include <algorithm>
#include <vector>
using namespace std;

typedef long double coord_t;         // coordinate type
typedef long double coord2_t;  // must be big enough to hold 2*max(|coordinate|)^2

struct Point {
    coord_t x, y;
    long double weight;
    int constraint_id;
    
    bool operator <(const Point &p) const {
        return x < p.x || (x == p.x && y < p.y);
    }
};
        
// Returns a list of points on the convex hull in counter-clockwise order.
// Note: the last point in the returned list is the same as the first one.
vector<Point> convex_hull(vector<Point> P);
        
#endif /* defined(__DistributedLPSolver__convex_hull__) */
