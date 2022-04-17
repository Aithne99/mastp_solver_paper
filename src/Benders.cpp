#include "Benders.h"

#include <cassert>
#include <iostream>
#include <set>

using namespace std;

void DFS(const Instance &instance, const vector<double> &x,
            set<pair<double, int>> &active, vector<bool> &visited, int face,
            vector<double> &w) {

    // this is how you do it with recursion, readably
    //visited[face] = true;

    //if (!active.empty())
    //{
    //    int edge = active.begin()->second;
    //    w[edge] += instance.face_area_[face];
    //} else
    //{
    //    assert(face == 0);
    //}

    //for (const auto &adj : instance.face_graph_[face]) {
    //    if (visited[adj.to_])
    //        continue;

    //    // this is always a 1-1 correspondence, I have no idea why this was treated as a for loop
    //    int edge = adj.circle_;
    //    pair<double, int> item(-x[edge], edge);
    //    if (adj.grows_)
    //    {
    //        active.insert(item);
    //    }
    //    else
    //    {
    //        active.erase(item);
    //    }    

    //    DFS(instance, x, active, visited, adj.to_, w);

    //    if (adj.grows_)
    //    {
    //        active.erase(item);
    //    }
    //    else
    //    {
    //        active.insert(item);
    //    }
    //}

    // and this is how you do it without needing 800000 stack frames
    for (auto e : instance.edge_to_circle_)
    {
        if (x[e] == 1)
        {
            for (auto f : instance.hitting_set_[e])
            {
                if (visited[f])
                    continue;
                w[e] += instance.face_area_[f];
                visited[f] = true;
            }
        }
    }
    for (auto e : instance.edge_to_circle_)
    {
        if (x[e] == 0)
        {
            for (auto f : instance.hitting_set_[e])
            {
                if (visited[f])
                    continue;
                w[e] += instance.face_area_[f];
                visited[f] = true;
            }
        }
    }
}

vector<double> SeparateBendersCut(const Instance &instance,
                                    const vector<double> &x, const double ub) {

    set<pair<double, int>> active;
    vector<double> w(instance.edges_.size());
    vector<bool> visited(instance.num_faces_);

    DFS(instance, x, active, visited, 0, w);

    double sum = 0.0;

    double ub_of_lb = 0.;
    for (int e = 0; e < instance.edges_.size(); e++) {
    ub_of_lb += w[e] * x[e];
    }
    if (ub_of_lb < ub) {
    return vector<double>();
    }
    double violation = -x[instance.edges_.size()];
    double obj = -violation;
    for (int i = 0; i < instance.edges_.size(); i++) {
    violation += x[i] * w[i];
    sum += w[i];
    ;
    }

    if (obj > 1.e-2)
    violation /= obj;

    if (violation < 1.e-4)
    return vector<double>();

    return w;
}