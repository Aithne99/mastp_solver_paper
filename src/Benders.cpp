#include "Benders.h"

#include <cassert>
#include <iostream>
#include <set>
#include "Timer.h"

using namespace std;

void DFS(const Instance &instance, const vector<double> &x,
            set<pair<double, int>> &active, vector<bool> &visited, int face,
            vector<double> &w) {
    //visited[face] = true;

    //if (!active.empty()) {
    //    int edge = active.begin()->second;
    //    w[edge] += instance.face_area_[face];
    //} else {
    //    assert(face == 0);
    //}

    //for (const auto &adj : instance.face_graph_[face]) {
    //    if (visited[adj.to_])
    //        continue;

    //    for (int edge : instance.circle_to_edges_[adj.circle_]) {
    //        pair<double, int> item(-x[edge], edge);
    //        if (adj.grows_) {
    //            active.insert(item);
    //        } else {
    //            active.erase(item);
    //        }
    //    }

    //    DFS(instance, x, active, visited, adj.to_, w);

    //    for (int edge : instance.circle_to_edges_[adj.circle_]) {
    //        pair<double, int> item(-x[edge], edge);
    //        if (adj.grows_) {
    //            active.erase(item);
    //        } else {
    //            active.insert(item);
    //        }
    //    }
    //}
    // and this is how you do it without needing 8000000 stack frames
    // edges which are present in the current solution are always weighed by all the faces which are part of their hitting set AND haven't been counted yet
    for (int e = 0; e < instance.edges_.size(); ++e)
    {
        if (x[e] == 1)
        {
            for (auto f : instance.hitting_set_[instance.edge_to_circle_[e]])
            {
                if (visited[f])
                    continue;
                w[e] += instance.face_area_[f];
                visited[f] = true;
            }
        }
    }
    // edges which are not present in the current solution are weighed by the new faces they would be bringing in
    for (int e = 0; e < instance.edges_.size(); ++e)
    {
        if (x[e] != 1)
        {
            for (auto f : instance.hitting_set_[instance.edge_to_circle_[e]])
            {
                if (visited[f])
                    continue;
                w[e] += instance.face_area_[f];
                //visited[f] = true; // faces which would be brought into the solution by multiple currently inactive edges should be counted towards the weight of all edges they belong to
            }
        }
    }
}

vector<double> SeparateBendersCut(const Instance &instance,
                                    const vector<double> &x, const double ub) {

    set<pair<double, int>> active;
    vector<double> w(instance.edges_.size());
    vector<bool> visited(instance.num_faces_);

    //Timer timer;
    //timer.Start();

    DFS(instance, x, active, visited, 0, w);
    //timer.Pause();
    //std::cout << "\n Time spent in Benders DFS: " << timer.Read();

    double sum = 0.0;

    double ub_of_lb = 0.;
    for (int e = 0; e < instance.edges_.size(); e++)
    {
        ub_of_lb += w[e] * x[e];
    }
    if (ub_of_lb < ub)
    {
        return vector<double>();
    }
    double violation = -x[instance.edges_.size()];
    double obj = -violation;
    for (int i = 0; i < instance.edges_.size(); i++)
    {
        violation += x[i] * w[i];
        sum += w[i];
    }

    if (obj > 1.e-2)
        violation /= obj;

    if (violation < 1.e-4)
        return vector<double>();

    return w;
}