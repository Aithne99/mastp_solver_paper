#ifdef PREPROCESS_ARRANGEMENT
#include "Arrangement.h"
#else
#include "SolverGurobi.h"
#include "Timer.h"
#endif
#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include "Instance.h"

using namespace std;

//std::random_device rd;
std::mt19937 g(30);


int main(int argc, char *argv[]) {
    Instance instance;

#ifdef PREPROCESS_ARRANGEMENT
    cin >> instance.n_;
    for (int i = 0; i < instance.n_; i++) {
        int x, y;
        cin >> x >> y;
        instance.points_.push_back({ x, y });
    }

    for (int i = 0; i < instance.n_; i++) {
        for (int j = i + 1; j < instance.n_; j++) {
            instance.edges_.push_back({ i, j });
        }
    }

    BuildArrangementData(instance);

    instance.PrintArrangementData(argc, argv);
#else
    instance.load(argc, argv);

    std::vector<double> areas(instance.edges_.size());
    for (int i = 0; i < instance.edges_.size(); ++i)
    {
        for (auto f : instance.hitting_set_[instance.edge_to_circle_[i]])
        {
            areas[i] += instance.face_area_[f];
        }
    }

    std::sort(instance.edge_order.begin(), instance.edge_order.end(), [&areas](int left, int right) { return areas[left] > areas[right]; });

    //std::sort(instance.edge_order.begin(), instance.edge_order.end(), [&instance](int left, int right) { return instance.hitting_set_[instance.edge_to_circle_[left]].size() < instance.hitting_set_[instance.edge_to_circle_[right]].size(); });

    for (auto e : instance.edge_order)
        std::cout << e << " ";

    std::cout << std::endl;

    Timer timer;
    timer.Start();

    SolverGurobi solver(instance);
//  solver.stats_.build_arrangement_time_ = timer.Read() / 1000.;
    solver.stats_.num_faces_ = instance.num_faces_;
    solver.BuildModel();
    solver.Run();
    timer.Pause();

    std::cout << '\a'; // beep
    std::cout << timer.Read();

    _sleep(INT_MAX);
#endif
}
