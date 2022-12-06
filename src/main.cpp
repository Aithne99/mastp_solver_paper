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
std::mt19937 g(10);


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

    Timer timer;
    timer.Start();

    SolverGurobi solver(instance);
//  solver.stats_.build_arrangement_time_ = timer.Read() / 1000.;
    solver.stats_.num_faces_ = instance.num_faces_;
    solver.BuildModel();
    solver.Run();
    timer.Pause();

    std::cout << '\a'; // beep

    _sleep(INT_MAX);
#endif
}
