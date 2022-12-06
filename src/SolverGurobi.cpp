#include "SolverGurobi.h"
#include "DisjointSetForest.h"

#include <iomanip>
#include <iostream>
#include <cmath>

using namespace std;

void SolverGurobi::BuildModel() {
    env_ = new GRBEnv();

    env_->set(GRB_StringParam_LogFile, "logfile.txt");

    model_ = new GRBModel(*env_);

    model_->getEnv().set(GRB_IntParam_DualReductions, 0);
    model_->getEnv().set(GRB_IntParam_PreCrush, 1);
  // Variables: [edges, faces]

#ifdef BENDERS
    const int num_vars = instance_.edges_.size() + 1;

    auto vars = model_->addVars(num_vars - 1, GRB_BINARY);
    model_->addVar(0, GRB_INFINITY, 1.f, GRB_CONTINUOUS);

    model_->update();
    // Warm start benders
    vector<double> x(num_vars);
    for (int k = 0; k < num_vars - 1 && false; k++) {
        x[k] = 1.;

        auto cut = SeparateBendersCut(instance_, x);

        vector<int> bc_ind(num_vars);
        vector<double> bc_val(num_vars);
        for (int i = 0; i < instance_.edges_.size(); i++) {
            bc_ind[i] = i;
            bc_val[i] = cut[i];
        }
        bc_ind.back() = instance_.edges_.size();
        bc_val.back() = -1.;

        int beg[1] = { 0 };
        char sense[1] = { 'L' };
        double rhs[1] = { 0. };
        
        GRBLinExpr* expr = new GRBLinExpr();
        expr->addTerms(bc_val.data(), vars, bc_val.size());

        model_->addConstrs(expr, sense, rhs, nullptr, num_vars);
        x[k] = 0.;
    }
#else
    const int num_vars = instance_.edges_.size() + instance_.num_faces_;
    vector<double> var_obj(num_vars);
    vector<char> var_type(num_vars);
    for (int i = 0; i < instance_.edges_.size(); i++) {
        var_type[i] = GRB_CONTINUOUS;
    }
    for (int i = 0; i < instance_.num_faces_; i++) {
        const int var_index = instance_.edges_.size() + i;
        var_type[var_index] = GRB_BINARY;
        var_obj[var_index] = instance_.face_area_[i];
    }

    model_->addVars(nullptr, nullptr, var_obj.data(), var_type.data(), nullptr, num_vars);
    model_->update();
    // Covering constraints
    vector<bool> visited(instance_.num_faces_);
    set<int> inclusion_set;

#ifdef PREPROCESS
    AddReducedCoveringConstraintsDFS(visited, 0, inclusion_set);
#else
    AddCoveringConstraintsDFS(visited, 0, inclusion_set);
#endif
#endif
    // Spanning tree constraints

    AddCardinalityConstraint();
}

void CutCallback::callback()
{
    if (/*where == GRB_CB_MIPNODE || */ where == GRB_CB_MIPSOL)
    {
        int node_depth;
        double nd = getDoubleInfo(GRB_CB_MIPSOL_NODCNT);
        node_depth = (int)nd;
        const Instance& instance = solver.instance_;
        auto& stats = solver.stats_;
        auto vars = solver.model_->getVars();
#ifdef BENDERS
        const int num_vars = instance.edges_.size() + 1;
#else
        const int num_vars = instance.edges_.size() + instance.num_faces_;
#endif
        double* xptr = getSolution(vars, num_vars);
        vector<double> x(xptr, xptr + num_vars);
        //if (wherefrom == CPX_CALLBACK_MIP_CUT_FEAS) {
        //    int sol_src = -1;
        //    if (CPXgetcallbacknodeinfo(env, wherefrom, 0,
        //        CPX_CALLBACK_INFO_LAZY_SOURCE, &sol_src)) {

        //        cerr << "GUROBI - CutCallback - Failed to get candidate solution source."
        //            << endl;
        //    }
        //}

        int num_cuts_added = 0;

        auto secs = SeparateSEC(instance, x);
        stats.num_sec_separations_++;
        if (node_depth == 0) {
            stats.num_sec_separations_root_++;
        }

        sort(secs.begin(), secs.end());
        reverse(secs.begin(), secs.end());
#ifndef MAT
        if (!secs.empty())
            secs.resize(1);
#endif
        for (const auto& sec : secs) {
            AddSEC(sec.second);
            num_cuts_added++;

            stats.num_sec_++;
            if (node_depth == 0) {
                stats.num_sec_root_++;
            }
        }

#ifdef SEPARATE_TRIANGLE_CUTS
        if (num_cuts_added == 0) {
            int num_triineq =
                SeparateTriangleInequalities(x);
            num_cuts_added++;

            stats.num_triineq_separations_++;
            stats.num_triineq_ += num_triineq;
            if (node_depth == 0) {
                stats.num_triineq_separations_root_++;
                stats.num_triineq_root_ += num_triineq;
            }
        }
#endif

#ifdef BENDERS

        bool integral = true;
        for (int e = 0; e < instance.edges_.size(); e++) {
            if (min(fabs(x[e] - 0), fabs(1 - x[e])) > 1.e-4) {
                integral = false;
                break;
            }
        }
        if (num_cuts_added == 0 || integral) {
            double early_stop_th = integral || solver.benders_early_stop_ub ==
                numeric_limits<double>::infinity()
                ? 0.
                : solver.benders_early_stop_ub;

            auto benders_cut = SeparateBendersCut(instance, x, early_stop_th);
            solver.stats_.num_benders_separations_++;
            if (node_depth == 0) {
                solver.stats_.num_benders_separations_root_++;
            }

            if (!benders_cut.empty()) {
                AddBendersCut(benders_cut);
                num_cuts_added++;
                //std::cout << "\r\n Benders cut: ";
                //for (auto coef : benders_cut)
                //    std::cout << coef << " ";
                //std::cout << "\r\n";
                double bvio = -x.back();
                for (int i = 0; i < instance.edges_.size(); i++) {
                    bvio += x[i] * benders_cut[i];
                }

                solver.stats_.num_benders_cuts_++;
                if (node_depth == 0) {
                    solver.stats_.num_benders_cuts_root_++;
                }
            }
        }
#endif

        if (num_cuts_added) {
            return;
        }

        if (node_depth == 0 && false) {
            //double dual_bound = 0.;
            //if (CPXgetcallbacknodeinfo(env, wherefrom, 0,
            //    CPX_CALLBACK_INFO_NODE_OBJVAL, &dual_bound)) {

            //    cerr << "GUROBI - CutCallback - Failed to get dual bound" << endl;
            //}

            //solver.stats_.objective_lower_bound_root_ = dual_bound;
        }
    }
#ifndef  MAT

    if (where == GRB_CB_MIPNODE)
    {
        if (getIntInfo(GRB_CB_MIPNODE_STATUS) != GRB_OPTIMAL)
        {
            return;
        }
        const Instance& instance = solver.instance_;

        auto vars = solver.model_->getVars();
        int numVars = solver.model_->get(GRB_IntAttr_NumVars);

        double* x = getNodeRel(vars, numVars);

        vector<int> tree = solver.ClosestTree(x);
        solver.ImproveTree(tree);

        if (tree.size() != instance.n_ - 1) {
            std::cout << "what" << '\n';
            return;
        }

        solver.stats_.num_runs_heuristic_++;

#ifdef BENDERS
            const int num_vars = instance.edges_.size() + 1;
#else
            const int num_vars = instance.edges_.size() + instance.num_faces_;
#endif
            fill(x, x + num_vars, 0);
            for (int i : tree) {
                x[i] = 1.0;
            }
            vector<bool> visited(instance.num_faces_);
#ifdef BENDERS
            double objval_p = 0;
            objval_p = solver.DfsComputeObj(visited, 0, 0, x);
            x[instance.edges_.size()] = objval_p;
#else
            solver.DfsComputeY(visited, 0, 0, x);

            double objval_p = 0;
            for (int i = 0; i < instance.num_faces_; i++) {
                if (x[i + instance.edges_.size()] == 1.) {
                    objval_p += instance.face_area_[i];
                }
            }
#endif

            //*useraction_p = CPX_CALLBACK_SET;

            setSolution(vars, x, num_vars);
            useSolution();
            double node_depth = getDoubleInfo(GRB_CB_MIPNODE_NODCNT);

            if (node_depth == 0) {
                solver.stats_.objective_upper_bound_root_ =
                    max(solver.stats_.objective_upper_bound_root_, objval_p);
            }

            solver.benders_early_stop_ub = min(solver.benders_early_stop_ub, objval_p);
        //}
    }
#endif
}

int SolverGurobi::ArcIndex(const int root, const int edge,
    const bool reversed) const {
    return instance_.edges_.size() * (2 * root + 1) + 2 * edge + reversed;
}

void SolverGurobi::AddImplicationConstraint(const int a, const int b) {
    const int a_var = instance_.edges_.size() + a;
    const int b_var = instance_.edges_.size() + b;
    GRBVar v[2];
    v[0] = model_->getVar(b_var);
    v[1] = model_->getVar(a_var);
    double val[2] = { 1.0, -1.0 };
    int beg[1] = { 0 };
    char sense[1] = { 'G' };

    GRBLinExpr* expr = new GRBLinExpr();
    expr->addTerms(val, v, 2);
    model_->addConstrs(expr, sense, nullptr, nullptr, 1);
}

void SolverGurobi::AddCoveringConstraint(const int face, const int edge) {
    const int face_var = instance_.edges_.size() + face;
    int ind[2] = { face_var, edge };
    double val[2] = { 1.0, -1.0 };
    char sense[1] = { 'G'};

    GRBVar v[2];
    v[0] = model_->getVar((int)face_var);
    v[1] = model_->getVar((int)edge);
    GRBLinExpr* expr = new GRBLinExpr();
    expr->addTerms(val, v, 2);
    model_->addConstrs(expr, sense, nullptr, nullptr, 1);
}

void SolverGurobi::AddArcEdgeCouplingConstraints() {
    //const int num_edges = instance_.edges_.size();
    //for (int r = 0; r < instance_.n_; r++) {
    //    for (int e = 0; e < num_edges; e++) {
    //        int ind[3] = { e, ArcIndex(r, e, false), ArcIndex(r, e, true) };
    //        double val[3] = { -1., 1., 1. };
    //        int beg[1] = { 0 };
    //        char sense[1] = { 'E' };
    //        double rhs[1] = { 0. };
    //        GRBVar v[3];
    //        v[0] = model_->getVar(ind[0]);
    //        v[1] = model_->getVar(ind[1]);
    //        v[2] = model_->getVar(ind[2]);
    //        GRBLinExpr* expr = new GRBLinExpr();
    //        expr->addTerms(val, v, 3);
    //        model_->addConstrs(expr, sense, rhs, nullptr, 1);
    //    }
    //}
}

void SolverGurobi::AddCardinalityConstraint() {
    const int num_edges = instance_.edges_.size();
    vector<int> card_ind(num_edges);
    vector<double> card_val(num_edges, 1.0);
    for (int i = 0; i < num_edges; i++) {
        card_ind[i] = i;
    }

    int beg[1] = { 0 };
    char sense[1] = { 'E' };
#ifdef MAT
    double rhs[1] = { double(instance_.n_) };
#else
    double rhs[1] = { double(instance_.n_ - 1) };
#endif
    auto vars = model_->getVars();
    GRBLinExpr *expr = new GRBLinExpr();
    expr->addTerms(card_val.data(), vars, card_val.size());
    model_->addConstrs(expr, sense, rhs, nullptr, 1);
#ifdef MAT  
    double one[1] = { 1 };
    double two[1] = { 2 };

    for (int i = 0; i < instance_.n_; ++i)
    {
        GRBLinExpr* expr2 = new GRBLinExpr();
        for (int j = 0; j < instance_.edges_.size(); ++j)
        {
            if (instance_.edges_[j].first == i || instance_.edges_[j].second == i)
            {
                expr2->addTerms(one, &vars[j], 1);
            }
        }
        model_->addConstrs(expr2, sense, two, nullptr, 1);
    }
#endif
}

void SolverGurobi::AddCoveringConstraintsDFS(vector<bool>& visited, const int face,
    set<int>& circles) {

    visited[face] = true;

    for (const int circle : circles) {
        for (const int edge : instance_.circle_to_edges_[circle]) {
            AddCoveringConstraint(face, edge);
        }
    }

    for (const auto& e : instance_.face_graph_[face]) {
        if (visited[e.to_])
            continue;

        if (e.grows_) {
            circles.insert(e.circle_);
        }
        else {
            circles.erase(e.circle_);
        }

        AddCoveringConstraintsDFS(visited, e.to_, circles);

        if (e.grows_) {
            circles.erase(e.circle_);
        }
        else {
            circles.insert(e.circle_);
        }
    }


    //for (int i = 0; i < instance_.edges_.size(); ++i)
    //{
    //    int circle = instance_.edge_to_circle_[i];
    //    for (auto f : instance_.hitting_set_[circle])
    //        AddCoveringConstraint(f, i);
    //}
}

void SolverGurobi::AddReducedCoveringConstraintsDFS(vector<bool>& visited,
    const int face,
    set<int>& circles) {
    visited[face] = true;

    int parent = -1;
    int extra_circle = 0;

    for (const auto& e : instance_.face_graph_[face]) {
        if (visited[e.to_])
            continue;

        if (e.grows_) {
            circles.insert(e.circle_);
        }
        else {
            circles.erase(e.circle_);

            const int new_w = instance_.circle_to_edges_[e.circle_].size();
            const int old_w = instance_.circle_to_edges_[extra_circle].size();

            if (parent == -1 || new_w < old_w) {
                parent = e.to_;
                extra_circle = e.circle_;
            }
        }

        AddReducedCoveringConstraintsDFS(visited, e.to_, circles);

        if (e.grows_) {
            circles.erase(e.circle_);
        }
        else {
            circles.insert(e.circle_);
        }
    }

    if (parent == -1) {
        for (const int circle : circles) {
            for (const int edge : instance_.circle_to_edges_[circle]) {
                AddCoveringConstraint(face, edge);
            }
        }
    }
    else {
        AddImplicationConstraint(parent, face);
        for (const int edge : instance_.circle_to_edges_[extra_circle]) {
            AddCoveringConstraint(face, edge);
        }
    }
}

void CutCallback::AddBendersCut(const vector<double>& opt_cut) {
    const auto& instance = solver.instance_;

    vector<int> xs = { 3, 5, 8, 18, 23, 28, 31, 47, 50, 55, 56 };
    double norm = 0.;
    for (const double& e : opt_cut)
        norm += e * e;
    norm += 1;
    norm = sqrt(norm);

    vector<int> ind;
    vector<double> val;
    double xseval = 0.;
    for (int i = 0; i < instance.edges_.size(); i++) {
        if (binary_search(xs.begin(), xs.end(), i))
            xseval += opt_cut[i];
        if (fabs(opt_cut[i]) < 1.e-9)
            continue;

        ind.push_back(i);
        val.push_back(opt_cut[i]);
    }
    ind.push_back(instance.edges_.size());
    val.push_back(-1.);

    GRBVar* v = solver.model_->getVars();
    vector<GRBVar> _v;
    for (auto i : ind)
    {
        _v.push_back(solver.model_->getVar(i));
    }

    GRBLinExpr lhs;
    lhs.addTerms(val.data(), _v.data(), val.size());
    char sense = 'L';
    double rhs = 0.;
    addLazy(lhs, sense, rhs);
}

void CutCallback::AddCutset(const int root, const vector<int>& arcs) {

    const Instance& instance = solver.instance_;

    vector<int> ind = arcs;
    for (auto& i : ind)
        i += (root * 2 + 1) * instance.edges_.size();
    vector<double> val(arcs.size(), 1.);

    GRBVar* v = solver.model_->getVars();
    vector<GRBVar> _v;
    int i = 0;
    for (int j = 0; j < solver.model_->get(GRB_IntAttr_NumVars); ++j)
    {
        if (ind[i] == j)
        {
            _v.push_back(v[j]);
            ++i;
        }
    }

    GRBLinExpr lhs;
    lhs.addTerms(val.data(), _v.data(), val.size());

    char sense = 'G';
    double rhs = 1.;

    addLazy(lhs, sense, rhs);
}

void CutCallback::AddSEC(const vector<int>& S) {
    if (S.empty())
        return;

    const Instance& instance = solver.instance_;

    vector<int> ind;
    vector<double> val;
    vector<int> xs = { 3, 5, 8, 18, 23, 28, 31, 47, 50, 55, 56 };
    double xs_sec_vio = 0.;
    for (int i = 0; i < instance.edges_.size(); i++) {
        const int a = instance.edges_[i].first;
        const int b = instance.edges_[i].second;
        if (!binary_search(S.begin(), S.end(), a))
            continue;
        if (!binary_search(S.begin(), S.end(), b))
            continue;

        ind.push_back(i);
        val.push_back(1.0);
    }
    const double rhs = S.size() - 1;
    xs_sec_vio -= rhs;
    for (int i = 0; i < ind.size(); i++) {
        if (binary_search(xs.begin(), xs.end(), ind[i])) {
            xs_sec_vio += val[i];
        }
    }

    GRBVar* v = solver.model_->getVars();
    vector<GRBVar> _v;
    int i = 0;
    for (auto i : ind)
    {
        _v.push_back(solver.model_->getVar(i));
    }

    GRBLinExpr lhs;
    lhs.addTerms(val.data(), _v.data(), val.size());

    char sense = 'L';

    addLazy(lhs, sense, rhs);
}

double SolverGurobi::DfsComputeObj(vector<bool>& visited, int face, int count,
    double* x) {
    double obj = 0.;

    //visited[face] = true;

    //if (count) {
    //    obj += instance_.face_area_[face];
    //}

    //for (const auto &e : instance_.face_graph_[face]) {
    //    if (visited[e.to_])
    //        continue;

    //    int delta = 0;
    //    for (int i : instance_.circle_to_edges_[e.circle_]) {
    //        if (x[i] > 0.) {
    //            delta += e.grows_ ? 1 : -1;
    //            break;
    //        }
    //    }
    //    obj += DfsComputeObj(visited, e.to_, count + delta, x);
    //}

    for (int e = 0; e < instance_.edges_.size(); ++e)
    {
        if (x[e])
        {
            for (int f : instance_.hitting_set_[instance_.edge_to_circle_[e]])
            {
                if (visited[f])
                    continue;
                visited[f] = true;
                obj += instance_.face_area_[f];
            }
        }
    }

    return obj;
}

void SolverGurobi::DfsComputeY(vector<bool>& visited, int face, int count,
    double* x) {

    for (auto e : instance_.edge_to_circle_)
    {
        if (x[e])
        {
            for (auto f : instance_.hitting_set_[e])
            {
                //if (visited[f])
                //    continue;
                //visited[f] = true;
                x[instance_.edges_.size() + f] = 1.0;
            }
        }
    }
}

vector<int> SolverGurobi::ClosestTree(double* x) {
    vector<int> tree;
    vector<int> ord(instance_.edges_.size());
    for (int i = 0; i < instance_.edges_.size(); i++) {
        ord[i] = i;
    }
    sort(ord.begin(), ord.end(), [x](int i, int j) { return x[i] > x[j]; });

    DisjointSetForest components;
    components.Reset(instance_.n_);
    for (int i : ord) {
        int a = instance_.edges_[i].first;
        int b = instance_.edges_[i].second;

        if (components.Find(a) == components.Find(b))
            continue;

        components.Join(a, b);

        tree.push_back(i);
    }

    return tree;
}

int FindBestFaceInTwoCircles(const Instance& instance, const vector<double>& x,
    const int c1, const int c2) {
    if (instance.hitting_set_[c1].size() > instance.hitting_set_[c2].size())
        return FindBestFaceInTwoCircles(instance, x, c2, c1);

    double best_value = numeric_limits<double>::infinity();
    int best_face = 0;

    for (const int f : instance.hitting_set_[c1]) {
        if (binary_search(instance.hitting_set_[c2].begin(),
            instance.hitting_set_[c2].end(), f)) {

            const double value = x[instance.edges_.size() + f];
            if (value < best_value) {
                best_value = value;
                best_face = f;
            }
        }
    }

    return best_face;
}

int CutCallback::SeparateTriangleInequalities(const vector<double>& x) {
    // The idea here is to fix one side of the triangle and compute
    // the best two faces for each vertex completing the triangle.

    const Instance& instance = solver.instance_;

    vector<pair<int, int>> best_face(instance.n_);
    vector<pair<int, int>> tip_edge(instance.n_);

    int num_vio_triang = 0;
    for (int k = 0; k < instance.edges_.size(); k++) {
        const int c_k = instance.edge_to_circle_[k];

        fill(best_face.begin(), best_face.end(), pair<int, int>(0, 0));

        const int r = instance.edges_[k].first;
        const int s = instance.edges_[k].second;

        for (int i = 0; i < instance.edges_.size(); i++) {
            if (i == k)
                continue;

            const int c_i = instance.edge_to_circle_[i];

            if (instance.edges_[i].first == r || instance.edges_[i].second == r) {
                const int t = instance.edges_[i].first == r ? instance.edges_[i].second
                    : instance.edges_[i].first;

                best_face[t].first = FindBestFaceInTwoCircles(instance, x, c_k, c_i);
                tip_edge[t].first = i;
            }

            if (instance.edges_[i].first == s || instance.edges_[i].second == s) {
                const int t = instance.edges_[i].first == s ? instance.edges_[i].second
                    : instance.edges_[i].first;

                best_face[t].second = FindBestFaceInTwoCircles(instance, x, c_k, c_i);
                tip_edge[t].second = i;
            }
        }

        for (int t = 0; t < instance.n_; t++) {
            if (best_face[t].first == 0 || best_face[t].second == 0)
                continue;

            const int i = tip_edge[t].first;
            const int j = tip_edge[t].second;
            const int f_a = best_face[t].first;
            const int f_b = best_face[t].second;
            const int var_f_a = instance.edges_.size() + f_a;
            const int var_f_b = instance.edges_.size() + f_b;

            const double violation = x[i] + x[j] + x[k] - x[var_f_a] - x[var_f_b];
            if (violation > 1.e-5) {
                num_vio_triang++;

                vector<int> ind = { i, j, k, var_f_a, var_f_b };
                vector<double> coef = { -1., -1., -1., 1., 1. };

                if (var_f_a == var_f_b) {
                    ind.pop_back();
                    coef.pop_back();
                    coef.back() += 1.0;
                    continue;
                }

                GRBVar* v = solver.model_->getVars();
                vector<GRBVar> _v;
                for (auto i : ind)
                {
                    _v.push_back(solver.model_->getVar(i));
                }
                GRBLinExpr lhs;
                lhs.addTerms(coef.data(), _v.data(), coef.size());

                char sense = 'G';
                double rhs = 0.;

                addLazy(lhs, sense, rhs);
            }
        }
    }

    return num_vio_triang;
}

void SolverGurobi::ImproveTree(vector<int>& tree) {
    vector<int> count(instance_.num_faces_);
    for (int edge : tree) {
        for (int face : instance_.hitting_set_[edge]) {
            count[face]++;
        }
    }

    DisjointSetForest dsf;

    while (true) {
        int best_remove = -1;
        int best_add = -1;
        double best_delta = 0.;
        for (int remove : tree) {
            double area_decrease = 0.;
            for (int face : instance_.hitting_set_[remove]) {
                count[face]--;
                if (count[face] == 0) {
                    area_decrease += instance_.face_area_[face];
                }
            }

            dsf.Reset(instance_.n_);
            for (int edge : tree) {
                if (edge != remove) {
                    dsf.Join(instance_.edges_[edge].first, instance_.edges_[edge].second);
                }
            }

            for (int add = 0; add < instance_.edges_.size(); add++) {
                if (add == remove)
                    continue;
                if (dsf.Find(instance_.edges_[add].first) ==
                    dsf.Find(instance_.edges_[add].second))
                    continue;

                double area_increase = 0;
                for (int face : instance_.hitting_set_[add]) {
                    if (count[face] == 0) {
                        area_increase += instance_.face_area_[face];
                    }
                }

                double delta = area_increase - area_decrease;
                if (delta < best_delta) {
                    best_delta = delta;
                    best_remove = remove;
                    best_add = add;
                }
            }

            for (int face : instance_.hitting_set_[remove]) {
                count[face]++;
            }
        }

        if (best_delta > -1.e-2)
            break;

        for (int& edge : tree) {
            if (edge == best_remove) {
                edge = best_add;
                break;
            }
        }
        for (int face : instance_.hitting_set_[best_remove]) {
            count[face]--;
        }
        for (int face : instance_.hitting_set_[best_add]) {
            count[face]++;
        }
    }
}

void SolverGurobi::SetSolverParameters() {
    model_->set(GRB_DoubleParam_TimeLimit, 3 * 60 * 60);

    model_->set(GRB_IntParam_Threads, 1);
    model_->set(GRB_IntParam_LazyConstraints, 1);
}

void SolverGurobi::Run() {
    SetSolverParameters();

    stats_.num_sec_ = 0;
    stats_.num_sec_separations_ = 0;
    stats_.num_sec_root_ = 0;
    stats_.num_sec_separations_root_ = 0;
    stats_.num_triineq_ = 0;
    stats_.num_triineq_separations_ = 0;
    stats_.num_triineq_root_ = 0;
    stats_.num_triineq_separations_root_ = 0;
    stats_.num_runs_heuristic_ = 0;
    stats_.num_vars_ = model_->get(GRB_IntAttr_NumVars);
    stats_.num_constraints_ = model_->get(GRB_IntAttr_NumConstrs);
    stats_.num_benders_cuts_ = 0;
    stats_.num_benders_cuts_root_ = 0;
    stats_.num_benders_separations_ = 0;
    stats_.num_benders_separations_root_ = 0;

    benders_early_stop_ub = numeric_limits<double>::infinity();

    //model_->write("./tmp/cpx.lp");

    //GRBsetintparam(env_, CPXPARAM_ScreenOutput, CPX_ON);

    CutCallback* cb = new CutCallback(*this);

    model_->setCallback(cb);

    // Solve root
    // Gurobi already provides root relaxation statistics
    //model_->set(GRB_DoubleParam_NodeLimit, 0);
    //    
    //model_->optimize();

    //stats_.runtime_root_ = model_->get(GRB_DoubleAttr_Runtime);

    //stats_.objective_upper_bound_root_ = model_->get(GRB_DoubleAttr_ObjVal);
    //stats_.objective_lower_bound_root_ = model_->get(GRB_DoubleAttr_ObjBound);

    //double* root_x =model_->get(GRB_DoubleAttr_X, model_->getVars(), model_->get(GRB_IntAttr_NumVars));
    //if (CPXgetx(env_, model_, root_x.data(), 0, root_x.size() - 1)) {
    //    cerr << "GUROBI: Can't get root x!" << endl;
    //    exit(1);
    //}

    //cout << endl
    //    << endl
    //    << "ROOT ENDED #########################################################"
    //    << endl
    //    << endl;

    //model_->reset();
    //model_->set(GRB_DoubleParam_NodeLimit, 9223372036800000000ll);

    model_->optimize();

    // Get final stats

    stats_.runtime_ = model_->get(GRB_DoubleAttr_Runtime);
    stats_.objective_upper_bound_ = model_->get(GRB_DoubleAttr_ObjVal);
    stats_.objective_lower_bound_ = model_->get(GRB_DoubleAttr_ObjBound);
    stats_.num_nodes_ = model_->get(GRB_DoubleAttr_NodeCount);
    //stats_.num_nodes_left_ = CPXgetnodeleftcnt(env_, model_);

    cout << fixed << setprecision(6);
    cout << "STATS:" << endl;
    cout << "{" << endl;
    cout << "\"objective_lower_bound_root\": "
        << stats_.objective_lower_bound_root_ << "," << endl;
    cout << "\"objective_upper_bound_root\": "
        << stats_.objective_upper_bound_root_ << "," << endl;
    cout << "\"runtime_root\": " << stats_.runtime_root_ << "," << endl;
    cout << "\"runtime\": " << stats_.runtime_ << "," << endl;
    cout << "\"num_nodes\": " << stats_.num_nodes_ << "," << endl;
    cout << "\"num_nodes_left\": " << stats_.num_nodes_left_ << "," << endl;
    cout << "\"num_sec\": " << stats_.num_sec_ << "," << endl;
    cout << "\"num_sec_separations\": " << stats_.num_sec_separations_ << ","
        << endl;
    cout << "\"num_sec_root\": " << stats_.num_sec_root_ << "," << endl;
    cout << "\"num_sec_separations_root\": " << stats_.num_sec_separations_root_
        << "," << endl;
    cout << "\"num_vars\": " << stats_.num_vars_ << "," << endl;
    cout << "\"num_constraints\": " << stats_.num_constraints_ << "," << endl;
    cout << "\"build_arrangement_time\": " << stats_.build_arrangement_time_
        << "," << endl;
    cout << "\"objective_lower_bound\": " << stats_.objective_lower_bound_ << ","
        << endl;
    cout << "\"objective_upper_bound\": " << stats_.objective_upper_bound_ << ","
        << endl;
    cout << "\"num_faces\": " << stats_.num_faces_ << "," << endl;
    cout << "\"num_triineq\": " << stats_.num_triineq_ << "," << endl;
    cout << "\"num_triineq_separations\": " << stats_.num_triineq_separations_
        << "," << endl;
    cout << "\"num_triineq_root\": " << stats_.num_triineq_root_ << "," << endl;
    cout << "\"num_triineq_separations_root\": "
        << stats_.num_triineq_separations_root_ << "," << endl;
    cout << "\"num_benders_separations\": " << stats_.num_benders_separations_
        << "," << endl;
    cout << "\"num_benders_separations_root\": "
        << stats_.num_benders_separations_root_ << "," << endl;
    cout << "\"num_benders_cuts\": " << stats_.num_benders_cuts_ << "," << endl;
    cout << "\"num_benders_cuts_root\": " << stats_.num_benders_cuts_root_ << ","
        << endl;
    cout << "\"num_runs_heuristic\": " << stats_.num_runs_heuristic_ << endl;

    cout << "}" << endl;

    auto vars = model_->getVars();
    for (int i = 0; i < instance_.edges_.size(); ++i)
        std::cout << vars[i].get(GRB_DoubleAttr_X) << " ";
}
