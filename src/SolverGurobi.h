#ifndef H_SOLVER_GUROBI_
#define H_SOLVER_GUROBI_

#include <vector>
#include <set>

#ifdef WIN32
#include "C:\gurobi912\win64\include\gurobi_c++.h"
#else
#include "gurobi_c++.h"
#endif

#include "Instance.h"
#include "SEC.h"
#include "Benders.h"

using namespace std;

struct SolverGurobi {
  Instance& instance_;

  GRBEnv* env_ = nullptr;
  GRBModel* model_ = nullptr;

  SolverGurobi(Instance& instance) : instance_(instance) {}

  void BuildModel();

  void Run();

  
  //todo: make this function private
  void DfsComputeY(vector<bool>& visited, int face, int count, double* x);
  double DfsComputeObj(vector<bool>& visited, int face, int count, double* x);
  vector<int> ClosestTree(double* x);
  void ImproveTree(vector<int>& tree);
  
  //Stats
  struct {
    double objective_lower_bound_root_;
    double objective_upper_bound_root_; 
    double runtime_root_;
    double runtime_;
    int num_nodes_;
    int num_nodes_left_;
    int num_sec_;
    int num_sec_separations_;
    int num_sec_root_;
    int num_sec_separations_root_;
    int num_vars_;
    int num_constraints_;
    double build_arrangement_time_;
    double objective_lower_bound_;
    double objective_upper_bound_;
    int num_faces_;
    int num_triineq_;
    int num_triineq_separations_;
    int num_triineq_root_;
    int num_triineq_separations_root_;
    int num_runs_heuristic_;
    int num_benders_separations_;
    int num_benders_separations_root_;
    int num_benders_cuts_;
    int num_benders_cuts_root_;
  } stats_;

  //
  double benders_early_stop_ub;

  int ArcIndex(const int root, const int edge, const bool reversed) const;

  private:

  void AddCoveringConstraintsDFS(vector<bool>& visited, const int face,
                                set<int>& circles);

  void AddReducedCoveringConstraintsDFS(vector<bool>& visited, const int face,
                                set<int>& circles);

  void AddCoveringConstraint(const int face, const int edge);
  void AddCardinalityConstraint();
  void AddImplicationConstraint(const int a, const int b);


  void AddArcEdgeCouplingConstraints();

  void SetSolverParameters();
};

class CutCallback : public GRBCallback
{
public:
    CutCallback(SolverGurobi& _solver) : solver(_solver) {}
    void callback();
    void AddSEC(const vector<int>& S);
    void AddBendersCut(const vector<double>& opt_cut);
    int SeparateTriangleInequalities(const vector<double>& x);
    void AddCutset(const int root, const vector<int>& arcs);
private:
    SolverGurobi& solver;
    //GRBVar* vars;
    //int numVars;
};
#endif