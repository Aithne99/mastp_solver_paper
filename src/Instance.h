#ifndef H_INSTANCE_
#define H_INSTANCE_

#include <utility>
#include <vector>
#include <iostream>
#include <algorithm>
#include <random>

using namespace std;

extern std::random_device rd;
extern std::mt19937 g;

struct Instance {
  int n_;
  std::vector<pair<int, int>> points_;
  vector<pair<int, int>> edges_;
  vector<int> edge_to_circle_;
  vector<vector<int>> circle_to_edges_;

  struct FaceAdjacency {
    int to_;
    int circle_;
    bool grows_;

    FaceAdjacency(int t = 0, int c = 0, bool g = false) : to_(t), circle_(c), grows_(g) {}
  };

  int num_faces_;
  std::vector<int> edge_order;
  vector<vector<FaceAdjacency>> face_graph_;
  vector<int> inclusion_set_size_;
  vector<double> face_area_;

  vector<vector<int>> hitting_set_;

  void PrintArrangementData(int argc, char* argv[])
  {
      std::cout << n_ << "\n";
      std::cout << points_.size() << "\n";
      for (auto p : points_)
      {
          std::cout << p.first << " " << p.second;
          std::cout << "\n";
      }
      std::cout << "\n";

      std::cout << edges_.size() << "\n";
      for (auto e : edges_)
      {
          std::cout << e.first << " " << e.second;
          std::cout << "\n";
      }

      std::cout << edge_to_circle_.size() << "\n";
      for (auto e : edge_to_circle_)
          std::cout << e << " ";
      std::cout << "\n";

      std::cout << circle_to_edges_.size() << "\n";
      for (auto e : circle_to_edges_)
      {
          std::cout << e.size() << "\n";
          for (auto c : e)
          {
              std::cout << c << " ";
          }
          std::cout << "\n";
      }

      std::cout << num_faces_ << "\n";

      std::cout << face_graph_.size() << "\n";
      for (auto g : face_graph_)
      {
          std::cout << g.size() << "\n";
          for (auto f : g)
          {
              std::cout << f.circle_ << " " << f.to_ << " " << f.grows_ << "\n";
          }
          std::cout << "\n";
      }

      std::cout << inclusion_set_size_.size() << "\n";
      for (auto e : inclusion_set_size_)
          std::cout << e << " ";
      std::cout << "\n";

      std::cout << face_area_.size() << "\n";
      for (auto e : face_area_)
          std::cout << e << " ";
      std::cout << "\n";

      std::cout << hitting_set_.size() << "\n";
      for (auto e : hitting_set_)
      {
          std::cout << e.size() << "\n";
          for (auto c : e)
          {
              std::cout << c << " ";
          }
          std::cout << "\n";
      }
  }
  void load(int argc, char* argv[])
  {
      std::cin >>  n_;
      int size;
      std::cin >> size;

      for (int i = 0; i < size; ++i)
      {
          int p, q;
          std::cin >> p >> q;
          points_.push_back(std::make_pair(p, q));
      }

      std::cin >> size;
      for (int i = 0; i < size; ++i)
      {
          int p, q;
          std::cin >> p >> q;
          edges_.push_back(std::make_pair(p, q));
      }

      for (int i = 0; i < size; ++i)
          edge_order.push_back(i);

      std::shuffle(edge_order.begin(), edge_order.end(), g);


      std::cin >> size;
      for (int i = 0; i < size; ++i)
      {
          int e;
          std::cin >> e;
          edge_to_circle_.push_back(e);
      }

      std::cin >> size;
      for (int i = 0; i < size; ++i)
      {
          int _size;
          int e;
          std::cin >> _size;
          std::vector<int> temp(_size);
          for (int j = 0; j < _size; ++j)
          {
              std::cin >> e;
              temp[j] = e;
          }
          circle_to_edges_.push_back(temp);
      }

      std::cin >> num_faces_;

      std::cin >> size;
      for (int i = 0; i < size; ++i)
      {
          int _size;
          int e, f;
          bool g;
          std::cin >> _size;
          std::vector<FaceAdjacency> temp(_size);
          for (int j = 0; j < _size; ++j)
          {
              std::cin >> f >> e >> g;
              temp[j] = (FaceAdjacency(e, f, g));
          }
          face_graph_.push_back(temp);
      }

      std::cin >> size;
      for (int i = 0; i < size; ++i)
      {
          int e;
          std::cin >> e;
          inclusion_set_size_.push_back(e);
      }

      std::cin >> size;
      for (int i = 0; i < size; ++i)
      {
          double e;
          std::cin >> e;
          face_area_.push_back(e);
      }

      std::cin >> size;
      for (int i = 0; i < size; ++i)
      {
          int _size;
          int e;
          std::cin >> _size;
          std::vector<int> temp(_size);
          for (int j = 0; j < _size; ++j)
          {
              std::cin >> e;
              temp[j] = e;
          }
          hitting_set_.push_back(temp);
      }
  }
};

#endif
