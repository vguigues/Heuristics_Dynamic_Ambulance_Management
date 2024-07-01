#ifndef _POLICY_TESTER_H
#define _POLICY_TESTER_H

#include "instance.h"
#include "main.h"
#include "solver.h"

class Instance;
class OSRMHelper;
class Travel;
class Solver;
class DistrictData;

using namespace std;
struct SecondStageWaitingTime {
  double waiting_on_scene;
  double waiting_to_hospital;
  double mean_total_time;
};

struct Stats {
  double mean_waiting_on_scene;
  double mean_waiting_on_scene_penalized;
  double mean_waiting_to_hospital;
  double max_waiting_on_scene;
  double max_waiting_on_scene_penalized;
  double max_waiting_to_hospital;
  double actual_std_dev;
  double actual_std_dev_avg;
  double std_dev;
  double std_dev_avg;
  double waiting_on_scene_q90;
  double waiting_to_hospital_q90;
  double mean_total;
  double max_total;
  double q90_total;

  Stats(vector<double> waiting_on_scene,
        vector<double> waiting_on_scene_penalized,
        vector<double> waiting_to_hospital);
};

class PolicyTester {
 public:
  PolicyTester(Instance &ins);
  ~PolicyTester();

  const string real_data1 =
      "calibration/Scenarios_For_Tests/"
      "Real_Data_Continuous_Scenarios/baseScenario.txt";

  GRBEnv env;
  Instance &ins;
  Travel &travel;
  vector<Ambulance> ambulances;
  double time;

  vector<vector<int>> block_indexes;
  vector<vector<double>> waiting_on_scene;
  vector<vector<double>> waiting_on_scene_penalized;
  vector<vector<double>> waiting_to_hospital;
  vector<vector<double>> calls_end;
  vector<vector<int>> which_ambulance;

  vector<vector<double>> lambda;
  vector<DistrictData> data;
  vector<double> run_times;
  vector<tuple<int, int, double>> solver_run_times;

  // const vector<string> policies{"queue", "forward", "priorities", "minmax",
  // 	"gen_forward", "gen_minmax","non_miopyc"};

  // const vector<string> policies{"non_miopyc"};
  // const vector<string> policies{"queue", "markov_prep", "preparedness",
  //                               "district", "ordered"};
  // const vector<string> policies{"cg", "enumerate"};
  const vector<string> policies{"forward", "priorities", "minmax",
                                "non_miopyc"};
  // const vector<string> policies{"preparedness", "district", "ordered"};
  // const vector<string> policies{"prep2"};
  // const vector<string> policies{"markov_prep"};
  const map<string, string> policies_names{
      {"queue", "CA"},         {"forward", "BM"},  {"forward_prep", "BMP"},
      {"priorities", "GHP1"},  {"minmax", "GHP2"}, {"non_miopyc", "NM"},
      {"preparedness", "PR1"}, {"prep2", "PR2"},   {"district", "DS"},
      {"ordered", "ORD"},      {"coverage", "CV"}, {"markov_prep", "MP"},
      {"queue_deficit", "QD"}};

  shared_ptr<Solver> get_solver(const string &policy, vector<Call> &calls,
                                vector<Ambulance> &ambulances, Travel &travel,
                                int g = 0, double time = 0.0);

  void run();

  void one_stage_old();
  void one_stage();
  void two_stage();
  void two_stage_tree();
  void two_stage_tree_new();

  double penalized_response_time(double response_time, double amb_type,
                                 double call_priority) {
    if (g_params.amb_setup == "us") {
      if (amb_type == -1) {
        return response_time * ins.penalties[call_priority];
      }
      return response_time * ins.penalties[call_priority] +
             1000 * ins.penalty_matrix[amb_type][call_priority];
    } else if (g_params.amb_setup == "rj") {
      return response_time * ins.penalties[call_priority];
    } else {
      fmt::print("Unkown amb setup: {}\n", g_params.amb_setup);
      exit(1);
    }
    return response_time * ins.penalties[call_priority];
  }
  vector<vector<vector<Call>>> load_scenario_tree();

  int set_next_event(vector<Call> &this_scenario, int &event_call,
                     int &index_call, vector<int> &queue);

  SecondStageWaitingTime get_waiting_time(int amb_id, Call &call,
                                          const string &policy,
                                          vector<int> &index_scenarios,
                                          int nearest_base_id,
                                          std::vector<int> queue, int sc,
                                          double time);

  vector<vector<Call>> get_future_scenarios(
      double time, vector<Call> &this_scenario, int nb_realizations, int T,
      vector<vector<vector<Call>>> &my_scenarios);

  SecondStageWaitingTime get_waiting_time_tree(
      int amb_id, int call_id, const string &policy,
      vector<vector<Call>> &future_scenarios, vector<Call> &this_scenario,
      double time, std::vector<int> &queue, int base);

  SecondStageWaitingTime get_waiting_time_tree_return(
      int amb_id, int base, const string &policy,
      vector<vector<Call>> &future_scenarios, vector<Call> &this_scenario,
      double time, std::vector<int> &queue, double min_time_to_free = 0);

  SecondStageWaitingTime get_waiting_time_return(
      int amb_id, const string &policy, vector<int> &index_scenarios,
      double time, std::vector<int> queue, int sc, int base);

  vector<int> get_nearest_base(vector<Call> &calls);

  bool can_answer(Ambulance &amb, Call &call);

  SecondStageWaitingTime get_waiting_time_model(
      int amb_id, int call_id, const string &policy,
      vector<vector<Call>> &future_scenarios, vector<Call> &this_scenario,
      double time, std::vector<int> &queue, int base);

  SecondStageWaitingTime get_waiting_time_model_return(
      int amb_id, int base, const string &policy,
      vector<vector<Call>> &future_scenarios, vector<Call> &this_scenario,
      double time, std::vector<int> &queue);

  double run_model(vector<Call> &scenario, int amb_id, int call_id, int base_id,
                   bool is_return);
};

#endif