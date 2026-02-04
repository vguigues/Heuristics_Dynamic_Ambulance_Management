#include <fmt/core.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <fstream>
#include <functional>
#include <memory>
#include <random>
#include <unordered_map>

#include "esma/data.hpp"
#include "esma/entities.hpp"
#include "esma/experiment.hpp"
#include "esma/param.hpp"
#include "esma/travel.hpp"

template <typename AmbulanceType = esma::AmbulanceTypeUS,
          typename CallType = esma::CallTypeUS>
void run_solver_sweep(esma::Data<AmbulanceType, CallType>& data,
                      const esma::Params& params, esma::Travel& travel,
                      const std::vector<int>& n_ambs,
                      const std::string& experiment_name,
                      std::vector<std::string> solvers = {
                          "dummy_queue", "best_myopic", "non_myopic", "ghp1",
                          "ghp2", "coverage", "centrality", "ordered",
                          "district", "preparedness", "prep2",
                          "dist_centrality"}) {
  auto scenarios = esma::generate_scenarios<AmbulanceType, CallType>(
      data, params, 1, params.t_simul_begin * si::seconds,
      params.t_simul_end * si::seconds, data.regions(), data.priorities());

  std::ofstream latex_alloc(
      fmt::format("../results/tables/{}.tex", experiment_name), std::ios::out);
  if (!latex_alloc) {
    throw std::ios_base::failure(fmt::format(
        "ERROR: Could not open LaTeX output file at ../results/tables/{}.tex",
        experiment_name));
  }
  std::ofstream latex_rt(
      fmt::format("../results/tables/{}_rt.tex", experiment_name),
      std::ios::out);
  if (!latex_rt) {
    throw std::ios_base::failure(
        fmt::format("ERROR: Could not open LaTeX output file at "
                    "../results/tables/{}_rt.tex",
                    experiment_name));
  }
  latex_alloc
      << std::unitbuf;  // flush each write so lines hit disk immediately
  latex_rt << std::unitbuf;
  latex_alloc << "% Auto-generated solver sweep results\n\n";
  latex_rt << "% Auto-generated solver sweep results\n\n";

  std::vector<int> amb_counts;
  amb_counts.reserve(n_ambs.size());
  for (int n : n_ambs) {
    if (n > 0) {
      amb_counts.push_back(n);
    }
  }

  struct Metrics {
    double mean;
    double q09;
    double max;
    bool valid;
  };
  struct SolverResults {
    std::string name;
    std::vector<Metrics> alloc;
    std::vector<Metrics> resp;
  };

  auto compute_results =
      [&](const std::vector<std::string>& solver_names,
          const std::function<void(const SolverResults&)>& on_row) {
        std::vector<SolverResults> out;
        out.reserve(solver_names.size());
        for (const auto& solver_name : solver_names) {
          SolverResults entry;
          entry.name = solver_name;
          entry.alloc.reserve(amb_counts.size());
          entry.resp.reserve(amb_counts.size());
          for (size_t idx = 0; idx < amb_counts.size(); ++idx) {
            int n = amb_counts[idx];
            if (solver_name == "dist_centrality_rollout") {
              entry.alloc.push_back({0.0, 0.0, 0.0, false});
              entry.resp.push_back({0.0, 0.0, 0.0, false});
              continue;
            }
            esma::Params params_n = params;
            params_n.nb_ambulances = n;
            esma::Data<AmbulanceType, CallType> data_n(params_n, travel);
            auto ambulances = data_n.ambulances();  // sized to current n
            esma::Simulator<AmbulanceType, CallType> sim(
                params_n, data_n, params.log_level, travel,
                params.t_simul_begin * si::seconds, solver_name, solver_name);
            sim.set_log_level(0);
            fmt::print("Begun {}\n", solver_name);
            esma::SimulationResult result = sim.run(
                params.t_simul_begin * si::seconds, scenarios[0], ambulances);
            result.write_results(
                fmt::format("../results/response_times/{}_{}_a{}.txt",
                            experiment_name, solver_name, n));
            result.write_summary(fmt::format("../results/summary_{}_{}_a{}.txt",
                                             experiment_name, solver_name, n));
            entry.alloc.push_back({result.mean_penalized_response_time(),
                                   result.q09_penalized_response_time(),
                                   result.max_penalized_response_time(), true});
            entry.resp.push_back({result.mean_response_time(),
                                  result.q09_response_time(),
                                  result.max_response_time(), true});
            fmt::print("Finished {}, a = {}\n", solver_name, n);
          }
          on_row(entry);
          out.push_back(std::move(entry));
        }
        return out;
      };

  std::vector<std::string> rollout_solvers;
  rollout_solvers.reserve(solvers.size());
  for (const auto& solver : solvers) {
    rollout_solvers.push_back(fmt::format("{}_rollout", solver));
  }
  std::vector<SolverResults> rollout_results =
      compute_results(rollout_solvers, [](const SolverResults&) {});

  std::unordered_map<int, size_t> amb_index;
  amb_index.reserve(amb_counts.size());
  for (size_t i = 0; i < amb_counts.size(); ++i) {
    amb_index[amb_counts[i]] = i;
  }

  struct TableStream {
    std::ofstream& latex;
    const std::vector<int>& counts;
    const std::unordered_map<int, size_t>& amb_index;
    std::string caption;
    std::string scale;
    bool use_response;

    void begin() {
      latex << "\\begin{table}[h!]\n";
      latex << fmt::format("\\scalebox{{{}}}{{\n", scale);
      latex << "\\centering\n";
      latex << "\\begin{tabular}{@{}l";
      for (size_t i = 0; i < counts.size(); ++i) {
        latex << "|ccc";
      }
      latex << "@{}}\n\\toprule\n";
      latex << "Number of Ambulances";
      for (size_t i = 0; i < counts.size(); ++i) {
        const char* bar = (i + 1 == counts.size()) ? "" : "|";
        latex << fmt::format(" & \\multicolumn{{3}}{{c{}}}{{{}}}", bar,
                             counts[i]);
      }
      latex << " \\\\\n\\midrule\n";
      latex << "Heuristic";
      for (size_t i = 0; i < counts.size(); ++i) {
        latex << " & Mean & Q0.9 & Max";
      }
      latex << " \\\\\n\\midrule\n";
    }

    void write_row(const SolverResults& row) {
      latex << row.name;
      for (int n : counts) {
        const auto& metrics = use_response ? row.resp[amb_index.at(n)]
                                           : row.alloc[amb_index.at(n)];
        if (!metrics.valid) {
          latex << " & - & - & -";
        } else {
          latex << fmt::format(" & {:.0f} & {:.0f} & {:.0f}", metrics.mean,
                               metrics.q09, metrics.max);
        }
      }
      latex << " \\\\\n";
    }

    void end() {
      latex << "\\bottomrule\n\\end{tabular}\n}\n";
      latex << fmt::format("\\caption{{{}}}\n", caption);
      latex << "\\end{table}\n\n";
    }
  };

  TableStream alloc_base{
      latex_alloc, amb_counts,
      amb_index,   "Allocation costs by solver and fleet size",
      "1",         false};
  TableStream rt_base{latex_rt,  amb_counts,
                      amb_index, "Response times by heuristic and fleet size",
                      "0.7",     true};
  TableStream alloc_rollout{
      latex_alloc, amb_counts,
      amb_index,   "Allocation costs by solver rollout and fleet size",
      "1",         false};
  TableStream rt_rollout{
      latex_rt,  amb_counts,
      amb_index, "Response times by heuristic rollout and fleet size",
      "0.7",     true};

  alloc_base.begin();
  rt_base.begin();
  compute_results(solvers, [&](const SolverResults& row) {
    alloc_base.write_row(row);
    rt_base.write_row(row);
  });
  alloc_base.end();
  rt_base.end();

  alloc_rollout.begin();
  rt_rollout.begin();
  compute_results(rollout_solvers, [&](const SolverResults& row) {
    alloc_rollout.write_row(row);
    rt_rollout.write_row(row);
  });
  alloc_rollout.end();
  rt_rollout.end();

  latex_alloc.close();
  latex_rt.close();
}

template <typename AmbulanceType = esma::AmbulanceTypeUS,
          typename CallType = esma::CallTypeUS>
void run_solver_sweep_no_rollout(esma::Data<AmbulanceType, CallType>& data,
                                 const esma::Params& params,
                                 esma::Travel& travel,
                                 const std::vector<int>& n_ambs,
                                 const std::string& experiment_name,
                                 std::vector<std::string> solvers = {
                                     "dummy_queue", "best_myopic", "non_myopic",
                                     "ghp1", "ghp2", "coverage", "centrality",
                                     "ordered", "district", "preparedness",
                                     "prep2", "dist_centrality"}) {
  auto scenarios = esma::generate_scenarios<AmbulanceType, CallType>(
      data, params, 1, params.t_simul_begin * si::seconds,
      params.t_simul_end * si::seconds, data.regions(), data.priorities());

  std::ofstream latex_alloc(
      fmt::format("../results/tables/{}.tex", experiment_name), std::ios::out);
  if (!latex_alloc) {
    throw std::ios_base::failure(fmt::format(
        "ERROR: Could not open LaTeX output file at ../results/tables/{}.tex",
        experiment_name));
  }
  std::ofstream latex_rt(
      fmt::format("../results/tables/{}_rt.tex", experiment_name),
      std::ios::out);
  if (!latex_rt) {
    throw std::ios_base::failure(
        fmt::format("ERROR: Could not open LaTeX output file at "
                    "../results/tables/{}_rt.tex",
                    experiment_name));
  }
  latex_alloc
      << std::unitbuf;  // flush each write so lines hit disk immediately
  latex_rt << std::unitbuf;
  latex_alloc << "% Auto-generated solver sweep results\n\n";
  latex_rt << "% Auto-generated solver sweep results\n\n";

  std::vector<int> amb_counts;
  amb_counts.reserve(n_ambs.size());
  for (int n : n_ambs) {
    if (n > 0) {
      amb_counts.push_back(n);
    }
  }

  struct Metrics {
    double mean;
    double q09;
    double max;
    bool valid;
  };
  struct SolverResults {
    std::string name;
    std::vector<Metrics> alloc;
    std::vector<Metrics> resp;
  };

  std::vector<SolverResults> results;
  results.reserve(solvers.size());
  for (const auto& solver_name : solvers) {
    SolverResults entry;
    entry.name = solver_name;
    entry.alloc.reserve(amb_counts.size());
    entry.resp.reserve(amb_counts.size());
    for (size_t idx = 0; idx < amb_counts.size(); ++idx) {
      int n = amb_counts[idx];
      if (solver_name == "dist_centrality_rollout") {
        entry.alloc.push_back({0.0, 0.0, 0.0, false});
        entry.resp.push_back({0.0, 0.0, 0.0, false});
        continue;
      }
      esma::Params params_n = params;
      params_n.nb_ambulances = n;
      esma::Data<AmbulanceType, CallType> data_n(params_n, travel);
      auto ambulances = data_n.ambulances();
      esma::Simulator<AmbulanceType, CallType> sim(
          params_n, data_n, params.log_level, travel,
          params.t_simul_begin * si::seconds, solver_name, solver_name);
      sim.set_log_level(0);
      fmt::print("Begun {}\n", solver_name);
      esma::SimulationResult result =
          sim.run(params.t_simul_begin * si::seconds, scenarios[0], ambulances);
      result.write_results(
          fmt::format("../results/response_times/{}_{}_a{}.txt",
                      experiment_name, solver_name, n));
      result.write_summary(fmt::format("../results/summary_{}_{}_a{}.txt",
                                       experiment_name, solver_name, n));
      entry.alloc.push_back({result.mean_penalized_response_time(),
                             result.q09_penalized_response_time(),
                             result.max_penalized_response_time(), true});
      entry.resp.push_back({result.mean_response_time(),
                            result.q09_response_time(),
                            result.max_response_time(), true});
      fmt::print("Finished {}, a = {}\n", solver_name, n);
    }
    results.push_back(std::move(entry));
  }

  const size_t split_index = amb_counts.size() > 5 ? 5 : amb_counts.size();
  std::vector<int> amb_left(amb_counts.begin(),
                            amb_counts.begin() + split_index);
  std::vector<int> amb_right(amb_counts.begin() + split_index,
                             amb_counts.end());
  std::unordered_map<int, size_t> amb_index;
  amb_index.reserve(amb_counts.size());
  for (size_t i = 0; i < amb_counts.size(); ++i) {
    amb_index[amb_counts[i]] = i;
  }

  auto write_table = [&](std::ofstream& latex, const std::string& caption,
                         const std::string& scale, bool add_hspace,
                         bool use_response) {
    auto write_block = [&](const std::vector<int>& counts, bool with_hspace,
                           bool add_center) {
      if (counts.empty()) {
        return;
      }
      latex << fmt::format("\\scalebox{{{}}}{{\n", scale);
      if (add_center) {
        latex << "\\centering\n";
      }
      if (with_hspace) {
        latex << "\\hspace{2cm}\n";
      }
      latex << "\\begin{tabular}{@{}l";
      for (size_t i = 0; i < counts.size(); ++i) {
        latex << "|ccc";
      }
      latex << "@{}}\n\\toprule\n";
      latex << "Number of Ambulances";
      for (size_t i = 0; i < counts.size(); ++i) {
        const char* bar = (i + 1 == counts.size()) ? "" : "|";
        latex << fmt::format(" & \\multicolumn{{3}}{{c{}}}{{{}}}", bar,
                             counts[i]);
      }
      latex << " \\\\\n\\midrule\n";
      latex << "Heuristic";
      for (size_t i = 0; i < counts.size(); ++i) {
        latex << " & Mean & Q0.9 & Max";
      }
      latex << " \\\\\n\\midrule\n";
      for (const auto& row : results) {
        latex << row.name;
        for (int n : counts) {
          const auto& metrics =
              use_response ? row.resp[amb_index[n]] : row.alloc[amb_index[n]];
          if (!metrics.valid) {
            latex << " & - & - & -";
          } else {
            latex << fmt::format(" & {:.0f} & {:.0f} & {:.0f}", metrics.mean,
                                 metrics.q09, metrics.max);
          }
        }
        latex << " \\\\\n";
      }
      latex << "\\bottomrule\n\\end{tabular}\n}\n";
    };

    latex << "\\begin{table}[h!]\n";
    write_block(amb_left, false, true);
    if (!amb_right.empty()) {
      write_block(amb_right, add_hspace, false);
    }
    latex << fmt::format("\\caption{{{}}}\n", caption);
    latex << "\\end{table}\n\n";
  };

  write_table(latex_alloc, "Allocation costs by solver and fleet size", "1",
              true, false);
  write_table(latex_rt, "Response times by heuristic and fleet size", "0.7",
              false, true);

  latex_alloc.close();
  latex_rt.close();
}

template <typename AmbulanceType = esma::AmbulanceTypeUS,
          typename CallType = esma::CallTypeUS>
void run_solver_sweep_random(esma::Data<AmbulanceType, CallType>& data,
                             const esma::Params& params, esma::Travel& travel,
                             const std::vector<int>& n_ambs,
                             const std::string& experiment_name,
                             std::vector<std::string> solvers = {
                                 "dummy_queue", "best_myopic", "non_myopic",
                                 "ghp1", "ghp2", "coverage", "centrality",
                                 "ordered", "district", "preparedness", "prep2",
                                 "dist_centrality"}) {
  (void)data;
  (void)travel;
  esma::Params random_params = params;
  random_params.time_random = true;

  esma::Travel random_travel(false);
  random_travel.set_params(random_params.intercept_time * si::seconds,
                           random_params.M, random_params.delta,
                           random_params.lambda_time);
  esma::Data<AmbulanceType, CallType> random_data(random_params, random_travel);

  run_solver_sweep(random_data, random_params, random_travel, n_ambs,
                   experiment_name, std::move(solvers));
}

template <typename AmbulanceType = esma::AmbulanceTypeUS,
          typename CallType = esma::CallTypeUS>
void test_toy_one_stage(esma::Data<AmbulanceType, CallType>& data,
                        const esma::Params& params, esma::Travel& travel) {
  std::vector<esma::Call<esma::CallTypeUS>> scenario =
      esma::generate_toy_scenario<esma::AmbulanceTypeUS, esma::CallTypeUS>(
          data, params);
  for (auto& call : scenario) {
    fmt::print("call {}, t {}, p {}, l {} {} h {} {} {} b {} {} {}\n",
               call.id(), call.time().value(), call.priority_index(),
               call.location().latitude, call.location().longitude,
               call.hospital().id(), call.hospital().location().latitude,
               call.hospital().location().longitude, call.cleaning_base().id(),
               call.cleaning_base().location().latitude,
               call.cleaning_base().location().longitude);
  }
  std::string dispatch_solver = "dummy_queue";
  std::string reassignment_solver = "dummy_queue";
  esma::Simulator<esma::AmbulanceTypeUS, esma::CallTypeUS> sim(
      params, data, travel, params.t_simul_begin * si::seconds, dispatch_solver,
      reassignment_solver);

  esma::SimulationResult result =
      sim.run(0 * si::seconds, scenario, data.ambulances());
  result.write_results("../results.txt");
  sim.write_trajectories("../trajectories.txt");
}

template <typename AmbulanceType = esma::AmbulanceTypeUS,
          typename CallType = esma::CallTypeUS>
void test_two_stage(esma::Data<AmbulanceType, CallType>& data,
                    const esma::Params& params, esma::Travel& travel,
                    std::string policy) {
  std::string dispatch_solver = policy;
  std::string reassignment_solver = policy;
  auto t0 = std::chrono::high_resolution_clock::now();
  esma::Simulator<esma::AmbulanceTypeUS, esma::CallTypeUS> sim_one_stage(
      params, data, params.log_level, travel,
      params.t_simul_begin * si::seconds, dispatch_solver, reassignment_solver);
  auto scenarios = esma::generate_scenarios<AmbulanceType, CallType>(
      data, params, 6, params.t_simul_begin * si::seconds,
      params.t_simul_end * si::seconds, data.regions(), data.priorities());
  auto amb_copy = data.ambulances();
  sim_one_stage.set_log_level(2);
  esma::SimulationResult result_one_stage =
      sim_one_stage.run(0 * si::seconds, scenarios[0], amb_copy);
  auto dt = std::chrono::high_resolution_clock::now() - t0;
  auto time_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(dt).count();
  fmt::print("Finished one stage in {} ms\n", time_ms);
  std::cin.get();
  result_one_stage.write_results(
      fmt::format("../results/results_one_stage_{}.txt", policy));
  result_one_stage.write_summary(
      fmt::format("summary_one_stage_{}.txt", policy));
  sim_one_stage.write_trajectories(
      fmt::format("../trajectories/trajectories_one_stage_{}.txt", policy));
  dispatch_solver = policy + "_rollout";
  reassignment_solver = policy + "_rollout";
  esma::Simulator<esma::AmbulanceTypeUS, esma::CallTypeUS> sim(
      params, data, params.log_level, travel,
      params.t_simul_begin * si::seconds, dispatch_solver, reassignment_solver);
  fmt::print("Running two stage\n");
  //   std::cin.get();
  sim.set_log_level(2);
  amb_copy = data.ambulances();
  t0 = std::chrono::high_resolution_clock::now();
  esma::SimulationResult result =
      sim.run(0 * si::seconds, scenarios[0], amb_copy,
              std::vector<esma::Call<CallType>>(), std::vector<int>(),
              std::vector<int>());
  dt = std::chrono::high_resolution_clock::now() - t0;
  double time_s =
      static_cast<double>(
          std::chrono::duration_cast<std::chrono::milliseconds>(dt).count()) /
      1000.0;
  result.write_results(fmt::format("../results_two_stage_{}.txt", policy));
  result.write_summary(
      fmt::format("summary_two_stage_{}.txt", dispatch_solver));
  sim.write_trajectories(
      fmt::format("../trajectories_two_stage_{}.txt", policy));

  fmt::print("Policy {} finished. Mean one_stage = {}, Mean two_stage = {}\n",
             policy, result_one_stage.mean_penalized_response_time(),
             result.mean_penalized_response_time());
  fmt::print("Time one_stage = {} ms, time two_stage = {} s\n", time_ms,
             time_s);
}

template <typename AmbulanceType = esma::AmbulanceTypeUS,
          typename CallType = esma::CallTypeUS>
void test_two_stage_random(esma::Data<AmbulanceType, CallType>& data,
                           const esma::Params& params, esma::Travel& travel,
                           std::string policy) {
  (void)data;
  (void)travel;
  esma::Params random_params = params;
  random_params.time_random = true;

  esma::Travel random_travel(false);
  random_travel.set_params(random_params.intercept_time * si::seconds,
                           random_params.M, random_params.delta,
                           random_params.lambda_time);
  esma::Data<AmbulanceType, CallType> random_data(random_params, random_travel);

  std::string dispatch_solver = policy;
  std::string reassignment_solver = policy;
  auto t0 = std::chrono::high_resolution_clock::now();
  esma::Simulator<esma::AmbulanceTypeUS, esma::CallTypeUS> sim_one_stage(
      random_params, random_data, random_params.log_level, random_travel,
      random_params.t_simul_begin * si::seconds, dispatch_solver,
      reassignment_solver);
  auto scenarios = esma::generate_scenarios<AmbulanceType, CallType>(
      random_data, random_params, 6, random_params.t_simul_begin * si::seconds,
      random_params.t_simul_end * si::seconds, random_data.regions(),
      random_data.priorities());
  auto amb_copy = random_data.ambulances();
  sim_one_stage.set_log_level(2);
  esma::SimulationResult result_one_stage =
      sim_one_stage.run(0 * si::seconds, scenarios[0], amb_copy);
  auto dt = std::chrono::high_resolution_clock::now() - t0;
  auto time_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(dt).count();
  fmt::print("Finished one stage in {} ms\n", time_ms);
  std::cin.get();
  result_one_stage.write_results(
      fmt::format("../results/results_one_stage_{}.txt", policy));
  result_one_stage.write_summary(
      fmt::format("summary_one_stage_{}.txt", policy));
  sim_one_stage.write_trajectories(
      fmt::format("../trajectories/trajectories_one_stage_{}.txt", policy));
  dispatch_solver = policy + "_rollout";
  reassignment_solver = policy + "_rollout";
  esma::Simulator<esma::AmbulanceTypeUS, esma::CallTypeUS> sim(
      random_params, random_data, random_params.log_level, random_travel,
      random_params.t_simul_begin * si::seconds, dispatch_solver,
      reassignment_solver);
  fmt::print("Running two stage\n");
  sim.set_log_level(4);

  amb_copy = random_data.ambulances();
  t0 = std::chrono::high_resolution_clock::now();
  esma::SimulationResult result =
      sim.run(0 * si::seconds, scenarios[0], amb_copy,
              std::vector<esma::Call<CallType>>(), std::vector<int>(),
              std::vector<int>());
  dt = std::chrono::high_resolution_clock::now() - t0;
  double time_s =
      static_cast<double>(
          std::chrono::duration_cast<std::chrono::milliseconds>(dt).count()) /
      1000.0;
  result.write_results(fmt::format("../results_two_stage_{}.txt", policy));
  result.write_summary(
      fmt::format("summary_two_stage_{}.txt", dispatch_solver));
  sim.write_trajectories(
      fmt::format("../trajectories_two_stage_{}.txt", policy));

  fmt::print("Policy {} finished. Mean one_stage = {}, Mean two_stage = {}\n",
             policy, result_one_stage.mean_penalized_response_time(),
             result.mean_penalized_response_time());
  fmt::print("Time one_stage = {} ms, time two_stage = {} s\n", time_ms,
             time_s);
}

template <typename AmbulanceType = esma::AmbulanceTypeUS,
          typename CallType = esma::CallTypeUS>
void test_quick_sweep(esma::Data<AmbulanceType, CallType>& data,
                      const esma::Params& params, esma::Travel& travel,
                      std::vector<int> amb_counts) {
  (void)data;
  (void)travel;
  esma::Params random_params = params;
  random_params.time_random = true;

  esma::Travel random_travel(false);
  random_travel.set_params(random_params.intercept_time * si::seconds,
                           random_params.M, random_params.delta,
                           random_params.lambda_time);
  esma::Data<AmbulanceType, CallType> random_data(random_params, random_travel);

  std::vector<std::string> solvers = {"prep2"};
  run_solver_sweep(random_data, random_params, random_travel, amb_counts,
                   "quick_sweep", solvers);
}

template <typename AmbulanceType = esma::AmbulanceTypeUS,
          typename CallType = esma::CallTypeUS>
void compare_travel_time_difference(esma::Data<AmbulanceType, CallType>& data,
                                    const esma::Params& params,
                                    size_t sample_count = 1000) {
  auto ambulances = data.ambulances();
  if (ambulances.empty()) {
    fmt::print("No ambulances available to compute average speed.\n");
    return;
  }

  bu_velocity_t avg_speed = 0.0 * si::meters_per_second;
  for (const auto& amb : ambulances) {
    avg_speed += amb.speed();
  }
  avg_speed /= static_cast<double>(ambulances.size());

  std::vector<esma::Location> coords;
  for (const auto& region_samples : data.region_samples()) {
    coords.insert(coords.end(), region_samples.begin(), region_samples.end());
  }
  if (coords.size() < sample_count) {
    const auto& centers = data.centers();
    coords.insert(coords.end(), centers.begin(), centers.end());
    for (const auto& base : data.waiting_bases()) {
      coords.push_back(base.location());
    }
    for (const auto& hospital : data.hospitals()) {
      coords.push_back(hospital.location());
    }
    for (const auto& cleaning : data.cleaning_bases()) {
      coords.push_back(cleaning.location());
    }
  }

  if (coords.empty()) {
    fmt::print("No coordinates available to compare travel times.\n");
    return;
  }

  std::mt19937 rng(42);
  std::shuffle(coords.begin(), coords.end(), rng);
  std::vector<std::pair<esma::Location, esma::Location>> pairs;
  pairs.reserve(sample_count);
  if (coords.size() >= 2) {
    for (size_t i = 0; i + 1 < coords.size() && pairs.size() < sample_count;
         i += 2) {
      if (!(coords[i] == coords[i + 1])) {
        pairs.emplace_back(coords[i], coords[i + 1]);
      }
    }
  }

  size_t used = pairs.size();
  if (used == 0) {
    fmt::print("No valid coordinate pairs available.\n");
    return;
  }
  if (used < sample_count) {
    fmt::print("Only {} coordinate pairs available; using all of them.\n",
               used);
  }

  esma::Travel deterministic_travel(true);
  esma::Travel random_travel(false);
  random_travel.set_params(params.intercept_time * si::seconds, params.M,
                           params.delta, params.lambda_time);

  double total_det = 0.0;
  double total_rand = 0.0;
  double total_diff = 0.0;
  // Print random params:
  fmt::print(
      "Random travel params: intercept_time = {}s, M = {}, delta = "
      "{}, lambda_time = {}, used = {}, coords.size() = {}\n",
      params.intercept_time, params.M, params.delta, params.lambda_time, used,
      coords.size());

  for (size_t i = 0; i < used; ++i) {
    random_travel.set_seed(static_cast<unsigned int>(i + 1));
    const auto& source = pairs[i].first;
    const auto& target = pairs[i].second;
    const auto det_time =
        deterministic_travel.travel_time(source, target, avg_speed).value();
    const auto rand_time =
        random_travel.travel_time(source, target, avg_speed).value();
    const double diff = rand_time - det_time;
    total_det += det_time;
    total_rand += rand_time;
    total_diff += diff;
  }

  const double denom = static_cast<double>(used);
  fmt::print("Travel time comparison over {} coordinate pairs:\n", used);
  fmt::print("mean_det = {:.2f}s, mean_rand = {:.2f}s, mean_diff = {:.2f}s\n",
             total_det / denom, total_rand / denom, total_diff / denom);
}

int main(int argc, char const* argv[]) {
  fmt::print(
      "Ambulance Management Simulator v0.1 - Developed by Victor Hugo "
      "Nascimento\n");
  const esma::Params params = esma::load_params(argc, argv);
  esma::print_params(params);

  esma::Travel travel(!params.time_random);
  esma::Data<esma::AmbulanceTypeUS, esma::CallTypeUS> data(params, travel);
  fmt::print("data.regions.size() = {}\n", data.regions().size());
  // test_toy_one_stage<esma::AmbulanceTypeUS, esma::CallTypeUS>(data, params,
  // travel);
  // compare_travel_time_difference(data, params, 1000);
  // std::cin.get();

  // test_quick_sweep(data, params, travel, {16});
  // test_two_stage<esma::AmbulanceTypeUS, esma::CallTypeUS>(
  //     data, params, travel, "markov_preparedness");

  // test_two_stage_random<esma::AmbulanceTypeUS, esma::CallTypeUS>(
  //     data, params, travel, "prep2");

  std::vector<int> amb_counts = {10, 12, 14, 16, 18};
  std::string sweep_out = "solver_sweep.tex";
  if (!amb_counts.empty()) {
    int first = amb_counts.front();
    int last = amb_counts.back();
    sweep_out = fmt::format("no_rollout_street_{}_{}", first, last);
  }

  run_solver_sweep_no_rollout<esma::AmbulanceTypeUS, esma::CallTypeUS>(
      data, params, travel, amb_counts, sweep_out,
      {"dummy_queue", "coverage", "centrality", "ordered", "preparedness",
       "prep2", "dist_centrality", "tipat", "markov_preparedness"});

  // std::vector<int> amb_counts = {26, 28, 30};
  // std::string sweep_out = "";
  // if (!amb_counts.empty()) {
  //   int first = amb_counts.front();
  //   int last = amb_counts.back();
  //   sweep_out = fmt::format("det_mp_tipat_{}_{}", first, last);
  // }
  // run_solver_sweep<esma::AmbulanceTypeUS, esma::CallTypeUS>(
  //     data, params, travel, amb_counts, sweep_out,
  //     {"markov_preparedness", "tipat"});

  // if (!amb_counts.empty()) {
  //   int first = amb_counts.front();
  //   int last = amb_counts.back();
  //   sweep_out = fmt::format("nm_rerun_random_{}_{}", first, last);
  // }

  // run_solver_sweep_random<esma::AmbulanceTypeUS, esma::CallTypeUS>(
  //     data, params, travel, amb_counts, sweep_out, {"prep2"});
  return 0;
}
