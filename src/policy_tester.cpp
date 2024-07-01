#include "../include/policy_tester.h"

#include "../include/cg.h"
#include "../include/data.h"
#include "../include/osrm_helper.h"
#include "../include/solver.h"
#include "../include/travel.h"

PolicyTester::PolicyTester(Instance &ins)
    : ins(ins),
      travel(ins.travel),
      ambulances(ins.ambulances),
      time(0),
      waiting_on_scene(
          vector<vector<double>>(ins.nb_scenarios, vector<double>())),
      waiting_on_scene_penalized(
          vector<vector<double>>(ins.nb_scenarios, vector<double>())),
      waiting_to_hospital(
          vector<vector<double>>(ins.nb_scenarios, vector<double>())),
      calls_end(vector<vector<double>>(ins.nb_scenarios, vector<double>())),
      which_ambulance(vector<vector<int>>(ins.nb_scenarios, vector<int>())),
      lambda(7, vector<double>()) {
  // for(int i = 0; i < ins.nb_scenarios; ++i){
  // 	waiting_on_scene[i] = vector<double>(ins.nb_calls[i], GRB_INFINITY);
  // 	waiting_on_scene_penalized[i] = vector<double>(ins.nb_calls[i],
  // GRB_INFINITY); 	waiting_to_hospital[i] = vector<double>(ins.nb_calls[i],
  // GRB_INFINITY); 	calls_end[i] = vector<double>(ins.nb_calls[i],
  // GRB_INFINITY); 	which_ambulance[i] = vector<int>(ins.nb_calls[i], -1);
  // }

  if (find(policies.begin(), policies.end(), "district") != policies.end()) {
    for (int g = 0; g < 7; ++g) {
      fmt::print("g = {}\n", g);
      DistrictManager dm(ins, travel, g, ambulances);
      lambda[g] = dm.lambda;
      data.push_back(dm.get_districts());
    }
  }
}

shared_ptr<Solver> PolicyTester::get_solver(const string &policy,
                                            vector<Call> &calls,
                                            vector<Ambulance> &ambulances,
                                            Travel &travel, int g,
                                            double time) {
  if (policy == "forward") {
    return static_pointer_cast<Solver>(
        make_shared<ForwardSolver>(env, calls, ambulances, ins, travel, time));
  } else if (policy == "queue") {
    return static_pointer_cast<Solver>(
        make_shared<QueueSolver>(env, calls, ambulances, ins, travel, time));
  } else if (policy == "priorities" || policy == "priorities_a") {
    auto solver = static_pointer_cast<Solver>(
        make_shared<PrioritySolver>(env, calls, ambulances, ins, travel, time));
    solver->travel.set_forward(policy == "priorities");
    return solver;
  } else if (policy == "minmax") {
    return static_pointer_cast<Solver>(
        make_shared<MinMaxSolver>(env, calls, ambulances, ins, travel, time));
  } else if (policy == "gen_forward") {
    return static_pointer_cast<Solver>(make_shared<GenForwardSolver>(
        env, calls, ambulances, ins, travel, time));
  } else if (policy == "gen_minmax") {
    return static_pointer_cast<Solver>(make_shared<GenMinMaxSolver>(
        env, calls, ambulances, ins, travel, time));
  } else if (policy == "cg") {
    return static_pointer_cast<Solver>(
        make_shared<CGSolver>(env, calls, ambulances, ins, travel, time));
  } else if (policy == "minmaxp" || policy == "minmaxp_a") {
    auto solver = static_pointer_cast<Solver>(
        make_shared<MinMaxPSolver>(env, calls, ambulances, ins, travel, time));
    solver->travel.set_forward(policy == "minmaxp");
    return solver;
  } else if (policy == "non_miopyc") {
    return static_pointer_cast<Solver>(make_shared<NonMiopycSolver>(
        env, calls, ambulances, ins, travel, time));
  } else if (policy == "preparedness") {
    return static_pointer_cast<Solver>(make_shared<PreparednessSolver>(
        env, calls, ambulances, ins, travel, g, time));
  } else if (policy == "prep2") {
    return static_pointer_cast<Solver>(
        make_shared<Prep2Solver>(env, calls, ambulances, ins, travel, g, time));
  } else if (policy == "ordered") {
    return static_pointer_cast<Solver>(make_shared<OrderedSolver>(
        env, calls, ambulances, ins, travel, g, time));
  } else if (policy == "district") {
    return static_pointer_cast<Solver>(make_shared<DistrictSolver>(
        env, calls, ambulances, ins, travel, g, data[g], lambda[g], time));
  } else if (policy == "coverage") {
    return static_pointer_cast<Solver>(make_shared<CoverageSolver>(
        env, calls, ambulances, ins, travel, g, time));
  } else if (policy == "forward_prep") {
    return static_pointer_cast<Solver>(make_shared<ForwardPrepSolver>(
        env, calls, ambulances, ins, travel, g, time));
  } else if (policy == "minmax_prep") {
    return static_pointer_cast<Solver>(make_shared<MinMaxPrepSolver>(
        env, calls, ambulances, ins, travel, g, time));
  } else if (policy == "priorities_prep") {
    return static_pointer_cast<Solver>(make_shared<PriorityPrepSolver>(
        env, calls, ambulances, ins, travel, g, time));
  } else if (policy == "enumerate") {
    return static_pointer_cast<Solver>(make_shared<EnumerateSolver>(
        env, calls, ambulances, ins, travel, g, time));
  } else if (policy == "markov_prep") {
    return static_pointer_cast<Solver>(
        make_shared<MCSolver>(env, calls, ambulances, ins, travel, g, time));
  } else if (policy == "queue_deficit") {
    return static_pointer_cast<Solver>(make_shared<QueueDeficitSolver>(
        env, calls, ambulances, ins, travel, g, time));
  }

  cout << "ERROR: Unknown Policy: " << policy << "\n";
  exit(1);
}

void PolicyTester::run() {
  // one_stage_old();
  fmt::print("One stage\n");
  one_stage();
  fmt::print("=============== END Simulate Policies ====================\n");
  // cin.get();
  // two_stage();
  // two_stage_tree();
  // fmt::print("Two stage\n");
  // two_stage_tree_new();
  // fmt::print("=============== END Two Stage ====================\n");
}

void PolicyTester::one_stage_old() {
  for (size_t s = 0; s < ins.calls.size(); ++s) {
    for (size_t i = 0; i < ins.calls[s].size(); ++i) {
      for (size_t a = 0; a < ambulances.size(); ++a) {
        ins.calls[s][i].ambulances.push_back(a);
      }
    }
  }
  int nb_scen = g_params.n_scenarios;

  vector<double> all_waiting_on_scene;
  vector<double> all_waiting_on_scene_penalized;
  vector<double> all_waiting_to_hospital;
  all_waiting_on_scene.reserve(ins.nb_scenarios * ins.calls[0].size());
  all_waiting_on_scene_penalized.reserve(ins.nb_scenarios *
                                         ins.calls[0].size());
  all_waiting_to_hospital.reserve(ins.nb_scenarios * ins.calls[0].size());

  for (auto &policy : policies) {
    int s = 0;
    all_waiting_on_scene.clear();
    all_waiting_on_scene_penalized.clear();
    all_waiting_to_hospital.clear();
    for (size_t s = 0; s < ins.calls.size(); ++s) {
      waiting_on_scene[s].clear();
      waiting_on_scene_penalized[s].clear();
    }
    run_times.clear();
    // fmt::print("Policy {}\n", policy);
    for (auto &scenario : ins.calls) {
      // fmt::print("Scenario {} ({} calls)\n", s, scenario.size());
      auto solver =
          get_solver(policy, scenario, ins.ambulances, travel, s / 100);
      solver->run();
      // solver->print_results();
      run_times.insert(run_times.end(), solver->run_times.begin(),
                       solver->run_times.end());
      for (size_t i = 0; i < scenario.size(); ++i) {
        waiting_on_scene[s].push_back(solver->waiting_on_scene[i]);
        waiting_on_scene_penalized[s].push_back(
            solver->waiting_on_scene_penalized[i]);
      }
      all_waiting_on_scene.insert(all_waiting_on_scene.end(),
                                  solver->waiting_on_scene.begin(),
                                  solver->waiting_on_scene.end());
      all_waiting_on_scene_penalized.insert(
          all_waiting_on_scene_penalized.end(),
          solver->waiting_on_scene_penalized.begin(),
          solver->waiting_on_scene_penalized.end());
      all_waiting_to_hospital.insert(all_waiting_to_hospital.end(),
                                     solver->waiting_to_hospital.begin(),
                                     solver->waiting_to_hospital.end());
      ++s;
    }

    Stats stats(all_waiting_on_scene, all_waiting_on_scene_penalized,
                all_waiting_to_hospital);
    ofstream policy_results(
        fmt::format("results/one_stage_{}_a{}_s{}.txt", policy,
                    ambulances.size(), ins.calls.size()),
        ios::out);
    // policy_results << g_params.n_scenarios << "\n";
    for (int k = 0; k < g_params.n_scenarios; ++k) {
      policy_results << waiting_on_scene[k].size() << "\n";
      for (size_t i = 0; i < waiting_on_scene[k].size(); ++i) {
        auto &call = ins.calls[k][i];
        double ws = waiting_on_scene[k][i];
        double wsp = waiting_on_scene_penalized[k][i];
        // int which_amb = which_ambulance[k][i];
        // int call_end = calls_end[k][i];
        // policy_results << call.id << " " << call.time << " " << call.priority
        // << " "; policy_results << call.id << " " << call.time << " " <<
        // call.priority << " "; policy_results << fmt::format("{:.6f} {:.6f} ",
        // call.location.first, call.location.second); policy_results << ws << "
        // " << wsp << " " << call_end << " " << which_amb << "\n";
        policy_results << ws << " " << wsp << "\n";
      }
    }
    policy_results.close();

    fmt::print(
        "Policy {} Mean time = {:.1f} Mean pen = {:.1f}, Max pen = {:.1f}\n",
        policy, stats.mean_waiting_on_scene,
        stats.mean_waiting_on_scene_penalized,
        stats.max_waiting_on_scene_penalized);
    double min_run_time = *min_element(run_times.begin(), run_times.end());
    double avg_run_time =
        accumulate(run_times.begin(), run_times.end(), 0.0) / run_times.size();
    double max_run_time = *max_element(run_times.begin(), run_times.end());
    fmt::print("Run times: {}\t{}\t{}\n", min_run_time, avg_run_time,
               max_run_time);
    // std::cin.get();
  }
}

void PolicyTester::one_stage() {
  for (size_t s = 0; s < ins.calls.size(); ++s) {
    for (size_t i = 0; i < ins.calls[s].size(); ++i) {
      for (size_t a = 0; a < ambulances.size(); ++a) {
        ins.calls[s][i].ambulances.push_back(a);
      }
    }
  }

  bool debug = false;

  size_t num_calls_t = 100;
  int nb_scen_simu = 100;
  double t0 = 36 * 1800;
  int T = 4;  // 4
  int nb_realizations = 5;
  vector<vector<vector<Call>>> my_scenarios(
      T, vector<vector<Call>>(nb_realizations, vector<Call>()));

  for (int scen = 0; scen < nb_realizations; ++scen) {
    size_t j = 0;
    for (int t = 0; t < T; ++t) {
      while (j < ins.calls[scen].size() &&
             ins.calls[scen][j].time <= t0 + (t + 1) * 1800) {  // 1800
        my_scenarios[t][scen].push_back(ins.calls[scen][j]);
        ++j;
      }
    }
  }

  // fmt::print("Joined scenarios\n");

  // for(int t = 0; t < T; ++t){
  // 	for(int scen = 0; scen < nb_realizations; ++scen){
  // 		my_scenarios[t][scen] =
  // vector<Call>(my_scenarios[t][scen].begin(),
  // 			my_scenarios[t][scen].begin()+2);
  // 	}
  // }

  default_random_engine gen;
  uniform_int_distribution<int> rand_scen(0, nb_realizations - 1);
  block_indexes = vector<vector<int>>(g_params.n_scenarios, vector<int>(T, 0));
  for (int simu = 0; simu < g_params.n_scenarios; ++simu) {
    for (int t = 0; t < T; ++t) {
      block_indexes[simu][t] = rand_scen(gen);
    }
    if (false) {
      fmt::print("s = {}. block_indexes = {}\n", simu, block_indexes[simu]);
    }
  }
  // FIXME
  // block_indexes[0] = {1,2,4,4};

  vector<double> all_waiting_on_scene;
  vector<double> all_waiting_on_scene_penalized;
  vector<double> all_waiting_to_hospital;
  all_waiting_on_scene.reserve(ins.nb_scenarios * ins.calls[0].size());
  all_waiting_on_scene_penalized.reserve(ins.nb_scenarios *
                                         ins.calls[0].size());
  all_waiting_to_hospital.reserve(ins.nb_scenarios * ins.calls[0].size());
  vector<string> final_log;
  string run_times_file_name =
      fmt::format("results/tables/run_times_one_stage_{}_a{}.txt",
                  g_params.amb_setup, g_params.n_ambulances);
  ofstream run_times_file(run_times_file_name, std::ios::out);
  run_times_file << "Heuristic\tmin (us)\tmean (us)\tq0.9 (us)\tmax (us)\n";
  for (auto &policy : policies) {
    if (policy == "cg" || policy == "enumerate") {
      continue;
    }
    // int s = 0;
    all_waiting_on_scene.clear();
    all_waiting_on_scene_penalized.clear();
    all_waiting_to_hospital.clear();
    run_times.clear();
    fmt::print("Policy {}\n", policy);
    // for(auto& scenario: ins.calls){
    for (int s = 0; s < g_params.n_scenarios; ++s) {
      vector<Call> this_scenario;
      vector<int> this_scenario_indexes;
      for (int t = 0; t < T; ++t) {
        uniform_int_distribution<int> rand_scen(0, nb_realizations - 1);
        int u = block_indexes[s][t];
        // int u = block_indexes[0][t];
        this_scenario_indexes.push_back(u);
        auto my_end = (my_scenarios[t][u].size() > num_calls_t)
                          ? my_scenarios[t][u].begin() + num_calls_t
                          : my_scenarios[t][u].end();
        this_scenario.insert(this_scenario.end(), my_scenarios[t][u].begin(),
                             my_end);
      }
      for (size_t i = 0; i < this_scenario.size(); ++i) {
        this_scenario[i].id = i;
      }

      waiting_on_scene[s] = vector<double>(this_scenario.size(), GRB_INFINITY);
      waiting_on_scene_penalized[s] =
          vector<double>(this_scenario.size(), GRB_INFINITY);
      waiting_to_hospital[s] =
          vector<double>(this_scenario.size(), GRB_INFINITY);
      calls_end[s] = vector<double>(this_scenario.size(), GRB_INFINITY);
      which_ambulance[s] = vector<int>(this_scenario.size(), -1);
      // fmt::print("Scenario {} ({} calls)\n", s, scenario.size());
      auto solver =
          get_solver(policy, this_scenario, ins.ambulances, travel, 4, 0);
      // solver->set_debug_mode(true);
      solver->run();
      if (debug) {
        solver->print_results();
        fmt::print("/\\ Results s =  {}\n", s);
        cin.get();
      }
      if (policy == "cg" || policy == "enumerate") {
        continue;
      }
      // solver->print_results();
      run_times.insert(run_times.end(), solver->run_times.begin(),
                       solver->run_times.end());
      for (size_t i = 0; i < this_scenario.size(); ++i) {
        waiting_on_scene[s][i] = solver->waiting_on_scene[i];
        waiting_on_scene_penalized[s][i] =
            solver->waiting_on_scene_penalized[i];
        waiting_to_hospital[s][i] = solver->waiting_to_hospital[i];
        calls_end[s][i] = solver->calls_end[i];
        which_ambulance[s][i] = solver->which_ambulance[i];
      }
      all_waiting_on_scene.insert(all_waiting_on_scene.end(),
                                  solver->waiting_on_scene.begin(),
                                  solver->waiting_on_scene.end());
      all_waiting_on_scene_penalized.insert(
          all_waiting_on_scene_penalized.end(),
          solver->waiting_on_scene_penalized.begin(),
          solver->waiting_on_scene_penalized.end());
      all_waiting_to_hospital.insert(all_waiting_to_hospital.end(),
                                     solver->waiting_to_hospital.begin(),
                                     solver->waiting_to_hospital.end());
    }
    string type_return = (g_params.best_base) ? "_bb" : "";
    ofstream policy_results(
        fmt::format("results/one_stage/{}_{}_{}_{}{}.txt", g_params.amb_setup,
                    policy, ambulances.size(), ins.calls.size(), type_return),
        ios::out);
    // policy_results << g_params.n_scenarios << "\n";
    for (int k = 0; k < g_params.n_scenarios; ++k) {
      policy_results << waiting_on_scene[k].size() << "\n";
      for (size_t i = 0; i < waiting_on_scene[k].size(); ++i) {
        auto &call = ins.calls[k][i];
        double ws = waiting_on_scene[k][i];
        double wsp = waiting_on_scene_penalized[k][i];
        int which_amb = which_ambulance[k][i];
        double call_end = calls_end[k][i];
        // policy_results << call.id << " " << call.time << " " << call.priority
        // << " "; policy_results << call.id << " " << call.time << " " <<
        // call.priority << " "; policy_results << fmt::format("{:.6f} {:.6f} ",
        // call.location.first, call.location.second); policy_results << ws << "
        // " << wsp << " " << call_end << " " << which_amb << "\n";
        policy_results << which_amb << " " << ws << " " << wsp << " "
                       << call_end << "\n";
      }
    }
    policy_results.close();

    Stats stats(all_waiting_on_scene, all_waiting_on_scene_penalized,
                all_waiting_to_hospital);
    final_log.push_back(fmt::format(
        "Policy {} Mean actual = {:.0f} Mean pen = {:.0f}, Max pen = {:.1f}\n",
        policy, stats.mean_waiting_on_scene,
        stats.mean_waiting_on_scene_penalized,
        stats.max_waiting_on_scene_penalized));
    // fmt::print(final_log.back());

    sort(run_times.begin(), run_times.end());
    double min_run_time = run_times[0];
    double avg_run_time =
        accumulate(run_times.begin(), run_times.end(), 0.0) / run_times.size();
    double max_run_time = run_times.back();
    size_t index = static_cast<size_t>(0.9 * (run_times.size() - 1));
    double q09_run_time = (run_times.size() > 0) ? run_times[index] : 0.0;
    run_times_file << fmt::format("{}\t{}\t{}\t{}\t{}\n", policy, min_run_time,
                                  avg_run_time, q09_run_time, max_run_time);
    final_log.push_back(fmt::format("Run times: {}\t{}\t{}\n", min_run_time,
                                    avg_run_time, max_run_time));
    // fmt::print(final_log.back());
    // std::cin.get();
  }

  run_times_file.close();
  for (auto &line : final_log) {
    fmt::print(line);
  }
}

void PolicyTester::two_stage() {
  vector<double> all_waiting_on_scene;
  vector<double> all_waiting_on_scene_penalized;
  vector<double> all_waiting_to_hospital;

  bool debug = false;

  for (size_t s = 0; s < ins.calls.size(); ++s) {
    for (size_t i = 0; i < ins.calls[s].size(); ++i) {
      for (size_t a = 0; a < ambulances.size(); ++a) {
        ins.calls[s][i].ambulances.push_back(a);
      }
    }
  }

  for (auto &policy : policies) {  // each policy
    int sc = 0;
    fmt::print("Policy: {}\n", policy);
    all_waiting_on_scene.clear();
    all_waiting_on_scene_penalized.clear();
    all_waiting_to_hospital.clear();
    run_times.clear();
    for (int i = 0; i < ins.nb_scenarios; ++i) {
      which_ambulance[i] = std::vector<int>(ins.nb_calls[i], -1);
      waiting_on_scene[i] = std::vector<double>(ins.nb_calls[i], GRB_INFINITY);
      waiting_on_scene_penalized[i] =
          std::vector<double>(ins.nb_calls[i], GRB_INFINITY);
      waiting_to_hospital[i] =
          std::vector<double>(ins.nb_calls[i], GRB_INFINITY);
    }

    // each scenario
    for (auto &scenario : ins.calls) {
      auto nearest_base = get_nearest_base(scenario);
      ambulances = ins.ambulances;
      size_t calls_attended = 0;
      int event_call = 1;
      int index_call = 0;
      time = scenario[0].time;
      int amb_finish = -1;

      std::vector<int> queue;
      queue.reserve(scenario.size());
      int iter_count = 0;
      while (calls_attended < scenario.size()) {
        if (debug)
          fmt::print("calls_attended {} index_call {}\n", calls_attended,
                     index_call);
        vector<int> index_scenarios(ins.calls.size(), 0);
        // Find future calls in each scenario.
        for (size_t j = 0; j < ins.calls.size(); ++j) {
          while (static_cast<size_t>(index_scenarios[j]) <
                     ins.calls[j].size() &&
                 ins.calls[j][index_scenarios[j]].time <= time) {
            ++index_scenarios[j];
          }
        }

        if (event_call == 1) {
          queue.push_back(index_call);
          std::vector<int> queue_aux;
          queue_aux.reserve(queue.size());
          bool increased_ca = false;
          for (size_t i = 0; i < queue.size(); ++i) {
            auto t0 = std::chrono::high_resolution_clock::now();
            auto &call = scenario[queue[i]];
            double min_time = GRB_INFINITY;
            int index_amb = -1;
            double min_waiting_on_scene = GRB_INFINITY;
            double min_waiting_to_hospital = GRB_INFINITY;
            if (debug) {
              cout << "Call " << call << "\n";
            }

            for (size_t k = 0; k < ambulances.size(); ++k) {
              auto &ambulance = ambulances[k];
              if (ambulance.type <= call.priority &&
                  ambulance.arrival_time_at_f_last_trip <= time) {
                auto result = get_waiting_time(
                    ambulance.id, call, policy, index_scenarios,
                    nearest_base[queue[i]], queue, sc, time);
                if (debug) {
                  cout << "\t" << ambulance << " " << result.mean_total_time;
                  cout << "\n";
                }
                if (result.mean_total_time < min_time) {
                  index_amb = k;
                  min_time = result.mean_total_time;
                }
              }
            }

            auto result =
                get_waiting_time(-1, call, policy, index_scenarios,
                                 nearest_base[queue[i]], queue, sc, time);

            if (debug) cout << "\tqueue " << result.mean_total_time << "\n";
            if (result.mean_total_time < min_time) {
              index_amb = -1;
              min_time = result.mean_total_time;
            }

            if (index_amb >= 0) {
              min_time = travel.get_response_time(ambulances[index_amb], call,
                                                  time, true);
              double waiting_on_scene_i = time + min_time - call.time;
              double waiting_to_hospital_i = ambulances[index_amb].answer_call(
                  call, travel, ins, time, min_time, nearest_base[queue[i]]);
              which_ambulance[sc][queue[i]] = index_amb;
              waiting_on_scene[sc][queue[i]] = waiting_on_scene_i;
              waiting_on_scene_penalized[sc][queue[i]] =
                  penalized_response_time(waiting_on_scene_i,
                                          ambulances[index_amb].type,
                                          call.priority);
              waiting_to_hospital[sc][queue[i]] = waiting_to_hospital_i;
              calls_attended++;
              increased_ca = true;
            } else {
              if (debug) cout << "queue\n";
              for (size_t a = 0; a < ambulances.size(); ++a) {
                if (ambulances[a].arrival_time_at_f_last_trip <= time) {
                  call.ambulances.erase(
                      remove(call.ambulances.begin(), call.ambulances.end(), a),
                      call.ambulances.end());
                }
              }
              queue_aux.push_back(queue[i]);
            }
            auto dt = std::chrono::high_resolution_clock::now();
            run_times.push_back(
                std::chrono::duration_cast<chrono::nanoseconds>(dt - t0)
                    .count() /
                pow(10, 9));
          }
          queue.clear();
          for (auto i : queue_aux) {
            queue.push_back(i);
          }

          if (increased_ca) {
            iter_count = 0;
          } else {
            iter_count++;
          }
        } else {
          // Event_call == 0, get ambulance that just finished.
          if (amb_finish == -1) {
            fmt::print("ERROR! event_call = {}, but no ambulance is returning",
                       event_call);
            cin.get();
          }
          auto &amb = ambulances[amb_finish];
          if (debug) {
            cout << "amb_finish = " << amb_finish << "\n";
            cout << "return " << amb << "\n";
          }
          std::vector<int> queue_aux;
          queue_aux.reserve(queue.size());

          int ind_call = -1;
          int ind_base = -1;
          double min_time = GRB_INFINITY;
          double min_waiting_on_scene = GRB_INFINITY;
          double min_waiting_to_hospital = GRB_INFINITY;
          for (size_t i = 0; i < queue.size(); ++i) {
            auto &call = scenario[queue[i]];
            if (can_answer(amb, call)) {
              auto result =
                  get_waiting_time(amb_finish, call, policy, index_scenarios,
                                   nearest_base[queue[i]], queue, sc, time);
              if (debug) {
                cout << "\tCall " << call << " " << result.mean_total_time
                     << "\n";
              }
              if (result.mean_total_time < min_time) {
                ind_call = queue[i];
                min_time = result.mean_total_time;
              }
            }
          }

          for (size_t b = 0; b < ins.bases.size(); ++b) {
            auto result = get_waiting_time_return(
                amb_finish, policy, index_scenarios, time, queue, sc, b);

            if (debug) {
              cout << "\tbase " << b << " " << result.mean_total_time << "\n";
            }

            if (result.mean_total_time < min_time) {
              ind_call = -1;
              ind_base = b;
              min_time = result.mean_total_time;
            }
          }

          if (ind_call >= 0) {
            auto &call = scenario[ind_call];
            double min_time = travel.get_response_time(amb, call, time, true);
            double waiting_on_scene_i = time + min_time - call.time;
            calls_attended++;
            iter_count = 0;
            double waiting_to_hospital_i = amb.answer_call(
                call, travel, ins, time, min_time, nearest_base[ind_call]);
            which_ambulance[sc][ind_call] = amb_finish;
            waiting_on_scene[sc][ind_call] = waiting_on_scene_i;
            waiting_on_scene_penalized[sc][ind_call] = penalized_response_time(
                waiting_on_scene_i, amb.type, call.priority);
            waiting_to_hospital[sc][ind_call] = waiting_to_hospital_i;
            queue.erase(remove(queue.begin(), queue.end(), ind_call),
                        queue.end());
          } else if (ind_base >= 0) {
            auto base = ins.bases[ind_base];
            amb.base_location = base;
            amb.arrival_time_at_b_last_trip =
                time +
                travel.travel_time(amb.free_location, amb.base_location, amb);

            fmt::print("Chosen Base {}, arrival = {}\n", ind_base,
                       amb.arrival_time_at_b_last_trip);

            for (auto i : queue) {
              scenario[i].ambulances.erase(
                  remove(scenario[i].ambulances.begin(),
                         scenario[i].ambulances.end(), amb_finish),
                  scenario[i].ambulances.end());
            }
            // iter_count++;
          }
        }
        amb_finish = set_next_event(scenario, event_call, index_call, queue);
        // if(debug)
        // 	cin.get();
      }
    }

    int k = 0;
    for (int s = 0; s < ins.nb_scenarios; ++s) {
      for (size_t i = 0; i < ins.calls[s].size(); ++i) {
        // fmt::print("{} {} {}\n",s,i,waiting_on_scene_penalized[s][i]);
        all_waiting_on_scene.push_back(waiting_on_scene[s][i]);
        all_waiting_on_scene_penalized.push_back(
            waiting_on_scene_penalized[s][i]);
        all_waiting_to_hospital.push_back(waiting_to_hospital[s][i]);
      }
    }

    Stats stats(all_waiting_on_scene, all_waiting_on_scene_penalized,
                all_waiting_to_hospital);
    fmt::print("Policy {} Mean pen = {:.1f}, Max pen = {:.1f}\n", policy,
               stats.mean_waiting_on_scene_penalized,
               stats.max_waiting_on_scene_penalized);
    double min_run_time = *min_element(run_times.begin(), run_times.end());
    double avg_run_time =
        accumulate(run_times.begin(), run_times.end(), 0.0) / run_times.size();
    double max_run_time = *max_element(run_times.begin(), run_times.end());
    fmt::print("Run times: {}\t{}\t{}\n", min_run_time, avg_run_time,
               max_run_time);

    std::string base_filename =
        g_params.instance.substr(g_params.instance.find_last_of("/\\") + 1);
    std::string::size_type const p(base_filename.find_last_of('.'));
    std::string file_without_extension = base_filename.substr(0, p);
    std::ofstream results(
        fmt::format("results/{}_{}_{}.txt", file_without_extension, policy,
                    ins.nb_scenarios),
        std::ios::out);
    for (int i = 0; i < ins.nb_scenarios; ++i) {
      results << ins.nb_calls[i] << " ";
      for (int j = 0; j < ins.nb_calls[i]; ++j) {
        results << waiting_on_scene[i][j] << " ";
      }
      results << "\n";
    }
    results.close();
    // std::cin.get();
  }
}

void PolicyTester::two_stage_tree_new() {
  vector<double> all_waiting_on_scene;
  vector<double> all_waiting_on_scene_penalized;
  vector<double> all_waiting_to_hospital;

  bool debug = false;

  for (size_t s = 0; s < ins.calls.size(); ++s) {
    for (size_t i = 0; i < ins.calls[s].size(); ++i) {
      for (size_t a = 0; a < ambulances.size(); ++a) {
        ins.calls[s][i].ambulances.push_back(a);
      }
    }
  }

  size_t num_calls_t = 100;  // calls by time slot
  size_t num_decisions =
      ambulances.size();  // number of closest ambulances available
  double t0 = ins.time_horizon[0] * ins.slot_duration;  // initial time of tests
  int T = ins.time_horizon.size();  // number of time horizons
  int nb_realizations = 5;
  vector<vector<vector<Call>>> my_scenarios(
      T, vector<vector<Call>>(nb_realizations, vector<Call>()));

  for (int scen = 0; scen < nb_realizations; ++scen) {
    size_t j = 0;
    for (int t = 0; t < T; ++t) {
      while (j < ins.calls[scen].size() &&
             ins.calls[scen][j].time <=
                 t0 + (t + 1) * ins.slot_duration) {  // 1800
        my_scenarios[t][scen].push_back(ins.calls[scen][j]);
        ++j;
      }
      if (num_calls_t < my_scenarios[t][scen].size()) {
        my_scenarios[t][scen].erase(my_scenarios[t][scen].begin() + num_calls_t,
                                    my_scenarios[t][scen].end());
      }
    }
  }

  default_random_engine gen;
  uniform_int_distribution<int> rand_scen(0, nb_realizations - 1);
  block_indexes = vector<vector<int>>(g_params.n_scenarios, vector<int>(T, 0));
  for (int simu = 0; simu < g_params.n_scenarios; ++simu) {
    for (int t = 0; t < T; ++t) {
      block_indexes[simu][t] = rand_scen(gen);
    }
    if (debug) {
      fmt::print("s = {}. block_indexes = {}\n", simu, block_indexes[simu]);
    }
  }

  // FIXME
  // block_indexes[0] = {1,2,4,4};

  uniform_real_distribution<double> u(0, 1);
  vector<string> out_lines;

  string run_times_file_name =
      fmt::format("results/tables/run_times_two_stage_{}_{}.txt",
                  g_params.amb_setup, g_params.n_ambulances);
  ofstream run_times_file(run_times_file_name, std::ios::out);
  run_times_file << "Heuristic\tmin (us)\tmean (us)\tq0.9 (us)\tmax (us)\n";

  for (auto &policy : policies) {
    g_params.h_use_fixed_bases = (policy == "ordered" || policy == "district");
    bool free_at_base = (policy == "ordered" || policy == "district");
    run_times.clear();
    all_waiting_on_scene.clear();
    all_waiting_on_scene_penalized.clear();
    all_waiting_to_hospital.clear();
    vector<double> run_times_calls;
    vector<double> run_times_bases;

    for (int simu = 0; simu < g_params.n_scenarios; ++simu) {
      // fmt::print("simu = {}, block = {}\n", simu, block_indexes[simu]);
      if (simu > 0 && simu % 25 == 0) {
        fmt::print("{} {} scenarios completed\n", policy, simu);
        // cin.get();
      }
      default_random_engine u_gen;
      vector<Call> this_scenario;
      vector<int> this_scenario_indexes;
      for (int t = 0; t < 4; ++t) {  // 4
        int u = block_indexes[simu][t];
        // int u = block_indexes[0][t];
        this_scenario_indexes.push_back(u);
        auto my_end = (my_scenarios[t][u].size() > num_calls_t)
                          ? my_scenarios[t][u].begin() + num_calls_t
                          : my_scenarios[t][u].end();
        this_scenario.insert(this_scenario.end(), my_scenarios[t][u].begin(),
                             my_end);
      }
      // fmt::print("This scenario:\n");
      for (size_t i = 0; i < this_scenario.size(); ++i) {
        this_scenario[i].id = i;
        this_scenario[i].ambulances.clear();
        for (size_t a = 0; a < ambulances.size(); ++a) {
          this_scenario[i].ambulances.push_back(a);
        }
        // cout << this_scenario[i] << "\n";
        // fmt::print("{},{},{},{}\n", this_scenario[i].id,
        // this_scenario[i].priority, this_scenario[i].location.first,
        // this_scenario[i].location.second);
      }
      // cin.get();

      waiting_on_scene[simu] =
          vector<double>(this_scenario.size(), GRB_INFINITY);
      waiting_on_scene_penalized[simu] =
          vector<double>(this_scenario.size(), GRB_INFINITY);
      waiting_to_hospital[simu] =
          vector<double>(this_scenario.size(), GRB_INFINITY);
      calls_end[simu] = vector<double>(this_scenario.size(), GRB_INFINITY);
      which_ambulance[simu] = vector<int>(this_scenario.size(), -1);

      auto nearest_base = get_nearest_base(this_scenario);
      ambulances = ins.ambulances;
      if (policy == "district") {
        for (size_t a = 0; a < ambulances.size(); ++a) {
          ambulances[a].base_location =
              ins.bases[data[4].ambulance_district[a]];
        }
      }
      size_t calls_attended = 0;
      int event_call = 1;
      int index_call = 0;
      time = this_scenario[0].time;
      int amb_finish = -1;

      std::vector<int> queue;
      queue.reserve(this_scenario.size());
      double current_sum = 0;
      int current_solved = 0;
      while (calls_attended < this_scenario.size()) {
        vector<vector<Call>> future_scenarios = get_future_scenarios(
            time, this_scenario, nb_realizations, T, my_scenarios);

        // if ambulance arrived at base it's now free to answer all calls.
        //  for(auto& amb: ambulances){
        //  	if(amb.arrival_time_at_b_last_trip <= time){
        //  		for(size_t i = 0; i < this_scenario.size(); ++i){
        //  			this_scenario[i].ambulances.push_back(amb.id);
        //  		}
        //  	}
        //  }
        //  travel.set_delay_param(u(u_gen));
        //   travel.set_delay_param(-1);
        if (debug) {
          fmt::print(
              "Time = {}, index_call = {}, ev = {}, calls_attended = {} "
              "/ {}, queue = {}, policy = {}\n",
              time, index_call, event_call, calls_attended,
              this_scenario.size(), queue, policy);
          vector<double> queue_times;
          for (auto i : queue) {
            queue_times.push_back(this_scenario[i].time);
          }
          fmt::print("Times: {}\n", queue_times);
          current_sum = 0;
          current_solved = 0;
          int end_index = (event_call == 0) ? index_call : index_call - 1;
          for (int i = 0; i <= end_index; ++i) {
            auto &call = this_scenario[i];
            // cout << call << "\t" << waiting_on_scene[simu][i] << "\t";
            // cout << waiting_on_scene_penalized[simu][i] << "\t";
            // cout << which_ambulance[simu][i] << "\n";
            if (waiting_on_scene_penalized[simu][i] < GRB_INFINITY) {
              current_sum += waiting_on_scene_penalized[simu][i];
              current_solved += 1;
            }
          }
          double current_mean =
              (current_solved == 0) ? 0 : current_sum / current_solved;
          fmt::print(
              "Current sum = {}, Current mean = {:0f}, Current solved = {}\n",
              current_sum, current_mean, current_solved);
        }

        if (event_call == 1) {
          auto t0 = std::chrono::high_resolution_clock::now();
          auto &call = this_scenario[index_call];
          double min_time = GRB_INFINITY;
          int index_amb = -1;
          double min_waiting_on_scene = GRB_INFINITY;
          double min_waiting_to_hospital = GRB_INFINITY;
          if (debug) {
            int min_b = -1;
            double min_d = GRB_INFINITY;
            for (size_t b = 0; b < ins.bases.size(); ++b) {
              double d = travel.lat_long_distance(call.location, ins.bases[b]);
              if (d < min_d) {
                min_d = d;
                min_b = b;
              }
            }
            cout << "Call " << call << " closest to base " << min_b << "\n";
          }
          for (size_t k = 0; k < ambulances.size(); ++k) {
            auto &ambulance = ambulances[k];
            if (can_answer(ambulance, call) &&
                ambulance.arrival_time_at_f_last_trip <= time) {
              if (debug) {
                string location_info;
                int min_b = -1;
                double min_d = GRB_INFINITY;
                for (size_t b = 0; b < ins.bases.size(); ++b) {
                  double d = travel.lat_long_distance(ambulance.base_location,
                                                      ins.bases[b]);
                  if (d < min_d) {
                    min_d = d;
                    min_b = b;
                  }
                }

                if (ambulance.arrival_time_at_b_last_trip <= time) {
                  location_info =
                      fmt::format(", currently at base {}\n", min_b);
                } else {
                  location_info = fmt::format(
                      ", currently towards base {}, arrival_time = {:.1f}\n",
                      min_b, ambulance.arrival_time_at_b_last_trip);
                }
                cout << "Running Amb " << ambulance.id << location_info;
              }
              auto result = get_waiting_time_tree(
                  ambulance.id, call.id, policy, future_scenarios,
                  this_scenario, time, queue, nearest_base[index_call]);
              if (debug) {
                double decision_total = current_sum + result.mean_total_time;
                cout << "\t" << ambulance << " => " << result.mean_total_time
                     << " | total = " << decision_total;
                cout << "\n";
              }
              if (result.mean_total_time < min_time) {
                index_amb = k;
                min_time = result.mean_total_time;
              }
            }
          }
          if (debug) {
            cout << "Testing queue for " << call << "\n";
          }
          auto result = get_waiting_time_tree(
              -1, call.id, policy, future_scenarios, this_scenario, time, queue,
              nearest_base[index_call]);
          if (debug) {
            double decision_total = current_sum + result.mean_total_time;
            cout << "\tQueue =>" << result.mean_total_time
                 << "| total = " << decision_total << "\n";
          }
          if (result.mean_total_time < min_time) {
            index_amb = -1;
            min_time = result.mean_total_time;
          }
          if (index_amb >= 0) {
            auto &best_amb = ambulances[index_amb];
            min_time = travel.get_response_time(best_amb, call, time, true);
            double waiting_on_scene_i = time + min_time - call.time;
            double waiting_to_hospital_i = ambulances[index_amb].answer_call(
                call, travel, ins, time, min_time, nearest_base[index_call]);
            if (free_at_base) {
              best_amb.free_location = best_amb.base_location;
              best_amb.arrival_time_at_f_last_trip =
                  best_amb.arrival_time_at_b_last_trip;
            }
            which_ambulance[simu][index_call] = index_amb;
            waiting_on_scene[simu][index_call] = waiting_on_scene_i;
            waiting_on_scene_penalized[simu][index_call] =
                penalized_response_time(waiting_on_scene_i,
                                        ambulances[index_amb].type,
                                        call.priority);
            if (debug) {
              Location current_location;
              string type_current_location = "-";
              if (best_amb.arrival_time_at_b_last_trip <= time) {
                type_current_location = "b";
                current_location = best_amb.base_location;
              } else if (best_amb.arrival_time_at_f_last_trip < time) {
                type_current_location = "r";
                current_location = travel.ambulance_position(best_amb, time);
              } else if (abs(best_amb.arrival_time_at_f_last_trip - time) <
                         g_params.EPS) {
                type_current_location = "h";
                current_location = best_amb.free_location;
              }
              cout << "Chosen " << best_amb
                   << fmt::format(
                          " | {:.6f} {:.6f}({})", current_location.first,
                          current_location.second, type_current_location)
                   << " pt ";
              cout << waiting_on_scene_penalized[simu][index_call]
                   << fmt::format(" call = ({}, {})", index_call, call.time)
                   << "\n";
            }
            waiting_to_hospital[simu][index_call] = waiting_to_hospital_i;
            calls_attended++;
            auto dt = std::chrono::high_resolution_clock::now();
            run_times_calls.push_back(
                std::chrono::duration_cast<chrono::nanoseconds>(dt - t0)
                    .count() /
                pow(10, 9));
          } else {
            bool call_already_in_queue =
                find(queue.begin(), queue.end(), index_call) != queue.end();
            if (!call_already_in_queue) {
              queue.push_back(index_call);
            }
            if (debug) {
              fmt::print("Chose queue\n");
              if (call_already_in_queue) {
                fmt::print(
                    "WEIRD: Call in queue was reevaluated and went to "
                    "queue again.\n");
                cin.get();
              }
            }
          }
        } else {
          if (amb_finish == -1) {
            fmt::print("WEIRD!!! event_call = 0 but no return\n");
            fmt::print("Runtimes.size() = {}\n", run_times.size());
            fmt::print("policy = {}, s = {}, time = {}\n", policy, simu, time);
            double min_run_time =
                *min_element(run_times.begin(), run_times.end());
            double avg_run_time =
                accumulate(run_times.begin(), run_times.end(), 0.0) /
                run_times.size();
            double max_run_time =
                *max_element(run_times.begin(), run_times.end());
            fmt::print("Run times: {}\t{}\t{}\n", min_run_time, avg_run_time,
                       max_run_time);
            cin.get();
          }
          auto &amb = ambulances[amb_finish];
          std::vector<int> queue_aux;
          queue_aux.reserve(queue.size());
          int ind_call = -1;
          int ind_base = -1;
          double min_time = GRB_INFINITY;
          double min_waiting_on_scene = GRB_INFINITY;
          double min_waiting_to_hospital = GRB_INFINITY;
          double min_pen_time = GRB_INFINITY;
          for (size_t i = 0; i < queue.size(); ++i) {
            auto t0 = std::chrono::high_resolution_clock::now();
            auto &call = this_scenario[queue[i]];
            if (debug && !can_answer(amb, call)) {
              cout << "Call " << call << " | Can't be answered by " << amb.id
                   << "\n";
            }
            if (can_answer(amb, call)) {
              if (debug) {
                cout << "Testing " << amb << " for " << call << "\n";
              }
              auto result = get_waiting_time_tree(
                  amb_finish, call.id, policy, future_scenarios, this_scenario,
                  time, queue, nearest_base[queue[i]]);

              if (debug) {
                double decision_total = current_sum + result.mean_total_time;
                cout << "\tAmb " << amb.id << " Call " << call << " => "
                     << result.mean_total_time;
                cout << "| total = " << decision_total << "\n";
              }
              if (result.mean_total_time < min_time) {
                ind_call = queue[i];
                min_time = result.mean_total_time;
                for (auto &amb_aux : ambulances) {
                  if (can_answer(amb_aux, call) && amb_aux.id != amb.id &&
                      amb_aux.arrival_time_at_f_last_trip > time &&
                      amb_aux.arrival_time_at_f_last_trip - time <
                          min_pen_time) {
                    min_pen_time = amb_aux.arrival_time_at_f_last_trip - time;
                  }
                }
              }
            }
            auto dt = std::chrono::high_resolution_clock::now();
            run_times_calls.push_back(
                std::chrono::duration_cast<chrono::nanoseconds>(dt - t0)
                    .count() /
                pow(10, 9));
          }

          if (!free_at_base) {
            min_pen_time = (min_pen_time == GRB_INFINITY) ? 0 : min_pen_time;
            auto t0 = std::chrono::high_resolution_clock::now();
            for (size_t b = 0; b < ins.bases.size(); ++b) {
              auto result = get_waiting_time_tree_return(
                  amb_finish, b, policy, future_scenarios, this_scenario, time,
                  queue, min_pen_time);
              if (debug) {
                double decision_total = current_sum + result.mean_total_time;
                cout << "\tAmb " << amb.id << " Base " << b << " => "
                     << result.mean_total_time
                     << " | total = " << decision_total << "\n";
              }
              if (result.mean_total_time < min_time) {
                min_time = result.mean_total_time;
                ind_call = -1;
                ind_base = b;
              }
            }
            auto dt = std::chrono::high_resolution_clock::now();
            run_times_bases.push_back(
                std::chrono::duration_cast<chrono::nanoseconds>(dt - t0)
                    .count() /
                pow(10, 9));
          }
          if (ind_call >= 0) {
            auto &call = this_scenario[ind_call];
            double min_time = travel.get_response_time(amb, call, time, true);
            double waiting_on_scene_i = time + min_time - call.time;
            calls_attended++;
            // iter_count = 0;
            double waiting_to_hospital_i = amb.answer_call(
                call, travel, ins, time, min_time, nearest_base[ind_call]);
            if (free_at_base) {
              amb.free_location = amb.base_location;
              amb.arrival_time_at_f_last_trip = amb.arrival_time_at_b_last_trip;
            }
            which_ambulance[simu][ind_call] = amb_finish;
            waiting_on_scene[simu][ind_call] = waiting_on_scene_i;
            waiting_on_scene_penalized[simu][ind_call] =
                penalized_response_time(waiting_on_scene_i, amb.type,
                                        call.priority);
            if (debug) {
              cout << "Chosen call " << call << " ";
              cout << waiting_on_scene_penalized[simu][ind_call] << "\n";
            }
            waiting_to_hospital[simu][ind_call] = waiting_to_hospital_i;
            queue.erase(remove(queue.begin(), queue.end(), ind_call),
                        queue.end());
          } else if (ind_base >= 0) {
            // cout << "Chose Return\n";
            auto base = ins.bases[ind_base];
            amb.base_location = base;
            if (!free_at_base) {
              amb.arrival_time_at_b_last_trip =
                  time +
                  travel.travel_time(amb.free_location, amb.base_location, amb);
            } else {
              amb.arrival_time_at_b_last_trip = time;
            }

            if (debug) {
              fmt::print("Chosen Base {}\n", ind_base);
            }
            // If ambulance must return to base, it can't answer any calls in
            // queue. if(debug){ 	fmt::print("Amb {} forbidden to answer
            // to: ", amb_finish);
            // }
            // for(auto i: queue){
            // 	if(debug){
            // 		fmt::print("{} ", this_scenario[i].id);
            // 	}
            // 	this_scenario[i].ambulances.erase(remove(
            // 		this_scenario[i].ambulances.begin(),
            // 		this_scenario[i].ambulances.end(), amb_finish),
            // 		this_scenario[i].ambulances.end());
            // }

            // if(debug){
            // 	fmt::print("\n");
            // }
          }
        }
        amb_finish =
            set_next_event(this_scenario, event_call, index_call, queue);
        if (debug) {
          double current_sum = 0;
          int current_solved = 0;
          for (size_t i = 0; i < (size_t)index_call; ++i) {
            auto &call = this_scenario[i];
            // cout << call << "\t" << waiting_on_scene[simu][i] << "\t";
            // cout << waiting_on_scene_penalized[simu][i] << "\t";
            // cout << which_ambulance[simu][i] << "\n";
            if (waiting_on_scene_penalized[simu][i] < GRB_INFINITY) {
              current_sum += waiting_on_scene_penalized[simu][i];
              current_solved += 1;
            }
          }
          double current_mean =
              (current_solved == 0) ? 0 : current_sum / current_solved;
          fmt::print(
              "Current sum = {}, Current mean = {:0f}, Current solved = {}\n",
              current_sum, current_mean, current_solved);
          cout << "==================================\n";
          // cin.get();
        }
      }
      if (debug) {
        cout << "Call\tWait On Scene(s)\tWait Penalized(s)\tAmb:\n";
      }
      bool has_inf = false;
      double sum_pen = accumulate(waiting_on_scene_penalized[simu].begin(),
                                  waiting_on_scene_penalized[simu].end(), 0.0);
      for (size_t i = 0; i < this_scenario.size(); ++i) {
        auto &call = this_scenario[i];
        if (debug) {
          cout << call << "\t" << waiting_on_scene[simu][i] << "\t";
          cout << waiting_on_scene_penalized[simu][i] << "\t";
          cout << which_ambulance[simu][i] << "\n";
        }

        if (waiting_on_scene_penalized[simu][i] > pow(10, 9) ||
            waiting_on_scene[simu][i] > pow(10, 9) ||
            waiting_on_scene_penalized[simu][i] == GRB_INFINITY ||
            waiting_on_scene[simu][i] == GRB_INFINITY ||
            which_ambulance[simu][i] == -1) {
          fmt::print("ERROR: Call {} at scenario {} unresolved\n", i, simu);
          fmt::print("Block indexes = {}\n", block_indexes[simu]);
          cin.get();
        }
        all_waiting_on_scene.push_back(waiting_on_scene[simu][i]);
        all_waiting_on_scene_penalized.push_back(
            waiting_on_scene_penalized[simu][i]);
        all_waiting_to_hospital.push_back(waiting_to_hospital[simu][i]);
      }
      // if(debug){
      // 	cin.get();
      // }
    }
    // cin.get();
    string type_return = (g_params.best_base) ? "_bb" : "";
    ofstream policy_results(
        fmt::format("results/two_stage/{}_{}_{}_{}{}.txt", g_params.amb_setup,
                    policy, ambulances.size(), ins.calls.size(), type_return),
        ios::out);
    // policy_results << g_params.n_scenarios << "\n";
    for (int k = 0; k < g_params.n_scenarios; ++k) {
      policy_results << waiting_on_scene[k].size() << "\n";
      for (size_t i = 0; i < waiting_on_scene[k].size(); ++i) {
        auto &call = ins.calls[k][i];
        double ws = waiting_on_scene[k][i];
        double wsp = waiting_on_scene_penalized[k][i];
        int which_amb = which_ambulance[k][i];
        double call_end = calls_end[k][i];
        // policy_results << call.id << " " << call.time << " " << call.priority
        // << " "; policy_results << call.id << " " << call.time << " " <<
        // call.priority << " "; policy_results << fmt::format("{:.6f} {:.6f} ",
        // call.location.first, call.location.second); policy_results << ws << "
        // " << wsp << " " << call_end << " " << which_amb << "\n";
        policy_results << which_amb << " " << ws << " " << wsp << " "
                       << call_end << "\n";
      }
    }
    policy_results.close();
    double sum_all_pen = accumulate(all_waiting_on_scene_penalized.begin(),
                                    all_waiting_on_scene_penalized.end(), 0.0);
    Stats stats(all_waiting_on_scene, all_waiting_on_scene_penalized,
                all_waiting_to_hospital);
    out_lines.push_back(fmt::format(
        "Policy {} Mean real {:.0f}, Max real {:.0f}, Mean pen = {:.0f}, Max "
        "pen = {:.0f}, real_std_dev {:.0f}, real_std_dev_avg {:.0f},  std_dev "
        "{:.0f}, std_dev_avg {:.0f}\n",
        policy, stats.mean_waiting_on_scene, stats.max_waiting_on_scene,
        stats.mean_waiting_on_scene_penalized,
        stats.max_waiting_on_scene_penalized, stats.actual_std_dev,
        stats.actual_std_dev_avg, stats.std_dev, stats.std_dev_avg));
    fmt::print("{}", out_lines.back());
    double min_run_time =
        *min_element(run_times_calls.begin(), run_times_calls.end());
    double avg_run_time =
        accumulate(run_times_calls.begin(), run_times_calls.end(), 0.0) /
        run_times_calls.size();
    double max_run_time =
        *max_element(run_times_calls.begin(), run_times_calls.end());
    out_lines.push_back(fmt::format("Run times: {}\t{}\t{}\n", min_run_time,
                                    avg_run_time, max_run_time));
    // fmt::print("{}", out_lines.back());
    // std::cin.get();
  }
  ofstream solver_run_times_file(
      fmt::format("results/run_times_{}.txt", ambulances.size()), ios::out);
  solver_run_times_file << solver_run_times.size() << "\n";
  for (auto &elem : solver_run_times) {
    int num_calls, index_solver;
    double run_time_elem;
    tie(num_calls, index_solver, run_time_elem) = elem;
    solver_run_times_file << num_calls << " " << index_solver << " "
                          << fmt::format("{:.6f}\n", run_time_elem);
  }
  solver_run_times_file.close();

  for (auto &line : out_lines) {
    fmt::print("{}", line);
  }
}

void PolicyTester::two_stage_tree() {
  vector<double> all_waiting_on_scene;
  vector<double> all_waiting_on_scene_penalized;
  vector<double> all_waiting_to_hospital;

  bool debug = false;

  for (size_t s = 0; s < ins.calls.size(); ++s) {
    for (size_t i = 0; i < ins.calls[s].size(); ++i) {
      for (size_t a = 0; a < ambulances.size(); ++a) {
        ins.calls[s][i].ambulances.push_back(a);
      }
    }
  }

  size_t num_calls_t = 100;  // calls by time slot
  size_t num_decisions =
      ambulances.size();  // number of closest ambulances available
  double t0 = ins.time_horizon[0] * ins.slot_duration;  // initial time of tests
  int T = ins.time_horizon.size();  // number of time horizons
  int nb_realizations = 5;
  vector<vector<vector<Call>>> my_scenarios(
      T, vector<vector<Call>>(nb_realizations, vector<Call>()));

  for (int scen = 0; scen < nb_realizations; ++scen) {
    size_t j = 0;
    for (int t = 0; t < T; ++t) {
      while (j < ins.calls[scen].size() &&
             ins.calls[scen][j].time <=
                 t0 + (t + 1) * ins.slot_duration) {  // 1800
        my_scenarios[t][scen].push_back(ins.calls[scen][j]);
        ++j;
      }
      if (num_calls_t < my_scenarios[t][scen].size()) {
        my_scenarios[t][scen].erase(my_scenarios[t][scen].begin() + num_calls_t,
                                    my_scenarios[t][scen].end());
      }
    }
  }

  default_random_engine gen;
  uniform_int_distribution<int> rand_scen(0, nb_realizations - 1);
  block_indexes = vector<vector<int>>(g_params.n_scenarios, vector<int>(T, 0));
  for (int simu = 0; simu < g_params.n_scenarios; ++simu) {
    for (int t = 0; t < T; ++t) {
      block_indexes[simu][t] = rand_scen(gen);
    }
    fmt::print("s = {}. block_indexes = {}\n", simu, block_indexes[simu]);
  }

  uniform_real_distribution<double> u(0, 1);
  vector<string> out_lines;
  bool has_cg = find(policies.begin(), policies.end(), "cg") != policies.end();
  for (auto &policy : policies) {
    // fmt::print("Policy: {}\n", policy);
    g_params.h_use_fixed_bases = (policy == "ordered" || policy == "district");
    bool free_at_base = (policy == "ordered" || policy == "district");
    run_times.clear();
    all_waiting_on_scene.clear();
    all_waiting_on_scene_penalized.clear();
    all_waiting_to_hospital.clear();
    for (int simu = 0; simu < g_params.n_scenarios; ++simu) {
      if (simu > 0) {
        fmt::print("{} {} scenarios completed\n", policy, simu);
      }
      default_random_engine u_gen;
      vector<Call> this_scenario;
      vector<int> this_scenario_indexes;
      for (int t = 0; t < 4; ++t) {  // 4
        int u = block_indexes[simu][t];
        // int u = 0;
        this_scenario_indexes.push_back(u);
        auto my_end = (my_scenarios[t][u].size() > num_calls_t)
                          ? my_scenarios[t][u].begin() + num_calls_t
                          : my_scenarios[t][u].end();
        this_scenario.insert(this_scenario.end(), my_scenarios[t][u].begin(),
                             my_end);
      }
      // fmt::print("This scenario:\n");
      for (size_t i = 0; i < this_scenario.size(); ++i) {
        this_scenario[i].id = i;
        this_scenario[i].ambulances.clear();
        for (size_t a = 0; a < ambulances.size(); ++a) {
          this_scenario[i].ambulances.push_back(a);
        }
        // cout << this_scenario[i] << "\n";
        // fmt::print("{},{},{},{}\n", this_scenario[i].id,
        // this_scenario[i].priority, this_scenario[i].location.first,
        // this_scenario[i].location.second);
      }

      waiting_on_scene[simu] =
          vector<double>(this_scenario.size(), GRB_INFINITY);
      waiting_on_scene_penalized[simu] =
          vector<double>(this_scenario.size(), GRB_INFINITY);
      waiting_to_hospital[simu] =
          vector<double>(this_scenario.size(), GRB_INFINITY);
      calls_end[simu] = vector<double>(this_scenario.size(), GRB_INFINITY);
      which_ambulance[simu] = vector<int>(this_scenario.size(), -1);

      auto nearest_base = get_nearest_base(this_scenario);
      ambulances = ins.ambulances;
      if (policy == "district") {
        for (size_t a = 0; a < ambulances.size(); ++a) {
          ambulances[a].base_location =
              ins.bases[data[4].ambulance_district[a]];
        }
      }
      size_t calls_attended = 0;
      int event_call = 1;
      int index_call = 0;
      time = this_scenario[0].time;
      int amb_finish = -1;

      std::vector<int> queue;
      queue.reserve(this_scenario.size());
      // int iter_count = 0;
      while (calls_attended < this_scenario.size()) {
        vector<vector<Call>> future_scenarios = get_future_scenarios(
            time, this_scenario, nb_realizations, T, my_scenarios);

        for (auto &amb : ambulances) {
          // if ambulance arrived at base it's now free to answer all calls.
          if (amb.arrival_time_at_b_last_trip <= time) {
            for (size_t i = 0; i < this_scenario.size(); ++i) {
              this_scenario[i].ambulances.push_back(amb.id);
            }
          }
        }
        // travel.set_delay_param(u(u_gen));
        //  travel.set_delay_param(-1);
        if (debug) {
          fmt::print(
              "Time = {}, index_call = {}, ev = {}, calls_attended = {} "
              "/ {}, queue = {}, policy = {}\n",
              time, index_call, event_call, calls_attended,
              this_scenario.size(), queue, policy);
          vector<double> queue_times;
          for (auto i : queue) {
            queue_times.push_back(this_scenario[i].time);
          }
          fmt::print("Times: {}\n", queue_times);
        }
        if (event_call == 1 || (event_call == 0 && free_at_base)) {
          if (event_call == 1 && time <= this_scenario.back().time) {
            queue.push_back(index_call);
          }
          std::vector<int> queue_aux;
          queue_aux.reserve(queue.size());
          for (size_t i = 0; i < queue.size(); ++i) {
            auto t0 = std::chrono::high_resolution_clock::now();
            auto &call = this_scenario[queue[i]];
            double min_time = GRB_INFINITY;
            int index_amb = -1;
            double min_waiting_on_scene = GRB_INFINITY;
            double min_waiting_to_hospital = GRB_INFINITY;

            if (debug) {
              cout << "Call " << call << "\n";
            }
            if (!has_cg) {  // if model not involved, run all possible decisions
              for (size_t k = 0; k < ambulances.size(); ++k) {
                auto &ambulance = ambulances[k];
                if (can_answer(ambulance, call) &&
                    ambulance.arrival_time_at_f_last_trip <= time) {
                  if (debug) {
                    cout << "Running Amb " << ambulance.id << " Type "
                         << ambulance.type << "\n";
                  }
                  auto result = get_waiting_time_tree(
                      ambulance.id, call.id, policy, future_scenarios,
                      this_scenario, time, queue, nearest_base[queue[i]]);
                  if (debug) {
                    cout << "\t" << ambulance << " => "
                         << result.mean_total_time;
                    cout << "\n";
                  }
                  if (result.mean_total_time < min_time) {
                    index_amb = k;
                    min_time = result.mean_total_time;
                  }
                }
              }
            } else {  // else, choose num_decisions nearest ambulances to test
                      // the second stage
              vector<pair<double, int>> nearest_ambulances;
              for (size_t k = 0; k < ambulances.size(); ++k) {
                auto &ambulance = ambulances[k];
                if (can_answer(ambulance, call) &&
                    ambulance.arrival_time_at_f_last_trip <= time) {
                  double resp_time = penalized_response_time(
                      travel.get_response_time(ambulance, call, time),
                      ambulance.type, call.priority);
                  nearest_ambulances.push_back(make_pair(resp_time, k));
                }
              }
              std::sort(nearest_ambulances.begin(), nearest_ambulances.end());
              for (size_t k = 0;
                   k < num_decisions && k < nearest_ambulances.size(); ++k) {
                auto &ambulance = ambulances[nearest_ambulances[k].second];
                if (debug) {
                  cout << "Testing " << ambulance << " for " << call << "\n";
                }
                auto result = (true) ? get_waiting_time_tree(
                                           ambulance.id, call.id, policy,
                                           future_scenarios, this_scenario,
                                           time, queue, nearest_base[queue[i]])
                                     : get_waiting_time_model(
                                           ambulance.id, call.id, policy,
                                           future_scenarios, this_scenario,
                                           time, queue, nearest_base[queue[i]]);
                if (debug) {
                  cout << "\t" << ambulance << " => " << result.mean_total_time;
                  cout << "\n";
                }
                if (result.mean_total_time < min_time) {
                  index_amb = nearest_ambulances[k].second;
                  min_time = result.mean_total_time;
                }
              }
            }
            if (debug) {
              cout << "Testing queue for " << call << "\n";
            }
            auto result =
                (policy != "cg")
                    ? get_waiting_time_tree(-1, call.id, policy,
                                            future_scenarios, this_scenario,
                                            time, queue, nearest_base[queue[i]])
                    : get_waiting_time_model(
                          -1, call.id, policy, future_scenarios, this_scenario,
                          time, queue, nearest_base[queue[i]]);
            if (debug) {
              cout << "\tQueue =>" << result.mean_total_time << "\n";
            }
            if (result.mean_total_time < min_time) {
              index_amb = -1;
              min_time = result.mean_total_time;
            }

            if (index_amb >= 0) {
              auto &best_amb = ambulances[index_amb];
              min_time = travel.get_response_time(best_amb, call, time, true);
              double waiting_on_scene_i = time + min_time - call.time;
              double waiting_to_hospital_i = ambulances[index_amb].answer_call(
                  call, travel, ins, time, min_time, nearest_base[queue[i]]);
              if (free_at_base) {
                best_amb.free_location = best_amb.base_location;
                best_amb.arrival_time_at_f_last_trip =
                    best_amb.arrival_time_at_b_last_trip;
              }
              which_ambulance[simu][queue[i]] = index_amb;
              waiting_on_scene[simu][queue[i]] = waiting_on_scene_i;
              waiting_on_scene_penalized[simu][queue[i]] =
                  penalized_response_time(waiting_on_scene_i,
                                          ambulances[index_amb].type,
                                          call.priority);
              if (debug) {
                Location current_location;
                string type_current_location = "-";
                if (best_amb.arrival_time_at_b_last_trip <= time) {
                  type_current_location = "b";
                  current_location = best_amb.base_location;
                } else if (best_amb.arrival_time_at_f_last_trip < time) {
                  type_current_location = "r";
                  current_location = travel.ambulance_position(best_amb, time);
                } else if (abs(best_amb.arrival_time_at_f_last_trip - time) <
                           g_params.EPS) {
                  type_current_location = "h";
                  current_location = best_amb.free_location;
                } else {
                  type_current_location = "c";
                  current_location = travel.position_between_origin_destination(
                      best_amb.last_origin_location, best_amb.call->location,
                      best_amb.departure_time, time,
                      best_amb.arrival_time_at_c_last_trip, best_amb);
                }
                cout << "Chosen " << best_amb
                     << fmt::format(
                            " | {:.6f} {:.6f}({})", current_location.first,
                            current_location.second, type_current_location)
                     << " pt ";
                cout << waiting_on_scene_penalized[simu][queue[i]] << "\n";
              }
              waiting_to_hospital[simu][queue[i]] = waiting_to_hospital_i;
              calls_attended++;
            } else {
              queue_aux.push_back(queue[i]);
              if (debug) {
                cout << "Chose Queue\n";
              }
            }
            auto dt = std::chrono::high_resolution_clock::now();
            run_times.push_back(
                std::chrono::duration_cast<chrono::nanoseconds>(dt - t0)
                    .count() /
                pow(10, 9));
          }
          queue.clear();
          for (auto i : queue_aux) {
            queue.push_back(i);
          }
        } else {
          // Event_call == 0, get ambulance that just finished.
          if (amb_finish == -1) {
            fmt::print("WEIRD!!! event_call = 0 but no return\n");
            fmt::print("Runtimes.size() = {}\n", run_times.size());
            fmt::print("policy = {}, s = {}, time = {}\n", policy, simu, time);
            double min_run_time =
                *min_element(run_times.begin(), run_times.end());
            double avg_run_time =
                accumulate(run_times.begin(), run_times.end(), 0.0) /
                run_times.size();
            double max_run_time =
                *max_element(run_times.begin(), run_times.end());
            fmt::print("Run times: {}\t{}\t{}\n", min_run_time, avg_run_time,
                       max_run_time);
            cin.get();
          }
          auto &amb = ambulances[amb_finish];
          std::vector<int> queue_aux;
          queue_aux.reserve(queue.size());
          if (debug) {
            cout << "Return " << amb << "\n";
          }
          int ind_call = -1;
          int ind_base = -1;
          double min_time = GRB_INFINITY;
          double min_waiting_on_scene = GRB_INFINITY;
          double min_waiting_to_hospital = GRB_INFINITY;
          for (size_t i = 0; i < queue.size(); ++i) {
            auto t0 = std::chrono::high_resolution_clock::now();
            auto &call = this_scenario[queue[i]];
            if (debug && !can_answer(amb, call)) {
              cout << "Call " << call << " | Can't be answered by " << amb.id
                   << "\n";
            }
            if (can_answer(amb, call)) {
              if (debug) {
                cout << "Testing " << amb << " for " << call << "\n";
              }
              auto result =
                  (true)
                      ? get_waiting_time_tree(
                            amb_finish, call.id, policy, future_scenarios,
                            this_scenario, time, queue, nearest_base[queue[i]])
                      : get_waiting_time_model(
                            amb_finish, call.id, policy, future_scenarios,
                            this_scenario, time, queue, nearest_base[queue[i]]);
              if (debug) {
                cout << "\tAmb " << amb.id << " Call " << call << " => "
                     << result.mean_total_time;
                cout << "\n";
              }
              if (result.mean_total_time < min_time) {
                ind_call = queue[i];
                min_time = result.mean_total_time;
              }
            }
            auto dt = std::chrono::high_resolution_clock::now();
            run_times.push_back(
                std::chrono::duration_cast<chrono::nanoseconds>(dt - t0)
                    .count() /
                pow(10, 9));
          }
          if (!free_at_base && !has_cg) {
            for (size_t b = 0; b < ins.bases.size(); ++b) {
              auto result = get_waiting_time_tree_return(
                  amb_finish, b, policy, future_scenarios, this_scenario, time,
                  queue);
              if (debug) {
                cout << "\tAmb " << amb.id << " Base " << b << " => "
                     << result.mean_total_time << "\n";
              }
              if (result.mean_total_time < min_time) {
                min_time = result.mean_total_time;
                ind_call = -1;
                ind_base = b;
              }
            }
          } else if (!free_at_base) {
            vector<pair<double, int>> nearest_bases;
            for (size_t b = 0; b < ins.bases.size(); ++b) {
              double return_time =
                  travel.travel_time(amb.free_location, ins.bases[b], amb);
              nearest_bases.push_back(make_pair(return_time, b));
            }
            std::sort(nearest_bases.begin(), nearest_bases.end());
            for (size_t k = 0; k < num_decisions && k < nearest_bases.size();
                 ++k) {
              int b = nearest_bases[k].second;
              if (debug) {
                cout << "Testing " << amb << " for base " << b << "\n";
              }
              auto result = (true)
                                ? get_waiting_time_tree_return(
                                      amb_finish, b, policy, future_scenarios,
                                      this_scenario, time, queue)
                                : get_waiting_time_model_return(
                                      amb_finish, b, policy, future_scenarios,
                                      this_scenario, time, queue);
              if (debug) {
                cout << "\tAmb " << amb.id << " Base " << b << " => "
                     << result.mean_total_time << "\n";
              }
              if (result.mean_total_time < min_time) {
                min_time = result.mean_total_time;
                ind_call = -1;
                ind_base = b;
              }
            }
          }

          if (ind_call >= 0) {
            auto &call = this_scenario[ind_call];
            double min_time = travel.get_response_time(amb, call, time, true);
            double waiting_on_scene_i = time + min_time - call.time;
            calls_attended++;
            // iter_count = 0;
            double waiting_to_hospital_i = amb.answer_call(
                call, travel, ins, time, min_time, nearest_base[ind_call]);
            if (free_at_base) {
              amb.free_location = amb.base_location;
              amb.arrival_time_at_f_last_trip = amb.arrival_time_at_b_last_trip;
            }
            which_ambulance[simu][ind_call] = amb_finish;
            waiting_on_scene[simu][ind_call] = waiting_on_scene_i;
            waiting_on_scene_penalized[simu][ind_call] =
                penalized_response_time(waiting_on_scene_i, amb.type,
                                        call.priority);
            if (debug) {
              cout << "Chosen call " << call << " ";
              cout << waiting_on_scene_penalized[simu][ind_call] << "\n";
            }
            waiting_to_hospital[simu][ind_call] = waiting_to_hospital_i;
            queue.erase(remove(queue.begin(), queue.end(), ind_call),
                        queue.end());
          } else if (ind_base >= 0) {
            // cout << "Chose Return\n";
            auto base = ins.bases[ind_base];
            amb.base_location = base;
            if (!free_at_base) {
              amb.arrival_time_at_b_last_trip =
                  time +
                  travel.travel_time(amb.free_location, amb.base_location, amb);
            } else {
              amb.arrival_time_at_b_last_trip = time;
            }

            if (debug) {
              fmt::print("Chosen Base {}\n", ind_base);
            }

            // If ambulance must return to base, it can't answer any calls in
            // queue.
            if (debug) {
              fmt::print("Amb {} forbidden to answer to: ", amb_finish);
            }
            for (auto i : queue) {
              if (debug) {
                fmt::print("{} ", this_scenario[i].id);
              }
              this_scenario[i].ambulances.erase(
                  remove(this_scenario[i].ambulances.begin(),
                         this_scenario[i].ambulances.end(), amb_finish),
                  this_scenario[i].ambulances.end());
            }

            if (debug) {
              fmt::print("\n");
            }
          }
        }
        amb_finish =
            set_next_event(this_scenario, event_call, index_call, queue);
        if (debug) {
          cout << "==================================\n";
          // cin.get();
        }
      }
      if (debug) {
        cout << "Call\tWait On Scene(s)\tWait Penalized(s)\tAmb:\n";
      }

      for (size_t i = 0; i < this_scenario.size(); ++i) {
        auto &call = this_scenario[i];
        if (debug) {
          cout << call << "\t" << waiting_on_scene[simu][i] << "\t";
          cout << waiting_on_scene_penalized[simu][i] << "\t";
          cout << which_ambulance[simu][i] << "\n";
        }
        all_waiting_on_scene.push_back(waiting_on_scene[simu][i]);
        all_waiting_on_scene_penalized.push_back(
            waiting_on_scene_penalized[simu][i]);
        all_waiting_to_hospital.push_back(waiting_to_hospital[simu][i]);
      }

      if (debug) {
        cin.get();
      }
    }
    // cin.get();
    ofstream policy_results(
        fmt::format("results/two_stage_{}_a{}_s{}.txt", policy,
                    ambulances.size(), ins.calls.size()),
        ios::out);
    policy_results << all_waiting_on_scene.size() << "\n";
    for (size_t k = 0; k < all_waiting_on_scene.size(); ++k) {
      double ws = all_waiting_on_scene[k];
      double wsp = all_waiting_on_scene_penalized[k];
      policy_results << which_ambulance[0][k] << " " << ws << "\t" << wsp
                     << "\n";
    }
    policy_results.close();
    Stats stats(all_waiting_on_scene, all_waiting_on_scene_penalized,
                all_waiting_to_hospital);
    out_lines.push_back(fmt::format(
        "Policy {} Mean real {:.1f}, Max real {:.1f}, Mean pen = {:.1f}, Max "
        "pen = {:.1f}, real_std_dev {:.1f}, real_std_dev_avg {:.1f},  std_dev "
        "{:.1f}, std_dev_avg {:.1f}\n",
        policy, stats.mean_waiting_on_scene, stats.max_waiting_on_scene,
        stats.mean_waiting_on_scene_penalized,
        stats.max_waiting_on_scene_penalized, stats.actual_std_dev,
        stats.actual_std_dev_avg, stats.std_dev, stats.std_dev_avg));
    fmt::print("{}", out_lines.back());
    double min_run_time = *min_element(run_times.begin(), run_times.end());
    double avg_run_time =
        accumulate(run_times.begin(), run_times.end(), 0.0) / run_times.size();
    double max_run_time = *max_element(run_times.begin(), run_times.end());
    out_lines.push_back(fmt::format("Run times: {}\t{}\t{}\n", min_run_time,
                                    avg_run_time, max_run_time));
    fmt::print("{}", out_lines.back());
    // std::cin.get();
  }

  ofstream solver_run_times_file(
      fmt::format("results/run_times_{}.txt", ambulances.size()), ios::out);
  solver_run_times_file << solver_run_times.size() << "\n";
  for (auto &elem : solver_run_times) {
    int num_calls, index_solver;
    double run_time_elem;
    tie(num_calls, index_solver, run_time_elem) = elem;
    solver_run_times_file << num_calls << " " << index_solver << " "
                          << fmt::format("{:.6f}\n", run_time_elem);
  }
  solver_run_times_file.close();

  for (auto &line : out_lines) {
    fmt::print("{}", line);
  }
}

SecondStageWaitingTime PolicyTester::get_waiting_time_tree(
    int amb_id, int call_id, const string &policy,
    vector<vector<Call>> &future_scenarios, vector<Call> &this_scenario,
    double time, std::vector<int> &queue, int base) {
  Call *call = NULL;
  if (call_id >= 0) {
    call = &this_scenario[call_id];
  }

  bool debug = false;

  vector<Ambulance> temp_ambulances = ambulances;
  double min_time = 0;
  double waiting_time = 0;
  double waiting_to_hospital = 0;
  double pen_waiting_time = 0;
  if (amb_id >= 0 && call_id >= 0) {
    Ambulance &amb = temp_ambulances[amb_id];
    min_time = travel.get_response_time(amb, *call, time, true);
    waiting_to_hospital =
        amb.answer_call(*call, travel, ins, time, min_time, base);
    waiting_time = time + min_time - call->time;
    pen_waiting_time =
        penalized_response_time(waiting_time, amb.type, call->priority);
    if (debug) {
      fmt::print(
          "Amb {}, Type {} dispatched to {}, pen_w = {} and will be "
          "free at {}\n",
          amb.id, amb.type, call_id, pen_waiting_time,
          amb.arrival_time_at_f_last_trip);
    }
  }

  double sum_total = 0;
  int num_total = 0;

  vector<int> busy_ambulances;
  if (amb_id == -1) {
    for (auto &amb : temp_ambulances) {
      if (can_answer(amb, *call) && amb.arrival_time_at_f_last_trip > time) {
        busy_ambulances.push_back(amb.id);
      }
    }

    if (busy_ambulances.size() == 0) {
      return SecondStageWaitingTime{0, 0, GRB_INFINITY};
    }
  }

  bool is_call_in_queue =
      find(queue.begin(), queue.end(), call_id) != queue.end();
  vector<int> queue_aux = queue;
  if (amb_id == -1 && !is_call_in_queue) {
    queue_aux.push_back(call_id);
  }

  for (size_t j = 0; j < future_scenarios.size(); ++j) {
    auto future_calls = future_scenarios[j];
    vector<Ambulance> aux_ambulances = temp_ambulances;

    if (!queue_aux.empty()) {
      // if queue isn't empty, calls of queue are added to future scenarios.
      for (int i = queue_aux.size() - 1; i >= 0; --i) {
        // If testing for call_id in queue or a call in queue is not the current
        // call.
        if (amb_id == -1 || queue_aux[i] != call_id) {
          future_calls.insert(future_calls.begin(),
                              this_scenario[queue_aux[i]]);
          if (amb_id == -1 && queue_aux[i] == call_id) {
            // If the current call must go to queue,
            // all free ambulances are excluded from attending it.
            // and will only be possible to attend it after arriving at the
            // base.
            for (auto &amb : aux_ambulances) {
              if (amb.arrival_time_at_f_last_trip <= time) {
                future_calls[0].ambulances.erase(
                    remove(future_calls[0].ambulances.begin(),
                           future_calls[0].ambulances.end(), amb.id),
                    future_calls[0].ambulances.end());
              }
            }
          }
        }
      }
    }

    if (!future_calls.empty()) {
      std::sort(future_calls.begin(), future_calls.end(),
                [&](const Call &lhs, const Call &rhs) {
                  return lhs.time < rhs.time;
                });
      for (size_t i = 0; i < future_calls.size(); ++i) {
        future_calls[i].id = i;
        // if(debug){
        // 	fmt::print("Call {}, time {}, ambs = {}\n", future_calls[i].id,
        // future_calls[i].time, 		future_calls[i].ambulances);
        // }
      }
      auto solver =
          get_solver(policy, future_calls, aux_ambulances, travel, 4, time);
      if (debug) {
        solver->set_debug_mode(true);
      }
      solver->run();

      if (policy == "cg" || policy == "enumerate") {
        // fmt::print("|C| = {}, index = {}, time = {}\n", future_calls.size(),
        // solver->index_solver,  solver->run_time);
        solver_run_times.push_back(make_tuple(
            future_calls.size(), solver->index_solver, solver->run_time));
      } else {
        double scenario_run_time =
            accumulate(solver->run_times.begin(), solver->run_times.end(), 0.0);
        int index_solver = -1;
        for (size_t i = 0; i < policies.size(); ++i) {
          if (policy == policies[i]) {
            index_solver = i;
          }
        }
        solver_run_times.push_back(
            make_tuple(future_calls.size(), index_solver, scenario_run_time));
      }

      double sum_scenario = 0;
      int num_scenario = 0;

      if (policy != "cg" && policy != "enumerate") {
        sum_scenario =
            accumulate(solver->waiting_on_scene_penalized.begin(),
                       solver->waiting_on_scene_penalized.end(), 0.0);
        // Total waiting time not really needed
        // for(size_t i = 0; i < solver->waiting_to_hospital.size(); ++i){
        // 	if(solver->which_ambulance[i] == -1){
        // 		fmt::print("ERROR: unsolved call {} {} w1 = {}\n", i,
        // policy); 		solver->print_results(); cin.get();
        // 	}
        // 	int amb_type = ambulances[solver->which_ambulance[i]].type;
        // 	int call_type = future_calls[i].priority;
        // 	sum_scenario +=
        // penalized_response_type((solver->waiting_to_hospital[i] +
        // solver->calls[i].time_at_hospital), amb_type, call_type);
        // }
        // sum_scenario += accumulate(solver->waiting_to_hospital.begin(),
        // 	solver->waiting_to_hospital.end(), 0.0);

        sum_total += sum_scenario;
        num_scenario = solver->waiting_on_scene_penalized.size();
        num_total += num_scenario;

        if (sum_scenario > pow(10, 9)) {
          fmt::print(
              "Warning! scenario {} has sum_scenario = {}, policy = {}\n", j,
              sum_scenario, policy);
          solver =
              get_solver(policy, future_calls, aux_ambulances, travel, 4, time);
          solver->set_debug_mode(true);
          solver->run();
          solver->print_results();
          cin.get();
        }

        if (debug) {
          solver->print_results();
          fmt::print("\t\tScenario {}, sum_scenario = {}\n", j, sum_scenario);
        }

      } else {
        if (amb_id >= 0 && call_id >= 0) {
          auto &amb = temp_ambulances[amb_id];
          // sum_total += solver->obj + penalized_response_time(
          //                                waiting_time + waiting_to_hospital +
          //                                    call->time_at_hospital,
          //                                amb.type, call->priority);
          sum_total += solver->obj;
        } else {
          sum_total += solver->obj;
        }
        ++num_total;
      }
    } else {
      // if(debug){
      // 	fmt::print("Entrou call empty!\n");
      // }
      if (amb_id >= 0) {
        // return
        // SecondStageWaitingTime{0,0,waiting_time*ins.penalty_matrix[ambulances[amb_id].type][call->priority]};
        return SecondStageWaitingTime{
            0, 0,
            penalized_response_time(waiting_time, ambulances[amb_id].type,
                                    call->priority)};
      } else {
        return {0, 0, pen_waiting_time};
      }
    }
  }

  double f_x0 = 0;
  bool testing_call_on_queue = amb_id == -1 && call_id != -1;
  if (!testing_call_on_queue) {
    auto &amb_in = ambulances[amb_id];
    f_x0 = penalized_response_time(waiting_time, amb_in.type,
                                   this_scenario[call_id].priority);
  }

  double mean_total = (sum_total / future_scenarios.size());
  // if(debug){
  // 	fmt::print("f_x0 = {}, mean_total = {}\n", f_x0, mean_total);
  // }

  // double result_total = (num_total == 0) ? f_x0 : mean_total + (f_x0 -
  // mean_total) / (num_total + 1);
  double result_total = (num_total == 0) ? 0 : f_x0 + mean_total;
  // fmt::print("\t\tresult_total = {}\n", result_total);
  return SecondStageWaitingTime{0, 0, result_total};
}

SecondStageWaitingTime PolicyTester::get_waiting_time_tree_return(
    int amb_id, int base, const string &policy,
    vector<vector<Call>> &future_scenarios, vector<Call> &this_scenario,
    double time, std::vector<int> &queue, double min_time_to_free) {
  bool debug = false;
  if (debug) {
    fmt::print("Entered debug mode\n");
  }

  vector<Ambulance> temp_ambulances = ambulances;
  auto &amb = temp_ambulances[amb_id];
  amb.base_location = ins.bases[base];
  double time_amb_to_base =
      travel.travel_time(amb.free_location, amb.base_location, amb);
  amb.arrival_time_at_b_last_trip = time + time_amb_to_base;

  double min_time = 0;
  double waiting_to_hospital = 0;
  double sum_total = 0;
  int num_total = 0;

  for (size_t j = 0; j < future_scenarios.size(); ++j) {
    auto future_calls = future_scenarios[j];
    vector<Ambulance> aux_ambulances = temp_ambulances;
    for (int i = queue.size() - 1; i >= 0; --i) {
      future_calls.insert(future_calls.begin(), this_scenario[queue[i]]);
      // By returning the current ambulance to base,
      // we forbid it to answer any calls in queue:
      future_calls[0].ambulances.erase(
          remove(future_calls[0].ambulances.begin(),
                 future_calls[0].ambulances.end(), amb_id),
          future_calls[0].ambulances.end());
    }
    if (debug) {
      fmt::print("Processing Scenarios {}: num_calls = {}\n", j,
                 future_calls.size());
    }
    if (!future_calls.empty()) {
      std::sort(future_calls.begin(), future_calls.end(),
                [&](const Call &lhs, const Call &rhs) {
                  return lhs.time < rhs.time;
                });
      for (size_t i = 0; i < future_calls.size(); ++i) {
        future_calls[i].id = i;
      }
      auto solver =
          get_solver(policy, future_calls, aux_ambulances, travel, 4, time);
      if (debug) {
        solver->set_debug_mode(true);
      }
      solver->run();
      if (debug) {
        solver->set_debug_mode(false);
      }
      if (policy == "enumerate" || policy == "cg") {
        int index_solver = (policy == "enumerate") ? 4 : solver->index_solver;
        solver_run_times.push_back(
            make_tuple(future_calls.size(), index_solver, solver->run_time));
      }
      double sum_scenario = 0;
      int num_scenario = 0;
      if (policy != "cg" && policy != "enumerate") {
        sum_scenario =
            accumulate(solver->waiting_on_scene_penalized.begin(),
                       solver->waiting_on_scene_penalized.end(), 0.0);
        // for(size_t i = 0; i < solver->waiting_to_hospital.size(); ++i){
        // 	int amb_type = ambulances[solver->which_ambulance[i]].type;
        // 	int call_type = future_calls[i].priority;
        // 	sum_scenario +=
        // penalized_response_time(solver->waiting_to_hospital[i] +
        // solver->calls[i].time_at_hospital, amb_type, call_type);
        // }
        sum_total += sum_scenario;
        // sum_total += accumulate(solver->waiting_to_hospital.begin(),
        // 	solver->waiting_to_hospital.end(), 0.0);
        num_scenario = solver->waiting_on_scene_penalized.size();
        num_total += num_scenario;

        if (sum_scenario > pow(10, 9)) {
          fmt::print(
              "Warning! scenario {} has sum_scenario = {}, policy = {}\n", j,
              sum_scenario, policy);
          solver =
              get_solver(policy, future_calls, aux_ambulances, travel, 4, time);
          solver->set_debug_mode(true);
          solver->run();
          solver->print_results();
          cin.get();
        }

        if (debug) {
          solver->print_results();
          fmt::print("\t\tScenario {}, sum_scenario = {}\n", j, sum_scenario);
        }
      } else {
        sum_total += solver->obj;
        ++num_total;
      }
    } else {
      // If no future calls, return to closest base
      return {0, 0, time_amb_to_base};
    }
  }
  double f_x0 = 0;
  double result_total =
      (num_total == 0) ? 0 : f_x0 + sum_total / future_scenarios.size();
  return SecondStageWaitingTime{0, 0, result_total};
}

vector<vector<Call>> PolicyTester::get_future_scenarios(
    double time, vector<Call> &this_scenario, int nb_realizations, int T,
    vector<vector<vector<Call>>> &my_scenarios) {
  int t_slot = 0;
  double t0 = ins.time_horizon[0] * ins.slot_duration;
  double slot_duration = 1800;  // T*1800 / 4
  using namespace fmt;

  if (time <= t0 + slot_duration) {
    t_slot = 0;
  } else if (time <= t0 + 2 * slot_duration) {
    t_slot = 1;
  } else if (time <= t0 + 3 * slot_duration) {
    t_slot = 2;
  } else {
    t_slot = 3;
  }

  size_t lb_scen = 0;
  while (lb_scen < this_scenario.size() &&
         this_scenario[lb_scen].time <= time) {
    ++lb_scen;
  }
  size_t ub_scen = lb_scen;
  while (ub_scen < this_scenario.size() &&
         this_scenario[ub_scen].time <= t0 + (t_slot + 1) * slot_duration) {
    ++ub_scen;
  }
  size_t nb_scen = 0;
  if (time >= t0 + 3 * slot_duration) {  // Last time slot, no future scenarios.
    nb_scen = 1;
    vector<vector<Call>> future_scenarios(nb_scen, vector<Call>());
    int ind = 0;
    for (size_t i = lb_scen; i < ub_scen; ++i) {
      future_scenarios[ind].push_back(this_scenario[i]);
    }
    return future_scenarios;
  }

  size_t num_calls_t = 100;           // former 20
  size_t num_future_scenarios = 100;  // former 15 - 25
  default_random_engine gen(1);
  // uniform_int_distribution<int> rand_scen(0,nb_realizations-1);

  // FIXME Remove
  //  vector<int> this_block_indexes =  {1,2,4,4};
  //  nb_scen = 1;
  //  vector<vector<Call>> future_scenarios(nb_scen, vector<Call>());
  //  for(size_t i = lb_scen; i < ub_scen; ++i){
  //  	future_scenarios[0].push_back(this_scenario[i]);
  //  }
  //  if(time < t0 + slot_duration){
  //  	for(int t = 1; t < T; ++t){
  //  		future_scenarios[0].insert(future_scenarios[0].end(),
  //  my_scenarios[t][this_block_indexes[t]].begin(),
  //  			my_scenarios[t][this_block_indexes[t]].end());
  //  	}
  //  }else if(time < t0 + 2*slot_duration){
  //  	future_scenarios[0].insert(future_scenarios[0].end(),
  //  my_scenarios[2][this_block_indexes[2]].begin(),
  //  			my_scenarios[2][this_block_indexes[2]].end());
  //  	future_scenarios[0].insert(future_scenarios[0].end(),
  //  my_scenarios[3][this_block_indexes[3]].begin(),
  //  			my_scenarios[3][this_block_indexes[3]].end());
  //  }else if(time < t0 + 3*slot_duration){
  //  	future_scenarios[0].insert(future_scenarios[0].end(),
  //  my_scenarios[3][this_block_indexes[3]].begin(),
  //  			my_scenarios[3][this_block_indexes[3]].end());
  //  }else{
  //  	return vector<vector<Call>>();
  //  }
  //  for(size_t i = 0; i < future_scenarios.size(); ++i){
  //  	future_scenarios[0][i].id = i;
  //  }
  //  return future_scenarios;

  if (time < t0 + slot_duration) {
    t_slot = 0;
    // nb_scen = 125;
    nb_scen = 25;
    uniform_int_distribution<int> rand_scen(0, nb_realizations - 1);
    vector<vector<Call>> future_scenarios(nb_scen, vector<Call>());
    int ind = 0;

    for (size_t ind = 0; ind < nb_scen; ++ind) {
      for (size_t i = lb_scen; i < ub_scen; ++i) {
        future_scenarios[ind].push_back(this_scenario[i]);
      }

      for (int t = 1; t < T; ++t) {
        int u = rand_scen(gen);
        auto &current_my_scenarios = my_scenarios[t][u];
        // auto& current_my_scenarios = my_scenarios[t][this_block_indexes[t]];
        auto my_end = (current_my_scenarios.size() > num_calls_t)
                          ? current_my_scenarios.begin() + num_calls_t
                          : current_my_scenarios.end();
        // auto my_end = my_scenarios[t][u].end();
        future_scenarios[ind].insert(future_scenarios[ind].end(),
                                     current_my_scenarios.begin(), my_end);
        // my_scenarios[t][u].begin(), my_scenarios[t][u].end());
      }
    }

    for (size_t i = 0; i < nb_scen; ++i) {
      for (size_t j = 0; j < future_scenarios[i].size(); ++j) {
        future_scenarios[i][j].id = j;
      }
    }
    future_scenarios.erase(
        future_scenarios.begin() + min(nb_scen, num_future_scenarios),
        future_scenarios.end());
    // fmt::print("future_scenarios\n");
    // for(size_t tscen = 0; tscen < nb_scen; ++tscen){
    // 	fmt::print("Future {}:\n", tscen);
    // 	for(auto& call: future_scenarios[tscen]){
    // 		cout << call <<  "\n";
    // 	}
    // 	cin.get();
    // }
    return future_scenarios;
  } else if (time < t0 + 2 * slot_duration) {
    t_slot = 1;
    nb_scen = 25;
    uniform_int_distribution<int> rand_scen(0, nb_realizations - 1);
    vector<vector<Call>> future_scenarios(nb_scen, vector<Call>());
    int ind = 0;
    for (int u = 0; u < nb_realizations; ++u) {
      for (int v = 0; v < nb_realizations; ++v) {
        for (size_t i = lb_scen; i < this_scenario.size(); ++i) {
          future_scenarios[ind].push_back(this_scenario[i]);
        }
        auto &current_my_scenarios = my_scenarios[2][u];
        // auto& current_my_scenarios = my_scenarios[2][this_block_indexes[1]];
        auto my_end = (current_my_scenarios.size() > num_calls_t)
                          ? current_my_scenarios.begin() + num_calls_t
                          : current_my_scenarios.end();
        future_scenarios[ind].insert(future_scenarios[ind].end(),
                                     current_my_scenarios.begin(), my_end);

        auto &current_my_scenarios_v = my_scenarios[3][v];
        // auto& current_my_scenarios_v =
        // my_scenarios[3][this_block_indexes[2]];
        my_end = (current_my_scenarios_v.size() > num_calls_t)
                     ? current_my_scenarios_v.begin() + num_calls_t
                     : current_my_scenarios_v.end();
        future_scenarios[ind].insert(future_scenarios[ind].end(),
                                     current_my_scenarios_v.begin(), my_end);
        ++ind;
      }
    }

    for (size_t i = 0; i < nb_scen; ++i) {
      for (size_t j = 0; j < future_scenarios[i].size(); ++j) {
        future_scenarios[i][j].id = j;
      }
    }
    future_scenarios.erase(
        future_scenarios.begin() + min(nb_scen, num_future_scenarios),
        future_scenarios.end());
    return future_scenarios;
  } else if (time < t0 + 3 * slot_duration) {
    t_slot = 2;
    nb_scen = 5;
    uniform_int_distribution<int> rand_scen(0, nb_realizations - 1);
    vector<vector<Call>> future_scenarios(nb_scen, vector<Call>());
    for (int u = 0; u < nb_realizations; ++u) {
      for (size_t i = lb_scen; i < ub_scen; ++i) {
        future_scenarios[u].push_back(this_scenario[i]);
      }
      auto &current_my_scenarios = my_scenarios[3][u];
      // auto& current_my_scenarios = my_scenarios[3][this_block_indexes[3]];
      auto my_end = (my_scenarios[3][u].size() > num_calls_t)
                        ? my_scenarios[3][u].begin() + num_calls_t
                        : my_scenarios[3][u].end();
      // auto my_end = my_scenarios[3][u].end();
      future_scenarios[u].insert(future_scenarios[u].end(),
                                 my_scenarios[3][u].begin(), my_end);
    }
    // future_scenarios[0].insert(future_scenarios[0].end(),
    // my_scenarios[3][block_indexes[0][3]].begin(),
    // 	my_scenarios[3][block_indexes[0][3]].end());
    for (size_t i = 0; i < nb_scen; ++i) {
      for (size_t j = 0; j < future_scenarios[i].size(); ++j) {
        future_scenarios[i][j].id = j;
      }
    }
    future_scenarios.erase(
        future_scenarios.begin() + min(nb_scen, num_future_scenarios),
        future_scenarios.end());
    return future_scenarios;
  } else {
    return vector<vector<Call>>();
  }
}

SecondStageWaitingTime PolicyTester::get_waiting_time(
    int amb_id, Call &call, const string &policy, vector<int> &index_scenarios,
    int nearest_base_id, std::vector<int> queue, int sc, double time) {
  vector<Ambulance> temp_ambulances = ambulances;
  double min_time = 0;
  double waiting_to_hospital = 0;
  if (amb_id >= 0) {
    Ambulance &amb = temp_ambulances[amb_id];
    min_time = travel.get_response_time(amb, call, time, true);
    waiting_to_hospital =
        amb.answer_call(call, travel, ins, time, min_time, nearest_base_id);
  }

  double sum_total = 0;
  int num_total = 0;

  vector<int> busy_ambulances;
  if (amb_id == -1) {
    for (auto &amb : temp_ambulances) {
      if (amb.arrival_time_at_f_last_trip > time) {
        busy_ambulances.push_back(amb.id);
      }
    }

    if (busy_ambulances.size() == 0) {
      return SecondStageWaitingTime{0, 0, GRB_INFINITY};
    }
  }

  for (size_t j = 0; j < index_scenarios.size(); ++j) {
    if (static_cast<size_t>(index_scenarios[j]) < ins.calls[j].size() &&
        ins.calls[j][index_scenarios[j]].time >= call.time) {
      vector<Ambulance> aux_ambulances = temp_ambulances;
      vector<Call> future_calls(ins.calls[j].begin() + index_scenarios[j],
                                ins.calls[j].end());

      for (auto i : queue) {
        if (amb_id == -1 || (i != call.id)) {
          future_calls.insert(future_calls.begin(), ins.calls[sc][i]);
          if (i == call.id) {
            // If the current call must go to queue,
            // all free ambulances are excluded.
            for (auto &amb : aux_ambulances) {
              if (amb.arrival_time_at_f_last_trip <= time) {
                future_calls[0].ambulances.erase(
                    remove(future_calls[0].ambulances.begin(),
                           future_calls[0].ambulances.end(), amb_id),
                    future_calls[0].ambulances.end());
              }
            }
          }
        }
      }

      for (size_t i = 0; i < future_calls.size(); ++i) {
        future_calls[i].id = i;
      }

      auto solver =
          get_solver(policy, future_calls, aux_ambulances, travel, time);
      solver->run();

      // if(amb_id == -1){
      // 	solver->print_results();
      // 	cout << "curr call " << penalized_response_time(min_time, -1,
      // call.priority) << "\n"; 	std::cin.get();
      // }

      if (policy != "cg" && policy != "enumerate") {
        sum_total += accumulate(solver->waiting_on_scene_penalized.begin(),
                                solver->waiting_on_scene_penalized.end(), 0.0);
        // sum_total += accumulate(solver->waiting_to_hospital.begin(),
        // 	solver->waiting_to_hospital.end(), 0.0);
        num_total += solver->waiting_on_scene.size();
      } else {
        sum_total += solver->obj;
        ++num_total;
      }
      // fmt::print("Sum total {}\n", sum_total);
      // std::cin.get();
    }
  }

  if (amb_id != -1) {
    sum_total += penalized_response_time(min_time, -1, call.priority);
    num_total += 1;
  }

  double result_total = (num_total == 0) ? 0 : sum_total / num_total;
  return SecondStageWaitingTime{(time + min_time) - call.time,
                                waiting_to_hospital, result_total};
}

SecondStageWaitingTime PolicyTester::get_waiting_time_return(
    int amb_id, const string &policy, vector<int> &index_scenarios, double time,
    std::vector<int> queue, int sc, int base) {
  vector<Ambulance> temp_ambulances = ambulances;

  auto &amb = temp_ambulances[amb_id];
  amb.base_location = ins.bases[base];
  amb.arrival_time_at_b_last_trip =
      time + travel.travel_time(amb.free_location, amb.base_location, amb);

  double min_time = 0;
  double waiting_to_hospital = 0;
  double sum_total = 0;
  int num_total = 0;
  for (size_t j = 0; j < index_scenarios.size(); ++j) {
    if (static_cast<size_t>(index_scenarios[j]) < ins.calls[j].size() &&
        ins.calls[j][index_scenarios[j]].time >= time) {
      vector<Ambulance> aux_ambulances = temp_ambulances;
      vector<Call> future_calls(ins.calls[j].begin() + index_scenarios[j],
                                ins.calls[j].end());

      for (auto i : queue) {
        future_calls.insert(future_calls.begin(), ins.calls[sc][i]);
        // By returning the current ambulance to base,
        // we forbid it to answer any calls in queue:
        future_calls[0].ambulances.erase(
            remove(future_calls[0].ambulances.begin(),
                   future_calls[0].ambulances.end(), amb_id),
            future_calls[0].ambulances.end());
      }

      for (size_t i = 0; i < future_calls.size(); ++i) {
        future_calls[i].id = i;
      }

      auto solver = get_solver(policy, future_calls, aux_ambulances, travel);
      solver->run();
      // solver->print_results();
      // std::cin.get();

      if (policy != "cg" && policy != "enumerate") {
        sum_total += accumulate(solver->waiting_on_scene_penalized.begin(),
                                solver->waiting_on_scene_penalized.end(), 0.0);
        // sum_total += accumulate(solver->waiting_to_hospital.begin(),
        // 	solver->waiting_to_hospital.end(), 0.0);
        num_total += solver->waiting_on_scene.size();
      } else {
        sum_total += solver->obj;
        ++num_total;
      }
      // fmt::print("Sum total {}\n", sum_total);
      // std::cin.get();
    }
  }

  double result_total = (num_total == 0) ? 0 : sum_total / num_total;
  return SecondStageWaitingTime{(time + min_time) - time, waiting_to_hospital,
                                result_total};
}

vector<int> PolicyTester::get_nearest_base(vector<Call> &calls) {
  vector<int> nearest_base;
  for (size_t i = 0; i < calls.size(); ++i) {
    auto &call = calls[i];
    double min_dist = GRB_INFINITY;
    int min_ind = -1;
    for (int b = 0; b < ins.nb_bases; ++b) {
      double dist = GRB_INFINITY;
      auto &base = ins.bases[b];
      if (call.clean_needed) {
        dist =
            travel.lat_long_distance(ins.cleaning_bases[call.cleaning], base);
      } else if (call.hosp_needed) {
        dist = travel.lat_long_distance(ins.hospitals[call.hospital], base);
      } else {
        dist = travel.lat_long_distance(call.location, base);
      }

      if (dist < min_dist) {
        min_dist = dist;
        min_ind = b;
      }
    }
    nearest_base.push_back(min_ind);
  }

  return nearest_base;
}

int PolicyTester::set_next_event(vector<Call> &this_scenario, int &event_call,
                                 int &index_call, std::vector<int> &queue) {
  vector<pair<double, int>> future_arrival_times;

  bool time_greater_last_call = time > this_scenario.back().time;
  for (size_t i = 0; i < ambulances.size(); ++i) {
    auto &amb = ambulances[i];
    if (amb.arrival_time_at_f_last_trip > time + g_params.EPS) {
      future_arrival_times.push_back(
          make_pair(amb.arrival_time_at_f_last_trip, i));
    }
  }

  int nb_calls = this_scenario.size();
  // fmt::print("Future Arrival times size = {}\n",
  // future_arrival_times.size());
  if (future_arrival_times.size() > 0) {
    auto min_elem =
        *min_element(future_arrival_times.begin(), future_arrival_times.end());

    double min_arrival_time = min_elem.first;
    if (index_call < nb_calls - 1 &&
        this_scenario[index_call + 1].time < min_arrival_time + g_params.EPS) {
      event_call = 1;
      index_call += 1;
      time = this_scenario[index_call].time;
    } else {
      event_call = 0;
      time = min_arrival_time;
      return min_elem.second;
    }
  } else {
    event_call = 1;
    if (time_greater_last_call) {
      index_call = queue[0];
      // fmt::print("WARNING: Reset detected in simulator: queue = {},
      // index_call = {} / {}\n", 	queue, index_call,
      // this_scenario.size());
    } else {
      index_call += 1;
      time = this_scenario[index_call].time;
    }
    return -1;
  }
  // else if(index_call < nb_calls - 1){
  // 	event_call = 1;
  // 	index_call += 1;
  // 	time = this_scenario[index_call].time;
  // }

  return -1;
}

Stats::Stats(vector<double> waiting_on_scene,
             vector<double> waiting_on_scene_penalized,
             vector<double> waiting_to_hospital) {
  assert(waiting_on_scene.size() == waiting_to_hospital.size() &&
         "Stats ERROR: Waiting times vectors must have the same size!");
  if (waiting_on_scene.size() == 0) {
    mean_waiting_on_scene = GRB_INFINITY;
    mean_waiting_on_scene_penalized = GRB_INFINITY;
    mean_waiting_to_hospital = GRB_INFINITY;
    max_waiting_on_scene = GRB_INFINITY;
    max_waiting_on_scene_penalized = GRB_INFINITY;
    max_waiting_to_hospital = GRB_INFINITY;
    std_dev = GRB_INFINITY;
    std_dev_avg = GRB_INFINITY;
    mean_total = GRB_INFINITY;
    max_total = GRB_INFINITY;
    waiting_on_scene_q90 = -1;
    waiting_to_hospital_q90 = -1;
    q90_total = -1;
  } else {
    mean_waiting_on_scene =
        accumulate(waiting_on_scene.begin(), waiting_on_scene.end(), 0.0) /
        waiting_on_scene.size();
    mean_waiting_on_scene_penalized =
        accumulate(waiting_on_scene_penalized.begin(),
                   waiting_on_scene_penalized.end(), 0.0) /
        waiting_on_scene_penalized.size();
    mean_waiting_to_hospital = accumulate(waiting_to_hospital.begin(),
                                          waiting_to_hospital.end(), 0.0) /
                               waiting_to_hospital.size();
    max_waiting_on_scene =
        *max_element(waiting_on_scene.begin(), waiting_on_scene.end());
    max_waiting_on_scene_penalized = *max_element(
        waiting_on_scene_penalized.begin(), waiting_on_scene_penalized.end());
    max_waiting_to_hospital =
        *max_element(waiting_to_hospital.begin(), waiting_to_hospital.end());
    vector<double> total = waiting_on_scene;
    transform(total.begin(), total.end(), waiting_to_hospital.begin(),
              total.begin(), plus<double>());
    mean_total = accumulate(total.begin(), total.end(), 0.0) / total.size();
    max_total = *max_element(total.begin(), total.end());
    sort(waiting_on_scene.begin(), waiting_on_scene.end());
    sort(waiting_to_hospital.begin(), waiting_to_hospital.end());
    sort(total.begin(), total.end());

    waiting_on_scene_q90 =
        waiting_on_scene[(int)floor(0.9 * waiting_on_scene.size())];
    waiting_to_hospital_q90 =
        waiting_to_hospital[(int)floor(0.9 * waiting_on_scene.size())];
    q90_total = total[(int)floor(0.9 * waiting_on_scene.size())];

    double sum_diffs = 0;
    double actual_sum_diffs = 0;
    for (size_t i = 0; i < waiting_on_scene.size(); ++i) {
      sum_diffs += pow(
          waiting_on_scene_penalized[i] - mean_waiting_on_scene_penalized, 2);
      actual_sum_diffs += pow(waiting_on_scene[i] - mean_waiting_on_scene, 2);
    }
    actual_std_dev = sqrt(actual_sum_diffs / waiting_on_scene.size());
    actual_std_dev_avg = actual_std_dev / sqrt(waiting_on_scene.size());
    std_dev = sqrt(sum_diffs / waiting_on_scene_penalized.size());
    std_dev_avg = std_dev / sqrt(waiting_on_scene_penalized.size());
  }
}

bool PolicyTester::can_answer(Ambulance &amb, Call &call) {
  bool amb_allowed = false;
  for (int amb_id : call.ambulances) {
    if (amb_id == amb.id) {
      amb_allowed = true;
      break;
    }
  }
  return amb_allowed &&
         (g_params.amb_setup == "us" || amb.type <= call.priority);
}

SecondStageWaitingTime PolicyTester::get_waiting_time_model(
    int amb_id, int call_id, const string &policy,
    vector<vector<Call>> &future_scenarios, vector<Call> &this_scenario,
    double time, std::vector<int> &queue, int base) {
  Call *call = NULL;
  if (call_id >= 0) {
    call = &this_scenario[call_id];
  }

  bool is_call_amb = amb_id != -1 && call_id != -1;
  bool is_queue = call_id != -1 && amb_id == -1;

  vector<int> busy_ambulances;
  vector<Ambulance> temp_ambulances = ambulances;
  double sum_total = 0;
  int num_total = 0;
  if (amb_id == -1) {
    for (auto &amb : temp_ambulances) {
      if (amb.arrival_time_at_f_last_trip > time) {
        busy_ambulances.push_back(amb.id);
      }
    }

    if (busy_ambulances.size() == 0) {
      return SecondStageWaitingTime{0, 0, GRB_INFINITY};
    }
  }

  for (size_t j = 0; j < future_scenarios.size(); ++j) {
    auto &future_calls = future_scenarios[j];
    vector<Ambulance> aux_ambulances = temp_ambulances;

    if (!queue.empty()) {
      for (int i = queue.size() - 1; i >= 0; --i) {
        if (amb_id == -1 ||
            queue[i] != call_id) {  // if call in queue is not the current call
          future_calls.insert(future_calls.begin(), this_scenario[queue[i]]);
          if (queue[i] == call_id) {
            // If the current call must go to queue,
            // all free ambulances are excluded from attending it.
            // and will only be possible to attend it after arriving at the
            // base.
            for (auto &amb : aux_ambulances) {
              if (amb.arrival_time_at_f_last_trip <= time) {
                future_calls[0].ambulances.erase(
                    remove(future_calls[0].ambulances.begin(),
                           future_calls[0].ambulances.end(), amb.id),
                    future_calls[0].ambulances.end());
              }
            }
          }
        }
      }
    }
    if (!future_calls.empty()) {
      vector<int> this_queue;
      for (size_t i = 0; i < future_calls.size(); ++i) {
        future_calls[i].id = i;
        if (future_calls[i].time < time) {
          this_queue.push_back(i);
        }
      }
      if (is_call_amb) {
        Data data(call, &temp_ambulances[amb_id], ins, travel, this_queue,
                  future_calls, temp_ambulances, time);
        CGCall cg(data, env);

        // fix amb_id to call
        auto &amb = temp_ambulances[amb_id];
        if (amb.arrival_time_at_b_last_trip <= time) {
          int b_amb = -1;
          double min_dist = GRB_INFINITY;
          for (size_t b = 0; b < ins.bases.size(); ++b) {
            double d =
                travel.lat_long_distance(amb.base_location, ins.bases[b]);
            if (d < min_dist) {
              min_dist = d;
              b_amb = b;
            }
          }

          cg.fix_abh(amb.type, b_amb, call->hospital, 1);
        } else if (amb.arrival_time_at_f_last_trip <= time) {
          int b_amb = -1;
          double min_dist = GRB_INFINITY;
          for (size_t b = 0; b < ins.bases.size(); ++b) {
            double d =
                travel.lat_long_distance(amb.base_location, ins.bases[b]);
            if (d < min_dist) {
              min_dist = d;
              b_amb = b;
            }
          }

          int l_amb = -1;
          min_dist = GRB_INFINITY;
          Location current_location = travel.ambulance_position(amb, time);
          for (int l = 0; l < data.num_locals; ++l) {
            double d =
                travel.lat_long_distance(current_location, ins.centers[l]);
            if (d < min_dist) {
              min_dist = d;
              l_amb = l;
            }
          }
          cg.fix_albh(amb.type, l_amb, b_amb, call->hospital, 1);
        }
        auto t0 = std::chrono::high_resolution_clock::now();
        cg.solve();
        auto dt = std::chrono::high_resolution_clock::now();
        double this_run_time =
            std::chrono::duration_cast<chrono::milliseconds>(dt - t0).count();
        solver_run_times.push_back(
            make_tuple(future_calls.size(), 1, this_run_time));
        sum_total += cg.model.get(GRB_DoubleAttr_ObjVal);
        num_total += 1;
      } else if (is_queue) {
        Data data(call, NULL, ins, travel, this_queue, future_calls,
                  temp_ambulances, time);
        CGCall cg(data, env);
        cg.fix_queue();
        auto t0 = std::chrono::high_resolution_clock::now();
        cg.solve();
        auto dt = std::chrono::high_resolution_clock::now();
        double this_run_time =
            std::chrono::duration_cast<chrono::milliseconds>(dt - t0).count();
        solver_run_times.push_back(
            make_tuple(future_calls.size(), 1, this_run_time));
        sum_total += cg.model.get(GRB_DoubleAttr_ObjVal);
        num_total += 1;
      }
    }
  }

  return {0, 0, (num_total == 0) ? 0 : sum_total / num_total};
}

SecondStageWaitingTime PolicyTester::get_waiting_time_model_return(
    int amb_id, int base, const string &policy,
    vector<vector<Call>> &future_scenarios, vector<Call> &this_scenario,
    double time, std::vector<int> &queue) {
  vector<Ambulance> temp_ambulances = ambulances;
  auto &amb = temp_ambulances[amb_id];
  amb.base_location = ins.bases[base];
  double time_amb_to_base =
      travel.travel_time(amb.free_location, amb.base_location, amb);
  amb.arrival_time_at_b_last_trip = time + time_amb_to_base;

  double min_time = 0;
  double waiting_to_hospital = 0;
  double sum_total = 0;
  int num_total = 0;

  for (size_t j = 0; j < future_scenarios.size(); ++j) {
    auto &future_calls = future_scenarios[j];

    vector<Ambulance> aux_ambulances = temp_ambulances;
    for (int i = queue.size() - 1; i >= 0; --i) {
      future_calls.insert(future_calls.begin(), this_scenario[queue[i]]);
      // By returning the current ambulance to base,
      // we forbid it to answer any calls in queue:
      future_calls[0].ambulances.erase(
          remove(future_calls[0].ambulances.begin(),
                 future_calls[0].ambulances.end(), amb_id),
          future_calls[0].ambulances.end());
    }

    if (!future_calls.empty()) {
      vector<int> this_queue;
      for (size_t i = 0; i < future_calls.size(); ++i) {
        future_calls[i].id = i;
        if (future_calls[i].time < time) {
          this_queue.push_back(i);
        }
      }

      Data data(NULL, &amb, ins, travel, this_queue, future_calls,
                temp_ambulances, time);
      // fmt::print("Loaded data\n");
      CGAmbulance cg(data, env);
      // fmt::print("Loaded model\n");
      double min_dist = GRB_INFINITY;
      int min_h = -1;
      for (size_t h = 0; h < ins.hospitals.size(); ++h) {
        double d =
            travel.lat_long_distance(amb.free_location, ins.hospitals[h]);
        if (d < min_dist) {
          min_dist = d;
          min_h = h;
        }
      }
      cg.fix_return(amb.type, min_h, base, 1);
      auto t0 = std::chrono::high_resolution_clock::now();
      cg.solve();
      auto dt = std::chrono::high_resolution_clock::now();
      double this_run_time =
          std::chrono::duration_cast<chrono::milliseconds>(dt - t0).count();
      solver_run_times.push_back(
          make_tuple(future_calls.size(), 2, this_run_time));
      sum_total += cg.model.get(GRB_DoubleAttr_ObjVal);
      num_total += 1;
    }
  }

  return {0, 0, (num_total == 0) ? 0 : sum_total / num_total};
}

double PolicyTester::run_model(vector<Call> &scenario, int amb_id, int call_id,
                               int base_id, bool is_return) {
  if (scenario.size() == 0) {
    fmt::print("Empty scenario\n");
    return 0;
  }
  if ((amb_id == -1 || base_id == -1) && is_return) {
    fmt::print(
        "Error: calling return problem without ambulance (amb_id is {}) "
        "or base (base_id is {})\n",
        amb_id, base_id);
    exit(1);
  } else if (!is_return && call_id == -1) {
    fmt::print("Error: calling Call problem without call (call_id is {})\n",
               call_id);
  }

  bool call_problem = !is_return && amb_id != -1 && call_id != -1;
  bool queue_problem = !is_return && amb_id == -1;

  if (is_return) {
    // find which variable to fix

  } else {
  }

  return 0;
}

PolicyTester::~PolicyTester() {}