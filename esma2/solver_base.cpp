#include <fmt/color.h>

#include <boost/math/distributions/poisson.hpp>

#include "../include/call_model.h"
#include "../include/cg.h"
#include "../include/data.h"
#include "../include/full_model_det.h"
#include "../include/future.h"
#include "../include/instance.h"
#include "../include/solver.h"
#include "../include/travel.h"

Solver::Solver(GRBEnv &env, vector<Call> &calls, vector<Ambulance> &ambulances,
               Instance &ins, Travel &travel, double time)
    : env(env),
      calls(calls),
      ambulances(ambulances),
      travel(travel),
      ins(ins),
      time(time),
      event_call(0),
      index_call(0),
      first_time(calls.front().time),
      begin_time(time),
      last_time(calls.back().time + 7200),
      is_prepared(false),
      waiting_on_scene(calls.size(), GRB_INFINITY),
      waiting_on_scene_penalized(calls.size(), GRB_INFINITY),
      waiting_to_hospital(calls.size(), GRB_INFINITY),
      calls_end(calls.size(), GRB_INFINITY),
      which_ambulance(calls.size(), -1),
      obj(0),
      released_amb(-1) {
  set_calls_nearest_bases();

  int times_size = 3;
  time_ahead = 6 * 3600;
  lambda_base = vector<vector<vector<double>>>(
      ins.nb_bases,
      vector<vector<double>>(times_size * ins.time_horizon.size(),
                             vector<double>(ins.nb_priorities, 0)));
  times = vector<double>(ins.time_horizon.size(), 0);

  vector<set<int>> regions_by_base(ins.bases.size(), set<int>());
  for (int b = 0; b < ins.nb_bases; ++b) {
    for (int c = 0; c < ins.nb_priorities; ++c) {
      for (size_t t_ind = 0; t_ind < lambda_base[b].size(); ++t_ind) {
        int t = (t_ind < ins.time_horizon.size())
                    ? ins.time_horizon[t_ind]
                    : ins.time_horizon.back() + t_ind + 1;
        int g_ind = ins.g0;
        if (t >= ins.nb_times) {
          if (g_ind < ins.nb_days - 1) {
            g_ind += 1;
            t -= ins.nb_times;
          }
        }
        lambda_base[b][t_ind][c] += ins.lambda_bases(t, g_ind, b, c);
      }
    }
  }
  // for (size_t t = 0; t < times.size(); ++t) {
  //   times[t] = ins.time_horizon[t] * 0.5;
  // }
  // double t0 = times.back();
  // for (size_t t = 0; t < ins.time_horizon.size(); ++t) {
  //   times.push_back(t0 + 0.5);
  //   t0 = times.back();
  // }
}

void Solver::prepare() {
  bool debug = debug_mode;

  int first_new_call = 0;
  queue.clear();
  while (static_cast<size_t>(first_new_call) < calls.size() &&
         calls[first_new_call].time <= time) {
    queue.push_back(first_new_call);
    ++first_new_call;
  }

  index_call = first_new_call;
  if (debug) {
    fmt::print("Begin time {}\n", time);
    for (auto &call : calls) {
      cout << call << "\n";
    }
    fmt::print("queue = {}, index_call = {}\n", queue, index_call);
  }
  vector<pair<double, int>> future_arrival_times;
  for (size_t i = 0; i < ambulances.size(); ++i) {
    auto &amb = ambulances[i];
    if (debug) {
      cout << amb << "\n";
    }
    if (amb.arrival_time_at_f_last_trip > time) {
      future_arrival_times.push_back(
          make_pair(amb.arrival_time_at_f_last_trip, amb.id));
    }
  }

  if (future_arrival_times.size() > 0) {
    auto min_arrival =
        min_element(future_arrival_times.begin(), future_arrival_times.end());
    if (debug) {
      fmt::print("min_arrival = {}\n", *min_arrival);
    }
    if (static_cast<size_t>(index_call) >= calls.size() ||
        min_arrival->first <= calls[index_call].time) {
      event_call = 0;
      released_amb = min_arrival->second;
      time = min_arrival->first;
    } else {
      event_call = 1;
      time = calls[index_call].time;
      released_amb = -1;
    }
  } else {
    event_call = 1;
    released_amb = -1;
    if (static_cast<size_t>(index_call) < calls.size()) {
      time = calls[index_call].time;
    }
  }

  if (static_cast<size_t>(index_call) >= calls.size() &&
      future_arrival_times.size() == 0) {
    event_call = 1;
    index_call = queue[0];
  }

  if (debug) {
    fmt::print(
        "PREPARE: time = {}, queue = {}, index_call = {}, event_call = "
        "{}, released = {}, |calls| = {}\n",
        time, queue, index_call, event_call, released_amb, calls.size());
  }
  // cin.get();
}

void Solver::set_debug_mode(bool debug) { debug_mode = debug; }

double Solver::penalized_response_time(double response_time, int amb_type,
                                       int call_type) {
  if (g_params.amb_setup == "us") {
    if (amb_type == -1) {
      return response_time * ins.penalties[call_type];
    }
    return response_time * ins.penalties[call_type] +
           g_params.penalty * ins.penalty_matrix[amb_type][call_type];
  } else if (g_params.amb_setup == "rj") {
    return response_time * ins.penalties[call_type];
  } else {
    fmt::print("Unkown amb setup: {}\n", g_params.amb_setup);
    exit(1);
  }
  return response_time * ins.penalties[call_type];
}

// Find best base to send ambulance amb_id
int Solver::find_best_base(int amb_id, bool debug) {
  int index_base = -1;
  int deficit_global_max = -INT_MAX;
  auto &amb = ambulances[amb_id];
  if (debug) {
    fmt::print("Testing amb {} {}, type {}\n", amb_id, amb.id, amb.type);
  }

  for (int c = amb.type; c < ins.nb_priorities; ++c) {
    vector<tuple<int, int, int>> index_bases;
    for (size_t b = 0; b < ins.bases.size(); ++b) {
      // Travel time from free_location to base b
      double travel_time =
          travel.travel_time(amb.free_location, ins.bases[b], amb);
      double time_arrival = amb.arrival_time_at_f_last_trip + travel_time;

      // Count number of ambulances per type at base b
      int nb_amb = 0;
      vector<int> nb_amb_by_type(ins.nb_types_ambulance, 0);
      for (size_t a = 0; a < ambulances.size(); ++a) {
        auto &amb2 = ambulances[a];
        if (amb2.arrival_time_at_b_last_trip <= time_arrival &&
            (g_params.amb_setup == "us" || amb2.type <= c)) {
          bool equal = amb2.base_location == ins.bases[b];
          double d = travel.lat_long_distance(amb2.base_location, ins.bases[b]);
          if (d < g_params.EPS) {  // d == 0
            ++nb_amb_by_type[amb2.type];
          }
        }
      }

      if (debug) {
        fmt::print(
            "Testing call_type = {}, base {}, time_arrival = {}, "
            "nb_ambs = {}\n",
            c, b, time_arrival, nb_amb_by_type);
        fmt::print("Lambdas:\n");
        for (size_t t = 0; t < times.size(); ++t) {
          fmt::print("\t{} = {}\n", t, lambda_base[b][t][c]);
        }
        fmt::print("times  = {}, time_arrival = {}\n", times,
                   time_arrival / 3600);
      }

      // double arrival_b = time_arrival / 3600;
      // size_t index = 0;
      // // Index of time slot that ambulance reaches the base
      // while (index < times.size() && times[index] <= arrival_b) {
      //   ++index;
      // }

      // double lambda = 0;
      // if (index >= times.size()) {  // if arrival after the end of window,
      // get
      //                               // lambda with the last time slot
      //   lambda += lambda_base[b][times.size() - 1][c] *
      //             (arrival_b - times[times.size() - 1]);
      //   if (debug) {
      //     fmt::print("index pre = {}, lam = {}, p = {}\n", index, lambda, c);
      //   }
      // } else {
      //   if (index == 0) {
      //     fmt::print("Weird: time {} before horizon {}\n", arrival_b,
      //     times[0]); cin.get();
      //   }
      //   // lambda for the time interval from arrival_b to the times[index]
      //   // immediately after arrival_b
      //   lambda += lambda_base[b][index - 1][c] * (times[index] - arrival_b);

      //   // total time is the arrival_time at base b plus the predefined
      //   // time_ahead (currently set to 5400s| 1h30m)
      //   double total_time = (time_arrival + time_ahead) / 3600;
      //   if (debug) {
      //     fmt::print("index pre = {}, lam = {}, p = {}\n", index, lambda, c);
      //   }

      //   // lambda for the time interval from times[index] to the index
      //   // immediately before arrival_b + time_ahead
      //   while (index + 1 < times.size() && times[index + 1] <= total_time) {
      //     lambda +=
      //         (times[index + 1] - times[index]) * lambda_base[b][index][c];
      //     ++index;
      //     if (debug) {
      //       fmt::print("\tindex mid = {}, cum_lam = {}, lam_ind = {}, p =
      //       {}\n",
      //                  index, lambda, lambda_base[b][index][c], c);
      //     }
      //   }
      //   // Lambda for the interval [times[index], total_time - times[index]]
      //   if (index < times.size()) {
      //     lambda += (total_time - times[index]) * lambda_base[b][index][c];
      //     if (debug) {
      //       fmt::print("index post = {}, cum_lam = {}, lam_ind = {}, p =
      //       {}\n",
      //                  index, lambda, lambda_base[b][index][c], c);
      //     }
      //   } else {
      //     // times.back() is the last element of times array
      //     lambda +=
      //         (total_time - times.back()) * lambda_base[b][times.size() -
      //         1][c];
      //   }
      // }

      // int nb_calls = poisson_quantile(0.9, lambda, debug);
      int nb_calls = 0;
      // TODO: Remove after artificial example experiment
      if (g_params.artificial_example) {
        nb_calls =
            ((b == 11 && time < 12 * 3600) || (b == 22 && time >= 12 * 3600))
                ? 4
                : 0;
      }
      nb_amb = 0;
      for (int a = c; a >= 0; --a) {
        nb_amb += nb_amb_by_type[a];
      }
      if (debug) {
        fmt::print("c = {}, total_th = {}, nb_calls = {}, nb_amb = {}\n", c,
                   (time_arrival + time_ahead) / 3600, nb_calls, nb_amb);
      }
      double deficit_c = nb_calls - nb_amb;
      index_bases.push_back(make_tuple(deficit_c, nb_calls, b));
    }
    sort(index_bases.begin(), index_bases.end(),
         std::greater<tuple<int, int, int>>());
    int deficit_max = std::get<0>(index_bases[0]);
    int nb_calls = std::get<1>(index_bases[0]);
    if (deficit_max > 0) {
      double min_time = GRB_INFINITY;
      int min_index = -1;
      size_t k = 0;
      while (k < index_bases.size() &&
             std::get<0>(index_bases[k]) == deficit_max) {
        auto base = ins.bases[std::get<2>(index_bases[k])];
        double travel_time = travel.travel_time(amb.free_location, base, amb);
        if (travel_time < min_time) {
          min_time = travel_time;
          min_index = k;
        }
        ++k;
      }

      return std::get<2>(index_bases[min_index]);
    }
  }
  return -1;
}

double Solver::poisson_quantile(double q, double lambda, bool debug) {
  if (lambda < 0.0001) {
    return 0;
  }
  int k = 0;
  double cum_prob = 0;
  double factorial = 0;
  while (cum_prob < q) {
    if (k == 0) {
      factorial = fac(k);
    } else {
      factorial *= k;
    }
    // if(debug){
    // 	fmt::print("Cum_prob = {}, k = {}, factorial  = {}\n", cum_prob, k,
    // factorial);
    // 	// cin.get();
    // }
    cum_prob += (exp(-lambda) * pow(lambda, k)) / factorial;

    ++k;
  }
  if (debug) {
    fmt::print("Quantile {} for lam {} = {}, cum_prob = {}\n", q, lambda, k - 1,
               cum_prob);
  }
  return k - 1;
}

unsigned long int Solver::fac(int n) {
  unsigned long int prod = 1;
  for (int i = n; i >= 1; --i) {
    prod *= i;
  }

  return prod;
}

int Solver::get_best_index_amb(vector<int> &index_ambs) {
  int index_amb = -1;
  sort(index_ambs.begin(), index_ambs.end(),
       [&](int i, int j) { return ambulances[i].type > ambulances[j].type; });

  for (int i : index_ambs) {
    auto &amb = ambulances[i];
    if (amb.arrival_time_at_f_last_trip <= time) {
      index_amb = i;
      break;
    }
  }

  return index_amb;
}

void DistrictSolver::prepare() {
  bool debug = false;
  int first_new_call = 0;
  while (static_cast<size_t>(first_new_call) < calls.size() &&
         calls[first_new_call].time <= time) {
    queue.push_back(first_new_call);
    ++first_new_call;
  }
  index_call = first_new_call;
  if (debug) {
    fmt::print("Begin time {}\n", time);
    for (auto &call : calls) {
      cout << call << "\n";
    }
    fmt::print("queue = {}, index_call = {}\n", queue, index_call);
  }
  vector<pair<double, int>> future_arrival_times;
  for (size_t i = 0; i < ambulances.size(); ++i) {
    auto &amb = ambulances[i];
    if (debug) {
      cout << amb << "\n";
    }
    if (amb.arrival_time_at_b_last_trip >= time) {
      future_arrival_times.push_back(
          make_pair(amb.arrival_time_at_f_last_trip, amb.id));
    }
  }
  if (future_arrival_times.size() > 0) {
    auto min_arrival =
        min_element(future_arrival_times.begin(), future_arrival_times.end());
    if (debug) {
      fmt::print("min_arrival = {}\n", *min_arrival);
    }
    if (static_cast<size_t>(index_call) >= calls.size() ||
        min_arrival->first <= calls[index_call].time) {
      event_call = 0;
      released_amb = min_arrival->second;
      time = min_arrival->first;
    } else {
      event_call = 1;
      time = calls[index_call].time;
      released_amb = -1;
    }
  } else {
    event_call = 1;
    released_amb = -1;
    if (static_cast<size_t>(index_call) < calls.size()) {
      time = calls[index_call].time;
    }
  }

  if (debug) {
    fmt::print(
        "PREPARE: time = {}, queue = {}, index_call = {}, event_call = "
        "{}, released = {}, |calls| = {}\n",
        time, queue, index_call, event_call, released_amb, calls.size());
  }
  // cin.get();
}

void Solver::set_calls_nearest_bases() {
  for (size_t i = 0; i < calls.size(); ++i) {
    auto &call = calls[i];
    double min_dist = GRB_INFINITY;
    int min_ind = -1;
    for (int b = 0; b < ins.nb_bases; ++b) {
      double dist = GRB_INFINITY;
      auto &base = ins.bases[b];
      if (call.clean_needed) {
        dist = travel.vehicle_distance(ins.cleaning_bases[call.cleaning], base);
      } else if (call.hosp_needed) {
        dist = travel.vehicle_distance(ins.hospitals[call.hospital], base);
      } else {
        dist = travel.vehicle_distance(call.location, base);
      }

      if (dist < min_dist) {
        min_dist = dist;
        min_ind = b;
      }
    }
    nearest_base.push_back(min_ind);
  }
}

void Solver::set_next_event(int &event_call, int &index_call) {
  vector<pair<double, int>> future_arrival_times;
  released_amb = -1;
  bool debug = debug_mode;

  bool no_new_calls = time > calls.back().time;

  double former_time = time;
  double former_event = event_call;
  double former_index = index_call;
  for (size_t i = 0; i < ambulances.size(); ++i) {
    auto &amb = ambulances[i];
    // if(debug){
    // 	cout << amb << "\n";
    // }
    if (amb.arrival_time_at_f_last_trip > time + g_params.EPS) {
      future_arrival_times.push_back(
          make_pair(amb.arrival_time_at_f_last_trip, amb.id));
    }
  }
  int nb_calls = calls.size();
  int first_new_call = 0;
  while (first_new_call < nb_calls &&
         calls[first_new_call].time < time + g_params.EPS) {
    if (abs(calls[first_new_call].time - time) < g_params.EPS &&
        index_call < first_new_call) {  // needed to test if two calls or more
                                        // calls arrive same time
      break;
    }
    ++first_new_call;
  }
  double first_new_call_time =
      (first_new_call < nb_calls) ? calls[first_new_call].time : GRB_INFINITY;
  if (debug) {
    fmt::print(
        "SNE: Setting next event for time {}, index {}, ev {}, queue {}\n",
        time, index_call, event_call, queue);
    fmt::print("SNE: future_arrivals = {}\n", future_arrival_times);
    fmt::print(
        "SNE: min_arrival_time = {}\n",
        min_element(future_arrival_times.begin(), future_arrival_times.end())
            ->first);
    fmt::print("SNE: First new time {}\n", first_new_call_time);
  }
  if (future_arrival_times.size() > 0) {  // there are busy ambulances
    double min_arrival_time;
    tie(min_arrival_time, released_amb) =
        *min_element(future_arrival_times.begin(), future_arrival_times.end());

    // if(index_call >= nb_calls || min_arrival_time <= calls[index_call].time){
    // 	event_call = 0;
    // 	time = min_arrival_time;
    // }else if(former_index < nb_calls && calls[index_call+1].time <
    // min_arrival_time){ 	event_call = 1; 	index_call += 1;
    // time = calls[index_call].time; }else{ 	event_call = 1; 	time =
    // calls[index_call].time;
    // }

    if (min_arrival_time < first_new_call_time + g_params.EPS) {
      event_call = 0;
      time = min_arrival_time;
    } else {
      event_call = 1;
      index_call = first_new_call;
      time = first_new_call_time;
      released_amb = -1;
    }
    // if(index_call < nb_calls && (time < calls[index_call].time || (event_call
    // == 0 & time <= calls[index_call].time)) && 	calls[index_call].time
    // <= min_arrival_time){ 	event_call = 1; 	time =
    // calls[index_call].time; 	released_amb = -1; }else if(index_call <
    // nb_calls-1 && calls[index_call+1].time <= min_arrival_time){ //next call
    // comes before next arrival 	event_call = 1; 	index_call += 1;
    // time = calls[index_call].time; }else{ 	event_call = 0; 	time =
    // min_arrival_time;
    // }
  } else {
    event_call = 1;
    released_amb = -1;
    if (no_new_calls) {
      index_call = queue[0];
    } else {
      index_call = first_new_call;
      time = first_new_call_time;
    }
    // fmt::print("Empty arrivals, index_call = {}, |calls| = {}\n", index_call,
    // calls.size()); if(index_call >= calls.size()){ 	print_results();
    // }
  }
  // else if(index_call < nb_calls - 1){
  // 	event_call = 1;
  // 	index_call += 1;
  // 	time = calls[index_call].time;
  // }
  if (debug) {
    fmt::print("SNE: Next event at time {} set to {}, index_call = {}/{}\n",
               time, event_call, index_call, calls.size());
    // cin.get();
  }
}

bool Solver::can_answer(Ambulance &amb, Call &call) {
  bool amb_allowed = false;
  // fmt::print("call_id = {}, ambulances = {}\n", call.id, call.ambulances);
  for (size_t i = 0; i < call.ambulances.size(); ++i) {
    int amb_id = call.ambulances[i];
    if (amb_id == amb.id) {
      amb_allowed = true;
      break;
    }
  }

  return amb_allowed &&
         (g_params.amb_setup == "us" || amb.type <= call.priority);
}

void Solver::print_results() {
  double average_scene = 0;
  double average_hospital = 0;
  double average_pen = 0;
  cout << "Call\tWait On Scene(s)\tWait Penalized(s)\tAmb:\tType_amb:\n";
  for (size_t i = 0; i < calls.size(); ++i) {
    auto &call = calls[i];
    cout << call << "\t" << waiting_on_scene[i] << "\t";
    cout << waiting_on_scene_penalized[i] << "\t\t" << which_ambulance[i]
         << "\t\t" << ambulances[which_ambulance[i]].type << "\n";

    average_scene += waiting_on_scene[i];
    average_pen += waiting_on_scene_penalized[i];
    average_hospital += waiting_to_hospital[i];
  }
  average_scene /= calls.size();
  average_hospital /= calls.size();
  average_pen /= calls.size();
  // cout << "Average: Real " << average_scene << " Pen " << average_pen <<
  // "\n";
  cout << "\t\tAverage: Pen " << average_pen << "\n";
  // cout << "Total: Real " << average_scene*calls.size() << "  Pen ";
  // cout << "Total: Real " << average_scene*calls.size() << "  Pen ";
  cout << "\t\tTotal pen " << average_pen * calls.size() << "\n";
}

pair<int, double> Solver::get_closest_call(Ambulance &amb) {
  double min_time = GRB_INFINITY;
  int min_call = -1;
  for (auto i : queue) {
    auto &call = calls[i];
    // double response_time = penalized_response_time(
    //     travel.get_response_time(amb, call, time), amb.type, call.priority);
    double response_time = travel.get_response_time(amb, call, time);
    if (can_answer(amb, call) && response_time < min_time) {
      min_time = response_time;
      min_call = i;
    }
  }

  return make_pair(min_call, min_time);
}

pair<int, double> Solver::get_oldest_call(Ambulance &amb) {
  if (queue.empty()) {
    fmt::print("Querying oldest call from empty queue\n");
    return make_pair(-1, GRB_INFINITY);
  }
  auto min_elem =
      min_element(queue.begin(), queue.end(), [&](const int i, const int j) {
        return can_answer(amb, calls[i]) && calls[i].time < calls[j].time;
      });
  int min_call = -1;
  if (min_elem == queue.end()) {
    return make_pair(-1, GRB_INFINITY);
  } else {
    min_call = *min_elem;
  }
  double pen_resp_time = penalized_response_time(
      travel.get_response_time(amb, calls[min_call], time), amb.type,
      calls[min_call].priority);
  return make_pair(min_call, pen_resp_time);
}

Solver::~Solver() {}

NonMiopycSolver::NonMiopycSolver(GRBEnv &env, vector<Call> &calls,
                                 vector<Ambulance> &ambulances, Instance &ins,
                                 Travel &travel, double time)
    : Solver(env, calls, ambulances, ins, travel, time) {
  prepare();
  travel.set_forward(true);
}

void NonMiopycSolver::run() {
  bool debug = debug_mode;

  std::vector<bool> is_allocated(calls.size(), false);
  std::vector<std::vector<double>> travel_times(
      calls.size(), std::vector<double>(ambulances.size(), GRB_INFINITY));
  std::vector<std::vector<double>> pen_travel_times(
      calls.size(), std::vector<double>(ambulances.size(), GRB_INFINITY));
  int max_index = -1;
  std::vector<double> min_times(calls.size(), GRB_INFINITY);
  std::vector<double> min_times_p(calls.size(), GRB_INFINITY);
  std::vector<int> index_amb(calls.size(), -1);
  std::vector<int> amb_type(calls.size(), -1);
  vector<vector<int>> best_ambs(calls.size(), vector<int>());
  travel.set_forward(true);

  if (debug) {
    fmt::print("Queue = {} ,begin time = {}\n", queue, begin_time);
  }
  vector<bool> is_call_on_queue(calls.size(), false);
  for (auto i : queue) {
    is_call_on_queue[i] = true;
  }

  size_t i = 0;

  default_random_engine u_gen;
  uniform_real_distribution<double> u(0, 1);
  // uniform_int_distribution<int> u(0, 100);

  while (i < calls.size()) {
    // travel.set_delay_param(u(u_gen));
    if (debug) {
      fmt::print("Begin iter for call {}\n", i);
    }
    if (!is_allocated[i]) {
      if (debug) {
        int num_attended = 0;
        for (auto val : is_allocated) {
          num_attended += val;
        }
        fmt::print("Index i = {}, max_index = {}, attended = {} / {}\n", i,
                   max_index, num_attended, calls.size());
        std::cout << "Call " << calls[i]
                  << fmt::format(" ambs = {}", calls[i].ambulances) << "\n";
        for (auto &amb : ambulances) {
          cout << amb << "\n";
        }
        cin.get();
      }

      auto t0 = std::chrono::high_resolution_clock::now();
      if (i == static_cast<size_t>(max_index + 1)) {
        if (debug) {
          fmt::print("Call {} is new.\n", i);
          std::cin.get();
        }
        for (size_t j = 0; j < ambulances.size(); ++j) {
          auto &amb = ambulances[j];
          if (can_answer(amb, calls[i])) {
            // travel_times[i][j] = travel.get_response_time(amb, calls[i],
            // calls[i].time) * ins.penalty_matrix[amb.type][calls[i].priority];
            // travel_times[i][j] =
            // penalized_response_time(travel.get_response_time(amb, calls[i],
            // calls[i].time), amb.type, calls[i].priority);
            double time_dispatch = max(begin_time, calls[i].time);
            travel_times[i][j] =
                travel.get_response_time(amb, calls[i], time_dispatch);
            pen_travel_times[i][j] = penalized_response_time(
                travel_times[i][j], amb.type, calls[i].priority);
            if (debug) {
              double total_pen =
                  (calls[i].time < begin_time)
                      ? begin_time - calls[i].time + travel_times[i][j]
                      : travel_times[i][j];
              total_pen = penalized_response_time(total_pen, amb.type,
                                                  calls[i].priority);
              fmt::print(
                  "\tamb {} type {}, call {} p {}, tt = {} pen_tt = {} "
                  "total = {}\n",
                  j, amb.type, i, calls[i].priority, travel_times[i][j],
                  pen_travel_times[i][j], total_pen);
            }
            if (pen_travel_times[i][j] < min_times_p[i]) {
              min_times[i] = travel_times[i][j];
              min_times_p[i] = pen_travel_times[i][j];
              best_ambs[i].clear();
              best_ambs[i].push_back(j);
            } else if (abs(pen_travel_times[i][j] - min_times_p[i]) <
                       g_params.EPS) {
              best_ambs[i].push_back(j);
            }
          }
        }
        max_index = i;
        index_amb[i] = get_best_index_amb(best_ambs[i]);
      }
      if (debug) {
        for (size_t k = i; k <= static_cast<size_t>(max_index); ++k) {
          if (k == i) {
            fmt::print("-> k{}: ", k);
          } else {
            fmt::print("   k{}: ", k);
          }

          for (size_t j = 0; j < ambulances.size(); ++j) {
            double tt = pen_travel_times[k][j];
            if (tt > INT_MAX) {
              fmt::print("inf\t\t");
            } else {
              fmt::print("{:.3f}\t\t", tt);
            }
          }
          fmt::print("\n");
        }
        fmt::print("Min_times_p = {}\n", min_times_p[i]);
        fmt::print("best_ambs = {}\n", best_ambs[i]);
        fmt::print("index_amb = {}, busy = {}\n", index_amb[i],
                   (true)
                       ? get_best_index_amb(best_ambs[i]) == -1
                       : ambulances[get_best_index_amb(best_ambs[i])]
                                 .arrival_time_at_f_last_trip > calls[i].time);
      }
      while (!is_allocated[i]) {
        fmt::print("Beginnning iteration {}. Ambulances:\n", i);
        for (size_t a = 0; a < ambulances.size(); ++a) {
          std::cout << "\t" << ambulances[a] << "\n";
        }

        int index = get_best_index_amb(best_ambs[i]);
        if (index >= 0) {
          index_amb[i] = index;
          auto &amb = ambulances[index_amb[i]];
          is_allocated[i] = true;
          auto &call = calls[i];
          double waiting_on_scene_i =
              (is_call_on_queue[i]) * (begin_time - call.time) + min_times[i];
          double real_min_time = min_times[i];
          // double waiting_to_hospital_i =
          // amb.answer_call(calls[i],travel,ins,calls[i].time,
          // 	travel_times[i][amb.id] /
          // ins.penalty_matrix[amb.type][call.priority], nearest_base[i]);
          double time_dispatch = max(begin_time, calls[i].time);
          int base_id =
              (g_params.best_base) ? find_best_base(amb.id) : nearest_base[i];
          double waiting_to_hospital_i = amb.answer_call(
              calls[i], travel, ins, time_dispatch, real_min_time, base_id);

          if (debug) {
            fmt::print(
                "Amb {} dispatched to call i {}, W_scene = {}, t_free = "
                "{:.1f}, travel_location_hospital = {}\n",
                index_amb[i], i, waiting_on_scene_i,
                amb.arrival_time_at_f_last_trip,
                travel.travel_time(calls[i].location,
                                   ins.hospitals[calls[i].hospital], amb));
          }

          waiting_on_scene[i] = waiting_on_scene_i;
          waiting_on_scene_penalized[i] = penalized_response_time(
              waiting_on_scene_i, amb.type, call.priority);
          waiting_to_hospital[i] = waiting_to_hospital_i;
          which_ambulance[i] = amb.id;
          calls_end[i] = calls[i].end;

          for (int k = i + 1; k <= max_index; ++k) {
            if (!is_allocated[k] && can_answer(amb, calls[k])) {
              // travel_times[k][index_amb[i]] = travel.get_response_time(
              // 	amb,calls[k],calls[k].time) *
              // ins.penalty_matrix[amb.type][calls[k].priority];
              double time_dispatch = max(begin_time, calls[k].time);
              travel_times[k][amb.id] =
                  travel.get_response_time(amb, calls[k], time_dispatch);
              pen_travel_times[k][amb.id] = penalized_response_time(
                  travel_times[k][amb.id], amb.type, calls[k].priority);
              min_times[k] = GRB_INFINITY;
              min_times_p[k] = GRB_INFINITY;
              best_ambs[k].clear();
              for (size_t j = 0; j < ambulances.size(); ++j) {
                auto &amb_j = ambulances[j];
                if (can_answer(amb_j, calls[k])) {
                  if (pen_travel_times[k][j] < min_times_p[k]) {
                    min_times[k] = travel_times[k][j];
                    min_times_p[k] = pen_travel_times[k][j];
                    best_ambs[k].clear();
                    best_ambs[k].push_back(j);
                  } else if (abs(pen_travel_times[k][j] - min_times_p[k]) <
                             g_params.EPS) {
                    best_ambs[k].push_back(j);
                  }
                }
              }
            }
            index_amb[k] = get_best_index_amb(best_ambs[k]);
          }
          if (debug) {
            if (debug) {
              fmt::print("Update post-allocation\n");
            }
            for (size_t k = i; k <= static_cast<size_t>(max_index); ++k) {
              if (k == i) {
                fmt::print("-> k{}: ", k);
              } else {
                fmt::print("   k{}: ", k);
              }
              for (size_t j = 0; j < ambulances.size(); ++j) {
                double tt = pen_travel_times[k][j];
                if (tt > INT_MAX) {
                  fmt::print("inf\t");
                } else {
                  fmt::print("{}\t", tt);
                }
              }
              fmt::print("\n");
            }
            fmt::print("Min_times_p = {}\n", min_times_p);
            fmt::print("best_ambs = {}\n", best_ambs[i]);
            fmt::print(
                "index_amb = {}, busy = {}\n", index_amb[i],
                (true) ? get_best_index_amb(best_ambs[i]) == -1
                       : ambulances[get_best_index_amb(best_ambs[i])]
                                 .arrival_time_at_f_last_trip > calls[i].time);
          }

        } else {
          if (debug) {
            fmt::print(
                "Absolute best ambulance for call {} is not available.\n", i);
          }
          size_t ind = 0;
          while (!is_allocated[i] && ind < best_ambs[i].size()) {
            int which_k = found_ambulance(
                best_ambs, max_index, i, best_ambs[i][ind], is_allocated,
                travel_times, pen_travel_times, min_times, min_times_p,
                is_call_on_queue, debug);
            if (debug) {
              fmt::print(
                  "Best call (which_k) for best amb {} is {}. New "
                  "max_index = {}\n",
                  best_ambs[i][ind], which_k, max_index);
              std::cin.get();
            }
            if (which_k >= 0) {
              is_allocated[which_k] = true;
              auto &amb_i = ambulances[best_ambs[i][ind]];
              auto &call_k = calls[which_k];
              double time_to_free = 0;
              double ttime_to_call = min_times[which_k];
              if (amb_i.arrival_time_at_f_last_trip > calls[which_k].time) {
                time_to_free = amb_i.arrival_time_at_f_last_trip - call_k.time;
                ttime_to_call = travel.travel_time(amb_i.free_location,
                                                   call_k.location, amb_i);
              }
              double waiting_on_scene_i =
                  (is_call_on_queue[which_k]) * (begin_time - call_k.time) +
                  min_times[which_k];
              double real_min_time = waiting_on_scene_i;
              double time_dispatch = max(begin_time, calls[which_k].time);
              int base_id = (g_params.best_base) ? find_best_base(amb_i.id)
                                                 : nearest_base[which_k];
              double waiting_to_hospital_i =
                  amb_i.answer_call(calls[which_k], travel, ins, time_dispatch,
                                    real_min_time, base_id);
              waiting_on_scene[which_k] = waiting_on_scene_i;
              waiting_on_scene_penalized[which_k] = penalized_response_time(
                  waiting_on_scene_i, amb_i.type, call_k.priority);
              waiting_to_hospital[which_k] = waiting_to_hospital_i;
              which_ambulance[which_k] = amb_i.id;
              calls_end[which_k] = calls[which_k].end;

              if (debug) {
                fmt::print(
                    "Amb {} dispatched to call which_k {}, W_scene = "
                    "{}, time_to_free = {}, ttime_to_call = {}, t_free "
                    "= {:.1f}, travel_location_hospital = {}\n",
                    amb_i.id, which_k, waiting_on_scene_i, time_to_free,
                    ttime_to_call, amb_i.arrival_time_at_f_last_trip,
                    travel.travel_time(calls[which_k].location,
                                       ins.hospitals[calls[which_k].hospital],
                                       amb_i));
                fmt::print("Allocated: ");
                for (size_t k = i; k < calls.size(); ++k) {
                  if (is_allocated[k]) {
                    fmt::print("{} ", k);
                  }
                }
                fmt::print("\n");
                std::cin.get();
              }

              for (int k = i; k <= max_index; ++k) {
                if (!is_allocated[k]) {
                  if (can_answer(amb_i, calls[k])) {
                    double time_dispatch = max(begin_time, calls[k].time);
                    travel_times[k][amb_i.id] = travel.get_response_time(
                        amb_i, calls[k], time_dispatch);
                    pen_travel_times[k][amb_i.id] =
                        penalized_response_time(travel_times[k][amb_i.id],
                                                amb_i.type, calls[k].priority);
                  }
                  min_times[k] = GRB_INFINITY;
                  min_times_p[k] = GRB_INFINITY;
                  best_ambs[k].clear();
                  for (size_t j = 0; j < ambulances.size(); ++j) {
                    if (can_answer(ambulances[j], calls[k])) {
                      if (pen_travel_times[k][j] < min_times_p[k]) {
                        min_times[k] = travel_times[k][j];
                        min_times_p[k] = pen_travel_times[k][j];
                        best_ambs[k].clear();
                        best_ambs[k].push_back(j);
                      } else if (abs(pen_travel_times[k][j] - min_times_p[k]) <
                                 g_params.EPS) {
                        best_ambs[k].push_back(j);
                      }
                    }
                  }
                }
              }
              if (debug) {
                fmt::print("Update post-which_k dispatch\n");
                for (size_t k = i; k <= static_cast<size_t>(max_index); ++k) {
                  if (k == i) {
                    fmt::print("-> k{}: ", k);
                  } else {
                    fmt::print("   k{}: ", k);
                  }

                  for (size_t j = 0; j < ambulances.size(); ++j) {
                    double tt = pen_travel_times[k][j];
                    if (tt > INT_MAX) {
                      fmt::print("inf\t");
                    } else {
                      fmt::print("{}\t", tt);
                    }
                  }
                  fmt::print("\n");
                }
                fmt::print("Min_times_p = {}\n", min_times_p[i]);
                fmt::print("best_ambs = {}\n", best_ambs[i]);
                std::cin.get();
              }
              ++ind;
            } else {
              is_allocated[i] = true;
              auto &amb_i = ambulances[best_ambs[i][ind]];
              double time_to_free = 0;
              double ttime_to_call = min_times[i];
              if (amb_i.arrival_time_at_f_last_trip > calls[i].time) {
                time_to_free =
                    amb_i.arrival_time_at_f_last_trip - calls[i].time;
                ttime_to_call = travel.travel_time(amb_i.free_location,
                                                   calls[i].location, amb_i);
              }
              auto &call_i = calls[i];
              double waiting_on_scene_i =
                  (is_call_on_queue[i]) * (begin_time - call_i.time) +
                  min_times[i];
              double real_min_time =
                  travel.get_response_time(amb_i, call_i, call_i.time);
              double time_dispatch = max(begin_time, calls[i].time);
              int base_id = (g_params.best_base) ? find_best_base(amb_i.id)
                                                 : nearest_base[i];
              double waiting_to_hospital_i = amb_i.answer_call(
                  calls[i], travel, ins, time_dispatch, real_min_time, base_id);
              waiting_on_scene[i] = waiting_on_scene_i;
              waiting_on_scene_penalized[i] = penalized_response_time(
                  waiting_on_scene_i, amb_i.type, call_i.priority);
              waiting_to_hospital[i] = waiting_to_hospital_i;
              which_ambulance[i] = amb_i.id;
              calls_end[i] = calls[i].end;

              if (debug) {
                fmt::print(
                    "Amb {} dispatched to call i {}, W_scene = {}, "
                    "time_to_free = {}, ttime_to_call = {}, t_free = "
                    "{:.1f}, travel_location_hospital = {}\n",
                    amb_i.id, i, waiting_on_scene_i, time_to_free,
                    ttime_to_call, amb_i.arrival_time_at_f_last_trip,
                    travel.travel_time(calls[i].location,
                                       ins.hospitals[calls[i].hospital],
                                       amb_i));
                std::cin.get();
              }

              for (int k = i; k <= max_index; ++k) {
                if (!is_allocated[k]) {
                  if (can_answer(amb_i, calls[k])) {
                    double time_dispatch = max(begin_time, calls[k].time);
                    travel_times[k][amb_i.id] = travel.get_response_time(
                        amb_i, calls[k], calls[k].time);
                    pen_travel_times[k][amb_i.id] =
                        penalized_response_time(travel_times[k][amb_i.id],
                                                amb_i.type, calls[k].priority);
                  }
                  min_times[k] = GRB_INFINITY;
                  min_times_p[k] = GRB_INFINITY;
                  best_ambs[k].clear();
                  for (size_t j = 0; j < ambulances.size(); ++j) {
                    if (can_answer(ambulances[j], calls[k])) {
                      if (pen_travel_times[k][j] < min_times_p[k]) {
                        min_times[k] = travel_times[k][j];
                        min_times_p[k] = pen_travel_times[k][j];
                        best_ambs[k].clear();
                        best_ambs[k].push_back(j);
                      } else if (abs(pen_travel_times[k][j] - min_times_p[k]) <
                                 g_params.EPS) {
                        best_ambs[k].push_back(j);
                      }
                    }
                  }
                }
              }

              if (debug) {
                fmt::print("Update post-i dispatch\n");
                for (size_t k = i; k <= static_cast<size_t>(max_index); ++k) {
                  if (k == i) {
                    fmt::print("-> k{}: ", k);
                  } else {
                    fmt::print("   k{}: ", k);
                  }

                  for (size_t j = 0; j < ambulances.size(); ++j) {
                    double tt = pen_travel_times[k][j];
                    if (tt > INT_MAX) {
                      fmt::print("inf\t");
                    } else {
                      fmt::print("{}\t", tt);
                    }
                  }
                  fmt::print("\n");
                }
                fmt::print("Min_times = {}\n", min_times_p[i]);
                fmt::print("best_ambs[{}] = {}\n", i, best_ambs[i]);
                std::cin.get();
              }
            }
          }
        }
      }
      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
    } else {
      ++i;
    }
    if (debug) {
      fmt::print("end_iter, max_index {}\n", max_index);
      // cin.get();
    }
  }

  if (debug) {
    print_results();
    cin.get();
  }
}

int NonMiopycSolver::found_ambulance(
    vector<vector<int>> &best_ambs, int &max_index, int i, int index_amb,
    vector<bool> &is_allocated, vector<vector<double>> &travel_times,
    vector<vector<double>> &pen_travel_times, vector<double> &min_times,
    vector<double> &min_times_p, vector<bool> &is_call_on_queue, bool debug) {
  int which_k = -1;
  size_t k = max_index + 1;
  auto &amb = ambulances[index_amb];

  if (debug) {
    fmt::print("Evaluating future calls for amb {}\n", index_amb);
  }

  while (k < calls.size() && calls[k].time <= amb.arrival_time_at_f_last_trip) {
    if (!is_allocated[k]) {
      for (size_t j = 0; j < ambulances.size(); ++j) {
        if (can_answer(ambulances[j], calls[k])) {
          double time_dispatch = max(begin_time, calls[k].time);
          travel_times[k][j] =
              travel.get_response_time(ambulances[j], calls[k], time_dispatch);
          pen_travel_times[k][j] = penalized_response_time(
              travel_times[k][j], ambulances[j].type, calls[k].priority);
        }
      }

      min_times[k] = GRB_INFINITY;
      min_times_p[k] = GRB_INFINITY;
      best_ambs[k].clear();
      for (size_t j = 0; j < ambulances.size(); ++j) {
        auto &amb_j = ambulances[j];
        if (can_answer(amb_j, calls[k])) {
          if (pen_travel_times[k][j] < min_times_p[k]) {
            min_times[k] = travel_times[k][j];
            min_times_p[k] = pen_travel_times[k][j];
            best_ambs[k].clear();
            best_ambs[k].push_back(j);
          } else if (abs(pen_travel_times[k][j] - min_times_p[k]) <
                     g_params.EPS) {
            best_ambs[k].push_back(j);
          }
        }
      }
    }
    ++k;
  }

  if (k > static_cast<size_t>(max_index) && k <= calls.size()) {
    max_index = k - 1;
  }

  double max_min = -1;
  double w_time = (is_call_on_queue[i]) * (begin_time - calls[i].time) +
                  travel_times[i][index_amb];
  double t1 = penalized_response_time(w_time, amb.type, calls[i].priority);
  k = i + 1;
  if (debug) {
    fmt::print("t1 = {}\n", t1);
  }
  while (k < calls.size() && calls[k].time <= amb.arrival_time_at_f_last_trip) {
    if (!is_allocated[k]) {
      w_time = (is_call_on_queue[k]) * (begin_time - calls[k].time) +
               travel_times[k][index_amb];
      double t2 = penalized_response_time(w_time, amb.type, calls[k].priority);
      if (debug) {
        fmt::print("\tt2 of call {} = {} , best_ambs = {}\n", k, t2,
                   best_ambs[k]);
      }
      if (abs(min_times_p[k] - pen_travel_times[k][index_amb]) < g_params.EPS &&
          t2 > max_min) {
        which_k = k;
        max_min = t2;
      }
    }
    ++k;
  }
  if (max_min < t1 + g_params.EPS) {
    which_k = -1;
  }
  if (debug) {
    fmt::print("max_min = {}, which_k = {}\n", max_min, which_k);
    std::cin.get();
  }
  return which_k;
}

void NonMiopycSolver::compute_min_times(int k, vector<double> &min_times,
                                        vector<int> &index_amb,
                                        vector<int> &amb_type) {
  double min_time = GRB_INFINITY;
  int best_amb = -1;
  int best_type = -1;

  vector<int> best_ambs;

  for (size_t j = 0; j < ambulances.size(); ++j) {
    auto &amb = ambulances[j];
    double time_amb = GRB_INFINITY;
    if (can_answer(amb, calls[k])) {
      if (amb.arrival_time_at_b_last_trip <= calls[k].time) {
        time_amb =
            travel.travel_time(amb.base_location, calls[k].location, amb);
      } else if (amb.arrival_time_at_f_last_trip <= calls[k].time) {
        Location current_location =
            travel.ambulance_position(amb, calls[k].time);
        time_amb = travel.travel_time(current_location, calls[k].location, amb);
      } else {
        double time_free = amb.arrival_time_at_f_last_trip - calls[k].time;
        double time_free_to_call =
            travel.travel_time(amb.free_location, calls[k].location, amb);
        time_amb = time_free + time_free_to_call;
      }
      if (time_amb < min_time) {
        min_time = time_amb;
        best_ambs.clear();
        best_ambs.push_back(j);
      } else if ((abs(time_amb - min_time) < g_params.EPS)) {
        best_ambs.push_back(j);
      }
    }
  }

  min_times[k] = min_time;
  index_amb[k] = get_best_index_amb(best_ambs);
  amb_type[k] = ambulances[index_amb[k]].type;
}

ForwardSolver::ForwardSolver(GRBEnv &env, vector<Call> &calls,
                             vector<Ambulance> &ambulances, Instance &ins,
                             Travel &travel, double time)
    : Solver(env, calls, ambulances, ins, travel, time) {
  prepare();
  travel.set_forward(true);
}

void ForwardSolver::run() {
  bool debug = debug_mode;

  if (debug) {
    fmt::print("Queue = {}, begin_time = {}\n", queue, begin_time);
  }
  vector<bool> is_call_on_queue(calls.size(), false);
  for (auto i : queue) {
    is_call_on_queue[i] = true;
  }

  travel.set_forward(true);
  default_random_engine u_gen;
  uniform_int_distribution<int> u(0, 1);
  for (size_t i = 0; i < calls.size(); ++i) {
    // travel.set_delay_param(u(u_gen));
    auto t0 = std::chrono::high_resolution_clock::now();
    double min_time = GRB_INFINITY;
    double min_time_p = GRB_INFINITY;
    int index_amb = -1;
    auto &call = calls[i];
    time = call.time;
    if (debug) {
      cout << call << fmt::format(" Ambulances: {}", call.ambulances) << "\n";
    }
    for (auto &amb : ambulances) {
      if (can_answer(amb, call)) {
        double time_dispatch = max(begin_time, call.time);
        double time_amb = travel.get_response_time(amb, call, time_dispatch);
        double pen_time_amb =
            penalized_response_time(time_amb, amb.type, call.priority);

        // double time_amb =
        // penalized_response_time(travel.get_response_time(amb, call, time),
        // amb.type, call.priority);
        if (debug) {
          fmt::print(
              "\tAmb {}, type {}. Time_dispatch = {},  Travel time = "
              "{}, Pen Travel time = {}\n",
              amb.id, amb.type, time_dispatch, time_amb, pen_time_amb);
        }

        if (pen_time_amb < min_time_p) {
          min_time = time_amb;
          min_time_p = pen_time_amb;
          index_amb = amb.id;
        } else if (abs(pen_time_amb - min_time_p) < g_params.EPS &&
                   index_amb != -1 && amb.type > ambulances[index_amb].type) {
          min_time = time_amb;
          min_time_p = pen_time_amb;
          index_amb = amb.id;
        }
      }
    }

    if (index_amb >= 0) {
      auto &amb = ambulances[index_amb];
      double real_min_time = travel.get_response_time(amb, call, time);
      double waiting_on_scene_i =
          (is_call_on_queue[i]) * (begin_time - call.time) + min_time;
      double time_dispatch = max(begin_time, calls[i].time);
      int base_id =
          (g_params.best_base) ? find_best_base(amb.id) : nearest_base[i];
      double waiting_to_hospital_i = amb.answer_call(
          call, travel, ins, time_dispatch, real_min_time, base_id);
      waiting_on_scene[i] = waiting_on_scene_i;
      waiting_on_scene_penalized[i] =
          penalized_response_time(waiting_on_scene_i, amb.type, call.priority);
      waiting_to_hospital[i] = waiting_to_hospital_i;
      which_ambulance[i] = index_amb;
      calls_end[i] = call.end;
      if (debug) {
        fmt::print("Amb {} dispatched to call {}, w_pen = {}\n", index_amb,
                   call.id, waiting_on_scene_penalized[i]);
      }
    } else {
      waiting_on_scene[i] = GRB_INFINITY;
      waiting_on_scene_penalized[i] = GRB_INFINITY;
      waiting_to_hospital[i] = GRB_INFINITY;
      calls_end[i] = GRB_INFINITY;
      which_ambulance[i] = -1;
    }
    auto dt = std::chrono::high_resolution_clock::now();
    run_times.push_back(
        std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
  }

  if (debug) {
    print_results();
    cin.get();
  }
}

QueueSolver::QueueSolver(GRBEnv &env, vector<Call> &calls,
                         vector<Ambulance> &ambulances, Instance &ins,
                         Travel &travel, double time)
    : Solver(env, calls, ambulances, ins, travel, time) {
  prepare();
  travel.set_forward(false);
}

void QueueSolver::run() {
  calls_attended = 0;
  travel.set_forward(false);
  bool debug = debug_mode;
  using fmt::print;
  default_random_engine u_gen;
  uniform_real_distribution<double> u(0, 1);
  // uniform_int_distribution<int> u(0, 100);
  if (debug) {
    fmt::print("Begin_time = {}, time = {}, queue = {}\n", begin_time, time,
               queue);
    fmt::print("Calls:\n");
    for (auto &call : calls) {
      cout << call << "\n";
    }
    fmt::print("Ambulances:\n");
    for (auto &amb : ambulances) {
      cout << amb << "\n";
    }
  }

  while (static_cast<size_t>(calls_attended) < calls.size()) {
    travel.set_delay_param(u(u_gen));
    // for(auto& amb: ambulances){
    // 	double arrival_b = amb.arrival_time_at_b_last_trip;
    // 	if(arrival_b <= time){
    // 		for(size_t c = 0; c < calls.size(); ++c){
    // 			calls[c].ambulances.push_back(amb.id);
    // 		}
    // 	}
    // }
    if (event_call == 1 && static_cast<size_t>(index_call) < calls.size()) {
      queue.push_back(index_call);
    }
    if (debug) {
      fmt::print(
          "QueueSolver: Time {}, index = {}, event = {}, queue = {}, "
          "att = {}/{}\n",
          time, index_call, event_call, queue, calls_attended, calls.size());
    }
    std::vector<int> queue_aux;
    for (size_t i = 0; i < queue.size(); ++i) {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto &call = calls[queue[i]];
      double min_time_p = GRB_INFINITY;
      int index_amb = -1;
      for (auto &amb : ambulances) {
        if (can_answer(amb, call)) {
          // double time_amb =
          // penalized_response_time(travel.get_response_time(amb,call,time),
          // amb.type, call.priority);
          double time_amb = travel.get_response_time(amb, call, time);
          double pen_time_amb =
              penalized_response_time(time_amb, amb.type, call.priority);
          // print("\tAmb {} {} {}\n", amb.id, time_amb, amb.type);
          if (debug) {
            fmt::print("\tAmb {} time {} pen_time {} {} {} | {} {}\n", amb.id,
                       time_amb, pen_time_amb, amb.base_location.first,
                       amb.base_location.second, call.location.first,
                       call.location.second);
          }
          if (pen_time_amb < min_time_p) {
            min_time_p = pen_time_amb;
            index_amb = amb.id;
          }
        }
      }
      if (index_amb >= 0) {
        auto &amb = ambulances[index_amb];
        if (debug) {
          fmt::print("Call {} answered by {}\n", call.id, amb.id);
        }
        double real_min_time = travel.get_response_time(amb, call, time);
        double waiting_on_scene_i = time + real_min_time - call.time;
        int base_id = (g_params.best_base) ? find_best_base(amb.id)
                                           : nearest_base[queue[i]];
        double waiting_to_hospital_i =
            amb.answer_call(call, travel, ins, time, real_min_time, base_id);
        waiting_on_scene[queue[i]] = waiting_on_scene_i;
        waiting_on_scene_penalized[queue[i]] = penalized_response_time(
            waiting_on_scene_i, amb.type, call.priority);
        waiting_to_hospital[queue[i]] = waiting_to_hospital_i;
        which_ambulance[queue[i]] = index_amb;
        calls_end[queue[i]] = call.end;
        calls_attended++;
      } else {
        queue_aux.push_back(queue[i]);
        if (debug) {
          fmt::print("Call {} back to queue\n", call.id);
        }
      }
      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
    }
    queue.clear();
    for (auto i : queue_aux) {
      queue.push_back(i);
    }
    set_next_event(event_call, index_call);
  }
  if (debug) {
    print_results();
    // cin.get();
  }
}

DummyQueueSolver::DummyQueueSolver(GRBEnv &env, vector<Call> &calls,
                                   vector<Ambulance> &ambulances, Instance &ins,
                                   Travel &travel, double time)
    : Solver(env, calls, ambulances, ins, travel, time) {
  prepare();
  travel.set_forward(false);
}

void DummyQueueSolver::run() {
  calls_attended = 0;
  travel.set_forward(false);
  bool debug = false;
  using fmt::print;
  default_random_engine u_gen;
  uniform_real_distribution<double> u(0, 1);
  // uniform_int_distribution<int> u(0, 100);
  if (debug) {
    fmt::print("Begin_time = {}, time = {}, queue = {}\n", begin_time, time,
               queue);
    fmt::print("Calls:\n");
    for (auto &call : calls) {
      cout << call << "\n";
    }
    fmt::print("Ambulances:\n");
    for (auto &amb : ambulances) {
      cout << amb << "\n";
    }
  }

  while (static_cast<size_t>(calls_attended) < calls.size()) {
    travel.set_delay_param(u(u_gen));
    if (event_call == 1 && static_cast<size_t>(index_call) < calls.size()) {
      queue.push_back(index_call);
    }
    if (debug) {
      fmt::print(
          "DummyQueueSolver: Time {}, index = {}, event = {}, queue = {}, "
          "att = {}/{}\n",
          time, index_call, event_call, queue, calls_attended, calls.size());
    }
    std::vector<int> queue_aux;
    for (size_t i = 0; i < queue.size(); ++i) {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto &call = calls[queue[i]];
      double min_time = GRB_INFINITY;
      int index_amb = -1;
      for (auto &amb : ambulances) {
        if (can_answer(amb, call)) {
          // double time_amb =
          // penalized_response_time(travel.get_response_time(amb,call,time),
          // amb.type, call.priority);
          double time_amb = travel.get_response_time(amb, call, time);
          double pen_time_amb =
              penalized_response_time(time_amb, amb.type, call.priority);
          // print("\tAmb {} {} {}\n", amb.id, time_amb, amb.type);
          if (time_amb < min_time) {
            min_time = time_amb;
            index_amb = amb.id;
          }
        }
      }

      if (index_amb >= 0) {
        auto &amb = ambulances[index_amb];
        if (debug) {
          fmt::print("Call {} answered by {}\n", call.id, amb.id);
        }
        double real_min_time = travel.get_response_time(amb, call, time);
        double waiting_on_scene_i = time + real_min_time - call.time;
        int base_id = (g_params.best_base) ? find_best_base(amb.id)
                                           : nearest_base[queue[i]];
        double waiting_to_hospital_i =
            amb.answer_call(call, travel, ins, time, real_min_time, base_id);
        waiting_on_scene[queue[i]] = waiting_on_scene_i;
        waiting_on_scene_penalized[queue[i]] = penalized_response_time(
            waiting_on_scene_i, amb.type, call.priority);
        waiting_to_hospital[queue[i]] = waiting_to_hospital_i;
        which_ambulance[queue[i]] = index_amb;
        calls_end[queue[i]] = call.end;
        calls_attended++;
      } else {
        queue_aux.push_back(queue[i]);
        if (debug) {
          fmt::print("Call {} back to queue\n", call.id);
        }
      }
      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
    }
    queue.clear();
    for (auto i : queue_aux) {
      queue.push_back(i);
    }
    set_next_event(event_call, index_call);
  }
  if (debug) {
    print_results();
    // cin.get();
  }
}

PrioritySolver::PrioritySolver(GRBEnv &env, vector<Call> &calls,
                               vector<Ambulance> &ambulances, Instance &ins,
                               Travel &travel, double time)
    : Solver(env, calls, ambulances, ins, travel, time) {
  // travel.set_forward(false);
  queue.reserve(calls.size());
  prepare();
}

void PrioritySolver::run() {
  int calls_attended = 0;
  bool debug = debug_mode;

  travel.set_forward(true);
  default_random_engine u_gen;
  uniform_real_distribution<double> u(0, 1);
  // uniform_int_distribution<int> u(0, 100);

  if (debug) {
    fmt::print("Begin time = {}, queue = {}\n", begin_time, queue);
    for (auto &call : calls) {
      cout << call << fmt::format(" | Ambs = {}\n", call.ambulances);
      // if(call.ambulances.size() == 0){
      // 	cout << "| Unreachable!!!\n";
      // 	cin.get();
      // }else{
      // 	cout << "\n";
      // }
    }

    for (auto &amb : ambulances) {
      cout << amb << "\n";
    }
  }

  while (static_cast<size_t>(calls_attended) < calls.size()) {
    // travel.set_delay_param(u(u_gen));
    if (event_call == 1 && static_cast<size_t>(index_call) < calls.size()) {
      queue.push_back(index_call);
    }

    std::vector<std::pair<double, int>> sorted_queue;
    sorted_queue.reserve(queue.size());
    for (size_t i = 0; i < queue.size(); ++i) {
      auto &call = calls[queue[i]];
      double pen_resp_time =
          penalized_response_time((time - call.time), -1, call.priority);
      sorted_queue.push_back(std::make_pair(pen_resp_time, queue[i]));
    }
    std::sort(sorted_queue.begin(), sorted_queue.end(),
              std::greater<std::pair<double, int>>());

    if (debug) {
      fmt::print(
          "Time {:.1f}, index_call = {}, event_call = {}, attended = "
          "{}/{}, sorted_queue = {}\n",
          time, index_call, event_call, calls_attended, calls.size(),
          sorted_queue);

      set<int> set_queue(queue.begin(), queue.end());
      if (set_queue.size() != queue.size()) {
        fmt::print("{}\n", fmt::format(fg(fmt::color(0xFF0000)),
                                       "Queue has duplicates!!"));
        cin.get();
      }
    }

    std::vector<int> queue_aux;
    for (size_t k = 0; k < sorted_queue.size(); ++k) {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto call_ind = sorted_queue[k].second;
      auto &call = calls[call_ind];
      double min_time = GRB_INFINITY;
      double min_time_p = GRB_INFINITY;
      int index_amb = -1;
      if (debug) {
        cout << "Call " << call << "\n";
      }
      vector<int> index_ambs;

      for (auto &amb : ambulances) {
        if (can_answer(amb, call)) {
          // double time_amb =
          // penalized_response_time(travel.get_response_time(amb,call,time));
          double time_amb = travel.get_response_time(amb, call, time);
          double pen_time_amb =
              penalized_response_time(time_amb, amb.type, call.priority);
          double alloc_cost = penalized_response_time(
              time + time_amb - call.time, amb.type, call.priority);
          if (debug) {
            fmt::print(
                "\tAmb {} t_b = {} t_f = {} -> response_time = {} alloc_cost = "
                "{}\n",
                amb.id, amb.arrival_time_at_b_last_trip,
                amb.arrival_time_at_f_last_trip, time_amb, alloc_cost);
          }

          if (pen_time_amb < min_time_p) {
            min_time = time_amb;
            min_time_p = pen_time_amb;
            index_ambs.clear();
            index_ambs.push_back(amb.id);
          } else if (abs(pen_time_amb - min_time_p) < g_params.EPS) {
            index_ambs.push_back(amb.id);
          }
        }
      }

      index_amb = get_best_index_amb(index_ambs);

      if (debug) {
        fmt::print("index_amb = {}, index_ambs = {}\n", index_amb, index_ambs);
      }

      if (index_amb >= 0) {
        auto &amb = ambulances[index_amb];
        double waiting_on_scene_i = time + min_time - call.time;
        int base_id = (g_params.best_base) ? find_best_base(amb.id)
                                           : nearest_base[call_ind];
        double waiting_to_hospital_i =
            amb.answer_call(call, travel, ins, time, min_time, base_id);
        waiting_on_scene[call_ind] = waiting_on_scene_i;
        waiting_on_scene_penalized[call_ind] = penalized_response_time(
            waiting_on_scene_i, amb.type, call.priority);
        waiting_to_hospital[call_ind] = waiting_to_hospital_i;
        which_ambulance[call_ind] = index_amb;
        calls_end[call_ind] = call.end;
        calls_attended++;

        if (debug) {
          string msg_dispatch = fmt::format(
              fg(fmt::color(0xFFFF00)),
              "Call {} answered by {}. w_scene = {}, free = {}\n", call_ind,
              index_amb, waiting_on_scene_penalized[call_ind],
              amb.arrival_time_at_f_last_trip);
          fmt::print("{}", msg_dispatch);
        }

      } else {
        queue_aux.push_back(sorted_queue[k].second);
      }
      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
    }
    queue.clear();
    for (auto i : queue_aux) {
      queue.push_back(i);
    }
    set_next_event(event_call, index_call);

    if (debug) {
      fmt::print("Next event = {}, index = {}, queue = {}\n", event_call,
                 index_call, queue);
      cin.get();
    }
  }
}

MinMaxPSolver::MinMaxPSolver(GRBEnv &env, vector<Call> &calls,
                             vector<Ambulance> &ambulances, Instance &ins,
                             Travel &travel, double time)
    : Solver(env, calls, ambulances, ins, travel, time) {
  // travel.set_forward(false);
  queue.reserve(calls.size());
}

void MinMaxPSolver::run() {
  int event_call = 1;
  time = calls[0].time;
  int index_call = 0;
  calls_attended = 0;
  // travel.set_forward(false);
  while (static_cast<size_t>(calls_attended) < calls.size()) {
    if (event_call == 1) {
      queue.push_back(index_call);
    }

    std::vector<double> min_times(queue.size(), GRB_INFINITY);
    std::vector<double> min_times_p(queue.size(), GRB_INFINITY);
    std::vector<std::vector<double>> travel_times(
        queue.size(), std::vector<double>(ambulances.size(), GRB_INFINITY));
    std::vector<int> index_ambs(queue.size(), -1);

    // For every call in queue, determine the best ambulance and the best
    // response time
    for (size_t i = 0; i < queue.size(); ++i) {
      auto &call = calls[queue[i]];
      for (auto &amb : ambulances) {
        int j = amb.id;
        if (can_answer(amb, call)) {
          travel_times[i][j] = travel.get_response_time(amb, call, time) *
                               ins.penalty_matrix[amb.type][call.priority];

          if (travel_times[i][j] < min_times[i]) {
            min_times[i] = travel_times[i][j];
            index_ambs[i] = j;
          } else if (abs(travel_times[i][j] - min_times[i]) < g_params.EPS &&
                     index_ambs[i] != -1 &&
                     amb.type > ambulances[index_ambs[i]].type) {
            min_times[i] = travel_times[i][j];
            index_ambs[i] = j;
          }
        }
      }
      min_times[i] += time - call.time;
      // FIXME: EMS setup depends on ambulance and call. How to adjust in this
      // case?
      min_times_p[i] = penalized_response_time(min_times[i], -1, call.priority);
    }

    std::vector<int> remaining_indexes(queue.size(), 0);
    for (size_t i = 0; i < queue.size(); ++i) {
      remaining_indexes[i] = i;
    }

    std::vector<int> queue_aux;
    int nb_treated = 0;
    int total_queue = queue.size();

    while (nb_treated < total_queue) {
      auto t0 = std::chrono::high_resolution_clock::now();
      // get time and index of the worst call, and the best ambulance to such
      // call
      auto max_it = std::max_element(min_times_p.begin(), min_times_p.end());
      int current_call = std::distance(min_times_p.begin(), max_it);
      double max_time = min_times[current_call];
      int call_ind = queue[remaining_indexes[current_call]];
      auto &call = calls[call_ind];
      int index_amb = index_ambs[remaining_indexes[current_call]];
      auto &best_amb = ambulances[index_amb];

      // if best ambulance is free, it should answer the call
      if (index_amb >= 0 && best_amb.arrival_time_at_f_last_trip <= time) {
        double waiting_on_scene_i = max_time;
        double waiting_to_hospital_i = ambulances[index_amb].answer_call(
            call, travel, ins, time, max_time - (time - call.time),
            nearest_base[call_ind]);
        waiting_on_scene[call_ind] = waiting_on_scene_i;
        waiting_on_scene_penalized[call_ind] = penalized_response_time(
            waiting_on_scene_i, ambulances[index_amb].type,
            calls[call_ind].priority);
        waiting_to_hospital[call_ind] = waiting_to_hospital_i;
        which_ambulance[call_ind] = index_amb;
        calls_end[call_ind] = call.end;

        calls_attended++;
        // update travel times
        for (int k = 0; k < total_queue - nb_treated; ++k) {
          int ind = remaining_indexes[k];
          if (k != current_call && index_amb == index_ambs[ind]) {
            double time_to_h = best_amb.arrival_time_at_f_last_trip - time;
            double time_from_h_to_c = travel.travel_time(
                best_amb.free_location, calls[queue[ind]].location, best_amb);
            travel_times[ind][best_amb.id] = time_to_h + time_from_h_to_c;
            auto min_it = std::min_element(travel_times[ind].begin(),
                                           travel_times[ind].end());
            double min_val = *min_it;
            int min_ind = std::distance(travel_times[ind].begin(), min_it);
            min_times[k] = min_val + time - calls[queue[ind]].time;
            // FIXME: EMS System
            min_times_p[k] = penalized_response_time(
                min_times[k], -1, calls[queue[ind]].priority);
            index_ambs[ind] = min_ind;
          }
        }
      } else {
        queue_aux.push_back(queue[remaining_indexes[current_call]]);
      }
      if (nb_treated + 1 < total_queue) {
        min_times_p.erase(max_it);
        remaining_indexes.erase(
            std::remove(remaining_indexes.begin(), remaining_indexes.end(),
                        remaining_indexes[current_call]),
            remaining_indexes.end());
        min_times.erase(std::remove(min_times.begin(), min_times.end(),
                                    min_times[current_call]),
                        min_times.end());
      }
      nb_treated++;
      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count() *
          queue.size());
    }

    queue.clear();
    for (auto i : queue_aux) {
      queue.push_back(i);
    }

    set_next_event(event_call, index_call);
  }
}

GenForwardSolver::GenForwardSolver(GRBEnv &env, vector<Call> &calls,
                                   vector<Ambulance> &ambulances, Instance &ins,
                                   Travel &travel, double time)
    : Solver(env, calls, ambulances, ins, travel, time) {
  g_params.extended_model = false;
  travel.set_forward(true);
  queue.reserve(calls.size());
}

void GenForwardSolver::run() {
  int index_call = 0;
  queue.clear();
  while (static_cast<size_t>(index_call) < calls.size()) {
    fmt::print("index_call {}\n", index_call);
    time = calls[index_call].time;
    queue.push_back(index_call);
    int current = index_call + 1;
    while (static_cast<size_t>(current) < calls.size() &&
           calls[current].time - time < g_params.EPS) {
      queue.push_back(current);
      current += 1;
    }
    fmt::print("queue size {}\n", queue.size());
    if (queue.size() > 1) {
      // call model
      // call answer_call to each call in queue
      // set index_call += |answered calls|
      try {
        CallModel cm(env, *this);
        cm.solve();
        for (size_t i = 0; i < queue.size(); ++i) {
          auto &best_amb = ambulances[cm.opt_which_amb[i]];
          double min_time = cm.arrival_times[i] - time;
          double waiting_to_hospital_i =
              best_amb.answer_call(calls[queue[i]], travel, ins, time, min_time,
                                   nearest_base[queue[i]]);
          auto &call = calls[queue[i]];
          waiting_on_scene[queue[i]] =
              cm.arrival_times[i] - calls[queue[i]].time;
          waiting_on_scene_penalized[queue[i]] = penalized_response_time(
              waiting_on_scene[queue[i]], best_amb.type, call.priority);
          waiting_to_hospital[queue[i]] = waiting_to_hospital_i;
          which_ambulance[queue[i]] = cm.opt_which_amb[i];
          calls_end[queue[i]] = calls[queue[i]].end;
        }
        index_call = current;
        // print_results();
      } catch (GRBException &ex) {
        std::cout << ex.getMessage() << "\n";
        exit(1);
      }
    } else {
      double min_time = GRB_INFINITY;
      int index_amb = -1;
      auto &call = calls[queue[0]];
      time = call.time;
      for (auto &amb : ambulances) {
        if (can_answer(amb, call)) {
          double time_amb = travel.get_response_time(amb, call, time) *
                            ins.penalty_matrix[amb.type][call.priority];
          if (time_amb < min_time) {
            min_time = time_amb;
            index_amb = amb.id;
          }
        }
      }
      if (index_amb >= 0) {
        auto &amb = ambulances[index_amb];
        double waiting_to_hospital_i = amb.answer_call(
            call, travel, ins, time,
            min_time / ins.penalty_matrix[amb.type][call.priority],
            nearest_base[queue[0]]);
        waiting_on_scene[queue[0]] = time - call.time + min_time;
        waiting_on_scene_penalized[queue[0]] = penalized_response_time(
            waiting_on_scene[queue[0]], amb.type, calls[queue[0]].priority);
        waiting_to_hospital[queue[0]] = waiting_to_hospital_i;
        which_ambulance[queue[0]] = index_amb;
        calls_end[queue[0]] = calls[queue[0]].end;
        index_call += 1;
      } else {
        waiting_on_scene[queue[0]] = GRB_INFINITY;
        waiting_on_scene_penalized[queue[0]] = GRB_INFINITY;
        waiting_to_hospital[queue[0]] = GRB_INFINITY;
        which_ambulance[queue[0]] = -1;
        calls_end[queue[0]] = GRB_INFINITY;
      }
    }
    queue.clear();
    fmt::print("==================\n");
  }
}

GenMinMaxSolver::GenMinMaxSolver(GRBEnv &env, vector<Call> &calls,
                                 vector<Ambulance> &ambulances, Instance &ins,
                                 Travel &travel, double time)
    : Solver(env, calls, ambulances, ins, travel, time) {
  g_params.extended_model = false;
  travel.set_forward(false);
  queue.reserve(calls.size());
}

void GenMinMaxSolver::run() {
  int event_call = 1;
  time = calls[0].time;
  int index_call = 0;
  calls_attended = 0;

  default_random_engine u_gen;
  uniform_real_distribution<double> u(0, 1);
  // uniform_int_distribution<int> u(0, 100);
  travel.set_delay_param(u(u_gen));

  while (static_cast<size_t>(calls_attended) < calls.size()) {
    if (event_call == 1) {
      queue.push_back(index_call);
    }

    int current = -1;
    if (event_call == 1) {
      current = index_call + 1;
      while (current < static_cast<int>(calls.size()) &&
             calls[current].time - time < g_params.EPS) {
        queue.push_back(current);
        current += 1;
      }
    }

    if (queue.size() > 1) {
      try {
        CallModel cm(env, *this);
        cm.solve();
        for (size_t i = 0; i < queue.size(); ++i) {
          auto &best_amb = ambulances[cm.opt_which_amb[i]];
          double min_time = cm.arrival_times[i] - time;
          double waiting_to_hospital_i =
              best_amb.answer_call(calls[queue[i]], travel, ins, time, min_time,
                                   nearest_base[queue[i]]);
          waiting_on_scene[queue[i]] =
              cm.arrival_times[i] - calls[queue[i]].time;
          auto &call = calls[queue[i]];
          waiting_on_scene_penalized[queue[i]] =
              penalized_response_time(waiting_on_scene[queue[i]], best_amb.type,
                                      calls[queue[i]].priority);
          waiting_to_hospital[queue[i]] = waiting_to_hospital_i;
          which_ambulance[queue[i]] = cm.opt_which_amb[i];
          calls_end[queue[i]] = calls[queue[i]].end;
          calls_attended++;
        }
        // cm.print_results();
        if (event_call == 1) index_call = current - 1;
        queue.clear();
        set_next_event(event_call, index_call);
        continue;
      } catch (GRBException &ex) {
        std::cout << ex.getMessage() << "\n";
        exit(1);
      }
    }

    std::vector<double> min_times(queue.size(), GRB_INFINITY);
    std::vector<std::vector<double>> travel_times(
        queue.size(), std::vector<double>(ambulances.size(), GRB_INFINITY));
    std::vector<int> index_ambs(queue.size(), -1);

    for (size_t i = 0; i < queue.size(); ++i) {
      auto &call = calls[queue[i]];
      for (auto &amb : ambulances) {
        int j = amb.id;
        if (can_answer(amb, call)) {
          travel_times[i][j] = travel.get_response_time(amb, call, time) *
                               ins.penalty_matrix[amb.type][call.priority];

          if (travel_times[i][j] < min_times[i]) {
            min_times[i] = travel_times[i][j];
            index_ambs[i] = j;
          }
        }
      }
      min_times[i] += time - call.time;
    }

    std::vector<int> remaining_indexes(queue.size(), 0);
    for (size_t i = 0; i < queue.size(); ++i) {
      remaining_indexes[i] = i;
    }

    std::vector<int> queue_aux;
    int nb_treated = 0;
    int total_queue = queue.size();

    while (nb_treated < total_queue) {
      auto max_it = std::max_element(min_times.begin(), min_times.end());
      double max_time = *max_it;
      int current_call = std::distance(min_times.begin(), max_it);
      auto &call = calls[queue[remaining_indexes[current_call]]];
      int index_amb = index_ambs[remaining_indexes[current_call]];
      auto &best_amb = ambulances[index_amb];

      if (index_amb >= 0 && best_amb.arrival_time_at_f_last_trip <= time) {
        int call_ind = queue[remaining_indexes[current_call]];
        double min_time = max_time - (time - call.time);
        double waiting_to_hospital_i = best_amb.answer_call(
            call, travel, ins, time, min_time, nearest_base[call_ind]);

        waiting_on_scene[call_ind] = max_time;
        waiting_on_scene_penalized[call_ind] =
            penalized_response_time(waiting_on_scene[call_ind], best_amb.type,
                                    calls[call_ind].priority);
        waiting_to_hospital[call_ind] = waiting_to_hospital_i;
        which_ambulance[call_ind] = index_amb;
        calls_end[call_ind] = call.end;

        calls_attended++;
        // update travel times
        for (int k = 0; k < total_queue - nb_treated; ++k) {
          int ind = remaining_indexes[k];
          if (k != current_call && index_amb == index_ambs[ind]) {
            double time_to_h = best_amb.arrival_time_at_f_last_trip - time;
            double time_from_h_to_c = travel.travel_time(
                best_amb.free_location, calls[queue[ind]].location, best_amb);
            travel_times[ind][best_amb.id] = time_to_h + time_from_h_to_c;
            auto min_it = std::min_element(travel_times[ind].begin(),
                                           travel_times[ind].end());
            double min_val = *min_it;
            int min_ind = std::distance(travel_times[ind].begin(), min_it);
            min_times[k] = min_val + time - calls[queue[ind]].time;
            index_ambs[ind] = min_ind;
          }
        }
      } else {
        queue_aux.push_back(queue[remaining_indexes[current_call]]);
      }
      if (nb_treated + 1 < total_queue) {
        // erase current_call from min_times e remaining_indexes
        min_times.erase(max_it);
        remaining_indexes.erase(
            std::remove(remaining_indexes.begin(), remaining_indexes.end(),
                        remaining_indexes[current_call]),
            remaining_indexes.end());
      }
      nb_treated++;
    }

    queue.clear();
    for (auto i : queue_aux) {
      queue.push_back(i);
    }
    set_next_event(event_call, index_call);
  }
}

CGSolver::CGSolver(GRBEnv &env, vector<Call> &calls,
                   vector<Ambulance> &ambulances, Instance &ins, Travel &travel,
                   double time)
    : Solver(env, calls, ambulances, ins, travel, time) {
  prepare();
  g_params.extended_model = false;
  travel.set_forward(false);
}

void CGSolver::run() {
  if (calls.size() == 0) {
    obj = 0;
    fmt::print("No calls\n");
    return;
  }

  if (event_call == 1) {
    Data data(this, &calls[index_call], NULL);
    auto t0 = std::chrono::high_resolution_clock::now();
    CGCall cg(data, env);
    cg.solve();
    auto dt = std::chrono::high_resolution_clock::now();
    double time_cg =
        std::chrono::duration_cast<chrono::milliseconds>(dt - t0).count();
    obj = cg.model.get(GRB_DoubleAttr_ObjVal);
    run_time = time_cg;
    index_solver = 1;
  } else {
    Data data(this, NULL, &ambulances[released_amb]);
    auto t0 = std::chrono::high_resolution_clock::now();
    CGAmbulance cg(data, env);
    cg.solve();
    auto dt = std::chrono::high_resolution_clock::now();
    double time_cg =
        std::chrono::duration_cast<chrono::milliseconds>(dt - t0).count();
    obj = cg.model.get(GRB_DoubleAttr_ObjVal);
    run_time = time_cg;
    index_solver = 2;
  }
}

int CGSolver::get_return_base(CGCall &cg, Data &data, int t, int c, int a,
                              int l1, int l, int h) {
  int end_t =
      data.get_time_slot(time + t * data.quantum + data.tao[t][c][a][l1][l][h]);
  fmt::print("End t = {} {}\n", end_t,
             time + t * data.quantum + data.tao[t][c][a][l1][l][h]);
  fmt::print("quantum = {}\n", data.quantum);
  for (int b = 0; b < data.num_bases; ++b) {
    auto name = cg.yt_ahb[end_t][a][h][b].get(GRB_StringAttr_VarName);
    auto val = cg.yt_ahb[end_t][a][h][b].get(GRB_DoubleAttr_X);
    fmt::print("{} = {}\n", name, val);
    if (val > 0.5) {
      return b;
    }
  }

  for (int c1 = 0; c1 < data.types_call; ++c1) {
    if (a <= c1) {
      for (int l2 = 0; l2 < data.num_locals; ++l2) {
        for (auto h1 : data.H[c1][l2]) {
          auto name =
              cg.xt_cahlh[end_t][c1][a][h][l2][h1].get(GRB_StringAttr_VarName);
          auto val = cg.xt_cahlh[end_t][c1][a][h][l2][h1].get(GRB_DoubleAttr_X);
          if (val > g_params.EPS) {
            fmt::print("{} {}\n", name, val);
            return ins.nearest_base_to_region[l2];
          }
        }
      }
    }
  }
  int ind_b = -1;
  double min_d = GRB_INFINITY;
  for (int b = 0; b < data.num_bases; ++b) {
    if (travel.vehicle_distance(ins.hospitals[h], ins.bases[b]) < min_d) {
      double d = travel.vehicle_distance(ins.hospitals[h], ins.bases[b]);
      if (d < min_d) {
        min_d = d;
        ind_b = b;
      }
    }
  }

  return ind_b;
}

MinMaxSolver::MinMaxSolver(GRBEnv &env, vector<Call> &calls,
                           vector<Ambulance> &ambulances, Instance &ins,
                           Travel &travel, double time)
    : Solver(env, calls, ambulances, ins, travel, time) {
  prepare();
  travel.set_forward(false);
}

void MinMaxSolver::run() {
  calls_attended = 0;
  travel.set_forward(true);

  bool debug = debug_mode;

  if (debug) {
    fmt::print("Calls:\n");
    for (auto &call : calls) {
      cout << call << "\n";
    }
  }

  default_random_engine u_gen;
  uniform_real_distribution<double> u(0, 1);
  // uniform_int_distribution<int> u(0, 100);

  while (static_cast<size_t>(calls_attended) < calls.size()) {
    // travel.set_delay_param(u(u_gen));
    if (event_call == 1 && (size_t)index_call < calls.size()) {
      queue.push_back(index_call);
    }

    if (debug) {
      fmt::print(
          "(Minmax) Time = {:.1f} , index_call = {}, event_call = {}, "
          "attended = {} / {}, queue = {}\n",
          time, index_call, event_call, calls_attended, calls.size(), queue);
    }

    std::vector<double> min_times_p(queue.size(), GRB_INFINITY);
    std::vector<std::vector<double>> travel_times(
        queue.size(), std::vector<double>(ambulances.size(), GRB_INFINITY));
    std::vector<std::vector<double>> pen_travel_times(
        queue.size(), std::vector<double>(ambulances.size(), GRB_INFINITY));
    vector<int> index_ambs(queue.size(), -1);
    std::vector<vector<int>> best_ambs_for_call(queue.size(), vector<int>());

    // For every call in queue, determine the best ambulance and the best
    // response time
    for (size_t i = 0; i < queue.size(); ++i) {
      auto &call = calls[queue[i]];
      if (debug) {
        cout << "Call " << call << "\n";
      }
      for (auto &amb : ambulances) {
        int j = amb.id;
        if (can_answer(amb, call)) {
          travel_times[i][j] = travel.get_response_time(amb, call, time);
          pen_travel_times[i][j] = penalized_response_time(
              travel_times[i][j] + time - call.time, amb.type, call.priority);
          if (debug) {
            bool busy = amb.arrival_time_at_f_last_trip > time;
            fmt::print(
                "\tamb {}, type {}, time {:.3f}, call_time {:.3f}, traveltime "
                "{:.3f}, pen_time {:.3f}, "
                "busy {}\n ",
                amb.id, amb.type, time, call.time, travel_times[i][j],
                pen_travel_times[i][j], busy ? "y" : "n");
          }
          if (pen_travel_times[i][j] < min_times_p[i]) {
            min_times_p[i] = pen_travel_times[i][j];
            best_ambs_for_call[i].clear();
            best_ambs_for_call[i].push_back(j);
          } else if (abs(pen_travel_times[i][j] - min_times_p[i]) <
                     g_params.EPS) {
            best_ambs_for_call[i].push_back(j);
          }
        }
      }
      index_ambs[i] = get_best_index_amb(best_ambs_for_call[i]);
      if (debug) {
        fmt::print(
            "Call {} best_ambs = {}, index_amb = {}, min_time_p = {:1f}\n",
            call.id, best_ambs_for_call[i], index_ambs[i], min_times_p[i]);
      }
    }

    std::vector<int> remaining_indexes(queue.size(), 0);
    for (size_t i = 0; i < queue.size(); ++i) {
      remaining_indexes[i] = i;
    }

    std::vector<int> queue_aux;
    int nb_treated = 0;
    int total_queue = queue.size();

    while (nb_treated < total_queue) {
      auto t0 = std::chrono::high_resolution_clock::now();
      // get time and index of the worst call, and the best ambulance to such
      // call
      if (debug) {
        fmt::print("Nb_treated = {}\n", nb_treated);
        fmt::print("Queue = {}\n", queue);
        fmt::print("Remaining indexes = {}\n", remaining_indexes);
        // fmt::print("Min_times = {}\n", min_times);
        fmt::print("Min_times_p = [{:0.2f}]\n", fmt::join(min_times_p, ", "));
      }

      double max_time = -1;
      int current_call = -1;
      int index_amb = -1;
      // Chooses the call with maximum penalized waiting time. If two calls
      // realize the max_time, chooses one that has an available ambulance, if
      // possible.
      for (auto k : remaining_indexes) {
        if (min_times_p[k] > max_time) {
          max_time = min_times_p[k];
          current_call = k;
          index_amb = index_ambs[current_call];
        } else if (abs(min_times_p[k] - max_time) < g_params.EPS &&
                   index_amb == -1) {
          index_amb = get_best_index_amb(best_ambs_for_call[k]);
          current_call = k;
        }
      }

      int call_ind = queue[current_call];
      auto &call = calls[call_ind];
      auto &best_amb = ambulances[index_amb];

      if (debug) {
        fmt::print(
            "Treating call {}, min_time_p = {:.1f}, remaining_indexes = "
            "{}, min_times = {}, current_call = {}\n",
            call_ind, min_times_p[current_call], remaining_indexes, min_times_p,
            current_call);
      }

      // Best free ambulance should answer the call
      if (index_amb >= 0) {
        if (which_ambulance[call_ind] >= 0) {
          std::string aux = fmt::format(fg(fmt::color(0xFF0000)),
                                        "Call {} was already answered by {}\n",
                                        call_ind, which_ambulance[call_ind]);
          fmt::print("{}", aux);
          cin.get();
        }
        double waiting_on_scene_i =
            time + travel_times[current_call][index_amb] - call.time;
        int base_id = (g_params.best_base) ? find_best_base(index_amb)
                                           : nearest_base[call_ind];
        double waiting_to_hospital_i = ambulances[index_amb].answer_call(
            call, travel, ins, time, waiting_on_scene_i - (time - call.time),
            base_id);
        waiting_on_scene[call_ind] = waiting_on_scene_i;
        waiting_on_scene_penalized[call_ind] = penalized_response_time(
            waiting_on_scene_i, best_amb.type, call.priority);
        waiting_to_hospital[call_ind] = waiting_to_hospital_i;
        which_ambulance[call_ind] = index_amb;
        calls_end[call_ind] = call.end;

        if (debug) {
          std::string dispatch_msg = fmt::format(
              fg(fmt::color(0xFFFF00)),
              "Call {} answered by {}. w_scene = {}, free = {}\n", call_ind,
              index_amb, waiting_on_scene_penalized[call_ind],
              best_amb.arrival_time_at_f_last_trip);
          fmt::print("{}", dispatch_msg);
        }

        calls_attended++;
        // update travel times
        if (debug) {
          fmt::print("Updated travel times:\n");
        }
        for (auto ind : remaining_indexes) {
          auto has_index_amb = find(best_ambs_for_call[ind].begin(),
                                    best_ambs_for_call[ind].end(),
                                    index_amb) != best_ambs_for_call[ind].end();
          if (ind != current_call && has_index_amb) {
            double time_to_h = best_amb.arrival_time_at_f_last_trip - time;
            double time_from_h_to_c = travel.travel_time(
                best_amb.free_location, calls[queue[ind]].location, best_amb);
            travel_times[ind][best_amb.id] = time_to_h + time_from_h_to_c;
            pen_travel_times[ind][best_amb.id] = penalized_response_time(
                time + travel_times[ind][best_amb.id] - calls[queue[ind]].time,
                best_amb.type, calls[queue[ind]].priority);
            double min_time_p = GRB_INFINITY;
            vector<int> best_ambs_for_k;
            for (size_t j = 0; j < ambulances.size(); ++j) {
              double pen_waiting_time = penalized_response_time(
                  travel_times[ind][j] + time - calls[queue[ind]].time,
                  ambulances[j].type, calls[queue[ind]].priority);
              pen_travel_times[ind][j] = pen_waiting_time;
              if (debug) {
                fmt::print("\tcall {} amb {}: allocation cost = {}\n",
                           queue[ind], j, pen_travel_times[ind][j]);
              }
              if (pen_waiting_time < min_time_p) {
                min_time_p = pen_waiting_time;
                best_ambs_for_k.clear();
                best_ambs_for_k.push_back(j);
              } else if (abs(pen_waiting_time - min_time_p) < g_params.EPS) {
                best_ambs_for_k.push_back(j);
              }
            }
            index_ambs[ind] = get_best_index_amb(best_ambs_for_k);
            min_times_p[ind] = min_time_p;
            if (debug) {
              fmt::print("\tCall {}: pen = {:.3f}, index_amb = {}\n",
                         queue[ind], min_times_p[ind], index_ambs[ind]);
            }
          }
        }
      } else {
        if (debug) {
          fmt::print("Call {} back to queue\n", call_ind);
        }
        queue_aux.push_back(queue[current_call]);
      }
      if (nb_treated < total_queue) {
        // min_times.erase(min_times.begin() + current_call);
        // min_times_p.erase(min_times_p.begin() + current_call);
        remaining_indexes.erase(remove(remaining_indexes.begin(),
                                       remaining_indexes.end(), current_call),
                                remaining_indexes.end());
      }
      nb_treated++;
      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
    }

    queue.clear();
    for (auto i : queue_aux) {
      queue.push_back(i);
    }
    int former_event = event_call;
    set_next_event(event_call, index_call);
    if (debug) {
      fmt::print("Next event = {}, time {:.1f}, index = {}, queue = {}\n",
                 event_call, time, index_call, queue);
      cin.get();
    }
  }
}

PreparednessSolver::PreparednessSolver(GRBEnv &env, vector<Call> &calls,
                                       vector<Ambulance> &ambulances,
                                       Instance &ins, Travel &travel, int g,
                                       double time)
    : Solver(env, calls, ambulances, ins, travel, time), g(g) {
  prepare();
  travel.set_forward(false);
}

void PreparednessSolver::run() {
  bool debug = false;
  vector<Ambulance> available;
  available.reserve(ins.nb_ambulances);
  vector<double> delta_threshold{300, 600, 300, 600};
  vector<double> lambda(ins.nb_regions, 0);
  for (int r = 0; r < ins.nb_regions; ++r) {
    double sum = 0;
    for (auto t : ins.time_horizon) {
      for (int p = 0; p < ins.nb_priorities; ++p) {
        sum += ins.lambda(t, g, r, p);
      }
    }
    lambda[r] = sum;
  }
  if (debug) {
    fmt::print("Lambda:\n");
    for (size_t r = 0; r < lambda.size(); ++r) {
      if (lambda[r] > 0) {
        fmt::print("lam[{}] = {}\n", r, lambda[r]);
      }
    }
  }

  double first_time = time;
  travel.set_forward(false);
  calls_attended = 0;
  default_random_engine u_gen;
  uniform_real_distribution<double> u(0, 1);
  while (static_cast<size_t>(calls_attended) < calls.size()) {
    travel.set_delay_param(u(u_gen));
    // for (auto &amb : ambulances) {
    //   if (amb.arrival_time_at_b_last_trip <= time) {
    //     for (auto &call : calls) {
    //       call.ambulances.push_back(amb.id);
    //     }
    //   }
    // }
    if (debug) {
      fmt::print(
          "Time = {:.1f}, event = {}, index_call = {}, released = {}, "
          "attended = {}/{}\n",
          time, event_call, index_call, released_amb, calls_attended,
          calls.size());
    }
    if (event_call == 1) {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto &call = calls[index_call];
      if (debug) {
        cout << call << ", threshold = " << delta_threshold[call.priority]
             << "\n";
      }
      available.clear();
      double min_time = GRB_INFINITY;
      int index_amb = -1;
      for (auto &amb : ambulances) {
        double response_time = travel.get_response_time(amb, call, time);
        double pen_resp_time = (response_time < GRB_INFINITY)
                                   ? penalized_response_time(
                                         response_time, amb.type, call.priority)
                                   : GRB_INFINITY;
        if (debug) {
          fmt::print("\tamb {} , time = {:.2f}, pen = {:.2f}\n", amb.id,
                     response_time, pen_resp_time);
        }
        if (response_time < delta_threshold[call.priority] &&
            amb.arrival_time_at_f_last_trip <= time && can_answer(amb, call)) {
          available.push_back(amb);
        }

        if (amb.arrival_time_at_f_last_trip <= time &&
            response_time < min_time && can_answer(amb, call)) {
          min_time = response_time;
          index_amb = amb.id;
        }
      }
      if (debug) {
        fmt::print("Available:\n");
        for (auto &amb : available) {
          cout << amb << "\n";
        }
        fmt::print("END AVAILABLE\n");
      }

      if (available.size() > 0) {
        vector<vector<double>> prep(
            available.size(), vector<double>(ins.nb_regions, GRB_INFINITY));
        vector<double> area_prep(available.size(), GRB_INFINITY);
        for (size_t i = 0; i < available.size(); ++i) {
          double min_prep = GRB_INFINITY;
          for (int j = 0; j < ins.nb_regions; ++j) {
            double sum = 0;
            for (auto &amb : available) {
              if (amb.id != available[i].id &&
                  amb.arrival_time_at_b_last_trip <= time) {
                sum += 1 / travel.travel_time(amb.base_location, ins.centers[j],
                                              amb);
              } else if (amb.id != available[i].id &&
                         amb.arrival_time_at_f_last_trip <= time) {
                auto current_location = travel.ambulance_position(amb, time);
                sum += 1 / travel.travel_time(current_location, ins.centers[j],
                                              amb);
              }
            }
            prep[i][j] = (1 / lambda[j]) * sum;
            if (prep[i][j] < min_prep) {
              min_prep = prep[i][j];
            }
          }
          area_prep[i] = min_prep;
        }
        if (debug) {
          fmt::print("Area prep = {}\n", area_prep);
        }

        auto max_prep = max_element(area_prep.begin(), area_prep.end());
        int index_amb = distance(area_prep.begin(), max_prep);
        if (debug) {
          fmt::print("Max prep amb = {}\n", available[index_amb].id);
        }
        auto &amb = ambulances[available[index_amb].id];
        double real_min_time = travel.get_response_time(amb, call, time);
        double waiting_on_scene_i = time + real_min_time - call.time;
        double waiting_to_hospital_i = amb.answer_call(
            call, travel, ins, time, real_min_time, nearest_base[call.id]);
        waiting_on_scene[call.id] = waiting_on_scene_i;
        waiting_on_scene_penalized[call.id] = penalized_response_time(
            waiting_on_scene_i, amb.type, call.priority);
        if (debug) {
          fmt::print("Time {} Call {} attended by {} (prep rule)\n", time,
                     call.id, amb.id);
        }
        // fmt::print("| {:2f} {:.2f}\n", waiting_on_scene[call.id],
        // 	waiting_on_scene_penalized[call.id]);
        // fmt::print("Call location {}, amb location {}\n", call.location,
        // 	ambulances[index_amb].base_location);
        waiting_to_hospital[call.id] = waiting_to_hospital_i;
        which_ambulance[call.id] = index_amb;
        calls_end[call.id] = call.end;
        ++calls_attended;
      } else if (index_amb >= 0) {
        auto &amb = ambulances[index_amb];
        double real_min_time = travel.get_response_time(amb, call, time);
        double waiting_on_scene_i = time + real_min_time - call.time;
        double waiting_to_hospital_i = ambulances[index_amb].answer_call(
            call, travel, ins, time, real_min_time, nearest_base[call.id]);
        waiting_on_scene[call.id] = waiting_on_scene_i;
        waiting_on_scene_penalized[call.id] = penalized_response_time(
            waiting_on_scene_i, amb.type, call.priority);
        if (debug) {
          fmt::print("Call {} attended by {} (closest)\n", call.id, index_amb);
        }
        // fmt::print("| {:2f} {:.2f}\n", waiting_on_scene[call.id],
        // 	waiting_on_scene_penalized[call.id]);
        // fmt::print("Call location {}, amb location {}\n", call.location,
        // 	ambulances[index_amb].base_location);
        waiting_to_hospital[call.id] = waiting_to_hospital_i;
        which_ambulance[call.id] = index_amb;
        calls_end[call.id] = call.end;
        ++calls_attended;
      } else {
        queue.push_back(index_call);
      }
      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
    } else {
      // dispatch ambulance to closest call in queue.
      if (released_amb == -1) {
        fmt::print(
            "Error: Not event call but no released amb also. queue = {}\n",
            queue);
        std::cin.get();
      }
      auto &amb = ambulances[released_amb];
      double min_time = GRB_INFINITY;
      int min_call = -1;
      auto t0 = std::chrono::high_resolution_clock::now();
      tie(min_call, min_time) = get_closest_call(amb);
      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
      if (min_call >= 0) {
        auto &call = calls[min_call];
        double real_min_time = travel.get_response_time(amb, call, time);
        double waiting_on_scene_i = time + real_min_time - call.time;
        double waiting_to_hospital_i = amb.answer_call(
            call, travel, ins, time, real_min_time, nearest_base[call.id]);
        waiting_on_scene[call.id] = waiting_on_scene_i;
        waiting_on_scene_penalized[call.id] = penalized_response_time(
            waiting_on_scene_i, amb.type, call.priority);
        if (debug) {
          fmt::print("Call {} attended by {} (return)\n", min_call,
                     released_amb);
        }
        // fmt::print("| {:2f} {:.2f}\n", waiting_on_scene[call.id],
        // 	waiting_on_scene_penalized[call.id]);
        // fmt::print("Call location {}, amb location {}\n", call.location,
        // 	ambulances[amb.id].base_location);
        waiting_to_hospital[call.id] = waiting_to_hospital_i;
        which_ambulance[call.id] = amb.id;
        calls_end[call.id] = call.end;
        ++calls_attended;
        queue.erase(remove(queue.begin(), queue.end(), min_call), queue.end());
      } else {
        vector<int> reloc_available;
        for (auto &amb : ambulances) {
          if (amb.arrival_time_at_f_last_trip <= time) {
            reloc_available.push_back(amb.id);
          }
        }
        auto &amb = ambulances[released_amb];
        int b = get_return_base(released_amb, lambda, reloc_available);
        if (debug) {
          fmt::print("Ambulance  {} returned to base {}\n", released_amb, b);
        }
        amb.arrival_time_at_b_last_trip =
            travel.travel_time(amb.free_location, ins.bases[b], amb);
        amb.base_location = ins.bases[b];
      }
    }
    set_next_event(event_call, index_call);
    if (debug) {
      cin.get();
    }
  }
  if (debug) {
    print_results();
    cin.get();
  }
}

int PreparednessSolver::get_return_base(int amb_id, vector<double> &lambda,
                                        vector<int> &available) {
  auto &amb = ambulances[amb_id];
  vector<pair<double, int>> time_to_base;
  time_to_base.reserve(ins.nb_bases);
  for (int b = 0; b < ins.nb_bases; ++b) {
    double travel_time =
        travel.travel_time(amb.free_location, ins.bases[b], amb);
    time_to_base.push_back(make_pair(travel_time, b));
  }
  sort(time_to_base.begin(), time_to_base.end());

  size_t k = 0;
  vector<pair<double, int>> prep_to_base;
  prep_to_base.reserve(5);
  while (k < 5 && k < time_to_base.size()) {
    int b = time_to_base[k].second;
    double min_prep = GRB_INFINITY;
    int min_r = -1;
    for (int r = 0; r < ins.nb_regions; ++r) {
      double sum = 0;
      for (auto i : available) {
        auto &other_amb = ambulances[i];
        if (i == amb_id) {
          sum += 1 / travel.travel_time(ins.bases[b], ins.centers[r], amb);
        } else {
          if (other_amb.arrival_time_at_b_last_trip <= time) {
            sum += 1 / travel.travel_time(other_amb.base_location,
                                          ins.centers[r], other_amb);
          } else {
            auto current_location = travel.ambulance_position(other_amb, time);
            sum += 1 / travel.travel_time(current_location, ins.centers[r],
                                          other_amb);
          }
        }
      }
      double prep = (1 / lambda[r]) * sum;
      if (prep < min_prep) {
        min_prep = prep;
      }
    }
    prep_to_base.push_back(make_pair(min_prep, b));
    ++k;
  }
  sort(prep_to_base.begin(), prep_to_base.end(), greater<pair<double, int>>());
  int result = -1;
  if (prep_to_base[0].first > g_params.min_preparedness) {
    double max_prep = prep_to_base[0].first;
    int best_b = prep_to_base[0].second;
    double min_time =
        travel.travel_time(amb.free_location, ins.bases[best_b], amb);
    for (size_t i = 1; i < prep_to_base.size(); ++i) {
      auto prep = prep_to_base[i].first;
      int b = prep_to_base[i].second;
      if (prep < max_prep) {
        break;
      } else {
        double travel_time =
            travel.travel_time(amb.free_location, ins.bases[b], amb);
        if (travel_time < min_time) {
          min_time = travel_time;
          best_b = b;
        }
      }
    }
    result = best_b;
  } else {
    result = prep_to_base[0].second;
  }
  return result;
}

Prep2Solver::Prep2Solver(GRBEnv &env, vector<Call> &calls,
                         vector<Ambulance> &ambulances, Instance &ins,
                         Travel &travel, int g, double time)
    : Solver(env, calls, ambulances, ins, travel, time), g(g) {
  prepare();
  travel.set_forward(false);
}

void Prep2Solver::run() {
  bool debug = debug_mode;
  vector<pair<double, int>> available;
  available.reserve(ins.nb_ambulances);
  vector<double> delta_threshold{300, 600, 300, 600};
  vector<double> lambda(ins.nb_regions, 0);

  if (debug) {
    fmt::print("Calls:\n");
    for (auto &call : calls) {
      cout << call << " " << fmt::format("{}\n", call.ambulances);
    }
    fmt::print("Ambs:\n");
    for (auto &amb : ambulances) {
      cout << amb << "\n";
    }
    cin.get();
  }

  for (int r = 0; r < ins.nb_regions; ++r) {
    double sum = 0;
    // TODO: switch to 0 when considering all time slots
    for (int t = 36; t < 40; ++t) {
      for (int p = 0; p < ins.nb_priorities; ++p) {
        sum += ins.lambda(t, g, r, p);
      }
    }
    lambda[r] = sum;
  }

  calls_attended = 0;
  travel.set_forward(false);
  default_random_engine u_gen;
  uniform_real_distribution<double> u(0, 1);

  while (static_cast<size_t>(calls_attended) < calls.size()) {
    available.clear();
    travel.set_delay_param(u(u_gen));
    // for(auto& amb: ambulances){
    // 	if(amb.arrival_time_at_b_last_trip <= time){
    // 		for(auto& call: calls){
    // 			call.ambulances.push_back(amb.id);
    // 		}
    // 	}
    // }

    if (debug) {
      fmt::print(
          "Prep2: Time {}, event {}, index {}, attended {}/{}, queue {}\n",
          time, event_call, index_call, calls_attended, calls.size(), queue);
    }

    if (event_call == 1) {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto &call = calls[index_call];
      double min_time = GRB_INFINITY;
      double p_min = 0;
      for (Ambulance &amb : ambulances) {
        if (debug) {
          fmt::print("Testing ambulance {} for call {}\n", amb.id, index_call);
        }
        double response_time = travel.get_response_time(amb, call, time);
        if (debug) {
          fmt::print("Response_time amb {}  call {} = {}\n", amb.id, index_call,
                     response_time);
        }
        double pen_response_time =
            penalized_response_time(response_time, amb.type, call.priority);
        if (can_answer(amb, call) && amb.arrival_time_at_f_last_trip <= time) {
          available.push_back(make_pair(response_time, amb.id));
        }

        if (can_answer(amb, call) && response_time < min_time) {
          min_time = response_time;
        }
      }
      std::sort(available.begin(), available.end());
      bool low_priority = (g_params.amb_setup == "rj" && call.priority > 1) ||
                          (g_params.amb_setup == "us" &&
                           (call.priority == 1 || call.priority == 3));
      // fmt::print("Available.size() = {}\n", available.size());
      if (available.size() > 1 &&
          low_priority) {  // p \in {1,3} -> low_priority
        bool index_av = -1;
        for (size_t i = 0; i < available.size(); ++i) {
          double travel_time;
          int amb_id;
          tie(travel_time, amb_id) = available[i];
          if (travel_time < delta_threshold[call.priority]) {
            double min_prep = GRB_INFINITY;
            for (int r = 0; r < ins.nb_regions; ++r) {
              double sum = 0;
              for (auto &av : available) {
                auto &amb = ambulances[av.second];
                if (amb.id != amb_id &&
                    amb.arrival_time_at_b_last_trip <= time) {
                  sum += 1 / travel.travel_time(amb.base_location,
                                                ins.centers[r], amb);
                } else if (amb.id != amb_id &&
                           amb.arrival_time_at_f_last_trip <= time) {
                  auto current_location = travel.ambulance_position(amb, time);
                  sum += 1 / travel.travel_time(current_location,
                                                ins.centers[r], amb);
                }
              }
              double prep = (1 / lambda[r]) *
                            sum;  // prep of region r when i is unavailable
              if (prep < min_prep) {
                min_prep = prep;
              }
            }
            if (min_prep > p_min) {
              p_min = min_prep;
              index_av = i;
            }
          } else {
            index_av = i;
            break;
          }
        }
        // dispatch index_amb;
        min_time = available[index_av].first;
        auto &index_amb = available[index_av].second;
        auto &amb = ambulances[index_amb];
        double real_min_time = travel.get_response_time(amb, call, time);
        double waiting_on_scene_i = time + real_min_time - call.time;
        double waiting_to_hospital_i = amb.answer_call(
            call, travel, ins, time, real_min_time, nearest_base[call.id]);
        waiting_on_scene[call.id] = waiting_on_scene_i;
        waiting_on_scene_penalized[call.id] = penalized_response_time(
            waiting_on_scene_i, amb.type, call.priority);
        // fmt::print("Call {} attended by {} (prep) ", call.id, index_amb);
        // fmt::print("| {:2f} {:.2f}\n", waiting_on_scene[call.id],
        // waiting_on_scene_penalized[call.id]);
        // fmt::print("Call location {}, amb location {}\n", call.location,
        // 	ambulances[index_amb].base_location);
        waiting_to_hospital[call.id] = waiting_to_hospital_i;
        which_ambulance[call.id] = index_amb;
        calls_end[call.id] = call.end;
        ++calls_attended;
      } else if (available.size() > 0) {
        // dispatch to available[0]
        min_time = available[0].first;
        auto &index_amb = available[0].second;
        auto &amb = ambulances[index_amb];
        double real_min_time = travel.get_response_time(amb, call, time);
        double waiting_on_scene_i = time + real_min_time - call.time;
        double waiting_to_hospital_i = ambulances[index_amb].answer_call(
            call, travel, ins, time, real_min_time, nearest_base[call.id]);
        waiting_on_scene[call.id] = waiting_on_scene_i;
        waiting_on_scene_penalized[call.id] = penalized_response_time(
            waiting_on_scene_i, amb.type, call.priority);

        // fmt::print("Call {} attended by {} (closest) ", call.id, index_amb);
        // fmt::print("| {:2f} {:.2f}\n", waiting_on_scene[call.id],
        // waiting_on_scene_penalized[call.id]);
        // fmt::print("Call location {}, amb location {}\n", call.location,
        // 	ambulances[index_amb].base_location);
        waiting_to_hospital[call.id] = waiting_to_hospital_i;
        which_ambulance[call.id] = index_amb;
        calls_end[call.id] = call.end;
        ++calls_attended;
      } else {
        // fmt::print("Call {} sent to queue\n", call.id);
        queue.push_back(index_call);
      }

      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
    } else {
      if (queue.size() > 0) {
        // dispatch to oldest call
        if (released_amb == -1) {
          fmt::print("Error: Not event call but no released amb also.\n");
          cin.get();
        }
        auto &amb = ambulances[released_amb];
        auto t0 = std::chrono::high_resolution_clock::now();
        int min_call;
        double min_time;
        tie(min_call, min_time) = get_oldest_call(amb);
        auto dt = std::chrono::high_resolution_clock::now();
        run_times.push_back(
            std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
        if (min_call >= 0) {
          auto &call = calls[min_call];
          double real_min_time = travel.get_response_time(amb, call, time);
          double waiting_on_scene_i = time + real_min_time - call.time;
          double waiting_to_hospital_i = amb.answer_call(
              call, travel, ins, time, real_min_time, nearest_base[call.id]);
          waiting_on_scene[call.id] = waiting_on_scene_i;
          waiting_on_scene_penalized[call.id] = penalized_response_time(
              waiting_on_scene_i, amb.type, call.priority);
          // fmt::print("Call {} attended by {} (return)\n", min_call,
          // released_amb); fmt::print("| {:2f} {:.2f}\n",
          // waiting_on_scene[call.id], waiting_on_scene_penalized[call.id]);
          // fmt::print("Call location {}, amb location {}\n", call.location,
          // 	ambulances[amb.id].base_location);
          waiting_to_hospital[call.id] = waiting_to_hospital_i;
          which_ambulance[call.id] = amb.id;
          calls_end[call.id] = call.end;
          ++calls_attended;
          queue.erase(remove(queue.begin(), queue.end(), min_call),
                      queue.end());
        }
      } else {
        vector<int> reloc_available;
        for (auto &amb : ambulances) {
          if (amb.arrival_time_at_f_last_trip <= time) {
            reloc_available.push_back(amb.id);
          }
        }
        auto &amb = ambulances[released_amb];
        int b = get_return_base(released_amb, lambda, reloc_available);
        // fmt::print("Ambulance  {} returned to base {}\n", released_amb, b);
        amb.arrival_time_at_b_last_trip =
            travel.travel_time(amb.free_location, ins.bases[b], amb);
        amb.base_location = ins.bases[b];
      }
    }
    set_next_event(event_call, index_call);
    // cin.get();
  }
  // if(debug){
  // 	print_results();
  // }
}

int Prep2Solver::get_return_base(int amb_id, vector<double> &lambda,
                                 vector<int> &available) {
  auto &amb = ambulances[amb_id];
  vector<pair<double, int>> time_to_base;
  time_to_base.reserve(ins.nb_bases);
  for (int b = 0; b < ins.nb_bases; ++b) {
    double travel_time =
        travel.travel_time(amb.free_location, ins.bases[b], amb);
    time_to_base.push_back(make_pair(travel_time, b));
  }
  sort(time_to_base.begin(), time_to_base.end());

  size_t k = 0;
  vector<pair<double, int>> prep_to_base;
  prep_to_base.reserve(5);
  while (k < 5 && k < time_to_base.size()) {
    int b = time_to_base[k].second;
    double min_prep = GRB_INFINITY;
    int min_r = -1;
    for (int r = 0; r < ins.nb_regions; ++r) {
      double sum = 0;
      for (auto i : available) {
        auto &other_amb = ambulances[i];
        if (i == amb_id) {
          sum += 1 / travel.travel_time(ins.bases[b], ins.centers[r], amb);
        } else {
          if (other_amb.arrival_time_at_b_last_trip <= time) {
            sum += 1 / travel.travel_time(other_amb.base_location,
                                          ins.centers[r], other_amb);
          } else {
            auto current_location = travel.ambulance_position(other_amb, time);
            sum += 1 / travel.travel_time(current_location, ins.centers[r],
                                          other_amb);
          }
        }
      }
      double prep = (1 / lambda[r]) * sum;
      if (prep < min_prep) {
        min_prep = prep;
      }
    }
    prep_to_base.push_back(make_pair(min_prep, b));
    ++k;
  }
  sort(prep_to_base.begin(), prep_to_base.end(), greater<pair<double, int>>());
  int result = -1;
  if (prep_to_base[0].first > g_params.min_preparedness) {
    double max_prep = prep_to_base[0].first;
    int best_b = prep_to_base[0].second;
    double min_time =
        travel.travel_time(amb.free_location, ins.bases[best_b], amb);
    for (size_t i = 1; i < prep_to_base.size(); ++i) {
      auto prep = prep_to_base[i].first;
      int b = prep_to_base[i].second;
      if (prep < max_prep) {
        break;
      } else {
        double travel_time =
            travel.travel_time(amb.free_location, ins.bases[b], amb);
        if (travel_time < min_time) {
          min_time = travel_time;
          best_b = b;
        }
      }
    }
    result = best_b;
  } else {
    result = prep_to_base[0].second;
  }
  return result;
}

OrderedSolver::OrderedSolver(GRBEnv &env, vector<Call> &calls,
                             vector<Ambulance> &ambulances, Instance &ins,
                             Travel &travel, int g, double time)
    : Solver(env, calls, ambulances, ins, travel, time), g(g) {
  prepare();
}

void OrderedSolver::prepare() {
  bool debug = false;
  int first_new_call = 0;
  while (static_cast<size_t>(first_new_call) < calls.size() &&
         calls[first_new_call].time <= time) {
    queue.push_back(first_new_call);
    ++first_new_call;
  }
  index_call = first_new_call;
  if (debug) {
    fmt::print("Time {}\n", time);
    for (auto &call : calls) {
      cout << call << "\n";
    }
    fmt::print("queue = {}, index_call = {}, calls.size() = {}\n", queue,
               index_call, calls.size());
  }
  vector<pair<double, int>> future_arrival_times;
  for (size_t i = 0; i < ambulances.size(); ++i) {
    auto &amb = ambulances[i];
    if (debug) {
      cout << amb << "\n";
    }
    if (amb.arrival_time_at_b_last_trip >= time) {
      future_arrival_times.push_back(
          make_pair(amb.arrival_time_at_f_last_trip, amb.id));
    }
  }
  if (future_arrival_times.size() > 0) {
    auto min_arrival =
        min_element(future_arrival_times.begin(), future_arrival_times.end());
    if (debug) {
      fmt::print("min_arrival = {}\n", *min_arrival);
    }
    if (static_cast<size_t>(index_call) >= calls.size() ||
        min_arrival->first <= calls[index_call].time) {
      event_call = 0;
      released_amb = min_arrival->second;
      time = min_arrival->first;
    } else {
      event_call = 1;
      time = calls[index_call].time;
      released_amb = -1;
    }
  } else {
    event_call = 1;
    released_amb = -1;
    if (static_cast<size_t>(index_call) < calls.size()) {
      time = calls[index_call].time;
    }
  }

  if (debug) {
    fmt::print(
        "Ord Prepare, time = {}, index = {}, event = {}, queue = {}, "
        "released = {}\n",
        time, index_call, event_call, queue, released_amb);
    cin.get();
  }
}

void OrderedSolver::run() {
  calls_attended = 0;
  travel.set_forward(false);
  g_params.h_use_fixed_bases = true;
  bool debug = false;

  vector<double> lambda(ins.nb_regions, 0);
  for (int r = 0; r < ins.nb_regions; ++r) {
    double sum = 0;
    // TODO: switch to all time slots
    for (auto t : ins.time_horizon) {
      for (int p = 0; p < ins.nb_priorities; ++p) {
        sum += ins.lambda(t, g, r, p);
      }
    }
    lambda[r] = sum;
  }

  default_random_engine u_gen;
  uniform_real_distribution<double> u(0, 1);
  // travel.set_delay_param(u(u_gen));

  auto ordered_list = get_ordered_table(lambda);
  if (debug) {
    fmt::print("ordered list = {}\n", ordered_list);
  }

  while (static_cast<size_t>(calls_attended) < calls.size()) {
    travel.set_delay_param(u(u_gen));
    for (auto &amb : ambulances) {
      if (amb.arrival_time_at_b_last_trip <= time) {
        for (auto &call : calls) {
          call.ambulances.push_back(amb.id);
        }
      }
    }

    if (event_call == 1) {
      queue.push_back(index_call);
    }
    if (debug) {
      fmt::print(
          "Ordered | Time {}, index = {}, event = {}, attended = {}/{}, "
          "queue = {}\n",
          time, index_call, event_call, calls_attended, calls.size(), queue);
    }

    vector<int> queue_aux;
    for (auto i : queue) {
      auto &call = calls[i];
      int index_amb = -1;
      double min_time = GRB_INFINITY;
      string type_assign;
      auto t0 = std::chrono::high_resolution_clock::now();
      if (call.priority == 0 || call.priority == 2) {
        // dispatch closest available
        for (auto &amb : ambulances) {
          if (can_answer(amb, call) &&
              amb.arrival_time_at_b_last_trip <= time) {
            double resp_time = penalized_response_time(
                travel.get_response_time(amb, call, time), amb.type,
                call.priority);
            if (resp_time < min_time) {
              min_time = resp_time;
              index_amb = amb.id;
              type_assign = "closest";
            }
          }
        }
      } else {
        // dispatch from table
        for (size_t i = 0; i < ordered_list.size(); ++i) {
          auto &amb = ambulances[ordered_list[i]];
          if (amb.arrival_time_at_b_last_trip <= time) {
            index_amb = ordered_list[i];
            min_time = penalized_response_time(
                travel.get_response_time(amb, call, time), amb.type,
                call.priority);
            type_assign = "table";
            break;
          }
        }
      }

      if (index_amb >= 0) {
        Ambulance &best_amb = ambulances[index_amb];
        double real_min_time = travel.get_response_time(best_amb, call, time);
        double waiting_on_scene_i = time + real_min_time - call.time;
        double waiting_to_hospital_i = best_amb.answer_call(
            call, travel, ins, time, real_min_time, nearest_base[call.id]);
        best_amb.free_location = best_amb.base_location;
        best_amb.arrival_time_at_f_last_trip =
            best_amb.arrival_time_at_b_last_trip;
        waiting_on_scene[call.id] = waiting_on_scene_i;
        waiting_on_scene_penalized[call.id] = penalized_response_time(
            waiting_on_scene_i, best_amb.type, call.priority);
        if (debug) {
          fmt::print("Call {} attended by {} ({})\n", call.id, index_amb,
                     type_assign);
        }
        // fmt::print("| {:2f} {:.2f}\n", waiting_on_scene[call.id],
        // waiting_on_scene_penalized[call.id]);
        // fmt::print("Call location {}, amb location {}\n", call.location,
        // 	ambulances[index_amb].base_location);
        waiting_to_hospital[call.id] = waiting_to_hospital_i;
        which_ambulance[call.id] = index_amb;
        calls_end[call.id] = call.end;
        ++calls_attended;
      } else {
        if (debug) {
          fmt::print("Call {} back to queue\n", call.id);
        }
        queue_aux.push_back(call.id);
      }
      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
    }
    queue.clear();
    for (auto i : queue_aux) {
      queue.push_back(i);
    }
    set_next_event(event_call, index_call);
  }
  if (debug) {
    print_results();
    cin.get();
  }
}

vector<int> OrderedSolver::get_ordered_table(vector<double> &lambda) {
  vector<vector<pair<double, int>>> table(ins.nb_regions,
                                          vector<pair<double, int>>());
  int n = ins.nb_regions;
  int m = ambulances.size();
  for (int i = 0; i < n; ++i) {
    vector<pair<double, int>> ordered_ambs;
    for (int j = 0; j < m; ++j) {
      auto &amb = ambulances[j];
      double travel_time =
          travel.travel_time(amb.base_location, ins.centers[i], amb);
      ordered_ambs.push_back(make_pair(travel_time, amb.id));
    }
    sort(ordered_ambs.begin(), ordered_ambs.end());
    table[i] = ordered_ambs;
  }

  double sum_lambda = accumulate(lambda.begin(), lambda.end(), 0.0);
  vector<double> z(n, 0);
  for (int i = 0; i < n; ++i) {
    z[i] = lambda[i] / sum_lambda;
  }

  vector<vector<double>> B(m, vector<double>(m, GRB_INFINITY));
  for (int j = 0; j < m; ++j) {
    for (int k = 0; k < m; ++k) {
      double sum = 0;
      for (int i = 0; i < n; ++i) {
        if (table[i][j].second == k) {
          sum += z[i];
        }
      }
      B[k][j] = sum;
    }
  }

  vector<int> indexes(m, 0);
  for (int i = 0; i < m; ++i) {
    indexes[i] = i;
  }
  sort(indexes.begin(), indexes.end(), [&](const int i, const int j) {
    // i  is less than j if column i is less than column j
    for (int k = 0; k < m; ++k) {
      if (B[k][i] < B[k][j]) {
        return true;
      } else if (abs(B[k][i] - B[k][j]) < 0.0001) {
        continue;
      } else {
        return false;
      }
    }
    return false;
  });
  return indexes;
}

CoverageSolver::CoverageSolver(GRBEnv &env, vector<Call> &calls,
                               vector<Ambulance> &ambulances, Instance &ins,
                               Travel &travel, int g, double time)
    : Solver(env, calls, ambulances, ins, travel, time), g(g), T(900), q(0.8) {
  prepare();
}

void CoverageSolver::run() {
  bool debug = false;
  vector<double> lambda(ins.nb_regions, 0);
  for (int r = 0; r < ins.nb_regions; ++r) {
    double sum = 0;
    for (int t = 36; t < 40; ++t) {
      for (int p = 0; p < ins.nb_priorities; ++p) {
        sum += ins.lambda(t, g, r, p);
      }
    }
    lambda[r] = sum;
  }

  calls_attended = 0;
  travel.set_forward(false);

  vector<int> k(ins.nb_regions, 0);

  default_random_engine u_gen;
  uniform_real_distribution<double> u(0, 1);

  while (static_cast<size_t>(calls_attended) < calls.size()) {
    travel.set_delay_param(u(u_gen));
    // for(auto& amb: ambulances){
    // 	if(amb.arrival_time_at_b_last_trip <= time){
    // 		for(auto& call: calls){
    // 			call.ambulances.push_back(amb.id);
    // 		}
    // 	}
    // }

    if (event_call == 1) {
      queue.push_back(index_call);
    }

    if (debug) {
      fmt::print(
          "CoverageSolver: Time {}, index = {}, event = {}, queue = {}, "
          "att = {}/{}\n",
          time, index_call, event_call, queue, calls_attended, calls.size());
    }

    vector<int> queue_aux;
    for (auto i : queue) {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto &call = calls[i];

      vector<int> A;
      vector<int> A_plus;
      vector<int> A_minus;
      A_plus.reserve(ambulances.size() / 2);
      A_minus.reserve(ambulances.size() / 2);

      vector<double> min_times(ambulances.size(), GRB_INFINITY);

      for (size_t a = 0; a < ambulances.size(); ++a) {
        auto &amb = ambulances[a];
        if (can_answer(amb, call) && amb.is_idle(time)) {
          A.push_back(a);
          double response_time = travel.get_response_time(amb, call, time);
          double pen_response_time =
              penalized_response_time(response_time, amb.type, call.priority);
          if (debug) {
            fmt::print("Amb {} response = {:.2f}, pen response = {:.2f}\n",
                       amb.id, response_time, pen_response_time);
          }
          min_times[a] = pen_response_time;
          if (response_time < T) {
            A_plus.push_back(a);
          } else {
            A_minus.push_back(a);
          }
        } else {
          if (debug) {
            fmt::print("Amb {} X\n", amb.id);
          }
        }
      }

      if (debug) {
        fmt::print("A = {}\n", A);
        fmt::print("A+ = {}\n", A_plus);
        fmt::print("A- = {}\n", A_minus);
      }

      for (int r = 0; r < ins.nb_regions; ++r) {
        k[r] = get_k(ins.centers[r], A);
      }
      if (debug) {
        for (int r = 0; r < ins.nb_regions; ++r) {
          if (k[r] > 0) {
            fmt::print("k[{}] = {} ", r, k[r]);
            if ((r + 1) % 7 == 0) {
              fmt::print("\n");
            }
          }
        }
      }

      int index_amb = -1;
      if (A_plus.size() > 0) {
        index_amb = get_index_amb(k, A_plus, lambda);
      } else if (A_minus.size() > 0) {
        index_amb = get_index_amb(k, A_minus, lambda);
      }
      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());

      if (index_amb >= 0) {
        auto &best_amb = ambulances[index_amb];
        double real_min_time = travel.get_response_time(best_amb, call, time);
        double waiting_on_scene_i = time + real_min_time - call.time;
        double waiting_to_hospital_i = best_amb.answer_call(
            call, travel, ins, time, real_min_time, nearest_base[call.id]);
        waiting_on_scene[call.id] = waiting_on_scene_i;
        waiting_on_scene_penalized[call.id] = penalized_response_time(
            waiting_on_scene_i, best_amb.type, call.priority);
        if (debug) {
          fmt::print("Call {} attended by {}\n", call.id, index_amb);
        }
        // fmt::print("| {:2f} {:.2f}\n", waiting_on_scene[call.id],
        // 	waiting_on_scene_penalized[call.id]);
        // fmt::print("Call location {}, amb location {}\n", call.location,
        // 	ambulances[index_amb].base_location);
        waiting_to_hospital[call.id] = waiting_to_hospital_i;
        which_ambulance[call.id] = best_amb.id;
        calls_end[call.id] = call.end;
        ++calls_attended;
      } else {
        if (debug) {
          fmt::print("Call {} back to queue\n", i);
        }
        queue_aux.push_back(i);
      }
    }
    queue.clear();
    for (auto i : queue_aux) {
      queue.push_back(i);
    }
    set_next_event(event_call, index_call);
    // if(debug){
    // 	cin.get();
    // }
  }
  if (debug) {
    print_results();
    // cin.get();
  }
}

int CoverageSolver::get_index_amb(vector<int> &k, vector<int> &A,
                                  vector<double> &lambda) {
  double min_sum = GRB_INFINITY;
  int min_ind = -1;

  for (auto a : A) {
    auto &amb = ambulances[a];
    Location curr_location;
    if (amb.arrival_time_at_b_last_trip <= time) {
      curr_location = amb.base_location;
    } else if (amb.arrival_time_at_f_last_trip <= time) {
      curr_location = travel.ambulance_position(amb, time);
    }
    double sum = 0;
    for (int r = 0; r < ins.nb_regions; ++r) {
      double time_to_r = travel.travel_time(curr_location, ins.centers[r], amb);
      if (time_to_r <= T) {
        sum += lambda[r] * (1 - q) * pow(q, k[r] - 1);
      }
    }
    if (sum < min_sum) {
      min_sum = sum;
      min_ind = a;
    }
  }
  return min_ind;
}

int CoverageSolver::get_k(Location &target, vector<int> &A) {
  int result = 0;
  for (auto a : A) {
    auto &amb = ambulances[a];
    Location amb_location = amb.get_current_location(ins, time);
    double travel_time = travel.travel_time(amb_location, target, amb);
    if (travel_time < T) {
      ++result;
    }
  }

  return result;
}

DistrictSolver::DistrictSolver(GRBEnv &env, vector<Call> &calls,
                               vector<Ambulance> &ambulances, Instance &ins,
                               Travel &travel, int g, DistrictData &a_data,
                               vector<double> &a_lambda, double time)
    : Solver(env, calls, ambulances, ins, travel, time),
      g(g),
      data(a_data),
      lambda(a_lambda) {
  prepare();
}

void DistrictSolver::run() {
  g_params.h_use_fixed_bases = true;
  // fmt::print("Districts final:\n");
  // print_districts(data);
  // closest_cross(district_data, lambda);

  for (size_t a = 0; a < ambulances.size(); ++a) {
    ambulances[a].base_location = ins.bases[data.ambulance_district[a]];
  }

  heuristic_cross();
  // print_results();
}

void DistrictSolver::closest_cross() {
  int event_call = 1;
  time = calls[0].time;
  index_call = 0;
  calls_attended = 0;
  travel.set_forward(false);

  auto &D = data.D;
  auto &amb_dists = data.ambulance_district;

  while (static_cast<size_t>(calls_attended) < calls.size()) {
    for (auto &amb : ambulances) {
      if (amb.arrival_time_at_b_last_trip <= time) {
        for (auto &call : calls) {
          call.ambulances.push_back(amb.id);
        }
      }
    }
    fmt::print(
        "Time = {}, index_call = {}, event_call = {}, queue.size() = {}\n",
        time, index_call, event_call, queue.size());
    if (event_call == 1) {
      queue.push_back(index_call);
    }
    vector<int> queue_aux;
    for (auto c : queue) {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto &call = calls[c];
      int call_dist = D[call.region];
      int index_amb_intra = -1;
      double min_time_intra = GRB_INFINITY;
      int index_amb_inter = -1;
      double min_time_inter = GRB_INFINITY;
      for (size_t i = 0; i < ambulances.size(); ++i) {
        auto &amb = ambulances[i];
        int amb_dist = amb_dists[i];
        double resp_time = penalized_response_time(
            travel.get_response_time(amb, call, time), amb.type, call.priority);
        if (call_dist == amb_dist && resp_time < min_time_intra) {
          min_time_intra = resp_time;
          index_amb_intra = i;
        }
        if (call_dist != amb_dist && resp_time < min_time_inter) {
          min_time_inter = resp_time;
          index_amb_inter = i;
        }
      }
      int index_amb = -1;
      double min_time = GRB_INFINITY;
      if (index_amb_intra >= 0) {
        index_amb = index_amb_intra;
      } else if (index_amb_inter >= 0) {
        index_amb = index_amb_inter;
      }

      if (index_amb >= 0) {
        auto &best_amb = ambulances[index_amb];
        double real_min_time = travel.get_response_time(best_amb, call, time);
        double waiting_on_scene_i = time + real_min_time - call.time;
        double waiting_to_hospital_i = best_amb.answer_call(
            call, travel, ins, time, real_min_time, nearest_base[call.id]);
        best_amb.arrival_time_at_f_last_trip =
            best_amb.arrival_time_at_b_last_trip;
        best_amb.free_location = best_amb.base_location;
        waiting_on_scene[call.id] = waiting_on_scene_i;
        waiting_on_scene_penalized[call.id] = penalized_response_time(
            waiting_on_scene_i, best_amb.type, call.priority);
        // fmt::print("Call {} attended by {}\n", call.id, index_amb);
        // fmt::print("| {:2f} {:.2f}\n", waiting_on_scene[call.id],
        // 	waiting_on_scene_penalized[call.id]);
        // fmt::print("Call location {}, amb location {}\n", call.location,
        // 	ambulances[index_amb].base_location);
        waiting_to_hospital[call.id] = waiting_to_hospital_i;
        which_ambulance[call.id] = best_amb.id;
        calls_end[call.id] = call.end;
        ++calls_attended;
      } else {
        queue_aux.push_back(c);
      }
      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
    }
    queue.clear();
    for (auto i : queue_aux) {
      queue.push_back(i);
    }
    set_next_event(event_call, index_call);
  }
}

void DistrictSolver::heuristic_cross() {
  calls_attended = 0;
  travel.set_forward(false);

  bool debug = false;

  auto &D = data.D;
  auto &amb_dists = data.ambulance_district;
  auto table_intra = get_ordered_table_intra();
  // fmt::print("Table intra {}\n",table_intra);
  auto table_inter = get_ordered_table_inter();
  // fmt::print("Table intra {}\n",table_inter);

  vector<int> call_dists(calls.size(), -1);
  for (auto &call : calls) {
    if (D[call.region] == -1) {
      int r = call.region;
      double min_time = GRB_INFINITY;
      int min_amb = -1;
      for (size_t i = 0; i < ambulances.size(); ++i) {
        double travel_time = travel.travel_time(
            ins.centers[r], ambulances[i].base_location, ambulances[i]);
        if (travel_time < min_time) {
          min_time = travel_time;
          min_amb = i;
        }
      }
      call_dists[call.id] = amb_dists[min_amb];
    } else {
      call_dists[call.id] = D[call.region];
    }
  }

  default_random_engine u_gen;
  uniform_real_distribution<double> u(0, 1);
  // travel.set_delay_param(u(u_gen));

  while (static_cast<size_t>(calls_attended) < calls.size()) {
    travel.set_delay_param(u(u_gen));
    if (event_call == 1) {
      queue.push_back(index_call);
    }
    if (debug) {
      fmt::print(
          "DistrictSolver: Time = {:.3f}, index = {}, event = {}, queue "
          "= {}, att = {}/{}\n",
          time, index_call, event_call, queue, calls_attended, calls.size());
    }
    vector<int> queue_aux;
    for (auto c : queue) {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto &call = calls[c];
      if (debug) {
        cout << call << "\n";
      }
      int call_dist = call_dists[c];
      int index_amb_intra = -1;
      double min_time_intra = GRB_INFINITY;
      int index_amb_inter = -1;
      double min_time_inter = GRB_INFINITY;

      if (call.priority == 0 ||
          call.priority == 2) {  // send closest ambulance intra then inter
        for (size_t i = 0; i < ambulances.size(); ++i) {
          auto &amb = ambulances[i];
          int amb_dist = amb_dists[i];
          double resp_time =
              penalized_response_time(travel.get_response_time(amb, call, time),
                                      amb.type, call.priority);
          if (call_dist == amb_dist && resp_time < min_time_intra &&
              can_answer(amb, call)) {
            min_time_intra = resp_time;
            index_amb_intra = i;
          }
        }

        if (index_amb_intra == -1) {
          for (size_t i = 0; i < ambulances.size(); ++i) {
            auto &amb = ambulances[i];
            int amb_dist = amb_dists[i];
            double resp_time = penalized_response_time(
                travel.get_response_time(amb, call, time), amb.type,
                call.priority);
            if (call_dist != amb_dist && resp_time < min_time_inter &&
                can_answer(amb, call)) {
              min_time_inter = resp_time;
              index_amb_inter = i;
            }
          }
        }

      } else if (call.priority == 1 ||
                 call.priority == 3) {  // Check table of ordered preferences
        for (auto a : table_intra[call_dist]) {
          auto &amb = ambulances[a];
          double resp_time =
              penalized_response_time(travel.get_response_time(amb, call, time),
                                      amb.type, call.priority);
          if (amb.arrival_time_at_b_last_trip <= time &&
              can_answer(amb, call)) {
            min_time_intra = resp_time;
            index_amb_intra = a;
            break;
          }
        }

        if (index_amb_intra == -1) {
          for (auto a : table_inter[call_dist]) {
            auto &amb = ambulances[a];
            double resp_time = penalized_response_time(
                travel.get_response_time(amb, call, time), amb.type,
                call.priority);
            if (amb.arrival_time_at_b_last_trip <= time &&
                can_answer(amb, call)) {
              min_time_inter = resp_time;
              index_amb_inter = a;
              break;
            }
          }
        }
      }

      int index_amb = -1;
      double min_time = GRB_INFINITY;
      if (index_amb_intra >= 0) {
        if (debug) {
          fmt::print("min_time intra\n");
        }
        index_amb = index_amb_intra;
      } else if (index_amb_inter >= 0) {
        if (debug) {
          fmt::print("min_time inter\n");
        }
        index_amb = index_amb_inter;
      }

      if (index_amb >= 0) {
        auto &best_amb = ambulances[index_amb];
        double real_min_time = travel.get_response_time(best_amb, call, time);
        double waiting_on_scene_i = time + real_min_time - call.time;
        double waiting_to_hospital_i = best_amb.answer_call(
            call, travel, ins, time, real_min_time, nearest_base[call.id]);
        best_amb.arrival_time_at_f_last_trip =
            best_amb.arrival_time_at_b_last_trip;
        best_amb.free_location = best_amb.base_location;
        waiting_on_scene[call.id] = waiting_on_scene_i;
        waiting_on_scene_penalized[call.id] = penalized_response_time(
            waiting_on_scene_i, best_amb.type, call.priority);
        if (debug) {
          fmt::print("Call {} attended by {} ", call.id, index_amb);
          fmt::print("| {:2f} {:.2f}\n", waiting_on_scene[call.id],
                     waiting_on_scene_penalized[call.id]);
          // fmt::print("Call location {}, amb location {}\n", call.location,
          // 	ambulances[index_amb].base_location);
        }
        waiting_to_hospital[call.id] = waiting_to_hospital_i;
        which_ambulance[call.id] = best_amb.id;
        calls_end[call.id] = call.end;
        ++calls_attended;
      } else {
        if (debug) {
          fmt::print("Call {} back to queue\n", c);
        }
        queue_aux.push_back(c);
      }
      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
    }
    queue.clear();
    for (auto i : queue_aux) {
      queue.push_back(i);
    }
    set_next_event(event_call, index_call);
    if (debug) {
      cin.get();
    }
  }
  if (debug) {
    print_results();
    cin.get();
  }
}

unordered_map<int, vector<int>> DistrictSolver::get_ordered_table_intra() {
  int n_dist = data.districts.size();
  unordered_map<int, vector<int>> intra_table;
  for (int d = 0; d < n_dist; ++d) {
    vector<vector<pair<double, int>>> table(data.districts[d].size(),
                                            vector<pair<double, int>>());
    vector<pair<double, int>> ordered_ambs;
    auto ambs_of_dist = data.ambulances_of_district(d);
    double sum_lambda = 0.0;
    int c = 0;
    for (auto r : data.districts[d]) {
      for (size_t i = 0; i < ambs_of_dist.size(); ++i) {
        auto &amb = ambulances[ambs_of_dist[i]];
        double travel_time =
            travel.travel_time(amb.base_location, ins.centers[r], amb);
        ordered_ambs.push_back(make_pair(travel_time, i));
      }
      sort(ordered_ambs.begin(), ordered_ambs.end());
      table[c++] = ordered_ambs;
      sum_lambda += lambda[r];
    }

    vector<double> z(data.districts[d].size(), 0);
    for (size_t i = 0; i < data.districts[d].size(); ++i) {
      int r = data.districts[d][i];
      z[i] = lambda[r] / sum_lambda;
    }
    int m = ambs_of_dist.size();
    int n = data.districts[d].size();
    vector<vector<double>> B(m, vector<double>(m, GRB_INFINITY));
    for (int j = 0; j < m; ++j) {
      for (int k = 0; k < m; ++k) {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
          if (table[i][j].second == k) {
            sum += z[i];
          }
        }
        B[k][j] = sum;
      }
    }

    vector<int> indexes(m, 0);
    for (int i = 0; i < m; ++i) {
      indexes[i] = i;
    }
    sort(indexes.begin(), indexes.end(), [&](const int i, const int j) {
      // i  is less than j if column i is less than column j
      for (int k = 0; k < m; ++k) {
        if (B[k][i] < B[k][j]) {
          return true;
        } else if (abs(B[k][i] - B[k][j]) < 0.0001) {
          continue;
        } else {
          return false;
        }
      }
      return false;
    });
    vector<int> amb_index(indexes.size(), 0);
    for (size_t i = 0; i < indexes.size(); ++i) {
      amb_index[i] = ambs_of_dist[indexes[i]];
    }
    intra_table[d] = amb_index;
  }
  return intra_table;
}

unordered_map<int, vector<int>> DistrictSolver::get_ordered_table_inter() {
  int n_dist = data.districts.size();
  unordered_map<int, vector<int>> inter_table;
  for (int d = 0; d < n_dist; ++d) {
    // get_number of ambulances outside district d
    vector<int> ambs_outside;
    for (size_t i = 0; i < ambulances.size(); ++i) {
      if (data.ambulance_district[i] != d) {
        ambs_outside.push_back(i);
      }
    }
    int m = ambs_outside.size();
    vector<pair<double, int>> ordered_ambs;
    vector<vector<pair<double, int>>> table(data.districts[d].size(),
                                            vector<pair<double, int>>());
    double sum_lambda = 0;
    int c = 0;
    for (auto r : data.districts[d]) {
      for (size_t i = 0; i < ambs_outside.size(); ++i) {
        auto &amb = ambulances[ambs_outside[i]];
        double travel_time =
            travel.travel_time(amb.base_location, ins.centers[r], amb);
        ordered_ambs.push_back(make_pair(travel_time, i));
      }
      sort(ordered_ambs.begin(), ordered_ambs.end());
      table[c++] = ordered_ambs;
      sum_lambda += lambda[r];
    }
    vector<double> z(data.districts[d].size(), 0);
    for (size_t i = 0; i < data.districts[d].size(); ++i) {
      int r = data.districts[d][i];
      z[i] = lambda[r] / sum_lambda;
    }

    int n = data.districts[d].size();
    vector<vector<double>> B(m, vector<double>(m, GRB_INFINITY));
    for (int j = 0; j < m; ++j) {
      for (int k = 0; k < m; ++k) {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
          if (table[i][j].second == k) {
            sum += z[i];
          }
        }
        B[k][j] = sum;
      }
    }

    vector<int> indexes(m, 0);
    for (int i = 0; i < m; ++i) {
      indexes[i] = i;
    }
    sort(indexes.begin(), indexes.end(), [&](const int i, const int j) {
      // i  is less than j if column i is less than column j
      for (int k = 0; k < m; ++k) {
        if (B[k][i] < B[k][j]) {
          return true;
        } else if (abs(B[k][i] - B[k][j]) < 0.0001) {
          continue;
        } else {
          return false;
        }
      }
      return false;
    });

    vector<int> amb_index(indexes.size(), 0);
    for (size_t i = 0; i < indexes.size(); ++i) {
      amb_index[i] = ambs_outside[indexes[i]];
    }
    inter_table[d] = amb_index;
  }

  return inter_table;
}

DistrictManager::DistrictManager(Instance &ins, Travel &travel, int g,
                                 vector<Ambulance> &ambulances)
    : ins(ins),
      travel(travel),
      g(g),
      tab_fac(ins.nb_ambulances + 1, 0),
      ambulances(ambulances) {
  lambda = vector<double>(ins.nb_regions, 0);
  for (int r = 0; r < ins.nb_regions; ++r) {
    double sum = 0;
    for (int t = 36; t < 40; ++t) {
      for (int p = 0; p < ins.nb_priorities; ++p) {
        sum += ins.lambda(t, g, r, p);
      }
    }
    lambda[r] = sum;
  }
}

DistrictData DistrictManager::get_districts() {
  DistrictData d = get_m_districts();
  auto &districts = d.districts;
  auto &D = d.D;
  int m = ins.nb_bases;
  auto W = get_workload(d);
  local_search(d, W);
  fmt::print("Finished swap\n");
  W = get_workload(d);
  fmt::print("W after swap: {}\n", W);
  DistrictData result = merge_step(d, W);
  fmt::print("Finished merge\n");
  return result;
}

DistrictData DistrictManager::get_m_districts() {
  std::default_random_engine gen(600);
  // poisson_distribution<int> dist(lambda[r]);
  unordered_map<int, vector<int>> districts;
  for (int i = 0; i < ins.nb_bases; ++i) {
    districts.insert(make_pair(i, vector<int>()));
  }
  vector<int> D(ins.nb_regions, -1);
  vector<int> ambulance_district(ambulances.size(), -1);
  // For each region with calls, compute the closest base
  for (int r = 0; r < ins.nb_regions; ++r) {
    poisson_distribution<int> dist(lambda[r]);
    int nb_calls = dist(gen);
    if (nb_calls > 0) {
      double min_d = GRB_INFINITY;
      int min_b = -1;
      for (auto &amb : ambulances) {
        double travel_time =
            travel.travel_time(amb.base_location, ins.centers[r], amb);
        if (travel_time < min_d) {
          min_d = travel_time;
          min_b = amb.id;
        }
      }
      districts[min_b].push_back({r});
      D[r] = min_b;
    }
  }
  std::vector<pair<double, int>> districts_by_lambda;
  if (ambulances.size() > ins.nb_bases) {
    for (int i = 0; i < ins.nb_bases; ++i) {
      ambulance_district[i] = i;
      double sum_lambda = 0.0;
      for (auto r : districts[i]) {
        sum_lambda += lambda[r];
      }

      districts_by_lambda.push_back(std::make_pair(sum_lambda, i));
    }
    std::sort(districts_by_lambda.begin(), districts_by_lambda.end(),
              std::greater<std::pair<double, int>>());
    for (size_t a = ins.nb_bases; a < ambulances.size(); ++a) {
      ambulance_district[a] =
          districts_by_lambda[a % districts_by_lambda.size()].second;
    }
  } else {
    for (int i = 0; i < ins.nb_bases; ++i) {
      ambulance_district[i] = i;
    }
  }

  return DistrictData{districts, D, ambulance_district};
}

unordered_map<int, double> DistrictManager::get_workload(
    DistrictData &district_data) {
  auto map_accum = [](double value,
                      const std::unordered_map<int, double>::value_type &p) {
    return value + p.second;
  };
  auto &districts = district_data.districts;
  auto &D = district_data.D;
  std::default_random_engine gen(600);
  // int n = ins.nb_regions;
  int m = districts.size();
  vector<int> dist_index(ins.nb_bases, 0);
  int ind = 0;
  vector<double> W(m, GRB_INFINITY);
  unordered_map<int, double> rho_til;
  unordered_map<int, double> rho_til_prev;
  for (auto &kv : districts) {
    int i = kv.first;
    auto &district = kv.second;
    double sum = 0;
    rho_til_prev[i] = rho_til[i] = get_mean_rate_of_district(district, i);
  }
  double r = accumulate(std::begin(rho_til), end(rho_til), 0.0, map_accum) / m;

  vector<vector<vector<int>>> G(m, vector<vector<int>>(m, vector<int>()));
  for (int k = 0; k < m; ++k) {
    for (int j = 0; j < m; ++j) {
      for (int r = 0; r < ins.nb_regions; ++r) {
        if (D[r] >= 0) {
          vector<pair<double, int>> amb_times;
          for (auto &amb : ambulances) {
            amb_times.push_back(make_pair(
                travel.travel_time(amb.base_location, ins.centers[r], amb),
                amb.id));
          }
          sort(amb_times.begin(), amb_times.end());
          if (amb_times[k].second == j) {
            G[k][j].push_back(r);
          }
        }
      }
    }
  }

  int n = 1;
  bool stop = false;
  while (!stop) {
    ++n;
    for (int j = 0; j < m; ++j) {
      double sum = 1;
      for (int k = 0; k < m; ++k) {
        for (auto g : G[k][j]) {
          double q = Q(rho_til_prev[j], j, m);
          double term = (k > 0) ? lambda[g] * q * pow(r, k) : lambda[g];
          sum += term;
        }
      }
      rho_til[j] = 1 - (1 / sum);
    }
    r = accumulate(rho_til.begin(), rho_til.end(), 0.0, map_accum) / m;
    double gamma = normalize(rho_til, r);
    for (size_t j = 0; j < rho_til.size(); ++j) {
      rho_til[j] = rho_til[j] / gamma;
    }
    double max_val = -1;
    for (size_t j = 0; j < rho_til.size(); ++j) {
      double val = abs(rho_til[j] - rho_til_prev[j]);
      if (val > max_val) {
        max_val = val;
      }
    }
    rho_til_prev = rho_til;
    if (max_val < g_params.EPS) {
      stop = true;
    }
  }
  return rho_til;
}

double DistrictManager::normalize(unordered_map<int, double> &rho_til,
                                  double r) {
  return accumulate(rho_til.begin(), rho_til.end(), 0.0,
                    [](double value,
                       const std::unordered_map<int, double>::value_type &p) {
                      return value + p.second;
                    }) /
         r * rho_til.size();
}

unsigned long int DistrictManager::fac(int n) {
  if (n > ins.nb_ambulances) {
    fmt::print("Error  calling factorial of {} > {} = nb of ambulances\n", n,
               ins.nb_ambulances);
    exit(1);
  }
  if (tab_fac[n] > 0) {
    return tab_fac[n];
  }
  unsigned long int prod = 1;
  for (int i = n; i >= 1; --i) {
    prod *= i;
  }

  tab_fac[n] = prod;
  return tab_fac[n];
}

double DistrictManager::Q(double rho, int j, int N) {
  double sum_upper = 0;
  auto fac_N = fac(N);
  for (int k = j; k < N; ++k) {
    sum_upper += (fac(N - k - 1) * (N - k) / fac(k - j)) * (pow(N, k) / fac_N) *
                 pow(rho, k - j);
  }
  double sum_lower = 0;
  for (int i = 0; i < N; ++i) {
    sum_lower += (pow(N, i) / fac(i)) * pow(rho, i);
  }
  double result =
      sum_upper / ((1 - rho) * sum_lower + pow(N, N) * pow(rho, N) / fac_N);
  return result;
  // return pow(rho,j)*(1-rho);
}

double DistrictManager::var(unordered_map<int, double> &rho) {
  double mean =
      accumulate(rho.begin(), rho.end(), 0.0,
                 [](double value,
                    const std::unordered_map<int, double>::value_type &p) {
                   return value + p.second;
                 }) /
      rho.size();
  double sum = 0;
  for (size_t i = 0; i < rho.size(); ++i) {
    sum += pow(rho[i] - mean, 2);
  }

  return sum / rho.size();
}

void DistrictManager::local_search(DistrictData &district_data,
                                   unordered_map<int, double> &W) {
  auto &districts = district_data.districts;
  auto &D = district_data.D;
  double var_w = var(W);
  int m = ins.nb_bases;
  vector<int> region_of_district(ins.nb_bases, -1);
  for (int b = 0; b < ins.nb_bases; ++b) {
    region_of_district[b] = ins.nearest_region_to_base[b];
  }
  for (int k = 0; k < m; ++k) {
    int rk = region_of_district[k];
    for (int r = 0; r < ins.nb_regions; ++r) {
      if (D[r] == k) {
        for (auto s : ins.neighbors[r]) {
          int ds = D[s];
          int rs = region_of_district[ds];
          if (r != s && D[s] >= 0 && D[s] != k && s != rk) {
            // do a swap
            swap(district_data, k, ds, r, s);
            auto W2 = get_workload(district_data);
            double var_w2 = var(W2);
            if (var_w2 < var_w) {
              return;
            } else {
              swap(district_data, k, ds, r, s);
            }
          }
        }
      }
    }
  }
}

void DistrictManager::swap(DistrictData &district_data, int dr, int ds, int r,
                           int s) {
  auto &districts = district_data.districts;
  auto &D = district_data.D;
  districts[dr].erase(remove(districts[dr].begin(), districts[dr].end(), r),
                      districts[dr].end());
  districts[dr].push_back(s);
  D[s] = dr;
  districts[ds].erase(remove(districts[ds].begin(), districts[ds].end(), s),
                      districts[ds].end());
  districts[ds].push_back(r);
  D[r] = ds;
}

double DistrictManager::total_aexc(DistrictData &district_data,
                                   unordered_map<int, double> &rho) {
  auto &districts = district_data.districts;
  auto &D = district_data.D;
  double result = 0;
  for (size_t i = 0; i < districts.size(); ++i) {
    result += aexc(district_data, rho, i);
  }
  return result;
}

double DistrictManager::aexc(DistrictData &district_data,
                             unordered_map<int, double> &rho, int i) {
  int n = ins.nb_regions;
  auto &D = district_data.D;
  int m = district_data.districts.size();
  double rho_mean =
      accumulate(rho.begin(), rho.end(), 0.0,
                 [](double value,
                    const std::unordered_map<int, double>::value_type &p) {
                   return value + p.second;
                 }) /
      rho.size();

  vector<int> ambs_service(n, 0);

  double time_standard = 1200;
  for (int r = 0; r < n; ++r) {
    int count_r = 0;
    for (auto &amb : ambulances) {
      double travel_time =
          travel.travel_time(amb.base_location, ins.centers[r], amb);
      count_r += 1 * (travel_time <= time_standard);
    }
    ambs_service[r] = count_r;
  }

  double result = 0;
  default_random_engine gen;
  for (int r = 0; r < n; ++r) {
    if (D[r] == i) {
      int num_servers = ambs_service[r];
      for (int j = 1; j <= m - 1; ++j) {
        int y = (num_servers >= j + 1);
        poisson_distribution<int> nc(lambda[i]);
        int h = nc(gen);
        // double h = lambda[i];
        result += (1 - rho_mean) * pow(rho_mean, j) * h * y * Q(rho_mean, j, m);
      }
    }
  }
  return result;
}

double DistrictManager::get_mean_rate_of_district(vector<int> &district,
                                                  int b) {
  double result = 0;
  Location base = ins.bases[b];
  double time_on_scene = 600;
  double time_at_hospital = 600;

  for (int r : district) {
    Location h = ins.hospitals[ins.nearest_hospital_to_region[r]];
    double sum_service_time = 0;
    for (auto &s : ins.samples[r]) {
      sum_service_time +=
          travel.region_service_time(base, s, h, ambulances[0].speed);
      sum_service_time += time_on_scene + time_at_hospital;
    }
    double mean_service_time = sum_service_time / ins.samples[r].size();
    result += lambda[r] * mean_service_time;
  }

  return result;
}

DistrictData DistrictManager::merge_step(DistrictData &district_data,
                                         unordered_map<int, double> &rho) {
  double old_aexc = total_aexc(district_data, rho);
  DistrictData old_districts = district_data;

  bool stop = false;
  DistrictData result = district_data;
  int merge_count = 0;
  do {
    auto &districts = result.districts;
    vector<pair<int, int>> P;

    for (auto kvi : districts) {
      int i = kvi.first;
      for (auto kvj : districts) {
        int j = kvj.first;
        if (i != j) {
          if (can_merge(districts, i, j)) {
            P.push_back(make_pair(i, j));
          }
        }
      }
    }
    auto this_rho = get_workload(result);
    double this_aexc = total_aexc(result, this_rho);
    double max_aexc = this_aexc;
    pair<int, int> max_merge{-1, -1};
    for (auto p : P) {
      int i = p.first;
      int j = p.second;
      fmt::print("merge_count = {} | merge\n", merge_count);
      DistrictData merged = merge(result, i, j);
      fmt::print("merge_count = {} | get_workload\n", merge_count);
      auto merged_rho = get_workload(merged);
      fmt::print("merge_count = {} | total_aexc\n", merge_count);
      double merged_aexc = total_aexc(merged, merged_rho);
      // if(merged_aexc - this_aexc > 0){
      // 	fmt::print("Pair {}, gain = {} (abs = {})\n", p, (merged_aexc -
      // this_aexc)*100/this_aexc, merged_aexc - this_aexc);
      // }
      if (merged_aexc > max_aexc) {
        max_aexc = merged_aexc;
        max_merge = make_pair(i, j);
      }
    }
    if (max_merge != make_pair(-1, -1)) {
      result = merge(result, max_merge.first, max_merge.second);
      ++merge_count;
    } else {
      stop = true;
    }
  } while (!stop);
  // fmt::print("Merge Count = {}\n", merge_count);
  return result;
}

DistrictData DistrictManager::merge(DistrictData &district_data, int di,
                                    int dj) {
  DistrictData result = district_data;
  auto &districts = result.districts;
  auto &D = result.D;

  // insert dj in di
  result.districts[di].insert(districts[di].end(), districts[dj].begin(),
                              districts[dj].end());
  // sort di
  sort(districts[di].begin(), districts[di].end());
  // erase dj
  result.districts.erase(dj);
  // Set dj entries of D vector to di
  for (int i = 0; i < ins.nb_regions; ++i) {
    if (D[i] == dj) {
      D[i] = di;
    }
  }
  // Set dj entries of ambulance_district vector to di
  for (size_t i = 0; i < ambulances.size(); ++i) {
    if (result.ambulance_district[i] == dj) {
      result.ambulance_district[i] = di;
    }
  }
  return result;
}

bool DistrictManager::can_merge(unordered_map<int, vector<int>> &districts,
                                int i, int j) {
  if (districts[i].size() == 0 || districts[j].size() == 0) {
    return true;
  }

  for (auto r : districts[i]) {
    for (auto s : districts[j]) {
      if (ins.adj_matrix[r][s]) {
        return true;
      }
    }
  }

  return false;
}

ForwardPrepSolver::ForwardPrepSolver(GRBEnv &env, vector<Call> &calls,
                                     vector<Ambulance> &ambulances,
                                     Instance &ins, Travel &travel, int g,
                                     double time)
    : Solver(env, calls, ambulances, ins, travel, time), g(g) {}

void ForwardPrepSolver::run() {
  travel.set_forward(true);

  vector<vector<double>> lambda(ins.nb_times, vector<double>());
  for (int t = 0; t < ins.nb_times; ++t) {
    lambda[t] = get_lambda(g, t);
  }

  const double time_slot = 1800;
  for (size_t i = 0; i < calls.size(); ++i) {
    auto t0 = std::chrono::high_resolution_clock::now();
    double min_time = GRB_INFINITY;
    int index_amb = -1;
    auto &call = calls[i];
    time = call.time;

    for (auto &amb : ambulances) {
      if (can_answer(amb, call)) {
        // double time_amb = travel.get_response_time(amb, call, time) *
        // ins.penalty_matrix[amb.type][call.priority];
        double time_amb = penalized_response_time(
            travel.get_response_time(amb, call, time), amb.type, call.priority);
        // std::cout << amb << " " << time_amb << "\n";
        if (time_amb < min_time) {
          min_time = time_amb;
          index_amb = amb.id;
        } else if (abs(time_amb - min_time) < g_params.EPS && index_amb != -1 &&
                   amb.type > ambulances[index_amb].type) {
          min_time = time_amb;
          index_amb = amb.id;
        }
      }
    }

    if (index_amb >= 0) {
      auto &best_amb = ambulances[index_amb];
      double waiting_on_scene_i = time + min_time - call.time;
      auto finish_service_data = get_finish_service_data(best_amb, call, time);
      int t_begin = finish_service_data.free_time / 1800;
      if (t_begin >= ins.nb_times) {
        t_begin = ins.nb_times - 1;
      }
      int base_return = get_return_base(index_amb, finish_service_data,
                                        lambda[t_begin], travel.forward);
      double waiting_to_hospital_i = ambulances[index_amb].answer_call(
          call, travel, ins, time, min_time, base_return);

      waiting_on_scene[i] = waiting_on_scene_i;
      waiting_on_scene_penalized[i] = penalized_response_time(
          waiting_on_scene_i, best_amb.type, call.priority);
      waiting_to_hospital[i] = waiting_to_hospital_i;
      which_ambulance[i] = index_amb;
      calls_end[i] = call.end;
    } else {
      waiting_on_scene[i] = GRB_INFINITY;
      waiting_on_scene_penalized[i] = GRB_INFINITY;
      waiting_to_hospital[i] = GRB_INFINITY;
      calls_end[i] = GRB_INFINITY;
      which_ambulance[i] = -1;
    }
    auto dt = std::chrono::high_resolution_clock::now();
    run_times.push_back(
        std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
  }
}

MinMaxPrepSolver::MinMaxPrepSolver(GRBEnv &env, vector<Call> &calls,
                                   vector<Ambulance> &ambulances, Instance &ins,
                                   Travel &travel, int g, double time)
    : Solver(env, calls, ambulances, ins, travel, time), g(g) {
  queue.reserve(calls.size());
}

void MinMaxPrepSolver::run() {
  event_call = 1;
  time = calls[0].time;
  index_call = 0;
  calls_attended = 0;
  travel.set_forward(true);

  vector<vector<double>> lambda(ins.nb_times, vector<double>());
  for (int t = 0; t < ins.nb_times; ++t) {
    lambda[t] = get_lambda(g, t);
  }

  while (static_cast<size_t>(calls_attended) < calls.size()) {
    if (event_call == 1) {
      queue.push_back(index_call);
    }

    std::vector<double> min_times(queue.size(), GRB_INFINITY);
    std::vector<std::vector<double>> travel_times(
        queue.size(), std::vector<double>(ambulances.size(), GRB_INFINITY));
    std::vector<int> index_ambs(queue.size(), -1);
    // fmt::print("time = {}, index_call = {}, event = {}, queue = {}\n",time,
    // index_call, event_call, queue); For every call in queue, determine the
    // best ambulance and the best response time
    for (size_t i = 0; i < queue.size(); ++i) {
      auto &call = calls[queue[i]];
      for (auto &amb : ambulances) {
        int j = amb.id;
        if (can_answer(amb, call)) {
          // travel_times[i][j] = travel.get_response_time(amb, call, time) *
          // ins.penalty_matrix[amb.type][call.priority];
          travel_times[i][j] =
              penalized_response_time(travel.get_response_time(amb, call, time),
                                      amb.type, call.priority);

          if (travel_times[i][j] < min_times[i]) {
            min_times[i] = travel_times[i][j];
            index_ambs[i] = j;
          } else if (abs(travel_times[i][j] - min_times[i]) < g_params.EPS &&
                     index_ambs[i] != -1 &&
                     amb.type > ambulances[index_ambs[i]].type) {
            min_times[i] = travel_times[i][j];
            index_ambs[i] = j;
          }
        }
      }
      min_times[i] += time - call.time;
    }

    std::vector<int> remaining_indexes(queue.size(), 0);
    for (size_t i = 0; i < queue.size(); ++i) {
      remaining_indexes[i] = i;
    }

    std::vector<int> queue_aux;
    int nb_treated = 0;
    int total_queue = queue.size();

    while (nb_treated < total_queue) {
      auto t0 = std::chrono::high_resolution_clock::now();
      // get time and index of the worst call, and the best ambulance to such
      // call
      auto max_it = std::max_element(min_times.begin(), min_times.end());
      double max_time = *max_it;
      int current_call = std::distance(min_times.begin(), max_it);
      int call_ind = queue[remaining_indexes[current_call]];
      auto &call = calls[call_ind];
      int index_amb = index_ambs[remaining_indexes[current_call]];
      auto &best_amb = ambulances[index_amb];

      // if best ambulance is free, it should answer the call
      if (index_amb >= 0 && best_amb.arrival_time_at_f_last_trip <= time) {
        Ambulance &best_amb = ambulances[index_amb];
        double waiting_on_scene_i = max_time;
        FinishServiceData finish_service_data(
            get_finish_service_data(best_amb, call, time));
        // fmt::print("\tfinish: time {}, location {}\n",
        // finish_service_data.free_time, finish_service_data.free_location);
        int t_begin = finish_service_data.free_time / 1800;
        if (t_begin >= ins.nb_times) {
          t_begin = ins.nb_times - 1;
        }
        int base_return = get_return_base(index_amb, finish_service_data,
                                          lambda[t_begin], true);
        // fmt::print("\tbase_return = {}, nearest_base = {}\n", base_return,
        // nearest_base[call.id]);
        double waiting_to_hospital_i =
            best_amb.answer_call(call, travel, ins, time,
                                 max_time - (time - call.time), base_return);
        waiting_on_scene[call_ind] = waiting_on_scene_i;
        waiting_on_scene_penalized[call_ind] = penalized_response_time(
            waiting_on_scene_i, best_amb.type, call.priority);
        waiting_to_hospital[call_ind] = waiting_to_hospital_i;
        which_ambulance[call_ind] = index_amb;
        calls_end[call_ind] = call.end;

        calls_attended++;
        // update travel times
        for (int k = 0; k < total_queue - nb_treated; ++k) {
          int ind = remaining_indexes[k];
          if (k != current_call && index_amb == index_ambs[ind]) {
            double time_to_h = best_amb.arrival_time_at_f_last_trip - time;
            double time_from_h_to_c = travel.travel_time(
                best_amb.free_location, calls[queue[ind]].location, best_amb);
            travel_times[ind][best_amb.id] = time_to_h + time_from_h_to_c;
            auto min_it = std::min_element(travel_times[ind].begin(),
                                           travel_times[ind].end());
            double min_val = *min_it;
            int min_ind = std::distance(travel_times[ind].begin(), min_it);
            min_times[k] = min_val + time - calls[queue[ind]].time;
            index_ambs[ind] = min_ind;
          }
        }
      } else {
        queue_aux.push_back(queue[remaining_indexes[current_call]]);
      }
      if (nb_treated + 1 < total_queue) {
        min_times.erase(max_it);
        remaining_indexes.erase(
            std::remove(remaining_indexes.begin(), remaining_indexes.end(),
                        remaining_indexes[current_call]),
            remaining_indexes.end());
      }
      nb_treated++;
      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
    }

    queue.clear();
    for (auto i : queue_aux) {
      queue.push_back(i);
    }

    set_next_event(event_call, index_call);
  }
}

PriorityPrepSolver::PriorityPrepSolver(GRBEnv &env, vector<Call> &calls,
                                       vector<Ambulance> &ambulances,
                                       Instance &ins, Travel &travel, int g,
                                       double time)
    : Solver(env, calls, ambulances, ins, travel, time), g(g) {
  queue.reserve(calls.size());
}

void PriorityPrepSolver::run() {
  event_call = 1;
  time = calls[0].time;
  index_call = 0;
  calls_attended = 0;

  travel.set_forward(false);

  vector<vector<double>> lambda(ins.nb_times, vector<double>());
  for (int t = 0; t < ins.nb_times; ++t) {
    lambda[t] = get_lambda(g, t);
  }

  while (static_cast<size_t>(calls_attended) < calls.size()) {
    if (event_call == 1) {
      queue.push_back(index_call);
    }

    std::vector<std::pair<double, int>> sorted_queue;
    if (queue.size() > 1) {
      sorted_queue.reserve(queue.size());
      for (size_t i = 0; i < queue.size(); ++i) {
        auto &call = calls[queue[i]];
        double pen_resp_time =
            penalized_response_time((time - call.time), -1, call.priority);
        sorted_queue.push_back(std::make_pair(pen_resp_time, queue[i]));
      }
      std::sort(sorted_queue.begin(), sorted_queue.end(),
                std::greater<std::pair<double, int>>());
    } else if (queue.size() == 1) {
      auto &call = calls[queue[0]];
      double pen_resp_time =
          penalized_response_time((time - call.time), -1, call.priority);
      sorted_queue.push_back(std::make_pair(pen_resp_time, queue[0]));
    }

    std::vector<int> queue_aux;
    for (size_t k = 0; k < sorted_queue.size(); ++k) {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto call_ind = sorted_queue[k].second;
      auto &call = calls[call_ind];
      double min_time = GRB_INFINITY;
      int index_amb = -1;

      for (auto &amb : ambulances) {
        if (can_answer(amb, call)) {
          // double time_amb = travel.get_response_time(amb,call,time) *
          // ins.penalty_matrix[amb.type][call.priority];
          double time_amb =
              penalized_response_time(travel.get_response_time(amb, call, time),
                                      amb.type, call.priority);

          if (time_amb < min_time) {
            min_time = time_amb;
            index_amb = amb.id;
          } else if (abs(time_amb - min_time) < g_params.EPS &&
                     index_amb != -1 && amb.type > ambulances[index_amb].type) {
            min_time = time_amb;
            index_amb = amb.id;
          }
        }
      }

      if (index_amb >= 0 &&
          ambulances[index_amb].arrival_time_at_f_last_trip <= time) {
        auto &best_amb = ambulances[index_amb];
        double waiting_on_scene_i = time + min_time - call.time;
        FinishServiceData finish_service_data(
            get_finish_service_data(best_amb, call, time));
        int t_begin = finish_service_data.free_time / 1800;
        if (t_begin >= ins.nb_times) {
          t_begin = ins.nb_times - 1;
        }
        int base_return = get_return_base(index_amb, finish_service_data,
                                          lambda[t_begin], true);
        double waiting_to_hospital_i = ambulances[index_amb].answer_call(
            call, travel, ins, time, min_time, base_return);
        waiting_on_scene[call_ind] = waiting_on_scene_i;
        waiting_on_scene_penalized[call_ind] = penalized_response_time(
            waiting_on_scene_i, best_amb.type, call.priority);
        waiting_to_hospital[call_ind] = waiting_to_hospital_i;
        which_ambulance[call_ind] = index_amb;
        calls_end[call_ind] = call.end;
        calls_attended++;
      } else {
        queue_aux.push_back(sorted_queue[k].second);
      }
      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
    }

    queue.clear();
    for (auto i : queue_aux) {
      queue.push_back(i);
    }
    set_next_event(event_call, index_call);
  }
}

FinishServiceData Solver::get_finish_service_data(Ambulance &amb, Call &call,
                                                  double time) {
  // double min_time = travel.get_response_time(amb,call, time) *
  // ins.penalty_matrix[amb.type][call.priority];
  double min_time = penalized_response_time(
      travel.get_response_time(amb, call, time), amb.type, call.priority);
  Location free_location = call.location;
  double free_time = time + min_time + call.time_on_scene;
  if (call.hosp_needed && call.clean_needed) {
    free_location = ins.cleaning_bases[call.cleaning];
    free_time +=
        travel.travel_time(call.location, ins.hospitals[call.hospital], amb) +
        call.time_at_hospital +
        travel.travel_time(ins.hospitals[call.hospital],
                           ins.cleaning_bases[call.cleaning], amb) +
        call.cleaning_time;
  } else if (call.hosp_needed) {
    free_location = ins.hospitals[call.hospital];
    free_time +=
        travel.travel_time(call.location, ins.hospitals[call.hospital], amb) +
        call.time_at_hospital;
  } else if (call.clean_needed) {
    free_location = ins.cleaning_bases[call.cleaning];
    free_time += travel.travel_time(call.location,
                                    ins.cleaning_bases[call.cleaning], amb) +
                 call.cleaning_time;
  }
  return {free_time, free_location};
}

int Solver::get_return_base(int amb_id, FinishServiceData &finish_service_data,
                            vector<double> &lambda, bool forward) {
  int result = -1;
  assert(amb_id >= 0 || "Invalid ambulance id in get_return_base function");
  double free_time = finish_service_data.free_time;
  Location &free_location = finish_service_data.free_location;
  const int max_size = 34;
  auto &amb = ambulances[amb_id];
  auto time_to_base = time_to_bases(free_location);
  sort(time_to_base.begin(), time_to_base.end());
  auto prep_to_base =
      prep_to_bases(time_to_base, amb_id, max_size, free_time, lambda, forward);
  sort(prep_to_base.begin(), prep_to_base.end(), greater<pair<double, int>>());
  if (prep_to_base[0].first > g_params.min_preparedness) {
    double max_prep = prep_to_base[0].first;
    int best_b = prep_to_base[0].second;
    double min_time =
        travel.travel_time(amb.free_location, ins.bases[best_b], amb);
    for (size_t i = 1; i < prep_to_base.size(); ++i) {
      auto prep = prep_to_base[i].first;
      int b = prep_to_base[i].second;
      if (prep < max_prep) {
        break;
      } else {
        double travel_time =
            travel.travel_time(amb.free_location, ins.bases[b], amb);
        if (travel_time < min_time) {
          min_time = travel_time;
          best_b = b;
        }
      }
    }
    result = best_b;
  } else {
    result = prep_to_base[0].second;
  }
  return result;
}

vector<double> Solver::get_lambda(int g, int t_begin) {
  vector<double> lambda(ins.nb_regions, 0.0);
  // fmt::print("T begin = {}\n", t_begin);
  for (int r = 0; r < ins.nb_regions; ++r) {
    double sum = 0.0;
    for (int t = t_begin; t < ins.nb_times; ++t) {
      for (int p = 0; p < ins.nb_priorities; ++p) {
        sum += ins.lambda(t, g, r, p);
      }
    }
    lambda[r] = sum;
  }
  return lambda;
}

vector<double> Solver::prep(vector<double> &lambda) {
  vector<double> result(ins.nb_regions, 0.0);
  for (int r = 0; r < ins.nb_regions; ++r) {
    double sum = 0;
    for (size_t i = 0; i < ambulances.size(); ++i) {
      auto &other_amb = ambulances[i];
      double pw = p_weight[other_amb.type];
      Location curr_location = null_location;
      double travel_time = 0.0;
      if (other_amb.arrival_time_at_b_last_trip <= time) {
        curr_location = other_amb.base_location;
      } else if (other_amb.arrival_time_at_f_last_trip <= time) {
        curr_location = travel.ambulance_position(other_amb, time);
      } else {
        curr_location = other_amb.free_location;
        travel_time = other_amb.arrival_time_at_f_last_trip - time;
      }
      travel_time +=
          travel.travel_time(curr_location, ins.centers[r], other_amb);
      sum += pw / travel_time;
    }
    double lam = (lambda[r] <= g_params.EPS) ? g_params.EPS : lambda[r];
    result[r] = (1 / lam) * sum;
  }

  return result;
}

vector<pair<double, int>> Solver::time_to_bases(Location &free_location) {
  vector<pair<double, int>> result;
  result.reserve(ins.nb_bases);
  for (int b = 0; b < ins.nb_bases; ++b) {
    result.push_back(make_pair(
        travel.travel_time(free_location, ins.bases[b], ambulances[0]), b));
  }

  return result;
}

SortableVector Solver::prep_to_bases(SortableVector &time_to_base, int amb_id,
                                     size_t max_size, double free_time,
                                     vector<double> &lambda, bool forward) {
  const double threshold = 300;
  Ambulance &amb = ambulances[amb_id];
  size_t k = 0;
  vector<pair<double, int>> prep_to_base;
  prep_to_base.reserve(max_size);
  while (k < max_size && k < time_to_base.size()) {
    int b = time_to_base[k].second;
    double min_prep = GRB_INFINITY;
    int min_r = -1;
    for (int r = 0; r < ins.nb_regions; ++r) {
      double sum = 0;
      for (size_t i = 0; i < ambulances.size(); ++i) {
        Ambulance &other_amb = ambulances[i];
        double pw = p_weight[other_amb.type];
        // double pw = 1;
        Location curr_location = null_location;
        double travel_time = 0;
        if (i == static_cast<size_t>(amb_id)) {
          curr_location = ins.bases[b];
        } else {
          if (other_amb.arrival_time_at_b_last_trip <= time) {
            curr_location = other_amb.base_location;
          } else if (other_amb.arrival_time_at_f_last_trip <= time) {
            curr_location = travel.ambulance_position(other_amb, time);
          } else if (forward) {
            travel_time += other_amb.arrival_time_at_f_last_trip - time;
            curr_location = other_amb.free_location;
          }
        }
        if (curr_location != null_location) {
          travel_time +=
              travel.travel_time(curr_location, ins.centers[r], other_amb);
          sum += pw / travel_time;
        }
      }
      double lam = (lambda[r] < g_params.EPS) ? g_params.EPS : lambda[r];
      double prep = (1 / lambda[r]) * sum;
      if (prep < min_prep) {
        min_prep = prep;
      }
    }
    prep_to_base.push_back(make_pair(min_prep, b));
    ++k;
  }

  return prep_to_base;
}

EnumerateSolver::EnumerateSolver(GRBEnv &env, vector<Call> &calls,
                                 vector<Ambulance> &ambulances, Instance &ins,
                                 Travel &travel, int g, double time)
    : Solver(env, calls, ambulances, ins, travel, time),
      model(env),
      y(calls.size()),
      routes_by_ambulance(ambulances.size(), set<int>()),
      routes_by_call(calls.size(), set<int>()) {
  x.reserve(ambulances.size() * 40);
  routes.reserve(ambulances.size() * 40);
  prepare();
}

void EnumerateSolver::load_model() {
  int index_route = 0;
  for (auto &route : routes) {
    x.push_back(model.addVar(0, 1, route.cost, GRB_BINARY,
                             fmt::format("x_{}", index_route++)));
  }
  const double penalty_factor = 20000;
  for (auto &call : calls) {
    double penalty =
        (call.priority % 2 == 0) ? 3 * penalty_factor : 2 * penalty_factor;
    y[call.id] =
        model.addVar(0, 1, penalty, GRB_BINARY, fmt::format("y_{}", call.id));
  }
  for (auto &amb : ambulances) {
    GRBLinExpr con = 0;
    for (auto r : routes_by_ambulance[amb.id]) {
      con += x[r];
    }
    model.addConstr(con, GRB_LESS_EQUAL, 1,
                    fmt::format("routes_amb_{}", amb.id));
  }

  for (auto &call : calls) {
    GRBLinExpr con = y[call.id];
    for (auto r : routes_by_call[call.id]) {
      con += x[r];
    }
    model.addConstr(con, GRB_EQUAL, 1, fmt::format("routes_call_{}", call.id));
  }

  model.update();
  model.set(GRB_IntParam_OutputFlag, 0);
}

vector<Route> EnumerateSolver::generate_routes(Ambulance &amb) {
  vector<Route> generated_routes;
  const int max_depth = 3;  // TODO: CHANGE
  double free_time = max(time, amb.arrival_time_at_f_last_trip);

  shared_ptr<EnumNode> root(nullptr);
  for (size_t i = 0; i < calls.size(); ++i) {
    auto &call = calls[i];
    if ((call.time >= time || is_call_on_queue[i]) && (can_answer(amb, call))) {
      insert_enum_node(amb, call, root);
    }
  }
  while (!node_queue.empty()) {
    EnumNode &current_node = node_queue.front();
    auto &current_call = current_node.call;
    auto route = expand_route(amb, current_node);
    if (debug) {
      string route_ids = "";
      for (auto &node : route.node_sequence) {
        route_ids += fmt::format("{} ", node.node_index);
      }
      fmt::print(
          "Extending node: amb {}, call {}, depth = {}, current route {}:\n",
          amb.id, current_call.id, current_node.depth, route_ids);
    }
    generated_routes.push_back(route);
    if (current_node.depth < max_depth - 1) {
      for (size_t i = 0; i < calls.size(); ++i) {
        auto &call = calls[i];
        if (call.time > current_call.time) {
          if (debug) {
            cout << "\t" << call << "\n";
          }
          insert_enum_node(amb, call, make_shared<EnumNode>(current_node));
        }
      }
    }
    node_queue.pop();
  }

  return generated_routes;
}

void EnumerateSolver::insert_enum_node(Ambulance &amb, Call &call,
                                       shared_ptr<EnumNode> parent) {
  if (!can_answer(amb, call) && amb.arrival_time_at_b_last_trip <= call.time) {
    return;
  }
  // bool debug = parent == nullptr && (amb.id == 2 || amb.id == 14) && call.id
  // == 2;
  bool debug = false;
  Location origin;
  double origin_time;
  if (parent != nullptr) {
    auto &former_call = parent->call;
    double free_time = parent->time_finish_service;
    if (debug) {
      fmt::print("Parent waiting time = {}\n", parent->waiting_time);
      fmt::print("Finish parent time  = {}\n", free_time);
    }
    Location free_location = ins.hospitals[former_call.hospital];
    Location nearest_base =
        (g_params.h_use_fixed_bases)
            ? amb.base_location
            : ins.bases[ins.nearest_base_to_hospital[former_call.hospital]];
    origin_time = GRB_INFINITY;
    if (free_time > call.time) {
      if (debug) {
        fmt::print("free_time > call.time\n");
      }
      origin = free_location;
      origin_time = free_time;
    } else {
      if (debug) {
        fmt::print("free_time <= call.time\n");
      }
      double ttime = travel.travel_time(free_location, nearest_base, amb);
      origin = travel.position_between_origin_destination(
          free_location, nearest_base, free_time, call.time, free_time + ttime,
          amb);
      origin_time = call.time;
    }
  } else {
    // get ambulance current_location
    Location location;
    origin_time = GRB_INFINITY;
    if (amb.arrival_time_at_b_last_trip <= time) {
      origin = amb.base_location;
      origin_time = call.time;
      if (debug) {
        fmt::print("Origin 1\n");
      }
    } else if (amb.arrival_time_at_f_last_trip <= time) {
      location = travel.ambulance_position(amb, time);
      if (call.time > time) {
        double ttime = travel.travel_time(location, amb.base_location, amb);
        origin = travel.position_between_origin_destination(
            location, amb.base_location, time, call.time, time + ttime, amb);
        origin_time = call.time;
      } else {
        origin = location;
        origin_time = time;
      }
      if (debug) {
        fmt::print("Current location at time {}: {}\n", time, location);
        fmt::print("Where amb will be at call.time {}: {}\n", call.time,
                   origin);
      }
      if (debug) {
        fmt::print("Origin 2\n");
      }
    } else {
      location = amb.free_location;
      if (amb.arrival_time_at_f_last_trip >= call.time) {
        origin = amb.free_location;
        origin_time = amb.arrival_time_at_f_last_trip;
        if (debug) {
          fmt::print("Origin 3\n");
        }
      } else {
        origin = travel.position_between_origin_destination(
            amb.free_location, amb.base_location,
            amb.arrival_time_at_f_last_trip, call.time,
            amb.arrival_time_at_b_last_trip, amb);
        origin_time = call.time;
        if (debug) {
          fmt::print("Origin 4\n");
        }
      }
    }
  }

  double response_time = travel.travel_time(origin, call.location, amb);
  if (response_time > 3000) {
    return;
  }
  double time_arrival_scene = origin_time + response_time;
  double time_leave_scene = time_arrival_scene + call.time_on_scene;
  double time_arrival_hospital =
      time_leave_scene +
      travel.travel_time(call.location, ins.hospitals[call.hospital], amb);
  double time_finish_service = time_arrival_hospital + call.time_at_hospital;
  double waiting_time = time_arrival_scene - call.time;
  if (debug) {
    fmt::print("========================\n");
    cout << origin.first << " " << origin.second << "\n";
    cout << call << "\n";
    cout << amb << "\n";
    fmt::print("Creating node for amb {} and call {} with waiting_time {}\n",
               amb.id, call.id, waiting_time);
    fmt::print("Response time {}\n", response_time);
    fmt::print("time arrival scene {}\n", time_arrival_scene);
    fmt::print("time leave hospital {}\n", time_finish_service);
    fmt::print("Begin time {}\n", begin_time);
    fmt::print("Waiting time {}\n", waiting_time);
    cin.get();
  }
  node_queue.push(EnumNode{call, waiting_time, parent, time_finish_service});
  if (debug) {
    fmt::print("Child waiting time = {}\n", node_queue.back().waiting_time);
  }
}

Route EnumerateSolver::expand_route(Ambulance &amb, EnumNode &route_node) {
  vector<PathNode> route_calls;
  shared_ptr<EnumNode> current_node = make_shared<EnumNode>(route_node);
  double cost = 0;
  bool debug = false;
  while (current_node != nullptr) {
    auto current_call = current_node->call;
    if (debug) {
      fmt::print("Amb 0 to call {}: {}\n", current_call.id,
                 current_node->waiting_time);
    }
    cost += penalized_response_time(current_node->waiting_time, amb.type,
                                    current_call.priority);
    route_calls.push_back({TypeNode::CALL_NODE, current_call.id});
    current_node = current_node->prev;
  }
  if (debug) {
    cin.get();
  }
  reverse(route_calls.begin(), route_calls.end());

  return Route{amb.id, route_calls, cost};
}

void EnumerateSolver::run() {
  is_call_on_queue = vector<bool>(calls.size(), false);
  for (auto i : queue) {
    is_call_on_queue[i] = true;
  }
  int first_index_call = index_call;
  auto t0 = std::chrono::high_resolution_clock::now();
  int route_ind = 0;
  for (auto &amb : ambulances) {
    auto amb_routes = generate_routes(amb);
    for (auto &route : amb_routes) {
      routes.push_back(route);
      for (size_t i = 0; i < route.node_sequence.size(); ++i) {
        auto &call = calls[route.node_sequence[i].node_index];
        routes_by_call[call.id].insert(route_ind);
      }
      routes_by_ambulance[amb.id].insert(route_ind);
      ++route_ind;
    }
  }
  auto dt = std::chrono::high_resolution_clock::now();
  double time_enum =
      std::chrono::duration_cast<chrono::milliseconds>(dt - t0).count();
  t0 = std::chrono::high_resolution_clock::now();
  load_model();
  // fmt::print("Loaded model\n");
  model.optimize();
  dt = std::chrono::high_resolution_clock::now();
  double time_model =
      std::chrono::duration_cast<chrono::milliseconds>(dt - t0).count();
  // fmt::print("Sizes: ambs = {}, calls = {}, routes = {}. Time enum = {}, time
  // model = {}\n", ambulances.size(), calls.size(), routes.size(),
  // time_enum, time_model);
  obj = model.get(GRB_DoubleAttr_ObjVal);
  run_time = time_enum + time_model;
  index_solver = 0;
}

MCSolver::MCSolver(GRBEnv &env, vector<Call> &calls,
                   vector<Ambulance> &ambulances, Instance &ins, Travel &travel,
                   int g, double time)
    : Solver(env, calls, ambulances, ins, travel, time), g(g) {
  prepare();
  read_preparedness();
}

pair<vector<vector<int>>, vector<int>> MCSolver::partition_ambulances(int a0) {
  vector<vector<int>> A_b(ins.bases.size(), vector<int>());
  vector<int> A_prime;
  for (auto &amb : ambulances) {
    if (amb.arrival_time_at_f_last_trip <= time && amb.id != a0) {
      double min_dist = GRB_INFINITY;
      int min_b = -1;
      for (size_t b = 0; b < ins.bases.size(); ++b) {
        double d = travel.lat_long_distance(amb.base_location, ins.bases[b]);
        if (d < min_dist) {
          min_dist = d;
          min_b = b;
        }
      }
      if (min_b >= 0) {
        A_b[min_b].push_back(amb.id);
      } else {
        fmt::print("Error: Ambulance {} is not busy and not in any base\n",
                   amb.id);
      }
    } else {
      A_prime.push_back(amb.id);
    }
  }

  return make_pair(A_b, A_prime);
}

void MCSolver::read_preparedness() { preparedness = ins.preparedness; }

bool MCSolver::solve_selection(int i0, bool debug) {
  GRBModel model(env);
  vector<vector<int>> A_b;
  vector<int> A_prime;

  auto queue_aux = queue;
  auto i0_iter = find(queue.begin(), queue.end(), i0);
  size_t i0_index = distance(queue.begin(), i0_iter);
  bool is_new_call = i0_index == queue.size();
  if (is_new_call) {
    queue_aux.push_back(i0);
    i0_index = queue.size();
  }
  // fmt::print("queue_aux.size() = {}\n", queue_aux.size());
  vector<double> gamma(queue_aux.size(), 0);
  vector<int> calls_sub_regions(queue_aux.size(), -1);
  for (size_t i = 0; i < queue_aux.size(); ++i) {
    auto &call = calls[queue_aux[i]];
    calls_sub_regions[i] = (g_params.n_ambulances == 40)
                               ? ins.map_region_sub_region[call.region]
                               : 0;
    // fmt::print("call {} / {}\n", queue_aux[i], calls.size());
    double waiting_time = time - call.time;
    gamma[i] =
        penalized_response_time(3600 * (waiting_time + 1), -1, call.priority);
    // gamma[i] = penalized_response_time(1000, -1, call.priority);
  }

  double GAMMA = 1000;
  std::tie(A_b, A_prime) = partition_ambulances();
  vector<int> ambulances_sub_regions(ambulances.size(), -1);
  vector<vector<GRBVar>> x(ambulances.size(), vector<GRBVar>(queue_aux.size()));
  for (size_t a = 0; a < ambulances.size(); ++a) {
    ambulances_sub_regions[a] =
        (g_params.n_ambulances == 40)
            ? ins.map_ambulance_which_sub_region[ambulances[a].id]
            : 0;
    for (size_t i = 0; i < queue_aux.size(); ++i) {
      if (can_answer(ambulances[a], calls[queue_aux[i]]) &&
          ambulances_sub_regions[a] == calls_sub_regions[i]) {
        x[a][i] = model.addVar(0, 1, 0, GRB_BINARY,
                               fmt::format("x_{}_{}", a, queue_aux[i]));
      }
    }
  }
  vector<vector<GRBVar>> y(A_prime.size(), vector<GRBVar>(ins.bases.size()));
  for (size_t a = 0; a < A_prime.size(); ++a) {
    auto &amb = ambulances[A_prime[a]];
    for (size_t b = 0; b < ins.bases.size(); ++b) {
      if (ambulances_sub_regions[amb.id] == ins.map_base_sub_region[b]) {
        y[a][b] = model.addVar(0, 1, 0, GRB_BINARY,
                               fmt::format("y_{}_{}", A_prime[a], b));
      }
    }
  }

  if (debug) {
    fmt::print("selection variables loaded\n");
  }

  GRBLinExpr fo;
  bool force_forward = true;
  for (size_t a = 0; a < ambulances.size(); ++a) {
    for (size_t i = 0; i < queue_aux.size(); ++i) {
      auto &amb = ambulances[a];
      auto &call = calls[queue_aux[i]];
      if (can_answer(amb, call) &&
          ambulances_sub_regions[a] == calls_sub_regions[i]) {
        double travel_time = penalized_response_time(
            travel.get_response_time(amb, call, time, force_forward), amb.type,
            call.priority);
        // double response_time = time + travel_time - call.time;
        double response_time = travel_time;
        fo += response_time * x[a][i];
      }
    }
  }

  for (size_t i = 0; i < queue_aux.size(); ++i) {
    fo += gamma[i];
    for (size_t a = 0; a < ambulances.size(); ++a) {
      if (can_answer(ambulances[a], calls[queue_aux[i]]) &&
          ambulances_sub_regions[a] == calls_sub_regions[i]) {
        fo -= gamma[i] * x[a][i];
      }
    }
  }

  for (size_t b = 0; b < ins.bases.size(); ++b) {
    for (auto a : A_b[b]) {
      auto &amb = ambulances[a];
      for (size_t i = 0; i < queue_aux.size(); ++i) {
        if (can_answer(amb, calls[queue_aux[i]]) &&
            ambulances_sub_regions[a] == calls_sub_regions[i]) {
          fo += GAMMA * s_minus(amb.type, b, time, A_b[b]) * x[a][i];
        }
      }
    }
  }

  for (size_t a = 0; a < A_prime.size(); ++a) {
    auto &amb = ambulances[A_prime[a]];
    for (size_t b = 0; b < ins.bases.size(); ++b) {
      if (ambulances_sub_regions[amb.id] == ins.map_base_sub_region[b]) {
        fo += GAMMA * s_plus(amb.type, b, time, A_b[b]) * y[a][b];
      }
    }
  }

  model.setObjective(fo, GRB_MINIMIZE);

  if (debug) {
    fmt::print("selection objective loaded\n");
  }

  for (size_t b = 0; b < ins.bases.size(); ++b) {
    for (auto a : A_b[b]) {
      GRBLinExpr exp;
      for (size_t i = 0; i < queue_aux.size(); ++i) {
        if (can_answer(ambulances[a], calls[queue_aux[i]]) &&
            ambulances_sub_regions[a] == calls_sub_regions[i]) {
          exp += x[a][i];
        }
      }
      model.addConstr(exp, GRB_LESS_EQUAL, 1,
                      fmt::format("dispatch_idle_amb_{}", a));
    }
  }

  for (size_t a = 0; a < A_prime.size(); ++a) {
    GRBLinExpr exp;
    for (size_t i = 0; i < queue_aux.size(); ++i) {
      if (can_answer(ambulances[a], calls[queue_aux[i]]) &&
          ambulances_sub_regions[a] == calls_sub_regions[i]) {
        exp += x[a][i];
      }
    }

    for (size_t b = 0; b < ins.bases.size(); ++b) {
      if (ambulances_sub_regions[A_prime[a]] == ins.map_base_sub_region[b]) {
        exp += y[a][b];
      }
    }

    model.addConstr(exp, GRB_EQUAL, 1,
                    fmt::format("dispatch_busy_amb_{}", A_prime[a]));
  }

  for (size_t i = 0; i < queue_aux.size(); ++i) {
    GRBLinExpr exp;
    for (size_t a = 0; a < ambulances.size(); ++a) {
      if (can_answer(ambulances[a], calls[queue_aux[i]]) &&
          ambulances_sub_regions[a] == calls_sub_regions[i]) {
        exp += x[a][i];
      }
    }
    model.addConstr(exp, GRB_LESS_EQUAL, 1,
                    fmt::format("assignment_call_{}", queue_aux[i]));
  }

  if (debug) {
    model.write("select.lp");
    fmt::print("Wrote model\n");
    // cin.get();
  }
  model.update();

  model.set(GRB_IntParam_OutputFlag, 0);
  model.optimize();

  auto status = model.get(GRB_IntAttr_Status);
  if (status != GRB_OPTIMAL) {
    fmt::print("Error: selection model failed with status {}\n", status);
  }
  double run_time = model.get(GRB_DoubleAttr_Runtime);
  run_times.push_back(run_time);
  int index_amb = -1;
  int count = 0;
  for (size_t a = 0; a < ambulances.size(); ++a) {
    auto &amb = ambulances[a];
    if (can_answer(ambulances[a], calls[queue_aux[i0_index]]) &&
        ambulances_sub_regions[a] == calls_sub_regions[i0_index]) {
      double val = x[a][i0_index].get(GRB_DoubleAttr_X);

      if (val > 0.5) {
        index_amb = amb.id;
        ++count;
        // break;
      }
    }
  }

  if (debug) {
    fmt::print("index_amb = {}\n", index_amb);
    // cin.get();
  }

  if (count > 1) {
    fmt::print("Error: call {} was answered by {} ambulances\n", i0, count);
    exit(100);
  }

  if (index_amb != -1) {
    auto &amb = ambulances[index_amb];
    auto &call = calls[i0];
    double real_min_time =
        travel.get_response_time(amb, call, time, force_forward);
    double waiting_on_scene_i = time + real_min_time - call.time;
    double waiting_to_hospital_i = amb.answer_call(
        call, travel, ins, time, real_min_time, nearest_base[i0]);
    waiting_on_scene[i0] = waiting_on_scene_i;
    waiting_on_scene_penalized[i0] =
        penalized_response_time(waiting_on_scene_i, amb.type, call.priority);
    waiting_to_hospital[i0] = waiting_to_hospital_i;
    which_ambulance[i0] = index_amb;
    calls_end[i0] = call.end;

    if (debug) {
      auto sub_region = ins.map_region_sub_region[call.region];
      auto amb_sub_region = ambulances_sub_regions[index_amb];
      fmt::print(
          "Call {}, sr {}, answered by {} sr {}, finish_service {:.1f}, "
          "Real_min_time "
          "= {:.1f}\n",
          call.id, sub_region, amb.id, amb_sub_region,
          amb.arrival_time_at_f_last_trip, real_min_time);
    }
    if (!is_new_call) {
      queue.erase(remove(queue.begin(), queue.end(), i0), queue.end());
    }
    return true;
  } else {
    if (debug) {
      fmt::print("Call {} went to queue\n", i0);
    }
    if (is_new_call) {
      queue.push_back(i0);
    }
    return false;
  }
  return false;
}

bool MCSolver::solve_reassignment(int a0, bool debug) {
  GRBModel model(env);
  vector<vector<int>> A_b;
  vector<int> A_prime;
  vector<double> gamma(queue.size(), 0);
  vector<int> calls_sub_regions(queue.size(), -1);
  for (size_t i = 0; i < queue.size(); ++i) {
    auto &call = calls[queue[i]];
    calls_sub_regions[i] = (g_params.n_ambulances == 40)
                               ? ins.map_region_sub_region[call.region]
                               : 0;
    double waiting_time = time - call.time;
    gamma[i] =
        penalized_response_time(3600 * (waiting_time + 1), -1, call.priority);
    // gamma[i] = penalized_response_time(1000, -1, call.priority);
  }
  double GAMMA = 1000;

  std::tie(A_b, A_prime) = partition_ambulances(a0);
  if (debug) {
    fmt::print("Reassigning {}\n", a0);
    for (size_t b = 0; b < ins.bases.size(); ++b) {
      fmt::print("ambs at base {}: {}\n", b, A_b[b]);
    }
    fmt::print("A_prime = {}\n", A_prime);
  }
  // debug = time > calls.back().time && A_prime.size() == 1;
  vector<int> ambulances_sub_regions(ambulances.size(), -1);
  vector<vector<GRBVar>> x(ambulances.size(), vector<GRBVar>(queue.size()));
  for (size_t a = 0; a < ambulances.size(); ++a) {
    ambulances_sub_regions[a] =
        (g_params.n_ambulances == 40)
            ? ins.map_ambulance_which_sub_region[ambulances[a].id]
            : 0;
    for (size_t i = 0; i < queue.size(); ++i) {
      if (can_answer(ambulances[a], calls[queue[i]]) &&
          ambulances_sub_regions[a] == calls_sub_regions[i]) {
        x[a][i] = model.addVar(0, 1, 0, GRB_BINARY,
                               fmt::format("x_{}_{}", a, queue[i]));
      }
    }
  }

  vector<vector<GRBVar>> y(A_prime.size(), vector<GRBVar>(ins.bases.size()));
  int a0_index = -1;
  for (size_t a = 0; a < A_prime.size(); ++a) {
    if (A_prime[a] == a0) {
      a0_index = a;
    }
    for (size_t b = 0; b < ins.bases.size(); ++b) {
      if (ambulances_sub_regions[A_prime[a]] == ins.map_base_sub_region[b]) {
        y[a][b] = model.addVar(0, 1, 0, GRB_BINARY,
                               fmt::format("y_{}_{}", A_prime[a], b));
      }
    }
  }

  if (debug) {
    fmt::print("reassignment variables loaded\n");
  }

  // if(debug){
  // 	for(size_t b = 0; b < ins.bases.size(); ++b){
  // 		fmt::print("A_{} = {}\n", b, A_b[b]);
  // 	}
  // 	fmt::print("A_prime = {}\n", A_prime);
  // 	fmt::print("a0 = {}, a0_index = {}\n", a0, a0_index);
  // }

  model.update();

  GRBLinExpr fo;
  for (size_t a = 0; a < ambulances.size(); ++a) {
    for (size_t i = 0; i < queue.size(); ++i) {
      auto &amb = ambulances[a];
      auto &call = calls[queue[i]];
      if (can_answer(amb, call) &&
          ambulances_sub_regions[a] == calls_sub_regions[i]) {
        double travel_time = penalized_response_time(
            travel.get_response_time(amb, call, time), amb.type, call.priority);
        // double response_time = time + travel_time - call.time;
        double response_time = travel_time;
        fo += response_time * x[a][i];
      }
    }
  }
  if (debug) {
    fmt::print("resp_times\n");
  }
  for (size_t i = 0; i < queue.size(); ++i) {
    fo += gamma[i];
    for (size_t a = 0; a < ambulances.size(); ++a) {
      if (can_answer(ambulances[a], calls[queue[i]]) &&
          ambulances_sub_regions[a] == calls_sub_regions[i]) {
        if (debug) {
          fmt::print("coeff {} = {}\n", x[a][i].get(GRB_StringAttr_VarName),
                     gamma[i]);
        }
        fo -= gamma[i] * x[a][i];
      }
    }
  }
  // if (debug) {
  //   fmt::print("gamma\n");
  // }

  for (size_t b = 0; b < ins.bases.size(); ++b) {
    for (auto a : A_b[b]) {
      auto &amb = ambulances[a];
      for (size_t i = 0; i < queue.size(); ++i) {
        if (can_answer(amb, calls[queue[i]]) &&
            ambulances_sub_regions[a] == calls_sub_regions[i]) {
          double sm = s_minus(amb.type, b, time, A_b[b]);
          // if (debug) {
          //   fmt::print("coeff - {} = {}, GAMMA = {}, sm = {}\n",
          //              x[a][i].get(GRB_StringAttr_VarName), GAMMA * sm,
          //              GAMMA, sm);
          // }
          fo += GAMMA * s_minus(amb.type, b, time, A_b[b]) * x[a][i];
        }
      }
    }
  }

  for (size_t a = 0; a < A_prime.size(); ++a) {
    auto &amb = ambulances[A_prime[a]];
    for (size_t b = 0; b < ins.bases.size(); ++b) {
      if (ambulances_sub_regions[A_prime[a]] == ins.map_base_sub_region[b]) {
        double sp = s_plus(amb.type, b, time, A_b[b]);
        // if (debug) {
        //   fmt::print("coeff + {} = {}, GAMMA = {}, sp = {}\n",
        //              y[a][b].get(GRB_StringAttr_VarName), GAMMA * sp, GAMMA,
        //              sp);
        // }
        fo += GAMMA * sp * y[a][b];
      }
    }
  }
  // if (debug) {
  //   fmt::print("GAMMA\n");
  // }

  model.setObjective(fo, GRB_MINIMIZE);

  if (debug) {
    fmt::print("reassignment objective loaded\n");
  }

  for (size_t b = 0; b < ins.bases.size(); ++b) {
    for (auto a : A_b[b]) {
      GRBLinExpr exp;
      for (size_t i = 0; i < queue.size(); ++i) {
        if (can_answer(ambulances[a], calls[queue[i]]) &&
            ambulances_sub_regions[a] == calls_sub_regions[i]) {
          exp += x[a][i];
        }
      }
      model.addConstr(exp, GRB_LESS_EQUAL, 1,
                      fmt::format("dispatch_idle_amb_{}", a));
    }
  }

  for (size_t a = 0; a < A_prime.size(); ++a) {
    GRBLinExpr exp;
    for (size_t i = 0; i < queue.size(); ++i) {
      if (can_answer(ambulances[A_prime[a]], calls[queue[i]]) &&
          ambulances_sub_regions[A_prime[a]] == calls_sub_regions[i]) {
        exp += x[A_prime[a]][i];
      }
    }

    for (size_t b = 0; b < ins.bases.size(); ++b) {
      if (ambulances_sub_regions[A_prime[a]] == ins.map_base_sub_region[b]) {
        exp += y[a][b];
      }
    }
    model.addConstr(exp, GRB_EQUAL, 1,
                    fmt::format("dispatch_busy_amb_{}", A_prime[a]));
  }

  for (size_t i = 0; i < queue.size(); ++i) {
    GRBLinExpr exp;
    for (size_t a = 0; a < ambulances.size(); ++a) {
      if (can_answer(ambulances[a], calls[queue[i]]) &&
          ambulances_sub_regions[a] == calls_sub_regions[i]) {
        exp += x[a][i];
      }
    }
    model.addConstr(exp, GRB_LESS_EQUAL, 1,
                    fmt::format("assignment_call_{}", queue[i]));
  }
  model.update();
  if (debug) {
    model.write("reassign.lp");
    // fmt::print("Wrote reassignment model.\n");
    // cin.get();
  }

  model.set(GRB_IntParam_OutputFlag, 0);
  model.optimize();

  auto status = model.get(GRB_IntAttr_Status);
  if (status != GRB_OPTIMAL) {
    fmt::print("Error: reassigment model failed with status {}\n", status);
  }

  double run_time = model.get(GRB_DoubleAttr_Runtime);
  run_times.push_back(run_time);
  int call_index = -1;
  int count = 0;
  for (size_t i = 0; i < queue.size(); ++i) {
    if (can_answer(ambulances[a0], calls[queue[i]]) &&
        ambulances_sub_regions[a0] == calls_sub_regions[i]) {
      double val = x[a0][i].get(GRB_DoubleAttr_X);
      if (val > 0.5) {
        call_index = queue[i];
        // fmt::print("({}, {}) ", a0_index, queue[i]);
        ++count;
        // break;
      }
    }
  }
  // fmt::print("\n");

  if (count > 1) {
    fmt::print("ERROR: amb {} answered {} calls.\n", a0, count);
    exit(100);
  }

  if (debug && call_index != -1) {
    fmt::print("call_index = {}\n", call_index);
    // cin.get();
  }

  if (call_index != -1) {
    auto &amb = ambulances[a0];
    auto &call = calls[call_index];
    double real_min_time = travel.get_response_time(amb, call, time);
    double waiting_on_scene_i = time + real_min_time - call.time;
    double waiting_to_hospital_i = amb.answer_call(
        call, travel, ins, time, real_min_time, nearest_base[call_index]);
    waiting_on_scene[call_index] = waiting_on_scene_i;
    waiting_on_scene_penalized[call_index] =
        penalized_response_time(waiting_on_scene_i, amb.type, call.priority);
    waiting_to_hospital[call_index] = waiting_to_hospital_i;
    which_ambulance[call_index] = amb.id;
    calls_end[call_index] = call.end;
    queue.erase(remove(queue.begin(), queue.end(), call_index), queue.end());
    if (debug) {
      auto sub_region = ins.map_region_sub_region[call.region];
      auto amb_sub_region = ambulances_sub_regions[a0];
      fmt::print(
          "Call {} sr {} answered by {} sr {}, finish_service {:.1f}. "
          "Real_mean_time = {:.1f}\n",
          call.id, sub_region, amb.id, amb_sub_region,
          amb.arrival_time_at_f_last_trip, real_min_time);
    }
    return true;
  } else {
    int index_base = -1;
    count = 0;
    for (size_t b = 0; b < ins.bases.size(); ++b) {
      if (ambulances_sub_regions[a0] == ins.map_base_sub_region[b]) {
        double val = y[a0_index][b].get(GRB_DoubleAttr_X);
        if (val > 0.5) {
          index_base = b;
          ++count;
        }
      }
    }
    if (debug) {
      fmt::print("Amb {} must return to {}\n", a0_index, index_base);
      fmt::print("Queue {}\n", queue);
      model.write("debug_reassignment.lp");
      GRBVar *vars = model.getVars();
      for (int i = 0; i < model.get(GRB_IntAttr_NumVars); ++i) {
        string name = vars[i].get(GRB_StringAttr_VarName);
        double val = vars[i].get(GRB_DoubleAttr_X);
        if (val > 0.001) {
          fmt::print("{} => {}\n", name, val);
        }
      }
      fmt::print("========================\n");
      for (size_t i = 0; i < queue.size(); ++i) {
        fmt::print("call {} sub_region = {}\n", calls[queue[i]].id,
                   calls_sub_regions[i]);
      }
      fmt::print("========================\n");

      // cin.get();
      delete[] vars;
      // cin.get();
    }

    if (count > 1) {
      fmt::print("ERROR: amb {} was reassigned to {} bases.\n", a0, count);
      exit(100);
    }

    if (index_base == -1) {
      fmt::print(
          "ERROR: amb {} that was not dispatched was not reassigned either.\n",
          a0);
      exit(50);
    } else {
      auto &amb = ambulances[a0];
      int b = index_base;
      double arrival_time =
          travel.travel_time(amb.free_location, ins.bases[b], amb);
      amb.reassign_ambulance(ins.bases[b], arrival_time);
    }
    return false;
  }
  return false;
}

void MCSolver::run() {
  calls_attended = 0;
  bool debug = false;
  using fmt::print;
  default_random_engine u_gen;
  travel.set_forward(true);
  uniform_real_distribution<double> u(0, 1);
  while (static_cast<size_t>(calls_attended) < calls.size()) {
    travel.set_delay_param(u(u_gen));
    if (debug) {
      auto sub_region =
          (index_call >= 0 && static_cast<size_t>(index_call) < calls.size())
              ? ins.map_region_sub_region[calls[index_call].region]
              : -1;
      fmt::print(
          "MarkovSolver: Time {:.1f}, index = {}, sub_region = {}, event = {}, "
          "queue = "
          "{}, att = {}/{}\n",
          time, index_call, sub_region, event_call, queue, calls_attended,
          calls.size());
    }

    if (event_call == 1) {
      auto &call = calls[index_call];
      bool call_attended = solve_selection(index_call, debug);
      if (call_attended) {
        ++calls_attended;
      }
    } else {
      int a0 = -1;
      for (auto &amb : ambulances) {
        if (amb.arrival_time_at_f_last_trip == time) {
          a0 = amb.id;
        }
      }

      if (a0 == -1) {
        fmt::print("ERROR: Event_call = 0 but no ambulance is free\n");
        exit(1);
      }

      bool call_attended = solve_reassignment(a0, debug);
      if (call_attended) {
        ++calls_attended;
      }
    }
    set_next_event(event_call, index_call);
    if (debug) {
      cin.get();
    }
  }
  bool has_unsolved_calls = false;
  for (auto &call : calls) {
    if (which_ambulance[call.id] == -1) {
      has_unsolved_calls = true;
      break;
    }
  }

  if (has_unsolved_calls) {
    fmt::print("Some calls were unsolved\n");
    print_results();
    cin.get();
  }

  if (debug) {
    print_results();
    // cin.get();
  }
}

int MCSolver::get_time_slot(double time) {
  if (time < ins.time_horizon[0] * 1800) {
    fmt::print("ERROR: time {} before time_horizon begin {}\n", time,
               ins.time_horizon[0] * 1800);
    exit(1);
  }

  for (size_t i = 0; i < ins.time_horizon.size() - 1; ++i) {
    double slot_begin = ins.time_horizon[i] * 1800;
    double slot_end = ins.time_horizon[i + 1] * 1800;

    if (slot_begin <= time && time < slot_end) {
      return i;
    }
  }

  return ins.time_horizon.size() - 1;
}

double MCSolver::s_minus(int amb_type, int b, double a_time,
                         vector<int> &supply_b) {
  int t = get_time_slot(a_time);
  vector<int> m(ins.nb_types_ambulance, 0);
  for (auto a : supply_b) {
    auto &amb = ambulances[a];
    ++m[amb.type];
  }
  bool setup_us = ins.nb_types_ambulance == 2;
  double prep_before = (setup_us) ? preparedness(b, t, m[0], m[1])
                                  : preparedness(b, t, m[0], m[1], m[2]);
  --m[amb_type];
  double prep_after = (setup_us) ? preparedness(b, t, m[0], m[1])
                                 : preparedness(b, t, m[0], m[1], m[2]);
  // fmt::print("after {} {} {} {} {} = {} | before {} {} {} {} {} = {}\n",
  // 	b,t,m[0],m[1], m[2], prep_after, b,t,m[0],m[1], m[2], prep_before);
  return prep_after - prep_before;
}

double MCSolver::s_plus(int amb_type, int b, double a_time,
                        vector<int> &supply_b) {
  int t = get_time_slot(a_time);
  vector<int> m(ins.nb_types_ambulance, 0);
  for (auto a : supply_b) {
    auto &amb = ambulances[a];
    ++m[amb.type];
  }
  bool setup_us = ins.nb_types_ambulance == 2;
  double prep_before = (setup_us) ? preparedness(b, t, m[0], m[1])
                                  : preparedness(b, t, m[0], m[1], m[2]);
  ++m[amb_type];
  double prep_after = (setup_us) ? preparedness(b, t, m[0], m[1])
                                 : preparedness(b, t, m[0], m[1], m[2]);
  // fmt::print("after {} {} {} {} {} = {} | before {} {} {} {} {} = {}\n",
  // 	b,t,m[0],m[1], m[2], prep_after, b,t,m[0],m[1], m[2], prep_before);
  return prep_after - prep_before;
}

QueueDeficitSolver::QueueDeficitSolver(GRBEnv &env, vector<Call> &calls,
                                       vector<Ambulance> &ambulances,
                                       Instance &ins, Travel &travel, int g,
                                       double time)
    : Solver(env, calls, ambulances, ins, travel, time), g(g) {
  prepare();
  travel.set_forward(false);
}

void QueueDeficitSolver::run() {
  calls_attended = 0;
  bool debug = false;
  using fmt::print;
  default_random_engine u_gen;
  uniform_real_distribution<double> u(0, 1);
  while (static_cast<size_t>(calls_attended) < calls.size()) {
    travel.set_delay_param(u(u_gen));
    for (auto &amb : ambulances) {
      double arrival_b = amb.arrival_time_at_b_last_trip;
      if (arrival_b <= time) {
        for (size_t c = 0; c < calls.size(); ++c) {
          calls[c].ambulances.push_back(amb.id);
        }
      }
    }
    if (event_call == 1) {
      queue.push_back(index_call);
    }
    if (debug) {
      fmt::print(
          "QueueDeficitSolver: Time {}, index = {}, event = {}, queue = "
          "{}, att = {}/{}\n",
          time, index_call, event_call, queue, calls_attended, calls.size());
    }
    std::vector<int> queue_aux;
    for (size_t i = 0; i < queue.size(); ++i) {
      auto t0 = std::chrono::high_resolution_clock::now();
      auto &call = calls[queue[i]];
      double min_time = GRB_INFINITY;
      int index_amb = -1;
      for (auto &amb : ambulances) {
        if (can_answer(amb, call)) {
          double time_amb =
              penalized_response_time(travel.get_response_time(amb, call, time),
                                      amb.type, call.priority);
          // print("\tAmb {} {} {}\n", amb.id, time_amb, amb.type);
          if (time_amb < min_time) {
            min_time = time_amb;
            index_amb = amb.id;
          }
        }
      }

      if (index_amb >= 0) {
        auto &amb = ambulances[index_amb];

        double real_min_time = travel.get_response_time(amb, call, time);
        double waiting_on_scene_i = time + real_min_time - call.time;
        int best_base = find_best_base(amb.id, debug);
        double waiting_to_hospital_i =
            amb.answer_call(call, travel, ins, time, real_min_time, best_base);
        waiting_on_scene[queue[i]] = waiting_on_scene_i;
        waiting_on_scene_penalized[queue[i]] = penalized_response_time(
            waiting_on_scene_i, amb.type, call.priority);
        waiting_to_hospital[queue[i]] = waiting_to_hospital_i;
        which_ambulance[queue[i]] = index_amb;
        calls_end[queue[i]] = call.end;
        if (debug) {
          fmt::print("Amb {} calls {} and then goes to base {}\n", amb.id,
                     call.id, best_base);
        }
        calls_attended++;
      } else {
        queue_aux.push_back(queue[i]);
        if (debug) {
          fmt::print("Call {} back to queue\n", call.id);
        }
      }
      auto dt = std::chrono::high_resolution_clock::now();
      run_times.push_back(
          std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
    }
    queue.clear();
    for (auto i : queue_aux) {
      queue.push_back(i);
    }
    if (debug) {
      cin.get();
    }
    set_next_event(event_call, index_call);
  }
  // print_results();
  if (debug) {
    print_results();
    cin.get();
  }
}

CentralitySolver::CentralitySolver(GRBEnv &env, vector<Call> &calls,
                                   vector<Ambulance> &ambulances, Instance &ins,
                                   Travel &travel, int g, double time)
    : Solver(env, calls, ambulances, ins, travel, time), g(g) {
  prepare();
  travel.set_forward(false);
  int count_hosp_or_clean = 0;
  for (auto &call : calls) {
    if (call.hosp_needed || call.clean_needed) {
      ++count_hosp_or_clean;
    }
  }
  transport_prob = static_cast<double>(count_hosp_or_clean) / calls.size();
}

std::vector<std::vector<double>> CentralitySolver::compute_fitness() {
  std::vector<double> centrality(queue.size(), -1.0);
  for (size_t i = 0; i < queue.size(); ++i) {
    auto call_i = calls[queue[i]];
    double sum_central = 0.0;
    for (size_t j = 0; j < queue.size(); ++j) {
      if (i != j) {
        auto call_j = calls[queue[j]];
        double time_to_call =
            travel.travel_time(call_i.location, call_j.location, ambulances[0]);
        sum_central += 1 / (1 + (time_to_call / 60));
      }
    }
    centrality[i] = sum_central;
  }
  if (debug_mode) {
    fmt::print("Centrality = {}\n", centrality);
  }
  std::vector<std::vector<double>> fitness(
      ambulances.size(), std::vector<double>(queue.size(), -1.0));

  for (size_t a = 0; a < ambulances.size(); ++a) {
    auto &amb = ambulances[a];
    for (size_t c = 0; c < queue.size(); ++c) {
      auto call = calls[queue[c]];
      double time_to_call = travel.get_response_time(amb, call, time, true);
      double central = (abs(centrality[c]) < g_params.EPS) ? 1 : centrality[c];
      fitness[a][c] =
          pow(central, 1 - transport_prob) / (1 + (time_to_call / 60));
    }
  }

  return fitness;
}

std::vector<int> CentralitySolver::get_assignment() {
  std::vector<int> assignment(queue.size(), -1);
  auto fitness = compute_fitness();
  std::vector<std::tuple<double, int, int>> assigns_by_fitness;
  if (debug_mode) {
    fmt::print("FITNESS:\n");
  }
  for (size_t a = 0; a < ambulances.size(); ++a) {
    for (size_t c = 0; c < queue.size(); ++c) {
      if (debug_mode) {
        fmt::print("\ta = {}, c = {}, fitness = {}\n", a, c, fitness[a][c]);
      }
      assigns_by_fitness.push_back(std::make_tuple(fitness[a][c], a, c));
    }
  }
  if (debug_mode) {
    fmt::print("========== fitness end ============\n");
  }
  std::sort(assigns_by_fitness.begin(), assigns_by_fitness.end(),
            [](const auto &a, const auto &b) {
              return std::get<0>(a) > std::get<0>(b);
            });
  double fit;
  int a, c;

  std::vector<bool> available(ambulances.size(), false);
  for (size_t a = 0; a < ambulances.size(); ++a) {
    if (ambulances[a].arrival_time_at_f_last_trip <= time) {
      available[a] = true;
    }
  }
  std::vector<bool> assigned(queue.size(), false);
  for (size_t i = 0; i < assigns_by_fitness.size(); ++i) {
    double fit = std::get<0>(assigns_by_fitness[i]);
    int a = std::get<1>(assigns_by_fitness[i]);
    int c = std::get<2>(assigns_by_fitness[i]);
    if (available[a] && !assigned[c]) {
      assignment[c] = a;
      available[a] = false;
      assigned[c] = true;
    }
  }

  for (size_t c = 0; c < queue.size(); ++c) {
    if (!assigned[c]) {
      double best_fit = -GRB_INFINITY;
      int index_amb = -1;
      for (size_t a = 0; a < ambulances.size(); ++a) {
        if (available[a] && fitness[a][c] > best_fit) {
          best_fit = fitness[a][c];
          index_amb = a;
        }
      }

      if (index_amb != -1) {
        assignment[c] = index_amb;
        assigned[c] = true;
        available[index_amb] = false;
      }
    }
  }

  return assignment;
}

void CentralitySolver::run() {
  // debug_mode = true;
  bool debug = debug_mode;
  calls_attended = 0;
  travel.set_forward(false);
  default_random_engine u_gen;
  uniform_real_distribution<double> u(0, 1);

  while (static_cast<size_t>(calls_attended) < calls.size()) {
    travel.set_delay_param(u(u_gen));
    if (debug) {
      fmt::print(
          "Centrality: Time {}, event {}, index {}, attended {}/{}, queue {}\n",
          time, event_call, index_call, calls_attended, calls.size(), queue);
    }
    auto t0 = std::chrono::high_resolution_clock::now();
    if (event_call == 1 && static_cast<size_t>(index_call) < calls.size()) {
      queue.push_back(index_call);
    }
    auto assignment = get_assignment();
    std::vector<int> queue_aux;
    if (debug) {
      fmt::print("Ambulances:\n");
      for (auto &amb : ambulances) {
        fmt::print(
            "amb {}, type {}, arr_f = {:.1f}/{:.1f}, arr_b = {:.1f}/{:.1f}\n",
            amb.id, amb.type, amb.arrival_time_at_f_last_trip, time,
            amb.arrival_time_at_b_last_trip, time);
      }
      fmt::print("Assignment:\n");
      for (size_t i = 0; i < queue.size(); ++i) {
        auto call = calls[queue[i]];
        if (assignment[i] >= 0) {
          auto &amb = ambulances[assignment[i]];
          fmt::print("Call {}, type {}, time {:.1f} assigned to {}\n", call.id,
                     call.priority, call.time, amb.id);
        } else {
          fmt::print("Call {}, type {}, time {:.1f} must go to the queue\n",
                     call.id, call.priority, call.time);
        }
      }
    }
    for (size_t i = 0; i < queue.size(); ++i) {
      if (assignment[i] >= 0) {
        auto &call = calls[queue[i]];
        auto &amb = ambulances[assignment[i]];
        if (debug) {
          fmt::print("Call {} type {}, time {:.1f}, answered by {}\n", call.id,
                     call.priority, call.time, amb.id);
        }
        double real_min_time = travel.get_response_time(amb, call, time);
        double waiting_to_hospital_i = amb.answer_call(
            call, travel, ins, time, real_min_time, nearest_base[call.id]);
        double waiting_on_scene_i = time + real_min_time - call.time;
        waiting_on_scene[call.id] = waiting_on_scene_i;
        waiting_on_scene_penalized[call.id] = penalized_response_time(
            waiting_on_scene_i, amb.type, call.priority);
        which_ambulance[call.id] = assignment[i];
        calls_end[call.id] = call.end;
        calls_attended++;
      } else {
        queue_aux.push_back(queue[i]);
        if (debug) {
          fmt::print("Call {} back to queue\n", queue[i]);
        }
      }
    }
    auto dt = std::chrono::high_resolution_clock::now();
    run_times.push_back(
        std::chrono::duration_cast<chrono::microseconds>(dt - t0).count());
    queue.clear();
    for (auto i : queue_aux) {
      queue.push_back(i);
    }
    set_next_event(event_call, index_call);
    if (debug) {
      std::cin.get();
    }
  }
  if (debug) {
    print_results();
    cin.get();
  }
}

// if(min_prep < g_params.min_preparedness){
// 	auto x = get_relocations_heuristic(reloc_available, min_prep_r, lambda);
// 	for(int k = 0; k < reloc_available.size(); ++k){
// 		for(int r = 0; r < ins.nb_regions; ++r){
// 			if(x[k][r]){
// 				auto& amb = ambulances[reloc_available[k]];
// 				Location former_base = amb.base_location;
// 				auto new_base =
// ins.bases[ins.nearest_base_to_region[r]]; 				double
// travel_time = GRB_INFINITY;
// if(amb.arrival_time_at_b_last_trip <= time){
// travel_time = travel.travel_time(former_base, new_base, amb); }else{
// auto current_location = travel.ambulance_position(amb, time); travel_time =
// travel.travel_time(current_location, new_base, amb);
// 				}
// 				amb.base_location = new_base;
// 				amb.arrival_time_at_b_last_trip = time +
// travel_time;
// 			}
// 		}
// 	}
// }

// vector<vector<bool>> Prep2Solver::get_relocations_exact(vector<int>&
// reloc_available){ 	GRBModel model(env); 	stringstream name; 	double
// max_r_time = 1200; 	GRBVar z = model.addVar(0, max_r_time, 1.0,
// GRB_CONTINUOUS, "z"); 	vector<vector<GRBVar>> x(reloc_available.size(),
// vector<GRBVar>()); 	vector<vector<pair<double, int>>>
// reach(reloc_available.size(), 		vector<pair<double, int>>());
// 	x.reserve(reloc_available.size());
// 	for(int k = 0; k < reloc_available.size(); ++k){
// 		for(int r = 0; r < ins.nb_regions; ++r){
// 			// Create variable
// 			name << "x_" << k << "_" << r;
// 			x[k].push_back(model.addVar(0,1,0,GRB_BINARY,name.str()));
// 			name.str("");

// 			// Calculate reach
// 			auto& amb = ambulances[reloc_available[k]];
// 			double travel_time = GRB_INFINITY;
// 			if(amb.arrival_time_at_b_last_trip <= time){
// 				travel_time =
// travel.travel_time(amb.base_location, ins.centers[r], amb);
// }else{ 				auto current_location =
// travel.ambulance_position(amb, time); travel_time =
// travel.travel_time(current_location, ins.centers[r], amb);
// 			}
// 			if(travel_time < max_r_time){
// 				reach[k].push_back(make_pair(travel_time, r));
// 			}
// 		}
// 	}

// 	for(int k = 0; k < reloc_available.size(); ++k){
// 		GRBLinExpr lb_z;
// 		for(auto& aux: reach[k]){
// 			lb_z += aux.first*x[k][aux.second];
// 		}

// 		name << "z_con_" << k;
// 		model.addConstr(z, GRB_GREATER_EQUAL, lb_z, name.str());
// 		name.str("");
// 	}

// 	for(int k = 0; k < reloc_available.size(); ++k){
// 		GRBLinExpr uniqueness;
// 		for(auto& aux: reach[k]){
// 			uniqueness += x[k][aux.second];
// 		}
// 		name << "unique_" << k;
// 		model.addConstr(uniqueness, GRB_LESS_EQUAL, 1, name.str());
// 		name.str("");
// 	}

// 	GRBLinExpr ub_ambs;
// 	for(int k = 0; k < reloc_available.size(); ++k){
// 		for(auto& aux: reach[k]){
// 			ub_ambs += aux.first*x[k][aux.second];
// 		}
// 	}

// 	name << "ub_ambs";
// 	model.addConstr(ub_ambs, GRB_LESS_EQUAL, reloc_available.size(),
// name.str()); 	name.str("");

// 	//last constraint

// 	model.optimize();
// 	vector<vector<bool>> result(ambulances.size(),
// vector<bool>(ins.nb_regions, false)); 	for(int k = 0; k <
// reloc_available.size(); ++k){ 		for(int r = 0; r <
// ins.nb_regions; ++r){ 			if(x[k][r].get(GRB_DoubleAttr_X)
// > 0.5){ 				result[k][r] = true;
// 			}
// 		}
// 	}

// 	return result;
// }

// vector<vector<bool>> Prep2Solver::get_relocations_heuristic(vector<int>&
// reloc_available, 	int min_prep_r, vector<double>& lambda){
// 	vector<vector<bool>> x(reloc_available.size(),
// vector<bool>(ins.nb_regions, false)); 	int n = 5, m = 3, iter = 0,
// max_iter = 50; 	do{ 		vector<pair<double, int>> closest_ambs;
// for(int k = 0; k < reloc_available.size(); ++k){ 			auto&
// amb = ambulances[reloc_available[k]]; 			double
// travel_time = GRB_INFINITY;
// if(amb.arrival_time_at_b_last_trip <= time){
// travel_time = travel.travel_time(amb.base_location, ins.centers[min_prep_r],
// amb); 			}else{ 				auto
// current_location = travel.ambulance_position(amb, time); travel_time =
// travel.travel_time(current_location, ins.centers[min_prep_r],
// amb);
// 			}
// 			closest_ambs.push_back(make_pair(travel_time, k));
// 		}
// 		sort(closest_ambs.begin(), closest_ambs.end());
// 		int k = 0;
// 		while(k < n && k < closest_ambs.size()){

// 			++k;
// 		}
// 	}while(min_prep_r < g_params.min_preparedness && iter < max_iter);

// 	return x;
// }

// ModelSolver::ModelSolver(GRBEnv& env, vector<Call>& calls,
// 	vector<Ambulance>& ambulances, Instance& ins, Travel& travel):
// 	Solver(env, calls,ambulances,ins,travel){
// 	g_params.extended_model = false;
// 	travel.set_forward(false);
// 	queue.reserve(calls.size());
// }

// void ModelSolver::run(){
// 	int event_call = 1;
// 	time = calls[0].time;
// 	int index_call = 0;
// 	int calls_attended = 0;

// 	using fmt::print;
// 	std::vector<std::set<int>> H(calls.size(), std::set<int>());

// 	for(int i = 0; i < calls.size(); ++i){
// 		std::vector<std::pair<double, int>> sorted_hospitals;
// 		sorted_hospitals.reserve(ins.nb_hospitals);
// 		for(int h = 0; h < ins.nb_hospitals; ++h){
// 			double d = travel.vehicle_distance(calls[i].location,
// ins.hospitals[h]);
// sorted_hospitals.push_back(std::make_pair(d,h));
// 		}
// 		std::sort(sorted_hospitals.begin(), sorted_hospitals.end());
// 		for(int k = 0; k < g_params.n_nearest_hospitals && k <
// ins.nb_hospitals; ++k){
// H[i].insert(sorted_hospitals[k].second);
// 		}
// 	}
// 	while(calls_attended < calls.size()){
// 		//Find future calls in each scenario.
// 		vector<int> index_scenarios(ins.calls.size(), 0);
// 		for(int j = 0; j < ins.calls.size(); ++j){
// 			while(index_scenarios[j] < ins.calls[j].size() &&
// 				ins.calls[j][index_scenarios[j]].time <= time){
// 				++index_scenarios[j];
// 			}
// 		}
// 		// fmt::print("Time {}, Queue Size {}\n",time, queue.size());
// 		if(event_call == 1){
// 			auto& call = calls[index_call];
// 			Data data(this, &call, NULL);
// 			double min_avg_cost = GRB_INFINITY;
// 			int best_amb = -1;
// 			int best_h = -1;
// 			//for each ambulance capable of answering call
// 			std::cout << "Call " << call << ":\n";
// 			for(int a = 0; a < ambulances.size(); ++a){
// 				auto& amb = ambulances[a];
// 				if(amb.arrival_time_at_f_last_trip <= time &&
// can_answer(amb,call)){ 					int l0 = -1, b0
// = -1;

// 					double min_d = GRB_INFINITY;
// 					for(int b = 0; b < ins.nb_bases; ++b){
// 						double d =
// travel.vehicle_distance(amb.base_location, ins.bases[b]);
// if(d < min_d){ 							min_d =
// d; 							b0 = b;
// 						}
// 					}

// 					if(amb.arrival_time_at_b_last_trip >
// time){ 						auto current_location =
// travel.ambulance_position(amb,time);
// l0 = data.get_location_index(current_location);
// 					}

// 					for(auto h: H[call.id]){
// 						//Average over all scenarios
// 						double sum_cost = 0;
// 						for(int j = 0; j <
// ins.calls.size(); ++j){
// std::vector<Call> future_calls(ins.calls[j].begin() +
// index_scenarios[j], ins.calls[j].end()); FutureCall fc(data,env,*this,
// amb.type, l0,b0,h,future_calls); fc.solve();
// sum_cost += fc.obj;
// 						}

// 						double avg_cost = sum_cost /
// ins.calls.size(); 						fmt::print("a{}
// l{} b{} h{} => {}\n",amb.type, l0,b0,h, avg_cost);

// 						if(avg_cost < min_avg_cost){
// 							min_avg_cost = avg_cost;
// 							best_amb = a;
// 							best_h = h;
// 						}
// 					}
// 				}
// 			}

// 			// Average cost of not answering over all scenarios
// 			double sum_cost = 0;
// 			for(int j = 0; j < ins.calls.size(); ++j){
// 				std::vector<Call>
// future_calls(ins.calls[j].begin() +
// index_scenarios[j], ins.calls[j].end());
// FutureCall fc(data,env,*this, -1, -1,-1,-1,future_calls);
// fc.solve(); sum_cost += fc.obj;
// 			}

// 			double avg_cost = sum_cost / ins.calls.size();
// 			fmt::print("queue => {}\n",avg_cost);
// 			if(avg_cost < min_avg_cost){
// 				min_avg_cost = avg_cost;
// 				best_amb = best_h = -1;
// 			}

// 			//Answer call with best average cost
// 			if(best_amb == -1){
// 				fmt::print("queue\n");
// 				queue.push_back(index_call);
// 			}else{
// 				Ambulance& amb = ambulances[best_amb];
// 				std::cout << "Attended by " << amb << "\n";
// 				call.hospital = best_h;
// 				double min_time =
// travel.get_response_time(amb,call,time); 				double
// waiting_on_scene_i = time + min_time - call.time;
// double waiting_to_hospital_i = amb.answer_call(call, travel, ins, time,
// min_time, 0); //all returns to base 0 temporarily
// waiting_on_scene[call.id] = waiting_on_scene_i;
// waiting_on_scene_penalized[call.id] =
// waiting_on_scene_i*ins.penalties[call.priority];
// waiting_to_hospital[call.id] = waiting_to_hospital_i;
// which_ambulance[call.id] = best_amb;
// calls_end[call.id] = call.end; calls_attended++;
// 			}

// 			fmt::print("==========================\n");
// 		}else{
// 			Ambulance* ambulance = NULL;
// 			for(int a = 0; a < ambulances.size(); ++a){
// 				if(ambulances[a].arrival_time_at_f_last_trip ==
// time){ 					ambulance = &ambulances[a];
// 				}
// 			}
// 			std::cout << "Ambulance " << *ambulance << ":\n";
// 			Data data(this, NULL, ambulance);
// 			// std::cout << "Returning " << *ambulance << "\n";
// 			//for each base
// 			double min_avg_cost = GRB_INFINITY;
// 			int best_b = -1;
// 			for(int b = 0; b < ins.nb_bases; ++b){
// 				int sum_b = 0;
// 				for(int a = 0; a < data.types_amb; ++a){
// 					sum_b += data.A0_ab[a][b];
// 				}

// 				if(sum_b < data.cap_bases[b]){
// 					double sum_cost = 0;
// 					for(int j = 0; j < ins.calls.size();
// ++j){ 						std::vector<Call>
// future_calls(ins.calls[j].begin() +
// index_scenarios[j], ins.calls[j].end()); FutureAmbulance fc(data,env,*this,
// NULL, -1, b, future_calls); fc.solve();
// 						// if(j == 1 && index_call ==
// calls.size() - 1 && time > 8410){
// 						// 	GRBVar* vars =
// fc.model.getVars();
// 						// 	for(int i = 0; i <
// fc.model.get(GRB_IntAttr_NumVars); ++i){
// 						// 		auto name =
// vars[i].get(GRB_StringAttr_VarName);
// 						// 		auto val =
// vars[i].get(GRB_DoubleAttr_X);
// 						// 		if(val >
// g_params.EPS){
// 						//
// fmt::print("{} = {}\n",name, val);
// 						// 		}
// 						// 	}
// 						// 	delete[] vars;

// 						// 	//
// fc.model.write(fmt::format("model_b{}.lp", b));
// 						// 	std::cin.get();
// 						// }
// 						sum_cost += fc.obj;
// 					}

// 					//average cost over all scenarios
// 					double avg_cost = sum_cost /
// ins.calls.size(); 					fmt::print("Return to b
// => {}\n", avg_cost); 					if(avg_cost <
// min_avg_cost){ 						min_avg_cost =
// avg_cost; 						best_b = b;
// 					}
// 				}
// 			}
// 			int call_ind = -1;
// 			int best_h = -1;
// 			for(auto i: queue){
// 				if(can_answer(*ambulance,calls[i])){
// 					for(auto h: H[calls[i].id]){
// 						double sum_cost = 0;
// 						for(int j = 0; j <
// ins.calls.size(); ++j){
// std::vector<Call> future_calls(ins.calls[j].begin() +
// index_scenarios[j], ins.calls[j].end()); FutureAmbulance fc(data,env,*this,
// &calls[i], h, -1, future_calls); fc.solve();
// 							// if(j == 1 &&
// index_call == calls.size() - 1 && time > 8410){
// 							// 	GRBVar* vars =
// fc.model.getVars();
// 							// 	for(int i = 0; i
// < fc.model.get(GRB_IntAttr_NumVars); ++i){
// 							// 		auto
// name = vars[i].get(GRB_StringAttr_VarName);
// 							// 		auto val
// = vars[i].get(GRB_DoubleAttr_X);
// 							// 		if(val >
// g_params.EPS){
// 							//
// fmt::print("{} = {}\n",name, val);
// 							// 		}
// 							// 	}
// 							// 	delete[] vars;

// 							// 	//
// fc.model.write(fmt::format("model_b{}.lp", b));
// 							// 	std::cin.get();
// 							// }
// 							sum_cost += fc.obj;
// 						}

// 						double avg_cost = sum_cost /
// ins.calls.size();
// fmt::print("Answer call {} => {}\n", i, avg_cost);
// if(avg_cost < min_avg_cost){
// min_avg_cost = avg_cost;
// best_b =  -1; 							call_ind
// = i; 							best_h = h;
// 						}
// 					}
// 				}
// 			}

// 			if(best_b == -1){
// 				//assign ambulance to call_ind and best_h
// 				auto& call = calls[call_ind];
// 				std::cout << "Attending " << call << "\n";
// 				double min_time =
// travel.get_response_time(*ambulance, call, time);
// double waiting_on_scene_i = time + min_time - call.time;
// double waiting_to_hospital_i = ambulance->answer_call(call, travel, ins,
// time, 					min_time, 0); //returns to base
// 0 temporarily

// 				waiting_on_scene[call.id] = waiting_on_scene_i;
// 				waiting_on_scene_penalized[call.id] =
// waiting_on_scene_i*ins.penalties[call.priority];
// waiting_to_hospital[call.id] = waiting_to_hospital_i;
// which_ambulance[call.id] = ambulance->id;
// calls_end[call.id] = call.end;
// 				//remove call_ind from queue
// 				queue.erase(std::remove(queue.begin(),
// queue.end(), call_ind), 					queue.end());
// calls_attended++; 			}else{
// fmt::print("returning to {}\n", best_b);
// ambulance->base_location = ins.bases[best_b];
// 				ambulance->arrival_time_at_b_last_trip = time +
// travel.travel_time( ambulance->free_location, ambulance->base_location,
// *ambulance);
// 			}
// 		}
// 		fmt::print("==========================\n");
// 		set_next_event(event_call, index_call);
// 	}
// }

// void NonMiopycSolver::run(){
// 	std::vector<bool> is_allocated(calls.size(), false);
// 	std::vector<std::vector<double>> travel_times(calls.size(),
// 			std::vector<double>(ambulances.size(), GRB_INFINITY));
// 	int max_index = -1;
// 	std::vector<double> min_times(calls.size(), GRB_INFINITY);
// 	std::vector<int> index_amb(calls.size(), -1);
// 	std::vector<int> amb_type(calls.size(), -1);
// 	travel.set_forward(true);

// 	size_t i = 0;

// 	while(i < calls.size()){
// 		if(!is_allocated[i]){
// 			auto t0 = std::chrono::high_resolution_clock::now();
// 			if(i == static_cast<size_t>(max_index + 1)){
// 				for(size_t j = 0; j < ambulances.size(); ++j){
// 					auto & amb = ambulances[j];
// 					if(can_answer(amb, calls[i])){
// 						// travel_times[i][j] =
// travel.get_response_time(amb, calls[i], calls[i].time) *
// ins.penalty_matrix[amb.type][calls[i].priority];
// travel_times[i][j] = penalized_response_time(travel.get_response_time(amb,
// calls[i], calls[i].time), amb.type, calls[i].priority);
// 					}
// 				}

// 				compute_min_times(static_cast<int>(i),
// min_times, index_amb, amb_type); 				max_index = i;
// 			}

// 			while(!is_allocated[i]){
// 				auto& amb = ambulances[index_amb[i]];

// 				if(amb.arrival_time_at_f_last_trip <=
// calls[i].time){ 					is_allocated[i] = true;
// double waiting_on_scene_i = min_times[i];
// auto& call = calls[i]; 					double
// real_min_time = travel.get_response_time(amb, call, call.time);
// 					// double waiting_to_hospital_i =
// amb.answer_call(calls[i],travel,ins,calls[i].time,
// 					// 	travel_times[i][amb.id] /
// ins.penalty_matrix[amb.type][call.priority], nearest_base[i]);
// double waiting_to_hospital_i =
// amb.answer_call(calls[i],travel,ins,calls[i].time,
// real_min_time, nearest_base[i]);

// 					waiting_on_scene[i] =
// waiting_on_scene_i;
// waiting_on_scene_penalized[i] = waiting_on_scene_i*
// 						ins.penalties[calls[i].priority];
// 					waiting_to_hospital[i] =
// waiting_to_hospital_i;
// which_ambulance[i] = amb.id; calls_end[i] = calls[i].end;

// 					for(int k = i+1; k <= max_index; ++k){
// 						if(!is_allocated[k]){
// 							if(can_answer(amb,
// calls[k])){
// 								//
// travel_times[k][index_amb[i]] = travel.get_response_time(
// 								//
// amb,calls[k],calls[k].time) *
// ins.penalty_matrix[amb.type][calls[k].priority];
// 								travel_times[k][index_amb[i]]
// = penalized_response_time(travel.get_response_time(
// 									amb,calls[k],calls[k].time),
// amb.type, calls[k].priority);
// compute_min_times(k,min_times, index_amb, amb_type);
// 							}
// 						}
// 					}
// 				}else{
// 					int k = max_index + 1;
// 					bool cont = true;
// 					while(cont){
// 						if(static_cast<size_t>(k) >=
// calls.size()){ 							cont =
// false; 						}else{
// if(calls[k].time > amb.arrival_time_at_f_last_trip){
// cont = false; 							}else{
// if(!is_allocated[k]){
// for(size_t j = 0; j < ambulances.size(); ++j){ auto& amb_j = ambulances[j];
// if(can_answer(amb_j,calls[k])){
// 											//
// travel_times[k][j] = travel.get_response_time(amb_j,
// 											//
// calls[k], calls[k].time) * ins.penalty_matrix[amb_j.type][calls[k].priority];
// 											travel_times[k][j]
// = penalized_response_time(travel.get_response_time(amb_j,
// calls[k], calls[k].time), amb_j.type, calls[k].priority);
// 										}
// 									}
// 									compute_min_times(k,min_times,index_amb,amb_type);
// 								}
// 								k += 1;
// 							}
// 						}
// 					}

// 					max_index = k - 1;
// 					double maxmin = -1;
// 					int which_k = -1;
// 					double t1 =
// ins.penalties[calls[i].priority]*min_times[i];
// k = i + 1; 					if(static_cast<size_t>(i) <
// calls.size()
// - 1){ 						auto& amb_i =
// ambulances[index_amb[i]]; 						bool
// contd = calls[k].time <= amb_i.arrival_time_at_f_last_trip;

// 						while(contd){
// 							if(!is_allocated[k]){
// 								double t2 =
// ins.penalties[calls[k].priority]*min_times[k];

// 								if(abs(min_times[k]
// - travel_times[k][index_amb[i]]) < g_params.EPS &&
// t2 > maxmin){
// which_k = k;
//                                 	maxmin = t2;
// 								}
// 							}

// 							if(static_cast<size_t>(k+1)
// >= calls.size()){
// contd = false; 							}else
// if(calls[k+1].time > amb_i.arrival_time_at_f_last_trip){
// contd = false; 							}else{
// k += 1;
// 							}
// 						}
// 					}

// 					if(maxmin > t1){
// 						is_allocated[which_k] = true;
// 						auto& amb_i =
// ambulances[index_amb[i]]; 						double
// waiting_on_scene_i = min_times[which_k];
// auto& call_k = calls[which_k];
// double real_min_time = travel.get_response_time(amb_i, call_k, call_k.time);
// double waiting_to_hospital_i = amb_i.answer_call(calls[which_k],
// travel,ins, calls[which_k].time, real_min_time, nearest_base[which_k]);
// 						waiting_on_scene[which_k] =
// waiting_on_scene_i;
// waiting_on_scene_penalized[which_k] = waiting_on_scene_i*
// 							ins.penalties[calls[which_k].priority];
// 						waiting_to_hospital[which_k] =
// waiting_to_hospital_i;
// which_ambulance[which_k] = amb_i.id;
// calls_end[which_k] = calls[which_k].end;

// 						for(int k = i; k <= max_index;
// ++k){ if(!is_allocated[k]){
// if(can_answer(amb_i,calls[k])){
// 									//
// travel_times[k][index_amb[i]] = travel.get_response_time(
// 									//
// amb_i, calls[k], calls[k].time) *
// ins.penalty_matrix[amb_i.type][calls[k].priority];
// 									travel_times[k][index_amb[i]]
// = penalized_response_time(travel.get_response_time(
// amb_i, calls[k], calls[k].time), amb_i.type, calls[k].priority);
// 								}
// 								compute_min_times(k,min_times,index_amb,amb_type);
// 							}
// 						}
// 					}else{
// 						is_allocated[i] = true;
// 						auto& amb_i =
// ambulances[index_amb[i]]; 						double
// waiting_on_scene_i = min_times[i]; auto& call_i = calls[i];
// double real_min_time = travel.get_response_time(amb_i, call_i, call_i.time);
// double waiting_to_hospital_i =
// amb_i.answer_call(calls[i],travel,ins,calls[i].time,
// real_min_time, nearest_base[i]);

// 						waiting_on_scene[i] =
// waiting_on_scene_i;
// waiting_on_scene_penalized[i] = waiting_on_scene_i*
// 							ins.penalties[calls[i].priority];
// 						waiting_to_hospital[i] =
// waiting_to_hospital_i;
// which_ambulance[i] = amb_i.id;
// calls_end[i] = calls[i].end;

// 						for(int k = i; k <= max_index;
// ++k){ if(!is_allocated[k]){
// if(can_answer(amb_i,calls[k])){
// 									//
// travel_times[k][index_amb[i]] = travel.get_response_time(
// 									//
// amb_i,calls[k],calls[k].time) *
// ins.penalty_matrix[amb_i.type][calls[k].priority];
// 									travel_times[k][index_amb[i]]
// = penalized_response_time(travel.get_response_time(
// 										amb_i,calls[k],calls[k].time),
// amb_i.type, calls[k].priority);
// 								}
// 								compute_min_times(k,min_times,index_amb,amb_type);
// 							}
// 						}
// 					}
// 				}
// 			}
// 			auto dt = std::chrono::high_resolution_clock::now();
// 			run_times.push_back(std::chrono::duration_cast<chrono::microseconds>(dt
// - 				t0).count()); 		}else{
// 			++i;
// 		}
// 	}
// }

// void NonMiopycSolver::compute_min_times(int k, vector<double>& min_times,
// 	vector<int>& index_amb, vector<int> & amb_type){
// 	double min_time = GRB_INFINITY;
// 	int best_amb = -1;
// 	int best_type = -1;

// 	for(size_t j = 0; j < ambulances.size(); ++j){
// 		auto& amb = ambulances[j];
// 		double time_amb = GRB_INFINITY;
// 		if(can_answer(amb, calls[k])){
// 			if(amb.arrival_time_at_b_last_trip <= calls[k].time){
// 				time_amb = travel.travel_time(amb.base_location,
// calls[k].location, amb); 			}else
// if(amb.arrival_time_at_f_last_trip <= calls[k].time){
// Location current_location = travel.ambulance_position(amb,calls[k].time);
// time_amb = travel.travel_time(current_location,calls[k].location, amb);
// }else{ 				double time_free =
// amb.arrival_time_at_f_last_trip-calls[k].time;
// double time_free_to_call = travel.travel_time(amb.free_location,
// calls[k].location, amb); 				time_amb = time_free +
// time_free_to_call;
// 			}

// 			if((time_amb < min_time) || (abs(time_amb - min_time) <
// g_params.EPS && 				(amb.type > best_type))) {
// best_amb = j; 				best_type = amb.type;
// min_time = time_amb;
// 			}
// 		}
// 	}

// 	min_times[k] = min_time;
// 	index_amb[k] = best_amb;
// 	amb_type[k]  = best_type;
// }