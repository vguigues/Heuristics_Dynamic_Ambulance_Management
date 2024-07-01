#include "../include/data.h"

#include "../include/solver.h"
#include "../include/travel.h"

Data::Data(Call* a_call, Ambulance* a_ambulance, Instance& a_ins,
           Travel& a_travel, vector<int>& a_queue, vector<Call>& a_calls,
           vector<Ambulance>& a_ambulances, double a_time)
    : solver(NULL),
      ins(a_ins),
      travel(a_travel),
      call(a_call),
      ambulance(a_ambulance),
      queue(a_queue),
      calls(a_calls),
      ambulances(a_ambulances),
      time(a_time) {
  load_data();
}

Data::Data(Solver* solver, Call* call, Ambulance* ambulance)
    : solver(solver),
      ins(solver->ins),
      travel(solver->travel),
      call(call),
      ambulance(ambulance),
      queue(solver->queue),
      calls(solver->calls),
      ambulances(solver->ambulances),
      time(solver->time) {
  load_data();
}

void Data::load_data() {
  path = ins.path;
  if (call != NULL) {
    type_call0 = call->priority;
  }

  const bool debug = false;

  num_hosps = ins.nb_hospitals;
  num_bases = ins.nb_bases;
  num_locals = ins.nb_regions;
  types_call = ins.nb_priorities;
  types_amb = ins.nb_types_ambulance;
  num_ambs = ins.nb_ambulances;
  total_locals = num_locals + num_bases + num_hosps;

  num_times = 26;
  quantum = (20 * 1800) / (num_times - 2);
  max_cost = -GRB_INFINITY;

  hospitals = new int[num_hosps];
  bases = new int[num_bases];
  cap_bases = new int[num_bases];
  is_hospital = new bool[total_locals];
  is_base = new bool[total_locals];

  for (int l1 = 0; l1 < total_locals; ++l1) {
    is_base[l1] = is_hospital[l1] = false;
  }

  for (int b = 0; b < num_bases; ++b) {
    bases[b] = num_locals + b;
    is_base[bases[b]] = true;
  }

  for (int h = 0; h < num_hosps; ++h) {
    hospitals[h] = num_locals + num_bases + h;
    is_hospital[hospitals[h]] = true;
  }
  for (int b = 0; b < num_bases; ++b) {
    cap_bases[b] = ins.cap_bases[b];
  }

  A = ins.A;
  locations = ins.centers;
  locations.insert(locations.end(), ins.bases.begin(), ins.bases.end());
  locations.insert(locations.end(), ins.hospitals.begin(), ins.hospitals.end());
  // set_regions();
  if (call != NULL) {
    local_call0 = call->region;
  }

  H = new std::set<int>*[types_call];
  for (int c = 0; c < types_call; ++c) {
    H[c] = new std::set<int>[num_locals];
    for (int l = 0; l < num_locals; ++l) {
      H[c][l].insert(ins.nearest_hospital_to_region[l]);
    }
  }

  C = new int*[types_call];
  for (int c = 0; c < types_call; ++c) {
    C[c] = new int[num_locals];
    for (int l = 0; l < num_locals; ++l) {
      C[c][l] = 0;
    }
  }

  for (size_t i = 0; i < queue.size(); ++i) {
    auto future_call = calls[queue[i]];
    if (call == NULL || future_call != *call) {
      int location_index = -1;
      location_index = future_call.region;
      C[future_call.priority][location_index] += 1;
      if (debug) {
        if (C[future_call.priority][location_index] > 0) {
          fmt::print("C0_cl: c = {}, l = {}, val = {}\n", future_call.priority,
                     location_index, C[future_call.priority][location_index]);
        }
      }
    }
  }

  A0_ab = new int*[types_amb];
  for (int a = 0; a < types_amb; ++a) {
    A0_ab[a] = new int[num_bases];
    for (int b = 0; b < num_bases; ++b) {
      A0_ab[a][b] = 0;
    }
  }

  A0_ah = new int*[types_amb];
  for (int a = 0; a < types_amb; ++a) {
    A0_ah[a] = new int[num_hosps];
    for (int h = 0; h < num_hosps; ++h) {
      A0_ah[a][h] = 0;
    }
  }

  A0_alb = new int**[types_amb];
  for (int a = 0; a < types_amb; ++a) {
    A0_alb[a] = new int*[total_locals];
    for (int l = 0; l < total_locals; ++l) {
      A0_alb[a][l] = new int[num_bases];
      for (int b = 0; b < num_bases; ++b) {
        A0_alb[a][l][b] = 0;
      }
    }
  }
  A0_calh = new int***[types_call];
  for (int c = 0; c < types_call; ++c) {
    A0_calh[c] = new int**[types_amb];
    for (int a = 0; a < types_amb; ++a) {
      A0_calh[c][a] = new int*[total_locals];
      for (int l = 0; l < total_locals; ++l) {
        A0_calh[c][a][l] = new int[num_hosps];
        for (int h = 0; h < num_hosps; ++h) {
          A0_calh[c][a][l][h] = 0;
        }
      }
    }
  }

  A0_callh = new int****[types_call];
  for (int c = 0; c < types_call; ++c) {
    A0_callh[c] = new int***[types_amb];
    for (int a = 0; a < types_amb; ++a) {
      A0_callh[c][a] = new int**[total_locals];
      for (int l1 = 0; l1 < total_locals; ++l1) {
        A0_callh[c][a][l1] = new int*[num_locals];
        for (int l = 0; l < num_locals; ++l) {
          A0_callh[c][a][l1][l] = new int[num_hosps];
          for (int h = 0; h < num_hosps; ++h) {
            A0_callh[c][a][l1][l][h] = 0;
          }
        }
      }
    }
  }

  for (auto& amb : ambulances) {
    if (amb.arrival_time_at_f_last_trip > time + g_params.EPS) {
      Location origin_location = amb.last_origin_location;
      Location current_location = travel.position_between_origin_destination(
          origin_location, amb.call->location, amb.departure_time, time,
          amb.arrival_time_at_c_last_trip, amb);
      double travel_time =
          travel.travel_time(origin_location, amb.call->location, amb);
      if (current_location == amb.call->location) {
        current_location = travel.position_between_origin_destination(
            amb.call->location, ins.hospitals[amb.call->hospital],
            amb.departure_time + travel_time, time,
            amb.arrival_time_at_f_last_trip - amb.call->time_at_hospital, amb);
        int l1 = get_location_index(current_location);
        A0_calh[amb.call->priority][amb.type][l1][amb.call->hospital] += 1;
        if (debug) {
          fmt::print("A0_calh: c = {}, a = {}, l = {}, h = {}\n",
                     amb.call->priority, amb.type, l1, amb.call->hospital);
        }
      } else {
        int l1 = get_location_index(current_location);
        int l = amb.call->region;
        int h = amb.call->hospital;
        A0_callh[amb.call->priority][amb.type][l1][l][h] += 1;
        if (debug) {
          fmt::print("A0_callh, c = {}, a = {}, l1 = {}, l = {}, h = {}\n",
                     amb.call->priority, amb.type, l1, l, h);
        }
      }
    } else {
      if (amb.arrival_time_at_b_last_trip <= time) {
        int base_ind = -1;
        for (int b = 0; b < num_bases; ++b) {
          if (travel.lat_long_distance(amb.base_location, ins.bases[b]) <
              g_params.EPS) {
            base_ind = b;
          }
        }
        A0_ab[amb.type][base_ind] += 1;
        if (debug) {
          fmt::print("A0_ab, a = {}, b = {}\n", amb.type, base_ind);
        }
      }else if(time > amb.arrival_time_at_f_last_trip /*abs(time - amb.arrival_time_at_f_last_trip) > g_params.EPS*/){
        Location current_location = travel.ambulance_position(amb, time);
        int loc_ind = get_location_index(current_location);
        int base_ind = -1;
        for (int b = 0; b < num_bases; ++b) {
          if (travel.lat_long_distance(amb.base_location, ins.bases[b]) <
              g_params.EPS) {
            base_ind = b;
            break;
          }
        }
        // amb_nodes_map.push_back(std::make_pair(LOC,loc_ind));
        A0_alb[amb.type][loc_ind][base_ind] += 1;
        if (debug) {
          fmt::print("A0_alb, a = {}, l = {}, b = {}\n", amb.type, loc_ind,
                     base_ind);
        }
      } else if (abs(time - amb.arrival_time_at_f_last_trip) < g_params.EPS) {
        double min_d = GRB_INFINITY;
        int hosp_ind = -1;
        for (int h = 0; h < num_bases; ++h) {
          double d =
              travel.lat_long_distance(amb.free_location, ins.hospitals[h]);
          if (d < min_d) {
            min_d = d;
            hosp_ind = h;
          }
        }
        A0_ah[amb.type][hosp_ind] += 1;
        if (debug) {
          fmt::print("A0_ah, a = {}, h = {}\n", amb.type, hosp_ind);
        }
      }
    }
  }

  if (ambulance != NULL && call == NULL) {
    int sum = 0;
    for (int h = 0; h < num_hosps; ++h) {
      sum += A0_ah[ambulance->type][h];
    }
    if (sum < 1) {
      // double min_d = GRB_INFINITY;
      // int min_h = -1;
      // for(int h = 0; h < num_hosps; ++h){
      // 	double d = travel.lat_long_distance(ambulance->free_location,
      // ins.hospitals[h]); 	if(d < min_d){ 		min_d = d;
      // min_h = h;
      // 	}
      // }
      fmt::print(
          "ERROR: Ambulance {} returning at time {} but no ambulance in vector "
          "A0_ah.\n",
          ambulance->id, time);
      for (auto& amb : ambulances) {
        cout << amb << "\n";
      }
      cout << "Returning " << *ambulance << "\n";
      cin.get();
      // A0_ah[ambulance->type][min_h] += 1;
      // cin.get();
    } else if (sum > 1) {
      fmt::print(
          "Unlikely: Ambulance {} returning at time {} while {} ambulances "
          "also returning\n",
          ambulance->id, time, sum - 1);
      for (auto& amb : ambulances) {
        cout << amb << "\n";
      }
    }
  }
  set_times();

  std::set<int>**** Lt_tab = new std::set<int>***[num_times];
  for (int tao = 0; tao < num_times; ++tao) {
    Lt_tab[tao] = new std::set<int>**[num_times];
    for (int t = 0; t < num_times; ++t) {
      Lt_tab[tao][t] = new std::set<int>*[types_amb];
      for (int a = 0; a < types_amb; ++a) {
        Lt_tab[tao][t][a] = new std::set<int>[num_bases];
      }
    }
  }

  for (int t = 0; t < num_times; ++t) {
    for (int a = 0; a < types_amb; ++a) {
      for (int b = 0; b < num_bases; ++b) {
        for (int h = 0; h < num_hosps; ++h) {
          int forecast_l1 = L(t, a, hospitals[h], b);
          if (forecast_l1 != bases[b]) {
            Lt_tab[0][t][a][b].insert(forecast_l1);
          }
        }
      }
    }
  }

  for (int tao = 1; tao < num_times; ++tao) {
    for (int t = 0; t < num_times; ++t) {
      for (int a = 0; a < types_amb; ++a) {
        for (int b = 0; b < num_bases; ++b) {
          for (auto l1 : Lt_tab[tao - 1][t][a][b]) {
            int forecast_l1 = L(t, a, l1, b);
            if (forecast_l1 != bases[b]) {
              Lt_tab[tao][t][a][b].insert(forecast_l1);
            }
          }
        }
      }
    }
  }

  L_tab = new std::set<int>**[num_times];
  for (int t = 0; t < num_times; ++t) {
    L_tab[t] = new std::set<int>*[types_amb];
    for (int a = 0; a < types_amb; ++a) {
      L_tab[t][a] = new std::set<int>[num_bases];
      for (int b = 0; b < num_bases; ++b) {
        for (int tao = 0; tao < num_times; ++tao) {
          for (auto l1 : Lt_tab[tao][t][a][b]) {
            L_tab[t][a][b].insert(l1);
          }
        }
      }
    }
  }

  for (int tao = 0; tao < num_times; ++tao) {
    for (int t = 0; t < num_times; ++t) {
      for (int a = 0; a < types_amb; ++a) {
        delete[] Lt_tab[tao][t][a];
      }
      delete[] Lt_tab[tao][t];
    }
    delete[] Lt_tab[tao];
  }
  delete[] Lt_tab;

  // if(debug){
  // 	for(int b = 0; b < num_bases; ++b){
  // 		fmt::print("L_tab b = {}: ", b);
  // 		for(auto l: L_tab[0][0][b]){
  // 			fmt::print("{} ", l);
  // 		}
  // 		fmt::print("\n");
  // 	}
  // }

  lambda = new int**[num_times - 1];
  for (int t = 0; t < num_times - 1; ++t) {
    lambda[t] = new int*[types_call];
    for (int c = 0; c < types_call; ++c) {
      lambda[t][c] = new int[num_locals];
      for (int l = 0; l < num_locals; ++l) {
        lambda[t][c][l] = 0;
      }
    }
  }

  for (size_t i = 0; i < calls.size(); ++i) {
    auto& call = calls[i];
    int t = get_time_slot(call.time);
    if (debug) {
      std::cout << "num_times = " << num_times << ", ";
      std::cout << "T = " << t << " (" << call.time << "), ";
      std::cout << "P = " << call.priority << ", ";
      std::cout << "R = " << call.region << "\n";
    }
    if (t != -1) {
      lambda[t][call.priority][call.region] += 1;
    }
  }

  if (debug) {
    if (call != NULL) {
      cout << "Input call " << *call << "\n";
      fmt::print("c0, l0 = {} {}\n", type_call0, local_call0);
    } else if (ambulance != NULL) {
      cout << "Input ambulance " << *ambulance << "\n";
    }
    for (auto& call : calls) {
      cout << call << "\n";
    }
    fmt::print("current time = {}, quantum = {}\n", time, quantum);
    fmt::print(
        "============================= End debug Data "
        "=============================\n");
    cin.get();
  }

  relax["x0_abh"] = true;
  relax["x0_ahh"] = true;
  relax["x0_albh"] = true;
  relax["y0_ahb"] = true;
  relax["y0_ahclh"] = true;
  relax["xt_cablh"] = true;
  relax["xt_cahlh"] = true;
  relax["xt_calblh"] = true;
  relax["yt_ahb"] = true;
  relax["ct_cl"] = true;
  relax["at_ab"] = true;
  relax["at_alb"] = true;
}

void Data::set_times() {
  tao = new double*****[num_times];
  for (int t = 0; t < num_times; ++t) {
    tao[t] = new double****[types_call];
    for (int c = 0; c < types_call; ++c) {
      double time_on_scene = 600 + ins.delay_scene_p[c];
      double time_at_hospital = 600 - ins.delay_hosp_p[c];
      tao[t][c] = new double***[types_amb];
      for (int a = 0; a < types_amb; ++a) {
        tao[t][c][a] = new double**[total_locals];
        for (int l1 = 0; l1 < total_locals; ++l1) {
          tao[t][c][a][l1] = new double*[num_locals];
          for (int l = 0; l < num_locals; ++l) {
            tao[t][c][a][l1][l] = new double[num_hosps];
            for (int h = 0; h < num_hosps; ++h) {
              tao[t][c][a][l1][l][h] =
                  travel.travel_time(locations[l1], locations[l],
                                     ambulances[0]) +
                  time_on_scene +
                  travel.travel_time(locations[l], ins.hospitals[h],
                                     ambulances[0]) +
                  time_at_hospital;
              if (tao[t][c][a][l1][l][h] > max_cost) {
                max_cost = tao[t][c][a][l1][l][h];
              }
            }
          }
        }
      }
    }
  }

  tao_0 = new double***[types_call];
  for (int c = 0; c < types_call; ++c) {
    tao_0[c] = new double**[types_amb];
    for (int a = 0; a < types_amb; ++a) {
      tao_0[c][a] = new double*[total_locals];
      for (int l1 = 0; l1 < total_locals; ++l1) {
        tao_0[c][a][l1] = new double[num_hosps];
        for (int h = 0; h < num_hosps; ++h) {
          // double time_on_scene = 600 + ins.delay_scene_p[c];
          double time_at_hospital = 600 - ins.delay_hosp_p[c];
          tao_0[c][a][l1][h] =
              travel.travel_time(locations[l1], ins.hospitals[h],
                                 ambulances[0]) +
              time_at_hospital;
        }
      }
    }
  }
}

int Data::L(int t, int a, int l1, int b) {
  Location base;
  if (b < num_bases) {
    base = ins.bases[b];
  } else {
    base = locations[b];
  }
  double time_l1_b = travel.travel_time(locations[l1], base, ambulances[0]);
  if (time_l1_b <= quantum) {
    return (b < num_bases) ? bases[b] : b;
  } else {
    double travel_time = travel.travel_time(locations[l1], base, ambulances[0]);
    Location current_location = travel.position_between_origin_destination(
        locations[l1], base, time + t * quantum, time + (t + 1) * quantum,
        time + t * quantum + travel_time, ambulances[0]);
    return get_location_index(current_location);
  }
}

int Data::get_time_slot(double a_time) {
  if (a_time <= time) {
    return -1;
  }
  int t = (int)((a_time - time) / quantum);
  return t;
}

int Data::get_location_index(Location& location) {
  // for(int b = 0; b < num_bases; ++b){
  // 	if(travel.lat_long_distance(location, ins.bases[b]) < g_params.EPS){
  // 		return bases[b];
  // 	}
  // }

  // for(int h = 0; h < num_hosps; ++h){
  // 	if(travel.lat_long_distance(location, ins.hospitals[h]) < g_params.EPS){
  // 		return hospitals[h];
  // 	}
  // }

  double min_dist = GRB_INFINITY;
  double ind_l = -1;
  for (int l = 0; l < num_locals; ++l) {
    double dist_to_region = travel.lat_long_distance(location, ins.centers[l]);
    if (dist_to_region < min_dist) {
      min_dist = dist_to_region;
      ind_l = l;
    }
  }

  return ind_l;
}

Data::~Data() {
  for (int c = 0; c < types_call; ++c) {
    delete[] H[c];
  }
  for (int c = 0; c < types_call; ++c) {
    delete[] C[c];
  }
  for (int a = 0; a < types_amb; ++a) {
    delete[] A0_ab[a];
    delete[] A0_ah[a];
  }
  for (int a = 0; a < types_amb; ++a) {
    for (int l = 0; l < total_locals; ++l) {
      delete[] A0_alb[a][l];
    }
    delete[] A0_alb[a];
  }
  for (int c = 0; c < types_call; ++c) {
    for (int a = 0; a < types_amb; ++a) {
      for (int l = 0; l < total_locals; ++l) {
        delete[] A0_calh[c][a][l];
      }
      delete[] A0_calh[c][a];
    }
    delete[] A0_calh[c];
  }
  for (int c = 0; c < types_call; ++c) {
    for (int a = 0; a < types_amb; ++a) {
      for (int l1 = 0; l1 < total_locals; ++l1) {
        for (int l = 0; l < num_locals; ++l) {
          delete[] A0_callh[c][a][l1][l];
        }
        delete[] A0_callh[c][a][l1];
      }
      delete[] A0_callh[c][a];
    }
    delete[] A0_callh[c];
  }
  for (int t = 0; t < num_times - 1; ++t) {
    for (int c = 0; c < types_call; ++c) {
      delete[] lambda[t][c];
    }
    delete[] lambda[t];
  }

  for (int t = 0; t < num_times; ++t) {
    for (int c = 0; c < types_call; ++c) {
      for (int a = 0; a < types_amb; ++a) {
        for (int l1 = 0; l1 < total_locals; ++l1) {
          for (int l = 0; l < num_locals; ++l) {
            delete[] tao[t][c][a][l1][l];
          }
          delete[] tao[t][c][a][l1];
        }
        delete[] tao[t][c][a];
      }
      delete[] tao[t][c];
    }
    delete[] tao[t];
  }
  for (int c = 0; c < types_call; ++c) {
    for (int a = 0; a < types_amb; ++a) {
      for (int l1 = 0; l1 < total_locals; ++l1) {
        delete[] tao_0[c][a][l1];
      }
      delete[] tao_0[c][a];
    }
    delete[] tao_0[c];
  }

  for (int t = 0; t < num_times; ++t) {
    for (int a = 0; a < types_amb; ++a) {
      delete[] L_tab[t][a];
    }
    delete[] L_tab[t];
  }
  delete[] L_tab;

  delete[] H;
  delete[] A0_alb;
  delete[] A0_calh;
  delete[] A0_callh;
  delete[] lambda;
  delete[] tao;
  delete[] tao_0;
  delete[] C;
  delete[] A0_ab;
  delete[] A0_ah;
  delete[] hospitals;
  delete[] bases;
  delete[] is_hospital;
  delete[] cap_bases;
  delete[] is_base;
}

// Data::Data(EventSimulator& sim, Call* call, Ambulance* ambulance):
// ins(sim.ins), sim(sim), 	debug(g_params.debug), call(call),
// ambulance(ambulance){

// 	path = ins.path;
// 	if(call != NULL)
// 		type_call0 = call->priority;

// 	num_hosps = ins.nb_hospitals;
// 	num_bases = ins.nb_bases;
// 	num_locals = g_params.n_regions;
// 	types_call = ins.nb_priorities;
// 	types_amb = ins.nb_types_ambulance;
// 	num_ambs = ins.nb_ambulances;
// 	total_locals = num_locals+num_bases+num_hosps;

// 	num_times = g_params.n_time_slots;
// 	quantum = (sim.last_time-sim.time)/(num_times-2);

// 	hospitals = new int[num_hosps];
// 	bases = new int[num_bases];
// 	cap_bases = new int[num_bases];
// 	is_hospital = new bool[total_locals];
// 	is_base = new bool[total_locals];

// 	for(int l1 = 0; l1 < total_locals; ++l1){
// 		is_base[l1] = is_hospital[l1] = false;
// 	}

// 	for(int b = 0; b < num_bases; ++b){
// 		bases[b] = num_locals+b;
// 		is_base[bases[b]] = true;
// 	}

// 	for(int h = 0; h < num_hosps; ++h){
// 		hospitals[h] = num_locals+num_bases+h;
// 		is_hospital[hospitals[h]] = true;
// 	}

// 	for(int b = 0; b < num_bases; ++b){
// 		cap_bases[b] = ins.cap_bases[b];
// 	}

// 	A = ins.A;

// 	set_regions();
// 	if(call != NULL){
// 		local_call0 = get_location_index(call->location);
// 	}

// 	H = new std::set<int>*[types_call];
// 	for(int c = 0; c < types_call; ++c){
// 		H[c] = new std::set<int>[num_locals];
// 		for(int l = 0; l < num_locals; ++l){
// 			std::vector<std::pair<double, int>> dists;
// 			for(int h = 0; h < num_hosps; ++h){
// 				double dist = sim.norm(ins.hospitals[h],
// locations[l]);
// dists.push_back(std::make_pair(dist,h));
// 			}
// 			std::sort(dists.begin(), dists.end());
// 			for(int k = 0; k <
// std::min(g_params.n_nearest_hospitals,num_hosps); ++k)
// 				H[c][l].insert(dists[k].second);
// 		}
// 	}

// 	C = new int*[types_call];
// 	for(int c = 0; c < types_call; ++c){
// 		C[c] = new int[num_locals];
// 		for(int l = 0; l < num_locals; ++l){
// 			C[c][l] = 0;
// 		}
// 	}

// 	for(auto call_queue: sim.queue){
// 		if(call == NULL  || *call_queue != *call){
// 			int location_index = -1;
// 			location_index =
// get_location_index(call_queue->location);
// 			C[call_queue->priority][location_index] += 1;
// 		}
// 	}

// 	A0_ab = new int*[types_amb];
// 	for(int a = 0; a < types_amb; ++a){
// 		A0_ab[a] = new int[num_bases];
// 		for(int b = 0; b < num_bases; ++b){
// 			A0_ab[a][b] = 0;
// 		}
// 	}

// 	A0_ah = new int*[types_amb];
// 	for(int a = 0; a < types_amb; ++a){
// 		A0_ah[a] = new int[num_hosps];
// 		for(int h = 0; h < num_hosps; ++h){
// 			A0_ah[a][h] = 0;
// 		}
// 	}

// 	A0_alb = new int**[types_amb];
// 	for(int a = 0; a < types_amb; ++a){
// 		A0_alb[a] = new int*[total_locals];
// 		for(int l = 0; l < total_locals; ++l){
// 			A0_alb[a][l] = new int[num_bases];
// 			for(int b = 0; b < num_bases; ++b){
// 				A0_alb[a][l][b] = 0;
// 			}
// 		}
// 	}
// 	A0_calh = new int***[types_call];
// 	for(int c = 0; c < types_call; ++c){
// 		A0_calh[c] = new int**[types_amb];
// 		for(int a = 0; a < types_amb; ++a){
// 			A0_calh[c][a] = new int*[total_locals];
// 			for(int l = 0; l < total_locals; ++l){
// 				A0_calh[c][a][l] = new int[num_hosps];
// 				for(int h = 0; h < num_hosps; ++h){
// 					A0_calh[c][a][l][h] = 0;
// 				}
// 			}
// 		}
// 	}

// 	A0_callh = new int****[types_call];
// 	for(int c = 0; c < types_call; ++c){
// 		A0_callh[c] = new int***[types_amb];
// 		for(int a = 0; a < types_amb; ++a){
// 			A0_callh[c][a] = new int**[total_locals];
// 			for(int l1 = 0; l1 < total_locals; ++l1){
// 				A0_callh[c][a][l1] = new int*[num_locals];
// 				for(int l = 0; l < num_locals; ++l){
// 					A0_callh[c][a][l1][l] = new
// int[num_hosps]; 					for(int h = 0; h <
// num_hosps; ++h){
// A0_callh[c][a][l1][l][h] = 0;
// 					}
// 				}
// 			}
// 		}
// 	}

// 	for(auto amb: sim.ambulances){
// 		if(ambulance != NULL && amb.id == ambulance->id){
// 			amb_nodes_map.push_back(std::make_pair(HOSP, -1));
// 			continue;
// 		}

// 		std::cout << "Amb " << amb.id << " t" << amb.type << ": ";

// if(amb.busy){
// 	std::cout << "busy at ";
// 	Location origin_location = amb.last_origin_location;
// 	Location current_location = sim.ambulance_busy_position(amb, sim.time);
// 	int location_index = get_location_index(current_location);
// 	int hosp_ind = -1;
// 	for(int h = 0; h < num_hosps; ++h){
// 		if(sim.norm(amb.hospital_location, ins.hospitals[h]) < EPS){
// 			hosp_ind = h;
// 			break;
// 		}
// 	}
// 			double departure_time = amb.departure_time;
// 			double time_to_call = sim.travel_time(amb.trips.back(),
// amb.call->location, 				amb.speed);
// std::cout << sim.show_location(current_location) << ", index = ";
// std::cout << location_index << ". Going to: ";
// if(departure_time + time_to_call <= sim.time){
// std::cout << "hosp " << sim.show_location(amb.hospital_location)
// << ", from "; 				std::cout <<
// sim.show_location(amb.last_origin_location) << ",
// "; 				std::cout << "h = " << hosp_ind << "\n";
// 				A0_calh[amb.call->priority][amb.type][location_index][hosp_ind]
// += 1; 			}else{ 				int
// call_location_index = get_location_index(amb.call->location);
// std::cout << "call " << sim.show_location(amb.call->location) << ", from ";
// std::cout << sim.show_location(amb.last_origin_location);
// std::cout << ", index = " << call_location_index << " and hosp ";
// std::cout << sim.show_location(amb.hospital_location) << ", h = ";
// std::cout << hosp_ind
// << "\n";
// 				A0_callh[amb.call->priority][amb.type][location_index][call_location_index]
// 					[hosp_ind] += 1;
// 			}
// 			amb_nodes_map.push_back(std::make_pair(-1,-1));
// 		}else{
// if(amb.arrival_time_at_b_last_trip <= sim.time){
// 	int base_ind = -1;
// 	for(int b = 0; b < num_bases; ++b){
// 		if(sim.norm(amb.base_location, ins.bases[b]) < EPS){
// 			base_ind = b;
// 		}
// 	}
// 	std::cout << "at base ";
// 	std::cout << sim.show_location(amb.base_location);
// 	std::cout << ", b = " << base_ind << "\n";
// 	A0_ab[amb.type][base_ind] += 1;
// 	amb_nodes_map.push_back(std::make_pair(BASE,base_ind));
// }else if(abs(sim.time - amb.arrival_time_at_h_last_trip) > EPS){
// 	Location current_location = sim.ambulance_position(amb, sim.time);
// 	int loc_ind = get_location_index(current_location);
// 	std::cout << "returning, at ";
// 	std::cout << sim.show_location(current_location) << ", index = ";
// 	std::cout << loc_ind << ", to b = ";
// 	int base_ind = -1;
// 	for(int b = 0; b < num_bases; ++b){
// 		if(sim.norm(amb.base_location, ins.bases[b]) < EPS){
// 			base_ind = b;
// 			break;
// 		}
// 	}
// 	amb_nodes_map.push_back(std::make_pair(LOC,loc_ind));
// 	A0_alb[amb.type][loc_ind][base_ind] += 1;
// 	std::cout << base_ind << "\n";
// }
// 		}
// 	}
// 	if(ambulance != NULL){
// 		int hosp_ind = -1;
// 		for(int h = 0; h < num_hosps; ++h){
// 			if(sim.norm(ambulance->hospital_location,
// sim.ins.hospitals[h]) < EPS){ 				hosp_ind = h;
// break;
// 			}
// 		}
// 		if(hosp_ind != -1){
// 			std::cout << "Amb " << ambulance->id << " t" <<
// ambulance->type << "at h = "; 			std::cout << hosp_ind <<
// "\n"; 			A0_ah[ambulance->type][hosp_ind] += 1;
// 			amb_nodes_map[ambulance->id].second = hosp_ind;
// 		}else{
// 			std::cout << "ERROR: ambulance is not at any
// hospital\n"; 			exit(1);
// 		}
// 	}

// 	// tao
// 	set_times();
// 	// lt_tab
// std::set<int>**** Lt_tab = new std::set<int>***[num_times];
// for(int tao = 0; tao < num_times; ++tao){
// 	Lt_tab[tao] = new std::set<int>**[num_times];
// 	for(int t = 0;t < num_times; ++t){
// 		Lt_tab[tao][t] = new std::set<int>*[types_amb];
// 		for(int a = 0; a < types_amb; ++a){
// 			Lt_tab[tao][t][a] =  new std::set<int>[num_bases];
// 		}
// 	}
// }

// for(int t = 0; t < num_times; ++t){
// 	for(int a = 0; a < types_amb; ++a){
// 		for(int b = 0; b < num_bases; ++b){
// 			for(int h = 0; h < num_hosps; ++h){
// 				int forecast_l1 = L(t, a, hospitals[h], b);
// 				if(forecast_l1 != bases[b]){
// 					Lt_tab[0][t][a][b].insert(forecast_l1);
// 				}
// 			}
// 		}
// 	}
// }

// for(int tao = 1; tao < num_times; ++tao){
// 	for(int t = 0; t < num_times; ++t){
// 		for(int a = 0; a < types_amb; ++a){
// 			for(int b = 0; b < num_bases; ++b){
// 				for(auto l1: Lt_tab[tao-1][t][a][b]){
// 					int forecast_l1 = L(t, a, l1, b);
// 					if(forecast_l1 != bases[b]){
// 						Lt_tab[tao][t][a][b].insert(forecast_l1);
// 					}
// 				}
// 			}
// 		}
// 	}
// }

// L_tab = new std::set<int>**[num_times];
// for(int t = 0; t < num_times; ++t){
// 	L_tab[t] = new std::set<int>*[types_amb];
// 	for(int a = 0; a < types_amb; ++a){
// 		L_tab[t][a] = new std::set<int>[num_bases];
// 		for(int b = 0; b < num_bases; ++b){
// 			for(int tao = 0; tao < num_times; ++tao){
// 				for(auto l1: Lt_tab[tao][t][a][b])
// 					L_tab[t][a][b].insert(l1);
// 			}
// 		}
// 	}
// }

// lambda = new int**[num_times-1];
// for(int t = 0; t < num_times-1; ++t){
// 	lambda[t] = new int*[types_call];
// 	for(int c = 0; c < types_call; ++c){
// 		lambda[t][c] = new int[num_locals];
// 		for(int l = 0; l < num_locals; ++l){
// 			lambda[t][c][l] = 0;
// 		}
// 	}
// }

// for(auto call: sim.calls){
// 	int t = get_time_slot(call.time);
// 	if(t != -1){
// 		point_t bg_call_location(call.location.first,
// call.location.second); 		int l_index = -1;
// 		// std::cout << "location: "<< call.location.first << " " <<
// call.location.second << "\n"; 		for(int l = 0; l < num_locals;
// ++l){ 			if(bg::covered_by(bg_call_location,regions[l])){
// l_index = l; 				break;
// 			}
// 		}
// 		// std::cout << "num_time " << num_times << "\n";
// 		// std::cout << "time " << t << " " << call.time << "\n";
// 		// std::cout << "priority " << call.priority << "\n";
// 		// std::cout << "Region " << l_index << "\n";
// 		lambda[t][call.priority][l_index] += 1;
// 	}
// }

// relax["x0_abh"] = true;
// relax["x0_ahh"] = true;
// relax["x0_albh"] = true;
// relax["y0_ahb"] = true;
// relax["y0_ahclh"] = true;
// relax["xt_cablh"] = true;
// relax["xt_cahlh"] = true;
// relax["xt_calblh"] = true;
// relax["yt_ahb"] = true;
// relax["ct_cl"] = true;
// relax["at_ab"] = true;
// relax["at_alb"] = true;
// }

// void Data::set_regions(){
// 	int n_regions = g_params.n_regions;
// 	int L = (int) sqrt(n_regions);

// 	double width = ins.x_max - ins.x_min;
// 	double height = ins.y_max - ins.y_min;

// 	double grid_width = width / L;
// 	double grid_height = height / L;

// 	double x_left_origin = ins.x_min;
// 	double x_right_origin = ins.x_min + grid_width;
// 	double y_bottom_origin = ins.y_min;
// 	double y_top_origin = ins.y_min + grid_height;

// 	double x_left,x_right;
// 	int k = 0;
// 	for(int i = 0; i < L; ++i){
// 		x_left = x_left_origin;
// 		x_right = x_right_origin;
// 		for(int j = 0; j < L; ++j){
// 			polygon_t poly;
// 			bg::append(poly.outer(),
// point_t(x_left,y_bottom_origin)); 			bg::append(poly.outer(),
// point_t(x_left,y_top_origin)); 			bg::append(poly.outer(),
// point_t(x_right,y_top_origin)); 			bg::append(poly.outer(),
// point_t(x_right,y_bottom_origin)); 			bg::append(poly.outer(),
// point_t(x_left,y_bottom_origin)); 			regions.push_back(poly);

// 			point_t center;
// 			bg::centroid(poly,center);
// 			locations.push_back(std::make_pair(bg::get<0>(center),
// 				bg::get<1>(center)));
// 			// std::cout << "Polygon:\n";
// 			// for(auto it = boost::begin(bg::exterior_ring(poly));
// 			// 	it != boost::end(bg::exterior_ring(poly));
// ++it){
// 			//     double x = bg::get<0>(*it);
// 			//     double y = bg::get<1>(*it);
// 			//     std::cout << "\t" << x << " " << y << "\n";
// 			//     //use the coordinates...
// 			// }
// 			// std::cout << "Center: " << locations[k].first << " ";
// 			// std::cout << locations[k].second << "\n";
// 			// std::cin.get();
// 			++k;
// 			x_left += grid_width;
// 			x_right += grid_width;
// 		}
// 		y_top_origin += grid_height;
// 		y_bottom_origin += grid_height;
// 	}

// 	for(auto base: ins.bases){
// 		locations.push_back(base);
// 	}

// 	for(auto hospital: ins.hospitals){
// 		locations.push_back(hospital);
// 	}
// }