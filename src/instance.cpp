#include "../include/instance.h"

#include "../include/ambulance.h"
#include "../include/call.h"
#include "../include/osrm_helper.h"
#include "../include/travel.h"

// Constructor for tests generating scenarios for all 7 days;
Instance::Instance() : path(""), travel(false) {
  nb_times = 48;
  slot_duration = 0.5 * 3600;
  nb_days = 7;
  nb_scenarios = nb_days * g_params.n_scenarios;
  nb_hospitals = g_params.n_hospitals;
  nb_bases = g_params.n_bases;
  nb_cleaning_bases = 4;
  int nb_land_types = 5;
  nb_regions = 76;
  nb_ambulances = g_params.n_ambulances;

  if (g_params.amb_setup == "rj") {
    nb_types_ambulance = 3;
    nb_priorities = 3;
    A.push_back({0});
    A.push_back({0, 1});
    A.push_back({0, 1, 2});
    penalties = {4, 2, 1};
  } else if (g_params.amb_setup == "us") {
    nb_types_ambulance = 2;
    nb_priorities = 4;
    double theta_1 = 1, theta_2 = 4;
    // type 0: (ALS preferred, HOT), type 1: (ALS preferred, COLD), type 2: ([no
    // pref], HOT), type 3: ([no pref], COLD)
    nb_priorities = 4;       // Urgent, Non-urgent
    nb_types_ambulance = 2;  // 0 ALS, 1 BLS
    penalty_matrix = {{0, 0, theta_1, theta_1}, {theta_2, theta_2, 0, 0}};

    penalties = {theta_2, theta_1, theta_2, theta_1};

    A.push_back({0, 1});
    A.push_back({0, 1});
    A.push_back({0, 1});
    A.push_back({0, 1});
  } else {
    fmt::print("ERROR: Unknown setup: {}\n", g_params.amb_setup);
    exit(1);
  }

  load_lambda();
  load_lambda_bases();
  load_regions();
  time_horizon = {36, 37, 38, 39};

  if (g_params.best_base) {
    modify_lambdas();
  }

  auto hosp_arq = ifstream("hospitals.txt", ios::in);
  int total_h;
  hosp_arq >> total_h;
  if (nb_hospitals > total_h) {
    fmt::print("nb_hospitals {} > total hospitals {}\n", nb_hospitals, total_h);
    exit(1);
  }

  for (int h = 0; h < nb_hospitals; ++h) {
    double lat, longi;
    hosp_arq >> lat >> longi;
    hospitals.push_back(make_pair(lat, longi));
  }
  hosp_arq.close();
  nearest_hospital_to_region = vector<int>(nb_regions, -1);
  auto regions_hospitals_table = travel.table_in_out(centers, hospitals);
  for (int r = 0; r < nb_regions; ++r) {
    double min_time = GRB_INFINITY;
    for (int h = 0; h < nb_hospitals; ++h) {
      // fmt::print("{} {} to {} {} = {}\n",
      // 	centers[r].first, centers[r].second,
      // 	hospitals[h].first, hospitals[h].second,
      // 	regions_hospitals_table[r][h]);
      if (regions_hospitals_table[r][h] < min_time) {
        min_time = regions_hospitals_table[r][h];
        nearest_hospital_to_region[r] = h;
      }
    }
  }

  fmt::print("Loaded Hospitals\n");

  auto base_arq = ifstream("bases.txt", ios::in);
  int total_b;
  base_arq >> total_b;
  if (nb_bases > total_b) {
    fmt::print("nb_bases {} > total bases {}\n", nb_bases, total_b);
    exit(1);
  }
  for (int b = 0; b < nb_bases; ++b) {
    double lat, longi;
    base_arq >> lat >> longi;
    bases.push_back(make_pair(lat, longi));
    cap_bases.push_back(nb_ambulances);
  }
  base_arq.close();
  // std::random_device dev;
  // std::mt19937 gen(dev());
  std::default_random_engine gen(600);
  std::default_random_engine gen_times;
  auto samples_arq = ifstream(
      fmt::format("{}/samples.dat", g_params.generator_folder), ios::in);
  samples = std::vector<std::vector<Location>>(
      nb_regions, std::vector<Location>(100, null_location));
  for (int i = 0; i < nb_regions * 100; ++i) {
    int r, s;
    double lat, longi;
    samples_arq >> r >> s >> longi >> lat;
    samples[r][s].first = lat;
    samples[r][s].second = longi;
  }

  nearest_base_to_region = std::vector<int>(nb_regions, -1);

  for (int r = 0; r < nb_regions; ++r) {
    double min_d = GRB_INFINITY;
    for (int b = 0; b < nb_bases; ++b) {
      double d = travel.lat_long_distance(centers[r], bases[b]);
      if (d < min_d) {
        min_d = d;
        nearest_base_to_region[r] = b;
      }
    }
  }

  nearest_base_to_hospital = std::vector<int>(nb_hospitals, -1);

  for (int h = 0; h < nb_hospitals; ++h) {
    double min_d = GRB_INFINITY;
    for (int b = 0; b < nb_bases; ++b) {
      double d = travel.lat_long_distance(hospitals[h], bases[b]);
      if (d < min_d) {
        min_d = d;
        nearest_base_to_hospital[h] = b;
      }
    }
  }

  nearest_region_to_base = std::vector<int>(nb_bases, -1);
  for (int b = 0; b < nb_bases; ++b) {
    double min_d = GRB_INFINITY;
    for (int r = 0; r < nb_regions; ++r) {
      double d = travel.lat_long_distance(bases[b], centers[r]);
      if (d < min_d) {
        min_d = d;
        nearest_region_to_base[b] = r;
      }
    }
  }
  auto &nbr = nearest_base_to_region;
  vector<set<int>> regions_by_base(bases.size(), set<int>());
  for (int b = 0; b < nb_bases; ++b) {
    for (int r = 0; r < nb_regions; ++r) {
      if (b == nbr[r]) {
        regions_by_base[b].insert(r);
        // lambda_base[b][t_ind][c] += ins.lambda(t,4,r,c);
      }
    }
  }

  samples_arq.close();
  boost::filesystem::create_directories("daily_scenarios");

  // Tests for a friday night
  ofstream daily_scenarios(
      fmt::format("daily_scenarios/{}_s{}_a{}.txt", g_params.instance_type,
                  g_params.n_scenarios, nb_ambulances),
      std::ios::out);
  // std::scientific(daily_scenarios);
  total_scenarios = max(5, g_params.n_scenarios);
  delay_scene_p = vector<double>(nb_priorities, 0.0);
  delay_hosp_p = vector<double>(nb_priorities, 0.0);
  g0 = 4;

  for (int s = 0; s < total_scenarios; ++s) {
    std::vector<Call> scenario;
    int id = 0;
    // 2h == 40, 3h == 44, 4h == 48
    for (auto t : time_horizon) {
      for (int r = 0; r < nb_regions; ++r) {
        for (int p = 0; p < nb_priorities; ++p) {
          double param = get_param(t, 4, r, p);
          poisson_distribution<int> dist(param);
          int nb_calls = dist(gen);
          double slot_fraction = slot_duration / nb_calls;
          for (int k = 0; k < nb_calls; ++k) {
            double begin_slot = t * slot_duration + 1;
            double end_slot = begin_slot + (k + 1) * slot_fraction - 1;
            std::uniform_int_distribution<int> int_dist(begin_slot, end_slot);
            std::uniform_int_distribution<int> sample_rand(0, 99);
            std::uniform_real_distribution<double> u_scene(0, 1);
            std::uniform_real_distribution<double> u_hospital(0, 1);
            std::exponential_distribution<double> scene_rand(
                1 / (-penalties[p] * 200 * log(1 - u_scene(gen_times))));
            std::exponential_distribution<double> hosp_rand(
                1 / (-penalties[p] * 200 * log(1 - u_hospital(gen_times))));
            double time = int_dist(gen);
            Location location = samples[r][sample_rand(gen)];
            if (delay_scene_p[p] == 0.0) {
              delay_scene_p[p] = scene_rand(gen_times);
            }
            if (delay_hosp_p[p] == 0.0) {
              delay_hosp_p[p] = min(600.0, hosp_rand(gen_times));
            }

            // double time_on_scene = 600 + delay_scene_p[p];
            // double time_at_hospital = 600 - delay_hosp_p[p];
            double time_on_scene = 600;
            double time_at_hospital = 600;
            double cleaning_time = 2400;
            int min_h = nearest_hospital_to_region[r];
            int cb = 0;
            scenario.push_back(Call(id++, time, location, min_h, cb, p, true,
                                    false, time_on_scene, time_at_hospital,
                                    cleaning_time));
            scenario.back().region = r;
            scenario.back().disc_time = t;
          }
        }
      }
    }
    // fmt::print("Delay scene {}, delay hosp {}\n", delay_scene_p,
    // delay_hosp_p);
    std::sort(scenario.begin(), scenario.end());
    for (size_t i = 0; i < scenario.size(); ++i) {
      scenario[i].id = i;
    }
    daily_scenarios << scenario.size() << "\n";
    for (auto &call : scenario) {
      // Time - Region index - Priority - Day - Time on scene - Lat - Long -
      // Time at hospital - TimeCleaningBase - Cleaning needed - Hospital needed
      // - index hospital
      daily_scenarios << call.time << " " << call.region << " ";
      daily_scenarios << call.priority << " " << 4 << " ";
      daily_scenarios << call.time_on_scene << " ";
      daily_scenarios << call.location.first << " " << call.location.second
                      << " ";
      daily_scenarios << call.time_at_hospital << " ";
      daily_scenarios << call.cleaning_time << " ";
      daily_scenarios << call.clean_needed << " " << call.hosp_needed << " ";
      daily_scenarios << call.hospital << " " << call.cleaning << "\n";
      // cout << call.time_on_scene << "\t" << call.time_at_hospital << "\n";
    }
    calls.push_back(scenario);
    nb_calls.push_back(scenario.size());
  }
  daily_scenarios.close();
  fmt::print("Loaded Calls, {} scenarios\n", calls.size());
  // print_scenario_stats();

  // std::uniform_int_distribution<int> amb_dist(0, nb_bases-1);
  amb_count_by_type = vector<int>(nb_types_ambulance, 0);
  // if nb_ambulances is not multiple of nb_types_ambulance, this map creates
  // more ambulances of less advanced types.
  vector<int> types_map;
  if (g_params.amb_setup == "rj") {
    types_map = {2, 1, 0};
  } else if (g_params.amb_setup == "us") {
    types_map = {1, 0};
  }
  for (int a = 0; a < nb_ambulances; ++a) {
    // int b = amb_dist(gen);
    Location base_location = bases[a % nb_bases];
    Location hospital_location = null_location;
    int amb_type = types_map[a % nb_types_ambulance];
    ++amb_count_by_type[amb_type];
    ambulances.push_back(Ambulance(a, base_location, hospital_location, -1, -1,
                                   0, (amb_type), 16.6667));
  }
}

// Instance::Instance(): path(g_params.instance), travel(false){
// 	nb_types_ambulance = 3;
// 	nb_priorities = 3;
// 	slot_duration = 3600;
// 	path = g_params.instance;
// 	nb_times = g_params.n_time_slots;
// 	nb_regions = g_params.n_regions;
// 	nb_hospitals = g_params.n_hospitals;
// 	nb_bases = g_params.n_bases;
// 	nb_cleaning_bases = g_params.n_cleaning_bases;
// 	nb_ambulances = g_params.n_ambulances;

// 	bool is_euclidian = path == simulated_rect_1 || path ==
// simulated_rect_2;

// 	penalties = {4,2,1};

// 	A.push_back({0});
// 	A.push_back({0,1});
// 	A.push_back({0,1,2});

// 	if(is_euclidian){
// 		load_euclidian();
// 	}else{
// 		load_real_data();
// 	}

// 	for(int a = 0; a < nb_ambulances; ++a){
// 		// int b = amb_dist(gen);
// 		Location base_location = bases[a % nb_bases];
// 		Location hospital_location = null_location;
// 		const double speed = 0.078; //approximate speed needed to cross
// square in 30m 		ambulances.push_back(Ambulance(a, base_location,
// hospital_location, -1, -1, 0, 			(a %
// nb_types_ambulance), speed));
// 	}
// 	// fmt::print("Ambulances:\n");
// 	// for(auto& amb: ambulances){
// 	// 	std::cout << amb << " " << amb.base_location.first << " " <<
// amb.base_location.second;
// 	// 	std::cout << "\n";
// 	// }
// }

void Instance::load_euclidian() {
  travel.euclidian = true;

  nb_hospitals = 4;
  nb_bases = 4;
  nb_cleaning_bases = 1;

  hospitals.push_back(make_pair(50, 0));
  hospitals.push_back(make_pair(100, 50));
  hospitals.push_back(make_pair(50, 100));
  hospitals.push_back(make_pair(0, 50));

  bases.push_back(make_pair(0, 0));
  bases.push_back(make_pair(100, 0));
  bases.push_back(make_pair(100, 100));
  bases.push_back(make_pair(0, 100));

  cleaning_bases.push_back(make_pair(50, 50));

  default_random_engine gen;
  auto calls_arq = ifstream(path, ios::in);
  int scenario_size = 0;
  double time, time_on_scene, lat, longi, time_at_hospital, time_at_cleaning;
  int priority, cleaning_needed, hospital_needed, index_hospital;
  int k = 0;
  while (k++ < g_params.n_scenarios) {
    calls_arq >> scenario_size;
    nb_calls.push_back(scenario_size);
    std::vector<Call> scenario;
    scenario.reserve(scenario_size);
    int id = 0;
    for (int i = 0; i < scenario_size; ++i) {
      calls_arq >> time >> lat >> longi >> priority >> cleaning_needed;
      calls_arq >> time_at_cleaning >> index_hospital >> time_on_scene;
      calls_arq >> hospital_needed >> time_at_hospital;
      priority -= 1;
      index_hospital -= 1;
      uniform_int_distribution<int> int_dist((time - 1) * slot_duration,
                                             time * slot_duration);

      Location call_location = make_pair(lat, longi);
      double min_d = GRB_INFINITY;
      for (int h = 0; h < nb_hospitals; ++h) {
        double d = travel.euclidian_distance(call_location, hospitals[h]);
        if (d < min_d) {
          min_d = d;
          index_hospital = h;
        }
      }

      int index_cleaning = 0;

      if (priority == 0) {
        std::uniform_int_distribution<int> cb_rand(0, 4);
        if (cb_rand(gen) == 0) {
          cleaning_needed = 1;
        }
      }

      scenario.push_back(Call(
          id++, time * slot_duration, call_location, index_hospital,
          index_cleaning, priority, hospital_needed, cleaning_needed,
          0.333 * slot_duration, 0.333 * slot_duration, 0.666 * slot_duration));
    }
    // fmt::print("scenario size {}\n", scenario.size());
    // scenario = std::vector<Call>(scenario.begin(), scenario.begin()+10);
    calls.push_back(scenario);
  }
  calls_arq.close();
  nb_scenarios = min(static_cast<int>(calls.size()), g_params.n_scenarios);
  calls = std::vector<std::vector<Call>>(calls.begin(),
                                         calls.begin() + nb_scenarios);
}

void Instance::load_real_data() {
  auto hosp_arq = ifstream("hospitals.txt", ios::in);
  int total_h;
  hosp_arq >> total_h;
  if (nb_hospitals > total_h) {
    fmt::print("nb_hospitals {} > total hospitals {}\n", nb_hospitals, total_h);
    exit(1);
  }
  for (int h = 0; h < nb_hospitals; ++h) {
    double lat, longi;
    hosp_arq >> lat >> longi;
    hospitals.push_back(make_pair(lat, longi));
  }
  hosp_arq.close();
  // fmt::print("Hospitals:\n");
  // for(auto h: hospitals){
  // 	fmt::print("{} {}\n",h.first, h.second);
  // }

  auto base_arq = ifstream("bases.txt", ios::in);
  int total_b;
  base_arq >> total_b;
  if (nb_bases > total_b) {
    fmt::print("nb_bases {} > total bases {}\n", nb_bases, total_b);
    exit(1);
  }
  for (int b = 0; b < nb_bases; ++b) {
    double lat, longi;
    base_arq >> lat >> longi;
    bases.push_back(make_pair(lat, longi));
    cap_bases.push_back((nb_ambulances / nb_bases) + 4);
  }
  base_arq.close();
  // fmt::print("Bases:\n");
  // for(auto h: bases){
  // 	fmt::print("{} {}\n",h.first, h.second);
  // }

  auto cleaning_arq = ifstream("cleaning.txt", ios::in);
  int total_cb;
  cleaning_arq >> total_cb;
  if (nb_cleaning_bases > total_cb) {
    fmt::print("cleaning_bases {} > total cleaning bases {}\n",
               nb_cleaning_bases, total_cb);
    exit(1);
  }
  for (int c = 0; c < nb_cleaning_bases; ++c) {
    double lat, longi;
    cleaning_arq >> lat >> longi;
    cleaning_bases.push_back(std::make_pair(lat, longi));
  }
  cleaning_arq.close();
  // fmt::print("Cleaning Bases:\n");
  // for(auto h: cleaning_bases){
  // 	fmt::print("{} {}\n",h.first, h.second);
  // }

  auto samples_arq = ifstream(
      fmt::format("{}/samples.dat", g_params.generator_folder), ios::in);
  std::vector<std::vector<Location>> samples(
      nb_regions, std::vector<Location>(100, null_location));
  for (int i = 0; i < nb_regions * 100; ++i) {
    int r, s;
    double lat, longi;
    samples_arq >> r >> s >> longi >> lat;
    samples[r][s].first = lat;
    samples[r][s].second = longi;
  }
  samples_arq.close();

  std::default_random_engine gen;
  auto calls_arq = ifstream(path, ios::in);
  int scenario_size = 0;
  // Time - Region index - Priority - Day - Time on scene - Lat - Long - Time at
  // hospital - TimeCleaningBase - Cleaning needed - Hospital needed - index
  // hospital
  double time, time_on_scene, lat, longi, time_at_hospital, time_at_cleaning;
  int region, priority, day, cleaning_needed, hospital_needed, index_hospital,
      index_cleaning;
  while (calls_arq >> scenario_size) {
    nb_calls.push_back(scenario_size);
    std::vector<Call> scenario;
    scenario.reserve(scenario_size);
    int id = 0;
    for (int i = 0; i < scenario_size; ++i) {
      calls_arq >> time >> region >> priority >> day >> time_on_scene;
      calls_arq >> lat >> longi >> time_at_hospital >> time_at_cleaning;
      calls_arq >> cleaning_needed >> hospital_needed >> index_hospital;
      // fmt::print("line: {} {} {} {} {} {} {} {} {} {} {}
      // {}\n",time,region,priority,
      // 	day,time_on_scene,lat,longi,time_at_hospital,time_at_cleaning,
      // cleaning_needed, 	hospital_needed, index_hospital);
      region -= 1;
      priority -= 1;
      std::uniform_int_distribution<int> int_dist((time - 1) * slot_duration,
                                                  time * slot_duration);
      std::uniform_int_distribution<int> sample_rand(0, 99);
      double min_d = GRB_INFINITY;
      int sample_id = sample_rand(gen);
      auto call_location = samples[region][sample_id];
      index_cleaning = -1;
      Location previous_location =
          (hospital_needed) ? hospitals[index_hospital] : call_location;
      for (int c = 0; c < nb_cleaning_bases; ++c) {
        double travel_time = travel.travel_time(
            previous_location, cleaning_bases[c], ambulances[0]);
        if (travel_time < min_d) {
          min_d = travel_time;
          index_cleaning = c;
        }
      }
      min_d = GRB_INFINITY;
      index_hospital = -1;
      if (g_params.h_random_hospital) {
        std::uniform_int_distribution<int> h_rand(0, nb_hospitals - 1);
        index_hospital = h_rand(gen);
      } else {
        for (int h = 0; h < nb_hospitals; ++h) {
          double travel_time =
              travel.travel_time(call_location, hospitals[h], ambulances[0]);
          if (travel_time < min_d) {
            min_d = travel_time;
            index_hospital = h;
          }
        }
      }

      if (priority == 0) {
        std::uniform_int_distribution<int> cb_rand(0, 4);
        if (cb_rand(gen) == 0) {
          cleaning_needed = 1;
        }
      }

      scenario.push_back(Call(
          id++, time * slot_duration, call_location, index_hospital,
          index_cleaning, priority, hospital_needed, cleaning_needed,
          0.333 * slot_duration, 0.333 * slot_duration, 0.666 * slot_duration));
      auto &call = scenario.back();
      call.region = region;
      // std::cout << call << " " << call.location.first << " " <<
      // call.location.second << "\n";
    }
    // fmt::print("=============================\n");
    calls.push_back(scenario);
  }
  calls_arq.close();
  nb_scenarios = min(static_cast<int>(calls.size()), g_params.n_scenarios);
  calls = std::vector<std::vector<Call>>(calls.begin(),
                                         calls.begin() + nb_scenarios);
}

Instance::Instance(std::string path) : path(path), travel(true) {
  std::ifstream arq(path, std::ios::in);
  std::cout << path << "\n";
  arq >> nb_scenarios;
  // arq >> nb_calls;
  arq >> nb_hospitals;
  arq >> nb_bases;
  arq >> nb_cleaning_bases;
  arq >> nb_types_ambulance;
  arq >> nb_ambulances;
  arq >> nb_priorities;
  arq >> x_min >> x_max >> y_min >> y_max;

  nb_calls = vector<int>(nb_scenarios, 0);
  std::cout << nb_scenarios << " ";
  std::cout << nb_hospitals << " " << nb_bases << " ";
  std::cout << nb_cleaning_bases << " " << nb_types_ambulance << " ";
  std::cout << nb_ambulances << " " << nb_priorities << "\n";
  std::cout << x_min << " " << x_max << " " << y_min << " " << y_max << "\n";

  A = vector<set<int>>(nb_priorities, set<int>());
  // int aux = 0;
  std::string aux_str;
  // arq.ignore();
  // for(int c = 0; c < nb_priorities; ++c){
  // 	std::getline(arq, aux_str);
  // 	std::istringstream ss(aux_str);
  // 	while(ss >> aux){
  // 		A[c].insert(aux);
  // 	}
  // }

  // for(int c = 0; c < nb_priorities; ++c){
  // 	std::cout << c << ": ";
  // 	for(auto a: A[c]){
  // 		std::cout << a << " ";
  // 	}
  // 	std::cout << "\n";
  // }
  // std::cin.get();

  double time, lat, longi, time_on_scene, time_at_hospital, cleaning_time;
  int priority, hospital, cleaning, hosp_needed, clean_needed;

  penalties.reserve(nb_priorities);
  for (int i = 0; i < nb_priorities; ++i) {
    double aux;
    arq >> aux;
    penalties.push_back(aux);
    std::cout << penalties.back() << " ";
  }
  std::cout << "\n";

  penalties = {4, 2, 1};

  for (int k = 0; k < nb_scenarios; ++k) {
    calls.push_back(vector<Call>());
    arq >> nb_calls[k];
    for (int i = 0; i < nb_calls[k]; ++i) {
      arq >> time >> lat >> longi >> hospital >> cleaning >> priority;
      arq >> hosp_needed >> clean_needed >> time_on_scene;
      arq >> time_at_hospital >> cleaning_time;
      Location call_location(lat, longi);
      // auto time_equipment = time_on_scene;
      calls.back().push_back(Call(
          i, time, call_location, hospital, cleaning, priority, hosp_needed,
          clean_needed, time_on_scene, time_at_hospital, cleaning_time));
      std::cout << time << " " << lat << " " << longi << " " << hospital << " ";
      std::cout << cleaning << " " << priority << " " << hosp_needed << " ";
      std::cout << clean_needed << " " << time_on_scene << " "
                << time_at_hospital;
      std::cout << " " << cleaning_time << "\n";
    }
  }

  for (int h = 0; h < nb_hospitals; ++h) {
    arq >> lat >> longi;
    hospitals.push_back(std::make_pair(lat, longi));
    std::cout << lat << " " << longi << "\n";
  }

  int cap;
  for (int b = 0; b < nb_bases; ++b) {
    arq >> lat >> longi >> cap;
    bases.push_back(std::make_pair(lat, longi));
    cap_bases.push_back(cap);
    std::cout << lat << " " << longi << " " << cap << "\n";
  }

  for (int cb = 0; cb < nb_cleaning_bases; ++cb) {
    arq >> lat >> longi;
    cleaning_bases.push_back(std::make_pair(lat, longi));
    std::cout << lat << " " << longi << "\n";
  }

  int type;
  double speed;
  for (int a = 0; a < nb_ambulances; ++a) {
    arq >> lat >> longi >> type >> speed;
    Location base_location(lat, longi);
    Location hospital_location = null_location;
    ambulances.push_back(
        Ambulance(a, base_location, hospital_location, -1, -1, 0, type, speed));
    std::cout << lat << " " << longi << " " << type << " " << speed << "\n";
  }

  if (path == "instances/no_myopic7.3.txt") {
    ambulances[2].free_location = std::make_pair(5, 10);
    double t_to_b0 = travel.travel_time(ambulances[2].free_location, bases[0],
                                        ambulances[2]);
    double t_to_b1 = travel.travel_time(ambulances[2].free_location, bases[1],
                                        ambulances[2]);
    if (t_to_b0 < t_to_b1) {
      ambulances[2].base_location = bases[0];
    } else {
      ambulances[2].base_location = bases[1];
    }
    ambulances[2].arrival_time_at_f_last_trip = 1;
    ambulances[2].arrival_time_at_b_last_trip =
        1 + travel.travel_time(ambulances[2].free_location,
                               ambulances[2].base_location, ambulances[2]);
    std::cout << ambulances[2] << "\n";
    ambulances[2].set_new_point(0, hospitals[0], TripType::TO_HOSPITAL);
    ambulances[2].set_new_point(1, hospitals[0], TripType::AT_HOSPITAL);
    if (t_to_b0 < t_to_b1) {
      ambulances[2].set_new_point(1 + t_to_b0, ambulances[2].base_location,
                                  TripType::TO_BASE);
    } else {
      ambulances[2].set_new_point(1 + t_to_b1, ambulances[2].base_location,
                                  TripType::TO_BASE);
    }

  } else if (path == "instances/no_myopic7.4.txt") {
    ambulances[1].free_location = hospitals[0];
    ambulances[1].base_location = bases[0];
    ambulances[1].arrival_time_at_f_last_trip = 3;
    double t_to_b0 =
        travel.travel_time(ambulances[1].free_location,
                           ambulances[1].base_location, ambulances[1]);
    ambulances[1].arrival_time_at_b_last_trip = 3 + t_to_b0;

    ambulances[1].set_new_point(0, hospitals[0], TripType::TO_HOSPITAL);
    ambulances[1].set_new_point(3, hospitals[0], TripType::AT_HOSPITAL);
    ambulances[1].set_new_point(3 + t_to_b0, ambulances[1].base_location,
                                TripType::TO_BASE);
  }

  // exit(1);
}

void Instance::print_scenario_stats() {
  for (size_t s = 0; s < calls.size(); ++s) {
    auto &scenario = calls[s];
    int g = s % 7;
    vector<double> lam_r_agg(nb_regions, 0.0);
    for (int r = 0; r < nb_regions; ++r) {
      double sum = 0.0;
      for (int t = 0; t < nb_times; ++t) {
        for (int p = 0; p < nb_priorities; ++p) {
          sum += lambda(t, g, r, p);
        }
      }
      lam_r_agg[r] = sum;
    }
    double total_lam = accumulate(lam_r_agg.begin(), lam_r_agg.end(), 0.0);
    double total_calls = scenario.size();
    vector<int> calls_by_region(nb_regions, 0);
    for (auto &call : scenario) {
      ++calls_by_region[call.region];
    }
    double avg_diff = 0.0;
    vector<double> diffs;
    for (int r = 0; r < nb_regions; ++r) {
      double frac_lam = lam_r_agg[r] / total_lam;
      double frac_calls = calls_by_region[r] / total_calls;
      diffs.push_back(abs(frac_calls - frac_lam));
      avg_diff += diffs.back();
      fmt::print("r = {}, pct_lam = {:.4f}, pct_calls = {:.4f}, diff {:.4f}\n",
                 r, frac_lam * 100, frac_calls * 100, diffs.back());
    }
    fmt::print("Scenario {}: avg_diff = {:.4f}, max_diff = {:.4f}\n", s,
               avg_diff / nb_regions, *max_element(diffs.begin(), diffs.end()));
    cin.get();
  }
}

double Instance::get_param(int t, int g, int r, int p) {
  const string ins_type = g_params.instance_type;
  if (ins_type == "one_priority") {
    double sum = 0.0;
    for (int p = 0; p < 4; ++p) {
      sum += lambda(t, g, r, p);
    }
    return sum;
  } else if (ins_type == "two_priorities") {
    if (p >= 2) {
      return (lambda(t, g, r, 2) + lambda(t, g, r, 3)) / 2;
    } else {
      return (lambda(t, g, r, 0) + lambda(t, g, r, 1)) / 2;
    }
  } else if (ins_type == "four_priorities") {
    return lambda(t, g, r, p);
  }

  fmt::print("ERROR: Unknown instance_type = {}\n", ins_type);
  exit(1);
}

double Instance::survival_function(double time) { return exp(time); }

double Instance::get_max_service_time(vector<Call> &calls) {
  double max_time = -GRB_INFINITY;
  for (auto &call : calls) {
    double avg_center_time = 0;
    for (auto &center : centers) {
      avg_center_time +=
          travel.travel_time(center, call.location, ambulances[0]);
    }
    avg_center_time /= centers.size();
  }

  return max_time;
}

void Instance::load_lambda() {
  const ulong samu_priorities = 3;
  bool is_setup_us = g_params.amb_setup == "us";
  xt::xarray<double>::shape_type shape = {(ulong)nb_times, (ulong)nb_days,
                                          (ulong)nb_regions,
                                          (ulong)samu_priorities + is_setup_us};
  lambda = xt::zeros<double>(shape);
  fmt::print("Shape p: {}\n", lambda.shape(3));
  auto lambda_arq = ifstream(fmt::format("{}/sorted_empirical_estimation.dat",
                                         g_params.generator_folder),
                             ios::in);
  ulong file_T, file_G, file_R, file_P;
  int nb_obs_total;
  lambda_arq >> file_T >> file_G >> file_R >> file_P >> nb_obs_total;
  if (file_T != (ulong)nb_times || file_G != (ulong)nb_days ||
      file_R != (ulong)nb_regions || file_P != samu_priorities) {
    fmt::print(
        "Lambdas shape ({},{},{},{}) different of settings ({},{},{},{})\n",
        file_T, file_G, file_R, file_P, nb_times, nb_days, nb_regions, 3);
    exit(1);
  }
  int t, g, r, p;
  double val;
  for (int i = 0; i < nb_times * nb_days * nb_regions * 3; ++i) {
    lambda_arq >> t >> g >> r >> p >> val;
    lambda(t, g, r, p + is_setup_us) = val;
  }
  lambda_arq.close();

  if (g_params.amb_setup == "us") {
    for (int t = 0; t < nb_times; ++t) {
      for (int g = 0; g < nb_days; ++g) {
        for (int r = 0; r < nb_regions; ++r) {
          lambda(t, g, r, 0) = lambda(t, g, r, 1);
        }
      }
    }
  }
  fmt::print("Loaded Lambda\n");
}

void Instance::load_lambda_bases() {
  const ulong max_bases = 34;
  const ulong samu_priorities = 3;
  bool is_setup_us = g_params.amb_setup == "us";
  xt::xarray<double>::shape_type shape = {(ulong)nb_times, (ulong)nb_days,
                                          max_bases,
                                          (ulong)samu_priorities + is_setup_us};
  lambda_bases = xt::zeros<double>(shape);

  auto lambda_arq =
      ifstream(fmt::format("calibration/Bases_voronoi/empirical_estimation.dat",
                           g_params.generator_folder),
               ios::in);
  ulong file_T, file_G, file_R, file_P;
  int nb_obs_total;
  lambda_arq >> file_T >> file_G >> file_R >> file_P >> nb_obs_total;
  if (file_T != (ulong)nb_times || file_G != (ulong)nb_days ||
      file_R != max_bases || file_P != (ulong)samu_priorities) {
    fmt::print(
        "Lambdas shape ({},{},{},{}) different of settings ({},{},{},{})\n",
        file_T, file_G, file_R, file_P, nb_times, nb_days, bases.size(), 3);
    exit(1);
  }
  int t, g, r, p;
  double val;
  for (int i = 0; i < nb_times * nb_days * 34 * 3; ++i) {
    lambda_arq >> t >> g >> r >> p >> val;
    lambda_bases(t, g, r, p + is_setup_us) = val;
  }
  lambda_arq.close();

  if (g_params.amb_setup == "us") {
    for (int t = 0; t < nb_times; ++t) {
      for (int g = 0; g < nb_days; ++g) {
        for (ulong r = 0; r < max_bases; ++r) {
          lambda(t, g, r, 0) = lambda(t, g, r, 1);
        }
      }
    }
  }
  fmt::print("Loaded Lambda Bases\n");
}

void Instance::load_regions() {
  int nb_land_types = 0;
  auto neighbors_arq = ifstream(
      fmt::format("{}/sorted_neighbors.dat", g_params.generator_folder),
      ios::in);
  centers = vector<Location>(nb_regions, null_location);
  neighbors = vector<vector<int>>(nb_regions, vector<int>());
  adj_matrix =
      vector<vector<bool>>(nb_regions, vector<bool>(nb_regions, false));

  string aux_str;
  int num_regions, num_regressors;
  std::getline(neighbors_arq, aux_str);
  std::istringstream ss(aux_str);
  ss >> num_regions >> num_regressors;
  if (num_regions != nb_regions) {
    fmt::print("Error: neighbors file has {} regions, but nb_regions = {}\n",
               num_regions, nb_regions);
    cin.get();
  }
  while (true) {
    int ind, terrain_type, s;
    double lat, longi, dist;
    std::getline(neighbors_arq, aux_str);
    if (aux_str == "END") {
      break;
    }
    std::istringstream ss(aux_str);
    ss >> ind >> lat >> longi;  // >> terrain_type;
    centers[ind] = make_pair(lat, longi);
    double regressor;
    for (int j = 0; j < nb_land_types; ++j) {
      ss >> regressor;
    }
    while (ss >> s >> dist) {
      neighbors[ind].push_back(s);
      adj_matrix[ind][s] = true;
    }
  }
  neighbors_arq.close();

  fmt::print("Loaded Regions\n");
}

void Instance::read_preparedness() {
  ulong nb, nt, nta, n0, n1, n2;
  string file_name = fmt::format("preps/preparedness_{}_{}a_{}b.txt",
                                 g_params.amb_setup, nb_ambulances, nb_bases);
  // fmt::print("file_name = {}\n",file_name);
  ifstream prep_file(file_name, std::ios::in);

  prep_file >> nb >> nt >> nta;
  if (nb != static_cast<ulong>(nb_bases)) {
    fmt::print(
        "Error: number of bases in preparedness file {} not compatible "
        "with number of bases in the simulator: {}",
        nb, nb_bases);
    exit(1);
  }
  if (nt != time_horizon.size()) {
    fmt::print(
        "Error: time horizon size in preparedness file {} not "
        "compatible with time horizon size in the simulator: {}",
        nt, time_horizon.size());
    exit(1);
  }

  if (nta != static_cast<ulong>(nb_types_ambulance)) {
    fmt::print(
        "Error: preparedness file has {} types of ambulance but "
        "instance has {}.\n",
        nta, nb_types_ambulance);
    exit(1);
  }

  if (nta == 2) {
    prep_file >> n0 >> n1;
    if (n0 != (ulong)amb_count_by_type[0] ||
        n1 != (ulong)amb_count_by_type[1]) {
      fmt::print(
          "Error: Supply in preparedness file ({}, {}) not compatible "
          "with supply of simulator ({}, {})\n",
          n0, n1, amb_count_by_type[0], amb_count_by_type[1]);
      exit(1);
    }
    preparedness = xt::zeros<double>({nb, nt, n0 + 1, n1 + 1});
    double val;
    int b, t, supply0, supply1;
    for (ulong k = 0; k < nb * nt * (n0 + 1) * (n1 + 1); ++k) {
      prep_file >> b >> t >> supply0 >> supply1 >> val;
      preparedness(b, t, supply0, supply1) = val;
    }
  } else if (nta == 3) {
    prep_file >> n0 >> n1 >> n2;
    if (n0 != (ulong)amb_count_by_type[0] ||
        n1 != (ulong)amb_count_by_type[1] ||
        n2 != (ulong)amb_count_by_type[2]) {
      fmt::print(
          "Error: Supply in preparedness file ({}, {}, {}) not "
          "compatible with supply of simulator ({}, {}, {})\n",
          n0, n1, n2, amb_count_by_type[0], amb_count_by_type[1],
          amb_count_by_type[2]);
      exit(1);
    }
    preparedness = xt::zeros<double>({nb, nt, n0 + 1, n1 + 1, n2 + 1});
    double val;
    int b, t, supply0, supply1, supply2;
    for (ulong k = 0; k < nb * nt * (n0 + 1) * (n1 + 1) * (n2 + 1); ++k) {
      prep_file >> b >> t >> supply0 >> supply1 >> supply2 >> val;
      preparedness(b, t, supply0, supply1, supply2) = val;
    }
  }
  prep_file.close();
}

void Instance::modify_lambdas() {
  vector<int> west = {22, 23, 32, 33};
  vector<int> south = {29, 30, 39, 40};

  // At the first time slot, south lambdas are multiplied by 10.
  for (int g = 0; g < nb_days; ++g) {
    for (int r : south) {
      for (int p = 0; p < nb_priorities; ++p) {
        lambda(time_horizon[0], g, r, p) *= 10;
      }
    }
    for (int r : west) {
      for (int p = 0; p < nb_priorities; ++p) {
        lambda(time_horizon[0], g, r, p) = 0;
      }
    }
  }

  // At the last time slot, west lambdas are multiplied by 10
  for (int g = 0; g < nb_days; ++g) {
    for (int r : west) {
      for (int p = 0; p < nb_priorities; ++p) {
        lambda(time_horizon.back(), g, r, p) *= 10;
      }
    }
    for (int r : south) {
      for (int p = 0; p < nb_priorities; ++p) {
        lambda(time_horizon[0], g, r, p) = 0;
      }
    }
  }
}

Instance::~Instance() {}

// auto scenarios_arq = ofstream(fmt::format("{}/scenarios.txt",
// g_params.generator_folder), 	ios::out);

// for(int s = 0; s < nb_scenarios; ++s){
// 	int id = 0;
// 	for(int t = 0; t < nb_times; ++t){
// 		std::uniform_int_distribution<int> int_dist(t*slot_duration,
// (t+1)*slot_duration); 		std::uniform_int_distribution<int>
// sample_rand(0,99); 		for(int g = 0; g < 1; ++g){
// for(int r = 0; r < nb_regions;
// ++r){ 				for(int p = 0; p < nb_priorities; ++p){
// poisson_distribution<int> dist(lambda(t,g,r,p));
// int this_nb_calls = dist(gen); nb_calls[s] += this_nb_calls;
// for(int c = 0; c < this_nb_calls; ++c){
// double time = int_dist(gen);
// 						// Location location =
// samples[r][sample_rand(gen)];
// Location location = centers[r]; double time_on_scene = 0.1*slot_duration;
// double time_at_hospital = 0.2*slot_duration;
// double cleaning_time = 1*slot_duration;
// double min_time = GRB_INFINITY;
// int min_h = nearest_hospital[r];
// int cb = 0;
// calls[s].push_back(Call(id++, time, location, min_h, cb, p, true, false,
// time_on_scene, time_at_hospital, cleaning_time)); calls[s].back().region = r;
// calls[s].back().disc_time = t;
// 					}
// 				}
// 			}
// 		}
// 	}
// 	std::sort(calls[s].begin(), calls[s].end());
// 	for(size_t i = 0; i < calls[s].size(); ++i){
// 		calls[s][i].id = i;
// 	}
// 	scenarios_arq << calls[s].size() << "\n";
// 	int this_g = 5;
// 	for(auto& call: calls[s]){
// 		scenarios_arq << call.disc_time << " ";
// 		scenarios_arq << call.region << " " << call.priority << " " <<
// this_g << " "; 		scenarios_arq <<
// call.time_on_scene/slot_duration << " " << call.location.first;
// scenarios_arq << " " << call.location.second << " ";
// scenarios_arq << call.time_at_hospital/slot_duration << " ";
// scenarios_arq << call.cleaning_time/slot_duration << " ";
// scenarios_arq << call.clean_needed
// << " " << call.hosp_needed <<
// "\n";
// 		// fmt::print("{} {} ({}, {}) {} {}\n", call.id, call.time,
// 		// 	call.location.first, call.location.second,
// call.hosp_needed,
// 		// 	call.clean_needed);
// 	}
// }
// scenarios_arq.close();
