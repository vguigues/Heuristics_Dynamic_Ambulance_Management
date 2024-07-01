#include "../include/travel.h"
#include "../include/instance.h"
#include "../include/solver.h"


Solver::Solver(GRBEnv& env, vector<Call>& calls, vector<FastEquipment>& equips, 
	vector<Ambulance>& ambulances, Instance& ins, Travel& travel): 
    env(env), calls(calls), equips(equips), ambulances(ambulances), 
	travel(travel), ins(ins), time(0), first_time(calls.front().time), 
	last_time(calls.back().time+7200), waiting_on_scene(calls.size(), GRB_INFINITY),
	waiting_on_scene_penalized(calls.size(),GRB_INFINITY),
	waiting_to_hospital(calls.size(), GRB_INFINITY),
	calls_end(calls.size(), GRB_INFINITY),
	which_ambulance(calls.size(), -1), obj(0){
   
    set_calls_nearest_bases();
}

void Solver::set_calls_nearest_bases(){
	for(int i = 0; i < calls.size(); ++i){
		auto& call = calls[i];
		double min_dist = GRB_INFINITY;
		int min_ind = -1;
		for(int b = 0; b < ins.nb_bases; ++b){
			double dist = GRB_INFINITY;
			auto& base = ins.bases[b];
			if(call.clean_needed){
				dist = travel.vehicle_distance(ins.cleaning_bases[call.cleaning], base);
			}else if(call.hosp_needed){
				dist = travel.vehicle_distance(ins.hospitals[call.hospital], base); 
			}else{
				dist = travel.vehicle_distance(call.location, base);
			}

			if(dist < min_dist){
				min_dist = dist;
				min_ind = b;
			}
		}
		nearest_base.push_back(min_ind);
	}
}


void Solver::set_next_event(int& event_call, int& index_call){
	vector<double> future_arrival_times;
	for(int i = 0; i < ambulances.size(); ++i){
		auto& amb = ambulances[i];
		if(amb.arrival_time_at_f_last_trip > time + g_params.EPS){
			future_arrival_times.push_back(amb.arrival_time_at_f_last_trip);
		}
	}

	int nb_calls = calls.size();
	if(future_arrival_times.size() > 0){
		double min_arrival_time = *min_element(future_arrival_times.begin(),
			future_arrival_times.end());

		if(index_call < nb_calls-1 && calls[index_call+1].time <= 
			min_arrival_time){
			event_call = 1;
			index_call += 1;
			time = calls[index_call].time;
		}else{
			event_call = 0;
			time = min_arrival_time;
		}
	}else if(index_call < nb_calls - 1){
		event_call = 1;
		index_call += 1;
		time = calls[index_call].time;
	}
}

bool Solver::can_answer(Ambulance& amb, Call& call){
	bool amb_allowed = false;
	for(int amb_id: call.ambulances){
		if(amb_id == amb.id){
			amb_allowed = true;
			break;
		}
	}
	return amb_allowed && amb.type <= call.priority;
}

void Solver::print_results(){
	double average_scene = 0;
	double average_hospital = 0;
	double average_pen = 0;
	cout << "Call\tWait On Scene(s)\tWait Penalized(s)\tAmb:\n";
	for(int i = 0; i < calls.size(); ++i){
		auto& call = calls[i];
		cout << call << "\t"<< waiting_on_scene[i] << "\t";
		cout << waiting_on_scene_penalized[i] << "\t" << which_ambulance[i] << "\n";

		average_scene += waiting_on_scene[i];
		average_pen += waiting_on_scene_penalized[i];
		average_hospital += waiting_to_hospital[i];
	}
	average_scene /= calls.size();
	average_hospital /= calls.size();
	average_pen /= calls.size();
	cout << "Average: Real " << average_scene << " Pen " << average_pen << "\n";
	cout << "Total: Real " << average_scene*calls.size() << "  Pen "; 
	cout << average_pen*calls.size() << "\n";
}


Solver::~Solver(){}

MinMaxSolver::MinMaxSolver(GRBEnv& env, vector<Call>& calls, vector<FastEquipment>& equips,
    vector<Ambulance>& ambulances, 
	Instance& ins, Travel& travel): Solver(env, calls, equips, ambulances,ins,travel){
	travel.set_forward(false);
	queue.reserve(calls.size());
}


void MinMaxSolver::run(){
    int event_call = 1;
	time = calls[0].time;
	int index_call = 0;
	int calls_attended = 0;
	travel.set_forward(true);

    vector<FastEquipment> motos;
    vector<FastEquipment> drones;

    for(auto& fe: equips){
        if(fe.type == FastType::MOTORCYCLE){
            motos.push_back(fe);
        }else if(fe.type == FastType::DRONE){
            drones.push_back(fe);
        }
    }

	while(calls_attended < calls.size()){
        if(event_call == 1){
			queue.push_back(index_call);
		}

        vector<double> min_times(queue.size(), GRB_INFINITY);
        vector<double> min_pen_times(queue.size(), GRB_INFINITY);
        vector<double> min_times_ambs(queue.size(), GRB_INFINITY);
        vector<double> min_times_motos(queue.size(), GRB_INFINITY);
        vector<double> min_times_drones(queue.size(), GRB_INFINITY);

		vector<vector<double>> travel_times(queue.size(), 
			vector<double>(ambulances.size(), GRB_INFINITY));
        vector<vector<double>> travel_times_motos(queue.size(), 
            vector<double>(motos.size(), GRB_INFINITY));
        vector<vector<double>> travel_times_drones(queue.size(), 
            vector<double>(drones.size(), GRB_INFINITY));

		vector<int> index_ambs(queue.size(), -1);
        vector<int> index_motos(queue.size(), -1);
        vector<int> index_drones(queue.size(), -1);


        for(int i = 0; i < queue.size(); ++i){
            auto& call = calls[queue[i]];
            if(call.priority == -1){
                int j = 0;
                for(auto& amb: ambulances){
                    travel_times[i][j] = travel.get_response_time(amb,call,time);
                    if(travel_times[i][j] < min_times[i]){
                        min_times_ambs[i] = travel_times[i][j];
                        index_ambs = j;
                    }else if(abs(travel_times[i][j] - min_times[i]) < g_params.EPS && 
                        index_ambs[i] != -1 && amb.type > ambulances[index_ambs[i]].type){
                        min_times_ambs[i] = travel_times[i][j];
                        index_ambs[i] = j;
                    }
                    ++j;
                }
                j = 0;
                for(auto& moto: motos){
                    travel_times_motos[i][j] = travel.get_response_time(moto, call, time);
                    if(travel_times_motos[i][j] < min_times[i]){
                        min_times_motos[i] = travel_times_motos[i][j];
                        index_motos[i] = j;
                    }
                    ++j
                }
                j = 0;
                for(auto& drone: drones){
                    travel_times_drones[i][j] = travel.get_response_time(drone, call, time);
                    if(travel_times_drones[i][j] < min_times[i]){
                        min_times_drones[i] = travel_times_drones[i][j];
                        index_drones[i] = j;
                    }
                    ++j
                }

                if(min_times_motos[i] < min_times_ambs[i] || min_times_drones < min_times_ambs){
                    min_times[i] = min_times_ambs[i];
                    if(min_times_motos[i] < min_times_drones[i]){
                        min_pen_times[i] = ins.survival_function(min_times_motos[i]) + 
                            ins.penalties[0]*(min_times_ambs[i] - min_times_motos[i]);
                        index_drones[i] = -1;
                    }else{
                        min_pen_times[i] = ins.survival_function(min_times_drones[i]) +
                            ins.penalties[0]*(min_times_ambs[i] - min_times_drones[i]);
                        index_motos[i] = -1;
                    }
                }else{
                    min_times[i] = min_times_ambs[i];
                    min_pen_times[i] = ins.survival_function(min_times_ambs[i]);
                    index_motos[i] = -1;
                    index_drones[i] = -1;
                }
            }else{
                for(auto& amb: ambulances){
					travel_times[i][j] = travel.get_response_time(amb,call,time);
                    if(travel_times[i][j] < min_times[i]){
                        min_times_ambs[i] = travel_times[i][j];
                        index_ambs = j;
                    }else if(abs(travel_times[i][j] - min_times[i]) < g_params.EPS && 
                        index_ambs[i] != -1 && amb.type > ambulances[index_ambs[i]].type){
                        min_times_ambs[i] = travel_times[i][j];
                        index_ambs[i] = j;
                    }
					if(travel_times[i][j] < min_times[i]){
						min_times[i] = travel_times[i][j];
						index_ambs[i] = j;
					}
                    ++j;
				}
				min_times[i] = time + min_times[i] - call.time;
				min_pen_times[i] = ins.penalties[call.priority]*min_times[i];
            }
        }

		std::vector<int> remaining_indexes(queue.size(),0);
		for(int i = 0; i < queue.size(); ++i){
			remaining_indexes[i] = i;
		}

		std::vector<int> queue_aux;
		int nb_treated = 0;
		int total_queue = queue.size();

		while(nb_treated < total_queue){
			auto t0 = std::chrono::high_resolution_clock::now();
			//get time and index of the worst call, and the best ambulance to such call
			auto max_it = std::max_element(min_pen_times.begin(), min_pen_times.end());
			double max_time = *max_it;
			int current_call = std::distance(min_times.begin(), max_it);
			int call_ind = queue[remaining_indexes[current_call]];
			auto& call = calls[call_ind];
			int index_amb = index_ambs[remaining_indexes[current_call]];
			int index_drone = index_drones[remaining_indexes[current_call]];
			int index_moto = index_motos[remaining_indexes[current_call]];
			auto& best_amb = ambulances[index_amb];


			//TODO: Ajustar despacho
			if(call.priority == -1){
				if(index_moto >= 0){
					auto& best_moto = motos[index_moto];
					if(best_moto.arrival_time_at_f_last_trip <= time){ // free

					}
				}else(index_drone >= 0){
					auto& best_drone = drones[index_drone];
					if(best_moto.arrival_time_at_f_last_trip <= time){ // free

					}
				}
				if(amb.arrival_time_at_f_last_trip <= time){

				}else{

				}
			}else{
				if(index_amb >= 0 && best_amb.arrival_time_at_f_last_trip <= time){
					double waiting_on_scene_i = max_time;
					double waiting_to_hospital_i = best_amb.answer_call(call, 
						travel, ins, time, max_time - (time - call.time), 
						nearest_base[call_ind]);
					waiting_on_scene[call_ind] = waiting_on_scene_i;
					waiting_on_scene_penalized[call_ind] = 
						waiting_on_scene_i*ins.penalties[call.priority];
					waiting_to_hospital[call_ind] = waiting_to_hospital_i;
					which_ambulance[call_ind] = index_amb;
					calls_end[call_ind] = call.end;
					calls_attended++;

					//update travel times
					for(int k = 0; k < total_queue-nb_treated; ++k){
						int ind = remaining_indexes[k];
						if(k != current_call && index_amb == index_ambs[ind]){
							double time_to_h = best_amb.arrival_time_at_f_last_trip -
								time;
							double time_from_h_to_c = travel.travel_time(best_amb.free_location,
								calls[queue[ind]].location);
							travel_times[ind][best_amb.id] = time_to_h + time_from_h_to_c;
							auto min_it = std::min_element(travel_times[ind].begin(),
								travel_times[ind].end());
							double min_val = *min_it;
							int min_ind = std::distance(travel_times[ind].begin(), min_it);
							min_times[k] = min_val + time - calls[queue[ind]].time;
							index_ambs[ind] = min_ind;
						}
					}
				}
			}
		}
    }
}