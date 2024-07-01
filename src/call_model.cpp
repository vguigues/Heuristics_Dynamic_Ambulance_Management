#include "../include/travel.h"
#include "../include/ambulance.h"
#include "../include/call.h"
#include "../include/instance.h"
#include "../include/solver.h"
#include "../include/call_model.h"


// CallModel::CallModel(GRBEnv& env, Simulator& sim): ins(sim.ins), 
// 	calls(calls), queue(queue), ambulances(ambulances), time(time), 
// 	travel(sim.osrm, ins), model(env){
// 	load_model();
// }

CallModel::CallModel(GRBEnv& env, Solver& solver): ins(solver.ins),
	calls(solver.calls), queue(solver.queue), ambulances(solver.ambulances), 
	time(solver.time), travel(solver.travel), model(env){
	load_model();
}

void CallModel::load_model(){
	n_calls = queue.size();
	opt_which_amb = std::vector<int>(n_calls, -1);
	arrival_times = std::vector<double>(n_calls, GRB_INFINITY);
	n_ambulances = ins.nb_ambulances;
	M = 80000;
	std::stringstream name;

	x_ik = new GRBVar*[n_calls];
	x_ijk = new GRBVar**[n_calls];
	z_ik = new GRBVar*[n_calls];
	
	for(int i = 0; i < n_calls; ++i){
		x_ik[i] = new GRBVar[n_ambulances];
		z_ik[i] = new GRBVar[n_ambulances];
		auto& call = calls[queue[i]];
		for(int k = 0; k < n_ambulances; ++k){
			auto& amb = ambulances[k];
			if(amb.type <= call.priority){
				name << "x_i" << i << "_k" << k;
				x_ik[i][k] = model.addVar(0,1,0,GRB_BINARY, name.str());
				name.str(""); 
				name << "z_i" << i << "_k" << k;
				z_ik[i][k] = model.addVar(0,1,0,GRB_BINARY, name.str());
				name.str("");
			}

		}
		x_ijk[i] = new GRBVar*[n_calls];
		for(int j = 0; j < n_calls; ++j){
			x_ijk[i][j] = new GRBVar[n_ambulances];
			if(i != j){
				auto& call_j = calls[queue[j]];
				for(int k = 0; k < n_ambulances; ++k){
					auto& amb = ambulances[k];
					if(amb.type <= call.priority && amb.type <= call_j.priority){
						name << "x_i" << i << "_j" << j << "_k" << k;
						x_ijk[i][j][k] = model.addVar(0,1,0,GRB_BINARY, name.str());
						name.str("");
					}
				}
			}
		}
	}

	y_ic = new GRBVar*[n_calls];
	y_ih = new GRBVar*[n_calls];
	y_ihc = new GRBVar**[n_calls];

	for(int i = 0; i < n_calls; ++i){
		y_ic[i] = new GRBVar[ins.nb_cleaning_bases];
		auto& call = calls[queue[i]];
		for(int cb = 0; cb < ins.nb_cleaning_bases; ++cb){
			name << "y_i" << i << "_c" << cb;
			if(g_params.extended_model){
				y_ic[i][cb] = model.addVar(0,1,0,GRB_BINARY, name.str());
			}else if(cb == call.cleaning && !call.hosp_needed && call.clean_needed){
				y_ic[i][cb] = model.addVar(1,1,0,GRB_BINARY, name.str());
			}else{
				y_ic[i][cb] = model.addVar(0,0,0,GRB_BINARY, name.str());
			}
			name.str("");
		}
		y_ih[i] = new GRBVar[ins.nb_hospitals];
		y_ihc[i] = new GRBVar*[ins.nb_hospitals];
		for(int h = 0; h < ins.nb_hospitals; ++h){
			name << "y_i" << i << "_h" << h;
			if(g_params.extended_model){
				y_ih[i][h] = model.addVar(0,1,0,GRB_BINARY, name.str());
			}else if(h == call.hospital && call.hosp_needed && !call.clean_needed){
				y_ih[i][h] = model.addVar(1,1,0,GRB_BINARY, name.str());
			}else{
				y_ih[i][h] = model.addVar(0,0,0,GRB_BINARY, name.str());
			}
			name.str("");
			y_ihc[i][h] = new GRBVar[ins.nb_cleaning_bases];
			for(int cb = 0; cb < ins.nb_cleaning_bases; ++cb){
				name << "y_i" << i << "_h" << h << "_c" << cb;
				if(g_params.extended_model){
					y_ihc[i][h][cb] = model.addVar(0,1,0,GRB_BINARY, name.str());
				}else if(h == call.hospital && cb == call.cleaning && 
					call.hosp_needed && call.clean_needed){
					y_ihc[i][h][cb] = model.addVar(1,1,0,GRB_BINARY, name.str());
				}else{
					y_ihc[i][h][cb] = model.addVar(0,0,0,GRB_BINARY, name.str());
				}
				name.str("");
			}
		}
	}

	t_i = new GRBVar[n_calls];
	for(int i = 0; i < n_calls; ++i){
		name << "t_" << i;
		t_i[i] = model.addVar(0,M,0, GRB_CONTINUOUS,name.str());
		name.str("");
	}

	M_p = new GRBVar[ins.nb_priorities];
	for(int p = 0; p < ins.nb_priorities; ++p){
		name << "M_p" << p;
		M_p[p] = model.addVar(0,M,ins.penalties[p],GRB_CONTINUOUS, name.str());
		// M_p[p] = model.addVar(0,M,1,GRB_CONTINUOUS, name.str());
		name.str("");
	}
	model.update();
	add_constraints();
}

void CallModel::add_constraints(){
	std::stringstream name;
	for(int k = 0; k < n_ambulances; ++k){
		auto& amb = ambulances[k];
		GRBLinExpr exp;
		for (int i = 0; i < n_calls; ++i){
			auto& call = calls[queue[i]];
			if(amb.type <= call.priority){
				exp += x_ik[i][k];
				
			}
		}
		name << "con_begin_ride_k" << k;
		model.addConstr(exp, GRB_LESS_EQUAL, 1, name.str());
		name.str("");
	}
	// fmt::print("Cons begin ride\n");
	// std::cin.get();

	for(int i = 0; i < n_calls; ++i){
		GRBLinExpr exp;
		auto& call = calls[queue[i]];
		for(int k = 0; k < n_ambulances; ++k){
			auto& amb = ambulances[k];
			if(amb.type <= call.priority){
				// fmt::print("i = {}, k = {}\n",i,k);
				exp += x_ik[i][k];
			}
		}

		for(int j = 0; j < n_calls; ++j){
			if(i != j){
				auto& call_j = calls[queue[j]];
				for(int k = 0; k < n_ambulances; ++k){
					auto& amb = ambulances[k];
					if(amb.type <= call.priority && amb.type <= call_j.priority){
						// fmt::print("i = {}, j = {}, k = {}\n",i, j,k);
						exp += x_ijk[j][i][k];
					}
				}
			}
		}
		name << "con_demand_" << i;
		try{
			model.addConstr(exp, GRB_EQUAL, 1, name.str());
		}catch(GRBException& ex){
			fmt::print("{} {}\n", ex.getMessage(), name.str());
			std::cin.get();
		}
		name.str("");
	}

	// fmt::print("Cons demand\n");


	for(int k = 0; k < n_ambulances; ++k){
		auto& amb = ambulances[k];
		for(int i = 0; i < n_calls; ++i){
			auto& call = calls[queue[i]];
			if(amb.type <= call.priority){
				GRBLinExpr exp;
				exp += x_ik[i][k];

				for(int j = 0; j < n_calls; ++j){
					auto& call_j = calls[queue[j]];
					if(j != i && amb.type <= call_j.priority){
						exp += x_ijk[j][i][k];
					}
				}

				exp -= z_ik[i][k];
				for(int j = 0; j < n_calls; ++j){
					auto& call_j = calls[queue[j]];
					if(j != i && amb.type <= call_j.priority){
						exp -= x_ijk[i][j][k];
					}
				}
				name << "flow_k" << k << "_i" << i;
				model.addConstr(exp, GRB_EQUAL, 0, name.str());
				name.str("");
			}
		}
	}

	// fmt::print("Cons flow\n");

	for(int i = 0; i < n_calls; ++i){
		auto& call = calls[queue[i]];
		if(!call.hosp_needed && call.clean_needed){
			GRBLinExpr exp;
			for(int c = 0; c < ins.nb_cleaning_bases; ++c){
				exp += y_ic[i][c];
			}
			name << "leave_cb_i" << i;
			model.addConstr(exp, GRB_EQUAL, 1, name.str());
			name.str("");
		}
	}

	

	for(int i = 0; i < n_calls; ++i){
		auto& call = calls[queue[i]];
		if(call.hosp_needed && !call.clean_needed){
			GRBLinExpr exp;
			for(int h = 0; h < ins.nb_hospitals; ++h){
				exp += y_ih[i][h];
			}
			name << "leave_h_i" << i;
			model.addConstr(exp, GRB_EQUAL, 1, name.str());
			name.str("");
		}
	}




	for(int i = 0; i < n_calls; ++i){
		auto& call = calls[queue[i]];
		if(call.hosp_needed && call.clean_needed){
			GRBLinExpr exp;
			for(int h = 0; h < ins.nb_hospitals; ++h){
				for(int c = 0; c < ins.nb_cleaning_bases; ++c){
					exp += y_ihc[i][h][c];
				}
			}
			name << "leave_h_cb_i" << i;
			model.addConstr(exp,GRB_EQUAL,1,name.str());
			name.str("");
		}
	}

	// fmt::print("Cons leave\n");
	// int c_h = calls.size() / ins.nb_hospitals + 1;
	// int n_h = 0;

	// for(int h = 0; h < ins.nb_hospitals; ++h){
	// 	GRBLinExpr exp;
	// 	for(int i = 0; i < n_calls; ++i){
	// 		auto& call = calls[queue[i]];
	// 		if(call.hosp_needed && !call.clean_needed){
	// 			exp += y_ih[i][h];
	// 		}
	// 	}

	// 	for(int i = 0; i < n_calls; ++i){
	// 		auto& call = calls[queue[i]];
	// 		for(int c = 0; c < ins.nb_cleaning_bases; ++c){
	// 			if(call.hosp_needed && call.clean_needed){
	// 				exp += y_ihc[i][h][c];
	// 			}
	// 		}
	// 	}
	// 	name << "hosp_cap_h" << h;
	// 	model.addConstr(exp, GRB_LESS_EQUAL, c_h-n_h, name.str());
	// 	name.str("");
	// }
	
	// fmt::print("Cons cap\n");
	for(int i = 0; i < n_calls; ++i){
		auto& call = calls[queue[i]];
		for(int k = 0; k < n_ambulances; ++k){
			auto& amb = ambulances[k];
			if(amb.type <= call.priority){
				double lhs;
				if(amb.arrival_time_at_b_last_trip <= time){
					lhs = time + travel.travel_time(amb.base_location, call.location, amb);
					// fmt::print("amb {} at base: {} {}\n",k,time, 
					// 	travel.travel_time(amb.base_location, call.location));
				}else if(amb.arrival_time_at_f_last_trip <= time){
					Location current_location = travel.ambulance_position(amb, time);
					lhs = time + travel.travel_time(current_location, call.location, amb);
					// fmt::print("amb {} returning: {} {}\n",k,time,
					// 	travel.travel_time(current_location, call.location));
				}else{
					Location free_location = amb.free_location;
					lhs = amb.arrival_time_at_f_last_trip + 
						travel.travel_time(free_location, call.location, amb);
					// fmt::print("amb {} busy: {} {}\n",k,time,
					// 	travel.travel_time(free_location, call.location));
				}
				GRBLinExpr rhs =  t_i[i] + M*(1-x_ik[i][k]);
				name << "time_first_i" << i << "_k" << k;
				model.addConstr(lhs, GRB_LESS_EQUAL, rhs, name.str());
				name.str("");
			}
		}
	}

	// fmt::print("Cons time first\n");

	// double speed = ambulances[0].speed; //assuming all speeds equal
	for(int i = 0; i < n_calls; ++i){
		auto& call_i = calls[queue[i]];
		if(!call_i.hosp_needed && !call_i.clean_needed){
			for(int j = 0; j < n_calls; ++j){
				if(i != j){
					auto& call_j = calls[queue[j]];
					double time_i_j = travel.travel_time(call_i.location, call_j.location,
						ambulances[0]);
					GRBLinExpr lhs = t_i[i] + call_i.time_on_scene + time_i_j;
					GRBLinExpr rhs = 1;
					for(int k = 0; k < n_ambulances; ++k){
						auto& amb = ambulances[k];
						int min_priority = min(call_i.priority, call_j.priority);
						if(amb.type <= min_priority){
							rhs -= x_ijk[i][j][k];
						}
					}

					name << "time_others_i" << i << "_j" << j;
					model.addConstr(lhs, GRB_LESS_EQUAL, t_i[j] + M*rhs, name.str());
					name.str("");
				}
			}
		}
	}

	// fmt::print("Cons time no hosp no cb \n");

	for(int i = 0; i < n_calls; ++i){
		auto& call_i = calls[queue[i]];
		if(!call_i.hosp_needed && call_i.clean_needed){
			for(int c = 0; c < ins.nb_cleaning_bases; ++c){
				for(int j = 0; j < n_calls; ++j){
					if(i != j){
						auto& call_j = calls[queue[j]];
						double time_ci_cb = travel.travel_time(call_i.location, 
							ins.cleaning_bases[c], ambulances[0]);
						double time_cb_cj = travel.travel_time(ins.cleaning_bases[c],
							call_j.location, ambulances[0]);
						GRBLinExpr lhs = t_i[i] + call_i.time_on_scene +
							y_ic[i][c]*(time_ci_cb + call_i.cleaning_time + time_cb_cj);
						GRBLinExpr rhs = 1;
						for(int k = 0; k < n_ambulances; ++k){
							auto& amb = ambulances[k];
							int min_priority = min(call_i.priority, call_j.priority);
							if(amb.type <= min_priority){
								rhs -= x_ijk[i][j][k];
							}
						}
						name << "time_others_i" << i << "_j" << j << "_c" << c;
						try{
							model.addConstr(lhs, GRB_LESS_EQUAL, t_i[j] + M*rhs, 
								name.str());
						}catch(GRBException& ex){
							fmt::print("{} {}\n",ex.getMessage(), name.str());
						}
						name.str("");
					}
				}
			}
		}
	}

	// fmt::print("Cons time no hosp yes cb \n");

	for(int i = 0; i < n_calls; ++i){
		auto& call_i = calls[queue[i]];
		if(call_i.hosp_needed && !call_i.clean_needed){
			for(int h = 0; h < ins.nb_hospitals; ++h){
				for(int j = 0; j < n_calls; ++j){
					if(i != j){
						auto& call_j = calls[queue[j]];
						double time_ci_h = travel.travel_time(call_i.location, 
							ins.hospitals[h], ambulances[0]);
						double time_h_cj = travel.travel_time(ins.hospitals[h], 
							call_j.location, ambulances[0]);
						GRBLinExpr lhs = t_i[i] + call_i.time_on_scene +
							y_ih[i][h]*(time_ci_h + call_i.time_at_hospital+ time_h_cj);
						GRBLinExpr rhs = 1;
						for(int k = 0; k < n_ambulances; ++k){
							auto& amb = ambulances[k];
							int min_priority = min(call_i.priority, call_j.priority);
							if(amb.type <= min_priority){
								rhs -= x_ijk[i][j][k];
							}
						}
						name << "time_others_i" << i << "_j" << j << "_h" << h;
						model.addConstr(lhs, GRB_LESS_EQUAL, t_i[j] + M*rhs, name.str());
						name.str("");
					}
				}
			}
		}
	}

	// fmt::print("Cons time yes hosp no cb \n");

	for(int i = 0; i < n_calls; ++i){
		auto& call_i = calls[queue[i]];
		if(call_i.hosp_needed && call_i.clean_needed){
			for(int h = 0; h < ins.nb_hospitals; ++h){
				for(int c = 0; c < ins.nb_cleaning_bases; ++c){
					for(int j = 0; j < n_calls; ++j){
						if(i != j){
							auto& call_j = calls[queue[j]];
							double time_ci_h = travel.travel_time(call_i.location,
								ins.hospitals[h], ambulances[0]);
							double time_h_cb = travel.travel_time(ins.hospitals[h],
								ins.cleaning_bases[c], ambulances[0]);
							double time_cb_cj = travel.travel_time(ins.cleaning_bases[c], 
								call_j.location, ambulances[0]);
							GRBLinExpr lhs = t_i[i] + call_i.time_on_scene + 
								y_ihc[i][h][c]*(time_ci_h + call_i.time_at_hospital +
								time_h_cb + call_i.cleaning_time + time_cb_cj);
							GRBLinExpr rhs = 1;
							for(int k = 0; k < n_ambulances; ++k){
								auto& amb = ambulances[k];
								int min_priority = min(call_i.priority, call_j.priority);
								if(amb.type <= min_priority){
									rhs -= x_ijk[i][j][k];
								}
							}
							name << "time_others_i" << i << "_j" << j << "_h" << h;
							name << "_c" << c;
							model.addConstr(lhs, GRB_LESS_EQUAL, t_i[j] + M*rhs, 
								name.str());
							name.str("");
						}
					}
				}
			}
		}
	}

	// fmt::print("Cons time yes hosp yes cb \n");
	// fmt::print("Cons time others\n");

	for(int p = 0; p < ins.nb_priorities; ++p){
		for(int i = 0; i < n_calls; ++i){
			auto& call = calls[queue[i]];
			if(call.priority == p){
				if(call.hosp_needed && call.clean_needed){
					for(int h = 0; h < ins.nb_hospitals; ++h){
						for(int c = 0; c < ins.nb_cleaning_bases; ++c){
							double time_c_h = travel.travel_time(call.location,
								ins.hospitals[h], ambulances[0]);
							GRBLinExpr rhs = t_i[i] + call.time_on_scene + 
								y_ihc[i][h][c]*time_c_h;
							name <<"M_con_p" << p << "_h" << h << "_c" << c << "_i" << i;
							model.addConstr(M_p[p], GRB_GREATER_EQUAL, rhs, name.str());
							name.str("");
						}
					}
				}else if(call.hosp_needed){
					for(int h = 0; h < ins.nb_hospitals; ++h){
						double time_c_h = travel.travel_time(call.location,
								ins.hospitals[h], ambulances[0]);
						GRBLinExpr rhs = t_i[i] + call.time_on_scene + 
								y_ih[i][h]*time_c_h;
						name << "M_con_p" << p << "_h" << h << "_i" << i;
						model.addConstr(M_p[p], GRB_GREATER_EQUAL, rhs, name.str());
						name.str("");
					}
				}else{
					name << "M_con_p" << p << "_i";
					//M3 >= t_i + time_on_scene?
					model.addConstr(M_p[p], GRB_GREATER_EQUAL, t_i[i], name.str());
					name.str("");
				}
			}
		}
	}

	// fmt::print("Cons M\n");


	model.update();
	// model.set(GRB_IntParam_Presolve, 2);
	model.set(GRB_IntParam_OutputFlag, 0);
	model.set(GRB_IntParam_NumericFocus, 3);
}


void CallModel::solve(){
	
	// model.write("opt_alloc.lp");
	// std::cin.get();
	model.optimize();

	// model.computeIIS();
	// model.write("debug.ilp");

	// GRBVar* vars = model.getVars();
	// for(int i = 0; i < model.get(GRB_IntAttr_NumVars); ++i){
	// 	auto name = vars[i].get(GRB_StringAttr_VarName);
	// 	auto val = vars[i].get(GRB_DoubleAttr_X);
	// 	if(val > g_params.EPS){
	// 		fmt::print("{} = {}\n", name,val);
	// 	}
	// }
	// delete[] vars;
	// std::cin.get();

	if(model.get(GRB_IntAttr_Status) == GRB_OPTIMAL){
		for(int k = 0; k < n_ambulances; ++k){
			auto& amb = ambulances[k];
			for(int i = 0; i < n_calls; ++i){
				auto& call_i = calls[queue[i]];
				if(amb.type <= call_i.priority){
					double val_xi = x_ik[i][k].get(GRB_DoubleAttr_X);
					double val_zi = z_ik[i][k].get(GRB_DoubleAttr_X);
					if(val_xi > 0.5 || val_zi > 0.5){
						opt_which_amb[i] = k;
					}else{
						for(int j = 0; j < n_calls; ++j){
							if(i != j){
								auto& call_j = calls[queue[j]];
								if(amb.type <= call_j.priority){
									double val_xij = x_ijk[i][j][k].get(
										GRB_DoubleAttr_X);
									if(val_xij > 0.5){
										opt_which_amb[i] = k;
									}
								}
							}
						}
					}
				}
			}
		}
		for(int i = 0; i < n_calls; ++i){
			arrival_times[i] = t_i[i].get(GRB_DoubleAttr_X);
		}
		// std::cin.get();
	}else{
		std::cout << "Model is Infeasible!! Check debug.ilp file for IIS.\n";
		model.set(GRB_IntParam_OutputFlag,1);
		model.computeIIS();
		model.write("debug.ilp");
		std::cout << "finished debug.ilp\n";
		exit(1);
	}
}

void CallModel::print_results(){
	int status = model.get(GRB_IntAttr_Status);
	if(status == GRB_OPTIMAL){

		std::cout << "OBJ = " << model.get(GRB_DoubleAttr_ObjVal) << "\n";
		auto vars = model.getVars();
		for(int i = 0; i < model.get(GRB_IntAttr_NumVars); ++i){
			double val = vars[i].get(GRB_DoubleAttr_X);
			if(val > 0.5){
				std::cout << vars[i].get(GRB_StringAttr_VarName) << " = ";
				std::cout << val << "\n";
			}
		}

		delete[] vars;
	}else{
		std::cout << "Model is Infeasible!! Check debug.ilp file for IIS.\n";
		model.computeIIS();
		model.write("debug.ilp");
	}
}


// int CallModel::ambulances_at_hospital(){
// 	int count = 0;
// 	for(int k = 0; k < n_ambulances; ++k){
// 		auto& amb = ambulances[k];
// 		if(amb.is_at_hospital(time, travel)){
// 			++count;
// 		}
// 	}
// 	return count;
// }


CallModel::~CallModel(){
	for(int i = 0; i < n_calls; ++i){
		delete[] x_ik[i]; 
		delete[] z_ik[i];
		delete[] y_ic[i];
		delete[] y_ih[i];
		for(int j = 0; j < n_calls; ++j){
			delete[] x_ijk[i][j];
		}
		delete[] x_ijk[i];
		for(int h = 0; h < ins.nb_hospitals; ++h){
			delete[] y_ihc[i][h];
		}
		delete[] y_ihc[i];
	}

	delete[] x_ik; delete[] x_ijk; delete[] z_ik; delete[] y_ic;delete[] y_ih; 
	delete [] y_ihc; delete[] t_i; delete[] M_p;
}