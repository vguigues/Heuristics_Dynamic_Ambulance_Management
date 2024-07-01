#include "../include/future.h"

FutureCall::FutureCall(Data & data, GRBEnv & env, Solver& solver, int a0, int l0, int b0,
	int h0, std::vector<Call>& calls): data(data), model(env), solver(solver), 
	type_call0(data.type_call0), local_call0(data.local_call0), a0(a0), l0(l0), b0(b0), 
	h0(h0), calls(calls){
    
    A0_ab = new int*[data.types_amb];
	A0_ah = new int*[data.types_amb];
	A0_alb = new int**[data.types_amb];
	A0_calh = new int***[data.types_call];
	A0_callh = new int****[data.types_call];
	C = new int*[data.types_call];


	lambda = new int**[data.num_times-1];
	for(int t = 0; t < data.num_times-1; ++t){
		lambda[t] = new int*[data.types_call];
		for(int c = 0; c < data.types_call; ++c){
			lambda[t][c] = new int[data.num_locals];
			for(int l = 0; l < data.num_locals; ++l){
				lambda[t][c][l] = 0;
			}
		}
	}

	for(size_t i = 1; i < calls.size(); ++i){
		auto& call = calls[i];
		int t = data.get_time_slot(call.time);
		if(t != -1){
			// std::cout << "num_times = " << num_times << ", ";
			// std::cout << "T = " << t << " (" << call.time << "), ";
			// std::cout << "P = " << call.priority << ", ";
			// std::cout << "R = " << call.region << "\n";
			lambda[t][call.priority][call.region] += 1;
		}
	}


	for(int a = 0; a < data.types_amb; ++a){
		A0_ab[a] = new int[data.num_bases];
		for(int b = 0; b < data.num_bases; ++b){
			A0_ab[a][b] = data.A0_ab[a][b];
			// if(A0_ab[a][b] > 0){
			// 	fmt::print("a{} b{}\n",a,b);
			// }
		}
		A0_ah[a] = new int[data.num_hosps];
		for(int h = 0; h < data.num_hosps; ++h){
			A0_ah[a][h] = data.A0_ah[a][h];
			if(A0_ah[a][h] > 0){
				fmt::print("=====> a{} h{}\n",a,h);
			}
		}
		A0_alb[a] = new int*[data.total_locals];
		for(int l1 = 0; l1 < data.total_locals; ++l1){
			A0_alb[a][l1] = new int[data.num_bases];
			for(int b = 0; b < data.num_bases; ++b){
				A0_alb[a][l1][b] = data.A0_alb[a][l1][b];
				// if(A0_alb[a][l1][b] > 0){
				// 	fmt::print("a{} l{} b{}\n",a,l1,b);
				// }
			}
		}
	}

	for(int c = 0; c < data.types_call; ++c){
		A0_calh[c] = new int**[data.types_amb];
		A0_callh[c] = new int***[data.types_amb];
		for(int a = 0; a < data.types_amb; ++a){
			A0_calh[c][a] = new int*[data.total_locals];
			for(int l1 = 0; l1 <  data.total_locals; ++l1){
				A0_calh[c][a][l1] = new int[data.num_hosps];
				for(int h = 0; h < data.num_hosps; ++h){
					A0_calh[c][a][l1][h] = data.A0_calh[c][a][l1][h];
					// if(A0_calh[c][a][l1][h] > 0){
					// 	fmt::print("c{} a{} l{} h{}\n",c,a,l1,h);
					// }
				}
			}
			A0_callh[c][a] = new int**[data.total_locals];
			for(int l1 = 0; l1 < data.total_locals; ++l1){
				A0_callh[c][a][l1] = new int*[data.num_locals];
				for(int l = 0; l < data.num_locals; ++l){
					A0_callh[c][a][l1][l] = new int[data.num_hosps];
					for(int h = 0; h < data.num_hosps; ++h){
						A0_callh[c][a][l1][l][h] = data.A0_callh[c][a][l1][l][h];
						// if(A0_callh[c][a][l1][l][h] > 0){
						// 	fmt::print("c{} a{} l'{} l{} h{}\n",c,a,l1,l,h);
						// }
					}
				}
			}
		}
		C[c] = new int[data.num_locals];
		for(int l = 0; l < data.num_locals; ++l){
			C[c][l] = data.C[c][l];
		}
	}

    std::stringstream name;

	xt_cablh = new GRBVar*****[data.num_times-1];
	for(int t = 0; t < data.num_times-1; ++t){
		xt_cablh[t] = new GRBVar****[data.types_call];
		for(int c = 0; c < data.types_call; ++c){
			xt_cablh[t][c] = new GRBVar***[data.types_amb];
			for(int a = 0; a < data.types_amb; ++a){
				xt_cablh[t][c][a] = new GRBVar**[data.num_bases];
				for(int b = 0; b < data.num_bases; ++b){
					xt_cablh[t][c][a][b] = new GRBVar*[data.num_locals];
					for(int l = 0; l < data.num_locals; ++l){
						xt_cablh[t][c][a][b][l] = new GRBVar[data.num_hosps];
						for(int h = 0; h < data.num_hosps; ++h){
							if(data.A[c].find(a) != data.A[c].end() &&
								data.H[c][l].find(h) != data.H[c][l].end()){
								name << "xt"<< t+1 <<"_c"<<c<<"_a"<<a<<"_b"<<b<<"_l"<<l<<"_h";
								name << h;
								xt_cablh[t][c][a][b][l][h] = model.addVar(0, 
									data.num_ambs, 0, (data.relax["xt_cablh"] ? 
									GRB_CONTINUOUS : GRB_INTEGER),
									name.str());
								name.str("");
							}
						}
					}
				}
			}
		}
	}


	xt_cahlh = new GRBVar*****[data.num_times-1];
	for(int t = 0; t < data.num_times-1; ++t){
		xt_cahlh[t] = new GRBVar****[data.types_call];
		for(int c = 0; c < data.types_call; ++c){
			xt_cahlh[t][c] = new GRBVar***[data.types_amb];
			for(int a = 0; a < data.types_amb; ++a){
				xt_cahlh[t][c][a] = new GRBVar**[data.num_hosps];
				for(int h1 = 0; h1 < data.num_hosps; ++h1){
					xt_cahlh[t][c][a][h1] = new GRBVar*[data.num_locals];
					for(int l = 0; l < data.num_locals; ++l){
						xt_cahlh[t][c][a][h1][l] = new GRBVar[data.num_hosps];
						for(int h = 0; h < data.num_hosps; ++h){
							if(data.A[c].find(a) != data.A[c].end() &&
								data.H[c][l].find(h) != data.H[c][l].end()){

								name << "xt"<< t+1 << "_c"<< c << "_a"<<a<<"_h'"<<h1;
								name << "_l"<< l <<"_h" <<h;
								xt_cahlh[t][c][a][h1][l][h] = model.addVar(0,
									data.num_ambs, 0, (data.relax["xt_cahlh"] ? 
									GRB_CONTINUOUS : GRB_INTEGER),
									name.str());
								name.str("");
							}
							
						}
					}
				}
			}
		}			
	}


	xt_calblh = new GRBVar******[data.num_times-1];
	for(int t = 0; t < data.num_times-1; ++t){
		xt_calblh[t] = new GRBVar*****[data.types_call];
		for(int c = 0; c < data.types_call; ++c){
			xt_calblh[t][c] = new GRBVar****[data.types_amb];
			for(int a = 0; a < data.types_amb; ++a){
				xt_calblh[t][c][a] = new GRBVar***[data.total_locals];
				for(int l1 = 0; l1 < data.total_locals; ++l1){
					xt_calblh[t][c][a][l1] = new GRBVar**[data.num_bases];
					for(int b = 0; b < data.num_bases; ++b){
						xt_calblh[t][c][a][l1][b] = new GRBVar*[data.num_locals];
						for(int l = 0; l < data.num_locals; ++l){
							xt_calblh[t][c][a][l1][b][l] = new GRBVar[data.num_hosps];
							for(int h = 0; h < data.num_hosps; ++h){
								if(data.A[c].find(a) != data.A[c].end() &&
									data.H[c][l].find(h) != data.H[c][l].end()
									&& data.L_tab[t+1][a][b].find(l1) !=
									data.L_tab[t+1][a][b].end()){
									name << "xt" << t+1 << "_c" << c << "_a" << a << "_l'";
									name << l1 << "_b" << b << "_l" << l << "_h" << h;
									xt_calblh[t][c][a][l1][b][l][h]=model.addVar(
										0, data.num_ambs, 0, (data.relax["xt_calblh"] 
										? GRB_CONTINUOUS: GRB_INTEGER),name.str());
									name.str("");
								}
							}
						}
					}
				}
			}
		}
	}

	yt_ahb = new GRBVar***[data.num_times-1];
	for(int t = 0; t < data.num_times-1; ++t){
		yt_ahb[t] = new GRBVar**[data.types_amb];
		for(int a = 0; a < data.types_amb; ++a){
			yt_ahb[t][a] = new GRBVar*[data.num_hosps];
			for(int h = 0; h < data.num_hosps; ++h){
				yt_ahb[t][a][h] = new GRBVar[data.num_bases];
				for(int b = 0; b < data.num_bases; ++b){
					name << "yt" << t+1 << "_a" << a << "_h" << h << "_b" << b; 
					yt_ahb[t][a][h][b] = model.addVar(0,data.num_ambs,0,
						(data.relax["yt_ahb"] ? GRB_CONTINUOUS : GRB_INTEGER),
						name.str());
					name.str("");
				}
			}
		}
	}

	Ct_cl = new GRBVar**[data.num_times];
	for(int t = 0; t < data.num_times; ++t){
		Ct_cl[t] = new GRBVar*[data.types_call];
		for(int c = 0; c < data.types_call; ++c){
			Ct_cl[t][c] = new GRBVar[data.num_locals];
			for(int l = 0; l < data.num_locals; ++l){
				name << "Ct" << t+1 << "_c" << c << "_l" << l; 
				Ct_cl[t][c][l] = model.addVar(0, GRB_INFINITY, 0,
					(data.relax["ct_cl"] ? GRB_CONTINUOUS : GRB_INTEGER), 
					name.str());
				name.str("");
			}
		}
	}

	At_ab = new GRBVar**[data.num_times];
	for(int t = 0; t < data.num_times; ++t){
		At_ab[t] = new GRBVar*[data.types_amb];
		for(int a = 0; a < data.types_amb; ++a){
			At_ab[t][a] = new GRBVar[data.num_bases];
			for(int b = 0; b < data.num_bases; ++b){
				name << "At" << t+1 << "_a" << a << "_b" << b;
				At_ab[t][a][b] = model.addVar(0,data.num_ambs, 0,
					(data.relax["at_ab"] ? GRB_CONTINUOUS : GRB_INTEGER), 
					name.str());
				name.str("");
			}
		}
	}

	At_alb = new GRBVar***[data.num_times];
	for(int t = 0; t < data.num_times; ++t){
		At_alb[t] = new GRBVar**[data.types_amb];
		for(int a = 0; a < data.types_amb; ++a){
			At_alb[t][a] = new GRBVar*[data.total_locals];
			for(int l = 0; l < data.total_locals; ++l){
				At_alb[t][a][l] = new GRBVar[data.num_bases];
				for(int b = 0; b < data.num_bases; ++b){
					if(data.L_tab[t][a][b].find(l) !=
						data.L_tab[t][a][b].end()){
						name << "At" << t+1 << "_a" << a << "_l" << l;
						name << "_b" << b;
						At_alb[t][a][l][b] = model.addVar(0,data.num_ambs, 0, 
							(data.relax["at_alb"] ? GRB_CONTINUOUS
							: GRB_INTEGER), name.str());
						name.str("");	
					}
				}
			}
		}
	}

	//Build the objective function. GRBLinExpr is a expression
	//built by the GRBVars and constants.
	GRBLinExpr fo;
	double factor = 1;
	if(type_call0 != -1 && l0 != -1){
		fo += f(0,type_call0,a0, l0, data.bases[b0], local_call0, h0, data, factor);
	}else if(type_call0 != -1 && b0 != -1){
		fo += f(0,type_call0,a0,data.bases[b0], local_call0, h0, data, factor);
	}

	for(int t = 1; t < data.num_times; ++t){
		for(int c = 0; c < data.types_call; ++c){
			for(int a = 0; a < data.types_amb; ++a){
				if(data.A[c].find(a) != data.A[c].end()){
					for(int l = 0; l < data.num_locals; ++l){
						for(int h = 0; h < data.num_hosps; ++h){
							if(data.H[c][l].find(h) != data.H[c][l].end()){
								for(int b = 0; b < data.num_bases; ++b){
									fo += f(t,c,a,data.bases[b],l,h, 
										data,factor)*
									xt_cablh[t-1][c][a][b][l][h];
								}

								for(int h1 = 0; h1 < data.num_hosps; ++h1){
									fo += f(t,c,a,data.hospitals[h1],l,h,data,factor)*
									xt_cahlh[t-1][c][a][h1][l][h];
								}

								for(int l1 = 0; l1 < data.total_locals; ++l1){
									for(int b = 0; b < data.num_bases; ++b){
										if(data.L_tab[t][a][b].find(l1) !=
											data.L_tab[t][a][b].end()){
											fo += f(t,c,a,l1,data.bases[b],
											l, h, data, factor)*
											xt_calblh[t-1][c][a][l1][b][l][h];
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	factor = data.max_cost;
	for(int t = 0; t < data.num_times; ++t){
		for(int c = 0; c < data.types_call; ++c){
			for(int l = 0; l < data.num_locals; ++l){
				fo += g_tcl(t+1,c,l,factor)*Ct_cl[t][c][l];
			}
		}
	}

	//Adds the objective function to the model
	model.setObjective(fo, GRB_MINIMIZE);


    add_bases_constraints_t0();
	add_hospitals_constraints_t0();
	add_locations_constraints_t0();
	add_queues_constraints_t0();

	add_bases_constraints();
	add_hospitals_constraints();
	add_locations_constraints();
	add_queues_constraints();
	add_ambs_location_constraints();
	add_base_cap_constraints();
}


void FutureCall::solve(){
    model.update();

    model.set(GRB_IntParam_OutputFlag,0);
	model.set(GRB_IntParam_Presolve, 0);
	model.set(GRB_IntParam_Threads, 1);

	model.optimize();
	int status = model.get(GRB_IntAttr_Status);
	obj = 0;
	if(status == GRB_OPTIMAL){
		obj = model.get(GRB_DoubleAttr_ObjVal);
		// GRBVar* vars = model.getVars();
		// for(int i = 0; i < model.get(GRB_IntAttr_NumVars); ++i){
		// 	auto name = vars[i].get(GRB_StringAttr_VarName);
		// 	auto val = vars[i].get(GRB_DoubleAttr_X);
		// 	if(val > g_params.EPS){
		// 		fmt::print("{} = {}\n",name,val);
		// 	}
		// }
		// delete[] vars;
		// std::cin.get();
	}else{
		std::cout << "ERROR: CALL MODEL INFEASIBLE!\n";

		fmt::print("{} {} {} {}\n", a0, l0, b0, h0);
		model.computeIIS();
		model.write("future.ilp");
		exit(1);
	}

}


void FutureCall::add_bases_constraints_t0(){
	std::stringstream name;
	GRBLinExpr con = 0;

	for(int a = 0; a < data.types_amb; ++a){
		for(int b = 0; b < data.num_bases; ++b){
			con += A0_ab[a][b];

			for(auto l1: data.L_tab[0][a][b]){
				if(data.L(0,a,l1,data.bases[b]) == data.bases[b]){
					con += A0_alb[a][l1][b];
				}
			}

			if(a == a0 && ((l0 == -1 && b0 == b) || (l0 != -1 && b0 == b &&
				data.L(0,a0,l0,data.bases[b0]) == data.bases[b0]))){
				con -= 1;
			}

			name << "flow_bases_a" << a << "_b" << b;
			model.addConstr(con, GRB_EQUAL, At_ab[0][a][b], name.str());
			name.str("");
			con = 0;
		}
	}

}

void FutureCall::add_hospitals_constraints_t0(){
	std::stringstream name;
	GRBLinExpr con = 0;

	// for(int a = 0; a < data.types_amb; ++a){
	// 	for(int h = 0; h < data.num_hosps; ++h){
	// 		name << "flow_hospitals_a" <<a << "_h" << h;
	// 		model.addConstr(con, GRB_EQUAL, data.A0_ah[a][h], name.str());
	// 		name.str("");
	// 		con = 0;
	// 	}
	// }
}

void FutureCall::add_locations_constraints_t0(){
	std::stringstream name;
	GRBLinExpr con = 0;

	for(int a = 0; a < data.types_amb; ++a){
		for(int b = 0; b < data.num_bases; ++b){
			for(auto l1: data.L_tab[0][a][b]){
				for(auto l2: data.L_tab[0][a][b]){
					if(data.L(0,a,l2,data.bases[b]) == l1){
						con += data.A0_alb[a][l2][b];
					}
				}
				if(a == a0 && l0 != -1 && data.L(0,a,l0,data.bases[b0]) == l1){
					con -= 1;
				}

				name << "flow_locals_a" << a << "_b" << b << "_l" << l1;
				model.addConstr(con, GRB_EQUAL, At_alb[0][a][l1][b], name.str());
				name.str("");
				con = 0;
			}
		}
	}
}

void FutureCall::add_queues_constraints_t0(){
	std::stringstream name;
	GRBLinExpr con = 0;

	for(int c = 0; c < data.types_call; ++c){
		for(int l = 0; l < data.num_locals; ++l){
			con += C[c][l];
			if(c == type_call0 && l == local_call0 && a0 == -1 && b0 == -1 && h0 == -1){
				con += 1;
			}

			name << "flow_queue_c" << c << "_l" << l;
			
			model.addConstr(con, GRB_EQUAL, Ct_cl[0][c][l], name.str());
			name.str("");
			con = 0;
		}
	}
	
}

void FutureCall::add_bases_constraints(){
	GRBLinExpr con;
	std::stringstream name;

	for(int t = 0; t < data.num_times-1; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			for(int b = 0; b < data.num_bases; ++b){
				con += At_ab[t][a][b];
				for(int h = 0; h < data.num_hosps; ++h){
					if(data.L(t,a, data.hospitals[h],data.bases[b]) == data.bases[b]){
						con += yt_ahb[t][a][h][b];
					}
				}
				for(int l1 = 0; l1 < data.total_locals; ++l1){
					if(data.L_tab[t+1][a][b].find(l1) != 
						data.L_tab[t+1][a][b].end() &&
						data.L(t,a,l1,data.bases[b]) == data.bases[b]){
						con += At_alb[t][a][l1][b];
						for(int c = 0; c < data.types_call; ++c){
							for(int l = 0; l < data.num_locals; ++l){
								for(auto h: data.H[c][l]){
									if(data.A[c].find(a) != data.A[c].end()){
										//If using xt0, put t+1 on 1st index
										con -= xt_calblh[t][c][a][l1][b][l][h];	
									}
								}
							}
						}
					}
				}
				for(int c = 0; c < data.types_call; ++c){
					for(int l = 0; l < data.num_locals; ++l){
						for(auto h: data.H[c][l]){
							if(data.A[c].find(a) != data.A[c].end()){
								//If using xt0, put t+1 on 1st index
								con -= xt_cablh[t][c][a][b][l][h];
							}
						}
					}
				}

				name << "flow_bases_t" << t+1 << "_a" << a << "_b" << b;
				model.addConstr(con, GRB_EQUAL, At_ab[t+1][a][b], name.str());
				name.str("");
				con = 0;
			} 
		}
	}
}

void FutureCall::add_hospitals_constraints(){
	GRBLinExpr	con;
	std::stringstream name;
	int quantum = data.quantum;


	int arrival_time = INT_MAX;
	if(a0 != -1 && b0 != -1 && h0 != -1){
		if(l0 != -1){
			arrival_time = data.tao[0][type_call0][a0][l0][data.local_call0][h0];
		}else{
			arrival_time = data.tao[0][type_call0][a0][data.bases[b0]][data.local_call0][h0];
		}
	}

	for(int t = 0; t < data.num_times-1; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			for(int h = 0; h < data.num_hosps; ++h){
				for(int c = 0; c < data.types_call; ++c){
					for(int l = 0; l < data.num_locals; ++l){
						for(int h1 = 0; h1 < data.num_hosps; ++h1){
							if(data.H[c][l].find(h1) !=
								data.H[c][l].end() &&
								data.A[c].find(a) != data.A[c].end()){
								//If using xt0, put t+1 on 1st index
								con -= xt_cahlh[t][c][a][h][l][h1];
							}
						}
					}
				}

				for(int b = 0; b < data.num_bases; ++b){
					con -= yt_ahb[t][a][h][b];
				}

				for(int c = 0; c < data.types_call; ++c){
					for(int l1 = 0; l1 < data.total_locals; ++l1){
						for(int l = 0; l < data.num_locals; ++l){
							if((t)*quantum <= data.tao[0][c][a][l1][l][h] 
							&& data.tao[0][c][a][l1][l][h] 
							<= (t+1)*quantum){
								con += A0_callh[c][a][l1][l][h];
							}	
						}
					}
				}

				for(int c = 0; c < data.types_call; ++c){
					for(int l1 = 0; l1 < data.total_locals; ++l1){
						if((t)*quantum <= data.tao_0[c][a][l1][h] 
						&& data.tao_0[c][a][l1][h] <= (t+1)*quantum){
							con += A0_calh[c][a][l1][h];
						}
					}
				}

				for(int c = 0; c < data.types_call; ++c){
					for(int l = 0; l < data.num_locals; ++l){
						if(data.H[c][l].find(h) != data.H[c][l].end()){
							for(int b = 0; b < data.num_bases; ++b){
								for(int t1 = 1; t1 < t+1; ++t1){
									if((t)*quantum <= t1*quantum+
										data.tao[t1][c][a][data.bases[b]][l][h] && 
										t1*quantum+
										data.tao[t1][c][a][data.bases[b]][l][h] 
										<= (t+1)*quantum && data.A[c].find(a) != 
										data.A[c].end()){
										con += xt_cablh[t1-1][c][a][b][l][h];
									}
								}
							}
						}
					}
				}

				for(int c = 0; c < data.types_call; ++c){
					for(int l = 0; l < data.num_locals; ++l){
						if(data.H[c][l].find(h) != data.H[c][l].end()){
							for(int h1 = 0; h1 < data.num_hosps; ++h1){
								for(int t1 = 1; t1 < t+1; ++t1){
									if((t)*quantum <= t1*quantum+
										data.tao[t1][c][a][data.hospitals[h1]][l][h] &&
										t1*quantum+data.tao[t1][c][a][data.hospitals[h1]][l][h] <= 
										(t+1)*quantum &&
										data.A[c].find(a) != data.A[c].end()){
										con += xt_cahlh[t1-1][c][a][h1][l][h];
									}
								}	
							}
						}
					}
				}

				if(a0 == a && h0 == h && t*quantum <= arrival_time && 
					arrival_time <= (t+1)*quantum){
					con += 1;
				}

				for(int c = 0; c < data.types_call; ++c){
					for(int l = 0; l < data.num_locals; ++l){
						if(data.H[c][l].find(h) != data.H[c][l].end()){
							for(int l1 = 0; l1 < data.total_locals; ++l1){
								for(int b = 0; b < data.num_bases; ++b){
									for(int t1 = 1; t1 < t+1; ++t1){
										if((t)*quantum <= t1*quantum+
											data.tao[t1][c][a][l1][l][h] &&
											t1*quantum+data.tao[t1][c][a][l1][l][h] <= (t+1)*quantum &&
											data.A[c].find(a) != data.A[c].end()
											&& data.L_tab[t1-1][a][b].find(l1) != 
											data.L_tab[t1-1][a][b].end()){
											con += xt_calblh[t1-1][c][a][l1][b][l][h];
										}
									}
								}
							}
						}
					}
				}

				name << "flow_hospitals_t" << t+1 << "_a" << a << "_h" << h;
				model.addConstr(con, GRB_EQUAL, 0,  name.str());
				name.str("");
				con = 0;
			}
		}
	}
}

void FutureCall::add_locations_constraints(){
	GRBLinExpr con;
	std::stringstream name;

	for(int t = 0; t < data.num_times-1; ++t){
		for(int a  = 0; a < data.types_amb; ++a){
			for(int b = 0; b < data.num_bases; ++b){
				for(auto l1: data.L_tab[t+1][a][b]){
					for(int h1 = 0; h1 < data.num_hosps; ++h1){
						if(data.L(t,a, data.hospitals[h1],data.bases[b]) == l1){
							con += yt_ahb[t][a][h1][b];
						}
					}

					for(auto l2: data.L_tab[t+1][a][b]){
						if(data.L(t,a,l2, data.bases[b]) == l1){
							con += At_alb[t][a][l2][b];
							for(int c = 0; c < data.types_call; ++c){
								for(int l = 0; l < data.num_locals; ++l){
									for(int h = 0; h < data.num_hosps; ++h){
										if(data.H[c][l].find(h) !=
											data.H[c][l].end() &&
											data.A[c].find(a)!=data.A[c].end()){
											//If using xt0, put t+1 on 1st index
											con -= xt_calblh[t][c][a][l2][b][l][h];	
										}
									}
								}
							}
						}
					}
					name << "flow_locals_t" << t+1 << "_" << a << "_"<<b<<"_" << l1;
					model.addConstr(con,GRB_EQUAL, At_alb[t+1][a][l1][b], name.str());
					name.str("");
					con = 0;
				}
			}
		}
	}
}

void FutureCall::add_queues_constraints(){
	GRBLinExpr con;
	std::stringstream name;

	for(int t = 0; t < data.num_times-1; ++t){
		for(int c = 0; c < data.types_call; ++c){
			for(int l = 0; l < data.num_locals; ++l){
				con += Ct_cl[t][c][l];
				con += lambda[t][c][l];
				for(auto a: data.A[c]){
					for(int b = 0; b < data.num_bases; ++b){
						for(auto h: data.H[c][l]){
							//If using xt0, put t+1 on 1st index
							con -= xt_cablh[t][c][a][b][l][h];
						}
					}
				}

				for(auto a: data.A[c]){
					for(int h1 = 0; h1 < data.num_hosps; ++h1){
						for(auto h: data.H[c][l]){
							//If using xt0, put t+1 on 1st index
							con -= xt_cahlh[t][c][a][h1][l][h];
						}
					}
				}

				for(auto a: data.A[c]){
					for(int b = 0; b < data.num_bases; ++b){
						for(auto l1: data.L_tab[t+1][a][b]){
							for(auto h: data.H[c][l]){
								//If using xt0, put t+1 on 1st index
								con -= xt_calblh[t][c][a][l1][b][l][h];
							}
						}
					}
				}
				name << "flow_queue_t" << t+1 << "_c" << c << "_l" << l;
				model.addConstr(con,GRB_EQUAL,Ct_cl[t+1][c][l], name.str());
				name.str("");
				con = 0;
			}
		}
	}
}

void FutureCall::add_ambs_location_constraints(){
	GRBLinExpr con;
	std::stringstream name;

	//Location constraints
	for(int t = 0; t < data.num_times-1; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			for(int b = 0; b < data.num_bases; ++b){
				for(int l1 = 0; l1 < data.total_locals; ++l1){
					for(int c = 0; c < data.types_call; ++c){
						for(int l = 0; l < data.num_locals; ++l){
							for(int h = 0; h < data.num_hosps; ++h){
								if(data.H[c][l].find(h) != data.H[c][l].end() &&
									data.A[c].find(a) != data.A[c].end() &&
									data.L_tab[t+1][a][b].find(l1) != 
									data.L_tab[t+1][a][b].end()){
									//If using xt0, put t+1 on 1st index
									con += xt_calblh[t][c][a][l1][b][l][h];
								}
							}
						}
					}
					if(data.L_tab[t][a][b].find(l1) != data.L_tab[t][a][b].end()){
						name << "dispatch_locals_t" << t+1 << "_a" << a << "_b";
						name << b << "_l" << l1;
						model.addConstr(con, GRB_LESS_EQUAL, At_alb[t][a][l1][b],
							name.str());
					}
					name.str("");
					con = 0;
				}
			}
		}
	}


	// Max of ambulance at bases
	for(int t = 0; t < data.num_times-1; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			for(int b = 0; b < data.num_bases; ++b){

				for(int c = 0; c < data.types_call; ++c){
					for(int l = 0; l < data.num_locals; ++l){
						for(int h = 0; h < data.num_hosps; ++h){
							if(data.H[c][l].find(h) != data.H[c][l].end() &&
								data.A[c].find(a) != data.A[c].end()){
								//If using xt0, put t+1 on 1st index
								con += xt_cablh[t][c][a][b][l][h];
							} 
						}
					}
				}

				name << "dispatch_bases_t" << t+1 << "_a" << a << "_b" << b;
				model.addConstr(con, GRB_LESS_EQUAL, At_ab[t][a][b], name.str());
				name.str("");
				con = 0;
			}
		}
	}
}

void FutureCall::add_base_cap_constraints(){
	GRBLinExpr con;
	std::stringstream name;

	for(int t = 0; t < data.num_times; ++t){
		for(int b = 0; b < data.num_bases; ++b){
			for(int a = 0; a < data.types_amb; ++a){
				con += At_ab[t][a][b];
			}
			name << "cap_bases_" << t+1 << "_" << b;
			model.addConstr(con,GRB_LESS_EQUAL,data.cap_bases[b], name.str());
			name.str("");
			con = 0;
		}
	}
}

FutureCall::~FutureCall(){
	for(int t = 0; t < data.num_times-1; ++t){
		for(int c = 0; c < data.types_call; ++c){
			for(int a = 0; a < data.types_amb; ++a){
				for(int b = 0; b < data.num_bases; ++b){
					for(int l = 0; l < data.num_locals; ++l){
						delete[] xt_cablh[t][c][a][b][l];
					}
					delete[] xt_cablh[t][c][a][b];
				}
				delete[] xt_cablh[t][c][a];
			}
			delete[] xt_cablh[t][c];
		}
		delete[] xt_cablh[t];
	}
	//If using xt0, loop until num_times instead of num_times-1
	for(int t = 0; t < data.num_times-1; ++t){
		for(int c = 0; c < data.types_call; ++c){
			for(int a = 0; a < data.types_amb; ++a){
				for(int h1 = 0; h1 < data.num_hosps; ++h1){
					for(int l = 0; l < data.num_locals; ++l){
						delete[] xt_cahlh[t][c][a][h1][l];
					}
					delete[] xt_cahlh[t][c][a][h1];
				}
				delete[] xt_cahlh[t][c][a];
			}
			delete[] xt_cahlh[t][c];
		}
		delete[] xt_cahlh[t];
	}
	//If using xt0, loop until num_times instead of num_times-1
	for(int t = 0; t < data.num_times-1; ++t){
		for(int c = 0; c < data.types_call; ++c){
			for(int a = 0; a < data.types_amb; ++a){
				for(int l1 = 0; l1 < data.total_locals; ++l1){
					for(int b = 0; b < data.num_bases; ++b){
						for(int l = 0; l < data.num_locals; ++l){
							delete[] xt_calblh[t][c][a][l1][b][l];
						}
						delete[] xt_calblh[t][c][a][l1][b];
					}
					delete[] xt_calblh[t][c][a][l1];
				}
				delete[] xt_calblh[t][c][a];
			}
			delete[] xt_calblh[t][c];
		}
		delete[] xt_calblh[t];
	}

	for(int t = 0; t < data.num_times-1; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			for(int h = 0; h < data.num_hosps; ++h){
				delete[] yt_ahb[t][a][h];
			}
			delete[] yt_ahb[t][a];
		}
		delete[] yt_ahb[t];
	}

	for(int t = 0; t < data.num_times; ++t){
		for(int c = 0; c < data.types_call; ++c){
			delete[] Ct_cl[t][c];
		}
		delete[] Ct_cl[t];
	}

	for(int t = 0; t < data.num_times; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			delete[] At_ab[t][a];
		}
		delete[] At_ab[t];
	}

	for(int t = 0; t < data.num_times; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			for(int l = 0; l < data.total_locals; ++l){
				delete[] At_alb[t][a][l];
			}
			delete[] At_alb[t][a];
		}
		delete[] At_alb[t];
	}



	for(int a = 0; a < data.types_amb; ++a){
		for(int l1 = 0; l1 < data.total_locals; ++l1){
			delete[] A0_alb[a][l1];
		}
		delete[] A0_ab[a];
		delete[] A0_ah[a];
		delete[] A0_alb[a];
	}

	for(int c = 0; c < data.types_call; ++c){
		for(int a = 0; a < data.types_amb; ++a){
			for(int l1 = 0; l1 <  data.total_locals; ++l1){
				delete[]  A0_calh[c][a][l1];
			}
			delete[]  A0_calh[c][a];
			for(int l1 = 0; l1 < data.total_locals; ++l1){
				for(int l = 0; l < data.num_locals; ++l){
					delete[] A0_callh[c][a][l1][l];
				}
				delete[]  A0_callh[c][a][l1];
			}
			delete[]  A0_callh[c][a];
		}
		delete[]  A0_calh[c];
		delete[]  A0_callh[c];
		delete[] C[c];
	}
	delete[] A0_ab; delete[] A0_ah; delete[] A0_alb; delete[] A0_calh;
	delete[] A0_callh; delete[] C;

	delete[] xt_cablh; delete[] xt_cahlh; delete[] xt_calblh;
	delete[] Ct_cl; delete[] At_ab; delete[] At_alb; delete[] yt_ahb;


	for(int t = 0; t < data.num_times-1; ++t){
		for(int c = 0; c < data.types_call; ++c){
			delete[] lambda[t][c];
		}
		delete[] lambda[t];
	}
	delete[] lambda;
}

//===========================================================================================


FutureAmbulance::FutureAmbulance(Data & data, GRBEnv & env, Solver& solver, 
	Call* call, int h_dest,int b_ret, std::vector<Call>& calls): data(data), 
	call(call), model(env), solver(solver), h_dest(h_dest), b_ret(b_ret), 
	calls(calls){

	Travel travel(false);
	a0 = data.ambulance->type;
	for(int h = 0; h < data.num_hosps; ++h){
		if(travel.euclidian_distance(data.ambulance->free_location, solver.ins.hospitals[h]) < 
			g_params.EPS){
			h0 = h;
			break;
		}
	}
	if(h_dest != -1 && call != NULL){
		c0 = call->priority;
		l0 = data.get_location_index(call->location);
	}else{
		h_dest = c0 = l0 = -1;
	}

	A0_ab = new int*[data.types_amb];
	A0_ah = new int*[data.types_amb];
	A0_alb = new int**[data.types_amb];
	A0_calh = new int***[data.types_call];
	A0_callh = new int****[data.types_call];
	C = new int*[data.types_call];


	lambda = new int**[data.num_times-1];
	for(int t = 0; t < data.num_times-1; ++t){
		lambda[t] = new int*[data.types_call];
		for(int c = 0; c < data.types_call; ++c){
			lambda[t][c] = new int[data.num_locals];
			for(int l = 0; l < data.num_locals; ++l){
				lambda[t][c][l] = 0;
			}
		}
	}

	for(size_t i = 1; i < calls.size(); ++i){
		auto& call = calls[i];
		int t = data.get_time_slot(call.time);
		if(t != -1){
			// std::cout << "num_times = " << num_times << ", ";
			// std::cout << "T = " << t << " (" << call.time << "), ";
			// std::cout << "P = " << call.priority << ", ";
			// std::cout << "R = " << call.region << "\n";
			lambda[t][call.priority][call.region] += 1;
		}
	}


	for(int a = 0; a < data.types_amb; ++a){
		A0_ab[a] = new int[data.num_bases];
		for(int b = 0; b < data.num_bases; ++b){
			A0_ab[a][b] = data.A0_ab[a][b];
		}
		A0_ah[a] = new int[data.num_hosps];
		for(int h = 0; h < data.num_hosps; ++h){
			A0_ah[a][h] = data.A0_ah[a][h];
		}
		A0_alb[a] = new int*[data.total_locals];
		for(int l1 = 0; l1 < data.total_locals; ++l1){
			A0_alb[a][l1] = new int[data.num_bases];
			for(int b = 0; b < data.num_bases; ++b){
				A0_alb[a][l1][b] = data.A0_alb[a][l1][b];
			}
		}
	}

	for(int c = 0; c < data.types_call; ++c){
		A0_calh[c] = new int**[data.types_amb];
		A0_callh[c] = new int***[data.types_amb];
		for(int a = 0; a < data.types_amb; ++a){
			A0_calh[c][a] = new int*[data.total_locals];
			for(int l1 = 0; l1 <  data.total_locals; ++l1){
				A0_calh[c][a][l1] = new int[data.num_hosps];
				for(int h = 0; h < data.num_hosps; ++h){
					A0_calh[c][a][l1][h] = data.A0_calh[c][a][l1][h];
				}
			}
			A0_callh[c][a] = new int**[data.total_locals];
			for(int l1 = 0; l1 < data.total_locals; ++l1){
				A0_callh[c][a][l1] = new int*[data.num_locals];
				for(int l = 0; l < data.num_locals; ++l){
					A0_callh[c][a][l1][l] = new int[data.num_hosps];
					for(int h = 0; h < data.num_hosps; ++h){
						A0_callh[c][a][l1][l][h] = data.A0_callh[c][a][l1][l][h];
					}
				}
			}
		}
		C[c] = new int[data.num_locals];
		for(int l = 0; l < data.num_locals; ++l){
			C[c][l] = data.C[c][l];
		}
	}

	std::stringstream name;

	xt_cablh = new GRBVar*****[data.num_times-1];
	for(int t = 0; t < data.num_times-1; ++t){
		xt_cablh[t] = new GRBVar****[data.types_call];
		for(int c = 0; c < data.types_call; ++c){
			xt_cablh[t][c] = new GRBVar***[data.types_amb];
			for(int a = 0; a < data.types_amb; ++a){
				xt_cablh[t][c][a] = new GRBVar**[data.num_bases];
				for(int b = 0; b < data.num_bases; ++b){
					xt_cablh[t][c][a][b] = new GRBVar*[data.num_locals];
					for(int l = 0; l <data.num_locals; ++l){
						xt_cablh[t][c][a][b][l] = new GRBVar[data.num_hosps];
						for(int h = 0; h < data.num_hosps; ++h){
							if(data.A[c].find(a) != data.A[c].end() &&
								data.H[c][l].find(h) != data.H[c][l].end()){
								name << "xt"<< t+1 <<"_c"<<c<<"_a"<<a<<"_b"<<b;
								name <<"_l"<<l<<"_h" << h;
								xt_cablh[t][c][a][b][l][h] = model.addVar(0,
									data.num_ambs,
									0, (data.relax["xt_cablh"] ? GRB_CONTINUOUS : GRB_INTEGER),
									name.str());
								name.str("");
							}
						}
					}
				}
			}
		}
	}

	xt_cahlh = new GRBVar*****[data.num_times-1];
	for(int t = 0; t < data.num_times-1; ++t){
		xt_cahlh[t] = new GRBVar****[data.types_call];
		for(int c = 0; c < data.types_call; ++c){
			xt_cahlh[t][c] = new GRBVar***[data.types_amb];
			for(int a = 0; a < data.types_amb; ++a){
				xt_cahlh[t][c][a] = new GRBVar**[data.num_hosps];
				for(int h1 = 0; h1 < data.num_hosps; ++h1){
					xt_cahlh[t][c][a][h1] = new GRBVar*[data.num_locals];
					for(int l = 0; l < data.num_locals; ++l){
						xt_cahlh[t][c][a][h1][l] = new GRBVar[data.num_hosps];
						for(int h = 0; h < data.num_hosps; ++h){
							if(data.A[c].find(a) != data.A[c].end() &&
								data.H[c][l].find(h) != data.H[c][l].end()){
								name << "xt"<< t+1 << "_c"<< c << "_a"<<a<<"_h'"<<h1;
								name << "_l"<< l <<"_h" <<h;
								xt_cahlh[t][c][a][h1][l][h] = model.addVar(0, data.num_ambs,
									0, (data.relax["xt_cahlh"] ? GRB_CONTINUOUS : GRB_INTEGER),
									name.str());
								name.str("");
							}
							
						}
					}
				}
			}
		}			
	}

	xt_calblh = new GRBVar******[data.num_times-1];
	for(int t = 0; t < data.num_times-1; ++t){
		xt_calblh[t] = new GRBVar*****[data.types_call];
		for(int c = 0; c < data.types_call; ++c){
			xt_calblh[t][c] = new GRBVar****[data.types_amb];
			for(int a = 0; a < data.types_amb; ++a){
				xt_calblh[t][c][a] = new GRBVar***[data.total_locals];
				for(int l1 = 0; l1 < data.total_locals; ++l1){
					xt_calblh[t][c][a][l1] = new GRBVar**[data.num_bases];
					for(int b = 0; b < data.num_bases; ++b){
						xt_calblh[t][c][a][l1][b] = new GRBVar*[data.num_locals];
						for(int l = 0; l < data.num_locals; ++l){
							xt_calblh[t][c][a][l1][b][l] = new GRBVar[data.num_hosps];
							for(int h = 0; h < data.num_hosps; ++h){
								if(data.A[c].find(a) != data.A[c].end() &&
									data.H[c][l].find(h) != data.H[c][l].end()
									&& data.L_tab[t+1][a][b].find(l1) !=
									data.L_tab[t+1][a][b].end()){
									name << "xt" << t+1 << "_c" << c << "_a" << a << "_l'";
									name << l1 << "_b" << b << "_l" << l << "_h" << h;
									xt_calblh[t][c][a][l1][b][l][h] = model.addVar(
										0,data.num_ambs, 0,
										(data.relax["xt_calblh"] ? GRB_CONTINUOUS : GRB_INTEGER),
										name.str());
									name.str("");
								}
							}
						}
					}
				}
			}
		}
	}

	yt_ahb = new GRBVar***[data.num_times-1];
	for(int t = 0; t < data.num_times-1; ++t){
		yt_ahb[t] = new GRBVar**[data.types_amb];
		for(int a = 0; a < data.types_amb; ++a){
			yt_ahb[t][a] = new GRBVar*[data.num_hosps];
			for(int h = 0; h < data.num_hosps; ++h){
				yt_ahb[t][a][h] = new GRBVar[data.num_bases];
				for(int b = 0; b < data.num_bases; ++b){
					name << "yt" << t+1 << "_a" << a << "_h" << h << "_b" << b; 
					yt_ahb[t][a][h][b] = model.addVar(0,data.num_ambs,0,
						(data.relax["yt_ahb"] ? GRB_CONTINUOUS : GRB_INTEGER),
						name.str());
					name.str("");
				}
			}
		}
	}

	Ct_cl = new GRBVar**[data.num_times];
	for(int t = 0; t < data.num_times; ++t){
		Ct_cl[t] = new GRBVar*[data.types_call];
		for(int c = 0; c < data.types_call; ++c){
			Ct_cl[t][c] = new GRBVar[data.num_locals];
			for(int l = 0; l < data.num_locals; ++l){
				name << "Ct" << t+1 << "_c" << c << "_l" << l; 
				Ct_cl[t][c][l] = model.addVar(0, GRB_INFINITY, 0,
					(data.relax["ct_cl"] ? GRB_CONTINUOUS : GRB_INTEGER), name.str());
				name.str("");
			}
		}
	}

	At_ab = new GRBVar**[data.num_times];
	for(int t = 0; t < data.num_times; ++t){
		At_ab[t] = new GRBVar*[data.types_amb];
		for(int a = 0; a < data.types_amb; ++a){
			At_ab[t][a] = new GRBVar[data.num_bases];
			for(int b = 0; b < data.num_bases; ++b){
				name << "At" << t+1 << "_a" << a << "_b" << b;
				At_ab[t][a][b] = model.addVar(0,data.num_ambs, 0,
					(data.relax["at_ab"] ? GRB_CONTINUOUS : GRB_INTEGER), name.str());
				name.str("");
			}
		}
	}

	At_alb = new GRBVar***[data.num_times];
	for(int t = 0; t < data.num_times; ++t){
		At_alb[t] = new GRBVar**[data.types_amb];
		for(int a = 0; a < data.types_amb; ++a){
			At_alb[t][a] = new GRBVar*[data.total_locals];
			for(int l = 0; l < data.total_locals; ++l){
				At_alb[t][a][l] = new GRBVar[data.num_bases];
				for(int b = 0; b < data.num_bases; ++b){
					if( data.L_tab[t][a][b].find(l) ==
						data.L_tab[t][a][b].end())
						continue;
					name << "At" << t+1 << "_a" << a << "_l" << l << "_b" << b;
					At_alb[t][a][l][b] = model.addVar(0,data.num_ambs,0,
						(data.relax["at_alb"] ? GRB_CONTINUOUS : GRB_INTEGER), name.str());
					name.str("");
				}
			}
		}
	}

	GRBLinExpr fo;
	double factor = 1;

	if(c0 != -1){
		fo += f(0, c0, data.ambulance->type, h0, l0, h_dest, data, factor);
	}

	factor = 1;
	for(int t = 1; t < data.num_times; ++t){
		for(int c = 0; c < data.types_call; ++c){
			for(int a = 0; a < data.types_amb; ++a){
				if(data.A[c].find(a) != data.A[c].end()){
					for(int l = 0; l < data.num_locals; ++l){
						for(int h = 0; h < data.num_hosps; ++h){
							if(data.H[c][l].find(h) != data.H[c][l].end()){
								for(int b = 0; b < data.num_bases; ++b){
									fo += f(t-1,c,a,data.bases[b],l,h, 
										data,factor)*
									xt_cablh[t-1][c][a][b][l][h];
								}

								for(int h1 = 0; h1 < data.num_hosps; ++h1){
									fo += f(t-1,c,a,data.hospitals[h1],l,
										h,data,factor)*
									xt_cahlh[t-1][c][a][h1][l][h];
								}

								for(int l1 = 0; l1 < data.total_locals; ++l1){
									for(int b = 0; b < data.num_bases; ++b){
										if(data.L_tab[t][a][b].find(l1) !=
											data.L_tab[t][a][b].end()){
											fo += f(t-1,c,a,l1,data.bases[b],
											l, h, data, 
											factor)*
											xt_calblh[t-1][c][a][l1][b][l][h];
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	// factor = 125;
	for(int t = 0; t < data.num_times; ++t){
		factor = data.max_cost;
		for(int c = 0; c < data.types_call; ++c){
			for(int l = 0; l < data.num_locals; ++l){
				fo += g_tcl(t+1,c,l,factor)*Ct_cl[t][c][l];
			}
		}
		factor = 0;
		for(int a = 0; a < data.types_amb; ++a){
			for(int b = 0; b < data.num_bases; ++b){
				fo += g_tab(t+1,a,b,factor)*At_ab[t][a][b];
				for(int l1 = 0; l1 < data.total_locals; ++l1){
					if(data.L_tab[t][a][b].find(l1) != 
						data.L_tab[t][a][b].end()){
						fo+= g_talb(t+1,a,l1,data.bases[b],data,factor)*
						At_alb[t][a][l1][b];
					}
				}
			}
		}
	}

	model.setObjective(fo,GRB_MINIMIZE);

	model.update();

	add_bases_constraints_t0();
	add_hospitals_constraints_t0();
	add_locations_constraints_t0();
	add_queues_constraints_t0();

	add_bases_constraints();
	add_hospitals_constraints();
	add_locations_constraints();
	add_queues_constraints();
	add_ambs_location_constraints();
	add_base_cap_constraints();

	model.update();	

}



void FutureAmbulance::solve(){
    model.update();

    model.set(GRB_IntParam_OutputFlag,0);
	model.set(GRB_IntParam_Presolve, 0);
	model.set(GRB_IntParam_Threads, 1);

	model.optimize();
	int status = model.get(GRB_IntAttr_Status);
	if(status == GRB_OPTIMAL){
		obj = model.get(GRB_DoubleAttr_ObjVal);
	}else{
		std::cout << "ERROR: AMBULANCE MODEL INFEASIBLE!\n";
		model.computeIIS();
		model.write("future.ilp");
		exit(1);
	}
}

void FutureAmbulance::add_bases_constraints_t0(){
	std::stringstream name;
	GRBLinExpr con = 0;

	//0.5 - Flow at bases, t = 0
	for(int a = 0; a < data.types_amb; ++a){
		for(int b = 0; b < data.num_bases; ++b){
			con += A0_ab[a][b];
			if(a0 == a && b_ret != -1 && b == b_ret &&
				data.L(0,a,data.hospitals[h0], data.bases[b_ret]) == data.bases[b_ret]){
				con += 1;
			}

			for(int l1 = 0; l1 < data.total_locals; ++l1){
				if(data.L(0,a,l1,data.bases[b]) == data.bases[b] &&
					(data.L_tab[0][a][b].find(l1) != data.L_tab[0][a][b].end() ||
						data.is_hospital[l1])){
					con += A0_alb[a][l1][b];
				}
			}

			name << "flow_bases_a" << a << "_b" << b;
			model.addConstr(con, GRB_EQUAL, At_ab[0][a][b], name.str());
			name.str("");
			con = 0;
		}
	}
}

void FutureAmbulance::add_hospitals_constraints_t0(){
	std::stringstream name;
	GRBLinExpr con = 0;
	
	// 0.6 - Flow at hospitals
	for(int a = 0; a < data.types_amb; ++a){
		for(int h = 0; h < data.num_hosps; ++h){
			if(h == h0 && a == a0){
				con = 1; //either the ambulance leave to base or to another call
				name << "flow_hospitals_a" << a << "_h" << h;
				model.addConstr(con, GRB_EQUAL, A0_ah[a][h], name.str());
				name.str("");
				con = 0;
			}else if(A0_ah[a][h] > 0){
				fmt::print("WEIRD ====> a{} h{} = {}, return is ({},{})\n",a,h,A0_ah[a][h],
					a0,h0);
				exit(1);
			}
		}
	}
}

void FutureAmbulance::add_locations_constraints_t0(){
	std::stringstream name;
	GRBLinExpr con = 0;
	// 0.7 - Fluxo at location between hospitals and bases
	for(int a = 0; a < data.types_amb; ++a){
		for(int b = 0; b < data.num_bases; ++b){
			for(auto l1: data.L_tab[0][a][b]){
				if(b_ret != -1 && data.L(0,a,data.hospitals[h0], data.bases[b_ret]) == l1){
					con += 1;
				}

				for(auto l2: data.L_tab[0][a][b]){
					if(data.L(0,a,l2, data.bases[b]) == l1){
						con += A0_alb[a][l2][b];
					}
				}
				name << "flow_locals_a" << a << "_b" << b << "_l" << l1;
				model.addConstr(con,GRB_EQUAL,At_alb[0][a][l1][b],
					name.str());
				name.str("");
				con = 0;
			}
		}
	}
}

void FutureAmbulance::add_queues_constraints_t0(){
	std::stringstream name;
	GRBLinExpr con = 0;

	// 0.8 - Flow at queues
	for(int c = 0; c < data.types_call; ++c){
		for(int l = 0; l < data.num_locals; ++l){
			con += C[c][l];
			if(l == l0 && c == c0 && b_ret == -1){
				con -= 1;
			}

			name << "flow_queues_c" << c << "_l" << l;
			model.addConstr(con, GRB_EQUAL, Ct_cl[0][c][l],name.str());
			name.str("");
			con = 0;
		}
	}
}

void FutureAmbulance::add_bases_constraints(){
	std::stringstream name;
	GRBLinExpr con = 0;

	// 0.9 - Flow at bases,  t > 0
	for(int t = 0; t < data.num_times-1; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			for(int b = 0; b < data.num_bases; ++b){
				con += At_ab[t][a][b];
				for(int h = 0; h < data.num_hosps; ++h){
					if(data.L(t,a, data.hospitals[h],data.bases[b]) == data.bases[b]){
						con += yt_ahb[t][a][h][b];
					}	
				}
				for(auto l1: data.L_tab[t+1][a][b]){
					if(data.L(t,a,l1,data.bases[b]) == data.bases[b]){
						con += At_alb[t][a][l1][b];
						for(int c = 0; c < data.types_call; ++c){
							for(int l = 0; l < data.num_locals; ++l){
								for(auto h: data.H[c][l]){
									if(data.A[c].find(a) != data.A[c].end()){
										//If using xt0, put t+1 on 1st index
										con -= xt_calblh[t][c][a][l1][b][l][h];	
									}
								}
							}
						}
					}
				}
				for(int c = 0; c < data.types_call; ++c){
					for(int l = 0; l < data.num_locals; ++l){
						for(auto h: data.H[c][l]){
							if(data.A[c].find(a) != data.A[c].end()){
								//If using xt0, put t+1 on 1st index
								con -= xt_cablh[t][c][a][b][l][h];

							}
						}
					}
				}

				name << "flow_bases_t" << t+1 << "_a" << a << "_b" << b;
				try{
					model.addConstr(con, GRB_EQUAL, At_ab[t+1][a][b], name.str());
				}catch(GRBException ex){
					std::cout << "t" << t+1 << "_a" << a << "_b" << b << "\n";
					std::cout << ex.getErrorCode() << ": ";
					std::cout << ex.getMessage() << "\n";
					exit(1);
				}
				name.str("");
				con = 0;
			} 
		}
	}
}

void FutureAmbulance::add_hospitals_constraints(){
	std::stringstream name;
	GRBLinExpr con = 0;


	int arrival_time = INT_MAX;
	if(call != NULL && b_ret == -1){
		arrival_time = data.tao[0][c0][a0][data.hospitals[h0]][l0][h_dest];
	}

	// 0.10 - Flow at hospitals, t > 0
	for(int t = 0; t < data.num_times-1; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			for(int h = 0; h < data.num_hosps; ++h){

				for(int c = 0; c < data.types_call; ++c){
					for(int l = 0; l < data.num_locals; ++l){
						for(auto h1: data.H[c][l]){
							if(data.A[c].find(a) != data.A[c].end()){
								con -= xt_cahlh[t][c][a][h][l][h1];
							}
						}
					}
				}

				for(int b = 0; b < data.num_bases; ++b){
					con -= yt_ahb[t][a][h][b];
				}

				for(int c = 0; c < data.types_call; ++c){
					for(int l1 = 0; l1 < data.total_locals; ++l1){
						for(int l = 0; l < data.num_locals; ++l){
							if((t)*data.quantum <= 
							 data.tao[0][c][a][l1][l][h] && 
							 data.tao[0][c][a][l1][l][h] <= (t+1)*data.quantum){
								con += A0_callh[c][a][l1][l][h];
							}	
						}
					}
				}

				for(int c = 0; c < data.types_call; ++c){
					for(int l1 = 0; l1 < data.total_locals; ++l1){
						if((t)*data.quantum <= data.tao_0[c][a][l1][h] && 
							data.tao_0[c][a][l1][h] <= (t+1)*data.quantum){
							con += A0_calh[c][a][l1][h];
						}
					}
				}

				if(a0 == a && b_ret == -1){
					for(int c = 0; c < data.types_call; ++c){
						for(int l = 0; l < data.num_locals; ++l){
							for(int h1 = 0; h1 < data.num_hosps; ++h1){
								if((t)*data.quantum <= arrival_time  && 
									arrival_time <= (t+1)*data.quantum && 
									c0 == c && l0 == l && h0 == h1 && 
									data.H[c][l].find(h) != data.H[c][l].end() &&
									data.A[c].find(a) != data.A[c].end()){
									con += 1;
								}
							}
						}
					}
				}

				for(int c = 0; c < data.types_call; ++c){
					for(int l = 0; l < data.num_locals; ++l){
						if(data.H[c][l].find(h) != data.H[c][l].end()){
							for(int b = 0; b < data.num_bases; ++b){
								for(int t1 = 1; t1 < t+1; ++t1){
									if((t)*data.quantum <= 
										t1*data.quantum+
										data.tao[t1][c][a][data.bases[b]][l][h] &&
										t1*data.quantum+data.tao[t1][c][a][data.bases[b]][l][h] 
										<= (t+1)*data.quantum &&
										data.A[c].find(a) != data.A[c].end()){
										con += xt_cablh[t1-1][c][a][b][l][h];
									}
								}
							}
						}
					}
				}

				for(int c = 0; c < data.types_call; ++c){
					for(int l = 0; l < data.num_locals; ++l){
						if(data.H[c][l].find(h) != data.H[c][l].end()){
							for(int h1 = 0; h1 < data.num_hosps; ++h1){
								for(int t1 = 1; t1 < t+1; ++t1){
									if((t)*data.quantum <= 
										t1*data.quantum+
										data.tao[t1][c][a][data.hospitals[h1]][l][h] &&
										t1*data.quantum+
										data.tao[t1][c][a][data.hospitals[h1]][l][h] 
										<= (t+1)*data.quantum &&
										data.A[c].find(a) != data.A[c].end()){
										con += xt_cahlh[t1-1][c][a][h1][l][h];
									}
								}	
							}
						}
					}
				}

				for(int c = 0; c < data.types_call; ++c){
					for(int l = 0; l < data.num_locals; ++l){
						if(data.H[c][l].find(h) != data.H[c][l].end()){
							for(int l1 = 0; l1 < data.total_locals; ++l1){
								for(int b = 0; b < data.num_bases; ++b){
									for(int t1 = 1; t1 < t+1; ++t1){
										if((t)*data.quantum <= 
											t1*data.quantum+data.tao[t1][c][a][l1][l][h] &&
											t1*data.quantum+data.tao[t1][c][a][l1][l][h] <= 
											(t+1)*data.quantum &&
											data.A[c].find(a) != data.A[c].end()
											&& data.L_tab[t1-1][a][b].find(l1) != 
											data.L_tab[t1-1][a][b].end()){
											con += xt_calblh[t1-1][c][a][l1][b][l][h];
										}
									}
								}
							}
						}
					}
				}

				name << "flow_hospitals_t" << t+1 << "_a" << a << "_h" << h;
				model.addConstr(con, GRB_EQUAL, 0, name.str());
				name.str("");
				con = 0;

			}
		}
	}
}

void FutureAmbulance::add_locations_constraints(){
	std::stringstream name;
	GRBLinExpr con = 0;


	// 0.11 - Flow at locations between hospitals and bases,  t > 0
	for(int t = 0; t < data.num_times-1; ++t){
		for(int a  = 0; a < data.types_amb; ++a){
			for(int b = 0; b < data.num_bases; ++b){
				for(auto l1: data.L_tab[t+1][a][b]){
					for(int h1 = 0; h1 < data.num_hosps; ++h1){
							if(data.L(t,a, data.hospitals[h1], data.bases[b]) == l1){
								con += yt_ahb[t][a][h1][b];
							}
						}

						for(auto l2: data.L_tab[t+1][a][b]){
							if(data.L(t,a,l2, data.bases[b]) == l1){
								con += At_alb[t][a][l2][b];
								for(int c = 0; c < data.types_call; ++c){
									for(int l = 0; l < data.num_locals; ++l){
										for(auto h: data.H[c][l]){
											if(data.A[c].find(a) != data.A[c].end()){
												con -= xt_calblh[t][c][a][l2][b][l][h];	
											}
										}
									}
								}
							}
						}
						name << "flow_locals_t" << t+1 << "_" << a << "_"<<b<<"_" << l1;
						model.addConstr(con,GRB_EQUAL, At_alb[t+1][a][l1][b], name.str());
						name.str("");
						con = 0;
				}
			}
		}
	}
}

void FutureAmbulance::add_queues_constraints(){
	std::stringstream name;
	GRBLinExpr con = 0;

	// 0.12 - Flow at queues, t > 0
	for(int t = 0; t < data.num_times-1; ++t){
		for(int c = 0; c < data.types_call; ++c){
			for(int l = 0; l < data.num_locals; ++l){
				con += Ct_cl[t][c][l];
				con += lambda[t][c][l];
				for(auto a: data.A[c]){
					for(int b = 0; b < data.num_bases; ++b){
						for(auto h: data.H[c][l]){
							con -= xt_cablh[t][c][a][b][l][h];
						}
					}
				}

				for(auto a: data.A[c]){
					for(int h1 = 0; h1 < data.num_hosps; ++h1){
						for(auto h: data.H[c][l]){
							con -= xt_cahlh[t][c][a][h1][l][h];
						}
					}
				}

				for(auto a: data.A[c]){
					for(int b = 0; b < data.num_bases; ++b){
						for(auto l1: data.L_tab[t+1][a][b]){
							for(auto h: data.H[c][l]){
								con -= xt_calblh[t][c][a][l1][b][l][h];
							}
						}
					}
				}
				name << "flow_queue_t" << t+1 << "_c" << c << "_l" << l;
				model.addConstr(con,GRB_EQUAL,Ct_cl[t+1][c][l], name.str());
				name.str("");
				con = 0;
			}
		}
	}
}

void FutureAmbulance::add_ambs_location_constraints(){
	GRBLinExpr con;
	std::stringstream name;

	//Location constraints
	for(int t = 0; t < data.num_times-1; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			for(int b = 0; b < data.num_bases; ++b){
				for(auto l1: data.L_tab[t+1][a][b]){
					for(int c = 0; c < data.types_call; ++c){
						for(int l = 0; l < data.num_locals; ++l){
							for(auto h: data.H[c][l]){
								if(data.A[c].find(a) != data.A[c].end()){
									//If using xt0, put t+1 on 1st index
									con += xt_calblh[t][c][a][l1][b][l][h];
								}
							}
						}
					}
					name << "dispatch_locals_t" << t+1 << "_a" << a << "_b";
					name << b << "_l" << l1;
					model.addConstr(con, GRB_LESS_EQUAL, At_alb[t][a][l1][b],
						name.str());
					name.str("");
					con = 0;
				}
			}
		}
	}


	// Max of ambulance at bases
	for(int t = 0; t < data.num_times-1; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			for(int b = 0; b < data.num_bases; ++b){
				for(int c = 0; c < data.types_call; ++c){
					for(int l = 0; l < data.num_locals; ++l){
						for(auto h: data.H[c][l]){
							if(data.A[c].find(a) != data.A[c].end()){
								//If using xt0, put t+1 on 1st index
								con += xt_cablh[t][c][a][b][l][h];
							}
						}
					}
				}

				name << "dispatch_bases_t" << t+1 << "_a" << a << "_b" << b;
				model.addConstr(con, GRB_LESS_EQUAL, At_ab[t][a][b], name.str());
				name.str("");
				con = 0;
			}
		}
	}
}

void FutureAmbulance::add_base_cap_constraints(){
	GRBLinExpr con;
	std::stringstream name;

	for(int t = 0; t < data.num_times; ++t){
		for(int b = 0; b < data.num_bases; ++b){
			for(int a = 0; a < data.types_amb; ++a){
				con += At_ab[t][a][b];
			}
			name << "cap_bases_" << t+1 << "_" << b;
			model.addConstr(con,GRB_LESS_EQUAL,data.cap_bases[b], name.str());
			name.str("");
			con = 0;
		}
	}
}

FutureAmbulance::~FutureAmbulance(){
	for(int t = 0; t < data.num_times-1; ++t){
		for(int c = 0; c < data.types_call; ++c){
			delete[] lambda[t][c];
		}
		delete[] lambda[t];
	}

	for(int t = 0; t < data.num_times-1; ++t){
		for(int c = 0; c < data.types_call; ++c){
			for(int a = 0; a < data.types_amb; ++a){
				for(int b = 0; b < data.num_bases; ++b){
					for(int l = 0; l < data.num_locals; ++l){
						delete[] xt_cablh[t][c][a][b][l];
					}
					delete[] xt_cablh[t][c][a][b];
				}
				delete[] xt_cablh[t][c][a];
			}
			delete[] xt_cablh[t][c];
		}
		delete[] xt_cablh[t];
	}


	for(int t = 0; t < data.num_times-1; ++t){
		for(int c = 0; c < data.types_call; ++c){
			for(int a = 0; a < data.types_amb; ++a){
				for(int h1 = 0; h1 < data.num_hosps; ++h1){
					for(int l = 0; l < data.num_locals; ++l){
						delete[] xt_cahlh[t][c][a][h1][l];
					}
					delete[] xt_cahlh[t][c][a][h1];
				}
				delete[] xt_cahlh[t][c][a];
			}
			delete[] xt_cahlh[t][c];
		}
		delete[] xt_cahlh[t];
	}
	//If using xt0, loop until num_times instead of num_times-1
	for(int t = 0; t < data.num_times-1; ++t){
		for(int c = 0; c < data.types_call; ++c){
			for(int a = 0; a < data.types_amb; ++a){
				for(int l1 = 0; l1 < data.total_locals; ++l1){
					for(int b = 0; b < data.num_bases; ++b){
						for(int l = 0; l < data.num_locals; ++l){
							delete[] xt_calblh[t][c][a][l1][b][l];
						}
						delete[] xt_calblh[t][c][a][l1][b];
					}
					delete[] xt_calblh[t][c][a][l1];
				}
				delete[] xt_calblh[t][c][a];
			}
			delete[] xt_calblh[t][c];
		}
		delete[] xt_calblh[t];
	}

	for(int t = 0; t < data.num_times-1; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			for(int h = 0; h < data.num_hosps; ++h){
				delete[] yt_ahb[t][a][h];
			}
			delete[] yt_ahb[t][a];
		}
		delete[] yt_ahb[t];
	}

	for(int t = 0; t < data.num_times; ++t){
		for(int c = 0; c < data.types_call; ++c){
			delete[] Ct_cl[t][c];
		}
		delete[] Ct_cl[t];
	}

	for(int t = 0; t < data.num_times; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			delete[] At_ab[t][a];
		}
		delete[] At_ab[t];
	}

	for(int t = 0; t < data.num_times; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			for(int l = 0; l < data.total_locals; ++l){
				delete[] At_alb[t][a][l];
			}
			delete[] At_alb[t][a];
		}
		delete[] At_alb[t];
	}



	for(int a = 0; a < data.types_amb; ++a){
		for(int l1 = 0; l1 < data.total_locals; ++l1){
			delete[] A0_alb[a][l1];
		}
		delete[] A0_ab[a];
		delete[] A0_ah[a];
		delete[] A0_alb[a];
	}

	for(int c = 0; c < data.types_call; ++c){
		for(int a = 0; a < data.types_amb; ++a){
			for(int l1 = 0; l1 <  data.total_locals; ++l1){
				delete[]  A0_calh[c][a][l1];
			}
			delete[]  A0_calh[c][a];
			for(int l1 = 0; l1 < data.total_locals; ++l1){
				for(int l = 0; l < data.num_locals; ++l){
					delete[] A0_callh[c][a][l1][l];
				}
				delete[]  A0_callh[c][a][l1];
			}
			delete[]  A0_callh[c][a];
		}
		delete[]  A0_calh[c];
		delete[]  A0_callh[c];
		delete[] C[c];
	}
	delete[] A0_ab; delete[] A0_ah; delete[] A0_alb; delete[] A0_calh;
	delete[] A0_callh; delete[] C;

	delete[] xt_cablh; delete[] xt_cahlh; delete[] xt_calblh;
	delete[] Ct_cl; delete[] At_ab; delete[] At_alb; delete[] yt_ahb;
	delete[] lambda;

}

