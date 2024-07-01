#include "../include/full_model_det.h"

CallModel1::CallModel1(Data & data, GRBEnv & env, Solver& solver):
	env(env), data(data), model(env), type_call0(data.type_call0), 
	local_call0(data.local_call0), solver(solver){

	src_type = -1;
	src_location = -1;
	amb_type = -1;
	loc_base = -1;
	dst_hosp = -1;

	A0_ab = new int*[data.types_amb];
	A0_ah = new int*[data.types_amb];
	A0_alb = new int**[data.types_amb];
	A0_calh = new int***[data.types_call];
	A0_callh = new int****[data.types_call];
	C = new int*[data.types_call];

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

	x0_abh = new GRBVar**[data.types_amb];
	for(int a = 0; a < data.types_amb; ++a){
		x0_abh[a] = new GRBVar*[data.num_bases];
		for(int b = 0; b < data.num_bases; ++b){
			x0_abh[a][b] = new GRBVar[data.num_hosps];
			for(int h = 0; h < data.num_hosps; ++h){
				if(data.A[type_call0].find(a) != data.A[type_call0].end() &&
					data.H[type_call0][local_call0].find(h) != 
					data.H[type_call0][local_call0].end()){

					name << "x0_a" << a << "_b" << b << "_h"  << h;
					x0_abh[a][b][h] = model.addVar(0, A0_ab[a][b], 0,
						(data.relax["x0_abh"] ? GRB_CONTINUOUS : GRB_INTEGER), 
						name.str());
					name.str("");
				}
			}
		}
	}

	x0_ahh = new GRBVar**[data.types_amb];
	for(int a = 0; a < data.types_amb; ++a){
		x0_ahh[a] = new GRBVar*[data.num_hosps];
		for(int h1 = 0; h1 < data.num_hosps; ++h1){
			x0_ahh[a][h1] = new GRBVar[data.num_hosps];
			for(int h = 0; h < data.num_hosps; ++h){
				if(data.A[type_call0].find(a) != data.A[type_call0].end() &&
					data.H[type_call0][local_call0].find(h) != 
					data.H[type_call0][local_call0].end()){

					name << "x0_a" << a << "_h'" << h1 << "_h" << h;
					x0_ahh[a][h1][h] = model.addVar(0, A0_ah[a][h1], 0,
						(data.relax["x0_ahh"] ? GRB_CONTINUOUS : GRB_INTEGER), 
						name.str());
					name.str("");
				}
				
			}
		}
	}

	x0_albh = new GRBVar***[data.types_amb];
	for(int a = 0; a < data.types_amb; ++a){
		x0_albh[a] = new GRBVar**[data.total_locals];
		for(int l = 0; l < data.total_locals; ++l){
			x0_albh[a][l] = new GRBVar*[data.num_bases];
			for(int b = 0; b < data.num_bases; ++b){
				x0_albh[a][l][b] = new GRBVar[data.num_hosps];
				for(int h = 0; h < data.num_hosps; ++h){
					if(data.A[type_call0].find(a) != data.A[type_call0].end() &&
						data.H[type_call0][local_call0].find(h) !=
						data.H[type_call0][local_call0].end() &&
						data.L_tab[0][a][b].find(l) != 
						data.L_tab[0][a][b].end()){

						name << "x0_a" << a << "_l" << l << "_b" << b << "_h" << h;
						x0_albh[a][l][b][h] = model.addVar(0, A0_alb[a][l][b],
							0, (data.relax["x0_albh"] ? GRB_CONTINUOUS : 
							GRB_INTEGER), name.str());
						name.str("");
					}
					
				}
			}
		}
	}

	y0_ahb = new GRBVar**[data.types_amb];
	for(int a = 0; a < data.types_amb; ++a){
		y0_ahb[a] = new GRBVar*[data.num_hosps];
		for(int h  = 0; h < data.num_hosps; ++h){
			y0_ahb[a][h] = new GRBVar[data.num_bases];
			for(int  b = 0; b < data.num_bases; ++b){
				name << "y0" << "_a" << a << "_h" << h << "_b" << b;
				y0_ahb[a][h][b] = model.addVar(0,data.num_ambs, 0,
					(data.relax["y0_ahb"] ? GRB_CONTINUOUS : GRB_INTEGER), 
					name.str());
				name.str("");
			}
		}
	}

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
	double factor = 0;
	for(int a = 0; a < data.types_amb; ++a){
		if(data.A[type_call0].find(a) == data.A[type_call0].end())
			continue;
		factor = 1;
		for(int b = 0; b < data.num_bases; ++b){
			for(int h = 0; h < data.num_hosps; ++h){
				if(data.H[type_call0][local_call0].find(h) != 
					data.H[type_call0][local_call0].end()){
					fo += f(0,type_call0, a, data.bases[b], local_call0, 
						h, data, factor) * 
						x0_abh[a][b][h];
				}
			}
		}

		for(int h1 = 0; h1 < data.num_hosps; ++h1){
			for(int h = 0; h < data.num_hosps; ++h){
				if(data.H[type_call0][local_call0].find(h) != 
					data.H[type_call0][local_call0].end()){
					fo += f(0,type_call0, a, data.hospitals[h1], local_call0, 
						h, data, factor)*
						x0_ahh[a][h1][h];	
				}
			}
		}

		for(int l1 = 0; l1 < data.total_locals; ++l1){
			for(int b = 0; b < data.num_bases; ++b){
				for(int h = 0; h < data.num_hosps; ++h){
					if(data.H[type_call0][local_call0].find(h) != 
					data.H[type_call0][local_call0].end() &&
					data.L_tab[0][a][b].find(l1) != data.L_tab[0][a][b].end()){
						fo+= f(0,type_call0,a,l1,data.bases[b],
							local_call0, h,data,factor)*
							x0_albh[a][l1][b][h];
					}
				}	
			}
		}
	}

	factor = 0;
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

		// for(int a = 0; a < data.types_amb; ++a){
		// 	for(int h = 0; h < data.num_hosps; ++h){
		// 		for(int b = 0; b < data.num_bases; ++b){
		// 			fo += data.times[data.hospitals[h]][data.bases[b]]*
		// 				yt_ahb[t-1][a][h][b];
		// 		}
		// 	}
		// }
	}

	for(int t = 0; t < data.num_times; ++t){
		for(int c = 0; c < data.types_call; ++c){
			factor = data.ins.penalties[c]*36000;
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


void CallModel1::solve(){
	model.update();
	model.set(GRB_IntParam_OutputFlag,0);
	model.set(GRB_IntParam_Presolve, 0);
	model.set(GRB_IntParam_Threads, 1);

	// std::cout << "Number of variables: " << model.get(GRB_IntAttr_NumVars) << "\n";
	model.optimize();

	// std::cout << "Z (DET): " << model.get(GRB_DoubleAttr_ObjVal) << "\n";

	// if(sim.events.top()->id == 10){
	// 	GRBVar* vars = model.getVars();
	// 	for(int i = 0; i < model.get(GRB_IntAttr_NumVars); ++i){
	// 		if(vars[i].get(GRB_DoubleAttr_X) > g_params.g_params.EPS){
	// 			std::cout << vars[i].get(GRB_StringAttr_VarName) << " = ";
	// 			std::cout << vars[i].get(GRB_DoubleAttr_X) << " (";
	// 			std::cout << vars[i].get(GRB_DoubleAttr_Obj) << ")\n";
	// 		}
	// 	}
	// 	delete[] vars;
	// }

	src_type = -1;
	src_location = -1;
	amb_type = -1;
	loc_base = -1;
	dst_hosp = -1;
	src_hosp = -1;
	base_return = -1;

	bool found = false;

	for(int a = 0; a < data.types_amb; ++a){
		for(int b = 0; b < data.num_bases; ++b){
			for(int h = 0; h < data.num_hosps; ++h){
				if(data.A[type_call0].find(a) != data.A[type_call0].end() &&
					data.H[type_call0][local_call0].find(h) !=
					data.H[type_call0][local_call0].end() && 
					x0_abh[a][b][h].get(GRB_DoubleAttr_X) > 0.5){
					src_type = BASE;
					src_location = b;
					dst_hosp = h;
					amb_type = a;
					found = true;
					// std::cout << x0_abh[a][b][h].get(GRB_StringAttr_VarName) << " ";
					// std::cout << x0_abh[a][b][h].get(GRB_DoubleAttr_X) << "\n";
					break;
				}
			}
			if(found) break;
		}
		if(found) break;
	}

	for(int a = 0; a < data.types_amb; ++a){
		for(int h1 = 0; h1 < data.num_hosps; ++h1){
			for(int h = 0; h < data.num_hosps; ++h){
				if(data.A[type_call0].find(a) != data.A[type_call0].end() &&
					data.H[type_call0][local_call0].find(h) !=
					data.H[type_call0][local_call0].end() && 
					x0_ahh[a][h1][h].get(GRB_DoubleAttr_X) > 0.5){
					src_type = HOSP;
					src_location = h1;
					dst_hosp = h;
					amb_type = a;
					found = true;
					// std::cout << x0_ahh[a][h1][h].get(GRB_StringAttr_VarName) << " ";
					// std::cout << x0_ahh[a][h1][h].get(GRB_DoubleAttr_X) << "\n";
					break;
				}
			}
			if(found) break;
		}
		if(found) break;
	}

	for(int a = 0; a < data.types_amb; ++a){
		for(int l1 = 0; l1 < data.total_locals; ++l1){
			for(int b = 0; b < data.num_bases; ++b){
				for(int h = 0; h < data.num_hosps; ++h){
					if(data.A[type_call0].find(a) != data.A[type_call0].end() &&
						data.H[type_call0][local_call0].find(h) !=
						data.H[type_call0][local_call0].end() &&
						data.L_tab[0][a][b].find(l1) != 
						data.L_tab[0][a][b].end() &&
						x0_albh[a][l1][b][h].get(GRB_DoubleAttr_X) > 0.5){
						src_type = LOC;
						src_location = l1;
						dst_hosp = h;
						amb_type = a;
						loc_base = b;

						// std::cout << x0_albh[a][l1][b][h].get(GRB_StringAttr_VarName) << " ";
						// std::cout << x0_albh[a][l1][b][h].get(GRB_DoubleAttr_X) << "\n";
						found = true;
					}
				}
				if(found) break;
			}
			if(found) break;
		}
		if(found) break;
	}
}


void CallModel1::add_bases_constraints_t0(){
	GRBLinExpr con;
	std::stringstream name;

	model.update();

	for(int a = 0; a < data.types_amb; ++a){
		for(int b = 0; b < data.num_bases; ++b){
			con += A0_ab[a][b];
			for(int h = 0; h < data.num_hosps; ++h){
				if(data.L(0,a, data.hospitals[h],data.bases[b]) == data.bases[b]){
					con += y0_ahb[a][h][b];
				}
			}
			for(int l1 = 0; l1 < data.total_locals; ++l1){
				if(data.L_tab[0][a][b].find(l1) != data.L_tab[0][a][b].end() &&
					data.L(0,a,l1,data.bases[b]) == data.bases[b]){
					con += A0_alb[a][l1][b];
					for(int h = 0; h < data.num_hosps; ++h){
						if(data.H[type_call0][local_call0].find(h) != 
							data.H[type_call0][local_call0].end() &&
							data.A[type_call0].find(a) != data.A[type_call0].end()){
							con -= x0_albh[a][l1][b][h];
						}
					}
				}
			}
			
			for(int h = 0; h < data.num_hosps; ++h){
				if(data.H[type_call0][local_call0].find(h) !=
					data.H[type_call0][local_call0].end() &&
					data.A[type_call0].find(a) != data.A[type_call0].end()){
					con -= x0_abh[a][b][h];
				}
			}

			name << "flow_base_a" << a  << "_b" << b;
			model.addConstr(con, GRB_EQUAL, At_ab[0][a][b], name.str());
			name.str("");
			con = 0;
		}
	}
}

void CallModel1::add_hospitals_constraints_t0(){
	GRBLinExpr con;
	std::stringstream name;


	for(int a = 0; a < data.types_amb; ++a){
		for(int h = 0; h < data.num_hosps; ++h){
			for(int h1 = 0; h1 < data.num_hosps; ++h1){
				if(data.H[type_call0][local_call0].find(h1) != 
					data.H[type_call0][local_call0].end() && 
					data.A[type_call0].find(a) != data.A[type_call0].end()){
					con += x0_ahh[a][h][h1];

				}
			}
			for(int b = 0; b < data.num_bases; ++b){
				con += y0_ahb[a][h][b];
			}

			name << "flow_hospital_a" << a << "_h" << h;
			model.addConstr(con, GRB_EQUAL, 
				A0_ah[a][h], name.str());
			name.str("");
			con = 0;
		}
	}
}

void CallModel1::add_locations_constraints_t0(){
	GRBLinExpr con;
	std::stringstream name;

	for(int a = 0; a < data.types_amb; ++a){
		for(int b = 0; b < data.num_bases; ++b){
			for(int l = 0; l < data.total_locals; ++l){
				if(data.L_tab[0][a][b].find(l) != data.L_tab[0][a][b].end()){
					for(int h = 0; h < data.num_hosps; ++h){
						if(data.L(0,a,data.hospitals[h],data.bases[b]) == l){
							con += y0_ahb[a][h][b];
						}
					}
					for(int l1 = 0; l1 < data.total_locals; ++l1){
						if(data.L_tab[0][a][b].find(l1) != 
							data.L_tab[0][a][b].end() &&
							data.L(0,a,l1, data.bases[b]) == l){
							con += A0_alb[a][l1][b];
							for(int h = 0; h < data.num_hosps; ++h){
								if(data.H[type_call0][local_call0].find(h) !=
									data.H[type_call0][local_call0].end() &&
									data.A[type_call0].find(a) != 
									data.A[type_call0].end()){
									con -= x0_albh[a][l1][b][h];
								}
							}
						}
					}
					name << "flow_locals_a" << a << "_b" << b  << "_l" << l;
					model.addConstr(con,GRB_EQUAL, 
						At_alb[0][a][l][b], name.str());

					name.str("");
					con = 0;
				}
			}
		}
	}
}

void CallModel1::add_queues_constraints_t0(){
	GRBLinExpr con;
	std::stringstream name;


	for(int c = 0; c < data.types_call; ++c){
		for(int l = 0; l < data.num_locals; ++l){
			con += C[c][l];
			if(c == type_call0 && l == local_call0){
				con += 1;
				for(int a = 0; a < data.types_amb; ++a){
					if(data.A[c].find(a) != data.A[c].end()){
						for(int b = 0; b < data.num_bases; ++b){
							for(int h = 0; h < data.num_hosps; ++h){
								if(data.H[c][l].find(h) !=
									data.H[c][l].end() &&
									data.A[type_call0].find(a) != 
									data.A[type_call0].end()){
									con -= x0_abh[a][b][h];
								}
							}
						}
					}
				}
				for(int a  = 0; a < data.types_amb; ++a){
					if(data.A[c].find(a) != data.A[c].end()){
						for(int h1 = 0; h1 < data.num_hosps; ++h1){
							for(int h = 0; h < data.num_hosps; ++h){
								if(data.H[c][l].find(h) !=
									data.H[c][l].end() &&
									data.A[type_call0].find(a) != 
									data.A[type_call0].end()){
									con -= x0_ahh[a][h1][h];
								}
							}
						}
					}
				}

				for(int a = 0; a < data.types_amb; ++a){
					if(data.A[c].find(a) != data.A[c].end()){
						for(int b = 0; b < data.num_bases; ++b){
							for(int l1 = 0; l1 < data.total_locals; ++l1){
								for(int h = 0; h < data.num_hosps; ++h){
									if(data.H[c][l].find(h) !=
										data.H[c][l].end() &&
										data.A[type_call0].find(a) != 
										data.A[type_call0].end() && 
										data.L_tab[0][a][b].find(l1) != 
										data.L_tab[0][a][b].end()){
										con -= x0_albh[a][l1][b][h];
									}
								}
							}
						}
					}
				}
			}

			name << "flow_queue_c" << c << "_l" << l;
			
			model.addConstr(con, GRB_EQUAL, 
				Ct_cl[0][c][l], name.str());
			name.str("");
			con = 0;
		}
	}
}

void CallModel1::add_bases_constraints(){
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
								for(int h = 0; h < data.num_hosps; ++h){
									if(data.H[c][l].find(h) !=
										data.H[c][l].end() &&
										data.A[c].find(a) != data.A[c].end()){
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
						for(int h = 0; h < data.num_hosps; ++h){
							if(data.H[c][l].find(h) !=
								data.H[c][l].end() &&
								data.A[c].find(a) != data.A[c].end()){
								//If using xt0, put t+1 on 1st index
								con -= xt_cablh[t][c][a][b][l][h];

							}
						}
					}
				}

				name << "flow_bases_t" << t+1 << "_a" << a << "_b" << b;
				try{
					model.addConstr(con, GRB_EQUAL, 
						At_ab[t+1][a][b], name.str());
				}catch(GRBException ex){
					std::cout << ex.getErrorCode() << "::: ";
					std::cout << ex.getMessage() << "\n";
					exit(1);
				}
				name.str("");
				con = 0;
			} 
		}
	}
}

void CallModel1::add_hospitals_constraints(){
	GRBLinExpr	con;
	std::stringstream name;
	int quantum = data.quantum;

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

				for(int b = 0; b < data.num_bases; ++b){
					if((t)*quantum <= data.tao[0][type_call0][a][data.bases[b]][local_call0][h] 
						&& data.tao[0][type_call0][a][data.bases[b]][local_call0][h] <= 
						(t+1)*quantum &&
						data.H[type_call0][local_call0].find(h) != 
						data.H[type_call0][local_call0].end() &&
						data.A[type_call0].find(a) != data.A[type_call0].end()){
						con += x0_abh[a][b][h];
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

				for(int h1 = 0; h1 < data.num_hosps; ++h1){
					if((t)*quantum <= 
						data.tao[0][type_call0][a][data.hospitals[h1]][local_call0][h] 
						&& data.tao[0][type_call0][a][data.hospitals[h1]][local_call0][h] 
						<= (t+1)*quantum &&
						data.H[type_call0][local_call0].find(h) != 
						data.H[type_call0][local_call0].end() &&
						data.A[type_call0].find(a) != data.A[type_call0].end()){
						con += x0_ahh[a][h1][h];
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

				for(int l1 = 0; l1 < data.total_locals; ++l1){
					for(int b = 0; b < data.num_bases; ++b){
						if((t)*quantum <= data.tao[0][type_call0][a][l1][local_call0][h] && 
							data.tao[0][type_call0][a][l1][local_call0][h] <= (t+1)*quantum &&
							data.H[type_call0][local_call0].find(h) != 
							data.H[type_call0][local_call0].end() &&
							data.A[type_call0].find(a) != 
							data.A[type_call0].end() &&
							data.L_tab[0][a][b].find(l1) != 
							data.L_tab[0][a][b].end()){
							con += x0_albh[a][l1][b][h];
						}
					}
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
void CallModel1::add_locations_constraints(){
	GRBLinExpr con;
	std::stringstream name;

	for(int t = 0; t < data.num_times-1; ++t){
		for(int a  = 0; a < data.types_amb; ++a){
			for(int b = 0; b < data.num_bases; ++b){
				for(int l1 = 0; l1 < data.total_locals; ++l1){
					if(data.L_tab[t+1][a][b].find(l1) !=
					 data.L_tab[t+1][a][b].end()){
						for(int h1 = 0; h1 < data.num_hosps; ++h1){
							if(data.L(t,a, data.hospitals[h1],data.bases[b]) == l1){
								con += yt_ahb[t][a][h1][b];
							}
						}

						for(int l2 = 0; l2 < data.total_locals; ++l2){
							if(data.L_tab[t+1][a][b].find(l2) != 
								data.L_tab[t+1][a][b].end() &&
								data.L(t,a,l2, data.bases[b]) == l1){
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
}
void CallModel1::add_queues_constraints(){
	GRBLinExpr con;
	std::stringstream name;

	for(int t = 0; t < data.num_times-1; ++t){
		for(int c = 0; c < data.types_call; ++c){
			for(int l = 0; l < data.num_locals; ++l){
				con += Ct_cl[t][c][l];
				con += data.lambda[t][c][l];
				for(int a = 0; a < data.types_amb; ++a){
					if(data.A[c].find(a) != data.A[c].end()){
						for(int b = 0; b < data.num_bases; ++b){
							for(int h = 0; h < data.num_hosps; ++h){
								if(data.H[c][l].find(h) != data.H[c][l].end() &&
									data.A[c].find(a) != data.A[c].end()){
									//If using xt0, put t+1 on 1st index
									con -= xt_cablh[t][c][a][b][l][h];
								}
							}
						}
					}
				}

				for(int a = 0; a < data.types_amb; ++a){
					if(data.A[c].find(a) != data.A[c].end()){
						for(int h1 = 0; h1 < data.num_hosps; ++h1){
							for(int h = 0; h < data.num_hosps; ++h){
								if(data.H[c][l].find(h) != data.H[c][l].end() &&
									data.A[c].find(a) != data.A[c].end()){
									//If using xt0, put t+1 on 1st index
									con -= xt_cahlh[t][c][a][h1][l][h];
								}
							}
						}
					}
				}

				for(int a = 0; a < data.types_amb; ++a){
					if(data.A[c].find(a) != data.A[c].end()){
						for(int b = 0; b < data.num_bases; ++b){
							for(int l1 = 0; l1 < data.total_locals; ++l1){
								for(int h = 0; h < data.num_hosps; ++h){
									if(data.H[c][l].find(h) != data.H[c][l].end() &&
										data.A[c].find(a) != data.A[c].end() &&
										data.L_tab[t+1][a][b].find(l1) != 
										data.L_tab[t+1][a][b].end()){
										//If using xt0, put t+1 on 1st index
										con -= xt_calblh[t][c][a][l1][b][l][h];
									}
								}
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
void CallModel1::add_ambs_location_constraints(){
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

	// model.update();
	// for(int a = 0; a < data.types_amb; ++a){
	// 	if(data.A[type_call0].find(a) == data.A[type_call0].end())
	// 		continue;
	// 	for(int b = 0; b < data.num_bases; ++b){
	// 		for(int h =  0; h < data.num_hosps; ++h){
	// 			if(data.H[type_call0][local_call0].find(h) ==
	// 				data.H[type_call0][local_call0].end())
	// 					continue;
	// 			con += x0_abh[a][b][h];
	// 		}
	// 	}

	// 	for(int h = 0; h <  data.num_hosps; ++h){
	// 		for(int h1 = 0; h1 < data.num_hosps;++h1){
	// 			if(data.H[type_call0][local_call0].find(h1) ==
	// 				data.H[type_call0][local_call0].end())
	// 				continue;

	// 			con += x0_ahh[a][h][h1];
	// 		}
	// 	}

	// 	for(int l1 = 0; l1 < data.total_locals; ++l1){
	// 		for(int b = 0; b < data.num_bases; ++b){
	// 			for(int h = 0; h < data.num_hosps; ++h){
	// 				if(data.H[type_call0][local_call0].find(h) ==
	// 				data.H[type_call0][local_call0].end() ||
	// 				data.L_tab[0][a][b].find(l1) == 
	// 				data.L_tab[0][a][b].end())
	// 					continue;

	// 				con += x0_albh[a][l1][b][h];					
	// 			}
	// 		}
	// 	}
	// }

	// name << "dispatch_call_c" << type_call0 << "_l" << local_call0;
	// model.addConstr(con,GRB_LESS_EQUAL,1,name.str());
	// name.str("");
}

void CallModel1::add_base_cap_constraints(){
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

void CallModel1::debug(){
	int temp = 0;
	for(int a  = 0; a < data.types_amb; ++a){
		for(int b = 0; b < data.num_bases; ++b){
			std::cout << "A0_a" << a << "_b" << b <<":  ";
			std::cout << A0_ab[a][b] << "\n";
			temp += A0_ab[a][b];
		}
	}
	std::cout << "TEMP: " << temp << "============\n";

	for(int a = 0; a < data.types_amb; ++a){
		for(int h = 0; h < data.num_hosps; ++h){
			std::cout << "A0_a" << a << "_h" << h << ": ";
			std::cout << A0_ah[a][h] << "\n";
			temp += A0_ah[a][h];
		}
	}

	for(int c = 0; c < data.types_call; ++c){
		for(int a = 0; a < data.types_amb; ++a){
			for(int l1 = 0; l1 < data.total_locals; ++l1){
				for(int h = 0; h < data.num_hosps; ++h){
					if(A0_calh[c][a][l1][h] != 0){
						std::cout << "A0_c" <<c << "_a" << a << "_l'" << l1;
						std::cout << "_h" << h << ": ";
						std::cout << A0_calh[c][a][l1][h] << "\n";
						temp += A0_calh[c][a][l1][h];
					}
				}
			}
		}
	}
	std::cout << "TEMP: " << temp << "==============\n";
	for(int c = 0; c < data.types_call; ++c){
		for(int a = 0; a < data.types_amb; ++a){
			for(int l1 = 0; l1 < data.total_locals; ++l1){
				for(int l = 0; l < data.num_locals; ++l){
					for(int h = 0; h < data.num_hosps; ++h){
						if(A0_callh[c][a][l1][l][h] != 0){
							std::cout << "A0_c" <<c << "_a" << a << "_l'" << l1;
							std::cout << "_l" << l << "_h" << h << ": ";
							std::cout << A0_callh[c][a][l1][l][h] << "\n";
							temp += A0_callh[c][a][l1][l][h];
						}
					}	
				}
			}
		}
	}
	std::cout << "TEMP: " << temp << "==============\n";

	std::cout << "INITIAL QUEUE:\n";
	for(int c = 0; c < data.types_call; ++c){
		for(int l = 0; l < data.num_locals; ++l){
			std::cout << C[c][l] << " ";
		}
		std::cout << "\n";
	}
}


CallModel1::~CallModel1(){
	for(int a = 0; a < data.types_amb; ++a){
		for(int b = 0; b < data.num_bases; ++b){
			delete[] x0_abh[a][b];
		}
		delete[] x0_abh[a];
	}
	for(int a = 0; a < data.types_amb; ++a){
		for(int h1 = 0; h1 < data.num_hosps; ++h1){
			delete[] x0_ahh[a][h1];
		}
		delete[] x0_ahh[a];
	}
	for(int a = 0; a < data.types_amb; ++a){
		for(int l = 0; l < data.total_locals; ++l){
			for(int b = 0; b < data.num_bases; ++b){
				delete[] x0_albh[a][l][b];
			}
			delete[] x0_albh[a][l];
		}
		delete[] x0_albh[a];
	}

	for(int a = 0; a < data.types_amb; ++a){
		for(int h  = 0; h < data.num_hosps; ++h){
			delete[] y0_ahb[a][h];
		}
		delete[] y0_ahb[a];
	}
	//If using xt0, loop until num_times instead of num_times-1
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

	delete[] x0_abh; delete[] x0_ahh; delete[] x0_albh; delete[] y0_ahb;
	delete[] xt_cablh; delete[] xt_cahlh; delete[] xt_calblh;
	delete[] Ct_cl; delete[] At_ab; delete[] At_alb; delete[] yt_ahb;
}


// //=================================================================


AmbulanceModel::AmbulanceModel(Data & data, GRBEnv & env, Solver& solver):
	data(data), model(env), solver(solver){
	std::stringstream name;
	

	A0_ab = new int*[data.types_amb];
	A0_ah = new int*[data.types_amb];
	A0_alb = new int**[data.types_amb];
	A0_calh = new int***[data.types_call];
	A0_callh = new int****[data.types_call];
	C = new int*[data.types_call];


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


	y0_ahb = new GRBVar**[data.types_amb];
	for(int a = 0; a < data.types_amb; ++a){
		y0_ahb[a] = new GRBVar*[data.num_hosps];
		for(int h = 0; h < data.num_hosps; ++h){
			y0_ahb[a][h] = new GRBVar[data.num_bases];
			for(int  b = 0; b < data.num_bases; ++b){
				name << "y0" << "_a" << a << "_h" << h << "_b" << b;
				y0_ahb[a][h][b] = model.addVar(0, A0_ah[a][h], 0,
					(data.relax["y0_ahb"] ? GRB_CONTINUOUS : GRB_INTEGER), 
					name.str());
				name.str("");
			}
		}
	}

	y0_ahclh = new GRBVar****[data.types_amb];
	for(int a = 0; a < data.types_amb; ++a){
		y0_ahclh[a] = new GRBVar***[data.num_hosps];
		for(int h1 = 0; h1 < data.num_hosps; ++h1){
			y0_ahclh[a][h1] = new GRBVar**[data.types_call];
			for(int c = 0; c < data.types_call; ++c){
				y0_ahclh[a][h1][c] = new GRBVar*[data.num_locals];
				for(int l  = 0; l < data.num_locals; ++l){
					y0_ahclh[a][h1][c][l] = new GRBVar[data.num_hosps];
					for(int h = 0; h < data.num_hosps; ++h){
						if(data.A[c].find(a) != data.A[c].end() &&
							data.H[c][l].find(h) != data.H[c][l].end()){
							name<<"y0_a" << a << "_h'"<<h1<<"_c"<<c<<"_l"<<l<<"_h"<<h;
							y0_ahclh[a][h1][c][l][h] = model.addVar(0,A0_ah[a][h], 0,
								(data.relax["y0_ahclh"] ? GRB_CONTINUOUS : GRB_INTEGER),
								name.str());
							name.str("");	
						}
					}
				}
			}
		}
	}

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
	double factor = 0;

	for(int a = 0; a < data.types_amb; ++a){
		for(int h = 0; h < data.num_hosps; ++h){
			// for(int b = 0; b < data.num_bases; ++b){
			// 	fo += data.times[data.hospitals[h]][data.bases[b]]*y0_ahb[a][h][b];
			// }
			for(int c = 0; c < data.types_call; ++c){
				for(int l = 0; l < data.num_locals; ++l){
					for(int h1 = 0; h1 < data.num_hosps; ++h1){
						if(data.A[c].find(a) != data.A[c].end() &&
							data.H[c][l].find(h1) != data.H[c][l].end()){
							fo += f(0,c,a,data.hospitals[h],l,h1,
								data,factor)*y0_ahclh[a][h][c][l][h1];
						}
					}
				}
			}
		}
	}

	factor = 0;
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

		// for(int a = 0; a < data.types_amb; ++a){
		// 	for(int h = 0; h < data.num_hosps; ++h){
		// 		for(int b = 0; b < data.num_bases; ++b){
		// 			fo += data.times[data.hospitals[h]][data.bases[b]]*
		// 				yt_ahb[t-1][a][h][b];
		// 		}
		// 	}
		// }
	}
	// factor = 125;
	for(int t = 0; t < data.num_times; ++t){
		for(int c = 0; c < data.types_call; ++c){
			factor = data.ins.penalties[c]*36000;
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


void AmbulanceModel::solve(){
	// if(data.debug){
		// model.write("amb_model.lp");
	// }
	model.set(GRB_IntParam_OutputFlag, 0);
	model.set(GRB_IntParam_Presolve, 0);
	model.set(GRB_IntParam_Threads, 1);
	// std::cout << "Number of variables: " << model.get(GRB_IntAttr_NumVars) << "\n";
	model.optimize();

	// GRBVar* vars = model.getVars();
	// for(int i = 0; i < model.get(GRB_IntAttr_NumVars); ++i){
	// 	if(vars[i].get(GRB_DoubleAttr_X) > g_params.EPS){
	// 		std::cout << vars[i].get(GRB_StringAttr_VarName) << " = ";
	// 		std::cout << vars[i].get(GRB_DoubleAttr_X) << "\n";
	// 	}
	// }
	// delete[] vars;
		
	std::cout << "STATUS " << model.get(GRB_IntAttr_Status) << "\n";
	std::cout << "Z (DET): " << model.get(GRB_DoubleAttr_ObjVal) << "\n";
	// std::cin.get();
	amb_type = -1;
	src_hosp = -1;
	base_return = -1;
	
	call_type = -1;
	call_location = -1;
	dst_hosp = -1;

	bool found = false;

	for(int a = 0; a < data.types_amb && !found; ++a){
		for(int h = 0; h < data.num_hosps && !found; ++h){
			for(int b = 0; b < data.num_bases && !found; ++b){
				if(y0_ahb[a][h][b].get(GRB_DoubleAttr_X) > g_params.EPS){
					amb_type = a;
					src_hosp = h;
					base_return = b;
					found = true;
					std::cout << y0_ahb[a][h][b].get(GRB_StringAttr_VarName) << "\n";
					break;
				}
			}
		}
	}

	for(int a = 0; a < data.types_amb && !found; ++a){
		for(int c = 0; c < data.types_call && !found; ++c){
			for(int h1 = 0; h1 < data.num_hosps  && !found; ++h1){
				for(int l = 0; l < data.num_locals && !found; ++l){
					for(int h = 0; h < data.num_hosps && !found; ++h){
						if(data.A[c].find(a) != data.A[c].end() &&
							data.H[c][l].find(h) != data.H[c][l].end() &&
							y0_ahclh[a][h1][c][l][h].get(GRB_DoubleAttr_X) > g_params.EPS){
							amb_type = a;
							call_type = c;
							call_location = l;
							src_hosp = h1;
							dst_hosp = h;
							found = true;
							std::cout << y0_ahclh[a][h1][c][l][h].get(GRB_StringAttr_VarName);
							std::cout << "\n";
							break;
						}
					}
				}
			}
		}
	}
	// GRBVar* vars = model.getVars();
	// std::cout << "OPTIMIZED! SOLUTION:\n";
	// for(int i = 0; i < model.get(GRB_IntAttr_NumVars); ++i){
	// 	if(vars[i].get(GRB_DoubleAttr_X) > g_params.EPS){
	// 		std::cout << vars[i].get(GRB_StringAttr_VarName) << " = ";
	// 		std::cout << vars[i].get(GRB_DoubleAttr_X) << "\n";
	// 	}
	// }
	// delete[] vars;
}

void AmbulanceModel::add_bases_constraints_t0(){
	std::stringstream name;
	GRBLinExpr con = 0;

	//0.5 - Flow at bases, t = 0
	for(int a = 0; a < data.types_amb; ++a){
		for(int b = 0; b < data.num_bases; ++b){
			con += A0_ab[a][b];
			for(int h = 0; h < data.num_hosps; ++h){
				if(data.L(0,a,data.hospitals[h], data.bases[b]) == data.bases[b]){
					con += y0_ahb[a][h][b];
				}
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

void AmbulanceModel::add_hospitals_constraints_t0(){
	std::stringstream name;
	GRBLinExpr con = 0;
	
	// 0.6 - Flow at hospitals
	for(int a = 0; a < data.types_amb; ++a){
		for(int h = 0; h < data.num_hosps; ++h){
			for(int b = 0; b < data.num_bases; ++b){
				con += y0_ahb[a][h][b];
			}

			for(int c = 0; c < data.types_call; ++c){
				for(int l = 0; l < data.num_locals; ++l){
					for(int h1 = 0; h1 < data.num_hosps; ++h1){
						if(data.H[c][l].find(h1) != data.H[c][l].end()){
							if(data.A[c].find(a) !=data. A[c].end())
								con += y0_ahclh[a][h][c][l][h1];
						}
					}
				}
			}

			name << "flow_hospitals_a" << a << "_h" << h;
			model.addConstr(con, GRB_EQUAL, A0_ah[a][h], name.str());
			name.str("");
			con = 0;
		}
	}
}

void AmbulanceModel::add_locations_constraints_t0(){
	std::stringstream name;
	GRBLinExpr con = 0;
	// 0.7 - Fluxo at location between hospitals and bases
	for(int a = 0; a < data.types_amb; ++a){
		for(int b = 0; b < data.num_bases; ++b){
			for(int l1 = 0; l1 < data.total_locals; ++l1){
				if(data.L_tab[0][a][b].find(l1) != data.L_tab[0][a][b].end()){
					for(int h = 0; h < data.num_hosps; ++h){
						if(data.L(0,a, data.hospitals[h],data.bases[b]) == l1){
							con += y0_ahb[a][h][b];
						}
					}

					for(int l2 = 0; l2 < data.total_locals; ++l2){
						if(data.L(0,a,l2, data.bases[b]) == l1 &&
							data.L_tab[0][a][b].find(l2) != 
							data.L_tab[0][a][b].end() ){
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
}


void AmbulanceModel::add_queues_constraints_t0(){
	std::stringstream name;
	GRBLinExpr con = 0;

	// 0.8 - Flow at queues
	for(int c = 0; c < data.types_call; ++c){
		for(int l = 0; l < data.num_locals; ++l){
			con += C[c][l];

			for(int h = 0; h < data.num_hosps; ++h){
				if(data.H[c][l].find(h) != data.H[c][l].end()){
					for(int a = 0; a < data.types_amb; ++a){
						for(int h1 = 0; h1 < data.num_hosps; ++h1){
							if(data.A[c].find(a) != data.A[c].end()){
								con -= y0_ahclh[a][h1][c][l][h];
							}
						}
					}
				}
			}

			name << "flow_queues_c" << c << "_l" << l;
			model.addConstr(con, GRB_EQUAL, Ct_cl[0][c][l],name.str());
			name.str("");
			con = 0;
		}
	}
}

void AmbulanceModel::add_bases_constraints(){
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
				for(int l1 = 0; l1 < data.total_locals; ++l1){
					if(data.L_tab[t+1][a][b].find(l1) != 
						data.L_tab[t+1][a][b].end() && 
						data.L(t,a,l1,data.bases[b]) == data.bases[b]){
						con += At_alb[t][a][l1][b];
						for(int c = 0; c < data.types_call; ++c){
							for(int l = 0; l < data.num_locals; ++l){
								for(int h = 0; h < data.num_hosps; ++h){
									if(data.H[c][l].find(h) !=
										data.H[c][l].end() &&
										data.A[c].find(a) != data.A[c].end()){
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
						for(int h = 0; h < data.num_hosps; ++h){
							if(data.H[c][l].find(h) !=
								data.H[c][l].end() &&
								data.A[c].find(a) != data.A[c].end()){
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


void AmbulanceModel::add_hospitals_constraints(){
	std::stringstream name;
	GRBLinExpr con = 0;

		// 0.10 - Flow at hospitals, t > 0
	for(int t = 0; t < data.num_times-1; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			for(int h = 0; h < data.num_hosps; ++h){

				for(int c = 0; c < data.types_call; ++c){
					for(int l = 0; l < data.num_locals; ++l){
						for(int h1 = 0; h1 < data.num_hosps; ++h1){
							if(data.H[c][l].find(h1) !=
								data.H[c][l].end() &&
								data.A[c].find(a) != data.A[c].end()){
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

				for(int c = 0; c < data.types_call; ++c){
					for(int l = 0; l < data.num_locals; ++l){
						for(int h1 = 0; h1 < data.num_hosps; ++h1){
							if((t)*data.quantum <= 
								data.tao[0][c][a][data.hospitals[h1]][l][h]  &&
								data.tao[0][c][a][data.hospitals[h1]][l][h] <= 
								(t+1)*data.quantum &&
								data.H[c][l].find(h) != data.H[c][l].end() &&
								data.A[c].find(a) != data.A[c].end()){
								con += y0_ahclh[a][h1][c][l][h];
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

void AmbulanceModel::add_locations_constraints(){
	std::stringstream name;
	GRBLinExpr con = 0;

	// 0.11 - Flow at locations between hospitals and bases,  t > 0
	for(int t = 0; t < data.num_times-1; ++t){
		for(int a  = 0; a < data.types_amb; ++a){
			for(int b = 0; b < data.num_bases; ++b){
				for(int l1 = 0; l1 < data.total_locals; ++l1){
					if(data.L_tab[t+1][a][b].find(l1) !=
						 data.L_tab[t+1][a][b].end()){
						for(int h1 = 0; h1 < data.num_hosps; ++h1){
							if(data.L(t,a, data.hospitals[h1], data.bases[b]) == l1){
								con += yt_ahb[t][a][h1][b];
							}
						}

						for(int l2 = 0; l2 < data.total_locals; ++l2){
							if(data.L_tab[t+1][a][b].find(l2) != 
								data.L_tab[t+1][a][b].end() &&
								data.L(t,a,l2, data.bases[b]) == l1){
								con += At_alb[t][a][l2][b];
								for(int c = 0; c < data.types_call; ++c){
									for(int l = 0; l < data.num_locals; ++l){
										for(int h = 0; h < data.num_hosps; ++h){
											if(data.H[c][l].find(h) !=
												data.H[c][l].end() &&
												data.A[c].find(a) != 
												data.A[c].end()){
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
}

void AmbulanceModel::add_queues_constraints(){
	std::stringstream name;
	GRBLinExpr con = 0;

	// 0.12 - Flow at queues, t > 0
	for(int t = 0; t < data.num_times-1; ++t){
		for(int c = 0; c < data.types_call; ++c){
			for(int l = 0; l < data.num_locals; ++l){
				con += Ct_cl[t][c][l];
				con += data.lambda[t][c][l];
				for(int a = 0; a < data.types_amb; ++a){
					if(data.A[c].find(a) != data.A[c].end()){
						for(int b = 0; b < data.num_bases; ++b){
							for(int h = 0; h < data.num_hosps; ++h){
								if(data.H[c][l].find(h) != data.H[c][l].end() &&
									data.A[c].find(a) != data.A[c].end()){
									//If using xt0, put t+1 on 1st index
									con -= xt_cablh[t][c][a][b][l][h];
								}
							}
						}
					}
				}

				for(int a = 0; a < data.types_amb; ++a){
					if(data.A[c].find(a) != data.A[c].end()){
						for(int h1 = 0; h1 < data.num_hosps; ++h1){
							for(int h = 0; h < data.num_hosps; ++h){
								if(data.H[c][l].find(h) != data.H[c][l].end() &&
									data.A[c].find(a) != data.A[c].end()){
									con -= xt_cahlh[t][c][a][h1][l][h];
								}
							}
						}
					}
				}

				for(int a = 0; a < data.types_amb; ++a){
					if(data.A[c].find(a) != data.A[c].end()){
						for(int b = 0; b < data.num_bases; ++b){
							for(int l1 = 0; l1 < data.total_locals; ++l1){
								for(int h = 0; h < data.num_hosps; ++h){
									if(data.H[c][l].find(h) != data.H[c][l].end() &&
										data.A[c].find(a) != data.A[c].end() &&
										data.L_tab[t+1][a][b].find(l1) != 
										data.L_tab[t+1][a][b].end()){
										//If using xt0, put t+1 on 1st index
										con -= xt_calblh[t][c][a][l1][b][l][h];
									}
								}
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

void AmbulanceModel::add_ambs_location_constraints(){
	std::stringstream name;
	GRBLinExpr con = 0;

	//Location constraints
	for(int t = 0; t < data.num_times-1; ++t){
		for(int a = 0; a < data.types_amb; ++a){
			for(int b = 0; b < data.num_bases; ++b){
				for(int l1 = 0; l1 < data.total_locals; ++l1){
					for(int c = 0; c < data.types_call; ++c){
						for(int l = 0; l < data.num_locals; ++l){
							for(int h = 0; h < data.num_hosps; ++h){
								if(data.H[c][l].find(h) != data.H[c][l].end() &&
									data.A[c].find(a) != data.A[c].end()  &&
									data.L_tab[t+1][a][b].find(l1) != 
									data.L_tab[t+1][a][b].end()){
									con += xt_calblh[t][c][a][l1][b][l][h];
								}
							}
						}
					}

					if(data.L_tab[t][a][b].find(l1) != data.L_tab[t][a][b].end()){
						name << "dispatch_locals_t" << t+1 << "_a" << a << "_b";
						name << b << "_l" << l1;
						model.addConstr(con, GRB_LESS_EQUAL,At_alb[t][a][l1][b],
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

void AmbulanceModel::add_base_cap_constraints(){
	std::stringstream name;
	GRBLinExpr con = 0;
	// 0.13 - Base capacities
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


AmbulanceModel::~AmbulanceModel(){
	for(int a = 0; a < data.types_amb; ++a){
		for(int h  = 0; h < data.num_hosps; ++h){
			delete[] y0_ahb[a][h];
		}
		delete[] y0_ahb[a];
	}

	for(int a = 0; a < data.types_amb; ++a){
		for(int h1 = 0; h1 < data.num_hosps; ++h1){
			for(int c = 0; c < data.types_call; ++c){
				for(int l  = 0; l < data.num_locals; ++l){
					delete[] y0_ahclh[a][h1][c][l];
				}
				delete[] y0_ahclh[a][h1][c];
			}
			delete[] y0_ahclh[a][h1];
		}
		delete[] y0_ahclh[a];
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

	delete[] y0_ahb; delete[] xt_cablh; delete[] xt_cahlh; delete[] xt_calblh; 
	delete[] y0_ahclh; delete[] Ct_cl; delete[] At_ab; delete[] At_alb; 
	delete[] yt_ahb;
}


	// if(model.get(GRB_IntAttr_Status) == 3){
	// 	model.write("ambulance.lp");
	// 	for(auto& amb: sim.ambulances){
	// 		std::cout << amb;
	// 	}

	// 	std::cout << "============================\n";
	// 	std::cout << "FILA DEBUG:\n";
	// 	for(int c = 0; c < data.types_call; ++c){
	// 		for(int l = 0; l < data.num_locals; ++l){
	// 			std::cout << C[c][l] << " ";
	// 		}
	// 		std::cout << "\n";
	// 	}
	// 	std::cout << "============================\n";
	// 	for(int c = 0; c < data.types_call; ++c){
	// 		for(int a = 0; a < data.types_amb; ++a){
	// 			for(int l1 = 0; l1 < data.total_locals; ++l1){
	// 				for(int h = 0; h < data.num_hosps; ++h){
	// 					if(A0_calh[c][a][l1][h] > 0){
	// 						std::cout << "A0_c" << c << "_a" << a;
	// 						std::cout << "_l'" << l1 << "_h" << h << " ";
	// 						std::cout << A0_calh[c][a][l1][h] << "\n"; 
	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	// 	std::cout << "============================\n";
	// 	for(int c = 0; c < data.types_call; ++c){
	// 		for(int a = 0; a < data.types_amb; ++a){
	// 			for(int l1 = 0; l1 < data.total_locals; ++l1){
	// 				for(int l = 0; l < data.num_locals; ++l){
	// 					for(int h = 0; h < data.num_hosps; ++h){
	// 						if(A0_callh[c][a][l1][l][h] > 0){
	// 							std::cout << "A0_c" << c << "_a" << a;
	// 							std::cout << "_l'" << l1 << "_l" << l;
	// 							std::cout << "_h" << h << " ";
	// 							std::cout <<A0_callh[c][a][l1][l][h]<<"\n";
	// 						}
	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	// 	std::cout << "============================\n";
	// 	for(int a = 0; a < data.types_amb; ++a){
	// 		for(int b = 0; b < data.num_bases; ++b){
	// 			if(A0_ab[a][b] > 0){
	// 				std::cout << "A0_a" << a << "_b" << b << " ";
	// 				std::cout << A0_ab[a][b] << "\n";
	// 			}
	// 		}
	// 	}
	// 	std::cout << "============================\n";
	// 	for(int a = 0; a < data.types_amb; ++a){
	// 		for(int h = 0; h < data.num_hosps; ++h){
	// 			if(A0_ah[a][h] > 0){
	// 				std::cout << "A0_a" << a << "_h" << h << " ";
	// 				std::cout << A0_ah[a][h] << "\n";
	// 			}
	// 		}
	// 	}
	// 	std::cout << "============================\n";
	// 	for(int a = 0; a < data.types_amb; ++a){
	// 		for(int l1 = 0; l1 < data.total_locals; ++l1){
	// 			for(int b = 0; b < data.num_bases; ++b){
	// 				if(A0_alb[a][l1][b] > 0){
	// 					std::cout << "A0_a" << a << "_l'" << l1;
	// 					std::cout << "_b" << b << " ";
	// 					std::cout << A0_alb[a][l1][b] << "\n";
	// 				}
	// 			}
				
	// 		}
	// 	}
	// 	model.computeIIS();
	// 	model.write("amb_t2.ilp");
	// 	std::cin.get();
	// 	std::cout << "";
	// }