#include "../include/generator_regressor.h"

using namespace std;

//Hardcoded test
GeneratorRegressor::GeneratorRegressor(GRBEnv& env): env(env){
	x_max = y_max = 10;
	n_x = n_y = 10;
	R = n_x*n_y;
	C = 1;
	T = 4;
	D = 7;

	which_group = vector<vector<int>>(D, 
		vector<int>(T, 0));

	for(int d = 0; d < 7; ++d){
		which_group[d][0] = 0;
		which_group[d][1] = 1;
		which_group[d][2] = 0;
		which_group[d][3] = 1;
	}

	// for(int d = 7; d < 15; ++d){
	// 	which_group[d][0] = 2;
	// 	which_group[d][1] = 3;
	// 	which_group[d][2] = 2;
	// 	which_group[d][3] = 3;
	// }

	vector<pair<int,int>> aux_group;
	for(int i = 0; i < 7; ++i){
		aux_group.push_back(make_pair(i,0));
	}
	groups.push_back(aux_group);
	aux_group.clear();
	for(int i = 0; i < 7; ++i){
		aux_group.push_back(make_pair(i,1));
	}
	groups.push_back(aux_group);
	aux_group.clear();
	for(int i = 0; i < 7; ++i){
		aux_group.push_back(make_pair(i,2));
	}
	groups[0].insert(groups[0].end(), aux_group.begin(), 
		aux_group.end());
	aux_group.clear();
	for(int i = 0; i < 7; ++i){
		aux_group.push_back(make_pair(i,3));
	}
	groups[1].insert(groups[1].end(), aux_group.begin(), 
		aux_group.end());
	aux_group.clear();
	// for(int i = 7; i < 15; ++i){
	// 	aux_group.push_back(make_pair(i,0));
	// }
	// groups.push_back(aux_group); //groups[2]
	// aux_group.clear();
	// for(int i = 7; i < 15; ++i){
	// 	aux_group.push_back(make_pair(i,1));
	// }
	// groups.push_back(aux_group); //groups[3]
	// aux_group.clear();
	// for(int i = 7; i < 15; ++i){
	// 	aux_group.push_back(make_pair(i,2));
	// }
	// groups[2].insert(groups[2].end(), aux_group.begin(), 
	// 	aux_group.end());
	// aux_group.clear();
	// for(int i = 7; i < 15; ++i){
	// 	aux_group.push_back(make_pair(i,3));
	// }
	// groups[3].insert(groups[3].end(), aux_group.begin(), 
	// 	aux_group.end());
	// aux_group.clear();

	fmt::print("Groups\n");
	for(int i = 0; i < groups.size(); ++i){
		fmt::print("{}: {}\n", i, groups[i]);
	}
	fmt::print("Which Group:\n");
	for(int d = 0; d < D; ++d){
		fmt::print("{}\n", which_group[d]);
	}
	
	nb_weeks = 52*100;
	nb_years = floor(nb_weeks/52);
	durations = vector<double>(T,6);
	nb_holidays_years = 8;
	// is_holidays = vector<pair<bool, int>>(nb_weeks*7, make_pair(false,-1));
	// vector<int> days_h;
	// for(int i = 0; i < 7; ++i){
	// 	days_h.push_back(i+1);
	// }
	// vector<int> index_h;
	// for(int i = 7; i < 15; ++i){
	// 	index_h.push_back(i+1);
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	for(int k = 0; k < days_h.size(); ++k){
	// 		is_holidays[year*52*7+days_h[k]] = make_pair(true,k+1);
	// 	}
	// }

	// for(int k = 0; k < days_h.size(); ++k){
	// 	if(nb_years*52*7 + days_h[k] <= nb_weeks*7){
	// 		is_holidays[nb_years*52*7+days_h[k]-1] = make_pair(true,k+1);
	// 	}
	// }

	nb_land_types = 2;
	nb_regressors = 1 + nb_land_types;

	beta_teorico = xt::zeros<double>({C,D,T,nb_regressors});
	regressors = xt::zeros<double>({nb_regressors, R});

	for(int d = 0; d < 7; ++d){
		beta_teorico(0,d,1,0) = 0.05;
		beta_teorico(0,d,3,0) = 0.05;

		beta_teorico(0,d,0,1) = 6;
		beta_teorico(0,d,1,1) = 18;
		beta_teorico(0,d,2,1) = 6;
		beta_teorico(0,d,3,1) = 18;
		
		beta_teorico(0,d,0,2) = 3;
		beta_teorico(0,d,1,2) = 6;
		beta_teorico(0,d,2,2) = 3;
		beta_teorico(0,d,3,2) = 6;
	}

	// for(int d = 7; d < 15; ++d){
	// 	beta_teorico(0,d,1,0) = 0.1;
	// 	beta_teorico(0,d,3,0) = 0.1;
		
	// 	beta_teorico(0,d,0,1) = 12;
	// 	beta_teorico(0,d,1,1) = 36;
	// 	beta_teorico(0,d,2,1) = 12;
	// 	beta_teorico(0,d,3,1) = 36;
		
	// 	beta_teorico(0,d,0,2) = 6;
	// 	beta_teorico(0,d,1,2) = 12;
	// 	beta_teorico(0,d,2,2) = 6;
	// 	beta_teorico(0,d,3,2) = 12;
	// }

	std::default_random_engine gen(600);
	std::uniform_real_distribution<double> rnd(0,1);
	for(int r = 0; r < R; ++r){
		int i = r / n_x;
		int j = r % n_x;
		if(i < 5 && j < 5){
			regressors(0, r) = 50 + 50*rnd(gen);
			regressors(1, r) = 0.5;
			regressors(2, r) = 0.25;
		}else if(i < 5 && j >= 5){
			regressors(0, r) = 50*rnd(gen);
			regressors(1, r) = 0.25;
			regressors(2, r) = 0.5;
		}else if(i >= 5 && j < 5){
			regressors(0, r) = 50*rnd(gen);
			regressors(1, r) = 0.25;
			regressors(2, r) = 0.5;
		}else{
			regressors(0, r) = 50 + 50*rnd(gen);
			regressors(1, r) = 0.5;
			regressors(2, r) = 0.25;
		}
	}

	sample_calls = xt::zeros<int>({C,D,T,R,static_cast<ulong>(7*nb_weeks)});
	nb_observations = xt::zeros<int>({C,D,T,R});
	nb_calls = xt::zeros<int>({C,D,T,R});

	int max_obs = 0;
	for(int index = 0; index < nb_weeks*7; ++index){ //each day in sample space
		int day = (index+1) % 7;
		if(day == 0){
			day = 6;
		}
		// if(is_holidays[index].first){
		// 	day = 6 + is_holidays[index].second;
		// }
		for(int c = 0; c < C; ++c){ //1
			for(int t = 0; t < T; ++t){ //4
				for(int r = 0; r < R; ++r){ // 100
					double rate = 0;
					for(int j = 0; j < nb_regressors; ++j){
						rate += beta_teorico(c,day,t,j)*regressors(j,r);
					}
					poisson_distribution<int> pd(rate);
					int this_nb_call = pd(gen);
					// fmt::print("Sample {} {} {} {} {}: {}\n", c,day,t,r, nb_observations(c,day,t,r),
					// 	this_nb_call);
					sample_calls(c,day,t,r,nb_observations(c,day,t,r)++) = this_nb_call;
					if(nb_observations(c,day,t,r) > max_obs){
						max_obs = nb_observations(c,day,t,r);
					}
					nb_calls(c,day,t,r) += this_nb_call;
				}
			}
		}
	}
	l_bounds = xt::zeros<double>(beta_teorico.shape());
	g_params.EPS = pow(10,-6);
	sigma = 0.5;
	beta_bar = 1;
	max_iter = 100;
	weight = 0;
	fmt::print("max obs = {}\n", max_obs);
	std::cout << "Initialized\n";
}

GeneratorRegressor::GeneratorRegressor(GRBEnv& env, std::string calls_path, 
		std::string neighbors_path, std::string info_path): env(env){
	// auto info_arq = ifstream(info_path, ios::in);
	// info_arq >> T >> G >> R >> P >> nb_land_types >> nb_holidays_years; 
	// slot_duration = 24 / T;
	// daily_obs = std::vector<int>(G, 0);
	// for(int g = 0; g < G; ++g){
	// 	info_arq >> daily_obs[g];
	// }
	// info_arq.close();

	// int max_obs = *max_element(daily_obs.begin(), daily_obs.end());

	// xt::xarray<double>::shape_type shape_d = {T,G,R,P};
	// xt::xarray<int>::shape_type shape_i = {T,G,R,P};
	// nb_observations = xt::zeros<int>(shape_i);

	// shape_i = {T,G,R,P,static_cast<ulong>(max_obs)};
	// sample_calls = xt::zeros<int>(shape_i);

	// shape_i = {T, nb_holidays_years, R, P};
	// nb_observations_holidays = xt::zeros<int>(shape_i);
	// shape_i = {T,nb_holidays_years,G,R,P};
	// nb_calls_holidays = xt::zeros<int>(shape_i);
	// shape_i = {T,G,R,P};
	// nb_calls_no_holidays = xt::zeros<int>(shape_i);

	// auto calls_arq = ifstream(calls_path, ios::in);
	// std::string aux_str;
	// is_holidays = std::vector<std::pair<bool, int>>(max_obs, make_pair(false,-1));
	// do{
	// 	std::getline(calls_arq, aux_str);
	// 	if(aux_str == "END"){
	// 		break;
	// 	}
	// 	std::istringstream ss(aux_str);
	// 	int t, g, r, p, j, h, val;
	// 	ss >> t >> g >> r >> p >> j >> val >> h;
	// 	sample_calls(t,g,r,p,j) = val;
	// 	if(h == -1){
	// 		nb_calls_no_holidays(t,g,r,p) += sample_calls(t,g,r,p,j);
	// 	}else{
	// 		is_holidays[j] = make_pair(true, h);
	// 		nb_observations_holidays(t,h,r,p) += 1;
	// 		nb_calls_holidays(t,h,g,r,p) += sample_calls(t,g,r,p,j);
	// 	}
	// 	nb_observations(t,g,r,p) += 1;
	// }while(true);

	// calls_arq.close();


	// neighbors = std::vector<vector<int>>(R, std::vector<int>());
	// xt::xarray<double>::shape_type dist_shape = {R, R};
	// regions = std::vector<Location>(R, null_location);
	// distance = xt::xarray<double>(dist_shape);
	// type = std::vector<int>(R,-1);
	// xt::xarray<double>::shape_type reg_shape = {R, nb_land_types};
	// regressors = xt::zeros<double>(reg_shape);
	// auto neighbors_arq = ifstream(neighbors_path, ios::in);
	// while(true){
	// 	int ind, terrain_type, s; 
	// 	double lat, longi, dist;
	// 	std::getline(neighbors_arq, aux_str);
	// 	if(aux_str == "END"){
	// 		break;
	// 	}
	// 	std::istringstream ss(aux_str);
	// 	ss >> ind >> lat >> longi >> terrain_type;
	// 	type[ind] = terrain_type;
	// 	regions[ind] = make_pair(lat, longi);
	// 	for(int j = 0; j < nb_land_types; ++j){
	// 		ss >> regressors(ind,j);
	// 	}
	// 	while(ss >> s >> dist){
	// 		distance(ind,s) = dist;
	// 		neighbors[ind].push_back(s);
	// 	}
	// }
	// neighbors_arq.close();
	std::cout << "Initialized\n";
}

double GeneratorRegressor::cross_validation(xt::xarray<double>& x_beta, 
	xt::xarray<double>& x_delta, double sigma, double beta_bar, 
		double proportion){
	double best_alpha = -1;
	// std::vector<double> alphas{0.001,0.01,0.05,0.1,0.5,1,2,5,10,50,100,1000};
	// // std::vector<double> alphas{0.001,0.5,1};
	// double max_likelihood = GRB_INFINITY;
	// int nb_obs = static_cast<int>(sample_calls.shape(4));
	// int nb_in_block = nb_obs*proportion;
	// xt::xarray<double>::shape_type shape = {T,G,R,P};

	// double beta_tilde = 1;
	// double beta_hat = 2;

	// for(double alpha: alphas){
	// 	double likelihood = 0;
	// 	for(int ind_cross = 0; ind_cross < floor(1/proportion); ++ind_cross){
	// 		nb_observations = nb_in_block*xt::ones<int>(shape);
	// 		nb_observations_holidays = xt::zeros<int>(nb_observations_holidays.shape());
	// 		nb_calls_holidays = xt::zeros<int>(nb_calls_holidays.shape());
	// 		nb_calls_no_holidays = xt::zeros<int>(nb_calls_no_holidays.shape());

	// 		for(int ind = ind_cross*nb_in_block; ind < (ind_cross+1)*nb_in_block; ++ind){
	// 			for(int t = 0; t < T; ++t){
	// 				for(int g = 0; g < G; ++g){
	// 					for(int r = 0; r < R; ++r){
	// 						for(int p = 0; p < P; ++p){
	// 							if(is_holidays[ind].first){
	// 								int k = is_holidays[ind].second;
	// 								nb_observations_holidays(t,k,r,p) += 1;
	// 								nb_calls_holidays(t,k,g,r,p) += sample_calls(t,g,r,p,ind);
	// 							}else{
	// 								nb_calls_no_holidays(t,g,r,p) += sample_calls(t,g,r,p,ind);
	// 							}
	// 						}
	// 					}
	// 				}
	// 			}
	// 		}

	// 		auto f_val = projected_gradient_armijo_feasible(x_beta, x_delta, alpha, 
	// 			sigma, beta_tilde, beta_hat);
	// 		xt::xarray<int> nb_observations_remaining = (nb_obs-nb_in_block)*xt::ones<int>(shape);
	// 		xt::xarray<int> nb_observations_holidays_remaining = xt::zeros<int>(
	// 			nb_observations_holidays.shape());
	// 		xt::xarray<int> nb_calls_holidays_remaining = xt::zeros<int>(nb_calls_holidays.shape());
	// 		xt::xarray<int> nb_calls_no_holidays_remaining = xt::zeros<int>(
	// 			nb_calls_no_holidays.shape());

	// 		for(int ind = 0; ind < ind_cross*nb_in_block; ++ind){
	// 			for(int t = 0; t < T; ++t){
	// 				for(int g = 0; g < G; ++g){
	// 					for(int r = 0; r < R; ++r){
	// 						for(int p = 0; p < P; ++p){
	// 							if(is_holidays[ind].first){
	// 								int k = is_holidays[ind].second;
	// 								nb_observations_holidays_remaining(t,k,r,p) += 1;
	// 								nb_calls_holidays_remaining(t,k,g,r,p) += sample_calls(t,g,r,p,ind);
	// 							}else{
	// 								nb_calls_no_holidays_remaining(t,g,r,p) += sample_calls(t,g,r,p,ind);
	// 							}
	// 						}
	// 					}
	// 				}
	// 			}
	// 		}

	// 		for(int ind = (ind_cross+1)*nb_in_block; ind < nb_obs; ++ind){
	// 			for(int t = 0; t < T; ++t){
	// 				for(int g = 0; g < G; ++g){
	// 					for(int r = 0; r < R; ++r){
	// 						for(int p = 0; p < P; ++p){
	// 							if(is_holidays[ind].first){
	// 								int k = is_holidays[ind].second;
	// 								nb_observations_holidays_remaining(t,k,r,p) += 1;
	// 								nb_calls_holidays_remaining(t,k,g,r,p) += sample_calls(t,g,r,p,ind);
	// 							}else{
	// 								nb_calls_no_holidays_remaining(t,g,r,p) += sample_calls(t,g,r,p,ind);
	// 							}
	// 						}
	// 					}
	// 				}
	// 			}
	// 		}

	// 		double f = 0;
	// 		xt::xarray<double> rates = xt::zeros<double>(shape);
	// 		for(int t = 0; t < T; ++t){
	// 			for(int g = 0; g < G; ++g){
	// 				for(int r = 0; r < R; ++r){
	// 					for(int p = 0; p < P; ++p){
	// 						for(int j = 0; j < nb_land_types; ++j){
	// 							rates(t,g,r,p) += x_beta(t,g,p,j)*regressors(r,j);
	// 						}
	// 						f += nb_observations_remaining(t,g,r,p)*rates(t,g,r,p) - 
	// 							nb_calls_no_holidays_remaining(t,g,r,p)*log(rates(t,g,r,p));
	// 					}
	// 				}
	// 			}
	// 		}
	// 		for(int t = 0; t < T; ++t){
	// 			for(int k = 0; k < nb_holidays_years; ++k){
	// 				for(int r = 0; r < R; ++r){
	// 					for(int p = 0; p < P; ++p){
	// 						f += nb_observations_holidays_remaining(t,k,r,p)*x_delta(t,k,r,p);
	// 						for(int g = 0; g < G; ++g){
	// 							f -= nb_calls_holidays_remaining(t,k,g,r,p)*log(rates(t,g,r,p)*x_delta(t,k,r,p));
	// 						}
	// 					}
	// 				}
	// 			}
	// 		}

	// 		likelihood += f;
	// 	}
	// 	fmt::print("Alpha = {}, likelihood = {}\n", alpha, likelihood);
	// 	if(likelihood < max_likelihood){
	// 		max_likelihood = likelihood;
	// 		best_alpha = alpha;
	// 	}
	// }
	
	// nb_observations = nb_obs*xt::ones<int>(shape);
	// nb_observations_holidays = xt::zeros<int>(nb_observations_holidays.shape());
	// nb_calls_holidays = xt::zeros<int>(nb_calls_holidays.shape());
	// nb_calls_no_holidays = xt::zeros<int>(nb_calls_no_holidays.shape());
	// for(int ind = 0; ind < nb_obs; ++ind){
	// 	for(int t = 0; t < T; ++t){
	// 		for(int g = 0; g < G; ++g){
	// 			for(int r = 0; r < R; ++r){
	// 				for(int p = 0; p < P; ++p){
	// 					if(is_holidays[ind].first){
	// 						int k = is_holidays[ind].second;
	// 						nb_observations_holidays(t,k,r,p) += 1;
	// 						nb_calls_holidays(t,k,g,r,p) += sample_calls(t,g,r,p,ind);
	// 					}else{
	// 						nb_calls_no_holidays(t,g,r,p) += sample_calls(t,g,r,p,ind);
	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	// auto f_val_best = projected_gradient_armijo_feasible(x_beta, x_delta, best_alpha, 
	// 	sigma, beta_tilde, beta_hat);

	return best_alpha;
}


void GeneratorRegressor::test(){
	double epsilon = g_params.EPS;
	xt::xarray<double> x_beta = 2*epsilon*xt::ones<double>(beta_teorico.shape());
	// x_beta = beta_teorico;
	auto f_val2 = projected_gradient_armijo_feasible(x_beta);
	fmt::print("Mean diff = {}\n", average_difference(x_beta));
	// write_params(x_beta);

	// std::cout << "Avg diff Feasible = " << average_difference(x_beta, x_delta) << "\n";
	
	// double proportion = 0.2;

	// x_beta = epsilon*xt::ones<double>(x_beta.shape());
	// x_delta = epsilon*xt::ones<double>(x_delta.shape());
	// double best_alpha = cross_validation(x_beta, x_delta,sigma, beta_bar, proportion);
	// std::cout << "Best alpha: " << best_alpha << "\n";
	// std::cout << "Avg diff Cross = " << average_difference(x_beta, x_delta) << "\n";
}


std::vector<double> GeneratorRegressor::projected_gradient_armijo_boundary(
	xt::xarray<double>& x_beta, xt::xarray<double>& x_delta, double alpha, 
	double sigma, double beta_bar){

	std::vector<double> f_val;
	// xt::xarray<double> z_beta = xt::zeros<double>(x_beta.shape());
	// xt::xarray<double> z_delta = xt::zeros<double>(x_delta.shape());
	// int k = 0;
	// double eps = 0.1;
	// f_val.reserve(max_iter);
	// while(k < max_iter){
	// 	double fold = oracle_objective_model2(x_beta, x_delta, alpha);
	// 	xt::xarray<double> gradient_beta, gradient_delta;
	// 	tie(gradient_beta, gradient_delta) = oracle_gradient_model2(x_beta, x_delta, 
	// 		alpha);
	// 	bool stop = false;
	// 	int j = 0;
	// 	double f = 0;
	// 	do{
	// 		xt::xarray<double> x_beta_aux = x_beta - (beta_bar/pow(2,j))*gradient_beta;
	// 		xt::xarray<double> x_delta_aux = x_delta - (beta_bar/pow(2,j))*gradient_delta;
	// 		xt::xarray<double> z_beta, z_delta;
	// 		tie(z_beta, z_delta) = projection_regressors(x_beta_aux, x_delta_aux, eps);
	// 		f = oracle_objective_model2(z_beta, z_delta, alpha);
	// 		double rhs = fold - sigma*(xt::sum(gradient_beta*(x_beta-z_beta))() + 
	// 			xt::sum(gradient_delta*(x_delta-z_delta))());
	// 		fmt::print("k = {}, j = {}, f = {}, rhs = {}\n", k, j, f, rhs);
	// 		if(f <= rhs){
	// 			stop = true;
	// 		}else{
	// 			++j;
	// 		}
	// 	}while(!stop);
	// 	f_val.push_back(f);
	// 	x_beta = z_beta;
	// 	x_delta = z_delta;
	// 	++k;
	// }

	return f_val;
}

std::vector<double> GeneratorRegressor::projected_gradient_armijo_feasible(
	xt::xarray<double>& x_beta){
	using xt::linalg::dot;
	xt::xarray<double> z_beta = xt::zeros<double>(x_beta.shape());
	double eps = g_params.EPS;
	int k = 0;
	std::vector<double> f_val;
	double b_param = 2;
	double beta_k = b_param;
	f_val.reserve(max_iter);
	int j = 0;
	while(k < max_iter){
		double fold = oracle_objective_model2(x_beta);
		xt::xarray<double> gradient_beta = oracle_gradient_model2(x_beta);
		xt::xarray<double> x_beta_aux = x_beta-beta_k*gradient_beta;
		try{
			z_beta = projection_regressors(x_beta_aux);
		}catch(GRBException& ex){
			fmt::print("{} {}\n",ex.getErrorCode(), ex.getMessage());
			cin.get();
		}
		bool stop = false;
		j = 0;
		xt::xarray<double> diff_aux = x_beta-z_beta;
		double rhs = mat_prod(gradient_beta, diff_aux);
		fmt::print("fold = {} rhs = {}\n", fold, rhs);
		double f = GRB_INFINITY;
		xt::xarray<double> z_aux_beta = xt::zeros<double>(x_beta.shape());
		while(!stop){
			z_aux_beta = x_beta + (1/pow(2,j))*(z_beta-x_beta);
			f = oracle_objective_model2(z_aux_beta);
			fmt::print("\t f = {} test {}\n", f, fold -f -(sigma/pow(2,j))*rhs);
			if(f <= fold-(sigma/pow(2,j))*rhs){
				stop = true;
			}else{
				++j;
			}
		}
		f_val.push_back(f);
		x_beta = z_aux_beta;
		if(k % 1 == 0){
			fmt::print("k = {}, f = {}, j = {}\n", k, f, j);
		}
		++k;
		beta_k = b_param / pow(2,j);
		// cin.get();
	}
	return f_val;
}

double GeneratorRegressor::mat_prod(xt::xarray<double>& a, xt::xarray<double>& b){
	double sum = 0;
	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int j = 0; j < nb_regressors; ++j){
					sum += a(c,d,t,j)*b(c,d,t,j);
				}
			}
		}
	}
	return sum;
}



xt::xarray<double> GeneratorRegressor::oracle_gradient_model2(xt::xarray<double>& x_beta){
	
	xt::xarray<double> gradient_beta = xt::zeros<double>(x_beta.shape());

	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int r = 0; r < R; ++r){
					double rates = 0;
					for(int j = 0; j < nb_regressors; ++j){
						rates += x_beta(c,d,t,j)*regressors(j,r);
					}
					for(int j = 0; j < nb_regressors; ++j){
						gradient_beta(c,d,t,j) += nb_observations(c,d,t,r)*
							regressors(j,r)-nb_calls(c,d,t,r)*regressors(j,r) / 
							rates;
					}
				}
			}
		}
	}

	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(auto& e1: groups[which_group[d][t]]){
					int d1 = e1.first;
					int t1 = e1.second;
					if((d1 != d) || (t1 != t)){
						for(int j = 0; j < nb_regressors; ++j){
							gradient_beta(c,d,t,j) += (2*weight / durations[t]) *
								((gradient_beta(c,d,t,j) / 
								durations[t]) - (gradient_beta(c,d1,t1,j) / 
								durations[t1]));
						}
					}
				}
			}
		}
	}

	return gradient_beta;
}


double GeneratorRegressor::oracle_objective_model2(xt::xarray<double>& x_beta){
	double f = 0;
	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int r = 0; r < R; ++r){
					double rates = 0;
					for(int j = 0; j < nb_regressors; ++j){
						rates += x_beta(c,d,t,j)*regressors(j,r);
					}
					f += nb_observations(c,d,t,r)*rates - 
						nb_calls(c,d,t,r)*log(rates);
				}
			}
		}
	}

	for(auto& group: groups){
		for(auto& e1: group){
			int d1 = e1.first;
			int t1 = e1.second;
			for(auto& e2: group){
				if(e1 != e2){
					int d2 = e2.first;
					int t2 = e2.second;
					for(int c = 0; c < C; ++c){
						for(int j = 0; j < nb_regressors; ++j){
							f += (weight/2) * pow(
								(x_beta(c,d1,t1,j)/durations[t1]) - 
								(x_beta(c,d2,t2,j)/durations[t2]),
								2);
						}
					}
				}
			}
		}	
	}
	return f;
}

xt::xarray<double> GeneratorRegressor::projection_regressors(xt::xarray<double>& x_beta){

	xt::xarray<GRBVar> y_beta(x_beta.shape());

	GRBModel model(env);
	stringstream name;

	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int j = 0; j < nb_regressors; ++j){
					name << "yb_" << c << "_" << d << "_" << t << "_" << j;
					double ub = (j == 0) ? 1 : pow(10,6);
					y_beta(c,d,t,j) = model.addVar(l_bounds(c,d,t,j), ub,
						0,GRB_CONTINUOUS, name.str());
					name.str("");
				}
			}
		}
	}

	GRBQuadExpr obj = 0;
	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int j = 0; j < nb_regressors; ++j){
					obj += 0.5*y_beta(c,d,t,j)*y_beta(c,d,t,j) -
						x_beta(c,d,t,j)*y_beta(c,d,t,j);
				}
			}
		}
	}
	try{
		model.setObjective(obj, GRB_MINIMIZE);
	}catch(GRBException& ex){
		cout << ex.getMessage() << "\n";
	}

	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int r = 0; r < R; ++r){
					GRBLinExpr con1 = 0;
					for(int j = 0; j < nb_regressors; ++j){
						con1 += y_beta(c,d,t,j)*regressors(j,r);
					}
					name << "con1_" << c << "_" << d << "_" << t << "_" << r;
					model.addConstr(con1, GRB_GREATER_EQUAL, g_params.EPS, name.str());
					name.str("");
					con1 = 0;
				}
			}
		}
	}
	model.update();
	// model.write("test.lp");
	model.set(GRB_IntParam_OutputFlag,0);
	model.set(GRB_IntParam_NumericFocus, 3);
	model.set(GRB_IntParam_DualReductions, 0);

	model.optimize();

	auto status = model.get(GRB_IntAttr_Status);

	xt::xarray<double> beta_val(y_beta.shape());

	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int j = 0; j < nb_regressors; ++j){
					beta_val(c,d,t,j) = y_beta(c,d,t,j).get(GRB_DoubleAttr_X);
				}
			}
		}
	}
	return beta_val;

}


double GeneratorRegressor::average_difference(xt::xarray<double>& x_beta){
	double mean_rate_beta = 0;
	vector<double> rates;
	vector<double> rates_est;
	for(int c = 0; c < C; ++c){
		for(int d = 0; d < D; ++d){
			for(int t = 0; t < T; ++t){
				for(int r = 0; r < R; ++r){
					double rate = 0;
					double rate_est = 0;
					for(int j = 0; j < nb_regressors; ++j){
						rate += beta_teorico(c,d,t,j)*regressors(j,r);
						rate_est += x_beta(c,d,t,j)*regressors(j,r);
					}
					rates.push_back(rate);
					rates_est.push_back(rate_est);
					mean_rate_beta += 100*(abs(rate-rate_est)/rate);
				}
			}
		}
	}

	return mean_rate_beta/(D*T*R*C);
}



void GeneratorRegressor::comp_wise_max(xt::xarray<double>& z ,xt::xarray<double>& a, double eps){
	// for(int t = 0; t < T; ++t){
	// 	for(int g = 0; g < G; ++g){
	// 		for(int r = 0; r < R; ++r){
	// 			for(int p = 0; p < P; ++p){
	// 				z(t,g,r,p) = max(a(t,g,r,p), eps);
	// 			}
	// 		}
	// 	}
	// }
}


bool GeneratorRegressor::is_neighbor(int r, int s){
	return r != s;
}


void GeneratorRegressor::write_params(xt::xarray<double>& x_beta){
	// ofstream out_file(fmt::format("{}/xRegressorT{}G{}I{}P{}K{}J{}_alpha{}.txt",
	// 	g_params.generator_folder,T,G,R,P, nb_holidays_years, nb_land_types,alpha), ios::out);

	// for(int t = 0; t < T; ++t){
	// 	for(int k = 0; k < nb_holidays_years; ++k){
	// 		for(int r = 0; r < R; ++r){
	// 			for(int p = 0; p < P; ++p){
	// 				out_file << t << " " << k << " " << r << " " << p;
	// 				out_file << " " << x_delta(t,k,r,p) << "\n";
	// 			}
	// 		}
	// 	}
	// }

	// for(int t = 0; t < T; ++t){
	// 	for(int g = 0; g < G; ++g){
	// 		for(int p = 0; p < P; ++p){
	// 			for(int j = 0; j < nb_land_types; ++j){
	// 				out_file << t << " " << g << " " << p << " " << j;
	// 				out_file << " " << x_beta(t,g,p,j) << "\n";
	// 			}
	// 		}
	// 	}
	// }

	// out_file << "END";
	// out_file.close();

	// ofstream plot(fmt::format("{}/params_plot_reg_{}.txt", g_params.generator_folder, alpha), 
	// 	ios::out);
	// for(int k = 0; k < nb_holidays_years; ++k){
	// 	plot << k << " ";
	// 	for(int t = 0; t < T; ++t){
	// 		double sum = 0;
	// 		for(int r = 0; r < R; ++r){
	// 			for(int p = 0; p < P; ++p){
	// 				sum += x_delta(t,k,r,p);
	// 			}
	// 		}
	// 		plot << sum << " ";
	// 	}
	// 	plot << "\n";
	// }
	// for(int g = 0; g < G; ++g){
	// 	plot  << g << " ";
	// 	for(int t = 0; t < T; ++t){
	// 		double sum = 0;
	// 		for(int p = 0; p < P; ++p){
	// 			for(int j = 0; j < nb_land_types; ++j){
	// 				sum += x_beta(t,g,p,j);
	// 			}
	// 		}
	// 		plot << sum << " ";
	// 	}
	// 	plot << "\n";
	// }
	// plot.close();
}


	// for(int j = 0; j < nb_weeks*nb_years*G; ++j){
	// 	is_holidays.push_back(make_pair(false,0));
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+1) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[((year)*nb_weeks)*G + day] = make_pair(true,0);
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+4) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[(11 + (year)*nb_weeks)*G + day] = make_pair(true,1);
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+6) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[(23 + (year)*nb_weeks)*G + day] = make_pair(true,2);
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+6) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[(30 + (year)*nb_weeks)*G + day] = make_pair(true,3);
	// }


	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+7) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[(40 + (year)*nb_weeks)*G + day] = make_pair(true,4);
	// }

	// for(int year = 0; year < nb_years; ++year){
	// 	int day = (year+3) % 7;
	// 	if(day == 0){
	// 		day = 7;
	// 	}
	// 	is_holidays[(50 + (year)*nb_weeks)*G + day] = make_pair(true,5);
	// }
