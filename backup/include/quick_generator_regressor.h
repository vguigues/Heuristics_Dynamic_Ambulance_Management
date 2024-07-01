#ifndef _QUICK_GENERATOR_REGRESSOR_H
#define _QUICK_GENERATOR_REGRESSOR_H

#include "main.h"


class QuickGeneratorRegressor{
public:
	GRBEnv& env;
	double x_max, y_max;
	int n_x, n_y;
	int nb_weeks, nb_years;
	double slot_duration;
	std::vector<double> duration;
	ulong T;
	ulong G;
	std::vector<int> daily_obs;
	ulong R;
	ulong P;
	ulong nb_holidays_years;
	ulong nb_land_types;
	std::vector<std::pair<bool, int>> is_holidays;
	std::vector<Location> regions;

	double alpha, sigma, beta, epsilon;
	int max_iter;

	xt::xarray<int> nb_observations;

	xt::xarray<int> nb_observations_holidays;
	xt::xarray<int> nb_calls_holidays;
	xt::xarray<int> nb_calls_no_holidays;

	xt::xarray<double> regressors;

	std::vector<std::vector<int>> neighbors;
	std::vector<int> type;
	xt::xarray<double> distance;
	xt::xarray<int> sample_calls;
	xt::xarray<double> beta_teorico;
	xt::xarray<double> delta_teorico;

	QuickGeneratorRegressor(GRBEnv& env, std::string calls_path = "calls.dat", 
		std::string neighbors_path = "neighbors.dat", std::string info_path = "info.dat");
	~QuickGeneratorRegressor() = default;

	std::pair<xt::xarray<double>,xt::xarray<double>>
	oracle_gradient_model2(xt::xarray<double>& x_beta, 
		xt::xarray<double>& x_delta, double alpha, int t, int p);
	double oracle_objective_model2(xt::xarray<double>& x_beta, 
		xt::xarray<double>& x_delta, double alpha, int t, int p);
	std::vector<double> projected_gradient_armijo_boundary(xt::xarray<double>& x_beta,
		xt::xarray<double>& x_delta, double alpha, double sigma, double beta_bar, int t, 
		int p);
	std::vector<double> projected_gradient_armijo_feasible(xt::xarray<double>& x_beta,
		xt::xarray<double>& x_delta, double alpha, double sigma, double beta_tilde, 
		double beta_hat, int t, int p);

	std::pair<xt::xarray<double>,xt::xarray<double>> 
	projection_regressors(xt::xarray<double>& x_beta, xt::xarray<double>& x_delta, int t, int p);


	double cross_validation(xt::xarray<double>& x_beta, xt::xarray<double>& x_delta, 
		double sigma, double beta_bar, double proportion);
	void test();
	bool is_neighbor(int r, int s);

	double average_difference(xt::xarray<double>& x_beta, xt::xarray<double>& x_delta);
	double max_difference(xt::xarray<double>& x_beta, xt::xarray<double>& x_delta);

	void comp_wise_max(xt::xarray<double>& z ,xt::xarray<double>& a);

	void write_params(xt::xarray<double>& x_beta, xt::xarray<double>& x_delta, double alpha);
};


#endif