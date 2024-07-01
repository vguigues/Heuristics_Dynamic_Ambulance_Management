#ifndef _GENERATOR_H
#define _GENERATOR_H

#include "main.h"


class Generator{
public:
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

	double alpha, sigma, beta, epsilon;
	int max_iter;

	std::vector<Location> regions;

	xt::xarray<int> nb_observations;
	xt::xarray<int> nb_calls;

	std::vector<std::vector<int>> neighbors;
	std::vector<int> type;
	xt::xarray<double> distance;
	xt::xarray<int> sample_calls; 
	xt::xarray<double> lambda_teorico;
	xt::xarray<double> lambda_teorico_agg;
	xt::xarray<double> calls_agg;
	ulong nb_holidays_years;
	ulong nb_land_types;
	// std::vector<std::pair<bool, int>> is_holidays;
	xt::xarray<double> regressors;

	Generator(std::string calls_path = "calls.dat", 
		std::string neighbors_path = "neighbors.dat", std::string info_path = "info.dat");
	~Generator() = default;

	xt::xarray<double> oracle_gradient_model1(xt::xarray<double>& lambda, double alpha);
	double oracle_objective_model1(xt::xarray<double>& lambda, double alpha);
	std::vector<double> projected_gradient_armijo_boundary(xt::xarray<double>& x, 
		double alpha, double sigma, double beta_bar);
	std::vector<double> projected_gradient_armijo_feasible(xt::xarray<double>& x, 
		double alpha, double sigma);


	double cross_validation(xt::xarray<double>& x, double sigma, double beta_bar, 
		double proportion);
	void test();
	bool is_neighbor(int r, int s);

	double average_difference(xt::xarray<double>& x);
	double max_difference(xt::xarray<double>& x);

	void comp_wise_max(xt::xarray<double>& z ,xt::xarray<double>& a, double eps);

	void write_params(xt::xarray<double>& x, double alpha);	
};


#endif