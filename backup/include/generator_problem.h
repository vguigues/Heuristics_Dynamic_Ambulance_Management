#ifndef _GENERATOR_PROBLEM_H
#define _GENERATOR_PROBLEM_H

#include "main.h"

class GeneratorProblem{
public:
    ulong T, R, C, D;
    double x_max, y_max;
    int n_x, n_y;
	int nb_weeks, nb_years;
    vector<double> durations;
    vector<pair<int, int>> groups;
    ulong nb_holidays_years;
    ulong nb_land_types;
    ulong nb_regressors;
    std::vector<std::pair<bool, int>> is_holidays;
    std::vector<Location> regions;
    
    xt::xarray<double> beta_teoricos;
    xt::xarray<int> sample_calls;
    xt::xarray<double> distance;
    std::vector<std::vector<int>> neighbors;
	std::vector<int> type;

    xt::xarray<int> nb_observations;
	xt::xarray<int> nb_observations_holidays;
	xt::xarray<int> nb_calls_holidays;
	xt::xarray<int> nb_calls_no_holidays;

	xt::xarray<double> regressors;


    GeneratorProblem();
    GeneratorProblem(std::string calls_path = "calls.dat", 
        std::string neighbors_path = "neighbors.dat", std::string info_path = "info.dat");

};


#endif
