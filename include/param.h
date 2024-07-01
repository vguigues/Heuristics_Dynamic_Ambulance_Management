#ifndef _PARAM_H
#define _PARAM_H

#include "boost/program_options.hpp"
#include <fstream>
#include <iostream>

class Param{
public:
	Param();
	Param(const Param& p);
	explicit Param(boost::program_options::variables_map vm);

	std::string instance, solver, generator_folder, instance_type, amb_setup;
	bool closest_base, best_base, debug, h_use_fixed_bases, h_forward, h_discard, h_order_priorities,
		h_order_time, extended_model, h_random_hospital;
	int n_nearest_hospitals, n_nearest_ambs, n_queue_calls_eval,n_scenarios,
		n_time_slots, n_regions, n_hospitals, n_bases, n_cleaning_bases, n_ambulances;
	double EPS = 0.001;
	double min_preparedness = 0.4;
	std::string osrm_map_path;


	Param& operator=(const Param& p);
	friend std::ostream& operator<<(std::ostream& out, const Param &s);
	~Param();
};

#endif