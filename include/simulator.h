#ifndef _SIMULATOR_H
#define _SIMULATOR_H


#include "instance.h"
#include "osrm_helper.h"
#include "call_model.h"

class Instance;
class ClassModel;


struct Stats{
	double mean_waiting_on_scene;
	double mean_waiting_to_hospital;
	double max_waiting_on_scene;
	double max_waiting_to_hospital;
	double waiting_on_scene_q90;
	double waiting_to_hospital_q90;
	double mean_total;
	double max_total;
	double q90_total;

	Stats(vector<double> waiting_on_scene, vector<double> waiting_to_hospital);
};

class Simulator{
public:
	Simulator(Instance& ins, OSRMHelper& osrm);
	~Simulator();

	GRBEnv env;
	Instance& ins;
	OSRMHelper& osrm;
	bool forward,discard;
	double time;

	std::string solver;

	std::vector<Call> calls;
	std::vector<Ambulance> ambulances;

	std::vector<int> nearest_base;

	std::vector<double> calls_end;
	std::vector<double> waiting_on_scene;
	std::vector<double> waiting_to_hospital;
	std::vector<double> waiting_on_scene_penalized;
	std::vector<double> waiting_to_hospital_penalized;	

	std::vector<int> which_ambulance;
	std::vector<int> queue;

	void run();
	void print_results();

	void forward_heuristic();
	void minmax_heuristic();
	void queue_heuristic();
	void priority_heuristic();

	void gen_forward_heuristic();
	void gen_minmax_heuristic();


	void set_next_event(int& event_call, int& index_call);

	void set_calls_destinations();
	void set_calls_nearest_bases();

	void reset();
	void set_calls(std::vector<Call>& calls);
	void set_solver(const std::string& solver);

	double get_response_time(Ambulance& amb, Call& call, double current_time);
	double travel_time(Location& a, Location& b, double speed);
	Location ambulance_position(Ambulance& amb, double t);
	double travel_time_from_position(Ambulance& amb,  Call& call);
	double norm(Location& a, Location& b);
	
};

#endif