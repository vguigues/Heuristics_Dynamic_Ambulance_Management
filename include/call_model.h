#ifndef _CALL_MODEL_H
#define _CALL_MODEL_H

#include "main.h"

using namespace std;

class Ambulance;
class Call;
class Instance;
class Simulator;
class Travel;
class Solver;

class CallModel{
public:
	// CallModel(GRBEnv& env, Simulator& sim);
	CallModel(GRBEnv& env, Solver& solver);
	~CallModel();

	Instance& ins;

	vector<Call>& calls;
	vector<int>& queue;
	vector<Ambulance>& ambulances;
	double time;
	Travel travel;

	GRBModel model;

	int n_calls;
	int n_ambulances;
	double M;

	GRBVar* M_p;
	GRBVar** x_ik;
	GRBVar*** x_ijk;
	GRBVar** z_ik;

	GRBVar** y_ic;
	GRBVar** y_ih;
	GRBVar*** y_ihc;
	GRBVar* t_i;

	std::vector<int> opt_which_amb;
	std::vector<double> arrival_times;

	void load_model();
	void add_constraints();
	void solve();
	void print_results();


	// int ambulances_at_hospital();
};

#endif