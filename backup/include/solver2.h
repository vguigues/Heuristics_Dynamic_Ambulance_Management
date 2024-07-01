#ifndef _SOLVER2_H
#define _SOLVER2_H

using namespace std;

#include "main.h"
#include "call.h"
#include "ambulance.h"
#include "fast_equipment.h"
#include "travel.h"


class Travel;

class Solver{
protected:
	Solver(GRBEnv& env, vector<Call>& calls, vector<FastEquipment>& equips, 
		vector<Ambulance>& ambulances, Instance& ins, Travel& travel);
public:
	virtual ~Solver();

	GRBEnv& env;
	vector<Call> calls;
	vector<FastEquipment> equips;
	vector<Ambulance> ambulances;
	Travel travel;
	Instance& ins;

	double time;
	double first_time;
	double last_time;

	vector<double> waiting_on_scene;
	vector<double> waiting_on_scene_penalized;
	vector<double> waiting_to_hospital;
	vector<double> calls_end;
	vector<int> which_ambulance;

	vector<int> queue;
	vector<int> nearest_base;
	double obj;

	std::vector<double> run_times;

	virtual void run() = 0;
	void set_next_event(int& event_call, int& index_call);
	void set_calls_nearest_bases();
	void print_results();
	bool can_answer(Ambulance& amb, Call& call);
};


class MinMaxSolver: public Solver{
public:
	MinMaxSolver(GRBEnv& env, vector<Call>& calls, vector<FastEquipment>& equips, 
		vector<Ambulance>& ambulances, Instance& ins, Travel& travel);

	virtual void run();	
};


#endif
