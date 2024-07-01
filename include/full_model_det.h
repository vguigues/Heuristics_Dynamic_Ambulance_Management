#ifndef FULL_MODEL_SIM
#define FULL_MODEL_SIM

#include "main.h"
#include "data.h"
#include "ambulance.h"
#include "gurobi_c++.h"
#include "solver.h"
#include "travel.h"

class CGSolver;
class Travel;

// int get_location(Data& data, int time, Ambulance & amb);

class CallModel1{
public:
	CallModel1(Data & data, GRBEnv & env, Solver& solver);
	~CallModel1();

	GRBEnv& env;
	Data & data;
	GRBModel model;
	int type_call0;
	int local_call0;
	Solver& solver;


	int src_type;
	int src_location;
	int amb_type;
	int loc_base;
	int base_return;
	int src_hosp;
	int dst_hosp;


	int** A0_ab;
	int** A0_ah;
	int*** A0_alb;
	int**** A0_calh;
	int***** A0_callh;
	int** C;


	GRBVar*** x0_abh;
	GRBVar*** x0_ahh;
	GRBVar**** x0_albh;
	GRBVar*** y0_ahb;

	GRBVar****** xt_cablh;
	GRBVar****** xt_cahlh;
	GRBVar******* xt_calblh;
	GRBVar**** yt_ahb;
	GRBVar*** Ct_cl;
	GRBVar*** At_ab;
	GRBVar**** At_alb;

	void add_bases_constraints_t0();
	void add_hospitals_constraints_t0();
	void add_locations_constraints_t0();
	void add_queues_constraints_t0();

	void add_bases_constraints();
	void add_hospitals_constraints();
	void add_locations_constraints();
	void add_queues_constraints();
	void add_ambs_location_constraints();
	void add_base_cap_constraints();


	void solve();

	// int is_at_hospital(Ambulance & amb);
	// int is_at_base(Ambulance & amb);
	// std::pair<int,int> is_at_location_base(Ambulance & amb);
	void debug();
	
};


class AmbulanceModel{
public:
	AmbulanceModel(Data & data, GRBEnv & env, Solver& solver);
	~AmbulanceModel();

	Data & data;
	GRBModel model;
	Solver& solver;

	int amb_type;
	int src_hosp;
	int base_return;
	
	int call_type;
	int call_location;
	int dst_hosp;


	int** A0_ab;
	int** A0_ah;
	int*** A0_alb;
	int**** A0_calh;
	int***** A0_callh;
	int** C;


	GRBVar*** y0_ahb;
	GRBVar***** y0_ahclh;


	GRBVar****** xt_cablh;
	GRBVar****** xt_cahlh;
	GRBVar******* xt_calblh;
	GRBVar**** yt_ahb;
	GRBVar*** Ct_cl;
	GRBVar*** At_ab;
	GRBVar**** At_alb;

	void add_bases_constraints_t0();
	void add_hospitals_constraints_t0();
	void add_locations_constraints_t0();
	void add_queues_constraints_t0();

	void add_bases_constraints();
	void add_hospitals_constraints();
	void add_locations_constraints();
	void add_queues_constraints();
	void add_ambs_location_constraints();
	void add_base_cap_constraints();

	void solve();

	// int is_at_hospital(Ambulance & amb);
	// int is_at_base(Ambulance & amb);
	// std::pair<int,int> is_at_location_base(Ambulance & amb);
	void debug();
};

#endif