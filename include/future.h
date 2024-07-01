#ifndef _FUTURE_H
#define _FUTURE_H

#include "main.h"
#include "data.h"
#include "solver.h"


class Solver;
class Data;

class FutureCall{
public:
    FutureCall(Data & data, GRBEnv & env, Solver& solver, int a0, int l0, int b0,
    	int h0, std::vector<Call>& calls);
    ~FutureCall();

	Data & data;
    GRBModel model;
	Solver& solver;
	int type_call0;
	int local_call0;
	int a0, l0, b0, h0;
	std::vector<Call>& calls;
    double obj;


	int** A0_ab;
	int** A0_ah;
	int*** A0_alb;
	int**** A0_calh;
	int***** A0_callh;
	int** C;
	int*** lambda;

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

};

class FutureAmbulance{
public:
    FutureAmbulance(Data & data, GRBEnv & env, Solver& solver, Call* call, int h_dest,
    	int b_ret, std::vector<Call>& calls);
    ~FutureAmbulance();

    Data & data;
	Call* call;
	GRBModel model;
	Solver& solver;
	int h_dest, b_ret;
	std::vector<Call>& calls;
	int c0, l0, h0, a0;
	double obj;


	int** A0_ab;
	int** A0_ah;
	int*** A0_alb;
	int**** A0_calh;
	int***** A0_callh;
	int** C;
	int*** lambda;

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
};


#endif