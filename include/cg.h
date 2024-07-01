#ifndef _CG_H
#define _CG_H


#include "main.h"
#include "data.h"
#include "solver.h"

class EventSimulator;
class Solver;

int get_time(int time);

class CGCall{
public:
	CGCall(Data& data, GRBEnv& env);
	~CGCall();

	Data& data;
	GRBEnv& env;
	GRBModel model;


	int columns_added;
	int type_call0;
	int local_call0;

	int src_type;
	int src_location;
	int amb_type;
	int loc_base;
	int dst_hosp;

	int iter;
	// int*** lambda; //moved to Data
	
	int** A0_ab; //clear
	int** A0_ah; //clear
	int*** A0_alb; //clear
	int**** A0_calh; //clear
	int***** A0_callh; //clear
	int** C; //clear

	// std::set<int>*** L_tab; //moved to Data

	std::vector<std::vector<int>> closest_l1_to_l;

	int**** vh_tall; //clear
	int** vc_tl; //clear
	int*** vb_tal; //clear


	GRBVar*** x0_abh; //clear
	GRBVar*** x0_ahh; //clear
	GRBVar**** x0_albh;   //clear
	GRBVar*** y0_ahb;  //clear

	GRBVar****** xt_cablh;  //clear
	GRBVar****** xt_cahlh;  //clear
	GRBVar******* xt_calblh; //clear
	GRBVar**** yt_ahb; //clear
	GRBVar*** Ct_cl; //clear
	GRBVar*** At_ab; //clear
	GRBVar**** At_alb; //clear


	GRBConstr** con_base_t0; //clear
	GRBConstr** con_hospital_t0; //clear
	GRBConstr*** con_location_t0; //clear
	GRBConstr** con_queue_t0; //clear


	GRBConstr*** con_base; //clear
	GRBConstr*** con_hospital; //clear
	GRBConstr**** con_location; //clear
	GRBConstr*** con_queue; //clear
	GRBConstr**** con_amb_location; //clear
	GRBConstr*** con_amb_base; //clear
	GRBConstr** con_cap_base; //clear


	double** beta0_ab; //clear
	double** alpha0_ah; //clear
	double*** psi0_alb; //clear
	double** phi0_cl; //clear
	double*** beta_tab; //clear
	double*** alpha_tah; //clear
	double**** psi_talb; //clear
	double*** phi_tcl; //clear
	double** ups_tb; //clear
	double**** theta_talb; //clear
	double*** gamma_tab; //clear



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


	int c_tl(int t, int l);
	int b_tal(int t, int a, int l);
	int h_tall(int t, int a, int l1, int l);

	void solve();

	void set_dual();
	void add_column(std::vector<int>& column);
	double sub_problem(int t, int a, int l, int l1, int c, int h, int b);
	bool pricing();
	void print_solution();


	void fix_abh(int a, int b, int h, double value){
		x0_abh[a][b][h].set(GRB_DoubleAttr_LB, value);
		x0_abh[a][b][h].set(GRB_DoubleAttr_UB, value);
	}

	void fix_ahh(int a, int h1, int h, double value){
		x0_ahh[a][h1][h].set(GRB_DoubleAttr_LB, value);
		x0_ahh[a][h1][h].set(GRB_DoubleAttr_UB, value);
	}

	void fix_albh(int a, int l, int b, int h, double value){
		x0_albh[a][l][b][h].set(GRB_DoubleAttr_LB, value);
		x0_albh[a][l][b][h].set(GRB_DoubleAttr_UB, value);
	}

	void fix_queue(){
		for(auto a: data.A[data.type_call0]){
			for(auto h: data.H[data.type_call0][data.local_call0]){
				for(int l = 0; l < data.num_locals; ++l){
					for(int b = 0; b < data.num_bases; ++b){
						fix_albh(a,l,b,h,0);
					}
				}

				for(int b = 0; b < data.num_bases; ++b){
					fix_abh(a,b,h,0);
				}

				for(int h1 = 0; h1 < data.num_hosps; ++h1){
					fix_ahh(a,h1,h,0);
				}
			}
		}
	}

	// int is_at_hospital(Ambulance & amb);
	// int is_at_base(Ambulance & amb);
	// std::pair<int,int> is_at_location_base(Ambulance & amb);
	
};


class CGAmbulance{
public:
	CGAmbulance(Data & data, GRBEnv & env);
	~CGAmbulance();

	Data& data;
	GRBEnv& env;
	GRBModel model;
	int columns_added;

	int amb_type;
	int src_hosp;
	int base_return;
	
	int call_type;
	int call_location;
	int dst_hosp;

	int iter;
	
	int** A0_ab; //clear
	int** A0_ah; //clear
	int*** A0_alb; //clear
	int**** A0_calh; //clear
	int***** A0_callh; //clear
	int** C; //clear

	std::vector<std::vector<int>> closest_l1_to_l;
	int**** vh_tall; //clear
	int** vc_tl; //clear
	int*** vb_tal; //clear

	GRBVar*** y0_ahb; //clear
	GRBVar***** y0_ahclh; //clear


	GRBVar****** xt_cablh; //clear
	GRBVar****** xt_cahlh; //clear
	GRBVar******* xt_calblh; //clear
	GRBVar**** yt_ahb; //clear
	GRBVar*** Ct_cl; //clear
	GRBVar*** At_ab; //clear
	GRBVar**** At_alb; //clear

	GRBConstr** con_base_t0; //clear
	GRBConstr** con_hospital_t0; //clear
	GRBConstr*** con_location_t0; //clear
	GRBConstr** con_queue_t0; //clear

	GRBConstr*** con_base; //clear
	GRBConstr*** con_hospital; //clear
	GRBConstr**** con_location; //clear
	GRBConstr*** con_queue; //clear
	GRBConstr**** con_amb_location; //clear
	GRBConstr*** con_amb_base; //clear
	GRBConstr** con_cap_base; //clear


	double** beta0_ab; //clear
	double** alpha0_ah; //clear
	double*** psi0_alb; //clear
	double** phi0_cl; //clear
	double*** beta_tab; //clear
	double*** alpha_tah; //clear
	double**** psi_talb; //clear
	double*** phi_tcl; //clear
	double** ups_tb; //clear
	double**** theta_talb; //clear
	double*** gamma_tab; //clear


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

	int c_tl(int t, int l);
	int b_tal(int t, int a, int l);
	int h_tall(int t, int a, int l1, int l);
	
	void solve();

	void fix_call(int a, int c, int h1, int l, int h, double value){
		y0_ahclh[a][c][h1][l][h].set(GRB_DoubleAttr_LB, value);
		y0_ahclh[a][c][h1][l][h].set(GRB_DoubleAttr_UB, value);
	}

	void fix_return(int a, int h, int b, double value){
		y0_ahb[a][h][b].set(GRB_DoubleAttr_LB, value);
		y0_ahb[a][h][b].set(GRB_DoubleAttr_UB, value);
	}

	void set_dual();
	void add_column(std::vector<int>& column);
	double sub_problem(int t, int a, int l, int l1, int c, int h, int b);
	bool pricing();
	void print_solution();



	// int is_at_hospital(Ambulance & amb);
	// int is_at_base(Ambulance & amb);
	// std::pair<int,int> is_at_location_base(Ambulance & amb);
};

typedef boost::multi_array<int, 2> Int2D;
typedef boost::multi_array<int, 3> Int3D;
typedef boost::multi_array<int, 4> Int4D;
typedef boost::multi_array<int, 5> Int5D;

typedef boost::multi_array<double, 2> Double2D;
typedef boost::multi_array<double, 3> Double3D;
typedef boost::multi_array<double, 4> Double4D;

typedef boost::multi_array<GRBVar, 2> Var2D;
typedef boost::multi_array<GRBVar, 3> Var3D;
typedef boost::multi_array<GRBVar, 4> Var4D;
typedef boost::multi_array<GRBVar, 5> Var5D;
typedef boost::multi_array<GRBVar, 6> Var6D;
typedef boost::multi_array<GRBVar, 7> Var7D;


typedef boost::multi_array<GRBConstr, 2> Con2D;
typedef boost::multi_array<GRBConstr, 3> Con3D;
typedef boost::multi_array<GRBConstr, 4> Con4D;


class SmallCGCall{
public:
	SmallCGCall(Data& data, GRBEnv& env, Solver& solver);
	~SmallCGCall();
	
	Data& data;
	GRBEnv& env;
	GRBModel model;
	Solver& solver;


	int src_type;
	int src_location;
	int amb_type;
	int loc_base;
	int dst_hosp;

	int type_call0;
	int local_call0;
	int iter;
	int columns_added;

	Int3D lambda;
	
	Int2D A0_ab;
	Int2D A0_ah;
	Int3D A0_alb;
	Int4D A0_calh;
	Int5D A0_callh;
	Int2D C;

	std::set<int>*** L_tab;

	std::vector<std::vector<int>> closest_l1_to_l;

	Int4D vh_tall;
	Int2D vc_tl;
	Int3D vb_tal;


	Var3D x0_abh;
	Var3D x0_ahh;
	Var4D x0_albh;
	Var3D y0_ahb;

	Var6D xt_cablh;
	Var6D xt_cahlh;
	Var7D xt_calblh;
	Var4D yt_ahb;
	Var3D Ct_cl;
	Var3D At_ab;
	Var4D At_alb;


	Con2D con_base_t0;
	Con2D con_hospital_t0;
	Con3D con_location_t0;
	Con2D con_queue_t0;


	Con3D con_base;
	Con3D con_hospital;
	Con4D con_location;
	Con3D con_queue;
	Con4D con_amb_location;
	Con3D con_amb_base;
	Con2D con_cap_base;


	double** beta0_ab;
	double** alpha0_ah;
	double*** psi0_alb;
	double** phi0_cl;
	double*** beta_tab;
	double*** alpha_tah;
	double**** psi_talb;
	double*** phi_tcl;
	double** ups_tb;
	double**** theta_talb;
	double*** gamma_tab;



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


	int c_tl(int t, int l);
	int b_tal(int t, int a, int l);
	int h_tall(int t, int a, int l1, int l);

	void solve();

	void set_dual();
	void add_column(std::vector<int>& column);
	double sub_problem(int t, int a, int l, int l1, int c, int h, int b);
	bool pricing();
	void print_solution();

	// int is_at_hospital(Ambulance & amb);
	// int is_at_base(Ambulance & amb);
	// std::pair<int,int> is_at_location_base(Ambulance & amb);

};


#endif