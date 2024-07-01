#ifndef _SOLVER_H
#define _SOLVER_H

using namespace std;

#include "main.h"
#include "call.h"
#include "ambulance.h"
#include "travel.h"

class Travel;
class Instance;
class CGCall;

typedef vector<pair<double, int>> SortableVector;
typedef struct{
	double free_time;
	Location free_location;
}FinishServiceData;

class Solver{
protected:
	Solver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, double time);
public:
	virtual ~Solver();

	GRBEnv& env;
	vector<Call> calls;
	vector<Ambulance> ambulances;
	Travel travel;
	Instance& ins;

	double time;
	double run_time;
	int index_solver;
	int event_call;
	int index_call;
	int calls_attended;
	double first_time; //time of first call
	double begin_time; // beginning time of simulation
	double last_time;
	bool is_prepared;
	bool debug_mode = false;

	vector<double> waiting_on_scene;
	vector<double> waiting_on_scene_penalized;
	vector<double> waiting_to_hospital;
	vector<double> calls_end;
	vector<int> which_ambulance;

	vector<double> times;
	vector<vector<vector<double>>> lambda_base;
	double time_ahead;

	vector<int> queue;
	vector<int> nearest_base;
	double obj;
	int released_amb;

	const vector<int> p_weight{4,2,1};

	std::vector<double> run_times;

	virtual void run() = 0;
	virtual void prepare();
	void set_next_event(int& event_call, int& index_call);
	void set_calls_nearest_bases();
	void print_results();
	bool can_answer(Ambulance& amb, Call& call);
	pair<int, double> get_closest_call(Ambulance& amb);
	pair<int, double> get_oldest_call(Ambulance& amb);
	vector<double> get_lambda(int g, int t_begin);
	vector<double> prep(vector<double>& lambda);
	SortableVector time_to_bases(Location& free_location);
	SortableVector prep_to_bases(SortableVector& time_to_base, int amb_id, size_t max_size, double free_time, vector<double>& lambda, bool forward);
	virtual int get_return_base(int amb_id, FinishServiceData& finish_service_data, vector<double>& lambda, bool forward);
	FinishServiceData get_finish_service_data(Ambulance& amb, Call& call, double time);
	double penalized_response_time(double response_time, int amb_type, int call_type);
	int find_best_base(int amb_id, bool debug = false);
	double poisson_quantile(double q, double lambda, bool debug = false);
	unsigned long int fac(int n);
	void set_debug_mode(bool debug);

	int get_best_index_amb(vector<int>& index_ambs);
};


class ForwardSolver: public Solver{
public:
	ForwardSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, double time);
	
	virtual void run();
};

class QueueSolver: public Solver{
public:
	QueueSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, double time);
	
	virtual void run();
};

class PreparednessSolver: public Solver{
public:
	PreparednessSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, int g, double time);

	int g;
	virtual void run();
	virtual int get_return_base(int amb_id, vector<double>& lambda, vector<int>&  available);
};


class ForwardPrepSolver: public Solver{
public:
	ForwardPrepSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, int g, double time);
	int g;
	
	virtual void run();
};


class MinMaxPrepSolver: public Solver{
public:
	MinMaxPrepSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, int g, double time);
	int g;
	
	virtual void run();
};

class PriorityPrepSolver: public Solver{
public:
	PriorityPrepSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, int g, double time);
	int g;
	virtual void run();
};

class Prep2Solver: public Solver{
public:
	Prep2Solver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, int g, double time);
	
	int g;
	virtual void run();
	vector<vector<bool>> get_relocations_exact(vector<int>& reloc_available);
	vector<vector<bool>> get_relocations_heuristic(vector<int>& reloc_available, int
		min_prep_r, vector<double>& lambda);
	virtual int get_return_base(int amb_id, vector<double>& lambda, vector<int>&  available);
};

class OrderedSolver: public Solver{
public:
	OrderedSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, int g, double time);

	int g;
	virtual void run();
	virtual void prepare();
	vector<int> get_ordered_table(vector<double>& lambda);
};



struct DistrictData{
	unordered_map<int, vector<int>> districts;
	vector<int> D; //District of region r
	vector<int> ambulance_district; // district of ambulance a 

	vector<int> ambulances_of_district(int d) const{
		vector<int> result;
		for(size_t i = 0; i < ambulance_district.size(); ++i){
			if(ambulance_district[i] == d){
				result.push_back(i);
			}
		}
		return result;
	}
};

class DistrictManager{
public:
	Instance& ins;
	Travel& travel;
	int g;
	vector<unsigned long int> tab_fac;
	vector<double> lambda;
	vector<Ambulance>& ambulances;

	DistrictManager(Instance& ins, Travel& travel, int g, vector<Ambulance>& ambulances);

	unsigned long int fac(int n);
	double total_aexc(DistrictData& district_data, unordered_map<int, double>& rho);
	double aexc(DistrictData& district_data, unordered_map<int, double>& rho, int i);
	DistrictData get_districts();
	DistrictData get_m_districts();
	unordered_map<int, double> get_workload(DistrictData& district_data);
	void local_search(DistrictData& district_data, unordered_map<int, double>& W);
	void swap(DistrictData& district_data, int dk, int ds, int r, int s);
	DistrictData merge_step(DistrictData& district_data, unordered_map<int, double>& rho);
	DistrictData merge(DistrictData& district_data, int di, int dj);
	bool can_merge(unordered_map<int,vector<int>>& districts, int i, int j);
	int get_K(vector<vector<int>>& districts);
	double normalize(unordered_map<int, double>& rho, double r);
	double Q(double rho, int j, int N);
	double get_mean_rate_of_district(vector<int>& district, int b);
	double var(unordered_map<int, double>& rho);

	void print_districts(DistrictData& data){
		fmt::print("Number of districts = {}\n", data.districts.size());
		for(auto& kv: data.districts){
			int i = kv.first;
			auto& district = kv.second;
			auto ambulances = data.ambulances_of_district(i);
			if(district.size() > 0){
				fmt::print("District {} regions: {}\n", i, district);
				fmt::print("District {} ambulances: {}\n", i, ambulances);
			}
		}
		fmt::print("D = {}\n", data.D);
		fmt::print("Ambulance districts = {}\n", data.ambulance_district);
	}
};


class DistrictSolver: public Solver{
public:
	DistrictSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, int g, DistrictData& a_data, vector<double>& a_lambda, double time);

	int g;
	DistrictData& data;
	vector<double>& lambda;


	virtual void run();
	virtual void prepare();

	void closest_cross();
	void heuristic_cross();

	unordered_map<int, vector<int>> get_ordered_table_intra();
	unordered_map<int, vector<int>> get_ordered_table_inter();

	void print_districts(DistrictData& data){
		fmt::print("Number of districts = {}\n", data.districts.size());
		for(auto& kv: data.districts){
			int i = kv.first;
			auto& district = kv.second;
			auto ambulances = data.ambulances_of_district(i);
			if(district.size() > 0){
				fmt::print("District {} regions: {}\n", i, district);
				fmt::print("District {} ambulances: {}\n", i, ambulances);
			}
		}
		fmt::print("D = {}\n", data.D);
		fmt::print("Ambulance districts = {}\n", data.ambulance_district);
	}

};

class CoverageSolver: public Solver{
public:
	CoverageSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, int g, double time);

	int g;
	double T;
	double q;
	virtual void run();
	int get_index_amb(vector<int>& k, vector<int>& A, vector<double>& lambda);
	int get_k(Location& loc_i, vector<int>& A);
};



class PrioritySolver: public Solver{
public:
	PrioritySolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, double time);

	virtual void run();
};

class MinMaxSolver: public Solver{
public:
	MinMaxSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, double time);

	virtual void run();	
};


class MinMaxPSolver: public Solver{
public:
	MinMaxPSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, double time);

	virtual void run();
};

class NonMiopycSolver: public Solver{
public:
	NonMiopycSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, double time);

	virtual void run();
	void compute_min_times(int k, vector<double>& min_times, vector<int>& index_amb, 
		vector<int> & amb_type);

	vector<int> get_best_busy_ambs(int k);

	int found_ambulance(vector<vector<int>>& best_ambs, int & max_index, int i, int index_amb, 
		vector<bool>& is_allocated, vector<vector<double>>& travel_times, vector<vector<double>>& pen_travel_times,
		vector<double>& min_times, vector<double>& min_times_p, vector<bool>& is_call_on_queue, bool debug);

	void compute_min_times(int i, int max_index, vector<vector<double>>& travel_times, 
		vector<double>& min_times, vector<vector<int>>& best_ambs);
};

class GenForwardSolver: public Solver{
public:
	GenForwardSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, double time);
	virtual void run();
};


class GenMinMaxSolver: public Solver{
public:
	GenMinMaxSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, double time);

	virtual void run();
};


class DeterministicSolver: public Solver{
public:
	DeterministicSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, double time);

	virtual void run();
};

class CGSolver: public Solver{
public:
	CGSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, double time);

	int get_return_base(CGCall& cg, Data& data, int t, int c, int a, int l1, 
		int l, int h);
	virtual void run();
};


enum TypeNode{
	BASE_NODE,
	HOSPITAL_NODE,
	CALL_NODE,
	RELOCATION_NODE
};

typedef struct{
	TypeNode node_type;
	int node_index;
}PathNode;

typedef struct{
	int ambulance_id;
	vector<PathNode> node_sequence;
	double cost;
}Route;


struct EnumNode{
	Call& call;
	double waiting_time;
	shared_ptr<EnumNode> prev;
	double time_finish_service;
	int depth;
	
	EnumNode(Call& a_call, double a_waiting_time, shared_ptr<EnumNode>& a_prev, double a_time_finish_service): call{a_call}, waiting_time{a_waiting_time}, prev{a_prev},
		time_finish_service{a_time_finish_service}{
		if(prev == nullptr){
			depth = 0;
		}else{
			depth = prev->depth+1;
		}
	}
};

class EnumerateSolver: public Solver{
public:
	EnumerateSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, int g, double time);
	
	bool debug = false;
	GRBModel model;
	vector<GRBVar> x;
	vector<GRBVar> y;
	vector<set<int>> routes_by_ambulance;
	vector<set<int>> routes_by_call;
	vector<Route> routes;
	std::queue<EnumNode> node_queue;
	vector<bool> is_call_on_queue;

	virtual void run();
	vector<Route> generate_routes(Ambulance& amb);
	Route expand_route(Ambulance& amb, EnumNode& route_node);
	void load_model();
	void insert_enum_node(Ambulance& amb, Call& call, shared_ptr<EnumNode> parent);
	double get_call_service_time(Call& call, Ambulance& amb);
	Route get_route(Ambulance& amb, EnumNode& route_node);
	double missed_call_penalty = 36000;
};

class MCSolver: public Solver{
public:
	MCSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel,  int g, double time);
	int g;
	xt::xarray<double> preparedness;
	std::vector<double> run_times_selection;

	int get_time_slot(double time);
	pair<vector<vector<int>>, vector<int>> partition_ambulances(int a0 = -1);
	bool solve_selection(int i0, bool debug);
	bool solve_reassignment(int a0, bool debug);
	void read_preparedness();

	double s_minus(int amb_type, int b, double a_time, vector<int>& supply_b);
	double s_plus(int amb_type, int b, double a_time, vector<int>& supply_b);

	virtual void run();
};


class QueueDeficitSolver: public Solver{
public: 
	QueueDeficitSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
		Instance& ins, Travel& travel, int g, double time);
	int g;
	bool debug = false;
	virtual void run();
};

#endif

// class ModelSolver: public Solver{
// public:
// 	ModelSolver(GRBEnv& env, vector<Call>& calls, vector<Ambulance>& ambulances, 
// 		Instance& ins, Travel& travel);

// 	virtual void run();
// };