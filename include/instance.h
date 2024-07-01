#ifndef _INSTANCE
#define _INSTANCE

#include "main.h"
#include "travel.h"

using namespace std;

class Call;
class Ambulance;
class Travel;
class OSRMHelper;

const string simulated_rect_1 = 
"calibration/Scenarios_For_Tests/Simulated_Data_Rectangle/simualtedRectangleUniform.txt";
const string simulated_rect_2 = 
"calibration/Scenarios_For_Tests/Simulated_Data_Rectangle/simualtedRectanglePoisson.txt";



class Instance{
public:
	Instance();
	explicit Instance(string path);
	~Instance();
	
	int nb_scenarios;
	vector<int> nb_calls;
	int nb_hospitals;
	int nb_bases;
	int nb_cleaning_bases;
	int nb_types_ambulance;
	int nb_ambulances;
	int nb_priorities;
	int nb_regions;
	int nb_times;
	int nb_days;
	int total_scenarios;
	vector<vector<double>> penalty_matrix;
	xt::xarray<double> preparedness;
	double max_service_time;

	string path;
	Travel travel;

	double x_min, x_max, y_min, y_max;
	double slot_duration;
	vector<set<int>> A; //for each priority c, set of ambulance types that can respond to c

	std::vector<std::vector<Location>> samples;
	vector<Location> centers;
	vector<Location> bases;
	vector<Location> cleaning_bases;
	vector<Location> hospitals;
	vector<int> cap_bases;
	vector<double> penalties;
	vector<int> nearest_base_to_region;
	vector<int> nearest_region_to_base;
	vector<int> nearest_hospital_to_region;
	vector<int> nearest_base_to_hospital;

	vector<vector<Call>> calls;
	vector<Ambulance> ambulances;
	vector<vector<vector<Call>>> scenarios_by_day;
	vector<vector<int>> neighbors;
	vector<vector<bool>> adj_matrix;
	vector<double> delay_scene_p;
	vector<double> delay_hosp_p;
	vector<int> amb_count_by_type;
	xt::xarray<double> lambda;
	xt::xarray<double> lambda_bases;
	int g0;
	vector<int> time_horizon;

	void load_euclidian();
	void load_real_data();


	void load_lambda();
	void load_lambda_bases();
	void load_regions();
	void read_preparedness();
	void modify_lambdas();

	void print_scenario_stats();
	double survival_function(double time);
	double get_param(int t, int g, int r, int p);
	double get_max_service_time(vector<Call>& calls);

	

};

#endif