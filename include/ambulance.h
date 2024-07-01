#ifndef _AMBULANCE_H
#define _AMBULANCE_H


#include "main.h"

class Simulator;
class Solver;
class Instance;
class Travel;
class Call;

using std::shared_ptr;

class Ambulance{
public:
	Ambulance();
	~Ambulance();

	Ambulance(int id, Location& base_location, Location& free_location,
			double arrival_time_at_c_last_trip,
			double arrival_time_at_f_last_trip,
			double arrival_time_at_b_last_trip, 
			int type, double speed);
	Call* call;
	int id;
	Location base_location;
	Location free_location;
	Location clean_location;
	double arrival_time_at_c_last_trip;
	double arrival_time_at_f_last_trip; //free
	double arrival_time_at_b_last_trip;
	int type;
	double speed;
	bool busy;
	Location last_origin_location;
	double departure_time;

	std::vector<double> times;
	std::vector<Location> trips;
	std::vector<TripType> trip_types;

	double answer_call(Call& call, Travel& travel, Instance& ins, double time, double min_time, 
		int base_id);

	void set_new_point(double time, Location& location, TripType trip_type);
	void reassign_ambulance(Location& base, double arrival_time);


	Location get_current_location(Instance& ins, double time);

	bool is_busy(double time);
	bool is_idle(double time);

	friend std::ostream& operator<<(std::ostream& out, const Ambulance& amb);
	
};


#endif