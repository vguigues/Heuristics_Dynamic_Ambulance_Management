#ifndef _FAST_EQUIPMENT_H
#define _FAST_EQUIPMENT_H

#include "main.h"

class Simulator;
class Solver;
class Instance;
class Travel;
class Call;

using std::shared_ptr;


enum class FastType{
    DRONE,
    MOTORCYCLE
};

class FastEquipment{
public:
    FastEquipment();
    ~FastEquipment();

    FastEquipment(int id, Location& base_location, Location& free_location,
			double arrival_time_at_c_last_trip,
			double arrival_time_at_f_last_trip,
			double arrival_time_at_b_last_trip, 
			FastType type, double speed);

    Call* call;
	int id;
	Location base_location;
	Location free_location;
	Location clean_location;
	double arrival_time_at_c_last_trip;
	double arrival_time_at_f_last_trip; //free
	double arrival_time_at_b_last_trip;
	FastType type;
	double speed;
	bool busy;
	Location last_origin_location;
	double departure_time;

	// int number_equipment;

    std::vector<double> times;
	std::vector<Location> trips;
	std::vector<TripType> trip_types;


    double answer_call(Call& call, Travel& travel, Instance& ins, double time, double min_time, 
		int nearest_base_id);

	void set_new_point(double time, Location& location, TripType trip_type);

	friend std::ostream& operator<<(std::ostream& out, const FastEquipment& fe);
};



#endif