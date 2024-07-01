#ifndef _CALL_H
#define _CALL_H

#include "main.h"

class Call{
public:
	Call();
	~Call();
	Call(int id, double time, Location& location, int hospital, int cleaning, int priority,
		bool to_hosp, bool to_clean, double time_on_scene, double time_on_hospital, double cleaning_time);
	int id;
	double time;
	Location location;
	int disc_time;
	int region;
	int hospital;
	int cleaning;
	int priority;
	bool hosp_needed;
	bool clean_needed;
	double time_equipment;
	double time_on_scene;
	double time_at_hospital;
	double cleaning_time;
	double end;
	bool first_aid;
	std::vector<int> ambulances;

	friend std::ostream& operator<<(std::ostream& out, const Call& call);

	bool operator==(const Call& c);
	bool operator!=(const Call& c);
	bool operator<(const Call& c);
};


#endif