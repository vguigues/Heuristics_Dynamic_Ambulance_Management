#include "../include/call.h"


Call::Call(): id(-1), time(-1), location(null_location), disc_time(-1), region(-1), hospital(-1),
	cleaning(false), priority(-1), hosp_needed(false), clean_needed(false), 
	time_equipment(-1), time_on_scene(-1), time_at_hospital(-1), cleaning_time(-1), 
	end(-1), first_aid(false){}

Call::Call(int id, double time, Location& location, int hospital, int cleaning, 
		int priority, bool hosp_needed, bool clean_needed,
		double time_on_scene, double time_at_hospital, double cleaning_time): id(id), 
		time(time), location(location), disc_time(-1), region(-1), hospital(hospital), 
		cleaning(cleaning), priority(priority), hosp_needed(hosp_needed), 
		clean_needed(clean_needed), time_equipment(600), 
		time_on_scene(time_on_scene), time_at_hospital(time_at_hospital), 
		cleaning_time(cleaning_time), end(-1), first_aid(false){
}



std::ostream& operator<<(std::ostream& out, const Call& call){
	out << "Id" << call.id << ", t" << call.time << ", ";
	// out << "Location: (" << call.location.first << ", " << call.location.second << "), ";
	out << "p" << call.priority << " r" << call.region << " h" << call.hospital;
	return out;
}


bool Call::operator==(const Call& c){
	return id == c.id;
}

bool Call::operator<(const Call& c){
	return time < c.time;
}

bool Call::operator!=(const Call& c){
	return id != c.id;
}


Call::~Call(){}