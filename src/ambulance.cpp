#include "../include/instance.h"
#include "../include/travel.h"
#include "../include/solver.h"
#include "../include/call.h"
#include "../include/ambulance.h"


Ambulance::Ambulance(int id, Location& base_location, Location& free_location,
		double arrival_time_at_c_last_trip, double arrival_time_at_f_last_trip, 
		double arrival_time_at_b_last_trip, int type, double speed): 
		call(NULL), id(id), base_location(base_location), 
		free_location(free_location), clean_location(null_location),
		arrival_time_at_c_last_trip(arrival_time_at_c_last_trip),
		arrival_time_at_f_last_trip(arrival_time_at_f_last_trip),
		arrival_time_at_b_last_trip(arrival_time_at_b_last_trip),
		type(type), speed(speed), busy(false), 
		last_origin_location(base_location), departure_time(0){
}

Location Ambulance::get_current_location(Instance& ins, double time){
	Travel& travel = ins.travel;
	if(arrival_time_at_b_last_trip <= time){
		return base_location;
	}else if(arrival_time_at_f_last_trip <= time){
		return travel.ambulance_position(*this, time);
	}else{
		double elapsed = departure_time + travel.travel_time(last_origin_location, 
			call->location, *this);
		if(elapsed > time){
			return travel.position_between_origin_destination(last_origin_location, 
				call->location, departure_time, time, arrival_time_at_c_last_trip, *this);
		}
	
		elapsed += call->time_on_scene;
		if(elapsed > time){
			return call->location;
		}
		double time_leave_call = elapsed;
		if(call->hosp_needed && call->clean_needed){
			//check if between call and hosp
			elapsed += travel.travel_time(call->location, ins.hospitals[call->hospital], *this);
			if(elapsed > time){
				return travel.position_between_origin_destination(call->location, 
					ins.hospitals[call->hospital], time_leave_call, time, arrival_time_at_f_last_trip - call->time_at_hospital, *this);
			}
			//check if at hosp
			elapsed += call->time_at_hospital;
			if(elapsed > time){
				return ins.hospitals[call->hospital];
			}
			//check if between hosp and clean
			double time_leave_hosp = elapsed;
			elapsed += travel.travel_time(ins.hospitals[call->hospital], 
				ins.cleaning_bases[call->cleaning], *this);
			if(elapsed > time){
				return travel.position_between_origin_destination(ins.hospitals[call->hospital], 
					ins.cleaning_bases[call->cleaning], time_leave_hosp, time, elapsed, *this);
			}
			///check if at clean
			elapsed += call->cleaning_time;
			if(elapsed > time){
				return ins.cleaning_bases[call->cleaning];
			}
		}else if(call->hosp_needed && !call->clean_needed){
			//check if between call and hosp
			elapsed += travel.travel_time(call->location, ins.hospitals[call->hospital], *this);
			if(elapsed > time){
				return travel.position_between_origin_destination(call->location, 
					ins.hospitals[call->hospital], time_leave_call, time, elapsed, *this);
			}
			//check if at hosp
			elapsed += call->time_at_hospital;
			if(elapsed > time){
				return ins.hospitals[call->hospital];
			}
		}else if(!call->hosp_needed && call->clean_needed){
			//check if between call and clean
			elapsed += travel.travel_time(call->location, ins.cleaning_bases[call->cleaning], *this);
			if(elapsed > time){
				return travel.position_between_origin_destination(call->location, 
					ins.cleaning_bases[call->cleaning], time_leave_call, time, elapsed, *this);
			}
			///check if at clean
			elapsed += call->cleaning_time;
			if(elapsed > time){
				return ins.cleaning_bases[call->cleaning];
			}
		}

		if(elapsed >= arrival_time_at_f_last_trip){
			return free_location;
		}
	}
	fmt::print("ERROR: could not find current location for amb {} at time {}\n", id, time);
	cin.get();
	return null_location;
}

bool Ambulance::is_busy(double time){ return arrival_time_at_f_last_trip > time; }
bool Ambulance::is_idle(double time){ return arrival_time_at_f_last_trip <= time; }

double Ambulance::answer_call(Call& call, Travel& travel, Instance& ins, double time, 
	double min_time, int base_id){
	this->call = &call;
	departure_time = time;

	if(arrival_time_at_b_last_trip <= time){
		last_origin_location = base_location;
		set_new_point(time, base_location, TripType::AT_BASE);
	}else if(arrival_time_at_f_last_trip < time){
		auto current_location = travel.ambulance_position(*this, time);
		last_origin_location = current_location;
		times.back() = time;
		trips.back() = current_location;
	}else{
		last_origin_location = free_location;
		times.erase(times.begin() + times.size() - 1);
		trips.erase(trips.begin() + trips.size() - 1);
		trip_types.erase(trip_types.begin() + trip_types.size()-2, 
			trip_types.end());
	}

	double time_arrival_scene = time + min_time;
	set_new_point(time_arrival_scene, call.location, TripType::TO_CALL);

	double time_leave_scene = time_arrival_scene + call.time_on_scene;
	set_new_point(time_leave_scene, call.location, TripType::AT_SCENE);

	double waiting_to_hospital_i = GRB_INFINITY;
	Location base;
	if(!g_params.h_use_fixed_bases){
		base = ins.bases[base_id];
	}else{
		base = base_location;
	}

	if(call.clean_needed && call.hosp_needed){
		auto& hospital = ins.hospitals[call.hospital];
		auto& cb = ins.cleaning_bases[call.cleaning];
		double t_from_c_to_h = travel.travel_time(call.location, hospital,*this);
		double t_from_h_to_cb = travel.travel_time(hospital, cb, *this);
		double t_from_cb_to_b = travel.travel_time(cb,base, *this);
		double time_leave_hospital = time_leave_scene + t_from_c_to_h +
			call.time_at_hospital;
		double time_arrival_cb = time_leave_hospital + t_from_h_to_cb;
		double time_leave_cb = time_arrival_cb + call.cleaning_time;
		double time_arrival_base = time_leave_cb + t_from_cb_to_b;

		call.end = time_leave_scene + t_from_c_to_h;
		waiting_to_hospital_i = call.time_on_scene + t_from_c_to_h;

		free_location = cb;
		base_location = base;
		arrival_time_at_f_last_trip = time_leave_cb;
		arrival_time_at_b_last_trip = time_arrival_base;

		set_new_point(time_leave_hospital - call.time_at_hospital, hospital,
			TripType::TO_HOSPITAL);
		set_new_point(time_leave_hospital, hospital, TripType::AT_HOSPITAL);

		set_new_point(time_arrival_cb, cb, TripType::TO_CB);
		set_new_point(time_leave_cb, cb, TripType::AT_CB);

		set_new_point(time_arrival_base, base_location, TripType::TO_BASE);
	}else if(call.clean_needed && !call.hosp_needed){
		auto& cb = ins.cleaning_bases[call.cleaning];
		double t_from_c_to_cb = travel.travel_time(call.location, cb, *this);
		double t_from_cb_to_b = travel.travel_time(cb,base, *this);
		double time_arrival_cb = time_leave_scene + t_from_c_to_cb;
		double time_leave_cb = time_arrival_cb + call.cleaning_time;
		double time_arrival_base = time_leave_cb + t_from_cb_to_b;

		call.end = time_leave_scene;
		waiting_to_hospital_i = call.time_on_scene;

		free_location = cb;
		base_location = base;
		arrival_time_at_f_last_trip = time_leave_cb;
		arrival_time_at_b_last_trip = time_arrival_base;


		set_new_point(time_arrival_cb, cb, TripType::TO_CB);
		set_new_point(time_leave_cb, cb, TripType::AT_CB);
		set_new_point(time_arrival_base, base_location, TripType::TO_BASE);

	}else if(!call.clean_needed && call.hosp_needed){
		auto& hospital = ins.hospitals[call.hospital];
		double t_from_c_to_h = travel.travel_time(call.location, hospital, *this);
		double t_from_h_to_b = travel.travel_time(hospital, base, *this);

		double time_leave_hospital = time_leave_scene + t_from_c_to_h +
			call.time_at_hospital;

		double time_arrival_base = time_leave_hospital + t_from_h_to_b;

		call.end = time_leave_scene + t_from_c_to_h;
		waiting_to_hospital_i = call.time_on_scene + t_from_c_to_h;

		free_location = hospital;
		base_location = base;
		arrival_time_at_f_last_trip = time_leave_hospital;
		arrival_time_at_b_last_trip = time_arrival_base;


		set_new_point(time_leave_hospital - call.time_at_hospital, hospital,
			TripType::TO_HOSPITAL);
		set_new_point(time_leave_hospital, hospital, TripType::AT_HOSPITAL);
		set_new_point(time_arrival_base, base_location, TripType::TO_BASE);

	}else if(!call.clean_needed && !call.hosp_needed){
		double t_from_c_to_b = travel.travel_time(call.location, base, *this);

		double time_arrival_base = time_leave_scene + t_from_c_to_b;

		call.end = time_leave_scene;
		waiting_to_hospital_i = call.time_on_scene;

		free_location = call.location;
		base_location = base;
		arrival_time_at_f_last_trip = time_leave_scene;
		arrival_time_at_b_last_trip = time_arrival_base;

		set_new_point(time_arrival_base, base_location, TripType::TO_BASE);
	}

	return waiting_to_hospital_i;
}

void Ambulance::reassign_ambulance(Location& base, double arrival_time){
	arrival_time_at_b_last_trip = arrival_time;
	base_location = base;
	times.pop_back();
	trips.pop_back();
	trip_types.pop_back();

	set_new_point(arrival_time, base, TripType::TO_BASE);
}




void Ambulance::set_new_point(double time, Location& location, TripType trip_type){
	// fmt::print("time location {} {} {}\n", time, location.first, location.second);
	times.push_back(time);
	trips.push_back(location);
	trip_types.push_back(trip_type);
}





std::ostream& operator<<(std::ostream& out, const Ambulance& amb){
	out << "Amb " << amb.id << ", type " << amb.type << " | ";
	out << "Arr c " << amb.arrival_time_at_c_last_trip << ", ";
	out << "Arr b " << amb.arrival_time_at_b_last_trip << ", Arr f ";
	out << amb.arrival_time_at_f_last_trip;
	return out;
}

Ambulance::~Ambulance(){}