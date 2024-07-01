#ifndef _TRAVEL_H
#define _TRAVEL_H

#include "main.h"
#include "osrm_helper.h"

class OSRMHelper;
class Instance;
class Call;
class Ambulance;
class FastEquipment;

class Travel {
 public:
  Travel(bool euclidian = true);

  OSRMHelper osrm;
  bool euclidian;
  bool forward;
  double delay_param = -1;

  void set_forward(bool a_forward);
  void set_delay_param(double u = -1);

  double get_delay(double travel_time, double v);
  double get_response_time(Ambulance& amb, Call& call, double current_time,
                           bool force_forward = false);
  double get_response_time(FastEquipment& fe, Call& call, double current_time,
                           bool force_forward = false);

  Location ambulance_position(Ambulance& amb, double t);
  Location equipment_position(FastEquipment& fe, double t);

  double travel_time(Location& a, Location& b, Ambulance& amb);
  double travel_time(Location& a, Location& b, FastEquipment& fe);

  double travel_time_from_position(Ambulance& amb, Call& call);
  double travel_time_from_position(FastEquipment& fe, Call& call);

  double vehicle_distance(Location& a, Location& b);
  double drone_distance(Location& a, Location& b);

  Location position_between_origin_destination(Location& a, Location& b,
                                               double t0, double t, double tf,
                                               Ambulance& amb);
  Location position_between_origin_destination(Location& a, Location& b,
                                               double t0, double t, double tf,
                                               FastEquipment& fe);

  double euclidian_travel_time(Location& a, Location& b, double speed,
                               bool free_flow = false);
  double lat_long_travel_time(Location& a, Location& b, double v,
                              bool free_flow = false);
  double street_travel_time(Location& a, Location& b, double v,
                            bool free_flow = false);

  double euclidian_distance(Location& a, Location& b);
  double lat_long_distance(Location& a, Location& b);
  double street_distance(Location& a, Location& b);

  double region_service_time(Location& base, Location& r, Location& hospital,
                             double speed);

  std::vector<std::vector<double>> table_in_out(std::vector<Location>& a,
                                                std::vector<Location>& b);
};

#endif