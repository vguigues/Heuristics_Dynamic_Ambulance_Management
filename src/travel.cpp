#include "../include/travel.h"

#include "../include/ambulance.h"
#include "../include/call.h"
#include "../include/fast_equipment.h"
#include "../include/osrm_helper.h"

Travel::Travel(bool euclidian)
    : osrm(), euclidian(euclidian), forward(false), delay_param(-1) {}

void Travel::set_forward(bool a_forward) { forward = a_forward; }

double Travel::vehicle_distance(Location &a, Location &b) {
  return lat_long_distance(a, b);
}

double Travel::drone_distance(Location &a, Location &b) {
  return euclidian_distance(a, b);
}

void Travel::set_delay_param(double u) { delay_param = u; }

double Travel::get_delay(double travel_time, double v) {
  if (delay_param < 0) {
    return 0;
  }
  std::default_random_engine gen(600);
  double lam = -1 * travel_time * log(1 - delay_param);
  std::exponential_distribution<double> exp(1 / lam);
  return exp(gen);
}

double Travel::get_response_time(Ambulance &amb, Call &call,
                                 double current_time, bool force_forward) {
  double time_this_ambulance = GRB_INFINITY;

  if (amb.arrival_time_at_b_last_trip <= current_time) {
    if (euclidian) {
      time_this_ambulance =
          euclidian_travel_time(amb.base_location, call.location, amb.speed);
    } else {
      time_this_ambulance =
          street_travel_time(amb.base_location, call.location, 3.6 * amb.speed);
    }
  } else if (amb.arrival_time_at_f_last_trip <= current_time) {
    Location current_location = ambulance_position(amb, current_time);
    if (euclidian) {
      time_this_ambulance =
          euclidian_travel_time(current_location, call.location, amb.speed);
    } else {
      time_this_ambulance =
          street_travel_time(current_location, call.location, 3.6 * amb.speed);
    }
  } else if (forward || force_forward) {
    double time_to_free = amb.arrival_time_at_f_last_trip - current_time;
    double time_free_to_call = GRB_INFINITY;
    if (euclidian) {
      time_free_to_call =
          euclidian_travel_time(amb.free_location, call.location, amb.speed);
    } else {
      time_free_to_call =
          street_travel_time(amb.free_location, call.location, 3.6 * amb.speed);
    }
    time_this_ambulance = time_to_free + time_free_to_call;
  }

  return time_this_ambulance;
}

double Travel::get_response_time(FastEquipment &fe, Call &call,
                                 double current_time, bool force_forward) {
  double time_this_ambulance = GRB_INFINITY;

  if (fe.arrival_time_at_b_last_trip <= current_time) {
    if (euclidian) {
      time_this_ambulance =
          euclidian_travel_time(fe.base_location, call.location, fe.speed);
    } else if (fe.type == FastType::DRONE) {
      time_this_ambulance =
          lat_long_travel_time(fe.base_location, call.location, 3.6 * fe.speed);
    } else {
      time_this_ambulance =
          street_travel_time(fe.base_location, call.location, 3.6 * fe.speed);
    }
  } else if (fe.arrival_time_at_f_last_trip <= current_time) {
    Location current_location = equipment_position(fe, current_time);
    if (euclidian) {
      time_this_ambulance =
          euclidian_travel_time(current_location, call.location, fe.speed);
    } else if (fe.type == FastType::DRONE) {
      time_this_ambulance =
          lat_long_travel_time(fe.base_location, call.location, 3.6 * fe.speed);
    } else {
      time_this_ambulance =
          street_travel_time(current_location, call.location, 3.6 * fe.speed);
    }
  } else if (forward || force_forward) {
    double time_to_free = fe.arrival_time_at_f_last_trip - current_time;
    double time_free_to_call = GRB_INFINITY;
    if (euclidian) {
      time_free_to_call =
          euclidian_travel_time(fe.free_location, call.location, fe.speed);
    } else if (fe.type == FastType::DRONE) {
      time_free_to_call =
          lat_long_travel_time(fe.free_location, call.location, 3.6 * fe.speed);
    } else {
      time_free_to_call =
          street_travel_time(fe.free_location, call.location, 3.6 * fe.speed);
    }
    time_this_ambulance = time_to_free + time_free_to_call;
  }

  return time_this_ambulance;
}

Location Travel::ambulance_position(Ambulance &amb, double t) {
  Location free_location = amb.free_location;
  double t0 = amb.arrival_time_at_f_last_trip;

  if (euclidian) {
    double d = euclidian_distance(free_location, amb.base_location);
    double ttravel = d / amb.speed;
    Location location = free_location;
    if (t < t0 + ttravel) {
      double x1, x2, y1, y2;
      std::tie(x1, x2) = free_location;
      std::tie(y1, y2) = amb.base_location;
      Location direction;
      direction.first = (y1 - x1) / d;
      direction.second = (y2 - x2) / d;
      location.first = x1 + amb.speed * direction.first * (t - t0);
      location.second = x2 + amb.speed * direction.second * (t - t0);
      return location;
    } else {
      return amb.base_location;
    }
  } else {
    const double R = 6371;
    const double radian = M_PI / 180;
    const double v = 60;

    double total_distance = lat_long_distance(free_location, amb.base_location);
    double current_total_time =
        amb.arrival_time_at_b_last_trip - amb.arrival_time_at_f_last_trip;
    double reduced_v = total_distance / current_total_time;

    double p11 = R * cos(radian * free_location.first) *
                 cos(radian * free_location.second);
    double p12 = R * cos(radian * free_location.first) *
                 sin(radian * free_location.second);
    double p13 = R * sin(radian * free_location.first);

    double p21 = R * cos(radian * amb.base_location.first) *
                 cos(radian * amb.base_location.second);
    double p22 = R * cos(radian * amb.base_location.first) *
                 sin(radian * amb.base_location.second);
    double p23 = R * sin(radian * amb.base_location.first);

    double d = sqrt(pow(p11 - p21, 2) + pow(p12 - p22, 2) + pow(p13 - p23, 2));
    double alpha = 2 * asin(d / (2 * R));
    double dearth = R * alpha;
    double ttravel = dearth / reduced_v;
    if (abs(current_total_time - ttravel) > 1) {
      fmt::print("WEIRD: current_total_time = {} != ttravel {}\n",
                 current_total_time, ttravel);
      std::cin.get();
    }
    double t0 = amb.arrival_time_at_f_last_trip;
    Location location = amb.free_location;
    if (t < t0 + ttravel) {
      double alpha0 = (t - t0) * reduced_v / R;
      double beta = sin(alpha - alpha0) / sin(alpha);
      double gamma =
          cos(alpha - alpha0) - sin(alpha - alpha0) * cos(alpha) / sin(alpha);
      double curr_p1 = beta * p11 + gamma * p21;
      double curr_p2 = beta * p12 + gamma * p22;
      double curr_p3 = beta * p13 + gamma * p23;
      location.first = (asin(curr_p3 / R) / M_PI) * 180;
      if (curr_p2 > -g_params.EPS) {
        location.second =
            (acos(curr_p1 / sqrt(pow(R, 2) - pow(curr_p3, 2))) / M_PI) * 180;
      } else {
        location.second =
            -(acos(curr_p1 / sqrt(pow(R, 2) - pow(curr_p3, 2))) / M_PI) * 180;
      }
      return location;
    } else {
      return amb.base_location;
    }
  }
}

Location Travel::equipment_position(FastEquipment &fe, double t) {
  Location free_location = fe.free_location;
  double t0 = fe.arrival_time_at_f_last_trip;

  if (euclidian) {
    double d = euclidian_distance(free_location, fe.base_location);
    double ttravel = d / fe.speed;
    Location location = free_location;
    if (t < t0 + ttravel) {
      double x1, x2, y1, y2;
      std::tie(x1, x2) = free_location;
      std::tie(y1, y2) = fe.base_location;
      Location direction;
      direction.first = (y1 - x1) / d;
      direction.second = (y2 - x2) / d;
      location.first = x1 + fe.speed * direction.first * (t - t0);
      location.second = x2 + fe.speed * direction.second * (t - t0);
      return location;
    } else {
      return fe.base_location;
    }
  } else {
    const double R = 6371;
    const double radian = M_PI / 180;
    const double v = 60;

    double p11 = R * cos(radian * free_location.first) *
                 cos(radian * free_location.second);
    double p12 = R * cos(radian * free_location.first) *
                 sin(radian * free_location.second);
    double p13 = R * sin(radian * free_location.first);

    double p21 = R * cos(radian * fe.base_location.first) *
                 cos(radian * fe.base_location.second);
    double p22 = R * cos(radian * fe.base_location.first) *
                 sin(radian * fe.base_location.second);
    double p23 = R * sin(radian * fe.base_location.first);

    double d = sqrt(pow(p11 - p21, 2) + pow(p12 - p22, 2) + pow(p13 - p23, 2));
    double alpha = 2 * asin(d / (2 * R));
    double dearth = R * alpha;
    double ttravel = dearth / v;
    double t0 = fe.arrival_time_at_f_last_trip;
    Location location = fe.free_location;
    if (t < t0 + ttravel) {
      double alpha0 = (t - t0) * v / R;
      double beta = sin(alpha - alpha0) / sin(alpha);
      double gamma =
          cos(alpha - alpha0) - sin(alpha - alpha0) * cos(alpha) / sin(alpha);
      double curr_p1 = beta * p11 + gamma * p21;
      double curr_p2 = beta * p12 + gamma * p22;
      double curr_p3 = beta * p13 + gamma * p23;
      location.first = (asin(curr_p3 / R) / M_PI) * 180;
      if (curr_p2 > -g_params.EPS) {
        location.second =
            (acos(curr_p1 / sqrt(pow(R, 2) - pow(curr_p3, 2))) / M_PI) * 180;
      } else {
        location.second =
            -(acos(curr_p1 / sqrt(pow(R, 2) - pow(curr_p3, 2))) / M_PI) * 180;
      }
      return location;
    } else {
      return fe.base_location;
    }
  }
}

double Travel::travel_time_from_position(Ambulance &amb, Call &call) {
  if (euclidian) {
    double d = euclidian_distance(amb.free_location, amb.base_location);
    double speed = amb.speed;
    double t = call.time;
    double ttravel = d / speed;
    double t0 = amb.arrival_time_at_f_last_trip;
    Location location = amb.free_location;
    if (t < t0 + ttravel) {
      double x1, x2, y1, y2;
      std::tie(x1, x2) = amb.free_location;
      std::tie(y1, y2) = amb.base_location;
      Location direction;
      direction.first = (y1 - x1) / d;
      direction.second = (y2 - x2) / d;
      location.first = x1 + amb.speed * direction.first * (t - t0);
      location.second = x2 + amb.speed * direction.second * (t - t0);
      return euclidian_travel_time(location, call.location, speed);
    } else {
      return euclidian_travel_time(amb.base_location, call.location, speed);
    }
  } else {
    Location location = ambulance_position(amb, call.time);
    return street_travel_time(location, call.location, 3.6 * amb.speed);
  }
}

double Travel::region_service_time(Location &base, Location &r,
                                   Location &hospital, double speed) {
  const double v = 3.6 * speed;
  double tbr = lat_long_travel_time(base, r, v);
  double trh = lat_long_travel_time(r, hospital, v);
  double thb = lat_long_travel_time(hospital, base, v);

  return tbr + trh + thb;
}

Location Travel::position_between_origin_destination(Location &a, Location &b,
                                                     double t0, double t,
                                                     double tf,
                                                     Ambulance &amb) {
  if (euclidian) {
    double d = euclidian_distance(a, b);
    double ttravel = d / amb.speed;
    Location location = a;
    if (t < t0 + ttravel) {
      double x1, x2, y1, y2;
      std::tie(x1, x2) = a;
      std::tie(y1, y2) = b;
      Location direction;
      direction.first = (y1 - x1) / d;
      direction.second = (y2 - x2) / d;
      location.first = x1 + amb.speed * direction.first * (t - t0);
      location.second = x2 + amb.speed * direction.second * (t - t0);
      return location;
    } else {
      return b;
    }
  } else {
    const double R = 6371;
    const double radian = M_PI / 180;
    const double v = 3.6 * amb.speed;

    if (abs(t - t0) < g_params.EPS) {
      return a;
    }

    double total_distance = lat_long_distance(a, b);
    double current_total_time = tf - t0;
    double reduced_v = total_distance / current_total_time;

    double p11 = R * cos(radian * a.first) * cos(radian * a.second);
    double p12 = R * cos(radian * a.first) * sin(radian * a.second);
    double p13 = R * sin(radian * a.first);

    double p21 = R * cos(radian * b.first) * cos(radian * b.second);
    double p22 = R * cos(radian * b.first) * sin(radian * b.second);
    double p23 = R * sin(radian * b.first);

    double d = sqrt(pow(p11 - p21, 2) + pow(p12 - p22, 2) + pow(p13 - p23, 2));
    double alpha = 2 * asin(d / (2 * R));
    double dearth = R * alpha;

    double ttravel = dearth / reduced_v;
    if (abs(current_total_time - ttravel) > 1) {
      fmt::print("WEIRD: current_total_time = {} != ttravel {}\n",
                 current_total_time, ttravel);
      std::cin.get();
    }
    Location location = a;
    if (t < t0 + ttravel) {
      double alpha0 = (t - t0) * reduced_v / R;
      double beta = sin(alpha - alpha0) / sin(alpha);
      double gamma =
          cos(alpha - alpha0) - sin(alpha - alpha0) * cos(alpha) / sin(alpha);
      double curr_p1 = beta * p11 + gamma * p21;
      double curr_p2 = beta * p12 + gamma * p22;
      double curr_p3 = beta * p13 + gamma * p23;
      location.first = (asin(curr_p3 / R) / M_PI) * 180;
      if (curr_p2 > -g_params.EPS) {
        location.second =
            (acos(curr_p1 / sqrt(pow(R, 2) - pow(curr_p3, 2))) / M_PI) * 180;
      } else {
        location.second =
            -(acos(curr_p1 / sqrt(pow(R, 2) - pow(curr_p3, 2))) / M_PI) * 180;
      }
      return location;
    } else {
      return b;
    }
  }
}

Location Travel::position_between_origin_destination(Location &a, Location &b,
                                                     double t0, double t,
                                                     double tf,
                                                     FastEquipment &fe) {
  if (euclidian) {
    double d = euclidian_distance(a, b);
    double speed = fe.speed;
    double ttravel = euclidian_travel_time(a, b, speed);
    Location location = a;
    if (t < t0 + ttravel) {
      double x1, x2, y1, y2;
      std::tie(x1, x2) = a;
      std::tie(y1, y2) = b;
      Location direction;
      direction.first = (y1 - x1) / d;
      direction.second = (y2 - x2) / d;
      location.first = x1 + fe.speed * direction.first * (t - t0);
      location.second = x2 + fe.speed * direction.second * (t - t0);
      return location;
    } else {
      return b;
    }
  } else {
    const double R = 6371;
    const double radian = M_PI / 180;
    const double v = 3.6 * fe.speed;

    if (abs(t - t0) < g_params.EPS) {
      return a;
    }

    double total_distance = lat_long_distance(a, b);
    double current_total_time = tf - t0;
    double reduced_v = total_distance / current_total_time;
    double p11 = R * cos(radian * a.first) * cos(radian * a.second);
    double p12 = R * cos(radian * a.first) * sin(radian * a.second);
    double p13 = R * sin(radian * a.first);

    double p21 = R * cos(radian * b.first) * cos(radian * b.second);
    double p22 = R * cos(radian * b.first) * sin(radian * b.second);
    double p23 = R * sin(radian * b.first);

    double d = sqrt(pow(p11 - p21, 2) + pow(p12 - p22, 2) + pow(p13 - p23, 2));
    double alpha = 2 * asin(d / (2 * R));
    double dearth = R * alpha;
    double ttravel = dearth / reduced_v;
    if (abs(current_total_time - ttravel) > 1) {
      fmt::print("WEIRD: current_total_time = {} != ttravel {}\n",
                 current_total_time, ttravel);
      std::cin.get();
    }
    Location location = a;
    if (t < t0 + ttravel) {
      double alpha0 = (t - t0) * reduced_v / R;
      double beta = sin(alpha - alpha0) / sin(alpha);
      double gamma =
          cos(alpha - alpha0) - sin(alpha - alpha0) * cos(alpha) / sin(alpha);
      double curr_p1 = beta * p11 + gamma * p21;
      double curr_p2 = beta * p12 + gamma * p22;
      double curr_p3 = beta * p13 + gamma * p23;
      location.first = (asin(curr_p3 / R) / M_PI) * 180;
      if (curr_p2 > -g_params.EPS) {
        location.second =
            (acos(curr_p1 / sqrt(pow(R, 2) - pow(curr_p3, 2))) / M_PI) * 180;
      } else {
        location.second =
            -(acos(curr_p1 / sqrt(pow(R, 2) - pow(curr_p3, 2))) / M_PI) * 180;
      }
      return location;
    } else {
      return b;
    }
  }
}

double Travel::travel_time(Location &a, Location &b, Ambulance &amb) {
  if (euclidian) {
    return euclidian_travel_time(a, b, amb.speed);
  } else {
    return street_travel_time(a, b, 3.6 * amb.speed);
  }
}

double Travel::travel_time(Location &a, Location &b, FastEquipment &fe) {
  if (euclidian) {
    return euclidian_travel_time(a, b, fe.speed);
  } else if (fe.type == FastType::DRONE) {
    return lat_long_travel_time(a, b, 3.6 * fe.speed);
  } else if (fe.type == FastType::MOTORCYCLE) {
    return street_travel_time(a, b, 3.6 * fe.speed);
  } else {
    return GRB_INFINITY;
  }
}

double Travel::euclidian_travel_time(Location &a, Location &b, double speed,
                                     bool free_flow) {
  return euclidian_distance(a, b) / speed;
}

double Travel::lat_long_travel_time(Location &a, Location &b, double v,
                                    bool free_flow) {
  double free_flow_travel_time = 3600 * (lat_long_distance(a, b) / v);
  if (free_flow) {
    // fmt::print("free flow\n");
    return free_flow_travel_time;
  }
  double delay =
      std::min(free_flow_travel_time, get_delay(free_flow_travel_time, v));
  // fmt::print("Free_flow = {}, Delay = {}, param {}\n", free_flow_travel_time,
  // delay, delay_param); return free_flow_travel_time +
  // std::max(free_flow_travel_time, get_delay(free_flow_travel_time));
  return free_flow_travel_time + delay;
}

double Travel::street_travel_time(Location &a, Location &b, double v,
                                  bool free_flow) {
  // TODO: Change for the streetmap time approach
  return lat_long_travel_time(a, b, v);
}

double Travel::euclidian_distance(Location &a, Location &b) {
  double x1, y1, x2, y2;
  std::tie(x1, y1) = a;
  std::tie(x2, y2) = b;
  return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

double Travel::lat_long_distance(Location &a, Location &b) {
  const double R = 6371;
  const double radian = M_PI / 180;

  double p11 = R * cos(radian * a.first) * cos(radian * a.second);
  double p12 = R * cos(radian * a.first) * sin(radian * a.second);
  double p13 = R * sin(radian * a.first);

  double p21 = R * cos(radian * b.first) * cos(radian * b.second);
  double p22 = R * cos(radian * b.first) * sin(radian * b.second);
  double p23 = R * sin(radian * b.first);

  double d = sqrt(pow(p11 - p21, 2) + pow(p12 - p22, 2) + pow(p13 - p23, 2));
  return 2 * R * asin(d / (2 * R));
}

double Travel::street_distance(Location &a, Location &b) {
  // TODO: Change to streetmap distances
  return lat_long_distance(a, b);
}

std::vector<std::vector<double>> Travel::table_in_out(
    std::vector<Location> &a, std::vector<Location> &b) {
  std::vector<std::vector<double>> result(
      a.size(), std::vector<double>(b.size(), GRB_INFINITY));

  for (size_t i = 0; i < a.size(); ++i) {
    for (size_t j = 0; j < b.size(); ++j) {
      result[i][j] = street_distance(a[i], b[j]);
    }
  }

  return result;
}
