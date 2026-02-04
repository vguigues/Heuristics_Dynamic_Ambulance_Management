#define BOOST_TEST_MODULE ESMA_TEST
#include <fmt/core.h>
#include <fmt/ranges.h>

#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <random>
#include <osrm/engine_config.hpp>
#include <osrm/json_container.hpp>
#include <osrm/osrm.hpp>
#include <osrm/route_parameters.hpp>

#include "esma/travel.hpp"

// Tolerância para comparações de ponto flutuante
const double tol = 1e-3;

using esma::Location;

namespace {
bool has_osrm_prefix_files(const std::filesystem::path& osrm_path) {
  const std::filesystem::path dir = osrm_path.parent_path();
  const std::string prefix = osrm_path.filename().string() + ".";
  if (!std::filesystem::exists(dir) || !std::filesystem::is_directory(dir)) {
    return false;
  }
  for (const auto& entry : std::filesystem::directory_iterator(dir)) {
    const auto name = entry.path().filename().string();
    if (name.rfind(prefix, 0) == 0) {
      return true;
    }
  }
  return false;
}

std::string resolve_osrm_path() {
  if (const char* env = std::getenv("ESMA_OSRM_DATA")) {
    return std::string(env);
  }
  const std::string candidates[] = {
      "data/rio_osrm/rio.osrm",
      "../data/rio_osrm/rio.osrm",
      "../../data/rio_osrm/rio.osrm",
  };
  for (const auto& candidate : candidates) {
    if (std::filesystem::exists(candidate) ||
        has_osrm_prefix_files(candidate)) {
      return candidate;
    }
  }
  return candidates[0];
}

void print_route_geometry(const Location& origin, const Location& destination) {
  std::string osrm_path = resolve_osrm_path();
  if (!std::filesystem::exists(osrm_path) &&
      !has_osrm_prefix_files(osrm_path)) {
    fmt::print("OSRM data not found near {}. Skipping geometry print.\n",
               osrm_path);
    return;
  }

  osrm::EngineConfig config;
  config.storage_config = {osrm_path};
  config.use_mmap = false;
  config.use_shared_memory = false;
  config.algorithm = osrm::EngineConfig::Algorithm::MLD;
  osrm::OSRM osrm_engine(config);

  osrm::RouteParameters params;
  params.coordinates.emplace_back(osrm::util::FloatLongitude{origin.longitude},
                                  osrm::util::FloatLatitude{origin.latitude});
  params.coordinates.emplace_back(
      osrm::util::FloatLongitude{destination.longitude},
      osrm::util::FloatLatitude{destination.latitude});
  params.geometries = osrm::RouteParameters::GeometriesType::GeoJSON;

  osrm::json::Object result;
  auto status = osrm_engine.Route(params, result);
  if (status != osrm::Status::Ok) {
    fmt::print("OSRM route failed; no geometry to print.\n");
    return;
  }

  const auto& routes =
      result.values.at("routes").get<osrm::json::Array>().values;
  if (routes.empty()) {
    fmt::print("OSRM route result empty.\n");
    return;
  }

  const auto& route = routes.front().get<osrm::json::Object>();
  const auto& geometry = route.values.at("geometry").get<osrm::json::Object>();
  const auto& coords =
      geometry.values.at("coordinates").get<osrm::json::Array>().values;

  fmt::print("Route geometry points: {}\n", coords.size());
  const size_t max_print = 29;
  for (size_t i = 0; i < coords.size() && i < max_print; ++i) {
    const auto& pair = coords[i].get<osrm::json::Array>().values;
    double lon = pair[0].get<osrm::json::Number>().value;
    double lat = pair[1].get<osrm::json::Number>().value;
    fmt::print("{}, {}, {}\n", i, lat, lon);
  }
  if (coords.size() > max_print) {
    fmt::print("... ({} more points)\n", coords.size() - max_print);
  }
}
}  // namespace

BOOST_AUTO_TEST_CASE(TestTravelDistance) {
  Location source = {-23.5505, -46.6333};  // São Paulo
  Location target = {-22.9068, -43.1729};  // Rio de Janeiro

  bu_length_t d = esma::travel_geodesic_distance(source, target);

  // Distância geodésica conhecida ~357 km
  BOOST_CHECK_CLOSE(d.value(), 360000.0, 1.0);  // tolerância 1%

  Location loc3 = {-22.907353412207307, -43.170207064908496};

  bu_length_t d2 = esma::travel_geodesic_distance(target, loc3);
  BOOST_CHECK_CLOSE(d2.value(), 283, 1.0);
}

BOOST_AUTO_TEST_CASE(TestTravelTime) {
  Location source = {-23.5505, -46.6333};
  Location target = {-22.9068, -43.1729};

  bu_velocity_t v = 25.0 * si::meters_per_second;  // 90 km/h

  bu_time_t t = esma::travel_geodesic_time(source, target, v);

  // Tempo esperado: distância / velocidade
  double expected_time =
      esma::travel_geodesic_distance(source, target).value() / 25.0;

  BOOST_CHECK_CLOSE(t.value(), expected_time, tol * 100);
}

BOOST_AUTO_TEST_CASE(TestTravelCurrentPosition) {
  Location source = {-23.5505, -46.6333};
  Location target = {-22.9068, -43.1729};

  bu_velocity_t v = 25.0 * si::meters_per_second;

  bu_time_t t0 = 0.0 * si::seconds;
  bu_time_t t = 1800.0 * si::seconds;  // meia hora

  Location pos =
      esma::travel_geodesic_current_position(source, target, v, t0, t);

  // A posição deve estar entre source e target
  BOOST_CHECK(pos.latitude > std::min(source.latitude, target.latitude) - tol);
  BOOST_CHECK(pos.latitude < std::max(source.latitude, target.latitude) + tol);
  BOOST_CHECK(pos.longitude >
              std::min(source.longitude, target.longitude) - tol);
  BOOST_CHECK(pos.longitude <
              std::max(source.longitude, target.longitude) + tol);
}

BOOST_AUTO_TEST_CASE(TestTravelRandomTime) {
  Location source = {-23.5505, -46.6333};
  Location target = {-22.9068, -43.1729};

  bu_velocity_t v = 25.0 * si::meters_per_second;

  unsigned int seed = 42;

  double t = esma::travel_geodesic_time(source, target, v).value();
  double random_t =
      esma::travel_geodesic_random_time(source, target, v, seed).value();

  fmt::print("t = {}, random_t = {}\n", t, random_t);
  // O tempo aleatório não pode ser negativo
  BOOST_CHECK(random_t > 0.0);

  // Deve ser próximo do tempo determinístico, dentro de algumas vezes a stddev
  double deterministic = esma::travel_geodesic_time(source, target, v).value();
  BOOST_CHECK(random_t > deterministic * 0.5);
  BOOST_CHECK(random_t < deterministic * 2.0);
}

BOOST_AUTO_TEST_CASE(TestTravelRandomCurrentPosition) {
  Location source = {-23.5505, -46.6333};
  Location target = {-22.9068, -43.1729};

  bu_velocity_t v = 25.0 * si::meters_per_second;

  bu_time_t t0 = 0.0 * si::seconds;
  bu_time_t t = 3600.0 * si::seconds;  // 1 hora

  unsigned int seed = 123;

  Location pos = esma::travel_geodesic_random_current_position(source, target,
                                                               v, t0, t, seed);

  // A posição deve ser plausível, entre source e target
  BOOST_CHECK(pos.latitude > std::min(source.latitude, target.latitude) - tol);
  BOOST_CHECK(pos.latitude < std::max(source.latitude, target.latitude) + tol);
  BOOST_CHECK(pos.longitude >
              std::min(source.longitude, target.longitude) - tol);
  BOOST_CHECK(pos.longitude <
              std::max(source.longitude, target.longitude) + tol);
}

BOOST_AUTO_TEST_CASE(TestStreetTravelDistanceAndTime) {
  Location origin = {-22.9068, -43.1729};       // Rio de Janeiro
  Location destination = {-22.9676, -43.2007};  // Copacabana -> Botafogo area
  fmt::print("Street travel test: origin=({}, {}), destination=({}, {})\n",
             origin.latitude, origin.longitude, destination.latitude,
             destination.longitude);

  bu_velocity_t v = 15.0 * si::meters_per_second;  // ~54 km/h

  bu_length_t d = esma::travel_street_distance(origin, destination);
  bu_time_t t = esma::travel_street_time(origin, destination, v);

  fmt::print("Street travel distance = {} m, time = {} s\n", d.value(),
             t.value());

  BOOST_CHECK(d.value() > 0.0);
  BOOST_CHECK(t.value() > 0.0);
  // Distance should be reasonable (< 1e6 m even if falling back to geodesic)
  BOOST_CHECK(d.value() < 1e6);

  print_route_geometry(origin, destination);
}

BOOST_AUTO_TEST_CASE(TestStreetTravelCurrentPosition) {
  Location origin = {-22.9068, -43.1729};       // Rio de Janeiro
  Location destination = {-22.9676, -43.2007};  // Copacabana -> Botafogo area
  bu_velocity_t v = 16.667 * si::meters_per_second;  // ~60 km/h

  esma::StreetTravel travel(true);
  bu_time_t t0 = 0.0 * si::seconds;

  bu_time_t total_time = travel.travel_time(origin, destination, v);
  BOOST_CHECK(total_time.value() > 0.0);

  Location pos_at_t0 = travel.current_position(origin, destination, v, t0, t0);
  BOOST_CHECK(pos_at_t0 == origin);

  bu_time_t t_before = -10.0 * si::seconds;
  Location pos_before =
      travel.current_position(origin, destination, v, t0, t_before);
  BOOST_CHECK(pos_before == origin);

  bu_time_t t_mid = total_time * 0.5;
  Location pos_mid = travel.current_position(origin, destination, v, t0, t_mid);
  BOOST_CHECK(pos_mid.latitude >= -90.0 && pos_mid.latitude <= 90.0);
  BOOST_CHECK(pos_mid.longitude >= -180.0 && pos_mid.longitude <= 180.0);
  fmt::print("Midpoint position: ({}, {}) at time {}/{}\n", pos_mid.latitude,
             pos_mid.longitude, t_mid.value(), total_time.value());

  bu_time_t t_after = total_time + 1.0 * si::seconds;
  Location pos_after =
      travel.current_position(origin, destination, v, t0, t_after);
  double dist_to_target =
      esma::travel_geodesic_distance(pos_after, destination).value();
  BOOST_CHECK(dist_to_target < 1000.0);
}

BOOST_AUTO_TEST_CASE(TestStreetTravelCoverageRio) {
  const std::string osrm_path = resolve_osrm_path();
  if (!std::filesystem::exists(osrm_path) &&
      !has_osrm_prefix_files(osrm_path)) {
    fmt::print("OSRM data not found near {}. Skipping coverage test.\n",
               osrm_path);
    return;
  }

  const double min_lat = -23.1;
  const double max_lat = -22.8;
  const double min_lon = -43.8;
  const double max_lon = -43.1;
  const int sample_count = 200;
  const double max_reasonable_m = 200000.0;

  std::mt19937 rng(123);
  std::uniform_real_distribution<double> lat_dist(min_lat, max_lat);
  std::uniform_real_distribution<double> lon_dist(min_lon, max_lon);

  esma::StreetTravel travel(true);
  int invalid_count = 0;
  int nan_count = 0;
  int zero_count = 0;
  int too_big_count = 0;

  for (int i = 0; i < sample_count; ++i) {
    Location from{lat_dist(rng), lon_dist(rng)};
    Location to{lat_dist(rng), lon_dist(rng)};
    bu_length_t dist = travel.travel_distance(from, to);
    double d = dist.value();
    bool is_nan = std::isnan(d);
    bool is_zero = d <= 0.0;
    bool is_too_big = d > max_reasonable_m;
    if (is_nan || is_zero || is_too_big) {
      ++invalid_count;
      nan_count += is_nan ? 1 : 0;
      zero_count += is_zero ? 1 : 0;
      too_big_count += is_too_big ? 1 : 0;
    }
  }

  fmt::print(
      "OSRM coverage test (Rio bbox): {} samples, invalid {}, nan {}, zero "
      "{}, too_big {}\n",
      sample_count, invalid_count, nan_count, zero_count, too_big_count);

  BOOST_CHECK(invalid_count < sample_count);
}
