-- Ambulance profile (based on car.lua, tuned for emergency response)

api_version = 4

Set = require('lib/set')
Sequence = require('lib/sequence')
WayHandlers = require("lib/way_handlers")
Relations = require("lib/relations")
TrafficSignal = require("lib/traffic_signal")
find_access_tag = require("lib/access").find_access_tag
limit = require("lib/maxspeed").limit
Utils = require("lib/utils")
Measure = require("lib/measure")

function setup()
  return {
    properties = {
      max_speed_for_map_matching      = 200/3.6, -- allow higher match speed for responses
      weight_name                     = 'duration',
      process_call_tagless_node      = false,
      u_turn_penalty                 = 5,   -- allow tighter rerouting
      continue_straight_at_waypoint  = true,
      use_turn_restrictions          = false, -- emergency vehicles can ignore
      left_hand_driving              = false,
      traffic_light_penalty          = 0,   -- assume lights can be bypassed
    },

    default_mode              = mode.driving,
    default_speed             = 15,
    oneway_handling           = false, -- allow contraflow when needed
    side_road_multiplier      = 1.0,
    turn_penalty              = 5,
    speed_reduction           = 0.9,
    turn_bias                 = 1.0,
    cardinal_directions       = false,

    vehicle_height = 2.6,
    vehicle_width = 2.2,
    vehicle_length = 6.0,
    vehicle_weight = 3500,

    suffix_list = {
      'N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'North', 'South', 'West', 'East', 'Nor', 'Sou', 'We', 'Ea'
    },

    barrier_whitelist = Set {
      'cattle_grid',
      'border_control',
      'toll_booth',
      'sally_port',
      'gate',
      'lift_gate',
      'no',
      'entrance',
      'height_restrictor',
      'arch',
      'emergency_access'
    },

    access_tag_whitelist = Set {
      'yes',
      'motorcar',
      'motor_vehicle',
      'vehicle',
      'permissive',
      'designated',
      'hov',
      'emergency'
    },

    access_tag_blacklist = Set {
      'no',
      'agricultural',
      'forestry',
      'psv',
      'customers',
      'private',
      'delivery',
      'destination'
    },

    service_access_tag_blacklist = Set {
        'private'
    },

    restricted_access_tag_list = Set {
      'private',
      'delivery',
      'destination',
      'customers',
    },

    access_tags_hierarchy = Sequence {
      'motorcar',
      'motor_vehicle',
      'vehicle',
      'access'
    },

    service_tag_forbidden = Set {
      -- keep empty so emergency lanes remain usable
    },

    restrictions = Sequence {
      'motorcar',
      'motor_vehicle',
      'vehicle'
    },

    classes = Sequence {
        'toll', 'motorway', 'ferry', 'restricted', 'tunnel'
    },

    excludable = Sequence {
        Set {'toll'},
        Set {'motorway'},
        Set {'ferry'}
    },

    avoid = Set {
      'area',
      'reversible',
      'impassable',
      'steps',
      'construction',
      'proposed'
    },

    speeds = Sequence {
      highway = {
        motorway        = 110,
        motorway_link   = 70,
        trunk           = 100,
        trunk_link      = 65,
        primary         = 85,
        primary_link    = 60,
        secondary       = 75,
        secondary_link  = 55,
        tertiary        = 60,
        tertiary_link   = 45,
        unclassified    = 45,
        residential     = 45,
        living_street   = 25,
        service         = 30
      }
    },

    service_penalties = {
      alley             = 0.7,
      parking           = 0.7,
      parking_aisle     = 0.7,
      driveway          = 0.7,
      ["drive-through"] = 0.7,
      ["drive-thru"] = 0.7
    },

    restricted_highway_whitelist = Set {
      'motorway',
      'motorway_link',
      'trunk',
      'trunk_link',
      'primary',
      'primary_link',
      'secondary',
      'secondary_link',
      'tertiary',
      'tertiary_link',
      'residential',
      'living_street',
      'unclassified',
      'service'
    },

    construction_whitelist = Set {
      'no',
      'widening',
      'minor',
    },

    route_speeds = {
      ferry = 8,
      shuttle_train = 15
    },

    bridge_speeds = {
      movable = 8
    },

    surface_speeds = {
      asphalt = nil,
      concrete = nil,
      ["concrete:plates"] = nil,
      ["concrete:lanes"] = nil,
      paved = nil,

      cement = 90,
      compacted = 90,
      fine_gravel = 80,

      paving_stones = 70,
      metal = 60,
      bricks = 60,

      grass = 45,
      wood = 40,
      sett = 45,
      grass_paver = 45,
      gravel = 45,
      unpaved = 45,
      ground = 40,
      dirt = 40,
      pebblestone = 40,
      tartan = 40,

      cobblestone = 35,
      clay = 35,

      earth = 25,
      stone = 25,
      rocky = 25,
      sand = 20,

      mud = 10
    },

    tracktype_speeds = {
      grade1 =  70,
      grade2 =  55,
      grade3 =  40,
      grade4 =  30,
      grade5 =  25
    },

    smoothness_speeds = {
      intermediate    =  90,
      bad             =  50,
      very_bad        =  30,
      horrible        =  15,
      very_horrible   =  10,
      impassable      =  0
    },

    maxspeed_table_default = {
      urban = 60,
      rural = 100,
      trunk = 120,
      motorway = 140
    },

    maxspeed_table = {
      ["none"] = 160
    },

    relation_types = Sequence {
      "route"
    },

    highway_turn_classification = {
    },

    access_turn_classification = {
    }
  }
end

function process_node(profile, node, result, relations)
  local access = find_access_tag(node, profile.access_tags_hierarchy)
  if access then
    if profile.access_tag_blacklist[access] and not profile.restricted_access_tag_list[access] then
      result.barrier = true
    end
  else
    local barrier = node:get_value_by_key("barrier")
    if barrier then
      local restricted_by_height = false
      if barrier == 'height_restrictor' then
         local maxheight = Measure.get_max_height(node:get_value_by_key("maxheight"), node)
         restricted_by_height = maxheight and maxheight < profile.vehicle_height
      end

      local bollard = node:get_value_by_key("bollard")
      local rising_bollard = bollard and "rising" == bollard

      local kerb = node:get_value_by_key("kerb")
      local highway = node:get_value_by_key("highway")
      local flat_kerb = kerb and ("lowered" == kerb or "flush" == kerb)
      local highway_crossing_kerb = barrier == "kerb" and highway and highway == "crossing"

      if not profile.barrier_whitelist[barrier]
                and not rising_bollard
                and not flat_kerb
                and not highway_crossing_kerb
                or restricted_by_height then
        result.barrier = true
      end
    end
  end

  result.traffic_lights = TrafficSignal.get_value(node)
end

function process_way(profile, way, result, relations)
  local data = {
    highway = way:get_value_by_key('highway'),
    bridge = way:get_value_by_key('bridge'),
    route = way:get_value_by_key('route')
  }

  if (not data.highway or data.highway == '') and
  (not data.route or data.route == '')
  then
    return
  end

  local handlers = Sequence {
    WayHandlers.default_mode,
    WayHandlers.blocked_ways,
    WayHandlers.avoid_ways,
    WayHandlers.handle_height,
    WayHandlers.handle_width,
    WayHandlers.handle_length,
    WayHandlers.handle_weight,
    WayHandlers.access,
    WayHandlers.oneway,
    WayHandlers.destinations,
    WayHandlers.ferries,
    WayHandlers.movables,
    WayHandlers.service,
    WayHandlers.hov,
    WayHandlers.speed,
    WayHandlers.maxspeed,
    WayHandlers.surface,
    WayHandlers.penalties,
    WayHandlers.classes,
    WayHandlers.turn_lanes,
    WayHandlers.classification,
    WayHandlers.roundabouts,
    WayHandlers.startpoint,
    WayHandlers.driving_side,
    WayHandlers.names,
    WayHandlers.weights,
    WayHandlers.way_classification_for_turn
  }

  WayHandlers.run(profile, way, result, data, handlers, relations)

  if profile.cardinal_directions then
      Relations.process_way_refs(way, relations, result)
  end
end

function process_turn(profile, turn)
  local turn_penalty = profile.turn_penalty
  local turn_bias = turn.is_left_hand_driving and 1. / profile.turn_bias or profile.turn_bias

  if turn.has_traffic_light then
      turn.duration = profile.properties.traffic_light_penalty
  end

  if turn.number_of_roads > 2 or turn.source_mode ~= turn.target_mode or turn.is_u_turn then
    if turn.angle >= 0 then
      turn.duration = turn.duration + turn_penalty / (1 + math.exp( -((13 / turn_bias) *  turn.angle/180 - 6.5*turn_bias)))
    else
      turn.duration = turn.duration + turn_penalty / (1 + math.exp( -((13 * turn_bias) * -turn.angle/180 - 6.5/turn_bias)))
    end

    if turn.is_u_turn then
      turn.duration = turn.duration + profile.properties.u_turn_penalty
    end
  end

  if profile.properties.weight_name == 'distance' then
     turn.weight = 0
  else
     turn.weight = turn.duration
  end

  if profile.properties.weight_name == 'routability' then
      if not turn.source_restricted and turn.target_restricted then
          turn.weight = constants.max_turn_weight
      end
  end
end

return {
  setup = setup,
  process_node = process_node,
  process_way = process_way,
  process_turn = process_turn
}
