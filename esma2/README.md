# ESMA

Emergency service simulation and optimization codebase.

**Requirements**

1. CMake 3.22+
2. C++17 compiler
3. Boost (components: unit_test_framework, program_options, filesystem, thread, system, iostreams, chrono, date_time, regex)
4. Gurobi (set `GUROBI_HOME`)
5. fmt, spdlog, cpptrace
6. OSRM library (either `/usr/local/lib/libosrm.a` + headers or an `osrm` CMake config)

**Build**

1. Set `GUROBI_HOME` if needed (default is `/opt/gurobi1203/linux64`):

```bash
export GUROBI_HOME=/path/to/gurobi/linux64
```

2. Configure and build:

```bash
cmake -S . -B build \
  -DGUROBI_HOME="${GUROBI_HOME}" \
  -DGUROBI_LIBRARY="${GUROBI_HOME}/lib/libgurobi120.so" \
  -DGUROBI_CXX_LIBRARY="${GUROBI_HOME}/src/build/libgurobi_c++.a"
cmake --build build
```

This generates the binary at `build/esma`.

Optional CMake parameters:

- `ESMA_ENABLE_TRAJECTORIES` (ON/OFF) enables writing of the ambulance trajectories in the simulation
- `ESMA_TRAVEL_STREET` (ON/OFF) enables street travel times/distances.
- `ESMA_TRAVEL_CARTESIAN` (ON/OFF) enables cartesian travel times/distances.

Street travel:

```bash
cmake -S . -B build \
  -DESMA_TRAVEL_STREET=ON \
  -DESMA_TRAVEL_CARTESIAN=OFF
cmake --build build
```

Geodesic travel (default):

```bash
cmake -S . -B build \
  -DESMA_TRAVEL_STREET=OFF \
  -DESMA_TRAVEL_CARTESIAN=OFF
cmake --build build
```

Cartesian travel:

```bash
cmake -S . -B build \
  -DESMA_TRAVEL_STREET=OFF \
  -DESMA_TRAVEL_CARTESIAN=ON
cmake --build build
```

Enable trajectories:

```bash
cmake -S . -B build -DESMA_ENABLE_TRAJECTORIES=ON
cmake --build build
```

Disable trajectories:

```bash
cmake -S . -B build -DESMA_ENABLE_TRAJECTORIES=OFF
cmake --build build
```

**Run**

1. From the repo root, run with the default config:

```bash
./build/esma -f ../default.cfg
```

2. Show available CLI options:

```bash
./build/esma --help
```

Config examples in the repo:

- `default.cfg`
- `default_hex.cfg`

**Tests**

```bash
cmake --build build --target run_tests
```
