#!/usr/bin/env bash
set -euo pipefail

# Configuration
GUROBI_HOME=${GUROBI_HOME:-/opt/gurobi1203/linux64}
BUILD_DIR=${BUILD_DIR:-build}
GENERATOR=${GENERATOR:-"Unix Makefiles"}

echo "Using GUROBI_HOME=${GUROBI_HOME}"
echo "Build dir: ${BUILD_DIR}"

cmake -S . -B "${BUILD_DIR}" -G "${GENERATOR}" \
  -DGUROBI_HOME="${GUROBI_HOME}" \
  -DGUROBI_LIBRARY="${GUROBI_HOME}/lib/libgurobi120.so" \
  -DGUROBI_CXX_LIBRARY="${GUROBI_HOME}/src/build/libgurobi_c++.a" \
  -DESMA_ENABLE_TRAJECTORIES=0

cmake --build "${BUILD_DIR}"
