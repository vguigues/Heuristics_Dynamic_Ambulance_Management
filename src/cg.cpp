#include "../include/cg.h"

int get_time(Data& data, int time) { return 1 + (int)time / data.quantum; }

CGCall::CGCall(Data& data, GRBEnv& env)
    : data(data),
      env(env),
      model(env),
      columns_added(0),
      type_call0(data.type_call0),
      local_call0(data.local_call0) {
  constexpr bool debug = false;
  if (debug) {
    std::cout << "DEBUG\n";
  }
  std::vector<std::pair<int, int>> intermediate;
  intermediate.reserve(data.total_locals);
  closest_l1_to_l.reserve(data.num_locals);
  Travel& travel = data.travel;
  // std::cout << "closest_l1_to_l:\n";
  for (int l = 0; l < data.num_locals; ++l) {
    closest_l1_to_l.push_back(std::vector<int>());
    closest_l1_to_l[l].reserve(data.total_locals);
    intermediate.clear();
    for (int l1 = 0; l1 < data.num_locals; ++l1) {
      auto dist_l1_l =
          travel.lat_long_distance(data.locations[l1], data.locations[l]);
      intermediate.push_back(std::make_pair(dist_l1_l, l1));
    }
    std::sort(intermediate.begin(), intermediate.end());
    // std::cout << l << ": ";
    for (auto& elem : intermediate) {
      // std::cout << elem.second << " ";
      closest_l1_to_l[l].push_back(elem.second);
    }
    // std::cout << "\n";
  }
  // std::cin.get();

  vh_tall = new int***[data.num_times];
  vc_tl = new int*[data.num_times];
  vb_tal = new int**[data.num_times];
  for (int t = 0; t < data.num_times; ++t) {
    vh_tall[t] = new int**[data.types_amb];
    vb_tal[t] = new int*[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      vh_tall[t][a] = new int*[data.num_locals];
      vb_tal[t][a] = new int[data.num_locals];
      for (int l1 = 0; l1 < data.num_locals; ++l1) {
        vh_tall[t][a][l1] = new int[data.num_locals];
        for (int l = 0; l < data.num_locals; ++l) {
          vh_tall[t][a][l1][l] = -1;
        }

        vb_tal[t][a][l1] = -1;
      }
    }

    vc_tl[t] = new int[data.num_locals];
    for (int l = 0; l < data.num_locals; ++l) {
      vc_tl[t][l] = -1;
    }
  }

  std::stringstream name;

  A0_ab = new int*[data.types_amb];
  A0_ah = new int*[data.types_amb];
  A0_alb = new int**[data.types_amb];
  A0_calh = new int***[data.types_call];
  A0_callh = new int****[data.types_call];
  C = new int*[data.types_call];

  for (int a = 0; a < data.types_amb; ++a) {
    A0_ab[a] = new int[data.num_bases];
    for (int b = 0; b < data.num_bases; ++b) {
      A0_ab[a][b] = data.A0_ab[a][b];
    }
    A0_ah[a] = new int[data.num_hosps];
    for (int h = 0; h < data.num_hosps; ++h) {
      A0_ah[a][h] = data.A0_ah[a][h];
    }
    A0_alb[a] = new int*[data.total_locals];
    for (int l1 = 0; l1 < data.total_locals; ++l1) {
      A0_alb[a][l1] = new int[data.num_bases];
      for (int b = 0; b < data.num_bases; ++b) {
        A0_alb[a][l1][b] = data.A0_alb[a][l1][b];
      }
    }
  }

  for (int c = 0; c < data.types_call; ++c) {
    A0_calh[c] = new int**[data.types_amb];
    A0_callh[c] = new int***[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      A0_calh[c][a] = new int*[data.total_locals];
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        A0_calh[c][a][l1] = new int[data.num_hosps];
        for (int h = 0; h < data.num_hosps; ++h) {
          A0_calh[c][a][l1][h] = data.A0_calh[c][a][l1][h];
        }
      }
      A0_callh[c][a] = new int**[data.total_locals];
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        A0_callh[c][a][l1] = new int*[data.num_locals];
        for (int l = 0; l < data.num_locals; ++l) {
          A0_callh[c][a][l1][l] = new int[data.num_hosps];
          for (int h = 0; h < data.num_hosps; ++h) {
            A0_callh[c][a][l1][l][h] = data.A0_callh[c][a][l1][l][h];
          }
        }
      }
    }
    C[c] = new int[data.num_locals];
    for (int l = 0; l < data.num_locals; ++l) {
      C[c][l] = data.C[c][l];
    }
  }

  x0_abh = new GRBVar**[data.types_amb];
  x0_ahh = new GRBVar**[data.types_amb];
  x0_albh = new GRBVar***[data.types_amb];
  y0_ahb = new GRBVar**[data.types_amb];

  for (auto a : data.A[type_call0]) {
    x0_abh[a] = new GRBVar*[data.num_bases];
    for (int b = 0; b < data.num_bases; ++b) {
      x0_abh[a][b] = new GRBVar[data.num_hosps];
      for (auto h : data.H[type_call0][local_call0]) {
        name << "x0_a" << a << "_b" << b << "_h" << h;
        x0_abh[a][b][h] =
            model.addVar(0, A0_ab[a][b], 0, GRB_CONTINUOUS, name.str());
        name.str("");
      }
    }

    x0_ahh[a] = new GRBVar*[data.num_hosps];
    for (int h1 = 0; h1 < data.num_hosps; ++h1) {
      x0_ahh[a][h1] = new GRBVar[data.num_hosps];
      for (auto h : data.H[type_call0][local_call0]) {
        name << "x0_a" << a << "_h'" << h1 << "_h" << h;
        x0_ahh[a][h1][h] =
            model.addVar(0, A0_ah[a][h1], 0, GRB_CONTINUOUS, name.str());
        name.str("");
      }
    }

    x0_albh[a] = new GRBVar**[data.total_locals];
    for (int l = 0; l < data.total_locals; ++l) {
      x0_albh[a][l] = new GRBVar*[data.num_bases];
      for (int b = 0; b < data.num_bases; ++b) {
        x0_albh[a][l][b] = new GRBVar[data.num_hosps];
        for (auto h : data.H[type_call0][local_call0]) {
          name << "x0_a" << a << "_l" << l << "_b" << b << "_h" << h;
          x0_albh[a][l][b][h] =
              model.addVar(0, A0_alb[a][l][b], 0, GRB_CONTINUOUS, name.str());
          name.str("");
        }
      }
    }

    y0_ahb[a] = new GRBVar*[data.num_hosps];
    for (int h = 0; h < data.num_hosps; ++h) {
      y0_ahb[a][h] = new GRBVar[data.num_bases];
      for (int b = 0; b < data.num_bases; ++b) {
        name << "y0" << "_a" << a << "_h" << h << "_b" << b;
        y0_ahb[a][h][b] =
            model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name.str());
        name.str("");
      }
    }
  }

  xt_cablh = new GRBVar*****[data.num_times - 1];
  xt_cahlh = new GRBVar*****[data.num_times - 1];
  xt_calblh = new GRBVar******[data.num_times - 1];
  yt_ahb = new GRBVar***[data.num_times - 1];
  Ct_cl = new GRBVar**[data.num_times];
  At_ab = new GRBVar**[data.num_times];
  At_alb = new GRBVar***[data.num_times];

  for (int t = 0; t < data.num_times - 1; ++t) {
    xt_cablh[t] = new GRBVar****[data.types_call];
    xt_cahlh[t] = new GRBVar****[data.types_call];
    xt_calblh[t] = new GRBVar*****[data.types_call];

    for (int c = 0; c < data.types_call; ++c) {
      xt_cablh[t][c] = new GRBVar***[data.types_amb];
      xt_cahlh[t][c] = new GRBVar***[data.types_amb];
      xt_calblh[t][c] = new GRBVar****[data.types_amb];
      for (int a = 0; a < data.types_amb; ++a) {
        xt_cablh[t][c][a] = new GRBVar**[data.num_bases];
        for (int b = 0; b < data.num_bases; ++b) {
          xt_cablh[t][c][a][b] = new GRBVar*[data.num_locals];
          for (int l = 0; l < data.num_locals; ++l) {
            xt_cablh[t][c][a][b][l] = new GRBVar[data.num_hosps];
            for (auto h : data.H[c][l]) {
              name << "xt" << t + 1 << "_c" << c << "_a" << a << "_b" << b;
              name << "_l" << l << "_h";
              name << h;
              xt_cablh[t][c][a][b][l][h] =
                  model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name.str());
              name.str("");
            }
          }
        }

        xt_cahlh[t][c][a] = new GRBVar**[data.num_hosps];
        for (int h1 = 0; h1 < data.num_hosps; ++h1) {
          xt_cahlh[t][c][a][h1] = new GRBVar*[data.num_locals];
          for (int l = 0; l < data.num_locals; ++l) {
            xt_cahlh[t][c][a][h1][l] = new GRBVar[data.num_hosps];
            for (int h = 0; h < data.num_hosps; ++h) {
              if (data.A[c].find(a) != data.A[c].end() &&
                  data.H[c][l].find(h) != data.H[c][l].end()) {
                name << "xt" << t + 1 << "_c" << c << "_a" << a << "_h'";
                name << h1 << "_l" << l << "_h" << h;
                xt_cahlh[t][c][a][h1][l][h] = model.addVar(
                    0, GRB_INFINITY, 0, GRB_CONTINUOUS, name.str());
                name.str("");
              }
            }
          }
        }

        xt_calblh[t][c][a] = new GRBVar***[data.total_locals];
        for (int l1 = 0; l1 < data.total_locals; ++l1) {
          xt_calblh[t][c][a][l1] = new GRBVar**[data.num_bases];
          for (int b = 0; b < data.num_bases; ++b) {
            xt_calblh[t][c][a][l1][b] = new GRBVar*[data.num_locals];
            for (int l = 0; l < data.num_locals; ++l) {
              xt_calblh[t][c][a][l1][b][l] = new GRBVar[data.num_hosps];
            }
          }
        }
      }
    }
    yt_ahb[t] = new GRBVar**[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      yt_ahb[t][a] = new GRBVar*[data.num_hosps];
      for (int h = 0; h < data.num_hosps; ++h) {
        yt_ahb[t][a][h] = new GRBVar[data.num_bases];
        for (int b = 0; b < data.num_bases; ++b) {
          name << "yt" << t + 1 << "_a" << a << "_h" << h << "_b" << b;
          yt_ahb[t][a][h][b] =
              model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name.str());
          name.str("");
        }
      }
    }
  }

  for (int t = 0; t < data.num_times; ++t) {
    Ct_cl[t] = new GRBVar*[data.types_call];
    for (int c = 0; c < data.types_call; ++c) {
      Ct_cl[t][c] = new GRBVar[data.num_locals];
      for (int l = 0; l < data.num_locals; ++l) {
        name << "Ct" << t + 1 << "_c" << c << "_l" << l;
        Ct_cl[t][c][l] =
            model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name.str());
        name.str("");
      }
    }
    At_ab[t] = new GRBVar*[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      At_ab[t][a] = new GRBVar[data.num_bases];
      for (int b = 0; b < data.num_bases; ++b) {
        name << "At" << t + 1 << "_a" << a << "_b" << b;
        At_ab[t][a][b] =
            model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name.str());
        name.str("");
      }
    }

    At_alb[t] = new GRBVar**[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      At_alb[t][a] = new GRBVar*[data.total_locals];
      for (int l = 0; l < data.total_locals; ++l) {
        At_alb[t][a][l] = new GRBVar[data.num_bases];
        for (int b = 0; b < data.num_bases; ++b) {
          // Increase L_tab to 0...T+1?
          if (data.L_tab[t][a][b].find(l) != data.L_tab[t][a][b].end()) {
            name << "At" << t + 1 << "_a" << a << "_l" << l;
            name << "_b" << b;
            At_alb[t][a][l][b] =
                model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name.str());
            name.str("");
          }
        }
      }
    }
  }

  // Build the objective function. GRBLinExpr is a expression
  // built by the GRBVars and constants.
  GRBLinExpr fo;
  double factor = 1;
  for (auto a : data.A[type_call0]) {
    factor = 1;
    for (int b = 0; b < data.num_bases; ++b) {
      for (auto h : data.H[type_call0][local_call0]) {
        fo += f(0, type_call0, a, data.bases[b], local_call0, h, data, factor) *
              x0_abh[a][b][h];
      }
    }

    for (int h1 = 0; h1 < data.num_hosps; ++h1) {
      for (auto h : data.H[type_call0][local_call0]) {
        fo += f(0, type_call0, a, data.hospitals[h1], local_call0, h, data,
                factor) *
              x0_ahh[a][h1][h];
      }
    }

    for (int b = 0; b < data.num_bases; ++b) {
      for (auto l1 : data.L_tab[0][a][b]) {
        for (auto h : data.H[type_call0][local_call0]) {
          fo += f(0, type_call0, a, l1, data.bases[b], local_call0, h, data,
                  factor) *
                x0_albh[a][l1][b][h];
        }
      }
    }
  }

  factor = 1;
  for (int t = 1; t < data.num_times; ++t) {
    for (int c = 0; c < data.types_call; ++c) {
      for (auto a : data.A[c]) {
        for (int l = 0; l < data.num_locals; ++l) {
          for (auto h : data.H[c][l]) {
            for (int b = 0; b < data.num_bases; ++b) {
              fo += f(t, c, a, data.bases[b], l, h, data, factor) *
                    xt_cablh[t - 1][c][a][b][l][h];
            }

            for (int h1 = 0; h1 < data.num_hosps; ++h1) {
              fo += f(t, c, a, data.hospitals[h1], l, h, data, factor) *
                    xt_cahlh[t - 1][c][a][h1][l][h];
            }

            // for(int l1 = 0; l1 < data.total_locals; ++l1){
            // 	for(int b = 0; b < data.num_bases; ++b){
            // 		if(data.L_tab[t][a][b].find(l1) !=
            // 			data.L_tab[t][a][b].end()){
            // 			fo += f(t-1,c,a,l1,data.bases[b],
            // 			l, h, data,
            // 			factor)*
            // 			xt_calblh[t-1][c][a][l1][b][l][h];
            // 		}
            // 	}
            // }
          }
        }
      }
    }
  }

  for (int t = 0; t < data.num_times; ++t) {
    for (int c = 0; c < data.types_call; ++c) {
      factor = data.ins.penalties[c] * 36000;
      for (int l = 0; l < data.num_locals; ++l) {
        fo += g_tcl(t + 1, c, l, factor) * Ct_cl[t][c][l];
      }
    }
    factor = 0;
    for (int a = 0; a < data.types_amb; ++a) {
      for (int b = 0; b < data.num_bases; ++b) {
        fo += g_tab(t + 1, a, b, factor) * At_ab[t][a][b];
        for (auto l1 : data.L_tab[t][a][b]) {
          fo += g_talb(t + 1, a, l1, data.bases[b], data, factor) *
                At_alb[t][a][l1][b];
        }
      }
    }
  }

  // Adds the objective function to the model
  model.setObjective(fo, GRB_MINIMIZE);

  con_base_t0 = new GRBConstr*[data.types_amb];
  for (int a = 0; a < data.types_amb; ++a) {
    con_base_t0[a] = new GRBConstr[data.num_bases];
  }

  con_hospital_t0 = new GRBConstr*[data.types_amb];
  for (int a = 0; a < data.types_amb; ++a) {
    con_hospital_t0[a] = new GRBConstr[data.num_hosps];
  }
  con_location_t0 = new GRBConstr**[data.types_amb];
  for (int a = 0; a < data.types_amb; ++a) {
    con_location_t0[a] = new GRBConstr*[data.total_locals];
    for (int l1 = 0; l1 < data.total_locals; ++l1) {
      con_location_t0[a][l1] = new GRBConstr[data.num_bases];
    }
  }
  con_queue_t0 = new GRBConstr*[data.types_call];
  for (int c = 0; c < data.types_call; ++c) {
    con_queue_t0[c] = new GRBConstr[data.num_locals];
  }

  con_base = new GRBConstr**[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    con_base[t] = new GRBConstr*[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      con_base[t][a] = new GRBConstr[data.num_bases];
    }
  }

  con_hospital = new GRBConstr**[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    con_hospital[t] = new GRBConstr*[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      con_hospital[t][a] = new GRBConstr[data.num_hosps];
    }
  }

  con_location = new GRBConstr***[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    con_location[t] = new GRBConstr**[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      con_location[t][a] = new GRBConstr*[data.total_locals];
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        con_location[t][a][l1] = new GRBConstr[data.num_bases];
      }
    }
  }

  con_queue = new GRBConstr**[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    con_queue[t] = new GRBConstr*[data.types_call];
    for (int c = 0; c < data.types_call; ++c) {
      con_queue[t][c] = new GRBConstr[data.num_locals];
    }
  }

  con_amb_location = new GRBConstr***[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    con_amb_location[t] = new GRBConstr**[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      con_amb_location[t][a] = new GRBConstr*[data.total_locals];
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        con_amb_location[t][a][l1] = new GRBConstr[data.num_bases];
      }
    }
  }

  con_amb_base = new GRBConstr**[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    con_amb_base[t] = new GRBConstr*[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      con_amb_base[t][a] = new GRBConstr[data.num_bases];
    }
  }

  con_cap_base = new GRBConstr*[data.num_times];
  for (int t = 0; t < data.num_times; ++t) {
    con_cap_base[t] = new GRBConstr[data.num_bases];
  }

  beta0_ab = new double*[data.types_amb];
  alpha0_ah = new double*[data.types_amb];
  psi0_alb = new double**[data.types_amb];
  for (int a = 0; a < data.types_amb; ++a) {
    beta0_ab[a] = new double[data.num_bases];
    for (int b = 0; b < data.num_bases; ++b) {
      beta0_ab[a][b] = 0;
    }

    alpha0_ah[a] = new double[data.num_hosps];
    for (int h = 0; h < data.num_hosps; ++h) {
      alpha0_ah[a][h] = 0;
    }

    psi0_alb[a] = new double*[data.total_locals];
    for (int l1 = 0; l1 < data.total_locals; ++l1) {
      psi0_alb[a][l1] = new double[data.num_bases];
      for (int b = 0; b < data.num_bases; ++b) {
        psi0_alb[a][l1][b] = 0;
      }
    }
  }

  phi0_cl = new double*[data.types_call];
  for (int c = 0; c < data.types_call; ++c) {
    phi0_cl[c] = new double[data.num_locals];
    for (int l = 0; l < data.num_locals; ++l) {
      phi0_cl[c][l] = 0;
    }
  }

  beta_tab = new double**[data.num_times - 1];
  alpha_tah = new double**[data.num_times - 1];
  psi_talb = new double***[data.num_times - 1];
  phi_tcl = new double**[data.num_times - 1];
  ups_tb = new double*[data.num_times - 1];
  theta_talb = new double***[data.num_times - 1];
  gamma_tab = new double**[data.num_times];

  for (int t = 0; t < data.num_times - 1; ++t) {
    beta_tab[t] = new double*[data.types_amb];
    alpha_tah[t] = new double*[data.types_amb];
    psi_talb[t] = new double**[data.types_amb];
    theta_talb[t] = new double**[data.types_amb];
    gamma_tab[t] = new double*[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      alpha_tah[t][a] = new double[data.num_hosps];
      psi_talb[t][a] = new double*[data.total_locals];
      theta_talb[t][a] = new double*[data.total_locals];
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        psi_talb[t][a][l1] = new double[data.num_bases];
        theta_talb[t][a][l1] = new double[data.num_bases];
        for (int b = 0; b < data.num_bases; ++b) {
          psi_talb[t][a][l1][b] = 0;
          theta_talb[t][a][l1][b] = 0;
        }
      }
      beta_tab[t][a] = new double[data.num_bases];
      gamma_tab[t][a] = new double[data.num_bases];

      for (int b = 0; b < data.num_bases; ++b) {
        gamma_tab[t][a][b] = 0;
        beta_tab[t][a][b] = 0;
      }
    }
    phi_tcl[t] = new double*[data.types_call];
    for (int c = 0; c < data.types_call; ++c) {
      phi_tcl[t][c] = new double[data.num_locals];
      for (int l = 0; l < data.num_locals; ++l) {
        phi_tcl[t][c][l] = 0;
      }
    }
    ups_tb[t] = new double[data.num_bases];
    for (int b = 0; b < data.num_bases; ++b) {
      ups_tb[t][b] = 0;
    }
  }
  gamma_tab[data.num_times - 1] = new double*[data.types_amb];
  for (int a = 0; a < data.types_amb; ++a) {
    gamma_tab[data.num_times - 1][a] = new double[data.num_bases];
    for (int b = 0; b < data.num_bases; ++b) {
      gamma_tab[data.num_times - 1][a][b] = 0;
    }
  }

  add_bases_constraints_t0();
  add_hospitals_constraints_t0();
  add_locations_constraints_t0();
  add_queues_constraints_t0();

  add_bases_constraints();
  add_hospitals_constraints();
  add_locations_constraints();
  add_queues_constraints();
  add_ambs_location_constraints();
  add_base_cap_constraints();

  model.update();
}

void CGCall::solve() {
  iter = 0;
  bool column_found = false;
  model.set(GRB_IntParam_OutputFlag, 0);
  model.set(GRB_IntParam_Presolve, 0);
  model.set(GRB_IntParam_Threads, 1);
  // std::cout << "Number of Variables: " << model.get(GRB_IntAttr_NumVars) <<
  // "\n";
  double sum_master_time = 0;
  double sum_pricing_time = 0;
  double sum_dual_time = 0;
  do {
    model.update();
    model.optimize();
    sum_master_time += model.get(GRB_DoubleAttr_Runtime);
    // double z =  model.get(GRB_DoubleAttr_ObjVal);
    // std::cout << "Iter " << iter << ", Z = " << z << "\n";
    std::chrono::high_resolution_clock::time_point t_dual =
        std::chrono::high_resolution_clock::now();
    set_dual();
    std::chrono::high_resolution_clock::time_point dt_dual =
        std::chrono::high_resolution_clock::now();
    sum_dual_time +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(dt_dual - t_dual)
            .count() /
        pow(10, 9);
    std::chrono::high_resolution_clock::time_point t0 =
        std::chrono::high_resolution_clock::now();
    column_found = pricing();
    std::chrono::high_resolution_clock::time_point dt =
        std::chrono::high_resolution_clock::now();
    sum_pricing_time +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(dt - t0).count() /
        pow(10, 9);
    // fmt::print("iter {} {}\n", iter,model.get(GRB_DoubleAttr_ObjVal));
    ++iter;
    if (iter > 10) {
      break;
    }
  } while (column_found);

  // std::cout << "Average Pricing Time: " << sum_pricing_time/iter << "\n";
  // std::cout << "Average Master Time: " << sum_master_time/iter << "\n";
  // std::cout << "Average Dual Time: " << sum_dual_time/iter << "\n";

  // model.optimize();
  // std::cout << "Columns Added: " << columns_added << "\n";

  src_type = -1;
  src_location = -1;
  amb_type = -1;
  loc_base = -1;

  // if(sim.events.top()->id == 14){
  // 	model.computeIIS();
  // 	model.write("id_14.ilp");
  // 	for(auto& amb: sim.ambulances){
  // 		std::cout << amb << "\n";
  // 	}
  // }

  if (model.get(GRB_IntAttr_Status) != 2) {
    fmt::print("Status {}\n", model.get(GRB_IntAttr_Status));
    model.computeIIS();
    model.write("cg.ilp");
    std::cin.get();
  }

  // if(model.get(GRB_DoubleAttr_ObjVal) > pow(10,6)){
  // 	for(auto& call: data.calls){
  // 		cout << call << "\n";
  // 	}
  // 	print_solution();
  // 	cin.get();
  // }

  // bool found = false;

  // for(int a = 0; a < data.types_amb; ++a){
  // 	for(int b = 0; b < data.num_bases; ++b){
  // 		for(int h = 0; h < data.num_hosps; ++h){
  // 			if(data.A[type_call0].find(a) !=
  // data.A[type_call0].end() &&
  // data.H[type_call0][local_call0].find(h) !=
  // 				data.H[type_call0][local_call0].end() &&
  // 				x0_abh[a][b][h].get(GRB_DoubleAttr_X) > 0.5){
  // 				src_type = BASE;
  // 				src_location = b;
  // 				dst_hosp = h;
  // 				amb_type = a;
  // 				found = true;
  // 				// std::cout <<
  // x0_abh[a][b][h].get(GRB_StringAttr_VarName) << " ";
  // 				// std::cout <<
  // x0_abh[a][b][h].get(GRB_DoubleAttr_X) << "\n";
  // break;
  // 			}
  // 		}
  // 		if(found) break;
  // 	}
  // 	if(found) break;
  // }

  // for(int a = 0; a < data.types_amb; ++a){
  // 	for(int h1 = 0; h1 < data.num_hosps; ++h1){
  // 		for(int h = 0; h < data.num_hosps; ++h){
  // 			if(data.A[type_call0].find(a) !=
  // data.A[type_call0].end() &&
  // data.H[type_call0][local_call0].find(h) !=
  // 				data.H[type_call0][local_call0].end() &&
  // 				x0_ahh[a][h1][h].get(GRB_DoubleAttr_X) > 0.5){
  // 				src_type = HOSP;
  // 				src_location = h1;
  // 				dst_hosp = h;
  // 				amb_type = a;
  // 				found = true;
  // 				// std::cout <<
  // x0_ahh[a][h1][h].get(GRB_StringAttr_VarName) << " ";
  // 				// std::cout <<
  // x0_ahh[a][h1][h].get(GRB_DoubleAttr_X) << "\n";
  // break;
  // 			}
  // 		}
  // 		if(found) break;
  // 	}
  // 	if(found) break;
  // }

  // for(int a = 0; a < data.types_amb; ++a){
  // 	for(int l1 = 0; l1 < data.total_locals; ++l1){
  // 		for(int b = 0; b < data.num_bases; ++b){
  // 			for(int h = 0; h < data.num_hosps; ++h){
  // 				if(data.A[type_call0].find(a) !=
  // data.A[type_call0].end() &&
  // data.H[type_call0][local_call0].find(h) !=
  // 					data.H[type_call0][local_call0].end() &&
  // 					data.L_tab[0][a][b].find(l1) !=
  // 					data.L_tab[0][a][b].end() &&
  // 					x0_albh[a][l1][b][h].get(GRB_DoubleAttr_X)
  // > 0.5){ 					src_type = LOC;
  // src_location = l1; 					dst_hosp = h;
  // amb_type = a; 					loc_base = b;
  // 					// std::cout <<
  // x0_albh[a][l1][b][h].get(GRB_StringAttr_VarName) << " ";
  // 					// std::cout <<
  // x0_albh[a][l1][b][h].get(GRB_DoubleAttr_X) << "\n";
  // found = true;
  // 				}
  // 			}
  // 			if(found) break;
  // 		}
  // 		if(found) break;
  // 	}
  // 	if(found) break;
  // }
}

void CGCall::set_dual() {
  // for(auto a: data.A[type_call0]){
  // 	for(int b = 0; b < data.num_bases; ++b){
  // 		beta0_ab[a][b] = con_base_t0[a][b].get(GRB_DoubleAttr_Pi);
  // 		// if(abs(beta0_ab[a][b]) > g_params.g_params.EPS){
  // 		// 	std::cout << "beta0_" << a << "_" << b << ": ";
  // 		// 	std::cout << beta0_ab[a][b] << "\n";
  // 		// }
  // 	}
  // }

  // for(int a = 0; a < data.types_amb; ++a){
  // 	for(int h = 0; h < data.num_hosps; ++h){
  // 		alpha0_ah[a][h] = con_hospital_t0[a][h].get(GRB_DoubleAttr_Pi);
  // 		// if(abs(alpha0_ah[a][h] )> g_params.g_params.EPS){
  // 		// 	std::cout << "alpha0_" << a << "_" << h << ": ";
  // 		// 	std::cout << alpha0_ah[a][h] << "\n";
  // 		// }
  // 	}
  // }

  // for(int a = 0; a < data.types_amb; ++a){
  // 	for(int b = 0; b < data.num_bases; ++b){
  // 		for(int l1 = 0; l1 < data.total_locals;++l1){
  // 			if(data.L_tab[0][a][b].find(l1) !=
  // data.L_tab[0][a][b].end()){ psi0_alb[a][l1][b] =
  // con_location_t0[a][l1][b].get(GRB_DoubleAttr_Pi);
  // 				// if(abs(psi0_alb[a][l1][b]) >
  // g_params.g_params.EPS){
  // 				// 	std::cout << "psi0_" << a << "_" << l1
  // << "_" << b << ": ";
  // 				// 	std::cout << psi0_alb[a][l1][b] << "\n";
  // 				// }
  // 			}
  // 		}
  // 	}
  // }

  // for(int c = 0; c < data.types_call; ++c){
  // 	for(int l = 0; l < data.num_locals; ++l){
  // 		phi0_cl[c][l] = con_queue_t0[c][l].get(GRB_DoubleAttr_Pi);
  // 		// if(abs(phi0_cl[c][l]) > g_params.g_params.EPS){
  // 		// 	std::cout << "phi0_" << c << "_" << l << ": ";
  // 		// 	std::cout << phi0_cl[c][l] << "\n";
  // 		// }
  // 	}
  // }
  // std::ofstream arq("dual_sols.txt");
  // arq << "Iter " << iter  << ":\n";
  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int b = 0; b < data.num_bases; ++b) {
        beta_tab[t][a][b] = con_base[t][a][b].get(GRB_DoubleAttr_Pi);
        gamma_tab[t][a][b] = con_amb_base[t][a][b].get(GRB_DoubleAttr_Pi);
        // if(gamma_tab[t][a][b] > g_params.EPS){
        // 	std::cout << "gamma_" << t << "_" << a << "_" << b << ": ";
        // 	std::cout << gamma_tab[t][a][b] << "\n";
        // }
      }

      for (int h = 0; h < data.num_hosps; ++h) {
        alpha_tah[t][a][h] = con_hospital[t][a][h].get(GRB_DoubleAttr_Pi);
        // if(abs(alpha_tah[t][a][h]) > g_params.EPS){
        // 	std::cout << "alpha_" << t+1 << "_" << a << "_" << h << ",";
        // 	std::cout << alpha_tah[t][a][h] << "\n";
        // }
      }

      for (int l1 = 0; l1 < data.num_locals; ++l1) {
        for (int b = 0; b < data.num_bases; ++b) {
          if (data.L_tab[1][a][b].find(l1) != data.L_tab[1][a][b].end()) {
            psi_talb[t][a][l1][b] =
                con_location[t][a][l1][b].get(GRB_DoubleAttr_Pi);

            // if(abs(psi_talb[t][a][l1][b]) > g_params.EPS){
            // 	arq << "psi_" << t+1 << "_" << a << "_" << l1 << "_";
            // 	arq << b << "," << psi_talb[t][a][l1][b] << "\n";
            // }
            if (l1 != data.bases[b]) {
              // if(abs(theta_talb[t][a][l1][b]-con_amb_location[t][a][l1][b].get(
              // 	GRB_DoubleAttr_Pi)) > g_params.EPS){
              // 		arq << "theta_" << t+1 << "_" << a << "_" << l1
              // << "_"; 		arq << b << "," <<
              // theta_talb[t][a][l1][b] << ","; 		arq
              // <<con_amb_location[t][a][l1][b].get(GRB_DoubleAttr_Pi);
              // arq <<
              // "\n";
              // }
              theta_talb[t][a][l1][b] =
                  con_amb_location[t][a][l1][b].get(GRB_DoubleAttr_Pi);
              // if(iter >= 4 && abs(theta_talb[t][a][l1][b]) > g_params.EPS){
              // 	arq << "theta_" << t+1 << "_" << a << "_" << l1 << "_";
              // 	arq << b << "," << theta_talb[t][a][l1][b] << "\n";
              // }
            }
          }
        }
      }
    }

    for (int b = 0; b < data.num_bases; ++b) {
      ups_tb[t][b] = con_cap_base[t][b].get(GRB_DoubleAttr_Pi);
      // if(abs(ups_tb[t][b]) > g_params.EPS){
      // 	std::cout << "ups_" << t+1 << "_" << b << ": " << ups_tb[t][b]
      // << "\n";
      // }
    }

    for (int c = 0; c < data.types_call; ++c) {
      for (int l = 0; l < data.num_locals; ++l) {
        // if(abs(phi_tcl[t][c][l] - con_queue[t][c][l].get(GRB_DoubleAttr_Pi))
        // > g_params.EPS){ 	arq << "phi_" << t+1 << "_" << c << "_" << l;
        // arq <<
        // "," <<  phi_tcl[t][c][l] << ","; 	arq <<
        // con_queue[t][c][l].get(GRB_DoubleAttr_Pi) << "\n";
        // }
        phi_tcl[t][c][l] = con_queue[t][c][l].get(GRB_DoubleAttr_Pi);
        // if(iter >= 4 && abs(phi_tcl[t][c][l]) > g_params.EPS){
        // 	std::cout << "phi_" << t+1 << "_" << c << "_" << l;
        // 	std::cout << "," << phi_tcl[t][c][l] << "\n";
        // }
      }
    }
  }

  // for(int a = 0; a < data.types_amb; ++a){
  // 	for(int b = 0; b < data.num_bases; ++b){
  // 		int t = data.num_times-1;
  // 		gamma_tab[t][a][b] =
  // con_amb_base[t][a][b].get(GRB_DoubleAttr_Pi);
  // if(abs(gamma_tab[t][a][b] > g_params.EPS){
  // std::cout << "gamma_" << t << "_" << a << "_" << b << ": ";
  // std::cout << gamma_tab[t][a][b] << "\n";
  // 		}
  // 	}
  // }
}

double CGCall::sub_problem(int t, int a, int l, int l1, int c, int h, int b) {
  int next_l = data.L(t + 1, a, l1, data.bases[b]);
  int arr_time =
      get_time(data, (t + 1) * data.quantum + data.tao[t][c][a][l1][l][h]);
  const int factor = 1;
  double reduced_cost = f(t + 1, c, a, l1, l, h, data, factor);
  if (next_l == data.bases[b]) {
    reduced_cost += beta_tab[t][a][b];
  }

  if (arr_time - 1 < data.num_times - 1) {
    reduced_cost -= alpha_tah[arr_time - 1][a][h];
  }

  if (data.L_tab[t + 1][a][b].find(next_l) != data.L_tab[t + 1][a][b].end()) {
    // std::cout << "psi!!!!\n";
    reduced_cost += psi_talb[t][a][next_l][b];
  }

  reduced_cost += phi_tcl[t][c][l];
  reduced_cost -= theta_talb[t][a][l1][b];

  // if(reduced_cost < -g_params.EPS){
  // 	std::cout << t+1 << " " << c << " " << a << " " << l1 << " ";
  // 	std::cout << b << " " << l << " " << h << ": ";
  // 	std::cout << reduced_cost << "\n";
  // }
  // if(reduced_cost < -g_params.EPS){
  // std::cout << t+1 << " " << c << " " << a << " " << l1 << " ";
  // std::cout << b << " " << l << " " << h << ": ";
  // std::cout << "next_l | base = " << next_l << " " << data.bases[b] << "\n";
  // std::cout << "arrival_time = " << arr_time << "\n";
  // std::cout << "tao = "<< data.tao[t+1][c][a][l1][l][h] << "\n";
  // std::cout << "f = " << f(t+1,c,a,l1,l,h,data,factor) << "\n";
  // std::cout << "beta = " << beta_tab[t][a][b] << "\n";
  // if(arr_time < data.num_times)
  // std::cout << "alpha = " << alpha_tah[arr_time-1][a][h] << "\n";
  // std::cout << "psi = " << psi_talb[t][a][next_l][b] << "\n";
  // std::cout << "phi = " << phi_tcl[t][c][l] << "\n";
  // std::cout << "theta = " << theta_talb[t][a][l1][b] << "\n";
  // std::cout << "reduced_cost = " << reduced_cost << "\n";
  // std::cin.get();
  // }

  return reduced_cost;
}

bool CGCall::pricing() {
  bool optimality_verified = true;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int l = 0; l < data.num_locals; ++l) {
        for (auto l1 : closest_l1_to_l[l]) {
          vh_tall[t][a][l1][l] = h_tall(t, a, l1, l);
        }
        vb_tal[t][a][l] = b_tal(t, a, l);
      }
    }

    for (int l = 0; l < data.num_locals; ++l) {
      vc_tl[t][l] = c_tl(t, l);
    }
  }

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int l = 0; l < data.num_locals; ++l) {
        bool negative_found = false;
        for (auto l1 : closest_l1_to_l[l]) {
          // int h = h_tall(t,a,l1,l);
          // int c = c_tl(t,l);
          // int b = b_tal(t,a,l1);
          int h = vh_tall[t][a][l1][l];
          int c = vc_tl[t][l];
          int b = vb_tal[t][a][l1];
          if (b != -1 && data.L_tab[t + 1][a][b].find(l1) !=
                             data.L_tab[t + 1][a][b].end()) {
            double reduced_cost = sub_problem(t, a, l, l1, c, h, b);
            if (reduced_cost < -g_params.EPS) {
              // std::cout << t+1 << " " << c << " " << a << " " << l1 << " " <<
              // b; std::cout << " " << l << " " << h << ": " << reduced_cost <<
              // "\n"; std::cin.get();
              negative_found = true;
              optimality_verified = false;
              std::vector<int> column{t, c, a, l1, b, l, h};
              add_column(column);
              ++columns_added;
            }
          }
          if (negative_found) break;
        }
      }
    }
  }
  // std::cout << "END PRICING\n";
  // std::cin.get();
  return !optimality_verified;
}

void CGCall::add_column(std::vector<int>& col) {
  GRBColumn column = GRBColumn();
  std::stringstream name;

  int t = col[0], c = col[1], a = col[2], l1 = col[3], b = col[4];
  int l = col[5], h = col[6];
  int t2 =
      get_time(data, (t + 1) * data.quantum + data.tao[t + 1][c][a][l1][l][h]);

  // name << "xt" << t+1 << "_c" << c << "_a" << a << "_l'";
  // name << l1 << "_b" << b << "_l" << l << "_h" << h;

  if (data.L_tab[t + 1][a][b].find(l1) != data.L_tab[t + 1][a][b].end() &&
      data.L(t, a, l1, data.bases[b]) == data.bases[b] &&
      data.H[c][l].find(h) != data.H[c][l].end() &&
      data.A[c].find(a) != data.A[c].end()) {
    // con_base[t][a][b]
    // std::cout << "-1 " << con_base[t][a][b].get(GRB_StringAttr_ConstrName) <<
    // "\n";
    column.addTerm(-1, con_base[t][a][b]);
  }

  if (t2 < data.num_times && data.H[c][l].find(h) != data.H[c][l].end() &&
      data.A[c].find(a) != data.A[c].end() &&
      data.L_tab[t][a][b].find(l1) != data.L_tab[t][a][b].end()) {
    // con_hospital[t2][a][h]
    // std::cout << "1 " <<
    // con_hospital[t2][a][h].get(GRB_StringAttr_ConstrName) << "\n";
    column.addTerm(1, con_hospital[t2 - 1][a][h]);
  }

  for (int l2 = 0; l2 < data.total_locals; ++l2) {
    if (data.L_tab[t + 1][a][b].find(l2) != data.L_tab[t + 1][a][b].end() &&
        data.L(t, a, l1, data.bases[b]) == l2 &&
        data.A[c].find(a) != data.A[c].end() &&
        data.H[c][l].find(h) != data.H[c][l].end()) {
      // con_location[t][a][l2][b]
      // std::cout << "-1 " <<
      // con_location[t][a][l2][b].get(GRB_StringAttr_ConstrName) << "\n";
      column.addTerm(-1, con_location[t][a][l2][b]);
    }
  }

  if (data.H[c][l].find(h) != data.H[c][l].end() &&
      data.A[c].find(a) != data.A[c].end() &&
      data.L_tab[t + 1][a][b].find(l1) != data.L_tab[t + 1][a][b].end()) {
    // con_queue[t][c][l]
    // std::cout << "-1 " << con_queue[t][c][l].get(GRB_StringAttr_ConstrName)
    // << "\n";
    column.addTerm(-1, con_queue[t][c][l]);
  }

  if (data.L_tab[t + 1][a][b].find(l1) != data.L_tab[t + 1][a][b].end() &&
      data.H[c][l].find(h) != data.H[c][l].end() &&
      data.A[c].find(a) != data.A[c].end()) {
    // std::cout << "1 " <<
    // con_amb_location[t][a][l1][b].get(GRB_StringAttr_ConstrName) << "\n";
    column.addTerm(1, con_amb_location[t][a][l1][b]);
  }
  // std::cout << "COLUMN: " << t << " " << c << " " << a << " " << l1 << " " <<
  // b << " "; std::cout << l << " " << h << "\n";
  try {
    name << "xt" << t + 1 << "_c" << c << "_a" << a << "_l'" << l1 << "_b" << b
         << "_l";
    name << l << "_h" << h;
    double factor = 1;
    double f_val = f(t + 1, c, a, l1, data.bases[b], l, h, data, factor);
    xt_calblh[t][c][a][l1][b][l][h] = model.addVar(
        0, GRB_INFINITY, f_val, GRB_CONTINUOUS, column, name.str());
  } catch (GRBException& ex) {
    std::cout << ex.getErrorCode() << ": " << ex.getMessage() << "\n";
    std::cin.get();
  }
}

void CGCall::add_bases_constraints_t0() {
  GRBLinExpr con;
  std::stringstream name;

  for (int a = 0; a < data.types_amb; ++a) {
    for (int b = 0; b < data.num_bases; ++b) {
      con += A0_ab[a][b];
      for (int h = 0; h < data.num_hosps; ++h) {
        if (data.L(0, a, data.hospitals[h], data.bases[b]) == data.bases[b]) {
          con += y0_ahb[a][h][b];
        }
      }
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        if (data.L_tab[0][a][b].find(l1) != data.L_tab[0][a][b].end() &&
            data.L(0, a, l1, data.bases[b]) == data.bases[b]) {
          con += A0_alb[a][l1][b];
          for (int h = 0; h < data.num_hosps; ++h) {
            if (data.H[type_call0][local_call0].find(h) !=
                    data.H[type_call0][local_call0].end() &&
                data.A[type_call0].find(a) != data.A[type_call0].end()) {
              con -= x0_albh[a][l1][b][h];
            }
          }
        }
      }

      for (int h = 0; h < data.num_hosps; ++h) {
        if (data.H[type_call0][local_call0].find(h) !=
                data.H[type_call0][local_call0].end() &&
            data.A[type_call0].find(a) != data.A[type_call0].end()) {
          con -= x0_abh[a][b][h];
        }
      }

      name << "flow_base_a" << a << "_b" << b;
      con_base_t0[a][b] =
          model.addConstr(con, GRB_EQUAL, At_ab[0][a][b], name.str());
      name.str("");
      con = 0;
    }
  }
}

void CGCall::add_hospitals_constraints_t0() {
  GRBLinExpr con;
  std::stringstream name;

  for (int a = 0; a < data.types_amb; ++a) {
    for (int h = 0; h < data.num_hosps; ++h) {
      for (int h1 = 0; h1 < data.num_hosps; ++h1) {
        if (data.H[type_call0][local_call0].find(h1) !=
                data.H[type_call0][local_call0].end() &&
            data.A[type_call0].find(a) != data.A[type_call0].end()) {
          con += x0_ahh[a][h][h1];
        }
      }
      for (int b = 0; b < data.num_bases; ++b) {
        con += y0_ahb[a][h][b];
      }

      name << "flow_hospital_a" << a << "_h" << h;
      con_hospital_t0[a][h] =
          model.addConstr(con, GRB_EQUAL, A0_ah[a][h], name.str());
      name.str("");
      con = 0;
    }
  }
}

void CGCall::add_locations_constraints_t0() {
  GRBLinExpr con;
  std::stringstream name;

  for (int a = 0; a < data.types_amb; ++a) {
    for (int b = 0; b < data.num_bases; ++b) {
      for (int l = 0; l < data.total_locals; ++l) {
        if (data.L_tab[0][a][b].find(l) != data.L_tab[0][a][b].end()) {
          for (int h = 0; h < data.num_hosps; ++h) {
            if (data.L(0, a, data.hospitals[h], data.bases[b]) == l) {
              con += y0_ahb[a][h][b];
            }
          }
          for (int l1 = 0; l1 < data.total_locals; ++l1) {
            if (data.L_tab[0][a][b].find(l1) != data.L_tab[0][a][b].end() &&
                data.L(0, a, l1, data.bases[b]) == l) {
              con += A0_alb[a][l1][b];
              for (int h = 0; h < data.num_hosps; ++h) {
                if (data.H[type_call0][local_call0].find(h) !=
                        data.H[type_call0][local_call0].end() &&
                    data.A[type_call0].find(a) != data.A[type_call0].end()) {
                  con -= x0_albh[a][l1][b][h];
                }
              }
            }
          }
          name << "flow_locals_a" << a << "_b" << b << "_l" << l;
          con_location_t0[a][l][b] =
              model.addConstr(con, GRB_EQUAL, At_alb[0][a][l][b], name.str());

          name.str("");
          con = 0;
        }
      }
    }
  }
}

void CGCall::add_queues_constraints_t0() {
  GRBLinExpr con;
  std::stringstream name;

  for (int c = 0; c < data.types_call; ++c) {
    for (int l = 0; l < data.num_locals; ++l) {
      con += C[c][l];
      if (c == type_call0 && l == local_call0) {
        con += 1;
        for (int a = 0; a < data.types_amb; ++a) {
          if (data.A[c].find(a) != data.A[c].end()) {
            for (int b = 0; b < data.num_bases; ++b) {
              for (int h = 0; h < data.num_hosps; ++h) {
                if (data.H[c][l].find(h) != data.H[c][l].end() &&
                    data.A[type_call0].find(a) != data.A[type_call0].end()) {
                  con -= x0_abh[a][b][h];
                }
              }
            }
          }
        }
        for (int a = 0; a < data.types_amb; ++a) {
          if (data.A[c].find(a) != data.A[c].end()) {
            for (int h1 = 0; h1 < data.num_hosps; ++h1) {
              for (int h = 0; h < data.num_hosps; ++h) {
                if (data.H[c][l].find(h) != data.H[c][l].end() &&
                    data.A[type_call0].find(a) != data.A[type_call0].end()) {
                  con -= x0_ahh[a][h1][h];
                }
              }
            }
          }
        }

        for (int a = 0; a < data.types_amb; ++a) {
          if (data.A[c].find(a) != data.A[c].end()) {
            for (int b = 0; b < data.num_bases; ++b) {
              for (int l1 = 0; l1 < data.total_locals; ++l1) {
                for (int h = 0; h < data.num_hosps; ++h) {
                  if (data.H[c][l].find(h) != data.H[c][l].end() &&
                      data.A[type_call0].find(a) != data.A[type_call0].end() &&
                      data.L_tab[0][a][b].find(l1) !=
                          data.L_tab[0][a][b].end()) {
                    con -= x0_albh[a][l1][b][h];
                  }
                }
              }
            }
          }
        }
      }

      name << "flow_queue_c" << c << "_l" << l;

      con_queue_t0[c][l] =
          model.addConstr(con, GRB_EQUAL, Ct_cl[0][c][l], name.str());
      name.str("");
      con = 0;
    }
  }
}

void CGCall::add_bases_constraints() {
  GRBLinExpr con;
  std::stringstream name;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int b = 0; b < data.num_bases; ++b) {
        con += At_ab[t][a][b];
        for (int h = 0; h < data.num_hosps; ++h) {
          if (data.L(t, a, data.hospitals[h], data.bases[b]) == data.bases[b]) {
            con += yt_ahb[t][a][h][b];
          }
        }
        for (int l1 = 0; l1 < data.total_locals; ++l1) {
          if (data.L_tab[t + 1][a][b].find(l1) !=
                  data.L_tab[t + 1][a][b].end() &&
              data.L(t, a, l1, data.bases[b]) == data.bases[b]) {
            con += At_alb[t][a][l1][b];
            // for(int c = 0; c < data.types_call; ++c){
            // 	for(int l = 0; l < data.num_locals; ++l){
            // 		for(int h = 0; h < data.num_hosps; ++h){
            // 			if(data.H[c][l].find(h) !=
            // 				data.H[c][l].end() &&
            // 				data.A[c].find(a) != data.A[c].end()){
            // 				//If using xt0, put t+1 on 1st index
            // 				con -= xt_calblh[t][c][a][l1][b][l][h];
            // 			}
            // 		}
            // 	}
            // }
          }
        }

        for (int c = 0; c < data.types_call; ++c) {
          for (int l = 0; l < data.num_locals; ++l) {
            for (int h = 0; h < data.num_hosps; ++h) {
              if (data.H[c][l].find(h) != data.H[c][l].end() &&
                  data.A[c].find(a) != data.A[c].end()) {
                // If using xt0, put t+1 on 1st index
                con -= xt_cablh[t][c][a][b][l][h];
              }
            }
          }
        }

        name << "flow_bases_t" << t + 1 << "_a" << a << "_b" << b;
        try {
          con_base[t][a][b] =
              model.addConstr(con, GRB_EQUAL, At_ab[t + 1][a][b], name.str());
        } catch (GRBException ex) {
          std::cout << ex.getErrorCode() << ": ";
          std::cout << ex.getMessage() << "\n";
          exit(1);
        }
        name.str("");
        con = 0;
      }
    }
  }
}

void CGCall::add_hospitals_constraints() {
  GRBLinExpr con;
  std::stringstream name;
  int quantum = data.quantum;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int h = 0; h < data.num_hosps; ++h) {
        for (int c = 0; c < data.types_call; ++c) {
          for (int l = 0; l < data.num_locals; ++l) {
            for (int h1 = 0; h1 < data.num_hosps; ++h1) {
              if (data.H[c][l].find(h1) != data.H[c][l].end() &&
                  data.A[c].find(a) != data.A[c].end()) {
                // If using xt0, put t+1 on 1st index
                con -= xt_cahlh[t][c][a][h][l][h1];
              }
            }
          }
        }

        for (int b = 0; b < data.num_bases; ++b) {
          con -= yt_ahb[t][a][h][b];
        }

        for (int c = 0; c < data.types_call; ++c) {
          for (int l1 = 0; l1 < data.total_locals; ++l1) {
            for (int l = 0; l < data.num_locals; ++l) {
              if ((t)*quantum <= data.tao[0][c][a][l1][l][h] &&
                  data.tao[0][c][a][l1][l][h] <= (t + 1) * quantum) {
                con += A0_callh[c][a][l1][l][h];
              }
            }
          }
        }

        for (int c = 0; c < data.types_call; ++c) {
          for (int l1 = 0; l1 < data.total_locals; ++l1) {
            if ((t)*quantum <= data.tao_0[c][a][l1][h] &&
                data.tao_0[c][a][l1][h] <= (t + 1) * quantum) {
              con += A0_calh[c][a][l1][h];
            }
          }
        }

        for (int b = 0; b < data.num_bases; ++b) {
          if ((t)*quantum <=
                  data.tao[0][type_call0][a][data.bases[b]][local_call0][h] &&
              data.tao[0][type_call0][a][data.bases[b]][local_call0][h] <=
                  (t + 1) * quantum &&
              data.H[type_call0][local_call0].find(h) !=
                  data.H[type_call0][local_call0].end() &&
              data.A[type_call0].find(a) != data.A[type_call0].end()) {
            con += x0_abh[a][b][h];
          }
        }

        for (int c = 0; c < data.types_call; ++c) {
          for (int l = 0; l < data.num_locals; ++l) {
            if (data.H[c][l].find(h) != data.H[c][l].end()) {
              for (int b = 0; b < data.num_bases; ++b) {
                for (int t1 = 1; t1 < t + 1; ++t1) {
                  if ((t)*quantum <=
                          t1 * quantum +
                              data.tao[t1][c][a][data.bases[b]][l][h] &&
                      t1 * quantum + data.tao[t1][c][a][data.bases[b]][l][h] <=
                          (t + 1) * quantum &&
                      data.A[c].find(a) != data.A[c].end()) {
                    con += xt_cablh[t1 - 1][c][a][b][l][h];
                  }
                }
              }
            }
          }
        }

        for (int h1 = 0; h1 < data.num_hosps; ++h1) {
          if ((t)*quantum <= data.tao[0][type_call0][a][data.hospitals[h1]]
                                     [local_call0][h] &&
              data.tao[0][type_call0][a][data.hospitals[h1]][local_call0][h] <=
                  (t + 1) * quantum &&
              data.H[type_call0][local_call0].find(h) !=
                  data.H[type_call0][local_call0].end() &&
              data.A[type_call0].find(a) != data.A[type_call0].end()) {
            con += x0_ahh[a][h1][h];
          }
        }

        for (int c = 0; c < data.types_call; ++c) {
          for (int l = 0; l < data.num_locals; ++l) {
            if (data.H[c][l].find(h) != data.H[c][l].end()) {
              for (int h1 = 0; h1 < data.num_hosps; ++h1) {
                for (int t1 = 1; t1 < t + 1; ++t1) {
                  if ((t)*quantum <=
                          t1 * quantum +
                              data.tao[t1][c][a][data.hospitals[h1]][l][h] &&
                      t1 * quantum +
                              data.tao[t1][c][a][data.hospitals[h1]][l][h] <=
                          (t + 1) * quantum &&
                      data.A[c].find(a) != data.A[c].end()) {
                    con += xt_cahlh[t1 - 1][c][a][h1][l][h];
                  }
                }
              }
            }
          }
        }

        // for(int c = 0; c < data.types_call; ++c){
        // 	for(int l = 0; l < data.num_locals; ++l){
        // 		if(data.H[c][l].find(h) != data.H[c][l].end()){
        // 			for(int l1 = 0; l1 < data.total_locals; ++l1){
        // 				for(int b = 0; b < data.num_bases; ++b){
        // 					for(int t1 = 1; t1 < t+1; ++t1){
        // 						if((t)*quantum <=
        // t1*quantum+
        // data.tao[t1][c][a][l1][l][h] &&
        // 							t1*quantum+data.tao[t1][c][a][l1][l][h]
        // <= (t+1)*quantum &&
        // data.A[c].find(a) != data.A[c].end()
        // 							&&
        // data.L_tab[t1-1][a][b].find(l1) !=
        // data.L_tab[t1-1][a][b].end()){
        // con
        // += xt_calblh[t1-1][c][a][l1][b][l][h];
        // 						}
        // 					}
        // 				}
        // 			}
        // 		}
        // 	}
        // }

        for (int b = 0; b < data.num_bases; ++b) {
          for (int l1 = 0; l1 < data.total_locals; ++l1) {
            if ((t)*quantum <= data.tao[0][type_call0][a][l1][local_call0][h] &&
                data.tao[0][type_call0][a][l1][local_call0][h] <=
                    (t + 1) * quantum &&
                data.H[type_call0][local_call0].find(h) !=
                    data.H[type_call0][local_call0].end() &&
                data.A[type_call0].find(a) != data.A[type_call0].end() &&
                data.L_tab[0][a][b].find(l1) != data.L_tab[0][a][b].end()) {
              con += x0_albh[a][l1][b][h];
            }
          }
        }

        name << "flow_hospitals_t" << t + 1 << "_a" << a << "_h" << h;
        con_hospital[t][a][h] = model.addConstr(con, GRB_EQUAL, 0, name.str());
        name.str("");
        con = 0;
      }
    }
  }
}

void CGCall::add_locations_constraints() {
  GRBLinExpr con;
  std::stringstream name;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int b = 0; b < data.num_bases; ++b) {
        for (int l1 = 0; l1 < data.total_locals; ++l1) {
          if (data.L_tab[t + 1][a][b].find(l1) !=
              data.L_tab[t + 1][a][b].end()) {
            for (int h1 = 0; h1 < data.num_hosps; ++h1) {
              if (data.L(t, a, data.hospitals[h1], data.bases[b]) == l1) {
                con += yt_ahb[t][a][h1][b];
              }
            }

            for (int l2 = 0; l2 < data.total_locals; ++l2) {
              if (data.L_tab[t + 1][a][b].find(l2) !=
                      data.L_tab[t + 1][a][b].end() &&
                  data.L(t, a, l2, data.bases[b]) == l1) {
                con += At_alb[t][a][l2][b];
                // for(int c = 0; c < data.types_call; ++c){
                // 	for(int l = 0; l < data.num_locals; ++l){
                // 		for(int h = 0; h < data.num_hosps; ++h){
                // 			if(data.H[c][l].find(h) !=
                // 				data.H[c][l].end() &&
                // 				data.A[c].find(a)!=data.A[c].end()){
                // 				//If using xt0, put t+1 on 1st
                // index 				con -=
                // xt_calblh[t][c][a][l2][b][l][h];
                // 			}
                // 		}
                // 	}
                // }
              }
            }

            name << "flow_locals_t" << t + 1 << "_" << a << "_" << b << "_"
                 << l1;
            con_location[t][a][l1][b] = model.addConstr(
                con, GRB_EQUAL, At_alb[t + 1][a][l1][b], name.str());
            name.str("");
            con = 0;
          }
        }
      }
    }
  }
}

void CGCall::add_queues_constraints() {
  GRBLinExpr con;
  std::stringstream name;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int c = 0; c < data.types_call; ++c) {
      for (int l = 0; l < data.num_locals; ++l) {
        con += Ct_cl[t][c][l];
        con += data.lambda[t][c][l];
        for (int a = 0; a < data.types_amb; ++a) {
          if (data.A[c].find(a) != data.A[c].end()) {
            for (int b = 0; b < data.num_bases; ++b) {
              for (int h = 0; h < data.num_hosps; ++h) {
                if (data.H[c][l].find(h) != data.H[c][l].end() &&
                    data.A[c].find(a) != data.A[c].end()) {
                  // If using xt0, put t+1 on 1st index
                  con -= xt_cablh[t][c][a][b][l][h];
                }
              }
            }
          }
        }

        for (int a = 0; a < data.types_amb; ++a) {
          if (data.A[c].find(a) != data.A[c].end()) {
            for (int h1 = 0; h1 < data.num_hosps; ++h1) {
              for (int h = 0; h < data.num_hosps; ++h) {
                if (data.H[c][l].find(h) != data.H[c][l].end() &&
                    data.A[c].find(a) != data.A[c].end()) {
                  // If using xt0, put t+1 on 1st index
                  con -= xt_cahlh[t][c][a][h1][l][h];
                }
              }
            }
          }
        }

        // for(int a = 0; a < data.types_amb; ++a){
        // 	if(data.A[c].find(a) != data.A[c].end()){
        // 		for(int b = 0; b < data.num_bases; ++b){
        // 			for(int l1 = 0; l1 < data.total_locals; ++l1){
        // 				for(int h = 0; h < data.num_hosps; ++h){
        // 					if(data.H[c][l].find(h) !=
        // data.H[c][l].end() &&
        // data.A[c].find(a) != data.A[c].end() &&
        // 						data.L_tab[t+1][a][b].find(l1)
        // != data.L_tab[t+1][a][b].end()){
        // 						//If using xt0, put t+1
        // on 1st index 						con -=
        // xt_calblh[t][c][a][l1][b][l][h];
        // 					}
        // 				}
        // 			}
        // 		}
        // 	}
        // }
        name << "flow_queue_t" << t + 1 << "_c" << c << "_l" << l;
        con_queue[t][c][l] =
            model.addConstr(con, GRB_EQUAL, Ct_cl[t + 1][c][l], name.str());
        name.str("");
        con = 0;
      }
    }
  }
}

void CGCall::add_ambs_location_constraints() {
  GRBLinExpr con = 0;
  std::stringstream name;

  // Location constraints
  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int b = 0; b < data.num_bases; ++b) {
        for (int l1 = 0; l1 < data.total_locals; ++l1) {
          if (data.L_tab[t + 1][a][b].find(l1) !=
              data.L_tab[t + 1][a][b].end()) {
            // for(int c = 0; c < data.types_call; ++c){
            // 	for(int l = 0; l < data.num_locals; ++l){
            // 		for(int h = 0; h < data.num_hosps; ++h){
            // 			if(data.H[c][l].find(h) != data.H[c][l].end() &&
            // 				data.A[c].find(a) != data.A[c].end() &&
            // 				data.L_tab[t][a][b].find(l1) !=
            // 				data.L_tab[t][a][b].end()){
            // 				//If using xt0, put t+1 on 1st index
            // 				con += xt_calblh[t][c][a][l1][b][l][h];
            // 			}
            // 		}
            // 	}
            // }
            name << "dispatch_locals_t" << t + 1 << "_a" << a << "_b";
            name << b << "_l" << l1;
            con_amb_location[t][a][l1][b] = model.addConstr(
                con, GRB_LESS_EQUAL, At_alb[t][a][l1][b], name.str());
            name.str("");
            con = 0;
          }
        }
      }
    }
  }

  // Max of ambulance at bases
  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int b = 0; b < data.num_bases; ++b) {
        for (int c = 0; c < data.types_call; ++c) {
          for (int l = 0; l < data.num_locals; ++l) {
            for (int h = 0; h < data.num_hosps; ++h) {
              if (data.H[c][l].find(h) != data.H[c][l].end() &&
                  data.A[c].find(a) != data.A[c].end()) {
                // If using xt0, put t+1 on 1st index
                con += xt_cablh[t][c][a][b][l][h];
              }
            }
          }
        }

        name << "dispatch_bases_t" << t + 1 << "_a" << a << "_b" << b;
        con_amb_base[t][a][b] =
            model.addConstr(con, GRB_LESS_EQUAL, At_ab[t][a][b], name.str());
        name.str("");
        con = 0;
      }
    }
  }
}

void CGCall::add_base_cap_constraints() {
  GRBLinExpr con;
  std::stringstream name;

  for (int t = 0; t < data.num_times; ++t) {
    for (int b = 0; b < data.num_bases; ++b) {
      for (int a = 0; a < data.types_amb; ++a) {
        con += At_ab[t][a][b];
      }
      name << "cap_bases_" << t + 1 << "_" << b;
      con_cap_base[t][b] =
          model.addConstr(con, GRB_LESS_EQUAL, data.cap_bases[b], name.str());
      name.str("");
      con = 0;
    }
  }
}

// int CGCall::get_location(Ambulance& amb){
// 	// Ambulance stopped
// 	if(amb.waiting_return || (!amb.busy && amb.base_ret == -1 &&
// !amb.waiting_return)){ 		return amb.location;
// 	}

// 	// Ambulance already at location
// 	if(amb.busy && amb.location == data.hospitals[amb.call_h_dest] ||
// 		!amb.busy && amb.base_ret != -1 && amb.location ==
// data.bases[amb.base_ret]){ 		return amb.location;
// 	}

// 	int departure_time = sim.get_time(amb.departure_time);
// 	int elapsed_time = (t0-departure_time)*data.quantum;
// 	int next = -1;
// 	std::vector<int> path;
// 	if(amb.busy){
// 		path = min_path(data, amb.origin, amb.call_location);
// 		std::vector<int> aux = min_path(data, amb.call_location,
// 			data.hospitals[amb.call_h_dest]);
// 		for(int i = 1; i < aux.size(); ++i){
// 			path.push_back(aux[i]);
// 		}
// 	}else{
// 		path = min_path(data, amb.origin, data.bases[amb.base_ret]);
// 	}

// 	if(path.size() == 1){
// 		next = path[0];
// 	}else if(path.size() == 2){
// 		if(data.times[path[0]][path[1]] < elapsed_time){
// 			next = path[1];
// 		}else{
// 			next = path[0];
// 		}
// 	}else{
// 		for(int i = 0; i < path.size()-1; ++i){
// 			if(data.times[path[i]][path[i+1]] < elapsed_time){
// 				next = path[i+1];
// 				elapsed_time -= data.times[path[i]][path[i+1]];
// 			}else{
// 				next = path[i];
// 				break;
// 			}
// 		}
// 	}
// 	return amb.location;
// }

// int CGCall::is_at_hospital(Ambulance & amb){
// 	for(int h = 0; h < data.num_hosps; ++h){
// 		if(amb.location == data.hospitals[h] && amb.waiting_return){
// 			return h;
// 		}
// 	}
// 	return -1;
// }

// int CGCall::is_at_base(Ambulance & amb){
// 	for(int b = 0; b < data.num_bases; ++b){
// 		if(amb.location == data.bases[b] && amb.base_ret == -1 &&
// !amb.waiting_return){ 			return b;
// 		}
// 	}
// 	return -1;
// }

// std::pair<int,int> CGCall::is_at_location_base(Ambulance & amb){
// 	if(amb.base_ret < 0)
// 		return std::make_pair(-1,-1);
// 	else
// 		return std::make_pair(amb.location,amb.base_ret);
// }

int CGCall::c_tl(int t, int l) {
  int arg_min = -1;
  double min_val = GRB_INFINITY;
  // std::cout << "phi_tcl " << t+1 << " " << l << ":\n";
  for (int c = 0; c < data.types_call; ++c) {
    // std::cout << c << ": phi = " << phi_tcl[t][c][l] << "\n";
    if (phi_tcl[t][c][l] < min_val) {
      min_val = phi_tcl[t][c][l];
      arg_min = c;
    }
  }
  // std::cout << "min = " << arg_min << "\n";
  // std::cin.get();
  return arg_min;
}

int CGCall::b_tal(int t, int a, int l) {
  int arg_min = -1;
  double min_val = GRB_INFINITY;
  // std::cout << "b_tal " << t+1 << " " << a << " " << l << ":\n";
  for (int b = 0; b < data.num_bases; ++b) {
    if (data.L_tab[t + 1][a][b].find(l) != data.L_tab[t + 1][a][b].end()) {
      double sum = 0;
      int next_l = data.L(t + 1, a, l, data.bases[b]);
      // std::cout << b << ": ";
      if (next_l == data.bases[b]) {
        sum += beta_tab[t][a][b];
        // std::cout << "beta = " << beta_tab[t][a][b] << " ";
      }

      if (data.L_tab[t + 1][a][b].find(next_l) !=
          data.L_tab[t + 1][a][b].end()) {
        sum += psi_talb[t][a][next_l][b];
        // std::cout << ", psi = " << psi_talb[t][a][next_l][b] << " ";
      }

      sum += theta_talb[t][a][l][b];
      // std::cout << ", theta = " << theta_talb[t][a][l][b] << ", ";
      // std::cout << "sum = " << sum << "\n";
      if (sum < min_val) {
        min_val = sum;
        arg_min = b;
      }
    }
  }
  // std::cout << "min = " << arg_min << "\n";
  // std::cin.get();
  return arg_min;
}

int CGCall::h_tall(int t, int a, int l1, int l) {
  int arg_min = -1;
  double min_val = GRB_INFINITY;
  int factor = 1;
  int current_time = (t + 1) * data.quantum;
  // std::cout << "h_tall " << t+1 << " "  << a << " " << l1 << " " << l <<
  // ":\n";
  for (auto h : data.H[0][l]) {
    int f_val = f(t + 1, 0, a, l1, l, h, data, factor);
    // std::cout << h << ": f = " << f_val << ", arr_time = ";
    int arrival_time =
        get_time(data, current_time + data.tao[t][0][a][l1][l][h]);
    double sum = f_val;
    if (arrival_time - 1 < data.num_times - 1) {
      // std::cout << ", alpha = " << alpha_tah[arrival_time-1][a][h] << "\n";
      sum -= alpha_tah[arrival_time - 1][a][h];
    } else {
      // std::cout << "\n";
    }

    if (sum < min_val) {
      min_val = sum;
      arg_min = h;
    }
  }
  // std::cout << "min = " << arg_min << "\n";
  // std::cin.get();
  return arg_min;
}

CGCall::~CGCall() {
  for (int a = 0; a < data.types_amb; ++a) {
    for (int l1 = 0; l1 < data.total_locals; ++l1) {
      delete[] A0_alb[a][l1];
    }
    delete[] A0_ab[a];
    delete[] A0_ah[a];
    delete[] A0_alb[a];
  }

  delete[] A0_ab;
  delete[] A0_ah;
  delete[] A0_alb;

  for (int c = 0; c < data.types_call; ++c) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        delete[] A0_calh[c][a][l1];
      }
      delete[] A0_calh[c][a];
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        for (int l = 0; l < data.num_locals; ++l) {
          delete[] A0_callh[c][a][l1][l];
        }
        delete[] A0_callh[c][a][l1];
      }
      delete[] A0_callh[c][a];
    }
    delete[] A0_calh[c];
    delete[] A0_callh[c];
    delete[] C[c];
  }

  delete[] A0_callh;
  delete[] A0_calh;
  delete[] C;

  for (int a = 0; a < data.types_amb; ++a) {
    for (int b = 0; b < data.num_bases; ++b) {
      delete[] x0_abh[a][b];
    }
    delete[] x0_abh[a];

    for (int h1 = 0; h1 < data.num_hosps; ++h1) {
      delete[] x0_ahh[a][h1];
    }
    delete[] x0_ahh[a];

    for (int l = 0; l < data.total_locals; ++l) {
      for (int b = 0; b < data.num_bases; ++b) {
        delete[] x0_albh[a][l][b];
      }
      delete[] x0_albh[a][l];
    }
    delete[] x0_albh[a];

    for (int h = 0; h < data.num_hosps; ++h) {
      delete[] y0_ahb[a][h];
    }
    delete[] y0_ahb[a];
  }
  delete[] x0_abh;
  delete[] x0_ahh;
  delete[] x0_albh;
  delete[] y0_ahb;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int c = 0; c < data.types_call; ++c) {
      for (int a = 0; a < data.types_amb; ++a) {
        for (int b = 0; b < data.num_bases; ++b) {
          for (int l = 0; l < data.num_locals; ++l) {
            delete[] xt_cablh[t][c][a][b][l];
          }
          delete[] xt_cablh[t][c][a][b];
        }
        delete[] xt_cablh[t][c][a];

        for (int h1 = 0; h1 < data.num_hosps; ++h1) {
          for (int l = 0; l < data.num_locals; ++l) {
            delete[] xt_cahlh[t][c][a][h1][l];
          }
          delete[] xt_cahlh[t][c][a][h1];
        }
        delete[] xt_cahlh[t][c][a];

        for (int l1 = 0; l1 < data.total_locals; ++l1) {
          for (int b = 0; b < data.num_bases; ++b) {
            for (int l = 0; l < data.num_locals; ++l) {
              delete[] xt_calblh[t][c][a][l1][b][l];
            }
            delete[] xt_calblh[t][c][a][l1][b];
          }
          delete[] xt_calblh[t][c][a][l1];
        }
        delete[] xt_calblh[t][c][a];
      }
      delete[] xt_cablh[t][c];
      delete[] xt_cahlh[t][c];
      delete[] xt_calblh[t][c];
    }
    delete[] xt_cablh[t];
    delete[] xt_cahlh[t];
    delete[] xt_calblh[t];

    for (int a = 0; a < data.types_amb; ++a) {
      for (int h = 0; h < data.num_hosps; ++h) {
        delete[] yt_ahb[t][a][h];
      }
      delete[] yt_ahb[t][a];
    }
    delete[] yt_ahb[t];
  }
  delete[] xt_cablh;
  delete[] xt_cahlh;
  delete[] xt_calblh;
  delete[] yt_ahb;

  for (int t = 0; t < data.num_times; ++t) {
    for (int c = 0; c < data.types_call; ++c) {
      delete[] Ct_cl[t][c];
    }
    delete[] Ct_cl[t];
    for (int a = 0; a < data.types_amb; ++a) {
      delete[] At_ab[t][a];
    }
    delete[] At_ab[t];

    for (int a = 0; a < data.types_amb; ++a) {
      for (int l = 0; l < data.total_locals; ++l) {
        delete[] At_alb[t][a][l];
      }
      delete[] At_alb[t][a];
    }
    delete[] At_alb[t];
  }
  delete[] Ct_cl;
  delete[] At_ab;
  delete[] At_alb;

  for (int a = 0; a < data.types_amb; ++a) {
    delete[] con_base_t0[a];
  }
  delete[] con_base_t0;

  for (int a = 0; a < data.types_amb; ++a) {
    delete[] con_hospital_t0[a];
  }
  delete[] con_hospital_t0;

  for (int a = 0; a < data.types_amb; ++a) {
    for (int l1 = 0; l1 < data.total_locals; ++l1) {
      delete[] con_location_t0[a][l1];
    }
    delete[] con_location_t0[a];
  }
  delete[] con_location_t0;

  for (int c = 0; c < data.types_call; ++c) {
    delete[] con_queue_t0[c];
  }
  delete[] con_queue_t0;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      delete[] con_base[t][a];
    }
    delete[] con_base[t];
  }
  delete[] con_base;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      delete[] con_hospital[t][a];
    }
    delete[] con_hospital[t];
  }
  delete[] con_hospital;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        delete[] con_location[t][a][l1];
      }
      delete[] con_location[t][a];
    }
    delete[] con_location[t];
  }
  delete[] con_location;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int c = 0; c < data.types_call; ++c) {
      delete[] con_queue[t][c];
    }
    delete[] con_queue[t];
  }
  delete[] con_queue;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        delete[] con_amb_location[t][a][l1];
      }
      delete[] con_amb_location[t][a];
    }
    delete[] con_amb_location[t];
  }
  delete[] con_amb_location;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      delete[] con_amb_base[t][a];
    }
    delete[] con_amb_base[t];
  }
  delete[] con_amb_base;

  for (int t = 0; t < data.num_times; ++t) {
    delete[] con_cap_base[t];
  }
  delete[] con_cap_base;

  for (int a = 0; a < data.types_amb; ++a) {
    delete[] beta0_ab[a];
    delete[] alpha0_ah[a];
    for (int l1 = 0; l1 < data.total_locals; ++l1) {
      delete[] psi0_alb[a][l1];
    }
    delete[] psi0_alb[a];
  }
  delete[] beta0_ab;
  delete[] alpha0_ah;
  delete[] psi0_alb;

  for (int c = 0; c < data.types_call; ++c) {
    delete[] phi0_cl[c];
  }
  delete[] phi0_cl;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        delete[] psi_talb[t][a][l1];
        delete[] theta_talb[t][a][l1];
      }
      delete[] alpha_tah[t][a];
      delete[] psi_talb[t][a];
      delete[] beta_tab[t][a];
      delete[] theta_talb[t][a];
      delete[] gamma_tab[t][a];
    }

    for (int c = 0; c < data.types_call; ++c) {
      delete[] phi_tcl[t][c];
    }
    delete[] beta_tab[t];
    delete[] alpha_tah[t];
    delete[] psi_talb[t];
    delete[] phi_tcl[t];
    delete[] ups_tb[t];
    delete[] theta_talb[t];
    delete[] gamma_tab[t];
  }
  delete[] beta_tab;
  delete[] alpha_tah;
  delete[] psi_talb;
  delete[] phi_tcl;
  delete[] ups_tb;
  delete[] theta_talb;
  delete[] gamma_tab;

  for (int t = 0; t < data.num_times; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int l1 = 0; l1 < data.num_locals; ++l1) {
        delete[] vh_tall[t][a][l1];
      }
      delete[] vh_tall[t][a];
      delete[] vb_tal[t][a];
    }

    delete[] vh_tall[t];
    delete[] vb_tal[t];
    delete[] vc_tl[t];
  }
  delete[] vh_tall;
  delete[] vc_tl;
  delete[] vb_tal;
}

//==============================================================================================

CGAmbulance::CGAmbulance(Data& data, GRBEnv& env)
    : data(data), env(env), model(env), columns_added(0) {
  std::stringstream name;

  std::vector<std::pair<int, int>> intermediate;
  intermediate.reserve(data.total_locals);
  closest_l1_to_l.reserve(data.num_locals);
  // std::cout << "closest_l1_to_l:\n";
  Travel& travel = data.travel;
  for (int l = 0; l < data.num_locals; ++l) {
    closest_l1_to_l.push_back(std::vector<int>());
    closest_l1_to_l[l].reserve(data.total_locals);
    intermediate.clear();
    for (int l1 = 0; l1 < data.num_locals; ++l1) {
      auto dist_l1_l =
          travel.lat_long_distance(data.locations[l1], data.locations[l]);
      intermediate.push_back(std::make_pair(dist_l1_l, l1));
    }
    std::sort(intermediate.begin(), intermediate.end());
    // std::cout << l << ": ";
    for (auto& elem : intermediate) {
      // std::cout << elem.second << " ";
      closest_l1_to_l[l].push_back(elem.second);
    }
    // std::cout << "\n";
  }
  // std::cin.get();

  vh_tall = new int***[data.num_times];
  vc_tl = new int*[data.num_times];
  vb_tal = new int**[data.num_times];
  for (int t = 0; t < data.num_times; ++t) {
    vh_tall[t] = new int**[data.types_amb];
    vb_tal[t] = new int*[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      vh_tall[t][a] = new int*[data.num_locals];
      vb_tal[t][a] = new int[data.num_locals];
      for (int l1 = 0; l1 < data.num_locals; ++l1) {
        vh_tall[t][a][l1] = new int[data.num_locals];
        for (int l = 0; l < data.num_locals; ++l) {
          vh_tall[t][a][l1][l] = -1;
        }

        vb_tal[t][a][l1] = -1;
      }
    }

    vc_tl[t] = new int[data.num_locals];
    for (int l = 0; l < data.num_locals; ++l) {
      vc_tl[t][l] = -1;
    }
  }

  A0_ab = new int*[data.types_amb];
  A0_ah = new int*[data.types_amb];
  A0_alb = new int**[data.types_amb];
  A0_calh = new int***[data.types_call];
  A0_callh = new int****[data.types_call];
  C = new int*[data.types_call];

  for (int a = 0; a < data.types_amb; ++a) {
    A0_ab[a] = new int[data.num_bases];
    for (int b = 0; b < data.num_bases; ++b) {
      A0_ab[a][b] = data.A0_ab[a][b];
    }
    A0_ah[a] = new int[data.num_hosps];
    for (int h = 0; h < data.num_hosps; ++h) {
      A0_ah[a][h] = data.A0_ah[a][h];
    }
    A0_alb[a] = new int*[data.total_locals];
    for (int l1 = 0; l1 < data.total_locals; ++l1) {
      A0_alb[a][l1] = new int[data.num_bases];
      for (int b = 0; b < data.num_bases; ++b) {
        A0_alb[a][l1][b] = data.A0_alb[a][l1][b];
      }
    }
  }

  for (int c = 0; c < data.types_call; ++c) {
    A0_calh[c] = new int**[data.types_amb];
    A0_callh[c] = new int***[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      A0_calh[c][a] = new int*[data.total_locals];
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        A0_calh[c][a][l1] = new int[data.num_hosps];
        for (int h = 0; h < data.num_hosps; ++h) {
          A0_calh[c][a][l1][h] = data.A0_calh[c][a][l1][h];
        }
      }
      A0_callh[c][a] = new int**[data.total_locals];
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        A0_callh[c][a][l1] = new int*[data.num_locals];
        for (int l = 0; l < data.num_locals; ++l) {
          A0_callh[c][a][l1][l] = new int[data.num_hosps];
          for (int h = 0; h < data.num_hosps; ++h) {
            A0_callh[c][a][l1][l][h] = data.A0_callh[c][a][l1][l][h];
          }
        }
      }
    }
    C[c] = new int[data.num_locals];
    for (int l = 0; l < data.num_locals; ++l) {
      C[c][l] = data.C[c][l];
    }
  }

  y0_ahb = new GRBVar**[data.types_amb];
  for (int a = 0; a < data.types_amb; ++a) {
    y0_ahb[a] = new GRBVar*[data.num_hosps];
    for (int h = 0; h < data.num_hosps; ++h) {
      y0_ahb[a][h] = new GRBVar[data.num_bases];
      for (int b = 0; b < data.num_bases; ++b) {
        name << "y0" << "_a" << a << "_h" << h << "_b" << b;
        y0_ahb[a][h][b] =
            model.addVar(0, A0_ah[a][h], 0, GRB_CONTINUOUS, name.str());
        name.str("");
      }
    }
  }

  y0_ahclh = new GRBVar****[data.types_amb];
  for (int a = 0; a < data.types_amb; ++a) {
    y0_ahclh[a] = new GRBVar***[data.num_hosps];
    for (int h1 = 0; h1 < data.num_hosps; ++h1) {
      y0_ahclh[a][h1] = new GRBVar**[data.types_call];
      for (int c = 0; c < data.types_call; ++c) {
        y0_ahclh[a][h1][c] = new GRBVar*[data.num_locals];
        for (int l = 0; l < data.num_locals; ++l) {
          y0_ahclh[a][h1][c][l] = new GRBVar[data.num_hosps];
          for (int h = 0; h < data.num_hosps; ++h) {
            if (data.A[c].find(a) != data.A[c].end() &&
                data.H[c][l].find(h) != data.H[c][l].end()) {
              name << "y0_a" << a << "_h'" << h1 << "_c" << c << "_l" << l
                   << "_h" << h;
              y0_ahclh[a][h1][c][l][h] =
                  model.addVar(0, A0_ah[a][h], 0, GRB_CONTINUOUS, name.str());
              name.str("");
            }
          }
        }
      }
    }
  }

  xt_cablh = new GRBVar*****[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    xt_cablh[t] = new GRBVar****[data.types_call];
    for (int c = 0; c < data.types_call; ++c) {
      xt_cablh[t][c] = new GRBVar***[data.types_amb];
      for (int a = 0; a < data.types_amb; ++a) {
        xt_cablh[t][c][a] = new GRBVar**[data.num_bases];
        for (int b = 0; b < data.num_bases; ++b) {
          xt_cablh[t][c][a][b] = new GRBVar*[data.num_locals];
          for (int l = 0; l < data.num_locals; ++l) {
            xt_cablh[t][c][a][b][l] = new GRBVar[data.num_hosps];
            for (int h = 0; h < data.num_hosps; ++h) {
              if (data.A[c].find(a) != data.A[c].end() &&
                  data.H[c][l].find(h) != data.H[c][l].end()) {
                name << "xt" << t + 1 << "_c" << c << "_a" << a << "_b" << b;
                name << "_l" << l << "_h" << h;
                xt_cablh[t][c][a][b][l][h] = model.addVar(
                    0, data.num_ambs, 0, GRB_CONTINUOUS, name.str());
                name.str("");
              }
            }
          }
        }
      }
    }
  }

  xt_cahlh = new GRBVar*****[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    xt_cahlh[t] = new GRBVar****[data.types_call];
    for (int c = 0; c < data.types_call; ++c) {
      xt_cahlh[t][c] = new GRBVar***[data.types_amb];
      for (int a = 0; a < data.types_amb; ++a) {
        xt_cahlh[t][c][a] = new GRBVar**[data.num_hosps];
        for (int h1 = 0; h1 < data.num_hosps; ++h1) {
          xt_cahlh[t][c][a][h1] = new GRBVar*[data.num_locals];
          for (int l = 0; l < data.num_locals; ++l) {
            xt_cahlh[t][c][a][h1][l] = new GRBVar[data.num_hosps];
            for (int h = 0; h < data.num_hosps; ++h) {
              if (data.A[c].find(a) != data.A[c].end() &&
                  data.H[c][l].find(h) != data.H[c][l].end()) {
                name << "xt" << t + 1 << "_c" << c << "_a" << a << "_h'" << h1;
                name << "_l" << l << "_h" << h;
                xt_cahlh[t][c][a][h1][l][h] = model.addVar(
                    0, data.num_ambs, 0, GRB_CONTINUOUS, name.str());
                name.str("");
              }
            }
          }
        }
      }
    }
  }

  xt_calblh = new GRBVar******[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    xt_calblh[t] = new GRBVar*****[data.types_call];
    for (int c = 0; c < data.types_call; ++c) {
      xt_calblh[t][c] = new GRBVar****[data.types_amb];
      for (int a = 0; a < data.types_amb; ++a) {
        xt_calblh[t][c][a] = new GRBVar***[data.total_locals];
        for (int l1 = 0; l1 < data.total_locals; ++l1) {
          xt_calblh[t][c][a][l1] = new GRBVar**[data.num_bases];
          for (int b = 0; b < data.num_bases; ++b) {
            xt_calblh[t][c][a][l1][b] = new GRBVar*[data.num_locals];
            for (int l = 0; l < data.num_locals; ++l) {
              xt_calblh[t][c][a][l1][b][l] = new GRBVar[data.num_hosps];
              // for(int h = 0; h < data.num_hosps; ++h){
              // 	if(data.A[c].find(a) != data.A[c].end() &&
              // 		data.H[c][l].find(h) != data.H[c][l].end()
              // 		&& data.L_tab[t+1][a][b].find(l1) !=
              // 		data.L_tab[t+1][a][b].end()){
              // 		name << "xt" << t+1 << "_c" << c << "_a" << a <<
              // "_l'"; 		name << l1 << "_b" << b << "_l" << l <<
              // "_h" << h; 		xt_calblh[t][c][a][l1][b][l][h] =
              // model.addVar( 			0,data.num_ambs, 0,
              // (data.relax["xt_calblh"] ? GRB_CONTINUOUS : GRB_INTEGER),
              // name.str()); 		name.str("");
              // 	}
              // }
            }
          }
        }
      }
    }
  }

  yt_ahb = new GRBVar***[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    yt_ahb[t] = new GRBVar**[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      yt_ahb[t][a] = new GRBVar*[data.num_hosps];
      for (int h = 0; h < data.num_hosps; ++h) {
        yt_ahb[t][a][h] = new GRBVar[data.num_bases];
        for (int b = 0; b < data.num_bases; ++b) {
          name << "yt" << t + 1 << "_a" << a << "_h" << h << "_b" << b;
          yt_ahb[t][a][h][b] =
              model.addVar(0, data.num_ambs, 0, GRB_CONTINUOUS, name.str());
          name.str("");
        }
      }
    }
  }

  Ct_cl = new GRBVar**[data.num_times];
  for (int t = 0; t < data.num_times; ++t) {
    Ct_cl[t] = new GRBVar*[data.types_call];
    for (int c = 0; c < data.types_call; ++c) {
      Ct_cl[t][c] = new GRBVar[data.num_locals];
      for (int l = 0; l < data.num_locals; ++l) {
        name << "Ct" << t + 1 << "_c" << c << "_l" << l;
        Ct_cl[t][c][l] =
            model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name.str());
        name.str("");
      }
    }
  }

  At_ab = new GRBVar**[data.num_times];
  for (int t = 0; t < data.num_times; ++t) {
    At_ab[t] = new GRBVar*[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      At_ab[t][a] = new GRBVar[data.num_bases];
      for (int b = 0; b < data.num_bases; ++b) {
        name << "At" << t + 1 << "_a" << a << "_b" << b;
        At_ab[t][a][b] =
            model.addVar(0, data.num_ambs, 0, GRB_CONTINUOUS, name.str());
        name.str("");
      }
    }
  }

  At_alb = new GRBVar***[data.num_times];
  for (int t = 0; t < data.num_times; ++t) {
    At_alb[t] = new GRBVar**[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      At_alb[t][a] = new GRBVar*[data.total_locals];
      for (int l = 0; l < data.total_locals; ++l) {
        At_alb[t][a][l] = new GRBVar[data.num_bases];
        for (int b = 0; b < data.num_bases; ++b) {
          if (data.L_tab[t][a][b].find(l) == data.L_tab[t][a][b].end())
            continue;
          name << "At" << t + 1 << "_a" << a << "_l" << l << "_b" << b;
          At_alb[t][a][l][b] =
              model.addVar(0, data.num_ambs, 0, GRB_CONTINUOUS, name.str());
          name.str("");
        }
      }
    }
  }

  GRBLinExpr fo;
  double factor = 1;

  for (int a = 0; a < data.types_amb; ++a) {
    for (int h = 0; h < data.num_hosps; ++h) {
      // for(int b = 0; b < data.num_bases; ++b){
      // 	fo +=
      // data.times[data.hospitals[h]][data.bases[b]]*y0_ahb[a][h][b];
      // }

      for (int c = 0; c < data.types_call; ++c) {
        for (int l = 0; l < data.num_locals; ++l) {
          for (int h1 = 0; h1 < data.num_hosps; ++h1) {
            if (data.A[c].find(a) != data.A[c].end() &&
                data.H[c][l].find(h1) != data.H[c][l].end()) {
              fo += f(0, c, a, data.hospitals[h], l, h1, data, factor) *
                    y0_ahclh[a][h][c][l][h1];
            }
          }
        }
      }
    }
  }

  for (int t = 1; t < data.num_times; ++t) {
    for (int c = 0; c < data.types_call; ++c) {
      for (int a = 0; a < data.types_amb; ++a) {
        if (data.A[c].find(a) != data.A[c].end()) {
          for (int l = 0; l < data.num_locals; ++l) {
            for (int h = 0; h < data.num_hosps; ++h) {
              if (data.H[c][l].find(h) != data.H[c][l].end()) {
                for (int b = 0; b < data.num_bases; ++b) {
                  fo += f(t - 1, c, a, data.bases[b], l, h, data, factor) *
                        xt_cablh[t - 1][c][a][b][l][h];
                }

                for (int h1 = 0; h1 < data.num_hosps; ++h1) {
                  fo += f(t - 1, c, a, data.hospitals[h1], l, h, data, factor) *
                        xt_cahlh[t - 1][c][a][h1][l][h];
                }

                // for(int l1 = 0; l1 < data.total_locals; ++l1){
                // 	for(int b = 0; b < data.num_bases; ++b){
                // 		if(data.L_tab[t][a][b].find(l1) !=
                // 			data.L_tab[t][a][b].end()){
                // 			fo += f(t-1,c,a,l1,data.bases[b],
                // 			l, h, data,
                // 			factor)*
                // 			xt_calblh[t-1][c][a][l1][b][l][h];
                // 		}
                // 	}
                // }
              }
            }
          }
        }
      }
    }

    // for(int a = 0; a < data.types_amb; ++a){
    // 	for(int h = 0; h < data.num_hosps; ++h){
    // 		for(int b = 0; b < data.num_bases; ++b){
    // 			fo += data.times[data.hospitals[h]][data.bases[b]]*
    // 				yt_ahb[t-1][a][h][b];
    // 		}
    // 	}
    // }
  }

  for (int t = 0; t < data.num_times; ++t) {
    for (int c = 0; c < data.types_call; ++c) {
      factor = data.ins.penalties[c] * 36000;
      for (int l = 0; l < data.num_locals; ++l) {
        fo += g_tcl(t + 1, c, l, factor) * Ct_cl[t][c][l];
      }
    }
    factor = 0;
    for (int a = 0; a < data.types_amb; ++a) {
      for (int b = 0; b < data.num_bases; ++b) {
        fo += g_tab(t + 1, a, b, factor) * At_ab[t][a][b];
        for (int l1 = 0; l1 < data.total_locals; ++l1) {
          if (data.L_tab[t][a][b].find(l1) != data.L_tab[t][a][b].end()) {
            fo += g_talb(t + 1, a, l1, data.bases[b], data, factor) *
                  At_alb[t][a][l1][b];
          }
        }
      }
    }
  }

  model.setObjective(fo, GRB_MINIMIZE);

  con_base_t0 = new GRBConstr*[data.types_amb];
  for (int a = 0; a < data.types_amb; ++a) {
    con_base_t0[a] = new GRBConstr[data.num_bases];
  }

  con_hospital_t0 = new GRBConstr*[data.types_amb];
  for (int a = 0; a < data.types_amb; ++a) {
    con_hospital_t0[a] = new GRBConstr[data.num_hosps];
  }
  con_location_t0 = new GRBConstr**[data.types_amb];
  for (int a = 0; a < data.types_amb; ++a) {
    con_location_t0[a] = new GRBConstr*[data.total_locals];
    for (int l1 = 0; l1 < data.total_locals; ++l1) {
      con_location_t0[a][l1] = new GRBConstr[data.num_bases];
    }
  }
  con_queue_t0 = new GRBConstr*[data.types_call];
  for (int c = 0; c < data.types_call; ++c) {
    con_queue_t0[c] = new GRBConstr[data.num_locals];
  }

  con_base = new GRBConstr**[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    con_base[t] = new GRBConstr*[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      con_base[t][a] = new GRBConstr[data.num_bases];
    }
  }

  con_hospital = new GRBConstr**[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    con_hospital[t] = new GRBConstr*[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      con_hospital[t][a] = new GRBConstr[data.num_hosps];
    }
  }

  con_location = new GRBConstr***[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    con_location[t] = new GRBConstr**[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      con_location[t][a] = new GRBConstr*[data.total_locals];
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        con_location[t][a][l1] = new GRBConstr[data.num_bases];
      }
    }
  }

  con_queue = new GRBConstr**[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    con_queue[t] = new GRBConstr*[data.types_call];
    for (int c = 0; c < data.types_call; ++c) {
      con_queue[t][c] = new GRBConstr[data.num_locals];
    }
  }

  con_amb_location = new GRBConstr***[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    con_amb_location[t] = new GRBConstr**[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      con_amb_location[t][a] = new GRBConstr*[data.total_locals];
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        con_amb_location[t][a][l1] = new GRBConstr[data.num_bases];
      }
    }
  }

  con_amb_base = new GRBConstr**[data.num_times - 1];
  for (int t = 0; t < data.num_times - 1; ++t) {
    con_amb_base[t] = new GRBConstr*[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      con_amb_base[t][a] = new GRBConstr[data.num_bases];
    }
  }

  con_cap_base = new GRBConstr*[data.num_times];
  for (int t = 0; t < data.num_times; ++t) {
    con_cap_base[t] = new GRBConstr[data.num_bases];
  }

  beta0_ab = new double*[data.types_amb];
  alpha0_ah = new double*[data.types_amb];
  psi0_alb = new double**[data.types_amb];
  for (int a = 0; a < data.types_amb; ++a) {
    beta0_ab[a] = new double[data.num_bases];
    for (int b = 0; b < data.num_bases; ++b) {
      beta0_ab[a][b] = 0;
    }

    alpha0_ah[a] = new double[data.num_hosps];
    for (int h = 0; h < data.num_hosps; ++h) {
      alpha0_ah[a][h] = 0;
    }

    psi0_alb[a] = new double*[data.total_locals];
    for (int l1 = 0; l1 < data.total_locals; ++l1) {
      psi0_alb[a][l1] = new double[data.num_bases];
      for (int b = 0; b < data.num_bases; ++b) {
        psi0_alb[a][l1][b] = 0;
      }
    }
  }

  phi0_cl = new double*[data.types_call];
  for (int c = 0; c < data.types_call; ++c) {
    phi0_cl[c] = new double[data.num_locals];
    for (int l = 0; l < data.num_locals; ++l) {
      phi0_cl[c][l] = 0;
    }
  }

  beta_tab = new double**[data.num_times - 1];
  alpha_tah = new double**[data.num_times - 1];
  psi_talb = new double***[data.num_times - 1];
  phi_tcl = new double**[data.num_times - 1];
  ups_tb = new double*[data.num_times - 1];
  theta_talb = new double***[data.num_times - 1];
  gamma_tab = new double**[data.num_times];

  for (int t = 0; t < data.num_times - 1; ++t) {
    beta_tab[t] = new double*[data.types_amb];
    alpha_tah[t] = new double*[data.types_amb];
    psi_talb[t] = new double**[data.types_amb];
    theta_talb[t] = new double**[data.types_amb];
    gamma_tab[t] = new double*[data.types_amb];
    for (int a = 0; a < data.types_amb; ++a) {
      alpha_tah[t][a] = new double[data.num_hosps];
      psi_talb[t][a] = new double*[data.total_locals];
      theta_talb[t][a] = new double*[data.total_locals];
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        psi_talb[t][a][l1] = new double[data.num_bases];
        theta_talb[t][a][l1] = new double[data.num_bases];
        for (int b = 0; b < data.num_bases; ++b) {
          psi_talb[t][a][l1][b] = 0;
          theta_talb[t][a][l1][b] = 0;
        }
      }
      beta_tab[t][a] = new double[data.num_bases];
      gamma_tab[t][a] = new double[data.num_bases];

      for (int b = 0; b < data.num_bases; ++b) {
        gamma_tab[t][a][b] = 0;
        beta_tab[t][a][b] = 0;
      }
    }
    phi_tcl[t] = new double*[data.types_call];
    for (int c = 0; c < data.types_call; ++c) {
      phi_tcl[t][c] = new double[data.num_locals];
      for (int l = 0; l < data.num_locals; ++l) {
        phi_tcl[t][c][l] = 0;
      }
    }
    ups_tb[t] = new double[data.num_bases];
    for (int b = 0; b < data.num_bases; ++b) {
      ups_tb[t][b] = 0;
    }
  }
  gamma_tab[data.num_times - 1] = new double*[data.types_amb];
  for (int a = 0; a < data.types_amb; ++a) {
    gamma_tab[data.num_times - 1][a] = new double[data.num_bases];
    for (int b = 0; b < data.num_bases; ++b) {
      gamma_tab[data.num_times - 1][a][b] = 0;
    }
  }

  add_bases_constraints_t0();
  add_hospitals_constraints_t0();
  add_locations_constraints_t0();
  add_queues_constraints_t0();

  add_bases_constraints();
  add_hospitals_constraints();
  add_locations_constraints();
  add_queues_constraints();
  add_ambs_location_constraints();
  add_base_cap_constraints();

  model.update();
}

void CGAmbulance::solve() {
  iter = 0;
  model.set(GRB_IntParam_OutputFlag, 0);
  model.set(GRB_IntParam_Presolve, 0);
  model.set(GRB_IntParam_Threads, 1);
  bool column_found = false;

  // std::cout << "Number of Variables: " << model.get(GRB_IntAttr_NumVars) <<
  // "\n";
  double sum_master_time = 0;
  double sum_pricing_time = 0;
  do {
    model.update();
    model.optimize();
    sum_master_time += model.get(GRB_DoubleAttr_Runtime);
    // double z =  model.get(GRB_DoubleAttr_ObjVal);
    // std::cout << "Iter " << iter << ", Z = " << z << "\n";
    set_dual();
    std::chrono::high_resolution_clock::time_point t0 =
        std::chrono::high_resolution_clock::now();
    column_found = pricing();
    std::chrono::high_resolution_clock::time_point dt =
        std::chrono::high_resolution_clock::now();
    sum_pricing_time +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(dt - t0).count() /
        pow(10, 9);
    ++iter;
  } while (column_found);

  // std::cout << "Columns Added: " << columns_added << "\n";
  // std::cout << "Average Pricing Time: " << sum_pricing_time/iter << "\n";
  // std::cout << "Average Master Time: " << sum_master_time/iter << "\n";
  // model.optimize();
  // model.write("cg_amb.lp");

  if (model.get(GRB_IntAttr_Status) != 2) {
    fmt::print("Status {}\n", model.get(GRB_IntAttr_Status));
    fmt::print("Num_iter: {}\n", iter);
    for (int a = 0; a < data.types_amb; ++a) {
      for (int h = 0; h < data.num_hosps; ++h) {
        if (data.A0_ah[a][h] > 0) {
          fmt::print("a = {}, h = {}, A0_ah = {}\n", a, h, data.A0_ah[a][h]);
        }
      }
    }
    for (int a = 0; a < data.types_amb; ++a) {
      for (int b = 0; b < data.num_bases; ++b) {
        if (data.A0_ab[a][b] > 0) {
          fmt::print("a = {}, b = {}, A0_ab = {}\n", a, b, data.A0_ab[a][b]);
        }
      }
    }

    for (int a = 0; a < data.types_amb; ++a) {
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        for (int b = 0; b < data.num_bases; ++b) {
          if (data.A0_alb[a][l1][b]) {
            fmt::print("a = {}, l = {}, b = {}, A0_alb = {}\n", a, l1, b,
                       data.A0_alb[a][l1][b]);
          }
        }
      }
    }

    for (int c = 0; c < data.types_call; ++c) {
      for (int a = 0; a < data.types_amb; ++a) {
        for (int l1 = 0; l1 < data.total_locals; ++l1) {
          for (int l = 0; l < data.num_locals; ++l) {
            for (int h = 0; h < data.num_hosps; ++h) {
              if (data.A0_callh[c][a][l1][l][h] > 0) {
                fmt::print(
                    "c = {}, a = {}, l1 = {}, l = {}, h = {}, A0_callh = {}\n",
                    c, a, l1, l, h, data.A0_callh[c][a][l1][l][h]);
              }
            }
          }
        }
      }
    }

    for (int c = 0; c < data.types_call; ++c) {
      for (int a = 0; a < data.types_amb; ++a) {
        for (int l1 = 0; l1 < data.total_locals; ++l1) {
          for (int h = 0; h < data.num_hosps; ++h) {
            if (data.A0_calh[c][a][l1][h] > 0) {
              fmt::print("c = {}, a = {}, l1 = {}, h = {}, A0_calh = {}\n", c,
                         a, l1, h, data.A0_calh[c][a][l1][h]);
            }
          }
        }
      }
    }

    for (int c = 0; c < data.types_call; ++c) {
      for (int l = 0; l < data.num_locals; ++l) {
        if (data.C[c][l] > 0) {
          fmt::print("c = {}, l = {}, C0_cl = {}\n");
        }
      }
    }

    for (auto& amb : data.ambulances) {
      cout << amb << "\n";
    }

    fmt::print("Queue = {}\n", data.queue);

    for (auto& call : data.calls) {
      cout << call << "\n";
    }

    fmt::print("Time = {}, quantum = {}\n", data.time, data.quantum);

    model.computeIIS();
    model.write("cg_amb.lp");
    model.write("cg_amb.ilp");
    std::cin.get();
  }

  // std::cin.get();

  amb_type = -1;
  src_hosp = -1;
  base_return = -1;

  call_type = -1;
  call_location = -1;
  dst_hosp = -1;

  // bool found = false;

  if (false) {
    for (auto& call : data.calls) {
      cout << call << "\n";
    }
    print_solution();
    cin.get();
  }
}

void CGAmbulance::set_dual() {
  // for(auto a: data.A[type_call0]){
  // 	for(int b = 0; b < data.num_bases; ++b){
  // 		beta0_ab[a][b] = con_base_t0[a][b].get(GRB_DoubleAttr_Pi);
  // 		// if(abs(beta0_ab[a][b]) > g_params.EPS){
  // 		// 	std::cout << "beta0_" << a << "_" << b << ": ";
  // 		// 	std::cout << beta0_ab[a][b] << "\n";
  // 		// }
  // 	}
  // }

  // for(int a = 0; a < data.types_amb; ++a){
  // 	for(int h = 0; h < data.num_hosps; ++h){
  // 		alpha0_ah[a][h] = con_hospital_t0[a][h].get(GRB_DoubleAttr_Pi);
  // 		// if(abs(alpha0_ah[a][h] )> g_params.EPS){
  // 		// 	std::cout << "alpha0_" << a << "_" << h << ": ";
  // 		// 	std::cout << alpha0_ah[a][h] << "\n";
  // 		// }
  // 	}
  // }

  // for(int a = 0; a < data.types_amb; ++a){
  // 	for(int b = 0; b < data.num_bases; ++b){
  // 		for(int l1 = 0; l1 < data.total_locals;++l1){
  // 			if(data.L_tab[0][a][b].find(l1) !=
  // data.L_tab[0][a][b].end()){ psi0_alb[a][l1][b] =
  // con_location_t0[a][l1][b].get(GRB_DoubleAttr_Pi);
  // 				// if(abs(psi0_alb[a][l1][b]) > g_params.EPS){
  // 				// 	std::cout << "psi0_" << a << "_" << l1
  // << "_" << b << ": ";
  // 				// 	std::cout << psi0_alb[a][l1][b] << "\n";
  // 				// }
  // 			}
  // 		}
  // 	}
  // }

  // for(int c = 0; c < data.types_call; ++c){
  // 	for(int l = 0; l < data.num_locals; ++l){
  // 		phi0_cl[c][l] = con_queue_t0[c][l].get(GRB_DoubleAttr_Pi);
  // 		// if(abs(phi0_cl[c][l]) > g_params.EPS){
  // 		// 	std::cout << "phi0_" << c << "_" << l << ": ";
  // 		// 	std::cout << phi0_cl[c][l] << "\n";
  // 		// }
  // 	}
  // }
  // std::ofstream arq("dual_sols.txt");
  // arq << "Iter " << iter  << ":\n";

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int b = 0; b < data.num_bases; ++b) {
        beta_tab[t][a][b] = con_base[t][a][b].get(GRB_DoubleAttr_Pi);
        gamma_tab[t][a][b] = con_amb_base[t][a][b].get(GRB_DoubleAttr_Pi);
        // if(gamma_tab[t][a][b] > g_params.EPS){
        // 	std::cout << "gamma_" << t << "_" << a << "_" << b << ": ";
        // 	std::cout << gamma_tab[t][a][b] << "\n";
        // }
      }

      for (int h = 0; h < data.num_hosps; ++h) {
        alpha_tah[t][a][h] = con_hospital[t][a][h].get(GRB_DoubleAttr_Pi);
        // if(abs(alpha_tah[t][a][h]) > g_params.EPS){
        // 	std::cout << "alpha_" << t+1 << "_" << a << "_" << h << ",";
        // 	std::cout << alpha_tah[t][a][h] << "\n";
        // }
      }

      for (int l1 = 0; l1 < data.num_locals; ++l1) {
        for (int b = 0; b < data.num_bases; ++b) {
          if (data.L_tab[1][a][b].find(l1) != data.L_tab[1][a][b].end()) {
            psi_talb[t][a][l1][b] =
                con_location[t][a][l1][b].get(GRB_DoubleAttr_Pi);
            // if(abs(psi_talb[t][a][l1][b]) > g_params.EPS){
            // 	arq << "psi_" << t+1 << "_" << a << "_" << l1 << "_";
            // 	arq << b << " = " << psi_talb[t][a][l1][b] << "\n";
            // }
            if (l1 != data.bases[b]) {
              theta_talb[t][a][l1][b] =
                  con_amb_location[t][a][l1][b].get(GRB_DoubleAttr_Pi);
              // if(abs(theta_talb[t][a][l1][b]) > g_params.EPS){
              // 	arq << "theta_" << t+1 << "_" << a << "_" << l1 << "_";
              // 	arq << b << "," << theta_talb[t][a][l1][b] << "\n";
              // }
            }
          }
        }
      }
    }

    for (int b = 0; b < data.num_bases; ++b) {
      ups_tb[t][b] = con_cap_base[t][b].get(GRB_DoubleAttr_Pi);
      // if(abs(ups_tb[t][b]) > g_params.EPS){
      // 	std::cout << "ups_" << t+1 << "_" << b << ": " << ups_tb[t][b]
      // << "\n";
      // }
    }

    for (int c = 0; c < data.types_call; ++c) {
      for (int l = 0; l < data.num_locals; ++l) {
        phi_tcl[t][c][l] = con_queue[t][c][l].get(GRB_DoubleAttr_Pi);
        // if(abs(phi_tcl[t][c][l]) > g_params.EPS){
        // 	std::cout << "phi_" << t+1 << "_" << c << "_" << l;
        // 	std::cout << "," << phi_tcl[t][c][l] << "\n";
        // }
      }
    }
  }

  // for(int a = 0; a < data.types_amb; ++a){
  // 	for(int b = 0; b < data.num_bases; ++b){
  // 		int t = data.num_times-1;
  // 		gamma_tab[t][a][b] =
  // con_amb_base[t][a][b].get(GRB_DoubleAttr_Pi);
  // if(abs(gamma_tab[t][a][b] > g_params.EPS){
  // std::cout << "gamma_" << t << "_" << a << "_" << b << ": ";
  // std::cout << gamma_tab[t][a][b] << "\n";
  // 		}
  // 	}
  // }
}

double CGAmbulance::sub_problem(int t, int a, int l, int l1, int c, int h,
                                int b) {
  int next_l = data.L(t + 1, a, l1, data.bases[b]);
  int arr_time =
      get_time(data, (t + 1) * data.quantum + data.tao[t][c][a][l1][l][h]);
  const int factor = 1;
  double reduced_cost = f(t + 1, c, a, l1, l, h, data, factor);
  if (next_l == data.bases[b]) {
    reduced_cost += beta_tab[t][a][b];
  }

  if (arr_time - 1 < data.num_times - 1) {
    reduced_cost -= alpha_tah[arr_time - 1][a][h];
  }

  if (data.L_tab[t + 1][a][b].find(next_l) != data.L_tab[t + 1][a][b].end()) {
    // std::cout << "psi!!!!\n";
    reduced_cost += psi_talb[t][a][next_l][b];
  }

  reduced_cost += phi_tcl[t][c][l];
  reduced_cost -= theta_talb[t][a][l1][b];

  // if(reduced_cost < -g_params.EPS){
  // 	std::cout << t+1 << " " << c << " " << a << " " << l1 << " ";
  // 	std::cout << b << " " << l << " " << h << ": ";
  // 	std::cout << reduced_cost << "\n";
  // }
  // if(reduced_cost < -g_params.EPS){
  // std::cout << t+1 << " " << c << " " << a << " " << l1 << " ";
  // std::cout << b << " " << l << " " << h << ": ";
  // std::cout << "next_l | base = " << next_l << " " << data.bases[b] << "\n";
  // std::cout << "arrival_time = " << arr_time << "\n";
  // std::cout << "tao = "<< data.tao[t+1][c][a][l1][l][h] << "\n";
  // std::cout << "f = " << f(t+1,c,a,l1,l,h,data,factor) << "\n";
  // std::cout << "beta = " << beta_tab[t][a][b] << "\n";
  // if(arr_time < data.num_times)
  // std::cout << "alpha = " << alpha_tah[arr_time-1][a][h] << "\n";
  // std::cout << "psi = " << psi_talb[t][a][next_l][b] << "\n";
  // std::cout << "phi = " << phi_tcl[t][c][l] << "\n";
  // std::cout << "theta = " << theta_talb[t][a][l1][b] << "\n";
  // std::cout << "reduced_cost = " << reduced_cost << "\n";
  // std::cin.get();
  // }

  return reduced_cost;
}

bool CGAmbulance::pricing() {
  bool optimality_verified = true;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int l = 0; l < data.num_locals; ++l) {
        for (auto l1 : closest_l1_to_l[l]) {
          vh_tall[t][a][l1][l] = h_tall(t, a, l1, l);
        }
        vb_tal[t][a][l] = b_tal(t, a, l);
      }
    }

    for (int l = 0; l < data.num_locals; ++l) {
      vc_tl[t][l] = c_tl(t, l);
    }
  }

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int l = 0; l < data.num_locals; ++l) {
        bool negative_found = false;
        for (auto l1 : closest_l1_to_l[l]) {
          // int h = h_tall(t,a,l1,l);
          // int c = c_tl(t,l);
          // int b = b_tal(t,a,l1);
          int h = vh_tall[t][a][l1][l];
          int c = vc_tl[t][l];
          int b = vb_tal[t][a][l1];
          if (b != -1 && data.L_tab[t + 1][a][b].find(l1) !=
                             data.L_tab[t + 1][a][b].end()) {
            double reduced_cost = sub_problem(t, a, l, l1, c, h, b);
            if (reduced_cost < -g_params.EPS) {
              // std::cout << t+1 << " " << c << " " << a << " " << l1 << " " <<
              // b; std::cout << " " << l << " " << h << ": " << reduced_cost <<
              // "\n"; std::cin.get();
              negative_found = true;
              optimality_verified = false;
              std::vector<int> column{t, c, a, l1, b, l, h};
              add_column(column);
              ++columns_added;
            }
          }
          if (negative_found) break;
        }
      }
    }
  }
  // std::cout << "END PRICING\n";
  // std::cin.get();
  return !optimality_verified;
}

void CGAmbulance::add_column(std::vector<int>& col) {
  GRBColumn column = GRBColumn();
  std::stringstream name;

  int t = col[0], c = col[1], a = col[2], l1 = col[3], b = col[4];
  int l = col[5], h = col[6];
  int t2 =
      get_time(data, (t + 1) * data.quantum + data.tao[t + 1][c][a][l1][l][h]);

  // name << "xt" << t+1 << "_c" << c << "_a" << a << "_l'";
  // name << l1 << "_b" << b << "_l" << l << "_h" << h;

  if (data.L_tab[t + 1][a][b].find(l1) != data.L_tab[t + 1][a][b].end() &&
      data.L(t, a, l1, data.bases[b]) == data.bases[b] &&
      data.H[c][l].find(h) != data.H[c][l].end() &&
      data.A[c].find(a) != data.A[c].end()) {
    // con_base[t][a][b]
    // std::cout << "-1 " << con_base[t][a][b].get(GRB_StringAttr_ConstrName) <<
    // "\n";
    column.addTerm(-1, con_base[t][a][b]);
  }

  if (t2 < data.num_times && data.H[c][l].find(h) != data.H[c][l].end() &&
      data.A[c].find(a) != data.A[c].end() &&
      data.L_tab[t][a][b].find(l1) != data.L_tab[t][a][b].end()) {
    // con_hospital[t2][a][h]
    // std::cout << "1 " <<
    // con_hospital[t2][a][h].get(GRB_StringAttr_ConstrName) << "\n";
    column.addTerm(1, con_hospital[t2 - 1][a][h]);
  }

  for (int l2 = 0; l2 < data.total_locals; ++l2) {
    if (data.L_tab[t + 1][a][b].find(l2) != data.L_tab[t + 1][a][b].end() &&
        data.L(t, a, l1, data.bases[b]) == l2 &&
        data.A[c].find(a) != data.A[c].end() &&
        data.H[c][l].find(h) != data.H[c][l].end()) {
      // con_location[t][a][l2][b]
      // std::cout << "-1 " <<
      // con_location[t][a][l2][b].get(GRB_StringAttr_ConstrName) << "\n";
      column.addTerm(-1, con_location[t][a][l2][b]);
    }
  }

  if (data.H[c][l].find(h) != data.H[c][l].end() &&
      data.A[c].find(a) != data.A[c].end() &&
      data.L_tab[t + 1][a][b].find(l1) != data.L_tab[t + 1][a][b].end()) {
    // con_queue[t][c][l]
    // std::cout << "-1 " << con_queue[t][c][l].get(GRB_StringAttr_ConstrName)
    // << "\n";
    column.addTerm(-1, con_queue[t][c][l]);
  }

  if (data.L_tab[t + 1][a][b].find(l1) != data.L_tab[t + 1][a][b].end() &&
      data.H[c][l].find(h) != data.H[c][l].end() &&
      data.A[c].find(a) != data.A[c].end()) {
    // std::cout << "1 " <<
    // con_amb_location[t][a][l1][b].get(GRB_StringAttr_ConstrName) << "\n";
    column.addTerm(1, con_amb_location[t][a][l1][b]);
  }
  // std::cout << "COLUMN: " << t << " " << c << " " << a << " " << l1 << " " <<
  // b << " "; std::cout << l << " " << h << "\n";
  try {
    name << "xt" << t + 1 << "_c" << c << "_a" << a << "_l'" << l1 << "_b" << b
         << "_l";
    name << l << "_h" << h;
    double factor = 1;
    double f_val = f(t + 1, c, a, l1, data.bases[b], l, h, data, factor);
    xt_calblh[t][c][a][l1][b][l][h] = model.addVar(
        0, GRB_INFINITY, f_val, GRB_CONTINUOUS, column, name.str());
  } catch (GRBException& ex) {
    std::cout << ex.getErrorCode() << ": " << ex.getMessage() << "\n";
    std::cin.get();
  }
}

void CGAmbulance::add_bases_constraints_t0() {
  std::stringstream name;
  GRBLinExpr con = 0;

  // 0.5 - Flow at bases, t = 0
  for (int a = 0; a < data.types_amb; ++a) {
    for (int b = 0; b < data.num_bases; ++b) {
      con += A0_ab[a][b];
      for (int h = 0; h < data.num_hosps; ++h) {
        if (data.L(0, a, data.hospitals[h], data.bases[b]) == data.bases[b]) {
          con += y0_ahb[a][h][b];
        }
      }

      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        if (data.L(0, a, l1, data.bases[b]) == data.bases[b] &&
            (data.L_tab[0][a][b].find(l1) != data.L_tab[0][a][b].end() ||
             data.is_hospital[l1])) {
          con += A0_alb[a][l1][b];
        }
      }

      name << "flow_bases_a" << a << "_b" << b;
      model.addConstr(con, GRB_EQUAL, At_ab[0][a][b], name.str());
      name.str("");
      con = 0;
    }
  }
}

void CGAmbulance::add_hospitals_constraints_t0() {
  std::stringstream name;
  GRBLinExpr con = 0;

  // 0.6 - Flow at hospitals
  for (int a = 0; a < data.types_amb; ++a) {
    for (int h = 0; h < data.num_hosps; ++h) {
      for (int b = 0; b < data.num_bases; ++b) {
        con += y0_ahb[a][h][b];
      }

      for (int c = 0; c < data.types_call; ++c) {
        for (int l = 0; l < data.num_locals; ++l) {
          for (int h1 = 0; h1 < data.num_hosps; ++h1) {
            if (data.H[c][l].find(h1) != data.H[c][l].end()) {
              if (data.A[c].find(a) != data.A[c].end())
                con += y0_ahclh[a][h][c][l][h1];
            }
          }
        }
      }

      name << "flow_hospitals_a" << a << "_h" << h;
      model.addConstr(con, GRB_EQUAL, A0_ah[a][h], name.str());
      name.str("");
      con = 0;
    }
  }
}

void CGAmbulance::add_locations_constraints_t0() {
  std::stringstream name;
  GRBLinExpr con = 0;
  // 0.7 - Fluxo at location between hospitals and bases
  for (int a = 0; a < data.types_amb; ++a) {
    for (int b = 0; b < data.num_bases; ++b) {
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        if (data.L_tab[0][a][b].find(l1) != data.L_tab[0][a][b].end()) {
          for (int h = 0; h < data.num_hosps; ++h) {
            if (data.L(0, a, data.hospitals[h], data.bases[b]) == l1) {
              con += y0_ahb[a][h][b];
            }
          }

          for (int l2 = 0; l2 < data.total_locals; ++l2) {
            if (data.L(0, a, l2, data.bases[b]) == l1 &&
                data.L_tab[0][a][b].find(l2) != data.L_tab[0][a][b].end()) {
              con += A0_alb[a][l2][b];
            }
          }
          name << "flow_locals_a" << a << "_b" << b << "_l" << l1;
          model.addConstr(con, GRB_EQUAL, At_alb[0][a][l1][b], name.str());
          name.str("");
          con = 0;
        }
      }
    }
  }
}

void CGAmbulance::add_queues_constraints_t0() {
  std::stringstream name;
  GRBLinExpr con = 0;

  // 0.8 - Flow at queues
  for (int c = 0; c < data.types_call; ++c) {
    for (int l = 0; l < data.num_locals; ++l) {
      con += C[c][l];

      for (int h = 0; h < data.num_hosps; ++h) {
        if (data.H[c][l].find(h) != data.H[c][l].end()) {
          for (int a = 0; a < data.types_amb; ++a) {
            for (int h1 = 0; h1 < data.num_hosps; ++h1) {
              if (data.A[c].find(a) != data.A[c].end()) {
                con -= y0_ahclh[a][h1][c][l][h];
              }
            }
          }
        }
      }

      name << "flow_queues_c" << c << "_l" << l;
      model.addConstr(con, GRB_EQUAL, Ct_cl[0][c][l], name.str());
      name.str("");
      con = 0;
    }
  }
}

void CGAmbulance::add_bases_constraints() {
  std::stringstream name;
  GRBLinExpr con = 0;

  // 0.9 - Flow at bases,  t > 0
  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int b = 0; b < data.num_bases; ++b) {
        con += At_ab[t][a][b];
        for (int h = 0; h < data.num_hosps; ++h) {
          if (data.L(t, a, data.hospitals[h], data.bases[b]) == data.bases[b]) {
            con += yt_ahb[t][a][h][b];
          }
        }
        for (int l1 = 0; l1 < data.total_locals; ++l1) {
          if (data.L_tab[t + 1][a][b].find(l1) !=
                  data.L_tab[t + 1][a][b].end() &&
              data.L(t, a, l1, data.bases[b]) == data.bases[b]) {
            con += At_alb[t][a][l1][b];
            // for(int c = 0; c < data.types_call; ++c){
            // 	for(int l = 0; l < data.num_locals; ++l){
            // 		for(int h = 0; h < data.num_hosps; ++h){
            // 			if(data.H[c][l].find(h) !=
            // 				data.H[c][l].end() &&
            // 				data.A[c].find(a) != data.A[c].end()){
            // 				//If using xt0, put t+1 on 1st index
            // 				con -= xt_calblh[t][c][a][l1][b][l][h];
            // 			}
            // 		}
            // 	}
            // }
          }
        }

        for (int c = 0; c < data.types_call; ++c) {
          for (int l = 0; l < data.num_locals; ++l) {
            for (int h = 0; h < data.num_hosps; ++h) {
              if (data.H[c][l].find(h) != data.H[c][l].end() &&
                  data.A[c].find(a) != data.A[c].end()) {
                // If using xt0, put t+1 on 1st index
                con -= xt_cablh[t][c][a][b][l][h];
              }
            }
          }
        }

        name << "flow_bases_t" << t + 1 << "_a" << a << "_b" << b;
        try {
          con_base[t][a][b] =
              model.addConstr(con, GRB_EQUAL, At_ab[t + 1][a][b], name.str());
        } catch (GRBException ex) {
          std::cout << ex.getErrorCode() << ": ";
          std::cout << ex.getMessage() << "\n";
          exit(1);
        }
        name.str("");
        con = 0;
      }
    }
  }
}

void CGAmbulance::add_hospitals_constraints() {
  GRBLinExpr con;
  std::stringstream name;
  int quantum = data.quantum;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int h = 0; h < data.num_hosps; ++h) {
        for (int c = 0; c < data.types_call; ++c) {
          for (int l = 0; l < data.num_locals; ++l) {
            for (int h1 = 0; h1 < data.num_hosps; ++h1) {
              if (data.H[c][l].find(h1) != data.H[c][l].end() &&
                  data.A[c].find(a) != data.A[c].end()) {
                // If using xt0, put t+1 on 1st index
                con -= xt_cahlh[t][c][a][h][l][h1];
              }
            }
          }
        }

        for (int b = 0; b < data.num_bases; ++b) {
          con -= yt_ahb[t][a][h][b];
        }

        for (int c = 0; c < data.types_call; ++c) {
          for (int l1 = 0; l1 < data.total_locals; ++l1) {
            for (int l = 0; l < data.num_locals; ++l) {
              if ((t)*quantum <= data.tao[0][c][a][l1][l][h] &&
                  data.tao[0][c][a][l1][l][h] <= (t + 1) * quantum) {
                con += A0_callh[c][a][l1][l][h];
              }
            }
          }
        }

        for (int c = 0; c < data.types_call; ++c) {
          for (int l1 = 0; l1 < data.total_locals; ++l1) {
            if ((t)*quantum <= data.tao_0[c][a][l1][h] &&
                data.tao_0[c][a][l1][h] <= (t + 1) * quantum) {
              con += A0_calh[c][a][l1][h];
            }
          }
        }

        for (int c = 0; c < data.types_call; ++c) {
          for (int l = 0; l < data.num_locals; ++l) {
            if (data.H[c][l].find(h) != data.H[c][l].end()) {
              for (int b = 0; b < data.num_bases; ++b) {
                for (int t1 = 1; t1 < t + 1; ++t1) {
                  if ((t)*quantum <=
                          t1 * quantum +
                              data.tao[t1][c][a][data.bases[b]][l][h] &&
                      t1 * quantum + data.tao[t1][c][a][data.bases[b]][l][h] <=
                          (t + 1) * quantum &&
                      data.A[c].find(a) != data.A[c].end()) {
                    con += xt_cablh[t1 - 1][c][a][b][l][h];
                  }
                }
              }
            }
          }
        }

        for (int c = 0; c < data.types_call; ++c) {
          for (int l = 0; l < data.num_locals; ++l) {
            for (int h1 = 0; h1 < data.num_hosps; ++h1) {
              if ((t)*data.quantum <=
                      data.tao[0][c][a][data.hospitals[h1]][l][h] &&
                  data.tao[0][c][a][data.hospitals[h1]][l][h] <=
                      (t + 1) * data.quantum &&
                  data.H[c][l].find(h) != data.H[c][l].end() &&
                  data.A[c].find(a) != data.A[c].end()) {
                con += y0_ahclh[a][h1][c][l][h];
              }
            }
          }
        }

        for (int c = 0; c < data.types_call; ++c) {
          for (int l = 0; l < data.num_locals; ++l) {
            if (data.H[c][l].find(h) != data.H[c][l].end()) {
              for (int h1 = 0; h1 < data.num_hosps; ++h1) {
                for (int t1 = 1; t1 < t + 1; ++t1) {
                  if ((t)*quantum <=
                          t1 * quantum +
                              data.tao[t1][c][a][data.hospitals[h1]][l][h] &&
                      t1 * quantum +
                              data.tao[t1][c][a][data.hospitals[h1]][l][h] <=
                          (t + 1) * quantum &&
                      data.A[c].find(a) != data.A[c].end()) {
                    con += xt_cahlh[t1 - 1][c][a][h1][l][h];
                  }
                }
              }
            }
          }
        }

        // for(int c = 0; c < data.types_call; ++c){
        // 	for(int l = 0; l < data.num_locals; ++l){
        // 		if(data.H[c][l].find(h) != data.H[c][l].end()){
        // 			for(int l1 = 0; l1 < data.total_locals; ++l1){
        // 				for(int b = 0; b < data.num_bases; ++b){
        // 					for(int t1 = 1; t1 < t+1; ++t1){
        // 						if((t)*quantum <=
        // t1*quantum+
        // data.tao[t1][c][a][l1][l][h] &&
        // 							t1*quantum+data.tao[t1][c][a][l1][l][h]
        // <= (t+1)*quantum &&
        // data.A[c].find(a) != data.A[c].end()
        // 							&&
        // data.L_tab[t1-1][a][b].find(l1) !=
        // data.L_tab[t1-1][a][b].end()){
        // con
        // += xt_calblh[t1-1][c][a][l1][b][l][h];
        // 						}
        // 					}
        // 				}
        // 			}
        // 		}
        // 	}
        // }

        name << "flow_hospitals_t" << t + 1 << "_a" << a << "_h" << h;
        con_hospital[t][a][h] = model.addConstr(con, GRB_EQUAL, 0, name.str());
        name.str("");
        con = 0;
      }
    }
  }
}

void CGAmbulance::add_locations_constraints() {
  GRBLinExpr con;
  std::stringstream name;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int b = 0; b < data.num_bases; ++b) {
        for (int l1 = 0; l1 < data.total_locals; ++l1) {
          if (data.L_tab[t + 1][a][b].find(l1) !=
              data.L_tab[t + 1][a][b].end()) {
            for (int h1 = 0; h1 < data.num_hosps; ++h1) {
              if (data.L(t, a, data.hospitals[h1], data.bases[b]) == l1) {
                con += yt_ahb[t][a][h1][b];
              }
            }

            for (int l2 = 0; l2 < data.total_locals; ++l2) {
              if (data.L_tab[t + 1][a][b].find(l2) !=
                      data.L_tab[t + 1][a][b].end() &&
                  data.L(t, a, l2, data.bases[b]) == l1) {
                con += At_alb[t][a][l2][b];
                // for(int c = 0; c < data.types_call; ++c){
                // 	for(int l = 0; l < data.num_locals; ++l){
                // 		for(int h = 0; h < data.num_hosps; ++h){
                // 			if(data.H[c][l].find(h) !=
                // 				data.H[c][l].end() &&
                // 				data.A[c].find(a)!=data.A[c].end()){
                // 				//If using xt0, put t+1 on 1st
                // index 				con -=
                // xt_calblh[t][c][a][l2][b][l][h];
                // 			}
                // 		}
                // 	}
                // }
              }
            }

            name << "flow_locals_t" << t + 1 << "_" << a << "_" << b << "_"
                 << l1;
            con_location[t][a][l1][b] = model.addConstr(
                con, GRB_EQUAL, At_alb[t + 1][a][l1][b], name.str());
            name.str("");
            con = 0;
          }
        }
      }
    }
  }
}

void CGAmbulance::add_queues_constraints() {
  GRBLinExpr con;
  std::stringstream name;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int c = 0; c < data.types_call; ++c) {
      for (int l = 0; l < data.num_locals; ++l) {
        con += Ct_cl[t][c][l];
        con += data.lambda[t][c][l];
        for (int a = 0; a < data.types_amb; ++a) {
          if (data.A[c].find(a) != data.A[c].end()) {
            for (int b = 0; b < data.num_bases; ++b) {
              for (int h = 0; h < data.num_hosps; ++h) {
                if (data.H[c][l].find(h) != data.H[c][l].end() &&
                    data.A[c].find(a) != data.A[c].end()) {
                  // If using xt0, put t+1 on 1st index
                  con -= xt_cablh[t][c][a][b][l][h];
                }
              }
            }
          }
        }

        for (int a = 0; a < data.types_amb; ++a) {
          if (data.A[c].find(a) != data.A[c].end()) {
            for (int h1 = 0; h1 < data.num_hosps; ++h1) {
              for (int h = 0; h < data.num_hosps; ++h) {
                if (data.H[c][l].find(h) != data.H[c][l].end() &&
                    data.A[c].find(a) != data.A[c].end()) {
                  // If using xt0, put t+1 on 1st index
                  con -= xt_cahlh[t][c][a][h1][l][h];
                }
              }
            }
          }
        }

        // for(int a = 0; a < data.types_amb; ++a){
        // 	if(data.A[c].find(a) != data.A[c].end()){
        // 		for(int b = 0; b < data.num_bases; ++b){
        // 			for(int l1 = 0; l1 < data.total_locals; ++l1){
        // 				for(int h = 0; h < data.num_hosps; ++h){
        // 					if(data.H[c][l].find(h) !=
        // data.H[c][l].end() &&
        // data.A[c].find(a) != data.A[c].end() &&
        // 						data.L_tab[t+1][a][b].find(l1)
        // != data.L_tab[t+1][a][b].end()){
        // 						//If using xt0, put t+1
        // on 1st index 						con -=
        // xt_calblh[t][c][a][l1][b][l][h];
        // 					}
        // 				}
        // 			}
        // 		}
        // 	}
        // }
        name << "flow_queue_t" << t + 1 << "_c" << c << "_l" << l;
        con_queue[t][c][l] =
            model.addConstr(con, GRB_EQUAL, Ct_cl[t + 1][c][l], name.str());
        name.str("");
        con = 0;
      }
    }
  }
}

void CGAmbulance::add_ambs_location_constraints() {
  GRBLinExpr con = 0;
  std::stringstream name;

  // Location constraints
  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int b = 0; b < data.num_bases; ++b) {
        for (int l1 = 0; l1 < data.total_locals; ++l1) {
          if (data.L_tab[t + 1][a][b].find(l1) !=
              data.L_tab[t + 1][a][b].end()) {
            // for(int c = 0; c < data.types_call; ++c){
            // 	for(int l = 0; l < data.num_locals; ++l){
            // 		for(int h = 0; h < data.num_hosps; ++h){
            // 			if(data.H[c][l].find(h) != data.H[c][l].end() &&
            // 				data.A[c].find(a) != data.A[c].end() &&
            // 				data.L_tab[t][a][b].find(l1) !=
            // 				data.L_tab[t][a][b].end()){
            // 				//If using xt0, put t+1 on 1st index
            // 				con += xt_calblh[t][c][a][l1][b][l][h];
            // 			}
            // 		}
            // 	}
            // }
            name << "dispatch_locals_t" << t + 1 << "_a" << a << "_b";
            name << b << "_l" << l1;
            con_amb_location[t][a][l1][b] = model.addConstr(
                con, GRB_LESS_EQUAL, At_alb[t][a][l1][b], name.str());
            name.str("");
            con = 0;
          }
        }
      }
    }
  }

  // Max of ambulance at bases
  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int b = 0; b < data.num_bases; ++b) {
        for (int c = 0; c < data.types_call; ++c) {
          for (int l = 0; l < data.num_locals; ++l) {
            for (int h = 0; h < data.num_hosps; ++h) {
              if (data.H[c][l].find(h) != data.H[c][l].end() &&
                  data.A[c].find(a) != data.A[c].end()) {
                // If using xt0, put t+1 on 1st index
                con += xt_cablh[t][c][a][b][l][h];
              }
            }
          }
        }

        name << "dispatch_bases_t" << t + 1 << "_a" << a << "_b" << b;
        con_amb_base[t][a][b] =
            model.addConstr(con, GRB_LESS_EQUAL, At_ab[t][a][b], name.str());
        name.str("");
        con = 0;
      }
    }
  }
}

void CGAmbulance::add_base_cap_constraints() {
  GRBLinExpr con;
  std::stringstream name;

  for (int t = 0; t < data.num_times; ++t) {
    for (int b = 0; b < data.num_bases; ++b) {
      for (int a = 0; a < data.types_amb; ++a) {
        con += At_ab[t][a][b];
      }
      name << "cap_bases_" << t + 1 << "_" << b;
      con_cap_base[t][b] =
          model.addConstr(con, GRB_LESS_EQUAL, data.cap_bases[b], name.str());
      name.str("");
      con = 0;
    }
  }
}

// int CGAmbulance::get_location(Ambulance& amb){
// 	// Ambulance stopped
// 	if(amb.waiting_return || (!amb.busy && amb.base_ret == -1 &&
// !amb.waiting_return)){ 		return amb.location;
// 	}

// 	// Ambulance already at location
// 	if(amb.busy && amb.location == data.hospitals[amb.call_h_dest] ||
// 		!amb.busy && amb.base_ret != -1 && amb.location ==
// data.bases[amb.base_ret]){ 		return amb.location;
// 	}

// 	int departure_time = sim.get_time(amb.departure_time);
// 	int elapsed_time = (t0-departure_time)*data.quantum;
// 	int next = -1;
// 	std::vector<int> path;
// 	if(amb.busy){
// 		path = min_path(data, amb.origin, amb.call_location);
// 		std::vector<int> aux = min_path(data, amb.call_location,
// 			data.hospitals[amb.call_h_dest]);
// 		for(int i = 1; i < aux.size(); ++i){
// 			path.push_back(aux[i]);
// 		}
// 	}else{
// 		path = min_path(data, amb.origin, data.bases[amb.base_ret]);
// 	}

// 	if(path.size() == 1){
// 		next = path[0];
// 	}else if(path.size() == 2){
// 		if(data.times[path[0]][path[1]] < elapsed_time){
// 			next = path[1];
// 		}else{
// 			next = path[0];
// 		}
// 	}else{
// 		for(int i = 0; i < path.size()-1; ++i){
// 			if(data.times[path[i]][path[i+1]] < elapsed_time){
// 				next = path[i+1];
// 				elapsed_time -= data.times[path[i]][path[i+1]];
// 			}else{
// 				next = path[i];
// 				break;
// 			}
// 		}
// 	}
// 	return amb.location;
// }

// int CGAmbulance::is_at_hospital(Ambulance & amb){
// 	for(int h = 0; h < data.num_hosps; ++h){
// 		if(amb.location == data.hospitals[h] && amb.waiting_return){
// 			return h;
// 		}
// 	}
// 	return -1;
// }

// int CGAmbulance::is_at_base(Ambulance & amb){
// 	for(int b = 0; b < data.num_bases; ++b){
// 		if(amb.location == data.bases[b] && amb.base_ret == -1 &&
// !amb.waiting_return){ 			return b;
// 		}
// 	}
// 	return -1;
// }

// std::pair<int,int> CGAmbulance::is_at_location_base(Ambulance & amb){
// 	if(amb.base_ret < 0)
// 		return std::make_pair(-1,-1);
// 	else
// 		return std::make_pair(amb.location,amb.base_ret);
// }

int CGAmbulance::c_tl(int t, int l) {
  int arg_min = -1;
  double min_val = GRB_INFINITY;
  // std::cout << "phi_tcl " << t+1 << " " << l << ":\n";
  for (int c = 0; c < data.types_call; ++c) {
    // std::cout << c << ": phi = " << phi_tcl[t][c][l] << "\n";
    if (phi_tcl[t][c][l] < min_val) {
      min_val = phi_tcl[t][c][l];
      arg_min = c;
    }
  }
  // std::cout << "min = " << arg_min << "\n";
  // std::cin.get();
  return arg_min;
}

int CGAmbulance::b_tal(int t, int a, int l) {
  int arg_min = -1;
  double min_val = GRB_INFINITY;
  // std::cout << "b_tal " << t+1 << " " << a << " " << l << ":\n";
  for (int b = 0; b < data.num_bases; ++b) {
    if (data.L_tab[t + 1][a][b].find(l) != data.L_tab[t + 1][a][b].end()) {
      double sum = 0;
      int next_l = data.L(t + 1, a, l, data.bases[b]);
      // std::cout << b << ": ";
      if (next_l == data.bases[b]) {
        sum += beta_tab[t][a][b];
        // std::cout << "beta = " << beta_tab[t][a][b] << " ";
      }

      if (data.L_tab[t + 1][a][b].find(next_l) !=
          data.L_tab[t + 1][a][b].end()) {
        sum += psi_talb[t][a][next_l][b];
        // std::cout << ", psi = " << psi_talb[t][a][next_l][b] << " ";
      }

      sum += theta_talb[t][a][l][b];
      // std::cout << ", theta = " << theta_talb[t][a][l][b] << ", ";
      // std::cout << "sum = " << sum << "\n";
      if (sum < min_val) {
        min_val = sum;
        arg_min = b;
      }
    }
  }
  // std::cout << "min = " << arg_min << "\n";
  // std::cin.get();
  return arg_min;
}

int CGAmbulance::h_tall(int t, int a, int l1, int l) {
  int arg_min = -1;
  double min_val = GRB_INFINITY;
  int factor = 1;
  int current_time = (t + 1) * data.quantum;
  // std::cout << "h_tall " << t+1 << " "  << a << " " << l1 << " " << l <<
  // ":\n";
  for (auto h : data.H[0][l]) {
    int f_val = f(t + 1, 0, a, l1, l, h, data, factor);
    // std::cout << h << ": f = " << f_val << ", arr_time = ";
    int arrival_time =
        get_time(data, current_time + data.tao[t][0][a][l1][l][h]);
    // std::cout << arrival_time;
    double sum = f_val;
    if (arrival_time - 1 < data.num_times - 1) {
      // std::cout << ", alpha = " << alpha_tah[arrival_time-1][a][h] << "\n";
      sum -= alpha_tah[arrival_time - 1][a][h];
    } else {
      // std::cout << "\n";
    }

    if (sum < min_val) {
      min_val = sum;
      arg_min = h;
    }
  }
  // std::cout << "min = " << arg_min << "\n";
  // std::cin.get();
  return arg_min;
}

CGAmbulance::~CGAmbulance() {
  for (int a = 0; a < data.types_amb; ++a) {
    for (int l1 = 0; l1 < data.total_locals; ++l1) {
      delete[] A0_alb[a][l1];
    }
    delete[] A0_ab[a];
    delete[] A0_ah[a];
    delete[] A0_alb[a];
  }

  delete[] A0_ab;
  delete[] A0_ah;
  delete[] A0_alb;

  for (int c = 0; c < data.types_call; ++c) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        delete[] A0_calh[c][a][l1];
      }
      delete[] A0_calh[c][a];
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        for (int l = 0; l < data.num_locals; ++l) {
          delete[] A0_callh[c][a][l1][l];
        }
        delete[] A0_callh[c][a][l1];
      }
      delete[] A0_callh[c][a];
    }
    delete[] A0_calh[c];
    delete[] A0_callh[c];
    delete[] C[c];
  }
  delete[] A0_callh;
  delete[] A0_calh;
  delete[] C;

  for (int a = 0; a < data.types_amb; ++a) {
    for (int h = 0; h < data.num_hosps; ++h) {
      delete[] y0_ahb[a][h];
      for (int c = 0; c < data.types_call; ++c) {
        for (int l = 0; l < data.num_locals; ++l) {
          delete[] y0_ahclh[a][h][c][l];
        }
        delete[] y0_ahclh[a][h][c];
      }
      delete[] y0_ahclh[a][h];
    }
    delete[] y0_ahb[a];
    delete[] y0_ahclh[a];
  }

  delete[] y0_ahb;
  delete[] y0_ahclh;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int c = 0; c < data.types_call; ++c) {
      for (int a = 0; a < data.types_amb; ++a) {
        for (int b = 0; b < data.num_bases; ++b) {
          for (int l = 0; l < data.num_locals; ++l) {
            delete[] xt_cablh[t][c][a][b][l];
          }
          delete[] xt_cablh[t][c][a][b];
        }
        delete[] xt_cablh[t][c][a];

        for (int h1 = 0; h1 < data.num_hosps; ++h1) {
          for (int l = 0; l < data.num_locals; ++l) {
            delete[] xt_cahlh[t][c][a][h1][l];
          }
          delete[] xt_cahlh[t][c][a][h1];
        }
        delete[] xt_cahlh[t][c][a];

        for (int l1 = 0; l1 < data.total_locals; ++l1) {
          for (int b = 0; b < data.num_bases; ++b) {
            for (int l = 0; l < data.num_locals; ++l) {
              delete[] xt_calblh[t][c][a][l1][b][l];
            }
            delete[] xt_calblh[t][c][a][l1][b];
          }
          delete[] xt_calblh[t][c][a][l1];
        }
        delete[] xt_calblh[t][c][a];
      }
      delete[] xt_cablh[t][c];
      delete[] xt_cahlh[t][c];
      delete[] xt_calblh[t][c];
    }
    delete[] xt_cablh[t];
    delete[] xt_cahlh[t];
    delete[] xt_calblh[t];

    for (int a = 0; a < data.types_amb; ++a) {
      for (int h = 0; h < data.num_hosps; ++h) {
        delete[] yt_ahb[t][a][h];
      }
      delete[] yt_ahb[t][a];
    }
    delete[] yt_ahb[t];
  }
  delete[] xt_cablh;
  delete[] xt_cahlh;
  delete[] xt_calblh;
  delete[] yt_ahb;
  for (int t = 0; t < data.num_times; ++t) {
    for (int c = 0; c < data.types_call; ++c) {
      delete[] Ct_cl[t][c];
    }
    delete[] Ct_cl[t];
    for (int a = 0; a < data.types_amb; ++a) {
      delete[] At_ab[t][a];
    }
    delete[] At_ab[t];

    for (int a = 0; a < data.types_amb; ++a) {
      for (int l = 0; l < data.total_locals; ++l) {
        delete[] At_alb[t][a][l];
      }
      delete[] At_alb[t][a];
    }
    delete[] At_alb[t];
  }
  delete[] Ct_cl;
  delete[] At_ab;
  delete[] At_alb;

  for (int a = 0; a < data.types_amb; ++a) {
    delete[] con_base_t0[a];
  }
  delete[] con_base_t0;

  for (int a = 0; a < data.types_amb; ++a) {
    delete[] con_hospital_t0[a];
  }
  delete[] con_hospital_t0;

  for (int a = 0; a < data.types_amb; ++a) {
    for (int l1 = 0; l1 < data.total_locals; ++l1) {
      delete[] con_location_t0[a][l1];
    }
    delete[] con_location_t0[a];
  }
  delete[] con_location_t0;

  for (int c = 0; c < data.types_call; ++c) {
    delete[] con_queue_t0[c];
  }
  delete[] con_queue_t0;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      delete[] con_base[t][a];
    }
    delete[] con_base[t];
  }
  delete[] con_base;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      delete[] con_hospital[t][a];
    }
    delete[] con_hospital[t];
  }
  delete[] con_hospital;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        delete[] con_location[t][a][l1];
      }
      delete[] con_location[t][a];
    }
    delete[] con_location[t];
  }
  delete[] con_location;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int c = 0; c < data.types_call; ++c) {
      delete[] con_queue[t][c];
    }
    delete[] con_queue[t];
  }
  delete[] con_queue;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        delete[] con_amb_location[t][a][l1];
      }
      delete[] con_amb_location[t][a];
    }
    delete[] con_amb_location[t];
  }
  delete[] con_amb_location;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      delete[] con_amb_base[t][a];
    }
    delete[] con_amb_base[t];
  }
  delete[] con_amb_base;

  for (int t = 0; t < data.num_times; ++t) {
    delete[] con_cap_base[t];
  }
  delete[] con_cap_base;

  for (int a = 0; a < data.types_amb; ++a) {
    delete[] beta0_ab[a];
    delete[] alpha0_ah[a];
    for (int l1 = 0; l1 < data.total_locals; ++l1) {
      delete[] psi0_alb[a][l1];
    }
    delete[] psi0_alb[a];
  }
  delete[] beta0_ab;
  delete[] alpha0_ah;
  delete[] psi0_alb;

  for (int c = 0; c < data.types_call; ++c) {
    delete[] phi0_cl[c];
  }
  delete[] phi0_cl;

  for (int t = 0; t < data.num_times - 1; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int l1 = 0; l1 < data.total_locals; ++l1) {
        delete[] psi_talb[t][a][l1];
        delete[] theta_talb[t][a][l1];
      }
      delete[] alpha_tah[t][a];
      delete[] psi_talb[t][a];
      delete[] theta_talb[t][a];
      delete[] beta_tab[t][a];
      delete[] gamma_tab[t][a];
    }

    for (int c = 0; c < data.types_call; ++c) {
      delete[] phi_tcl[t][c];
    }
    delete[] beta_tab[t];
    delete[] alpha_tah[t];
    delete[] psi_talb[t];
    delete[] phi_tcl[t];
    delete[] ups_tb[t];
    delete[] theta_talb[t];
    delete[] gamma_tab[t];
  }
  delete[] beta_tab;
  delete[] alpha_tah;
  delete[] psi_talb;
  delete[] phi_tcl;
  delete[] ups_tb;
  delete[] theta_talb;
  delete[] gamma_tab;

  for (int t = 0; t < data.num_times; ++t) {
    for (int a = 0; a < data.types_amb; ++a) {
      for (int l1 = 0; l1 < data.num_locals; ++l1) {
        delete[] vh_tall[t][a][l1];
      }
      delete[] vh_tall[t][a];
      delete[] vb_tal[t][a];
    }

    delete[] vh_tall[t];
    delete[] vb_tal[t];
    delete[] vc_tl[t];
  }
  delete[] vh_tall;
  delete[] vc_tl;
  delete[] vb_tal;
}

void CGCall::print_solution() {
  int num_vars = model.get(GRB_IntAttr_NumVars);
  GRBVar* vars = model.getVars();

  for (int i = 0; i < num_vars; ++i) {
    string name = vars[i].get(GRB_StringAttr_VarName);
    double val = vars[i].get(GRB_DoubleAttr_X);
    if (val > 0.0001 && name.rfind("At", 0) != 0) {
      fmt::print("{} = {}\n", vars[i].get(GRB_StringAttr_VarName), val);
    }
  }
  fmt::print("time = {}, quantum = {}\n", data.time, data.quantum);
}

void CGAmbulance::print_solution() {
  int num_vars = model.get(GRB_IntAttr_NumVars);
  GRBVar* vars = model.getVars();

  for (int i = 0; i < num_vars; ++i) {
    string name = vars[i].get(GRB_StringAttr_VarName);
    double val = vars[i].get(GRB_DoubleAttr_X);
    if (val > 0.0001 && name.rfind("At", 0) != 0) {
      fmt::print("{} = {}\n", vars[i].get(GRB_StringAttr_VarName), val);
    }
  }
  fmt::print("time = {}, quantum = {}\n", data.time, data.quantum);
  delete[] vars;
}