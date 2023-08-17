module IdentificationMap barrelDIRC_epid {
  set InputArray TrackMerger/tracks
  set OutputArray tracks

    add EfficiencyFormula {-11} {-11} {
      (eta< -1.64 || eta>=  1.25 || pt <    0.44 || pt >=    3.00) * ( 0.00 ) +
      ( -1.64 <= eta && eta <  -1.35) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      ( -1.64 <= eta && eta <  -1.35) * (   0.70 <= pt && pt <    0.95) * (0.819313) +
      ( -1.64 <= eta && eta <  -1.35) * (   0.95 <= pt && pt <    1.21) * (0.721698) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.21 <= pt && pt <    1.46) * (0.658191) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.46 <= pt && pt <    1.72) * (0.618107) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.72 <= pt && pt <    1.98) * (0.590891) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.98 <= pt && pt <    2.23) * (0.572019) +
      ( -1.64 <= eta && eta <  -1.35) * (   2.23 <= pt && pt <    2.49) * (0.558412) +
      ( -1.64 <= eta && eta <  -1.35) * (   2.49 <= pt && pt <    2.74) * (0.547921) +
      ( -1.64 <= eta && eta <  -1.35) * (   2.74 <= pt && pt <    3.00) * (0.540264) +
      ( -1.35 <= eta && eta <  -1.06) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      ( -1.35 <= eta && eta <  -1.06) * (   0.70 <= pt && pt <    0.95) * (0.904531) +
      ( -1.35 <= eta && eta <  -1.06) * (   0.95 <= pt && pt <    1.21) * (0.811067) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.21 <= pt && pt <    1.46) * (0.737047) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.46 <= pt && pt <    1.72) * (0.681233) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.72 <= pt && pt <    1.98) * (0.643418) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.98 <= pt && pt <    2.23) * (0.614759) +
      ( -1.35 <= eta && eta <  -1.06) * (   2.23 <= pt && pt <    2.49) * (0.593967) +
      ( -1.35 <= eta && eta <  -1.06) * (   2.49 <= pt && pt <    2.74) * (0.578366) +
      ( -1.35 <= eta && eta <  -1.06) * (   2.74 <= pt && pt <    3.00) * (0.566124) +
      ( -1.06 <= eta && eta <  -0.77) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      ( -1.06 <= eta && eta <  -0.77) * (   0.70 <= pt && pt <    0.95) * (0.964452) +
      ( -1.06 <= eta && eta <  -0.77) * (   0.95 <= pt && pt <    1.21) * (0.889668) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.21 <= pt && pt <    1.46) * (0.814844) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.46 <= pt && pt <    1.72) * (0.752748) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.72 <= pt && pt <    1.98) * (0.703903) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.98 <= pt && pt <    2.23) * (0.666318) +
      ( -1.06 <= eta && eta <  -0.77) * (   2.23 <= pt && pt <    2.49) * (0.637320) +
      ( -1.06 <= eta && eta <  -0.77) * (   2.49 <= pt && pt <    2.74) * (0.614699) +
      ( -1.06 <= eta && eta <  -0.77) * (   2.74 <= pt && pt <    3.00) * (0.597641) +
      ( -0.77 <= eta && eta <  -0.48) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      ( -0.77 <= eta && eta <  -0.48) * (   0.70 <= pt && pt <    0.95) * (0.991715) +
      ( -0.77 <= eta && eta <  -0.48) * (   0.95 <= pt && pt <    1.21) * (0.949781) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.21 <= pt && pt <    1.46) * (0.887818) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.46 <= pt && pt <    1.72) * (0.827077) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.72 <= pt && pt <    1.98) * (0.774428) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.98 <= pt && pt <    2.23) * (0.730866) +
      ( -0.77 <= eta && eta <  -0.48) * (   2.23 <= pt && pt <    2.49) * (0.695460) +
      ( -0.77 <= eta && eta <  -0.48) * (   2.49 <= pt && pt <    2.74) * (0.666119) +
      ( -0.77 <= eta && eta <  -0.48) * (   2.74 <= pt && pt <    3.00) * (0.642302) +
      ( -0.48 <= eta && eta <  -0.20) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      ( -0.48 <= eta && eta <  -0.20) * (   0.70 <= pt && pt <    0.95) * (0.998349) +
      ( -0.48 <= eta && eta <  -0.20) * (   0.95 <= pt && pt <    1.21) * (0.979229) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.21 <= pt && pt <    1.46) * (0.936414) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.46 <= pt && pt <    1.72) * (0.880768) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.72 <= pt && pt <    1.98) * (0.828355) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.98 <= pt && pt <    2.23) * (0.780618) +
      ( -0.48 <= eta && eta <  -0.20) * (   2.23 <= pt && pt <    2.49) * (0.739694) +
      ( -0.48 <= eta && eta <  -0.20) * (   2.49 <= pt && pt <    2.74) * (0.706492) +
      ( -0.48 <= eta && eta <  -0.20) * (   2.74 <= pt && pt <    3.00) * (0.678903) +
      ( -0.20 <= eta && eta <   0.09) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      ( -0.20 <= eta && eta <   0.09) * (   0.70 <= pt && pt <    0.95) * (0.999752) +
      ( -0.20 <= eta && eta <   0.09) * (   0.95 <= pt && pt <    1.21) * (0.992859) +
      ( -0.20 <= eta && eta <   0.09) * (   1.21 <= pt && pt <    1.46) * (0.965321) +
      ( -0.20 <= eta && eta <   0.09) * (   1.46 <= pt && pt <    1.72) * (0.920646) +
      ( -0.20 <= eta && eta <   0.09) * (   1.72 <= pt && pt <    1.98) * (0.871939) +
      ( -0.20 <= eta && eta <   0.09) * (   1.98 <= pt && pt <    2.23) * (0.822796) +
      ( -0.20 <= eta && eta <   0.09) * (   2.23 <= pt && pt <    2.49) * (0.779514) +
      ( -0.20 <= eta && eta <   0.09) * (   2.49 <= pt && pt <    2.74) * (0.743131) +
      ( -0.20 <= eta && eta <   0.09) * (   2.74 <= pt && pt <    3.00) * (0.711441) +
      (  0.09 <= eta && eta <   0.38) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      (  0.09 <= eta && eta <   0.38) * (   0.70 <= pt && pt <    0.95) * (0.999599) +
      (  0.09 <= eta && eta <   0.38) * (   0.95 <= pt && pt <    1.21) * (0.994663) +
      (  0.09 <= eta && eta <   0.38) * (   1.21 <= pt && pt <    1.46) * (0.971414) +
      (  0.09 <= eta && eta <   0.38) * (   1.46 <= pt && pt <    1.72) * (0.930049) +
      (  0.09 <= eta && eta <   0.38) * (   1.72 <= pt && pt <    1.98) * (0.881110) +
      (  0.09 <= eta && eta <   0.38) * (   1.98 <= pt && pt <    2.23) * (0.831849) +
      (  0.09 <= eta && eta <   0.38) * (   2.23 <= pt && pt <    2.49) * (0.788014) +
      (  0.09 <= eta && eta <   0.38) * (   2.49 <= pt && pt <    2.74) * (0.750547) +
      (  0.09 <= eta && eta <   0.38) * (   2.74 <= pt && pt <    3.00) * (0.719558) +
      (  0.38 <= eta && eta <   0.67) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      (  0.38 <= eta && eta <   0.67) * (   0.70 <= pt && pt <    0.95) * (0.999137) +
      (  0.38 <= eta && eta <   0.67) * (   0.95 <= pt && pt <    1.21) * (0.985442) +
      (  0.38 <= eta && eta <   0.67) * (   1.21 <= pt && pt <    1.46) * (0.945774) +
      (  0.38 <= eta && eta <   0.67) * (   1.46 <= pt && pt <    1.72) * (0.893651) +
      (  0.38 <= eta && eta <   0.67) * (   1.72 <= pt && pt <    1.98) * (0.840180) +
      (  0.38 <= eta && eta <   0.67) * (   1.98 <= pt && pt <    2.23) * (0.792935) +
      (  0.38 <= eta && eta <   0.67) * (   2.23 <= pt && pt <    2.49) * (0.751024) +
      (  0.38 <= eta && eta <   0.67) * (   2.49 <= pt && pt <    2.74) * (0.715834) +
      (  0.38 <= eta && eta <   0.67) * (   2.74 <= pt && pt <    3.00) * (0.686635) +
      (  0.67 <= eta && eta <   0.96) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      (  0.67 <= eta && eta <   0.96) * (   0.70 <= pt && pt <    0.95) * (0.994712) +
      (  0.67 <= eta && eta <   0.96) * (   0.95 <= pt && pt <    1.21) * (0.961722) +
      (  0.67 <= eta && eta <   0.96) * (   1.21 <= pt && pt <    1.46) * (0.903052) +
      (  0.67 <= eta && eta <   0.96) * (   1.46 <= pt && pt <    1.72) * (0.842425) +
      (  0.67 <= eta && eta <   0.96) * (   1.72 <= pt && pt <    1.98) * (0.789064) +
      (  0.67 <= eta && eta <   0.96) * (   1.98 <= pt && pt <    2.23) * (0.743371) +
      (  0.67 <= eta && eta <   0.96) * (   2.23 <= pt && pt <    2.49) * (0.706600) +
      (  0.67 <= eta && eta <   0.96) * (   2.49 <= pt && pt <    2.74) * (0.675897) +
      (  0.67 <= eta && eta <   0.96) * (   2.74 <= pt && pt <    3.00) * (0.650807) +
      (  0.96 <= eta && eta <   1.25) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      (  0.96 <= eta && eta <   1.25) * (   0.70 <= pt && pt <    0.95) * (0.976783) +
      (  0.96 <= eta && eta <   1.25) * (   0.95 <= pt && pt <    1.21) * (0.913388) +
      (  0.96 <= eta && eta <   1.25) * (   1.21 <= pt && pt <    1.46) * (0.842033) +
      (  0.96 <= eta && eta <   1.25) * (   1.46 <= pt && pt <    1.72) * (0.780163) +
      (  0.96 <= eta && eta <   1.25) * (   1.72 <= pt && pt <    1.98) * (0.730589) +
      (  0.96 <= eta && eta <   1.25) * (   1.98 <= pt && pt <    2.23) * (0.690807) +
      (  0.96 <= eta && eta <   1.25) * (   2.23 <= pt && pt <    2.49) * (0.659667) +
      (  0.96 <= eta && eta <   1.25) * (   2.49 <= pt && pt <    2.74) * (0.634500) +
      (  0.96 <= eta && eta <   1.25) * (   2.74 <= pt && pt <    3.00) * (0.614794)
    }

    add EfficiencyFormula {-11} {211} {
      (eta< -1.64 || eta>=  1.25 || pt <    0.44 || pt >=    3.00) * ( 0.00 ) +
      ( -1.64 <= eta && eta <  -1.35) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      ( -1.64 <= eta && eta <  -1.35) * (   0.70 <= pt && pt <    0.95) * (0.180687) +
      ( -1.64 <= eta && eta <  -1.35) * (   0.95 <= pt && pt <    1.21) * (0.278302) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.21 <= pt && pt <    1.46) * (0.341809) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.46 <= pt && pt <    1.72) * (0.381893) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.72 <= pt && pt <    1.98) * (0.409109) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.98 <= pt && pt <    2.23) * (0.427981) +
      ( -1.64 <= eta && eta <  -1.35) * (   2.23 <= pt && pt <    2.49) * (0.441588) +
      ( -1.64 <= eta && eta <  -1.35) * (   2.49 <= pt && pt <    2.74) * (0.452079) +
      ( -1.64 <= eta && eta <  -1.35) * (   2.74 <= pt && pt <    3.00) * (0.459736) +
      ( -1.35 <= eta && eta <  -1.06) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      ( -1.35 <= eta && eta <  -1.06) * (   0.70 <= pt && pt <    0.95) * (0.095469) +
      ( -1.35 <= eta && eta <  -1.06) * (   0.95 <= pt && pt <    1.21) * (0.188933) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.21 <= pt && pt <    1.46) * (0.262953) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.46 <= pt && pt <    1.72) * (0.318767) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.72 <= pt && pt <    1.98) * (0.356582) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.98 <= pt && pt <    2.23) * (0.385241) +
      ( -1.35 <= eta && eta <  -1.06) * (   2.23 <= pt && pt <    2.49) * (0.406033) +
      ( -1.35 <= eta && eta <  -1.06) * (   2.49 <= pt && pt <    2.74) * (0.421634) +
      ( -1.35 <= eta && eta <  -1.06) * (   2.74 <= pt && pt <    3.00) * (0.433876) +
      ( -1.06 <= eta && eta <  -0.77) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      ( -1.06 <= eta && eta <  -0.77) * (   0.70 <= pt && pt <    0.95) * (0.035548) +
      ( -1.06 <= eta && eta <  -0.77) * (   0.95 <= pt && pt <    1.21) * (0.110332) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.21 <= pt && pt <    1.46) * (0.185156) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.46 <= pt && pt <    1.72) * (0.247252) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.72 <= pt && pt <    1.98) * (0.296097) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.98 <= pt && pt <    2.23) * (0.333682) +
      ( -1.06 <= eta && eta <  -0.77) * (   2.23 <= pt && pt <    2.49) * (0.362680) +
      ( -1.06 <= eta && eta <  -0.77) * (   2.49 <= pt && pt <    2.74) * (0.385301) +
      ( -1.06 <= eta && eta <  -0.77) * (   2.74 <= pt && pt <    3.00) * (0.402359) +
      ( -0.77 <= eta && eta <  -0.48) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      ( -0.77 <= eta && eta <  -0.48) * (   0.70 <= pt && pt <    0.95) * (0.008285) +
      ( -0.77 <= eta && eta <  -0.48) * (   0.95 <= pt && pt <    1.21) * (0.050219) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.21 <= pt && pt <    1.46) * (0.112182) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.46 <= pt && pt <    1.72) * (0.172923) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.72 <= pt && pt <    1.98) * (0.225572) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.98 <= pt && pt <    2.23) * (0.269134) +
      ( -0.77 <= eta && eta <  -0.48) * (   2.23 <= pt && pt <    2.49) * (0.304540) +
      ( -0.77 <= eta && eta <  -0.48) * (   2.49 <= pt && pt <    2.74) * (0.333881) +
      ( -0.77 <= eta && eta <  -0.48) * (   2.74 <= pt && pt <    3.00) * (0.357698) +
      ( -0.48 <= eta && eta <  -0.20) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      ( -0.48 <= eta && eta <  -0.20) * (   0.70 <= pt && pt <    0.95) * (0.001651) +
      ( -0.48 <= eta && eta <  -0.20) * (   0.95 <= pt && pt <    1.21) * (0.020771) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.21 <= pt && pt <    1.46) * (0.063586) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.46 <= pt && pt <    1.72) * (0.119232) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.72 <= pt && pt <    1.98) * (0.171645) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.98 <= pt && pt <    2.23) * (0.219382) +
      ( -0.48 <= eta && eta <  -0.20) * (   2.23 <= pt && pt <    2.49) * (0.260306) +
      ( -0.48 <= eta && eta <  -0.20) * (   2.49 <= pt && pt <    2.74) * (0.293508) +
      ( -0.48 <= eta && eta <  -0.20) * (   2.74 <= pt && pt <    3.00) * (0.321097) +
      ( -0.20 <= eta && eta <   0.09) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      ( -0.20 <= eta && eta <   0.09) * (   0.70 <= pt && pt <    0.95) * (0.000248) +
      ( -0.20 <= eta && eta <   0.09) * (   0.95 <= pt && pt <    1.21) * (0.007141) +
      ( -0.20 <= eta && eta <   0.09) * (   1.21 <= pt && pt <    1.46) * (0.034679) +
      ( -0.20 <= eta && eta <   0.09) * (   1.46 <= pt && pt <    1.72) * (0.079354) +
      ( -0.20 <= eta && eta <   0.09) * (   1.72 <= pt && pt <    1.98) * (0.128061) +
      ( -0.20 <= eta && eta <   0.09) * (   1.98 <= pt && pt <    2.23) * (0.177204) +
      ( -0.20 <= eta && eta <   0.09) * (   2.23 <= pt && pt <    2.49) * (0.220486) +
      ( -0.20 <= eta && eta <   0.09) * (   2.49 <= pt && pt <    2.74) * (0.256869) +
      ( -0.20 <= eta && eta <   0.09) * (   2.74 <= pt && pt <    3.00) * (0.288559) +
      (  0.09 <= eta && eta <   0.38) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      (  0.09 <= eta && eta <   0.38) * (   0.70 <= pt && pt <    0.95) * (0.000401) +
      (  0.09 <= eta && eta <   0.38) * (   0.95 <= pt && pt <    1.21) * (0.005337) +
      (  0.09 <= eta && eta <   0.38) * (   1.21 <= pt && pt <    1.46) * (0.028586) +
      (  0.09 <= eta && eta <   0.38) * (   1.46 <= pt && pt <    1.72) * (0.069951) +
      (  0.09 <= eta && eta <   0.38) * (   1.72 <= pt && pt <    1.98) * (0.118890) +
      (  0.09 <= eta && eta <   0.38) * (   1.98 <= pt && pt <    2.23) * (0.168151) +
      (  0.09 <= eta && eta <   0.38) * (   2.23 <= pt && pt <    2.49) * (0.211986) +
      (  0.09 <= eta && eta <   0.38) * (   2.49 <= pt && pt <    2.74) * (0.249453) +
      (  0.09 <= eta && eta <   0.38) * (   2.74 <= pt && pt <    3.00) * (0.280442) +
      (  0.38 <= eta && eta <   0.67) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      (  0.38 <= eta && eta <   0.67) * (   0.70 <= pt && pt <    0.95) * (0.000863) +
      (  0.38 <= eta && eta <   0.67) * (   0.95 <= pt && pt <    1.21) * (0.014558) +
      (  0.38 <= eta && eta <   0.67) * (   1.21 <= pt && pt <    1.46) * (0.054226) +
      (  0.38 <= eta && eta <   0.67) * (   1.46 <= pt && pt <    1.72) * (0.106349) +
      (  0.38 <= eta && eta <   0.67) * (   1.72 <= pt && pt <    1.98) * (0.159820) +
      (  0.38 <= eta && eta <   0.67) * (   1.98 <= pt && pt <    2.23) * (0.207065) +
      (  0.38 <= eta && eta <   0.67) * (   2.23 <= pt && pt <    2.49) * (0.248976) +
      (  0.38 <= eta && eta <   0.67) * (   2.49 <= pt && pt <    2.74) * (0.284166) +
      (  0.38 <= eta && eta <   0.67) * (   2.74 <= pt && pt <    3.00) * (0.313365) +
      (  0.67 <= eta && eta <   0.96) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      (  0.67 <= eta && eta <   0.96) * (   0.70 <= pt && pt <    0.95) * (0.005288) +
      (  0.67 <= eta && eta <   0.96) * (   0.95 <= pt && pt <    1.21) * (0.038278) +
      (  0.67 <= eta && eta <   0.96) * (   1.21 <= pt && pt <    1.46) * (0.096948) +
      (  0.67 <= eta && eta <   0.96) * (   1.46 <= pt && pt <    1.72) * (0.157575) +
      (  0.67 <= eta && eta <   0.96) * (   1.72 <= pt && pt <    1.98) * (0.210936) +
      (  0.67 <= eta && eta <   0.96) * (   1.98 <= pt && pt <    2.23) * (0.256629) +
      (  0.67 <= eta && eta <   0.96) * (   2.23 <= pt && pt <    2.49) * (0.293400) +
      (  0.67 <= eta && eta <   0.96) * (   2.49 <= pt && pt <    2.74) * (0.324103) +
      (  0.67 <= eta && eta <   0.96) * (   2.74 <= pt && pt <    3.00) * (0.349193) +
      (  0.96 <= eta && eta <   1.25) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      (  0.96 <= eta && eta <   1.25) * (   0.70 <= pt && pt <    0.95) * (0.023217) +
      (  0.96 <= eta && eta <   1.25) * (   0.95 <= pt && pt <    1.21) * (0.086612) +
      (  0.96 <= eta && eta <   1.25) * (   1.21 <= pt && pt <    1.46) * (0.157967) +
      (  0.96 <= eta && eta <   1.25) * (   1.46 <= pt && pt <    1.72) * (0.219837) +
      (  0.96 <= eta && eta <   1.25) * (   1.72 <= pt && pt <    1.98) * (0.269411) +
      (  0.96 <= eta && eta <   1.25) * (   1.98 <= pt && pt <    2.23) * (0.309193) +
      (  0.96 <= eta && eta <   1.25) * (   2.23 <= pt && pt <    2.49) * (0.340333) +
      (  0.96 <= eta && eta <   1.25) * (   2.49 <= pt && pt <    2.74) * (0.365500) +
      (  0.96 <= eta && eta <   1.25) * (   2.74 <= pt && pt <    3.00) * (0.385206)
    }

    add EfficiencyFormula {211} {-11} {
      (eta< -1.64 || eta>=  1.25 || pt <    0.44 || pt >=    3.00) * ( 0.00 ) +
      ( -1.64 <= eta && eta <  -1.35) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      ( -1.64 <= eta && eta <  -1.35) * (   0.70 <= pt && pt <    0.95) * (0.180687) +
      ( -1.64 <= eta && eta <  -1.35) * (   0.95 <= pt && pt <    1.21) * (0.278302) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.21 <= pt && pt <    1.46) * (0.341809) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.46 <= pt && pt <    1.72) * (0.381893) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.72 <= pt && pt <    1.98) * (0.409109) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.98 <= pt && pt <    2.23) * (0.427981) +
      ( -1.64 <= eta && eta <  -1.35) * (   2.23 <= pt && pt <    2.49) * (0.441588) +
      ( -1.64 <= eta && eta <  -1.35) * (   2.49 <= pt && pt <    2.74) * (0.452079) +
      ( -1.64 <= eta && eta <  -1.35) * (   2.74 <= pt && pt <    3.00) * (0.459736) +
      ( -1.35 <= eta && eta <  -1.06) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      ( -1.35 <= eta && eta <  -1.06) * (   0.70 <= pt && pt <    0.95) * (0.095469) +
      ( -1.35 <= eta && eta <  -1.06) * (   0.95 <= pt && pt <    1.21) * (0.188933) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.21 <= pt && pt <    1.46) * (0.262953) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.46 <= pt && pt <    1.72) * (0.318767) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.72 <= pt && pt <    1.98) * (0.356582) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.98 <= pt && pt <    2.23) * (0.385241) +
      ( -1.35 <= eta && eta <  -1.06) * (   2.23 <= pt && pt <    2.49) * (0.406033) +
      ( -1.35 <= eta && eta <  -1.06) * (   2.49 <= pt && pt <    2.74) * (0.421634) +
      ( -1.35 <= eta && eta <  -1.06) * (   2.74 <= pt && pt <    3.00) * (0.433876) +
      ( -1.06 <= eta && eta <  -0.77) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      ( -1.06 <= eta && eta <  -0.77) * (   0.70 <= pt && pt <    0.95) * (0.035548) +
      ( -1.06 <= eta && eta <  -0.77) * (   0.95 <= pt && pt <    1.21) * (0.110332) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.21 <= pt && pt <    1.46) * (0.185156) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.46 <= pt && pt <    1.72) * (0.247252) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.72 <= pt && pt <    1.98) * (0.296097) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.98 <= pt && pt <    2.23) * (0.333682) +
      ( -1.06 <= eta && eta <  -0.77) * (   2.23 <= pt && pt <    2.49) * (0.362680) +
      ( -1.06 <= eta && eta <  -0.77) * (   2.49 <= pt && pt <    2.74) * (0.385301) +
      ( -1.06 <= eta && eta <  -0.77) * (   2.74 <= pt && pt <    3.00) * (0.402359) +
      ( -0.77 <= eta && eta <  -0.48) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      ( -0.77 <= eta && eta <  -0.48) * (   0.70 <= pt && pt <    0.95) * (0.008285) +
      ( -0.77 <= eta && eta <  -0.48) * (   0.95 <= pt && pt <    1.21) * (0.050219) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.21 <= pt && pt <    1.46) * (0.112182) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.46 <= pt && pt <    1.72) * (0.172923) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.72 <= pt && pt <    1.98) * (0.225572) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.98 <= pt && pt <    2.23) * (0.269134) +
      ( -0.77 <= eta && eta <  -0.48) * (   2.23 <= pt && pt <    2.49) * (0.304540) +
      ( -0.77 <= eta && eta <  -0.48) * (   2.49 <= pt && pt <    2.74) * (0.333881) +
      ( -0.77 <= eta && eta <  -0.48) * (   2.74 <= pt && pt <    3.00) * (0.357698) +
      ( -0.48 <= eta && eta <  -0.20) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      ( -0.48 <= eta && eta <  -0.20) * (   0.70 <= pt && pt <    0.95) * (0.001651) +
      ( -0.48 <= eta && eta <  -0.20) * (   0.95 <= pt && pt <    1.21) * (0.020771) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.21 <= pt && pt <    1.46) * (0.063586) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.46 <= pt && pt <    1.72) * (0.119232) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.72 <= pt && pt <    1.98) * (0.171645) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.98 <= pt && pt <    2.23) * (0.219382) +
      ( -0.48 <= eta && eta <  -0.20) * (   2.23 <= pt && pt <    2.49) * (0.260306) +
      ( -0.48 <= eta && eta <  -0.20) * (   2.49 <= pt && pt <    2.74) * (0.293508) +
      ( -0.48 <= eta && eta <  -0.20) * (   2.74 <= pt && pt <    3.00) * (0.321097) +
      ( -0.20 <= eta && eta <   0.09) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      ( -0.20 <= eta && eta <   0.09) * (   0.70 <= pt && pt <    0.95) * (0.000248) +
      ( -0.20 <= eta && eta <   0.09) * (   0.95 <= pt && pt <    1.21) * (0.007141) +
      ( -0.20 <= eta && eta <   0.09) * (   1.21 <= pt && pt <    1.46) * (0.034679) +
      ( -0.20 <= eta && eta <   0.09) * (   1.46 <= pt && pt <    1.72) * (0.079354) +
      ( -0.20 <= eta && eta <   0.09) * (   1.72 <= pt && pt <    1.98) * (0.128061) +
      ( -0.20 <= eta && eta <   0.09) * (   1.98 <= pt && pt <    2.23) * (0.177204) +
      ( -0.20 <= eta && eta <   0.09) * (   2.23 <= pt && pt <    2.49) * (0.220486) +
      ( -0.20 <= eta && eta <   0.09) * (   2.49 <= pt && pt <    2.74) * (0.256869) +
      ( -0.20 <= eta && eta <   0.09) * (   2.74 <= pt && pt <    3.00) * (0.288559) +
      (  0.09 <= eta && eta <   0.38) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      (  0.09 <= eta && eta <   0.38) * (   0.70 <= pt && pt <    0.95) * (0.000401) +
      (  0.09 <= eta && eta <   0.38) * (   0.95 <= pt && pt <    1.21) * (0.005337) +
      (  0.09 <= eta && eta <   0.38) * (   1.21 <= pt && pt <    1.46) * (0.028586) +
      (  0.09 <= eta && eta <   0.38) * (   1.46 <= pt && pt <    1.72) * (0.069951) +
      (  0.09 <= eta && eta <   0.38) * (   1.72 <= pt && pt <    1.98) * (0.118890) +
      (  0.09 <= eta && eta <   0.38) * (   1.98 <= pt && pt <    2.23) * (0.168151) +
      (  0.09 <= eta && eta <   0.38) * (   2.23 <= pt && pt <    2.49) * (0.211986) +
      (  0.09 <= eta && eta <   0.38) * (   2.49 <= pt && pt <    2.74) * (0.249453) +
      (  0.09 <= eta && eta <   0.38) * (   2.74 <= pt && pt <    3.00) * (0.280442) +
      (  0.38 <= eta && eta <   0.67) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      (  0.38 <= eta && eta <   0.67) * (   0.70 <= pt && pt <    0.95) * (0.000863) +
      (  0.38 <= eta && eta <   0.67) * (   0.95 <= pt && pt <    1.21) * (0.014558) +
      (  0.38 <= eta && eta <   0.67) * (   1.21 <= pt && pt <    1.46) * (0.054226) +
      (  0.38 <= eta && eta <   0.67) * (   1.46 <= pt && pt <    1.72) * (0.106349) +
      (  0.38 <= eta && eta <   0.67) * (   1.72 <= pt && pt <    1.98) * (0.159820) +
      (  0.38 <= eta && eta <   0.67) * (   1.98 <= pt && pt <    2.23) * (0.207065) +
      (  0.38 <= eta && eta <   0.67) * (   2.23 <= pt && pt <    2.49) * (0.248976) +
      (  0.38 <= eta && eta <   0.67) * (   2.49 <= pt && pt <    2.74) * (0.284166) +
      (  0.38 <= eta && eta <   0.67) * (   2.74 <= pt && pt <    3.00) * (0.313365) +
      (  0.67 <= eta && eta <   0.96) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      (  0.67 <= eta && eta <   0.96) * (   0.70 <= pt && pt <    0.95) * (0.005288) +
      (  0.67 <= eta && eta <   0.96) * (   0.95 <= pt && pt <    1.21) * (0.038278) +
      (  0.67 <= eta && eta <   0.96) * (   1.21 <= pt && pt <    1.46) * (0.096948) +
      (  0.67 <= eta && eta <   0.96) * (   1.46 <= pt && pt <    1.72) * (0.157575) +
      (  0.67 <= eta && eta <   0.96) * (   1.72 <= pt && pt <    1.98) * (0.210936) +
      (  0.67 <= eta && eta <   0.96) * (   1.98 <= pt && pt <    2.23) * (0.256629) +
      (  0.67 <= eta && eta <   0.96) * (   2.23 <= pt && pt <    2.49) * (0.293400) +
      (  0.67 <= eta && eta <   0.96) * (   2.49 <= pt && pt <    2.74) * (0.324103) +
      (  0.67 <= eta && eta <   0.96) * (   2.74 <= pt && pt <    3.00) * (0.349193) +
      (  0.96 <= eta && eta <   1.25) * (   0.44 <= pt && pt <    0.70) * (0.000000) +
      (  0.96 <= eta && eta <   1.25) * (   0.70 <= pt && pt <    0.95) * (0.023217) +
      (  0.96 <= eta && eta <   1.25) * (   0.95 <= pt && pt <    1.21) * (0.086612) +
      (  0.96 <= eta && eta <   1.25) * (   1.21 <= pt && pt <    1.46) * (0.157967) +
      (  0.96 <= eta && eta <   1.25) * (   1.46 <= pt && pt <    1.72) * (0.219837) +
      (  0.96 <= eta && eta <   1.25) * (   1.72 <= pt && pt <    1.98) * (0.269411) +
      (  0.96 <= eta && eta <   1.25) * (   1.98 <= pt && pt <    2.23) * (0.309193) +
      (  0.96 <= eta && eta <   1.25) * (   2.23 <= pt && pt <    2.49) * (0.340333) +
      (  0.96 <= eta && eta <   1.25) * (   2.49 <= pt && pt <    2.74) * (0.365500) +
      (  0.96 <= eta && eta <   1.25) * (   2.74 <= pt && pt <    3.00) * (0.385206)
    }

    add EfficiencyFormula {211} {211} {
      (eta< -1.64 || eta>=  1.25 || pt <    0.44 || pt >=    3.00) * ( 0.00 ) +
      ( -1.64 <= eta && eta <  -1.35) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      ( -1.64 <= eta && eta <  -1.35) * (   0.70 <= pt && pt <    0.95) * (0.819313) +
      ( -1.64 <= eta && eta <  -1.35) * (   0.95 <= pt && pt <    1.21) * (0.721698) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.21 <= pt && pt <    1.46) * (0.658191) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.46 <= pt && pt <    1.72) * (0.618107) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.72 <= pt && pt <    1.98) * (0.590891) +
      ( -1.64 <= eta && eta <  -1.35) * (   1.98 <= pt && pt <    2.23) * (0.572019) +
      ( -1.64 <= eta && eta <  -1.35) * (   2.23 <= pt && pt <    2.49) * (0.558412) +
      ( -1.64 <= eta && eta <  -1.35) * (   2.49 <= pt && pt <    2.74) * (0.547921) +
      ( -1.64 <= eta && eta <  -1.35) * (   2.74 <= pt && pt <    3.00) * (0.540264) +
      ( -1.35 <= eta && eta <  -1.06) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      ( -1.35 <= eta && eta <  -1.06) * (   0.70 <= pt && pt <    0.95) * (0.904531) +
      ( -1.35 <= eta && eta <  -1.06) * (   0.95 <= pt && pt <    1.21) * (0.811067) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.21 <= pt && pt <    1.46) * (0.737047) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.46 <= pt && pt <    1.72) * (0.681233) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.72 <= pt && pt <    1.98) * (0.643418) +
      ( -1.35 <= eta && eta <  -1.06) * (   1.98 <= pt && pt <    2.23) * (0.614759) +
      ( -1.35 <= eta && eta <  -1.06) * (   2.23 <= pt && pt <    2.49) * (0.593967) +
      ( -1.35 <= eta && eta <  -1.06) * (   2.49 <= pt && pt <    2.74) * (0.578366) +
      ( -1.35 <= eta && eta <  -1.06) * (   2.74 <= pt && pt <    3.00) * (0.566124) +
      ( -1.06 <= eta && eta <  -0.77) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      ( -1.06 <= eta && eta <  -0.77) * (   0.70 <= pt && pt <    0.95) * (0.964452) +
      ( -1.06 <= eta && eta <  -0.77) * (   0.95 <= pt && pt <    1.21) * (0.889668) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.21 <= pt && pt <    1.46) * (0.814844) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.46 <= pt && pt <    1.72) * (0.752748) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.72 <= pt && pt <    1.98) * (0.703903) +
      ( -1.06 <= eta && eta <  -0.77) * (   1.98 <= pt && pt <    2.23) * (0.666318) +
      ( -1.06 <= eta && eta <  -0.77) * (   2.23 <= pt && pt <    2.49) * (0.637320) +
      ( -1.06 <= eta && eta <  -0.77) * (   2.49 <= pt && pt <    2.74) * (0.614699) +
      ( -1.06 <= eta && eta <  -0.77) * (   2.74 <= pt && pt <    3.00) * (0.597641) +
      ( -0.77 <= eta && eta <  -0.48) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      ( -0.77 <= eta && eta <  -0.48) * (   0.70 <= pt && pt <    0.95) * (0.991715) +
      ( -0.77 <= eta && eta <  -0.48) * (   0.95 <= pt && pt <    1.21) * (0.949781) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.21 <= pt && pt <    1.46) * (0.887818) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.46 <= pt && pt <    1.72) * (0.827077) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.72 <= pt && pt <    1.98) * (0.774428) +
      ( -0.77 <= eta && eta <  -0.48) * (   1.98 <= pt && pt <    2.23) * (0.730866) +
      ( -0.77 <= eta && eta <  -0.48) * (   2.23 <= pt && pt <    2.49) * (0.695460) +
      ( -0.77 <= eta && eta <  -0.48) * (   2.49 <= pt && pt <    2.74) * (0.666119) +
      ( -0.77 <= eta && eta <  -0.48) * (   2.74 <= pt && pt <    3.00) * (0.642302) +
      ( -0.48 <= eta && eta <  -0.20) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      ( -0.48 <= eta && eta <  -0.20) * (   0.70 <= pt && pt <    0.95) * (0.998349) +
      ( -0.48 <= eta && eta <  -0.20) * (   0.95 <= pt && pt <    1.21) * (0.979229) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.21 <= pt && pt <    1.46) * (0.936414) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.46 <= pt && pt <    1.72) * (0.880768) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.72 <= pt && pt <    1.98) * (0.828355) +
      ( -0.48 <= eta && eta <  -0.20) * (   1.98 <= pt && pt <    2.23) * (0.780618) +
      ( -0.48 <= eta && eta <  -0.20) * (   2.23 <= pt && pt <    2.49) * (0.739694) +
      ( -0.48 <= eta && eta <  -0.20) * (   2.49 <= pt && pt <    2.74) * (0.706492) +
      ( -0.48 <= eta && eta <  -0.20) * (   2.74 <= pt && pt <    3.00) * (0.678903) +
      ( -0.20 <= eta && eta <   0.09) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      ( -0.20 <= eta && eta <   0.09) * (   0.70 <= pt && pt <    0.95) * (0.999752) +
      ( -0.20 <= eta && eta <   0.09) * (   0.95 <= pt && pt <    1.21) * (0.992859) +
      ( -0.20 <= eta && eta <   0.09) * (   1.21 <= pt && pt <    1.46) * (0.965321) +
      ( -0.20 <= eta && eta <   0.09) * (   1.46 <= pt && pt <    1.72) * (0.920646) +
      ( -0.20 <= eta && eta <   0.09) * (   1.72 <= pt && pt <    1.98) * (0.871939) +
      ( -0.20 <= eta && eta <   0.09) * (   1.98 <= pt && pt <    2.23) * (0.822796) +
      ( -0.20 <= eta && eta <   0.09) * (   2.23 <= pt && pt <    2.49) * (0.779514) +
      ( -0.20 <= eta && eta <   0.09) * (   2.49 <= pt && pt <    2.74) * (0.743131) +
      ( -0.20 <= eta && eta <   0.09) * (   2.74 <= pt && pt <    3.00) * (0.711441) +
      (  0.09 <= eta && eta <   0.38) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      (  0.09 <= eta && eta <   0.38) * (   0.70 <= pt && pt <    0.95) * (0.999599) +
      (  0.09 <= eta && eta <   0.38) * (   0.95 <= pt && pt <    1.21) * (0.994663) +
      (  0.09 <= eta && eta <   0.38) * (   1.21 <= pt && pt <    1.46) * (0.971414) +
      (  0.09 <= eta && eta <   0.38) * (   1.46 <= pt && pt <    1.72) * (0.930049) +
      (  0.09 <= eta && eta <   0.38) * (   1.72 <= pt && pt <    1.98) * (0.881110) +
      (  0.09 <= eta && eta <   0.38) * (   1.98 <= pt && pt <    2.23) * (0.831849) +
      (  0.09 <= eta && eta <   0.38) * (   2.23 <= pt && pt <    2.49) * (0.788014) +
      (  0.09 <= eta && eta <   0.38) * (   2.49 <= pt && pt <    2.74) * (0.750547) +
      (  0.09 <= eta && eta <   0.38) * (   2.74 <= pt && pt <    3.00) * (0.719558) +
      (  0.38 <= eta && eta <   0.67) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      (  0.38 <= eta && eta <   0.67) * (   0.70 <= pt && pt <    0.95) * (0.999137) +
      (  0.38 <= eta && eta <   0.67) * (   0.95 <= pt && pt <    1.21) * (0.985442) +
      (  0.38 <= eta && eta <   0.67) * (   1.21 <= pt && pt <    1.46) * (0.945774) +
      (  0.38 <= eta && eta <   0.67) * (   1.46 <= pt && pt <    1.72) * (0.893651) +
      (  0.38 <= eta && eta <   0.67) * (   1.72 <= pt && pt <    1.98) * (0.840180) +
      (  0.38 <= eta && eta <   0.67) * (   1.98 <= pt && pt <    2.23) * (0.792935) +
      (  0.38 <= eta && eta <   0.67) * (   2.23 <= pt && pt <    2.49) * (0.751024) +
      (  0.38 <= eta && eta <   0.67) * (   2.49 <= pt && pt <    2.74) * (0.715834) +
      (  0.38 <= eta && eta <   0.67) * (   2.74 <= pt && pt <    3.00) * (0.686635) +
      (  0.67 <= eta && eta <   0.96) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      (  0.67 <= eta && eta <   0.96) * (   0.70 <= pt && pt <    0.95) * (0.994712) +
      (  0.67 <= eta && eta <   0.96) * (   0.95 <= pt && pt <    1.21) * (0.961722) +
      (  0.67 <= eta && eta <   0.96) * (   1.21 <= pt && pt <    1.46) * (0.903052) +
      (  0.67 <= eta && eta <   0.96) * (   1.46 <= pt && pt <    1.72) * (0.842425) +
      (  0.67 <= eta && eta <   0.96) * (   1.72 <= pt && pt <    1.98) * (0.789064) +
      (  0.67 <= eta && eta <   0.96) * (   1.98 <= pt && pt <    2.23) * (0.743371) +
      (  0.67 <= eta && eta <   0.96) * (   2.23 <= pt && pt <    2.49) * (0.706600) +
      (  0.67 <= eta && eta <   0.96) * (   2.49 <= pt && pt <    2.74) * (0.675897) +
      (  0.67 <= eta && eta <   0.96) * (   2.74 <= pt && pt <    3.00) * (0.650807) +
      (  0.96 <= eta && eta <   1.25) * (   0.44 <= pt && pt <    0.70) * (1.000000) +
      (  0.96 <= eta && eta <   1.25) * (   0.70 <= pt && pt <    0.95) * (0.976783) +
      (  0.96 <= eta && eta <   1.25) * (   0.95 <= pt && pt <    1.21) * (0.913388) +
      (  0.96 <= eta && eta <   1.25) * (   1.21 <= pt && pt <    1.46) * (0.842033) +
      (  0.96 <= eta && eta <   1.25) * (   1.46 <= pt && pt <    1.72) * (0.780163) +
      (  0.96 <= eta && eta <   1.25) * (   1.72 <= pt && pt <    1.98) * (0.730589) +
      (  0.96 <= eta && eta <   1.25) * (   1.98 <= pt && pt <    2.23) * (0.690807) +
      (  0.96 <= eta && eta <   1.25) * (   2.23 <= pt && pt <    2.49) * (0.659667) +
      (  0.96 <= eta && eta <   1.25) * (   2.49 <= pt && pt <    2.74) * (0.634500) +
      (  0.96 <= eta && eta <   1.25) * (   2.74 <= pt && pt <    3.00) * (0.614794)
    }

  add EfficiencyFormula {0} {0} { 0.00 }
}