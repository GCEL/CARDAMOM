pars_lims = {0:  [1e-3,  0.1],     # Decomp rate [1e-5, 0.01]
1:  [0.43,   0.48],    # GPP to resp fraction [~0.46]
2:  [0.75,   1.5],     # GSI sens leaf growth [1.0, 1.025]
3:  [0.10,   1.0],     # NPP belowground allocation exponential parameter [0.01, 1.00]
4:  [1e-3,   2.0],     # GSI max leaf turnover [1e-5, 0.2]
5:  [1e-3,   1e-1],    # TOR roots [0.0001, 0.01]
6:  [1e-3,   1e-1],    # TOR litter [0.0001, 0.01]
7:  [1e-7,   1e-4],    # TOR SOM [1e-7, 0.001]
8:  [0.01,   0.20],    # T factor (Q10) [0.018,  0.08]
9:  [7,      25],      # PNUE [7, 20]
10: [1e-3,   1.0],     # GSI max labile turnover [1e-6, 0.2]
11: [230,    290],     # GSI min T (K) [225, 330]
12: [250,    300],     # GSI max T (K) [225, 330]
13: [3600,   20000],   # GSI min photoperiod (sec) [3600, 36000]
14: [35,     55],      # Leaf Mass per Area [20, 60]
15: [20,     100],     # initial labile pool size [1, 1000]
16: [20,     100],     # initial foliar pool size [1, 1000]
17: [40,     2000],    # initial root pool size [1, 1000]
18: [40,     2000],    # initial litter pool size [1, 10000]
19: [10000,  40000],   # GSI max photoperiod (sec) [3600, 64800]
20: [100,    3000],    # GSI min VPD (Pa) [1, 5500]
21: [1000,   5000],    # GSI max VPD (Pa) [1, 5500]
22: [1e-3,   0.5],     # critical GPP for LAI growth [1e-10, 0.30]
23: [0.96,   1.00],    # GSI sens for leaf senescenece [0.96, 1.00]
24: [0.5,    3.0],     # GSI growing stage/step [0.50, 1.5]
25: [1.0,    2.0],     # Initial GSI [1.0, 2.0]
26: [500,    1500],    # DM min lim for grazing (kg.DM.ha-1)
27: [1500,   3000],    # DM min lim for cutting (kg.DM.ha-1)
28: [0.25,   0.75],    # leaf:stem allocation [0.05, 0.75]
29: [19000,  21000],   # initial SOM pool size [5000, 10000] (UK) 19000, 21000
30: [0.015,  0.035],   # livestock demand in DM (1-3% of animal weight)
31: [0.01,   0.10],    # Post-grazing labile loss (fraction)
32: [0.50,   0.90],    # Post-cutting labile loss (fraction)
33: [0.1,    1.0]}     # min DM removal for grazing instance to occur (g.C.m-2.w-1)

for k, v in pars_lims.items():
    print(f"PI%parmin({k + 1}) = {v[0]}d0 \nPI%parmax({k + 1}) = {v[1]}d0 \n")
    # break
