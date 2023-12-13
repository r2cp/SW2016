This replication file contains data and programs for the paper

"Core Inflation and Trend Inflation"

by James Stock and Mark Watson

Revised January 2016

Data:  Data are contained in the subdirectory "data"

Programs:  All calculations are carried out using MATLAB (R2015b)

Matlab progams for the project are in the "matlab" directory.
These use a variety of utility programs that are in the directory "matlab_utility_programs" -- make sure this directory is in your MATLAB path.

UCSV Models:

The full-sample UCSV for the three univariate series are computed in the program:
ucsv_dp.m

The full-sample 17-component multivariate model is computed in 
mucsv_tvp_17c.m:  Multivariate 17-component model

The full-sample 3-component multivariate model is computed in 
mucsv_tvp_3c.m: Multivariate 3-component model

The recursively-estimated UCSV models used for POOS forecasting and shown in Table 2(b) and figure 9 are computed in 
ucsv_dp_poos.m
mucsv_tvp_17c_poos.m
mucsv_tvp_3c_poos.m

Constant-alpha versions of the component models are computed in 
mucsv_constAlpha_17c.m
mucsv_constAlpha_3c.m
mucsv_constAlpha_17c_poos.m
mucsv_constAlpha_3c_poos.m

Results reported in Paper:
 Table 1:  Tab_1.m
 Table 2:  Figure_8_Table2a.m and Figure_9_Table2b.m
 Table 3:  Tab_3.m
 
 Figure 1: Figure_1_2.m
 Figure 2: Figure_1_2.m
 Figure 3: Figure_3_4.m
 Figure 4: Figure_3_4.m
 Figure 5: Figure_5_6.m
 Figure 6: Figure_5_6.m
 Figure 7: Figure_7.m
 Figure 8: Figure_8_Table2a.m
 Figure 9: Figure_9_Table2b.m
 
 
 Results Reported in the online appendix:
 
 Tables A.1 and A.2:  AppendixTabs_A1_A2.m
 Tables A.3 - A.7: AppendixTabs_A3_A4_A5_A6_A7.m
 Tables A.8 - A.12:  AppendixTabs_A8_A9_A10_A11_A12.m
 Table A.13:
  The monthly recursive esimates are computed in
   ucsv_dp_poos_monthly.m
   mucsv_tvp_17c_poos_monthly.m
   mucsv_tvp_3c_poos_monthly.m
  The monthly values reported in the table are from AppendixTab_13.  The quarterly values are from Tab_3.m (these were reported in the paper).
 
 Appendix Figures A.1 - A.17: Appendix_mucsv17_figures.m
 Appendix Figures A.18 - A.20: Appendix_mucsv3_figures.m