# Analyzing Core Inflation and Trend Inflation: A Guatemalan Perspective

Replication files the short research paper for the Macroeconometrics class. This file outlines the main files used to produce the paper. The original codebase is based on the replication files by Stock and Watson (2016), and is contained in the `main` branch. For the corresponding results of the short paper, check the other git branches.

## Baseline Results

 The git branch used for this purpose is `gt/baseline`.  

- `pcomp_data_calendar_m_and_q_gt.m`: Main script to load Guatemalan time series and disaggregated data. 
- `ucsv_dp_gt_M.m`: Monthly univariate UCSVO model estimation for the CPI, CPIxE and CPIxFE. 
- `ucsv_dp_gt_Q.m`: Quarterly univariate UCSVO model estimation for the CPI, CPIxE and CPIxFE. 
- `mucsv_tvp_12c_gt_M.m` Monthly multivariate UCSVO model estimation for the 12 major CPI divisions. 

## UCSVO Model Variants 

In the corresponding git branches (`gt/no-jumps`, `gt/constant_vol_eps` and `gt/constant_vol_dtau`), the file `ucsv_dp_gt_M.m` estimates a monthly univariate UCSVO model estimation for the CPI, CPIxE and CPIxFE with the following feature modifications: 

- `ucsv_outlier_const_vol_e.m`: Function to estimate the univariate UCSVO with constant volatility of the transitory component. 
- `ucsv_outlier_const_vol_dtau.m`: Function to estimate the univariate UCSVO with constant volatility of the permanent or trend component. 
- `ucsv_outlier_nojumps.m`: Function to estimate the univariate UCSVO without the outlier adjusting feature of $s_t$.
- `Figure_1_2_gtM.m`: Script that produces the graph with the main objects from the univariate UCSVO model. 
- `Figure_2_compare_tau_gtM.m`: Script that produces the comparison graphs across model variants for the main objects of the UCSVO model. 
- `Figure_3_4_gt.m`: Script that produces the graph of the main multivariate UCSVO objects. 
 

## Out-of-Sample Forecasting Exercise 

In the branch `gt/poos-exercise`, the following files are used to estimate the UCSVO in a rolling-window fashion to perform the out-of-sample forecasting exercise: 

- `ucsv_dp_poos_monthly_gt_baseline.m`: Estimates the original fully-featured UCSVO model for the forecasting exercise. 
- `ucsv_dp_poos_monthly_gt_const_vol_dtau.m`: Estimates the UCSVO model with constant volatility of the permanent component for the forecasting exercise. 
- `ucsv_dp_poos_monthly_gt_const_vol_e.m`: Estimates the UCSVO model with constant volatility of the transitory component for the forecasting exercise. 
- `ucsv_dp_poos_monthly_gt_nojumps.m`: Estimates the UCSVO model without the outlier adjusting feature. 
- `Table_3_gt.m`: It produces the forecasting exercise results, based on the same exercise by Stock and Watson. 