% Appendix Tables A.1 and A.2
% 11/4/2015, MWW
%


clear all;
small = 1.0e-10;
big = 1.0e+6;                  
rng(663436);
 
  % -- File Directories  
  outdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/out/';
  figdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/fig/';
  matdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/mat/';

  % -- Read in Data --- 
  load_data = 1;  % 1 if reloading data from Excel, etc 
  mtoq_agg = 3;   % Temporal aggregation indicator of monthly to quarterly data
  pcomp_data_calendar_m_and_q;
  % Data Series Used
  dp_agg = dp_agg_q;
  dp_agg_xfe = dp_agg_xfe_q;
  dp_agg_xe = dp_agg_xe_q;
  dp_disagg = dp_disagg_q;
  calvec = calvec_q;
  dnobs = dnobs_q;
  calds = calds_q;
  nper = 4;
  outlabel = 'Q3';
  
  labvec = {'Headline Inflation';'Core XFE';'Core XE'};
  namevec = {'dp_agg';'dp_xfe';'dp_xe'};
  dp = [dp_agg dp_agg_xfe dp_agg_xe];
  
  outfile_name = [outdir 'AppendixTabs_A1_A2.out'];
  fileID = fopen(outfile_name,'w');
  for iseries = 1:size(dp,2);
    tmp = char(namevec(iseries));
    ulabel = [tmp '_' outlabel]; 
    fprintf(fileID,['Carrying out Calculations for ' ulabel '\n']);
    str_tmp = [matdir 'g_eps_prior_post' ulabel]; load(str_tmp);
    str_tmp = [matdir 'g_dtau_prior_post' ulabel]; load(str_tmp);
    str_tmp = [matdir 'ps_mean_pct' ulabel]; load(str_tmp);
    p = 1-ps_mean_pct;
    fprintf(fileID,'\n Gamma_EPS \n');
    fprintf(fileID,'gamma,prior,post \n');
    prtmat_comma(g_eps_prior_post,fileID,'%5.2f','\n');
    fprintf(fileID,'\n Gamma_dtau \n');
    fprintf(fileID,'gamma,prior,post \n');
    prtmat_comma(g_dtau_prior_post,fileID,'%5.2f','\n');
    fprintf(fileID,'\n Outlier Prob \n');
    fprintf(fileID,'16,50,67 \n');
    prtmat_comma(p(:,[4 3 2]),fileID,'%5.2f','\n');
  end;
  
