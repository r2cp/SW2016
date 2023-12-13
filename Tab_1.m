% Summary Calculations 
% May 17, 2015
%


clear all;
small = 1.0e-10;
big = 1.0e+6;                  
rng(63761);
 
  % -- File Directories  
  outdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/out/';
  figdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/fig/';
  matdir = '/Users/mwatson/Dropbox/p_comp/revision/matlab/mat/';
  
  
   % -- Read in Data --- 
  load_data = 1;  % 1 if reloading data from Excel, etc 
  mtoq_agg = 3;   % Temporal aggregation indicator of monthly to quarterly data
  pcomp_data_calendar_m_and_q;
  % Data Series Use
  share_avg = share_avg_q;
  calvec = calvec_q;
  dnobs = dnobs_q;
  calds = calds_q;
  nper = 4;
  label_suffix = '_Q3';
  
  
  % Table 1: Shares
  first_date = [1960 1];
  last_date = calds(end,:);
  ismpl1 = smpl(calvec,first_date,last_date,nper);
  ismpl2 = smpl(calvec,first_date,[1979 nper],nper);
  ismpl3 = smpl(calvec,[1980 1],[1999 nper],nper);
  ismpl4 = smpl(calvec,[2000 1],last_date,nper);
  ismpl = [ismpl1 ismpl2 ismpl3 ismpl4];
  tmp_mat = NaN*zeros(17,4);
  for i = 1:4;
      tmp = share_avg(ismpl(:,i)==1,:);
      tmp_mat(:,i) = (mean(tmp))';
  end;
  outfile_name = [outdir 'Table1_Shares' label_suffix '.out'];
  fileID = fopen(outfile_name,'w');
  for i = 1:17;
      tmpstr = char(labelvec_disagg(i));
      fprintf(fileID,[tmpstr ',']);
      prtmat_comma(tmp_mat(i,:),fileID,'%5.3f','\n');
  end;
  type(outfile_name);
 