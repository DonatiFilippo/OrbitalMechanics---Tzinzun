function [w_dep_datenum, w_fb_datenum, w_arr_datenum, t_opt_datenum] = datedata(dat1, dat2, dat3, t_opt_sol)
%% FUNCTION DESCRIPTION
%
%--------------------------------------------------------------------------
% DESCRIPTION:
% This function converts time windows and optimized solution dates from 
% Modified Julian Date 2000 (MJD2000) format to MATLAB's serial date number format (datenum). 
% This conversion facilitates plotting and date manipulations in MATLAB.
%
%--------------------------------------------------------------------------
% INPUTS:
%   dat1       [nx1]   Array of departure dates for the first leg [mjd2000]
%   dat2       [mx1]   Array of fly-by dates for the second leg   [mjd2000]
%   dat3       [px1]   Array of arrival dates for the third leg   [mjd2000]
%   t_opt_sol  [1x3]   Optimized solution dates [mjd2000] (departure, fly-by, arrival)
%
%--------------------------------------------------------------------------
% OUTPUTS:
%   w_dep_datenum  [nx1]   Serial date numbers for the departure window
%   w_fb_datenum   [mx1]   Serial date numbers for the fly-by window
%   w_arr_datenum  [px1]   Serial date numbers for the arrival window
%   t_opt_datenum  [1x3]   Serial date numbers for the optimized solution (departure, fly-by, arrival)
%
%--------------------------------------------------------------------------
% Group number : 27
%
% Created and maintained by : 
%   Azevedo Da Silva Esteban
%   Gavidia Pantoja Maria Paulina
%   Donati Filippo 
%   Domenichelli Eleonora
%
%--------------------------------------------------------------------------


% Prealocar las matrices de salida para eficiencia
w_dep_datenum = zeros(size(dat1));
w_fb_datenum = zeros(size(dat2));
w_arr_datenum = zeros(size(dat3));
t_opt_datenum = zeros(size(t_opt_sol));

for i = 1:length(dat1)
   w_dep_dates = mjd20002date(dat1(i));
   w_dep_datenum(i) = datenum(w_dep_dates);
end

for i = 1:length(dat2)
   w_fb_dates = mjd20002date(dat2(i));
   w_fb_datenum(i) = datenum(w_fb_dates);
end

for i = 1:length(dat3)
   w_arr_dates = mjd20002date(dat3(i));
   w_arr_datenum(i) = datenum(w_arr_dates);
end

for i = 1:length(t_opt_sol)
   t_opt_dates = mjd20002date(t_opt_sol(i));
   t_opt_datenum(i) = datenum(t_opt_dates);
end

end
