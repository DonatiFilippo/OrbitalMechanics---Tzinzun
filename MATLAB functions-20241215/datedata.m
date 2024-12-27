function [w_dep_datenum, w_fb_datenum, w_arr_datenum, t_opt_datenum] = datedata(dat1, dat2, dat3, t_opt_sol)

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