function [J0, UT] = J0_computation(jd, date)
%
% J0_computation - Compute the Julian Day at 0 UT, and the Universal Time
% of a given Date.
%
% PROTOTYPE:
%   [J0, UT] = J0_computation(jd, date)
%
% DESCRIPTION:
%   The function calculate the Julian Day number at 0 UT and the Universal
%   Time of a given Date, if provided with the Julian Day Number of the
%   same Date
%
% INPUT:
%   jd [1x1]        Date in Julian Day. The JD (Julian day) count is from 0 
%                   at 12:00 noon, 1 January -4712 (4713 BC), 
%                   Julian proleptic calendar. The corresponding date 
%                   in Gregorian calendar is 12:00 noon, 24 November -4713.
%
%   date [1x6]      Date in the Gregorian calendar, as a 6-elements vector
%                   [year, month, day, hour, minute, second]. For dates 
%                   before 1582, the resulting date components are valid 
%                   only in the Gregorian proleptic calendar. This is based
%                   on the Gregorian calendar but extended to cover dates 
%                   before its introduction. Date must be after 12:00 noon,
%                   24 November -4713.
%
% OUTPUT:
%   J0 [1x1]        Julian Date Number at 0 UT, of a given Date
%
%   UT [1x1]        Universal Time of a given Date                     [hr]
%
% CONTRIBUTORS:
%   Azevedo Da Silva Esteban
%   Gavidia Pantoja Maria Paulina
%   Donati Filippo 
%   Domenichelli Eleonora
%
%-------------------------------------------------------------------------

% Variable Extraction
hours = date(4);
minutes = date(5);
seconds = date(6);

% Universal Time Computation
UT = hours + minutes/60 + seconds/3600;

% Julian Day Number at 0 UT Computation
J0 = jd - UT/24;
end