function [val1 val2] = simplevals(theval,numafter);
%% takes any input and returns the value to the left of the decimal place
%% as one variable and numafter values are the decimal place. These values
%% can then be appended to the end of a figure name, in order to read out
%% the minimum and maximum values in the figure etc. 

if nargin<2;numafter=2;end;
sc = 10^numafter;
theval = round(theval*sc)/sc;

val1 = floor(theval*sign(theval))*sign(theval);
val2 = abs(round((theval-floor(theval*sign(theval))*sign(theval))*sc));
