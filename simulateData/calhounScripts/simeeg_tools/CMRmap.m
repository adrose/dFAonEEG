function map = CMRmap(varargin)
%% CMRmap is a Colormap with Mapping to Grayscale for Publications.
%% Obtained from http://www.mathworks.com/matlabcentral/fileexchange/2662

%% Call this just as you call any other colormap.

if isempty(varargin);varargin={128};end;

tempmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];

x = 1:8/(varargin{:}-1):9;											% 64 color levels instead of 9
x1 = 1:9;
%x1 = logspace(0,log10(9),9);
for i = 1:3
	map(:,i) = spline(x1,tempmap(:,i),x)';	% spline fit intermediate values
end
map = abs(map/max(max(map)));		% eliminate spurious values outside of range 

