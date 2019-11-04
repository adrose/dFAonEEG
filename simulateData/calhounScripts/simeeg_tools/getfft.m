%
% [data freq] = getfft(data,rate,freq,dim,N)
%
% GETFFT returns the fft and the frequency array of 'data.'  The Goertzel
% algorithm is used if the number of frequencies is small enough.
%
% data: data array to be ffted
% rate: sampling rate of data, defaults to 1 Hz
% freq: maximum frequency to keep, defaults to nyquist limit for real data
%       alternately, an vector of frequencies to keep
%  dim: dimension along which to FFT the data, defaults to first nonzero
%    N: length of fft to take, defaults to size(data,dim)
%
% data: ffted data
% freq: list of frequencies
%
% See also: FFT, GOERTZEL, FREQPARSE
%
% Written by Bill Winter, March 2006, modified March 2008 
function [x freq] = getfft(x,rate,freq,dim,N)
siz = size(x);
if nargin < 2 || isempty(rate), rate = 1; end
if nargin < 3 || isempty(freq), freq = inf; end
if nargin < 4 || isempty(dim), dim = find(siz > 1,1); end
if nargin < 5 || isempty(N), N = siz(dim); end
[freq ind] = freqparse(rate,freq,N,x);
if length(ind) < log2(N), x = goertzel(trunkpad(x,N,dim),ind,dim)/N;
else                                                % full fft
    x = fft(x,N,dim);                               % fft data
    in(1:length(siz)) = {':'};in{dim} = ind;
    x = x(in{:})/N;                                 % wanted frequencies
end
% x = x*rate/N;