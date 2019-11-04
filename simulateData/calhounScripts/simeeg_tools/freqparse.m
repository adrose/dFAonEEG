%
% [freq ind] = freqparse(rate,freq,N,flag)
%
% FREQPARSE processes a frequency vector.  Frequencies are modulated to the
% range [-rate/2 rate/2].
% 
% If N is a positive integer, frequencies are fit to discrete values
%    If freq is a scalar, it is treated as the maximum desired frequency.
%       If flag is real, only positive frequencies are retained.
%
% If N is empty, frequencies are treated as continuous
%
% inputs -
% rate: sampling rate of the data
% freq: desired frequencies; defaults to all
%       scalar: maximum frequency
%       vector: list of frequencies
%    N: length of fft; default is inf for continuous
% flag: 1 for real data, i for imaginary data; default is 1
%
% outputs -
% freq: processed frequencies
%  ind: index of fft for discrete frequency fourier transforms

% Written by Bill Winter, June 2008
function [freq ind] = freqparse(rate,freq,N,x)
if nargin < 2 || isempty(freq), freq = inf; end     % default: all freqs
if nargin < 3 || any(N ~= floor(N)), N = []; end    % default: not discrete
if nargin < 4, x = 1; end                           % default: real data
x = ~any(imag(x(:)));
freq = freq(:)';
if isempty(N) || ~isfinite(N) || N <= 0             % arbitrary frequencies
    freq = unique(mod(freq(:)',rate));              % unique frequencies
    ind = nan(size(freq));
else                                                % discrete frequencies
    ind = 1 + unique(mod(round(freq(:)'*N/rate),N));% unique frequencies
    if isscalar(freq)                               % make frequency vector
        ind(ind == 1) = N;
        freq = min(ceil((N+1)/2),ind);              % compare to nyquist
        ind = 1:N;
        if x, ind(ind > freq) = [];                 % if data is complex,
        else  ind((ind>freq)&(ind<N-freq+2)) = [];  % keep negative freq's
        end
    end
    freq = (ind - 1)*rate/N;
end
freq(freq > rate/2) = freq(freq > rate/2) - rate;   % negative frequencies