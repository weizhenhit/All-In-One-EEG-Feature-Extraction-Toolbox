function [A,f,tt] = func_hhspectrum(imf,t,l,aff)

% [A,f,tt] = HHSPECTRUM(imf,t,l,aff) computes the Hilbert-Huang spectrum
%
% inputs:
% 	- imf : matrix with one IMF per row
%   - t   : time instants
%   - l   : estimation parameter for instfreq
%   - aff : if 1, displays the computation evolution
%
% outputs:
%   - A   : amplitudes of IMF rows
%   - f   : instantaneous frequencies
%   - tt  : truncated time instants
%
% calls:
%   - hilbert  : computes the analytic signal
%   - instfreq : computes the instantaneous frequency

if nargin < 2

  t=1:size(imf,2);

end

if nargin < 3

  l=1;

end

if nargin < 4

  aff = 0;

end

lt=length(t);

tt=t((l+1):(lt-l));

for i=1:(size(imf,1)-1)

  an(i,:)=hilbert(imf(i,:)')';      % 调用hilbert函数得到第i个IMF分量的解析表示an(i,:)
  f(i,:)=func_instfreq(an(i,:)',tt,l)';  % 基于解析表示计算imfi的瞬时频率
  A=abs(an(:,l+1:end-l));   
  % disp(size(A));
  if aff

    disp(['mode ',int2str(i),' trait']);

  end

end
