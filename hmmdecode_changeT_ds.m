function [pStates,pSeq, fs, bs, s] = hmmdecode_changeT_ds(seq,T,e,changeT)
%   1st pos can be other than 1
%   See also HMMGENERATE, HMMESTIMATE, HMMVITERBI, HMMTRAIN.

%   Reference: Biological Sequence Analysis, Durbin, Eddy, Krogh, and
%   Mitchison, Cambridge University Press, 1998.  

%   Copyright 1993-2011 The MathWorks, Inc.


if nargin <= 3
    
    changeT = true;
    
end

if ~changeT
    
    assert(size(T,2)==1);
    
end

numStates = size(e,1);

numSymbols = size(e,2);


% add extra symbols to start to make algorithm cleaner at f0 and b0
seq = [numSymbols+1, seq ];
L = length(seq);

% This is what we'd like to do but it is numerically unstable
% warnState = warning('off');
% logTR = log(tr);
% logE = log(e);
% warning(warnState);
% f = zeros(numStates,L);
% f(1,1) = 1;
% % for count = 2:L
%     for state = 1:numStates
%         f(state,count) = logE(state,seq(count)) + log(sum( exp(f(:,count-1) + logTR(:,state))));
%     end
% end
% f = exp(f);

% so we introduce a scaling factor
fs = zeros(numStates,L);
fs(1,1) = 1;  % assume that we start in state 1.
s = zeros(1,L);
s(1) = 1;
len = (numStates-1)/2;
%first treat as the second, from the null states

tr = wrap_T(T(:,1),len);
for count = 2:L
    %tr = wrap_T(T(count-1),len);
    if changeT
        tr(1,1) = 1-sum(T(:,count-1));
        tr(1,2) = T(1,count-1);
        tr(1,2+len) = T(2,count-1);
    end
    for state = 1:numStates
        fs(state,count) = e(state,seq(count)) .* (sum(fs(:,count-1) .*tr(:,state)));
    end
    % scale factor normalizes sum(fs,count) to be 1. 
    s(count) =  sum(fs(:,count));
    fs(:,count) =  fs(:,count)./s(count);
end
bs = ones(numStates,L);
for count = L-1:-1:1
    if changeT
        tr(1,1) = 1-sum(T(:,count));
        tr(1,2) = T(1,count);
        tr(1,2+len) = T(2,count);
    end
    for state = 1:numStates
        bs(state,count) = (1/s(count+1)) * sum( tr(state,:)'.* bs(:,count+1) .* e(:,seq(count+1))); 
    end
end



    




pSeq = sum(log(s));
pStates = fs.*bs;

% get rid of the column that we stuck in to deal with the f0 and b0 
pStates(:,1) = [];

end



