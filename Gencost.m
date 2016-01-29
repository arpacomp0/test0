function Cost = Gencost(xx,alph,beta,gamma)
% Quadratic Generation Cost Function
% Cost is in units of $/hour
% Uses only first three components of argument xx, which
% are the per unit values of active power setpoint for each generator
%
PG_MW=1000*xx(1:3); % convert from pu to MWs, using 1000 MVA base
%
Cost=sum(alph+beta.*PG_MW+gamma.*PG_MW.*PG_MW);
%
end

