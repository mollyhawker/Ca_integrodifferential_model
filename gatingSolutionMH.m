function [gate] =gatingSolutionMH(rate,alpha,tdiff,G0)
%% Brady (1972) model.
% Integrals have been simplified to use a Riemann sum approximation
cumulative_rate=cumsum(rate(1:end-1).*tdiff);

% Single integral in the solution of gating variable model (Brady, 1972)
cumulative_exp=exp(cumulative_rate);

% Double integral in the solution of gating variable model (Brady, 1972)
integral_2=sum(cumulative_exp.*-alpha(1:end-1).*tdiff);

% full solution
gate=(cumulative_exp(end)^(-1))*(G0 - integral_2);

end