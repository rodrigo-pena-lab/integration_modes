
function gsyn = biexp(t,gpeak,tau_rise,tau_decay)

t_peak = ((tau_rise*tau_decay)/(tau_decay-tau_rise))*log(tau_decay/tau_rise); %ms

%Synapse conductance
%Biexponential(tau_rise and tau_decay)
a = exp(-t./tau_rise);
b = exp(-t./tau_decay);
c = exp(-t_peak/tau_rise);
d = exp(-t_peak/tau_decay);
gsyn = gpeak*(exp(a)-exp(b))/(exp(c)-exp(d));

end
