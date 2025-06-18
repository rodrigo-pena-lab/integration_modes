
function trace_PSCs = rand_PSCs(times_PSCs_pos,All_trace,Isyn_param)

%Generate Isyn with variability in kinetics (tau rise and decay)
tau_rise = Isyn_param.tau_rise;
tau_decay = Isyn_param.tau_decay;

for i = 1:length(times_PSCs_pos)
    
    %Kinetic variability
    var_rise = tau_rise;
    var_decay = tau_decay;
    
    gsyn = biexp(Isyn_param.t,Isyn_param.gpeak,var_rise,var_decay);
    Isyn = gsyn*(Isyn_param.VHold-Isyn_param.Erev); %pA    
    yy = times_PSCs_pos(i) + length(Isyn)-1;
    
    if yy <= length(All_trace)
        xx = All_trace(times_PSCs_pos(i):yy);
        All_trace(times_PSCs_pos(i):yy) = xx + Isyn;
    else
        xx = All_trace(times_PSCs_pos(i):end);
        All_trace(times_PSCs_pos(i):end) = xx + Isyn(1:length(xx));
    end
end

trace_PSCs = All_trace;

end