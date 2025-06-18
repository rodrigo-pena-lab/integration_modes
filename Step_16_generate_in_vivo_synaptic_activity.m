
clear
clc

%Add synthetic minis
%%Acquisition rate 20 KHz. Simulations with time intervals = 0.05 ms
acq_rate = 20000; %Hz
delta_t = 1000/acq_rate; %ms

%Set Holding voltage
VHold = -30; %mV

%Generate synaptic conductance (AMPA)
t = delta_t:delta_t:200; %ms
gpeak = 0.1; %nS
tau_rise = 1; %ms
tau_decay = 10; %ms
gsyn_AMPA = biexp(t,gpeak,tau_rise,tau_decay);

%Generate AMPA EPSC current
Erev_AMPA = 0; %mV
Isyn_AMPA = gsyn_AMPA*(VHold-Erev_AMPA); %pA

%Save all data in structure for sparse activity
AMPA_param.t = t;
AMPA_param.gpeak = gpeak;
AMPA_param.tau_rise = tau_rise;
AMPA_param.tau_decay = tau_decay;
AMPA_param.VHold = VHold;
AMPA_param.Erev = Erev_AMPA;

%For synchrony activity
AMPA_param_2 = AMPA_param;
AMPA_param_2.gpeak = 4*gpeak;

%Generate synaptic conductance (GABAA)
gpeak = 0.1; %nS
tau_rise = 2; %ms
tau_decay = 20; %ms
gsyn_GABAA = biexp(t,gpeak,tau_rise,tau_decay); %nS

%Generate GABAA IPSC current
Erev_GABAA = -60; %mV
Isyn_GABAA = gsyn_GABAA*(VHold-Erev_GABAA); %pA

%Save all data in structure
GABA_param.t = t;
GABA_param.gpeak = gpeak;
GABA_param.tau_rise = tau_rise;
GABA_param.tau_decay = tau_decay;
GABA_param.VHold = VHold;
GABA_param.Erev = Erev_GABAA;

%Generate trace with several EPSCs distributed evenly randomly
%times
num_PSCs = 200; %Number of generated PSCs
length_all_sim = 1; %Total time of simulation in s
Tsim = delta_t:delta_t:length_all_sim*1000; %Total simulation time in ms

cell_traces = cell(1,100);
cell_traces_2 = cell(1,100);
for p = 1:100

    %mEPSCs only
    times_PSCs = length_all_sim*1000*rand(1,num_PSCs); %Generate evenly random times where PSCs are
    times_PSCs_pos = times_PSCs*acq_rate/1000; %Vector element where PSC onset
    times_PSCs_pos = ceil(times_PSCs_pos);

    %Sparse activity
    EPSCs_pos = sort(times_PSCs_pos);
    All_trace = zeros(1,length(Tsim)); %Total simulation time in ms
    trace_EPSCs = rand_PSCs(EPSCs_pos,All_trace,AMPA_param);

    %Synchrony activity
    even_pos = 1:4:num_PSCs;
    EPSC_pos_2 = EPSCs_pos(even_pos);
    trace_EPSCs_2 = rand_PSCs(EPSC_pos_2,All_trace,AMPA_param_2);


    %mIPSCs only
    times_PSCs = length_all_sim*1000*rand(1,num_PSCs); %Generate evenly random times where PSCs are
    times_PSCs_pos = times_PSCs*acq_rate/1000; %Vector element where PSC onset
    times_PSCs_pos = ceil(times_PSCs_pos);
    IPSCs_pos = sort(times_PSCs_pos);

    All_trace = zeros(1,length(Tsim)); %Total simulation time in ms
    trace_IPSCs = rand_PSCs(times_PSCs_pos,All_trace,GABA_param);

    %Both mEPSCs and mIPSCs in same trace
    all_long_biphasic_traces = -(trace_EPSCs + trace_IPSCs);
    all_long_biphasic_traces_2 = -(trace_EPSCs_2 + trace_IPSCs);

    cell_traces{p} = all_long_biphasic_traces;
    cell_traces_2{p} = all_long_biphasic_traces_2;

end

save sparse_synchro_act.mat cell_traces cell_traces_2

figure()
plot(Tsim,all_long_biphasic_traces,'k','LineWidth',2)
box off
hold on
plot(Tsim,all_long_biphasic_traces_2,'r','LineWidth',2)
box off
legend({'Sparse' 'Synchronous'})
legend boxoff
ylabel('Current (a.u)')
xlabel('ms')
set(gca, 'FontSize', 20)

%Resize plot size
set(gcf,'units','inches','position',[0 0 16 4])

%Save vector figure in pdf
exportgraphics(gcf,'Step_4A_invivo.pdf','ContentType','vector');

