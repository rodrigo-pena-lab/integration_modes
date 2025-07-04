
clear
clc

Tsim = 400; %ms
VHold = -50;

%%Acquisition rate 20 KHz. Simulations with time intervals = 0.05 ms
acq_rate = 20000; %Hz
delta_t = 1000/acq_rate; %ms

%Generate synaptic conductance (AMPA)
t = delta_t:delta_t:200; %ms
gpeak = 0.04; % nS
tau_rise = 1; %ms
tau_decay = 10; %ms
gsyn_AMPA = biexp(t,gpeak,tau_rise,tau_decay);

%Generate AMPA EPSC current
Erev_AMPA = 0; %mV
Isyn_AMPA = gsyn_AMPA*(VHold-Erev_AMPA); %pA
Isyn_AMPA = Isyn_AMPA(1:1000); %50 ms trace
trace_EPSCs = -Isyn_AMPA;

%Multiple 3 EPSPs
mult_EPSP_trace = zeros(1,3*length(trace_EPSCs));
mult_EPSP_trace(1:length(trace_EPSCs)) = trace_EPSCs;
mult_EPSP_trace(201:length(trace_EPSCs)+200) = trace_EPSCs + mult_EPSP_trace(201:length(trace_EPSCs)+200) ; %10 ms after
mult_EPSP_trace(401:length(trace_EPSCs)+400) = trace_EPSCs + mult_EPSP_trace(401:length(trace_EPSCs)+400); %20 ms after

%%
% Specific capacitance = 1 microF/cm2
Cm = 1;
dt = 0.05;               % time step for forward euler method (ms)
loop  = ceil(Tsim/dt);   % no. of iterations of euler


%% Condition 1
%Conductance mS/cm2
gL = 0.3;
gak = 0;
gcaT = 0;
gnaT = 0;

%Reversal potential mV
eNa = 50;
eK =-77;
eCa = 50;
eL = VHold; %Leak reversal potential

%Injected current (EPSP train)
ton = 200; %ms
Iinj_ini = zeros(loop,1);
pos_on = fix(ton/dt);
Iinj_ini(pos_on+1:pos_on+length(mult_EPSP_trace)) = mult_EPSP_trace; %For synaptic current
Iinj = Iinj_ini;

% Initializing variable vectors
t = (1:loop)*dt;
V = zeros(loop,1);
nka = zeros(loop,1); %A-type K
lka = zeros(loop,1);
mcaT = zeros(loop,1); %T-type Ca
hcaT = zeros(loop,1);
mna = zeros(loop,1); %Transient sodium
hna = zeros(loop,1);

% Set initial values for the variables
V(1)=VHold;
nka(1) = 0; %A-type K
lka(1) = 1;
mcaT(1) = 0; %T-type Ca
hcaT(1) = 1;
mna(1) = 0; %Transient sodium
hna(1) = 1;

%Time constants
taunka = 0.2; %ms %A-type K
taulka = 5;
taumcaT = 0.5; %ms %T-type Ca
tauhcaT = 10;
mtauna = 0.05; %Transient sodium
htauna = 0.5;

%Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    %A-type potassium act/inact
    nka(i+1) = nka(i) +dt*((ninf_ka(V(i)) - nka(i))/taunka);
    lka(i+1) = lka(i) + dt*((linf_ka(V(i)) - lka(i))/taulka);

    %T-type calcium act/inact
    mcaT(i+1) = mcaT(i) + dt*((minf_caT(V(i)) - mcaT(i))/taumcaT);
    hcaT(i+1) = hcaT(i) + dt*((hinf_caT(V(i)) - hcaT(i))/tauhcaT);

    %Transient sodium act/inact
    mna(i+1) = mna(i) + dt*((minf_na(V(i)) - mna(i))/mtauna);
    hna(i+1) = hna(i) + dt*((hinf_na(V(i)) - hna(i))/htauna);

end
V1 = V;


%% Condition 2
%Conductance mS/cm2
gak = 10;

%Injected current (EPSP train)
Iinj = Iinj_ini + 5.7;

% Initializing variable vectors
V = zeros(loop,1);
nka = zeros(loop,1); %A-type K
lka = zeros(loop,1);
mcaT = zeros(loop,1); %T-type Ca
hcaT = zeros(loop,1);
mna = zeros(loop,1); %Transient sodium
hna = zeros(loop,1);

% Set initial values for the variables
V(1) = VHold;
nka(1) = 0; %A-type K
lka(1) = 1;
mcaT(1) = 0; %T-type Ca
hcaT(1) = 1;
mna(1) = 0; %Transient sodium
hna(1) = 1;

%Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    %A-type potassium act/inact
    nka(i+1) = nka(i) +dt*((ninf_ka(V(i)) - nka(i))/taunka);
    lka(i+1) = lka(i) + dt*((linf_ka(V(i)) - lka(i))/taulka);

    %T-type calcium act/inact
    mcaT(i+1) = mcaT(i) + dt*((minf_caT(V(i)) - mcaT(i))/taumcaT);
    hcaT(i+1) = hcaT(i) + dt*((hinf_caT(V(i)) - hcaT(i))/tauhcaT);

    %Transient sodium act/inact
    mna(i+1) = mna(i) + dt*((minf_na(V(i)) - mna(i))/mtauna);
    hna(i+1) = hna(i) + dt*((hinf_na(V(i)) - hna(i))/htauna);

end
V2 = V;



%% Condition 3
%Injected current (EPSP train)
Iinj = Iinj_ini + 16.6;

% Initializing variable vectors
V = zeros(loop,1);
nka = zeros(loop,1); %A-type K
lka = zeros(loop,1);
mcaT = zeros(loop,1); %T-type Ca
hcaT = zeros(loop,1);
mna = zeros(loop,1); %Transient sodium
hna = zeros(loop,1);

% Set initial values for the variables
V(1) = VHold;
nka(1) = 0; %A-type K
lka(1) = 1;
mcaT(1) = 0; %T-type Ca
hcaT(1) = 1;
mna(1) = 0; %Transient sodium
hna(1) = 1;

%Time constants ms
taulka = 10000; %K inact

%Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    %A-type potassium act/inact
    nka(i+1) = nka(i) +dt*((ninf_ka(V(i)) - nka(i))/taunka);
    lka(i+1) = lka(i) + dt*((linf_ka(V(i)) - lka(i))/taulka);

    %T-type calcium act/inact
    mcaT(i+1) = mcaT(i) + dt*((minf_caT(V(i)) - mcaT(i))/taumcaT);
    hcaT(i+1) = hcaT(i) + dt*((hinf_caT(V(i)) - hcaT(i))/tauhcaT);

    %Transient sodium act/inact
    mna(i+1) = mna(i) + dt*((minf_na(V(i)) - mna(i))/mtauna);
    hna(i+1) = hna(i) + dt*((hinf_na(V(i)) - hna(i))/htauna);

end
V3 = V;


figure

%Subtract baseline
V1 = V1 - mean(V1(fix(150/dt):fix(155/dt)));
V2 = V2 - mean(V2(fix(150/dt):fix(155/dt)));
V3 = V3 - mean(V3(fix(150/dt):fix(155/dt)));

subplot(2,3,1)
plot(t,V1,'k',t,V2,'b',t,V3,'r','LineWidth',2);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
xlim([190 280])
ylim([0 15])
box off
legend({'Passive' 'Inactivating' 'Non-inactivating'})
legend boxoff
title('Potassium')
box off
set(gca, 'FontSize', 20)

%Normalize
V1 = V1/max(V1(pos_on:pos_on+200));
V2 = V2/max(V2(pos_on:pos_on+200));
V3 = V3/max(V3(pos_on:pos_on+200));

subplot(2,3,4)
plot(t,V1,'k',t,V2,'b',t,V3,'r','LineWidth',2);
xlabel('Time (ms)');
ylabel('Normalized');
xlim([190 280])
ylim([-0.1 2.2])
box off
set(gca, 'FontSize', 20)





%% Sodium
%% Condition 1
%Conductance mS/cm2
gL = 0.3;
gak = 0;
gcaT = 0;
gnaT = 0;

% Initializing variable vectors
t = (1:loop)*dt;
V = zeros(loop,1);
nka = zeros(loop,1); %A-type K
lka = zeros(loop,1);
mcaT = zeros(loop,1); %T-type Ca
hcaT = zeros(loop,1);
mna = zeros(loop,1); %Transient sodium
hna = zeros(loop,1);

% Set initial values for the variables
V(1)=VHold;
nka(1) = 0; %A-type K
lka(1) = 1;
mcaT(1) = 0; %T-type Ca
hcaT(1) = 1;
mna(1) = 0; %Transient sodium
hna(1) = 1;

%Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    %A-type potassium act/inact
    nka(i+1) = nka(i) +dt*((ninf_ka(V(i)) - nka(i))/taunka);
    lka(i+1) = lka(i) + dt*((linf_ka(V(i)) - lka(i))/taulka);

    %T-type calcium act/inact
    mcaT(i+1) = mcaT(i) + dt*((minf_caT(V(i)) - mcaT(i))/taumcaT);
    hcaT(i+1) = hcaT(i) + dt*((hinf_caT(V(i)) - hcaT(i))/tauhcaT);

    %Transient sodium act/inact
    mna(i+1) = mna(i) + dt*((minf_na(V(i)) - mna(i))/mtauna);
    hna(i+1) = hna(i) + dt*((hinf_na(V(i)) - hna(i))/htauna);

end
V1 = V;


%% Condition 2
%Conductance mS/cm2
gnaT = 0.5;

%Injected current (EPSP train)
Iinj = Iinj_ini;

% Initializing variable vectors
V = zeros(loop,1);
nka = zeros(loop,1); %A-type K
lka = zeros(loop,1);
mcaT = zeros(loop,1); %T-type Ca
hcaT = zeros(loop,1);
mna = zeros(loop,1); %Transient sodium
hna = zeros(loop,1);

% Set initial values for the variables
V(1) = VHold;
nka(1) = 0; %A-type K
lka(1) = 1;
mcaT(1) = 0; %T-type Ca
hcaT(1) = 1;
mna(1) = 0; %Transient sodium
hna(1) = 1;

%Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    %A-type potassium act/inact
    nka(i+1) = nka(i) +dt*((ninf_ka(V(i)) - nka(i))/taunka);
    lka(i+1) = lka(i) + dt*((linf_ka(V(i)) - lka(i))/taulka);

    %T-type calcium act/inact
    mcaT(i+1) = mcaT(i) + dt*((minf_caT(V(i)) - mcaT(i))/taumcaT);
    hcaT(i+1) = hcaT(i) + dt*((hinf_caT(V(i)) - hcaT(i))/tauhcaT);

    %Transient sodium act/inact
    mna(i+1) = mna(i) + dt*((minf_na(V(i)) - mna(i))/mtauna);
    hna(i+1) = hna(i) + dt*((hinf_na(V(i)) - hna(i))/htauna);

end
V2 = V;



%% Condition 3
%Injected current (EPSP train)
Iinj = Iinj_ini;

% Initializing variable vectors
V = zeros(loop,1);
nka = zeros(loop,1); %A-type K
lka = zeros(loop,1);
mcaT = zeros(loop,1); %T-type Ca
hcaT = zeros(loop,1);
mna = zeros(loop,1); %Transient sodium
hna = zeros(loop,1);

% Set initial values for the variables
V(1) = VHold;
nka(1) = 0; %A-type K
lka(1) = 1;
mcaT(1) = 0; %T-type Ca
hcaT(1) = 1;
mna(1) = 0; %Transient sodium
hna(1) = 1;

%Time constants ms
htauna = 10000; %Na inact

%Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    %A-type potassium act/inact
    nka(i+1) = nka(i) +dt*((ninf_ka(V(i)) - nka(i))/taunka);
    lka(i+1) = lka(i) + dt*((linf_ka(V(i)) - lka(i))/taulka);

    %T-type calcium act/inact
    mcaT(i+1) = mcaT(i) + dt*((minf_caT(V(i)) - mcaT(i))/taumcaT);
    hcaT(i+1) = hcaT(i) + dt*((hinf_caT(V(i)) - hcaT(i))/tauhcaT);

    %Transient sodium act/inact
    mna(i+1) = mna(i) + dt*((minf_na(V(i)) - mna(i))/mtauna);
    hna(i+1) = hna(i) + dt*((hinf_na(V(i)) - hna(i))/htauna);

end
V3 = V;



%Subtract baseline
V1 = V1 - mean(V1(fix(150/dt):fix(155/dt)));
V2 = V2 - mean(V2(fix(150/dt):fix(155/dt)));
V3 = V3 - mean(V3(fix(150/dt):fix(155/dt)));

subplot(2,3,2)
plot(t,V1,'k',t,V2,'b',t,V3,'r','LineWidth',2);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
xlim([190 280])
ylim([0 15])
box off
title('Sodium')
set(gca, 'FontSize', 20)

%Normalize
V1 = V1/max(V1(pos_on:pos_on+200));
V2 = V2/max(V2(pos_on:pos_on+200));
V3 = V3/max(V3(pos_on:pos_on+200));

subplot(2,3,5)
plot(t,V1,'k',t,V2,'b',t,V3,'r','LineWidth',2);
xlabel('Time (ms)');
ylabel('Normalized');
xlim([190 280])
ylim([-0.1 2.2])
box off
set(gca, 'FontSize', 20)





%% Calcium
%% Condition 1
%Conductance mS/cm2
gL = 0.3;
gak = 0;
gcaT = 0;
gnaT = 0;

% Initializing variable vectors
t = (1:loop)*dt;
V = zeros(loop,1);
nka = zeros(loop,1); %A-type K
lka = zeros(loop,1);
mcaT = zeros(loop,1); %T-type Ca
hcaT = zeros(loop,1);
mna = zeros(loop,1); %Transient sodium
hna = zeros(loop,1);

% Set initial values for the variables
V(1)=VHold;
nka(1) = 0; %A-type K
lka(1) = 1;
mcaT(1) = 0; %T-type Ca
hcaT(1) = 1;
mna(1) = 0; %Transient sodium
hna(1) = 1;

%Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    %A-type potassium act/inact
    nka(i+1) = nka(i) +dt*((ninf_ka(V(i)) - nka(i))/taunka);
    lka(i+1) = lka(i) + dt*((linf_ka(V(i)) - lka(i))/taulka);

    %T-type calcium act/inact
    mcaT(i+1) = mcaT(i) + dt*((minf_caT(V(i)) - mcaT(i))/taumcaT);
    hcaT(i+1) = hcaT(i) + dt*((hinf_caT(V(i)) - hcaT(i))/tauhcaT);

    %Transient sodium act/inact
    mna(i+1) = mna(i) + dt*((minf_na(V(i)) - mna(i))/mtauna);
    hna(i+1) = hna(i) + dt*((hinf_na(V(i)) - hna(i))/htauna);

end
V1 = V;


%% Condition 2
%Conductance mS/cm2
gcaT = 0.1;

%Injected current (EPSP train)
Iinj = Iinj_ini - 0.85;

% Initializing variable vectors
V = zeros(loop,1);
nka = zeros(loop,1); %A-type K
lka = zeros(loop,1);
mcaT = zeros(loop,1); %T-type Ca
hcaT = zeros(loop,1);
mna = zeros(loop,1); %Transient sodium variables
hna = zeros(loop,1);

% Set initial values for the variables
V(1) = VHold;
nka(1) = 0; %A-type K
lka(1) = 1;
mcaT(1) = 0; %T-type Ca
hcaT(1) = 1;
mna(1) = 0; %Transient sodium variables
hna(1) = 1;

%Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    %A-type potassium act/inact
    nka(i+1) = nka(i) +dt*((ninf_ka(V(i)) - nka(i))/taunka);
    lka(i+1) = lka(i) + dt*((linf_ka(V(i)) - lka(i))/taulka);

    %T-type calcium act/inact
    mcaT(i+1) = mcaT(i) + dt*((minf_caT(V(i)) - mcaT(i))/taumcaT);
    hcaT(i+1) = hcaT(i) + dt*((hinf_caT(V(i)) - hcaT(i))/tauhcaT);

    %Transient sodium act/inact
    mna(i+1) = mna(i) + dt*((minf_na(V(i)) - mna(i))/mtauna);
    hna(i+1) = hna(i) + dt*((hinf_na(V(i)) - hna(i))/htauna);

end
V2 = V;



%% Condition 3
%Injected current (EPSP train)
Iinj = Iinj_ini - 2.475;

% Initializing variable vectors
V = zeros(loop,1);
nka = zeros(loop,1); %A-type K
lka = zeros(loop,1);
mcaT = zeros(loop,1); %T-type Ca
hcaT = zeros(loop,1);
mna = zeros(loop,1); %Transient sodium
hna = zeros(loop,1);

% Set initial values for the variables
V(1) = VHold;
nka(1) = 0; %A-type K
lka(1) = 1;
mcaT(1) = 0; %T-type Ca
hcaT(1) = 1;
mna(1) = 0; %Transient sodium
hna(1) = 1;

%Time constants ms
tauhcaT = 10000; %Ca inact

%Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    %A-type potassium act/inact
    nka(i+1) = nka(i) +dt*((ninf_ka(V(i)) - nka(i))/taunka);
    lka(i+1) = lka(i) + dt*((linf_ka(V(i)) - lka(i))/taulka);

    %T-type calcium act/inact
    mcaT(i+1) = mcaT(i) + dt*((minf_caT(V(i)) - mcaT(i))/taumcaT);
    hcaT(i+1) = hcaT(i) + dt*((hinf_caT(V(i)) - hcaT(i))/tauhcaT);

    %Transient sodium act/inact
    mna(i+1) = mna(i) + dt*((minf_na(V(i)) - mna(i))/mtauna);
    hna(i+1) = hna(i) + dt*((hinf_na(V(i)) - hna(i))/htauna);

end
V3 = V;



%Subtract baseline
V1 = V1 - mean(V1(fix(150/dt):fix(155/dt)));
V2 = V2 - mean(V2(fix(150/dt):fix(155/dt)));
V3 = V3 - mean(V3(fix(150/dt):fix(155/dt)));

subplot(2,3,3)
plot(t,V1,'k',t,V2,'b',t,V3,'r','LineWidth',2);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
xlim([190 280])
ylim([0 15])
box off
title('Calcium')
set(gca, 'FontSize', 20)

%Normalize
V1 = V1/max(V1(pos_on:pos_on+200));
V2 = V2/max(V2(pos_on:pos_on+200));
V3 = V3/max(V3(pos_on:pos_on+200));

subplot(2,3,6)
plot(t,V1,'k',t,V2,'b',t,V3,'r','LineWidth',2);
xlabel('Time (ms)');
ylabel('Normalized');
xlim([190 280])
ylim([-0.1 2.2])
box off
set(gca, 'FontSize', 20)


%Resize plot size
set(gcf,'units','inches','position',[0 0 16 8])

%Save vector figure in pdf
exportgraphics(gcf,'Step_15.pdf','ContentType','vector');







%A-type K variables
function n1ka = ninf_ka(V)
n1ka = 1/(1 + exp(-(V+50)/12));
end
function l1ka = linf_ka(V)
l1ka = 1/(1 + exp((V+60)/12));
end

%T-type Ca variables
function m1ca = minf_caT(V)
m1ca = 1/(1 + exp(-(V+40)/10));
end
function h1ca = hinf_caT(V)
h1ca = 1/(1 + exp((V+60)/10));
end

%Transient sodium variables
function n1na = minf_na(V)
n1na = 1/(1 + exp(-(V+28)/5));
end
function l1na = hinf_na(V)
l1na = 1/(1 + exp((V+62)/5));
end

