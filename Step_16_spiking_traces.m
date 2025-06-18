
clear
clc

load sparse_synchro_act.mat

all_rel_pass = [];
all_rel_k = [];
all_rel_Ca = [];
all_rel_Na = [];

%Select p from 1 to 100
p = 4;

%Load synaptic background signal (at 20 KHz)
all_long_biphasic_traces = cell_traces{p};
all_long_biphasic_traces_2 = cell_traces_2{p};

Iamp = 5; %microA/cm2
freq_in = 4; %Hz

Tsim = 1000; %ms
VHold = -80;

%Reversal potential
eNa = 50; %mV
eK =-80;
eCa = 50;
eL = VHold; %Leak reversal potential

%%
% Specific capacitance = 1 microF/cm2
Cm = 1;
dt = 0.05;               % time step for forward euler method (ms)
loop  = ceil(Tsim/dt);   % no. of iterations of euler
t = (1:loop)*dt;



%% Sparse activity


%% Condition 1
%Conductance mS/cm2
gL = 0.3;
gak = 0;
gcaT = 0;
gnaT = 0;

% Initializing variable vectors
V = zeros(loop,1);
nka = zeros(loop,1); %A-type K
lka = zeros(loop,1);
mcaT = zeros(loop,1); %T-type Ca
hcaT = zeros(loop,1);
mna = zeros(loop,1); %Transient sodium variables
hna = zeros(loop,1);

% Set initial values for the variables
V(1)=VHold;
nka(1) = 0; %A-type K
lka(1) = 1;
mcaT(1) = 0; %T-type Ca
hcaT(1) = 1;
mna(1) = 0; %Transient sodium variables
hna(1) = 1;

%Time constants
taunka = 0.2; %ms %A-type K
taulka = 5;
taumcaT = 0.5; %ms %T-type Ca
tauhcaT = 10;
mtauna = 0.05; %Transient sodium variables
htauna = 0.5;

%LIF parameters
Vth = -40; %mV
Erest = -60; %mV

tspk1 = zeros(loop,1);
%Injected current (constant)
Iinj = sin(2*pi*freq_in*t/1000);
Iinj = Iinj + 1;
Iinj = Iamp*Iinj + 0.5*all_long_biphasic_traces + 5;

% Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    if V(i+1) >= Vth     % elicit spikes
        V(i+1) = Erest;
        tspk1(i) = 1;
    end

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
gL = 0.3;
gak = 1;
gcaT = 0;
gnaT = 0;

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

% Euler method
tspk2 = zeros(loop,1);
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    if V(i+1) >= Vth     % elicit spikes
        V(i+1) = Erest;
        tspk2(i) = 1;
    end

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
%Conductance mS/cm2
gL = 0.3;
gak = 0;
gcaT = 0.1;
gnaT = 0;

% Initializing variable vectors
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

% Euler method
tspk3 = zeros(loop,1);
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    if V(i+1) >= Vth     % elicit spikes
        V(i+1) = Erest;
        tspk3(i) = 1;
    end

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




%% Condition 4
%Conductance mS/cm2
gL = 0.3;
gak = 0;
gcaT = 0; 
gnaT = 50; 

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

% Euler method
tspk4 = zeros(loop,1);
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    if V(i+1) >= Vth     % elicit spikes
        V(i+1) = Erest;
        tspk4(i) = 1;
    end

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
V4 = V;



V1_sparse = V1;
V2_sparse = V2;
V3_sparse = V3;
V4_sparse = V4;

tspk1_sparse = tspk1;
tspk2_sparse = tspk2;
tspk3_sparse = tspk3;
tspk4_sparse = tspk4;

%Plot spikes
figure
subplot(2,2,1)
plot(t,V1,'k');
spks = find(tspk1==1);
for k = 1:length(spks)
    line([t(spks(k)) t(spks(k))],[-50 0],'color','k')
end
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Trace w/o (Sparse)');
box off
freq_pass = length(spks);

subplot(2,2,2)
plot(t,V2,'b');
spks = find(tspk2==1);
for k = 1:length(spks)
    line([t(spks(k)) t(spks(k))],[-50 0],'color','b')
end
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('K');
box off
freq_K = length(spks);

subplot(2,2,3)
plot(t,V3,'b');
spks = find(tspk3==1);
for k = 1:length(spks)
    line([t(spks(k)) t(spks(k))],[-50 0],'color','b')
end
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Ca');
box off
freq_Ca = length(spks);

subplot(2,2,4)
plot(t,V4,'b');
spks = find(tspk4==1);
for k = 1:length(spks)
    line([t(spks(k)) t(spks(k))],[-50 0],'color','b')
end
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Na');
box off
freq_Na = length(spks);


%Resize plot size
set(gcf,'units','inches','position',[0 0 16 8])

%Save vector figure in pdf
exportgraphics(gcf,'Step_4A_spiking_sparse.pdf','ContentType','vector');





%% Synchronous

%% Condition 1
%Conductance mS/cm2
gL = 0.3;
gak = 0;
gcaT = 0;
gnaT = 0;

% Initializing variable vectors
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

tspk1 = zeros(loop,1);

%Injected current (constant)
Iinj = sin(2*pi*freq_in*t/1000);
Iinj = Iinj + 1;
Iinj = Iamp*Iinj + 0.5*all_long_biphasic_traces_2 + 5;

% Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    if V(i+1) >= Vth     % elicit spikes
        V(i+1) = Erest;
        tspk1(i) = 1;
    end

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
gL = 0.3;
gak = 1;
gcaT = 0;
gnaT = 0;

% Initializing variable vectors
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

% Euler method
tspk2 = zeros(loop,1);
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    if V(i+1) >= Vth     % elicit spikes
        V(i+1) = Erest;
        tspk2(i) = 1;
    end

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
%Conductance mS/cm2
gL = 0.3;
gak = 0;
gcaT = 0.1; 
gnaT = 0;

% Initializing variable vectors
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

% Euler method
tspk3 = zeros(loop,1);
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    if V(i+1) >= Vth     % elicit spikes
        V(i+1) = Erest;
        tspk3(i) = 1;
    end

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




%% Condition 4
%Conductance mS/cm2
gL = 0.3;
gak = 0;
gcaT = 0; 
gnaT = 50; 

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
V(1) = VHold;
nka(1) = 0; %A-type K
lka(1) = 1;
mcaT(1) = 0; %T-type Ca
hcaT(1) = 1;
mna(1) = 0; %Transient sodium
hna(1) = 1;

% Euler method
tspk4 = zeros(loop,1);
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i)); %Transient sodium

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + icaT + inaT + Iinj(i));

    if V(i+1) >= Vth     % elicit spikes
        V(i+1) = Erest;
        tspk4(i) = 1;
    end

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
V4 = V;


V1_sync = V1;
V2_sync = V2;
V3_sync = V3;
V4_sync = V4;

tspk1_sync = tspk1;
tspk2_sync = tspk2;
tspk3_sync = tspk3;
tspk4_sync = tspk4;



%Plot spikes
figure
subplot(2,2,1)
plot(t,V1,'k');
spks = find(tspk1==1);
for k = 1:length(spks)
    line([t(spks(k)) t(spks(k))],[-50 0],'color','k')
end
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Trace w/o (Synch)');
box off
freq_pass_2 = length(spks);

subplot(2,2,2)
plot(t,V2,'b');
spks = find(tspk2==1);
for k = 1:length(spks)
    line([t(spks(k)) t(spks(k))],[-50 0],'color','b')
end
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('K');
box off
freq_K_2 = length(spks);

subplot(2,2,3)
plot(t,V3,'b');
spks = find(tspk3==1);
for k = 1:length(spks)
    line([t(spks(k)) t(spks(k))],[-50 0],'color','b')
end
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Ca');
box off
freq_Ca_2 = length(spks);

subplot(2,2,4)
plot(t,V4,'b');
spks = find(tspk4==1);
for k = 1:length(spks)
    line([t(spks(k)) t(spks(k))],[-50 0],'color','b')
end
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Na');
box off
freq_Na_2 = length(spks);


%Frequency difference between sparse vs synchronous
diff_pass = freq_pass_2 - freq_pass;
diff_k = freq_K_2 - freq_K;
diff_Ca = freq_Ca_2 - freq_Ca;
diff_Na = freq_Na_2 - freq_Na;

%Relative change
rel_pass = diff_pass*100/freq_pass;
rel_k = diff_k*100/freq_K;
rel_Ca = diff_Ca*100/freq_Ca;
rel_Na = diff_Na*100/freq_Na;


%Resize plot size
set(gcf,'units','inches','position',[0 0 16 8])

%Save vector figure in pdf
exportgraphics(gcf,'Step_4A_spiking_sync.pdf','ContentType','vector');








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
