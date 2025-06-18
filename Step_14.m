
clear
clc

Tsim = 400; %ms
VHold = -50;

%%Acquisition rate 20 KHz. Simulations with time intervals = 0.05 ms
acq_rate = 20000; %Hz
delta_t = 1000/acq_rate; %ms

%Generate synaptic conductance (AMPA)
t = delta_t:delta_t:200; %ms
gpeak = 0.04; %nS
tau_rise = 1; %ms
tau_decay = 10; %ms
gsyn_AMPA = biexp(t,gpeak,tau_rise,tau_decay);

%Generate AMPA EPSC current
Erev_AMPA = 0; %mV
Isyn_AMPA = gsyn_AMPA*(VHold-Erev_AMPA); %pA
Isyn_AMPA = Isyn_AMPA(1:1000); %50 ms trace
trace_EPSCs = -Isyn_AMPA;


%%
% Specific capacitance = 1 microF/cm2
Cm = 1;
dt = 0.05;               % time step for forward euler method (ms)
loop  = ceil(Tsim/dt);   % no. of iterations of euler


%% Condition 1
%Conductance mS/cm2
gL = 0.3;
gak = 0;
gnaT = 0;

%Reversal potential
eNa = 50; %mV
eK =-77;
eCa = 50;
eL = VHold; %Leak reversal potential

%Injected current (constant)
%Current amplitude, duration of simulation, onset of current,
%offset of current
I = 1.5; %microA/cm2
ton = 200; %ms
toff  = 250; %ms
Iinj = zeros(loop,1);
pos_on = fix(ton/dt);
pos_off = fix(toff/dt);
Iinj(pos_on:pos_off) = I;

% Initializing variable vectors
V = zeros(loop,1);
nka = zeros(loop,1); %A-type K
lka = zeros(loop,1);
mna = zeros(loop,1); %Transient sodium
hna = zeros(loop,1);

% Set initial values for the variables
V(1) = VHold;
nka(1) = 0; %A-type K
lka(1) = 1;
mna(1) = 0; %Transient sodium
hna(1) = 1;

%Time constants ms
taunka = 0.2; %A-type K
taulka = 5;
mtauna = 0.05; %Transient sodium
htauna = 0.5;

%Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i));

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + inaT + Iinj(i));

    %A-type potassium act/inact
    nka(i+1) = nka(i) +dt*((ninf_ka(V(i)) - nka(i))/taunka);
    lka(i+1) = lka(i) + dt*((linf_ka(V(i)) - lka(i))/taulka);

    %Transient sodium act/inact
    mna(i+1) = mna(i) + dt*((minf_na(V(i)) - mna(i))/mtauna);
    hna(i+1) = hna(i) + dt*((hinf_na(V(i)) - hna(i))/htauna);

end
V1 = V;



%% Condition 2
%Conductance mS/cm2
gak = 10;
gnaT = 0;

% Initializing variable vectors
t = (1:loop)*dt;
V = zeros(loop,1);
nka = zeros(loop,1); %A-type K
lka = zeros(loop,1);
mna = zeros(loop,1); %Transient sodium variables
hna = zeros(loop,1);

% Set initial values for the variables
V(1) = VHold;
nka(1) = 0; %A-type K
lka(1) = 1;
mna(1) = 0; %Transient sodium
hna(1) = 1;

%Time constants ms
taunka = 0.2; %A-type K
taulka = 5;
mtauna = 0.1; %Transient sodium
htauna = 0.4;

%Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i));

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + inaT + Iinj(i));

    %A-type potassium act/inact
    nka(i+1) = nka(i) +dt*((ninf_ka(V(i)) - nka(i))/taunka);
    lka(i+1) = lka(i) + dt*((linf_ka(V(i)) - lka(i))/taulka);

    %Transient sodium act/inact
    mna(i+1) = mna(i) + dt*((minf_na(V(i)) - mna(i))/mtauna);
    hna(i+1) = hna(i) + dt*((hinf_na(V(i)) - hna(i))/htauna);

end
V2 = V;


figure
subplot(3,3,1)
plot(t,V1,'k',t,V2,'b','LineWidth',2);
xlabel('Time (ms)');
ylabel('mV');
title('K');
legend({'W/o' 'With'})
legend boxoff
xlim([180 300])
box off
set(gca, 'FontSize', 20)

%Subtract baseline
V1_2 = V1 - mean(V1(fix(150/dt):fix(155/dt)));
V2 = V2 - mean(V2(fix(150/dt):fix(155/dt)));

subplot(3,3,4)
plot(t,V1_2,'k',t,V2,'b','LineWidth',2);
xlabel('Time (ms)');
ylabel('mV');
xlim([180 300])
box off
set(gca, 'FontSize', 20)

%Normalize
V1_2 = V1_2/max(V1_2(pos_on:pos_off));
V2 = V2/max(V2(pos_on:pos_off));

subplot(3,3,7)
plot(t,V1_2,'k',t,V2,'b','LineWidth',2);
xlabel('Time (ms)');
ylabel('Normalized');
ylim([-0.1 1.1])
xlim([195 210])
box off
set(gca, 'FontSize', 20)




%% Condition 3
%Conductance mS/cm2
gak = 0; 
gnaT = 50;

% Initializing variable vector
V = zeros(loop,1);
nka = zeros(loop,1); %A-type K
lka = zeros(loop,1);
mna = zeros(loop,1); %Transient sodium variables
hna = zeros(loop,1);

% Set initial values for the variables
V(1) = VHold;
nka(1) = 0; %A-type K
lka(1) = 1;
mna(1) = 0; %Transient sodium variables
hna(1) = 1;

%Time constants
taunka = 0.2; %A-type K
taulka = 5;
mtauna = 0.1; %Transient sodium
htauna = 0.4;

%Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    ikA = gak*nka(i)^4*lka(i)*(eK-(V(i))); %A-type potassium
    inaT = gnaT*mna(i)*mna(i)*hna(i)*(eNa-V(i));

    V(i+1) = V(i) + (dt/Cm)*(iL + ikA + inaT + Iinj(i));

    %A-type potassium act/inact
    nka(i+1) = nka(i) +dt*((ninf_ka(V(i)) - nka(i))/taunka);
    lka(i+1) = lka(i) + dt*((linf_ka(V(i)) - lka(i))/taulka);

    %Transient sodium act/inact
    mna(i+1) = mna(i) + dt*((minf_na(V(i)) - mna(i))/mtauna);
    hna(i+1) = hna(i) + dt*((hinf_na(V(i)) - hna(i))/htauna);

end
V3 = V;


subplot(3,3,2)
plot(t,V1,'k',t,V3,'r','LineWidth',2);
xlabel('Time (ms)');
ylabel('mV');
title('Na');
legend({'W/o' 'With'})
legend boxoff
xlim([180 300])
box off
set(gca, 'FontSize', 20)

%Subtract baseline
V1 = V1 - mean(V1(fix(150/dt):fix(155/dt)));
V3 = V3 - mean(V3(fix(150/dt):fix(155/dt)));

subplot(3,3,5)
plot(t,V1,'k',t,V3,'r','LineWidth',2);
xlabel('Time (ms)');
ylabel('mV');
xlim([180 300])
box off
set(gca, 'FontSize', 20)

%Normalize
V1 = V1/max(V1(pos_on:pos_off));
V3 = V3/max(V3(pos_on:pos_off));

subplot(3,3,8)
plot(t,V1,'k',t,V3,'r','LineWidth',2);
xlabel('Time (ms)');
ylabel('Normalized');
ylim([-0.1 1.1])
xlim([195 210])
box off
set(gca, 'FontSize', 20)



%% Calcium
%% Condition 1
%Conductance mS/cm2
gcaT = 0;

% Initializing variable vectors
t = (1:loop)*dt;
V = zeros(loop,1);
mcaT = zeros(loop,1); %T-type Ca
hcaT = zeros(loop,1);

% Set initial values for the variables
V(1) = VHold;
mcaT(1) = 0; %T-type Ca
hcaT(1) = 1;

%Time constants ms
taumcaT = 0.5; %T-type Ca
tauhcaT = 10;

%Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium

    V(i+1) = V(i) + (dt/Cm)*(iL + icaT + Iinj(i));

    %T-type calcium act/inact
    mcaT(i+1) = mcaT(i) + dt*((minf_caT(V(i)) - mcaT(i))/taumcaT);
    hcaT(i+1) = hcaT(i) + dt*((hinf_caT(V(i)) - hcaT(i))/tauhcaT);

end
V1 = V;


%% Condition 2
%Conductance mS/cm2
gcaT = 0.2;

% Initializing variable vectors
V = zeros(loop,1);
mcaT = zeros(loop,1); %T-type Ca
hcaT = zeros(loop,1);

% Set initial values for the variables
V(1)=VHold;
mcaT(1) = 0; %T-type Ca
hcaT(1) = 1;

%Time constants ms
taumcaT = 0.5; %T-type Ca
tauhcaT = 10;

%Euler method
for i=1:loop-1

    %Currents
    iL = gL*(eL-V(i));
    icaT = gcaT*mcaT(i)*mcaT(i)*hcaT(i)*(eCa-(V(i))); %T-type calcium

    V(i+1) = V(i) + (dt/Cm)*(iL + icaT + Iinj(i));

    %T-type calcium act/inact
    mcaT(i+1) = mcaT(i) + dt*((minf_caT(V(i)) - mcaT(i))/taumcaT);
    hcaT(i+1) = hcaT(i) + dt*((hinf_caT(V(i)) - hcaT(i))/tauhcaT);

end
V2 = V;

%Plots
subplot(3,3,3)
plot(t,V1,'k',t,V2,'g','LineWidth',2);
xlabel('Time (ms)');
ylabel('mV');
legend({'W/o' 'With'})
legend boxoff
xlim([180 300])
box off
title('Ca')
set(gca, 'FontSize', 20)

%Subtract baseline
V1_2 = V1 - mean(V1(fix(150/dt):fix(155/dt)));
V2 = V2 - mean(V2(fix(150/dt):fix(155/dt)));

subplot(3,3,6)
plot(t,V1_2,'k',t,V2,'g','LineWidth',2);
xlabel('Time (ms)');
ylabel('mV');
xlim([180 300])
box off
set(gca, 'FontSize', 20)

%Normalize
V1_2 = V1_2/max(V1_2(pos_on:pos_off));
V2 = V2/max(V2(pos_on:pos_off));

subplot(3,3,9)
plot(t,V1_2,'k',t,V2,'g','LineWidth',2);
xlabel('Time (ms)');
ylabel('Normalized');
ylim([-0.1 1.1])
xlim([195 210])
box off
set(gca, 'FontSize', 20)

%Resize plot size
set(gcf,'units','inches','position',[0 0 16 8])

%Save vector figure in pdf
exportgraphics(gcf,'Step_3A.pdf','ContentType','vector');




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

