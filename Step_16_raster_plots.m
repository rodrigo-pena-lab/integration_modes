
clear
clc

load sparse_synchro_act.mat

all_spks_sparse_pass = [];
all_spks_sparse_k = [];
all_spks_sparse_Ca = [];
all_spks_sparse_Na = [];

all_spks_sync_pass = [];
all_spks_sync_k = [];
all_spks_sync_Ca = [];
all_spks_sync_Na = [];


Iamp = 5; %microA/cm2
freq_in = 4; %Hz

Tsim = 1000; %ms
VHold = -80;

%Reversal potential
eNa = 50; %mV
eK =-80;
eCa = 50;
eL = VHold; %Leak reversal potential

% Specific capacitance = 1 microF/cm2
Cm = 1;
dt = 0.05;               % time step for forward euler method (ms)
loop  = ceil(Tsim/dt);   % no. of iterations of euler
t = (1:loop)*dt;

%Time constants ms
taunka = 0.2; %A-type K
taulka = 5;
taumcaT = 0.5; %T-type Ca
tauhcaT = 10;
mtauna = 0.05; %Transient sodium
htauna = 0.5;

%LIF parameters
Vth = -40; %mV
Erest = -60; %mV

total_neurons = 100;
for p = 1:total_neurons

    %Load synaptic background signal (at 20 KHz)
    all_long_biphasic_traces = cell_traces{p};
    all_long_biphasic_traces_2 = cell_traces_2{p};


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

    tspk1 = zeros(loop,1);

    %Injected current (oscillatory)
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
    m = zeros(loop,1);
    h = zeros(loop,1);
    n = zeros(loop,1);
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
    V(1) = VHold;
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
    V(1)=VHold;
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

    all_spks_sparse_pass = [all_spks_sparse_pass; tspk1_sparse'];
    all_spks_sparse_k = [all_spks_sparse_k; tspk2_sparse'];
    all_spks_sparse_Ca = [all_spks_sparse_Ca; tspk3_sparse'];
    all_spks_sparse_Na = [all_spks_sparse_Na; tspk4_sparse'];








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
    V(1) = VHold;
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


    all_spks_sync_pass = [all_spks_sync_pass; tspk1_sync'];
    all_spks_sync_k = [all_spks_sync_k; tspk2_sync'];
    all_spks_sync_Ca = [all_spks_sync_Ca; tspk3_sync'];
    all_spks_sync_Na = [all_spks_sync_Na; tspk4_sync'];

end


%% Passive raster plots
%PSTH data
sum_sparse = sum(all_spks_sparse_pass);
sum_sync = sum(all_spks_sync_pass);
psth_sparse = zeros(1,100);
psth_sync = zeros(1,100);
for i = 1:100 %10 ms bin (200 points)
    psth_sparse(i) = sum(sum_sparse(200*(i-1)+1:200*i));
    psth_sync(i) = sum(sum_sync(200*(i-1)+1:200*i));
end
psth_sparse = psth_sparse/10;
psth_sync = psth_sync/10;

%Plot raster plot
figure()
subplot(4,4,1)
hold on
for i = 1:total_neurons
    tsp = all_spks_sparse_pass(i,:);
    tsp = find(tsp==1);
    plot(tsp*dt,i*ones(1,length(tsp)),'.b','linewidth',1)
end
xlabel('Spikes time (ms)')
ylabel('Neuron id')
ylim([0 (total_neurons+1)])
title('Passive Sparse')
xlim([200 450])

subplot(4,4,5)
hold on
for i = 1:total_neurons
    tsp = all_spks_sync_pass(i,:);
    tsp = find(tsp==1);
    plot(tsp*dt,i*ones(1,length(tsp)),'.r','linewidth',1)
end
xlabel('Spikes time (ms)')
ylabel('Neuron id')
ylim([0 (total_neurons+1)])
title('Passive Synchronous')
xlim([200 450])

subplot(4,4,9)
hold on
b2 = bar(10*(1:length(psth_sync)),psth_sync,'r');
b2.FaceAlpha = 0.5;
b1 = bar(10*(1:length(psth_sparse)),psth_sparse,'b');
b1.FaceAlpha = 0.5;
box off
xlim([0 1000])
xlabel('ms')
ylabel('Spiking freq (ms-1)')
legend({'Sync' 'Sparse'})
legend boxoff
xlim([200 450])
ylim([0 15])

subplot(4,4,13)
hold on
bar(10*(1:length(psth_sync)),100*(psth_sync-psth_sparse)./psth_sparse,'k');
box off
xlim([0 1000])
xlabel('ms')
ylabel('% diff spiking')
xlim([200 450])
ylim([-50 200])



%% Potassium raster plots
%PSTH data
sum_sparse = sum(all_spks_sparse_k);
sum_sync = sum(all_spks_sync_k);
psth_sparse = zeros(1,100);
psth_sync = zeros(1,100);
for i = 1:100 %10 ms bin (200 points)
    psth_sparse(i) = sum(sum_sparse(200*(i-1)+1:200*i));
    psth_sync(i) = sum(sum_sync(200*(i-1)+1:200*i));
end
psth_sparse = psth_sparse/10;
psth_sync = psth_sync/10;

%Plot raster plot
subplot(4,4,2)
hold on
for i = 1:total_neurons
    tsp = all_spks_sparse_k(i,:);
    tsp = find(tsp==1);
    plot(tsp*dt,i*ones(1,length(tsp)),'.b','linewidth',1)
end
xlabel('Spikes time (ms)')
ylabel('Neuron id')
ylim([0 (total_neurons+1)])
title('K Sparse')
xlim([200 450])

subplot(4,4,6)
hold on
for i = 1:total_neurons
    tsp = all_spks_sync_k(i,:);
    tsp = find(tsp==1);
    plot(tsp*dt,i*ones(1,length(tsp)),'.r','linewidth',1)
end
xlabel('Spikes time (ms)')
ylabel('Neuron id')
ylim([0 (total_neurons+1)])
title('K Synchronous')
xlim([200 450])

subplot(4,4,10)
hold on
b2 = bar(10*(1:length(psth_sync)),psth_sync,'r');
b2.FaceAlpha = 0.5;
b1 = bar(10*(1:length(psth_sparse)),psth_sparse,'b');
b1.FaceAlpha = 0.5;
box off
xlim([0 1000])
xlabel('ms')
xlim([200 450])
ylim([0 15])

subplot(4,4,14)
hold on
bar(10*(1:length(psth_sync)),100*(psth_sync-psth_sparse)./psth_sparse,'k');
box off
xlim([0 1000])
xlabel('ms')
ylabel('% diff spiking')
xlim([200 450])
ylim([-50 200])




%% Calcium raster plots
%PSTH data
sum_sparse = sum(all_spks_sparse_Ca);
sum_sync = sum(all_spks_sync_Ca);
psth_sparse = zeros(1,100);
psth_sync = zeros(1,100);
for i = 1:100 %10 ms bin (200 points)
    psth_sparse(i) = sum(sum_sparse(200*(i-1)+1:200*i));
    psth_sync(i) = sum(sum_sync(200*(i-1)+1:200*i));
end
psth_sparse = psth_sparse/10;
psth_sync = psth_sync/10;

%Plot raster plot
subplot(4,4,3)
hold on
for i = 1:total_neurons
    tsp = all_spks_sparse_Ca(i,:);
    tsp = find(tsp==1);
    plot(tsp*dt,i*ones(1,length(tsp)),'.b','linewidth',1)
end
xlabel('Spikes time (ms)')
ylabel('Neuron id')
ylim([0 (total_neurons+1)])
title('Ca Sparse')
xlim([200 450])

subplot(4,4,7)
hold on
for i = 1:total_neurons
    tsp = all_spks_sync_Ca(i,:);
    tsp = find(tsp==1);
    plot(tsp*dt,i*ones(1,length(tsp)),'.r','linewidth',1)
end
xlabel('Spikes time (ms)')
ylabel('Neuron id')
ylim([0 (total_neurons+1)])
title('Ca Synchronous')
xlim([200 450])

subplot(4,4,11)
hold on
b2 = bar(10*(1:length(psth_sync)),psth_sync,'r');
b2.FaceAlpha = 0.5;
b1 = bar(10*(1:length(psth_sparse)),psth_sparse,'b');
b1.FaceAlpha = 0.5;
box off
xlim([0 1000])
xlabel('ms')
ylabel('Spiking freq (ms-1)')
xlim([200 450])

subplot(4,4,15)
hold on
bar(10*(1:length(psth_sync)),100*(psth_sync-psth_sparse)./psth_sparse,'k');
box off
xlim([0 1000])
xlabel('ms')
ylabel('% diff spiking')
xlim([200 450])
ylim([-50 200])




%% Sodium raster plots
%PSTH data
sum_sparse = sum(all_spks_sparse_Na);
sum_sync = sum(all_spks_sync_Na);
psth_sparse = zeros(1,100);
psth_sync = zeros(1,100);
for i = 1:100 %10 ms bin (200 points)
    psth_sparse(i) = sum(sum_sparse(200*(i-1)+1:200*i));
    psth_sync(i) = sum(sum_sync(200*(i-1)+1:200*i));
end
psth_sparse = psth_sparse/10;
psth_sync = psth_sync/10;

%Plot raster plot
subplot(4,4,4)
hold on
for i = 1:total_neurons
    tsp = all_spks_sparse_Na(i,:);
    tsp = find(tsp==1);
    plot(tsp*dt,i*ones(1,length(tsp)),'.b','linewidth',1)
end
xlabel('Spikes time (ms)')
ylabel('Neuron id')
ylim([0 (total_neurons+1)])
title('Na Sparse')
xlim([200 450])

subplot(4,4,8)
hold on
for i = 1:total_neurons
    tsp = all_spks_sync_Na(i,:);
    tsp = find(tsp==1);
    plot(tsp*dt,i*ones(1,length(tsp)),'.r','linewidth',1)
end
xlabel('Spikes time (ms)')
ylabel('Neuron id')
ylim([0 (total_neurons+1)])
title('Na Synchronous')
xlim([200 450])

subplot(4,4,12)
hold on
b2 = bar(10*(1:length(psth_sync)),psth_sync,'r');
b2.FaceAlpha = 0.5;
b1 = bar(10*(1:length(psth_sparse)),psth_sparse,'b');
b1.FaceAlpha = 0.5;
box off
xlim([0 1000])
xlabel('ms')
ylabel('Spiking freq (ms-1)')
xlim([200 450])

subplot(4,4,16)
hold on
bar(10*(1:length(psth_sync)),100*(psth_sync-psth_sparse)./psth_sparse,'k');
box off
xlim([0 1000])
xlabel('ms')
ylabel('% diff spiking')
xlim([200 450])
ylim([-50 200])



%Resize plot size
set(gcf,'units','inches','position',[0 0 16 8])

%Save vector figure in pdf
exportgraphics(gcf,'Step_4A_raster.pdf','ContentType','vector');





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
