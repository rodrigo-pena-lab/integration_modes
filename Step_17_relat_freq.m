
clear
clc

tic

load sparse_synchro_act.mat

%Current amplitude, frequency
Iamp = 5; %microA/cm2
freq_in = 4; %Hz

Tsim = 1000; %ms
VHold = -80;

%%
% Specific capacitance = 1 microF/cm2
Cm = 1;
dt = 0.05;               % time step for forward euler method (ms)
loop  = ceil(Tsim/dt);   % no. of iterations of euler
t = (1:loop)*dt;

%Reversal potential
eNa = 50; %mV
eK =-80;
eCa = 50;
eL = VHold; %Leak reversal potential

%Time constants ms
taunka = 0.2; %A-type K
taulka = 5;
taumcaT = 0.5; %T-type Ca
tauhcaT = 10;
mtauna = 0.05; %Transient sodium
htauna = 0.5;

%Conductance
gL = 0.3;

%LIF parameters
Vth = -40; %mV
Erest = -60; %mV

mat3d = zeros(3,3);
contk = 0;
for j = 0:0.01:1 %K conductance
    contk = contk+1;
    contca = 0;
    for k = 0.1:-0.05:0 %Ca conductance
        contca = contca + 1;
        contna = 0;
        for l = 0:0.5:50 %Sodium conductance
            contna = contna+1;

            all_rel_k = zeros(1,100);
            parfor p = 1:100

                %Load synaptic background signal (at 20 KHz)
                all_long_biphasic_traces = cell_traces{p};
                all_long_biphasic_traces_2 = cell_traces_2{p};

                %Conductance
                gak = j;
                gcaT = k;
                gnaT = l;

                %% Sparse activity
                %Injected current (oscillatory)
                Iinj = sin(2*pi*freq_in*t/1000);
                Iinj = Iinj + 1;
                Iinj = Iamp*Iinj + 0.5*all_long_biphasic_traces + 5;

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
                spks = find(tspk2==1);
                freq_K = length(spks);






                %% Synchronous
                %Injected current (oscillatory)
                Iinj = sin(2*pi*freq_in*t/1000);
                Iinj = Iinj + 1;
                Iinj = Iamp*Iinj + 0.5*all_long_biphasic_traces_2 + 5;

                %% Condition 2
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
                spks = find(tspk2==1);
                freq_K_2 = length(spks);


                %Frequency difference between sparse vs synchronous
                diff_k = freq_K_2 - freq_K;

                %Relative change
                rel_k = diff_k*100/freq_K;
                all_rel_k(p) = rel_k;

            end

            avg_k = mean(all_rel_k);
            mat3d(contk,contca,contna) = avg_k;

        end
    end
end

save mat3d_Ca_2.mat mat3d

toc

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
