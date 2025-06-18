
clear
clc
 
V = -120:1:50;

%A-type K variables
n1ka = 1./(1 + exp(-(V+50)/12));
l1ka = 1./(1 + exp((V+60)/12));

figure()
subplot(1,3,1)
plot(V,n1ka,V,l1ka,'LineWidth',2)
xlabel('Voltage')
legend({'Activation' 'Inactivation'})
legend boxoff
box off
xlim([-120 50])
title('Potassium')


%T-type Ca variables
m1ca = 1./(1 + exp(-(V+40)/10));
h1ca = 1./(1 + exp((V+60)/10));

subplot(1,3,2)
plot(V,m1ca,V,h1ca,'LineWidth',2)
xlabel('Voltage')
legend({'Activation' 'Inactivation'})
legend boxoff
box off
xlim([-120 50])
title('Calcium')


%Transient sodium variables
n1na = 1./(1 + exp(-(V+28)/5));
l1na = 1./(1 + exp((V+62)/5));

subplot(1,3,3)
plot(V,n1na,V,l1na,'LineWidth',2)
xlabel('Voltage')
legend({'Activation' 'Inactivation'})
legend boxoff
box off
xlim([-120 50])
title('Sodium')


results = [V;n1ka;V;l1ka];

csvwrite('raw/results.csv',results)