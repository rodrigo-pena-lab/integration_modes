
clear
clc
load mat3d_Ca.mat

C = mat3d(:,1,:);
C = permute(C,[1 3 2]);
C = reshape(C,[],size(C,2),1);
C = flip(C);
C = flip(C,2);

%plot
figure()
subplot(1,3,3)
contourf(C)
hcb = colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = '% Relat Freq';
set(colorTitleHandle ,'String',titleString);
xlabel('g_{K}  mS/cm2')
ylabel('g_{Na} mS/cm2 ')
clim([10 40]);
xcond = 0:20:100;
xticks(xcond)
xticklabels(0.01*xcond)
xcond = 0:20:100;
yticks(xcond)
yticklabels(0.5*xcond)
title('g_{Ca} = 0.1')
set(gca, 'FontSize', 20)

C = mat3d(:,2,:);
C = permute(C,[1 3 2]);
C = reshape(C,[],size(C,2),1);
C = flip(C);
C = flip(C,2);

subplot(1,3,2)
contourf(C)
hcb = colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = '% Relat Freq';
set(colorTitleHandle ,'String',titleString);
xlabel('g_{K}  mS/cm2')
ylabel('g_{Na} mS/cm2 ')
clim([10 40]);
xcond = 0:20:100;
xticks(xcond)
xticklabels(0.01*xcond)
xcond = 0:20:100;
yticks(xcond)
yticklabels(0.5*xcond)
title('g_{Ca} = 0.05')
set(gca, 'FontSize', 20)

C = mat3d(:,3,:);
C = permute(C,[1 3 2]);
C = reshape(C,[],size(C,2),1);
C = flip(C);
C = flip(C,2);

subplot(1,3,1)
contourf(C)
hcb = colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = '% Relat Freq';
set(colorTitleHandle ,'String',titleString);
xlabel('g_{K}  mS/cm2')
ylabel('g_{Na} mS/cm2 ')
clim([10 40]);
xcond = 0:20:100;
xticks(xcond)
xticklabels(0.01*xcond)
xcond = 0:20:100;
yticks(xcond)
yticklabels(0.5*xcond)
title('g_{Ca} = 0')
set(gca, 'FontSize', 20)
colormap('hot');

%Resize plot size
set(gcf,'units','inches','position',[0 0 16 4])

%Save vector figure in pdf
exportgraphics(gcf,'Step_17_heatmaps.pdf','ContentType','vector');

