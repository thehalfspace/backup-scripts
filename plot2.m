% for 024.mat, idx1 = 154; idx2 = 155 + 200

idx1 = dynamic_it(154);
idx2 = qs_it(155)+200;
C = log10(abs(SLIPVEL(121:240, idx1:idx2)));

S = abs(STRESS(121:241, idx1:idx2));

xx = (time(idx1:idx2)- time(idx1));
yy = -1*FaultX(121:240)/1e3;
%yy = -1*FaultX/1e3;

figure();
subplot(3,1,1);
imagesc(xx,yy,C);
colorbar;
caxis([-4 2]);
title('Log(Slip rate (m/s))');
xlabel('Time (sec)');
ylabel('Depth (km)');

figure(gcf);

%figure();
subplot(3,1,2);
%plot(xx, log10(STRESS(216,idx1:idx2)), 'Linewidth', 2);
plot(xx, log10(abs(SLIPVEL(216,idx1:idx2))), 'Linewidth', 2);
title('Sliprate at 4.9 km depth')
xlabel('Time (sec)')
ylabel('Stress (MPa)')
ylim([-4.5 2])
xlim([-10 80])

subplot(3,1,3);
plot(xx, log10(abs(SLIPVEL(191,idx1:idx2))), 'Linewidth', 2);
title('Sliprate at 10 km depth')
xlabel('Time (sec)')
ylabel('Log Sliprate (m/s)')
ylim([-4.5 2])
xlim([-10 80])
figure(gcf)


%{
figure(); hold on;
for i = idx1:40:idx2
    plot(STRESS(:,i) + 0.1*(i-idx1), FaultX/1000, 'b-', 'Linewidth', 1);
    ylim([-24 0])
end
title('Stress profile for one seismic cycle')
xlabel('Shear stress')
ylabel('depth (km)')
hold off;
figure(gcf);
%}

% MFD
clear Mw; clear nMw; clear nMw2;
for i=3:length(slipvel_event(1,:))
    for j=239:481
        if slipvel_event(j,i)>=1e-3
            s = del_event(j,i) - del_event(j,i-1);
            
            Mo = 32*s;
            
            Mw(j,i) = (2/3)*(25 + log10(Mo)) - 10.7;
            
        end
    end
end

nMw2 = real(nonzeros(Mw));
nMw = nMw2(find(nMw2>5));
[counts, bins] = hist(nMw, 20);
scatter(bins, log10(counts), 'filled');
xlabel('Magnitude');
ylabel('Log of no. of earthquakes')
figure(gcf)


% Plot sliprate
 clear y_axis; clear x_axis;
for i=2:length(slipvel_event(1,:))
    [x_axis(i), idx] = max(slipvel_event(1:240,i));
     y_axis(i) = FaultX(idx)/1000;  
end

set(0,'DefaultAxesFontSize',14);
% histogram
figure();
hold on;
[counts, bins] = hist(y_axis, 20); 
barh(bins, counts, 'Facecolor', [0.4, 0.4, 1], 'EdgeColor', 'none'); 
a1 = plot(1:80, -0*ones(80,1), 'k--', 'Linewidth', 0.01);
M1 = 'Fault Zone Depth = Entire Domain';
hold off;
ylim([-18 0]); 
title('Hypocenter distribution of earthquakes')
xlabel('Number of earthquakes')
ylabel('Depth (km)')
legend(a1, M1)
box on
figure(gcf)

% estimate magnitude
k = 1;
mag1 = zeros(length(delfsec),1);
for i=1:max(size(Vfsec(1,:)))-10
    slr = Vfsec(213,i);
    [x_axis, idx] = max(delfsec(213,i));
    
    if slr <= 0.2 && slr>0.1
        mag1(k) = x_axis; k =k+1;
        
    end
end

mag1 = nonzeros(mag1);
mag1 = diff(mag1);


% estimate moment = integral mu.slip^2.element size
G = mu(1,1);
for j = 1: length(Vfsec(1,:))
    moment = 0;
for i = 237:481
    
    if Vfsec(i,j) > 1e-3
        moment = moment + G*delfsec(i,2000)^2 * dxe; 
    end
    
end
Mag1(j) = (2/3)*(log10(moment)) - 10.7;
end