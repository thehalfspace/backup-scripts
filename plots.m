set(0, 'DefaultFigureVisible', 'on')

yr2sec = 31536000;

% Plot initial stresses
figure(); hold on
a1 = plot(Seff/1e6, FaultX/1000, 'Linewidth', 1.5); M1 = 'Normal Stress';
a2 = plot(tauo/1e6, FaultX/1000, 'Linewidth', 1.5); M2 = 'Shear Stress';
a3 = plot((cca-ccb)*1e3 + 10, FaultX/1000, 'Linewidth', 1.5); M3 = 'a-b';
legend([a1; a2; a3], M1, M2, M3);
title('Effective normal stress and initial shear stress')
xlabel('Stress (MPa)')
ylabel('Depth (Km)')
%xlim([50e6 150e6])

plot(1:50, repmat(-1.8,1,50), 'k--')
plot(1:50, repmat(-12,1,50), 'k--')
plot(repmat(30,1,25), -24:0, 'k--')
plot(repmat(10,1,25), -24:0, 'k--')
ylim([-24/distN 0/distN])
hold off
figure(gcf)


% Nucleation size
G = repmat(RHO2*VS2^2, FaultNglob,1);

G(find(abs(FaultX)<=LX - THICKX)) = RHO2*VS2^2;

h_RA = (2/pi)*(G.*xLf.*ccb)./(Seff.*((ccb - cca).^2));
plot(h_RA/1000, FaultX/1000, 'Linewidth', 1.5);
title('Nucleation size')
xlabel('h_{RA} (km)')
ylabel('Depth (km)')
ylim([-24/distN 0e3/distN])
figure(gcf)

% Plot friction paramters
figure(); hold on
a1 = plot((cca-ccb), FaultX/1000, '-o', 'Linewidth', 1.5); M1 = 'a-b';
a2 = plot(cca, FaultX/1000, '-o', 'Linewidth', 1.5); M2 = 'a';
a3 = plot(ccb, FaultX/1000, '-o', 'Linewidth', 1.5); M3 = 'b';
%a4 = plot(stress_event(:,3), FaultX/1000, 'Linewidth', 1.5); M4 = 'stress';
title('Rate and state friction parameters a and b')
xlabel('Scaled friction and stress')
ylabel('Depth')
legend([a1; a2; a3], M1, M2, M3);
ylim([-24/distN 0e3/distN])
hold off
figure(gcf)

%   Plot shear stress as a function of time
figure();
plot(time(1:it-1)/yr2sec, Tauloc1, 'Linewidth', 2);
title('Shear stress (D_c = 8 mm)')
xlabel('Time (years)')
ylabel('tau - (fo*Sigma_{eff})')
%figname = sprintf('a_bShearStress00%d.jpg',iter001);
%saveas(gcf, figname)
figure(gcf)

set(0,'DefaultAxesFontSize',14);

%   Plot slip rate as a function of time
figure1 = figure; 
plot(time(1:it)/yr2sec, log10(SLIPVEL(191,:)), 'Linewidth', 2);
title('Slip rate at 5.6 km depth')
xlabel('Time (years)')
ylabel('Log of Slip rate (m/s)')
%xlim([230 680])
ylim([-4 2])
% Create textarrow
annotation(figure1,'textarrow',[0.671945701357466 0.763574660633489],...
    [0.845070422535211 0.77599530516432],'String',{'Mw ~ 7'},'LineWidth',2,...
    'FontSize',18,...
    'FontName','Helvetica Neue');
% Create textarrow
annotation(figure1,'textarrow',[0.315610859728507 0.329185520361991],...
    [0.605633802816901 0.401408450704225],'String',{'Mw ~ 4'},'LineWidth',2,...
    'FontSize',18,...
    'FontName','Helvetica Neue');
box on;
figure(gcf)

% Plot cumulative slip
% create a default color map ranging from red to light pink
c1 = [255, 196, 0]/255; % 
c2 = [204, 102, 0]/255; % 
c3 = [204, 204, 0]/255; % 
c4 = [253, 253, 0]/255;  % 
c5 = [204, 102, 0]/255; % 

lt = 3;
red = [255, 162, 0]/255;%[0.85, 0.325, 0.098];
pink = [153, 76, 0]/255;
colors_p = [linspace(red(1),pink(1),lt)',...
            linspace(red(2),pink(2),lt)',...
            linspace(red(3),pink(3),lt)'];

c1 = colors_p(1,:);
c2 = colors_p(2,:);
c3 = colors_p(3,:);
%c4 = colors_p(4,:);
%c5 = colors_p(5,:);

figure();
hold on;
for i=1:length(delf5yr(1,:))
    plot(delf5yr(:,i),FaultX/1000, 'Color', [0, 0.447, 0.741], 'LineWidth', 1)
end

for i=1:length(delfsec(1,:))
    [sliprate, idx] = max(Vfsec(1:480,i));
    x_axis = delfsec(311:end,i);
    y_axis = FaultX(311:end)/1000;
    
    if sliprate <= 0.01
        plot(x_axis, y_axis, 'Color', c3, 'LineWidth', 1)
    elseif sliprate <=0.1
        plot(x_axis, y_axis, 'Color', c2, 'LineWidth', 1)
    %elseif sliprate <=1
    %    plot(x_axis, y_axis, 'Color', c2, 'LineWidth', 1)
    %elseif sliprate <=1
    %    plot(x_axis, y_axis, 'Color', c2, 'LineWidth', 1)
    else
        plot(x_axis, y_axis, 'Color', c1, 'LineWidth', 1)
    end
    %plot(delfsec(155:end,i),FaultX(155:end)/1000, 'r', 'LineWidth', 1)
end

hold off;

%h1 = findobj('Color', c1); M1 = 'Slip rate <= 0.01';
%h2 = findobj('Color', c2); M2 = 'Slip rate <= 0.05';
%h3 = findobj('Color', c3); M3 = 'Slip rate <= 0.1';
%h4 = findobj('Color', c4); M4 = 'Slip rate <= 0.5';
%h5 = findobj('Color', c5); M5 = 'Slip rate > 0.5';

%legend([h1; h2; h3; h4; h5], M1, M2, M3, M4, M5);

title('Fault Zone depth = 9.6 km')
xlabel('Accumulated slip (m)')
ylabel('Depth (km)')
ylim([-24 0])
%xlim([12 35])
box on;
figure(gcf)

figure();
hold on;
plot(1:10, 5.*ones(10,1), 'Color', c1, 'Linewidth', 2.0);
plot(1:10, 6.*ones(10,1), 'Color', c2, 'Linewidth', 2.0);
plot(1:10, 7.*ones(10,1), 'Color', c3, 'Linewidth', 2.0);
plot(1:10, 9.*ones(10,1), 'Color', [0, 0.447, 0.741], 'Linewidth', 2.0);
plot(1:10, 3.*ones(10,1), 'k', 'Linewidth', 2.0);
plot(1:10, 10.*ones(10,1), 'k', 'Linewidth', 2.0);
hold off; figure(gcf)

%%{
subplot(1,2,2);
figure();
plot(SLIPVEL(:,4130), FaultX/1000, 'Linewidth', 1.0); 
ylim([-24 0])
title('Slip velocity at the onset of event')
xlabel('slip velocity')
ylabel('Depth (km)')
%}
figure(gcf)


% Plot slip velocity event
figure(); hold on;
for i = 2:62
    plot(slipvel_event(:,i), FaultX/1000, 'o', 'Linewidth', 1);
    xlim([0 0.01]); ylim([-24 0])
end
title('Slip at the onset of nucleation')
xlabel('Slip velocity')
ylabel('depth (km)')
hold off;
figure(gcf);

% Plot stress event
figure(); hold on;
for i = 2:62
    plot(stress_event(:,i) + 4*i, FaultX/1000, 'b-', 'Linewidth', 1);
    ylim([-24 0])
end
title('Stress profile for each event')
xlabel('Shear stress')
ylabel('depth (km)')
hold off;
figure(gcf);


%   Plot slip rate as a function of time
idx1 = dynamic_it(246);
idx2 = qs_it(249);

figure();
subplot(2,1,1);
plot(time(idx1:idx2)/yr2sec, log10(SLIPVEL(215,idx1:idx2)), 'Linewidth', 2);
title('Slip rate as a function of time (10 km depth)')
xlabel('Time (years)')
ylabel('Log of Slip rate (m/s)')

figure(gcf)

%idx1 = dynamic_it(173);
%idx2 = dynamic_it(176);
% Plot slip velocity
subplot(2,1,2);
figure();
k=0;
hold on;
for i= idx1:10:idx2    %66770:200:100662
    
    %for j = 1:239
        %if (SLIPVEL(j,i)) > 0.5
            k=k+1;
            plot(SLIPVEL(1:240,i)+k, FaultX(1:240)/1000, 'Linewidth', 1.0);
            plot(zeros(1,240)+k, FaultX(1:240)/1000, '--');
        %break;
    %end
    %end
end
hold off;
title('Slip velocity vs Depth every 50 iterations')
xlabel('Slip velocity (m/s)')
ylabel('depth (km)')
ylim([-24 0])
%xlim([0 900e-9])
figure(gcf)

figure();
plot(SLIPVEL(:,idx1+500), FaultX/1000, 'Linewidth', 2.0);
title('Slip velocity vs Depth (Time = 346.2 years)')
xlabel('Slip velocity (m/s)')
ylabel('depth (km)')
ylim([-16 0])
figure(gcf)


% Plot stress in time
figure();
%hold on;
tt = 5920; i = 1;
plot(1:50, repmat(-1.6,1,50), 'k--')
plot(1:50, repmat(-12,1,50), 'k--')
a1 = plot(1:50, repmat(-5.6,1,50), 'b--', 'LineWidth', 2); M1 = 'Fault Zone';
plot(repmat(30,1,25), -24:0, 'k--')
for tt =8000:11000 %[4500, 5000, 5500, 6000, 6400]
    %subplot(1,5,i);
    hold on;
    plot(1:40, repmat(-1.6,1,40), 'k--','LineWidth', 1)
    plot(1:40, repmat(-12,1,40), 'k--','LineWidth', 1)
    a1 = plot(1:40, repmat(-5.6,1,40), 'b--', 'LineWidth', 2); M1 = 'Fault Zone';
    plot(repmat(30,1,25), -24:0, 'k--','LineWidth', 1)
    plot(STRESS(:,tt),FaultX/1000, 'LineWidth', 1.5)
    i = i+1;
    %title(['Time =', num2str(time(tt)), 's, or ', num2str(time(tt)/yr2sec),'years'])
    ylim([-24 0])
    xlim([0 40])
    figure(gcf)
    pause(0.01)
end
legend([a1], M1)
title('Shear Stress during seismic event (t = 131.9733 y) every 400 iteratiopns/timesteps')
xlabel('Shear Stress (MPa)')
ylabel('depth (km)')
hold off;
