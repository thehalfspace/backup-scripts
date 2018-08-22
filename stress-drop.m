% Vfmax for the entire simulation
Vfmax = zeros(it,1);
idx = zeros(it,1);

for eo = 1:it-1
    Vfmax(eo) = max(abs(SLIPVEL(:,eo)));
end

% Plot stress drop for Mw 7.6 earthquake
indx = 5983;
stress_drop = tauafter(:,indx) - taubefore(:,indx);
figure()
hold on;
plot(stress_drop, FaultX/1e3, 'Linewidth', 1.5);
a1 = plot(-20:20, -8.*ones(41,1), 'k--'); M1 = 'Fault Zone';
xlabel("Stress Drop (MPa)")
ylabel("Depth (km)")
title("Stress drop for Mw 3.7 earthquake")
ylim([-24/distN 0e3/distN])
legend(a1, M1)
hold off;
figure(gcf);


% Plot sliprate to look at rupture propagation
ind1 = 70786; ind2 = 72866;
slr = SLIPVEL(:,ind1:200:ind2);
figure()
hold on;
plot(slr, FaultX/1e3, 'Linewidth', 1.5);
a1 = plot(-5:15, -8.*ones(21,1), 'k--'); M1 = 'Fault Zone';
xlabel("Sliprate (MPa)")
ylabel("Depth (km)")
title("Rupture propagation for Mw 7.6 earthquake")
ylim([-24/distN 0e3/distN])
legend(a1, M1)
hold off;
figure(gcf);


