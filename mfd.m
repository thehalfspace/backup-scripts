% idx 441 = 4 km depth
%for p = 300:481

clear moment; clear Mw; slip = zeros(181);
%% Calculate slip on each element with time, and then integrate with depth
for i = 2:328
      
    idx1 = dynamic_it(i);
    idx2 = qs_it(i+1);

    mo = 0; j = 1;
    net_area = 0; thres = 0;
    for depth = 300:481
        for k = idx1:idx2
            slip(j) = slip(j) + SLIPVEL(depth, k)*(time(k) - time(k-1));
        end
        j = j + 1;
    end
    
    thres = max(slip);
    net_area = sum(slip(slip> 0.01*thres)*dxe^2);
       
    Moment(i) = mu(1,1)*net_area;
end
Moment = Moment';


%% Moment calculation 2
for i = 2:328
      
    idx1 = dynamic_it(i);
    idx2 = qs_it(i+1);

    mo = 0;
    
    for k = idx1:idx2
        slip = sum(abs(SLIPVEL(300:481, k)))*(time(k) - time(k-1));
        mo = mo + slip*dxe^2;
    end
       
    Moment(i) = mu(1,1)*mo;
end
Moment = Moment';

Mw = (2/3)*log10(Moment.*1e7) - 10.7;

[counts, bins] = hist(Mw);
figure()
scatter(bins, log10(counts), 'filled');
xlabel('Magnitude');
ylabel('Log of no. of earthquakes')
figure(gcf)