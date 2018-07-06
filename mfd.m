% idx 441 = 4 km depth
%for p = 300:481

clear init; clear endd; clear slp;
    
j = 1; k = 1; flag = 1; p = 325;
slip = 0;
%init = 0; endd = 0;
    
    
for i = 1: length(SLIPVEL(1,:))
    
    if SLIPVEL(p, i) >= 1e-3
        if flag == 1 
            %slip = slip + SLIPVEL(441, i)*(time(i) - time(i-1));
            init(j) = i; j = j+ 1;
            flag = 0;
        end
    end
    
    if SLIPVEL(p, i)< 1e-3
        if flag == 0
            endd(k) = i; k = k+1;
            flag = 1;
        end
    end
end
init = init';
endd = endd';

%plot(time, SLIPVEL(p,:))


%sv = SLIPVEL(SLIPVEL > 1e-3)

for i = 2: length(endd)
    idx1 = init(i);
    idx2 = endd(i);
    
    slip = 0;
    for j = idx1:idx2
        slip = slip + SLIPVEL(p, j)*(time(j) - time(j-1));
    end
    
    slp(i) = slip;
end

%end
%plot(1:k-1, slp); figure(gcf)

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

[counts, bins] = hist(Mw, 10);
figure()
scatter(bins, log10(counts), 'filled');
xlabel('Magnitude');
ylabel('Log of no. of earthquakes')
figure(gcf)