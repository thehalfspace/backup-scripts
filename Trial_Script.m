function [TN,TNX,TNXX,TNXXX,TN4X]=TN_BASIS(T,JJ)
%This evaluates the JJ-th basis function at the point whose
% trigonometric argument is T. The corresponding argument for
% the Tn(x) as rational functions of x is
%  x = L cos(T)
% Needs to be independenty checked
TN = cos(JJ*T);  PT = - JJ * sin(JJ*T);
PTT = -JJ*JJ*TN;  PTTT  = -JJ*JJ*PT;  P4T  = -JJ*JJ*PTT ;
% Now transform from T-derivatives to Y-derivatives
S = sin(T);   C = cos(T);
TNX = - PT / S ;
TNXX = (S*PTT - C*PT)/(S*S*S) ;
TNXXX =(PT*(-S*S-3*C*C) + PTT*3*C*S - PTTT*S*S ) / S^5 ;
TN4X=(PT*(-15*C*C*C-9*C*S*S) + PTT*(4*S^4+15*C*C*S) -6*C*S*S*PTTT+S*S*S*P4T)/(S^7);



%NRsearch_NEW(fo(j),Vo(j),cca(j),ccb(j),Seff(j),tauNR(j),tauo(j),psi1(j),FaultZ(j),FaultVFree(j))

%%%%%%%%%%%%%%%%%%
% Input Parameters
fo = 0.6;
Vo = 1e-6;
cca = 0.0135;
ccb = 0.01;
Seff = 4.34e6;
tau = 0;
tauo = 3e6;
psi = -16.0682;
FaultZ = 4.18e6;
FaultVFree = 9.9586
%%%%%%%%%%%%%%%%%%


Vw = 10^10;
fact = 1+(Vo/Vw)*exp(-psi);

cA = cca*Seff;
eps = 0.001*cA;
k=0;
delta=inf;

while (abs(delta) > eps)
    fa = fact*tau/(Seff*cca);
    help = -(fo+ccb*psi)/cca;
    help1 = exp(help+fa);
    help2 = exp(help-fa);   
    Vf = Vo*(help1 - help2);
    Vfprime = fact*(Vo/cA)*(help1 + help2);
    delta = (FaultZ*FaultVFree - FaultZ*Vf + tauo - tau)/(1 + FaultZ*Vfprime);
    tau = tau + delta;
    k = k+1;
    if (abs(delta) > 10^10 || k == 1000)
        k
        Vf
        delta
        tau
        %Vfrpos
        error('N-R search fails to converge');
    end
end

%{
fa = fact*tau/(Seff*cca);
help = -(fo+ccb*psi)/cca;
help1 = exp(help+fa);
help2 = exp(help-fa);
Vf = Vo*(help1 - help2);
%}