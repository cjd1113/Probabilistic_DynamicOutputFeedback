
function [A,B1,B2,C1,C2,D11,D12,D21] = physicalcoords(M,C,K,PD,PC,Mm)

%Places the SMD system into standard, physical state space form.

A = [zeros(size(M,1)) eye(size(M,1));...
    -M\K -M\C];

B1 = [zeros(size(M,1),size(PD,2));... %disturbance input
    M\PD]; 

B2 = [zeros(size(M,1),size(PC,2));... %control input
    M\PC];

C1 = [Mm zeros(size(Mm,1),size(Mm,2));...
    zeros(size(Mm,1),size(Mm,2)), Mm]; %perf. output
C2 = C1; %measurement output

D11 = zeros(size(C1,1),size(B1,2));
D12 = zeros(size(C1,1),size(B2,2));
D21 = zeros(size(C2,1),size(B2,2));

return