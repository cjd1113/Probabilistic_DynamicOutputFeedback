
function [Ccc,V,D] = addDAMPING(Mcc,Kcc,zeta_r)

Nwn = size(Kcc,1);

%Solve the generalized eigenvalue problem
[V,D] = eigs(Kcc,Mcc,Nwn);

MM = V'*Mcc*V;
[m,n] = size(V);
tol = 10*eps;

%tolerance defined as an order of magnitude greater than machine epsilon to
%take care of numerical roundoff errors and ease computational complexity
for i = 1:m
    for j = 1:n
        if MM(i,j) < tol
            MM(i,j) = 0;
        end
    end
end


Mr = zeros(size(MM));
for i = 1:n
    Mr(i,i) = MM(i,i);
end

CC = 2*zeta_r*Mr*sqrt(D); %Add zeta_r damping to each mode

C = (Mcc*V/MM)*CC*(MM\(V')*Mcc); %Transformation to get back to physical coordinates - See Craig1981book page 450

Ccc = C;

end

