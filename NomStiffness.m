
%This function generates the uncertain interconnection stiffness matrix for
%an Euler Bernoulli beam element, where we have restricted this uncertain
%object to be +/- p_uncert% around the nominal value, E_nom

%Prepared by: Chris D'Angelo
%Date: May 21, 2018

function [K11nom,K12nom,K21nom,K22nom] = NomStiffness(width,height,length,E_nom)


E = E_nom;
L = length;
I = 1/12*width*height^3;

Ke = (E*I)/(L^3) * [12 6*L -12 6*L;...
    6*L 4*L^2 -6*L 2*L^2;...
    -12 -6*L 12 -6*L;...
    6*L 2*L^2 -6*L 4*L^2];

%Partition the stiffness element so that we have the block diagonal and
%cross terms for coupling interface forces / moments and displacements /
%rotations

K11nom = Ke(1:2,1:2);
K12nom = Ke(1:2,3:4);
K21nom = Ke(3:4,1:2);
K22nom = Ke(3:4,3:4);

end