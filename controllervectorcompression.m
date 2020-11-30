

function [ControllerParameters, zeroimaginarypartindices] = ...
    controllervectorcompression(Avectors,Bvectors,Cvectors,Dvectors,numrealeigs)

%%This function removes the inherent redundancy of complex conjugate
%%eigenvalue pairs.  Although there are 2 real and 2 complex parts to every
%%complex conjugate pair, we know that we only need 2 real numbers to
%%characterize each eigenvalue.  Therefore, we can perform this compression
%%to greatly reduce the number of GA parameters that we have.

%Prepared by: Chris D'Angelo
%Date: August 24, 2018

Avec_complexpart = Avectors(numrealeigs+1:end,:);
Bvec_complexpart = Bvectors(numrealeigs+1:end,:);
Cvec_complexpart = Cvectors(:,numrealeigs+1:end);
Dvec_complexpart = Dvectors;

%These are the "parsed" complex conjugate parts of the controller vectors.
%We have effectively halved their dimension since our eigenvalues /
%eigenvectors for complex numbers come in conjugate pairs.  
Avec_complexpartreduced = Avec_complexpart(1:2:end,:);
Bvec_complexpartreduced = Bvec_complexpart(1:2:end,:);
Cvec_complexpartreduced = Cvec_complexpart(:,1:2:end);
Dvec_complexpartreduced = Dvec_complexpart;

%Now, we need to recombine with the real-valued eigenvalues.

Avec_realcomplexreduced = [Avectors(1:numrealeigs,:);Avec_complexpartreduced];
Bvec_realcomplexreduced = [Bvectors(1:numrealeigs,:);Bvec_complexpartreduced];
Cvec_realcomplexreduced = [Cvectors(:,1:numrealeigs),Cvec_complexpartreduced];
Cvec_realcomplexreduced = Cvec_realcomplexreduced'; %note that we transpose the output vector to get it into a form that is more amenable to vectorization.  This is all just a lot of bookkeeping.
Dvec_realcomplexreduced = Dvec_complexpartreduced;


%Now, we need to vectorize these controller variables.

%These look like:
%[Areal; Acomplex]
%[Brealcolumns; Bcomplexcolumns]
%[Crealrows; Ccomplexrows];
Avectorized = Avec_realcomplexreduced(:);
Bvectorized = Bvec_realcomplexreduced(:);
Cvectorized = Cvec_realcomplexreduced(:);
Dvectorized = Dvec_realcomplexreduced(:); %Note that we will omit / set this term to zero.  We restrict our controller to be strictly proper.

ControllerParameters = [Avectorized;Bvectorized;Cvectorized];

%We now index where the zero-valued imaginary parts are located.
zeroimaginarypartindices = find(ControllerParameters == 0); %This will be declared as a global variable, useful for controller matrix reconstruction


end