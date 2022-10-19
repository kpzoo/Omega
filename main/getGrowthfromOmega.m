% Compute growth rate from omega estimates
function romega = getGrowthfromOmega(A, B, omega)

% Assumptions and notes
% - inputs refer to transformation under exp growth/decline
% - A is window length scalar, B is vector of omega values

% Input to Lambert W function
C = -B.*exp(-B);

% Start by assuming all solutions on principal Lambert branch
romega = A*(B + lambertw(C));

% Determine solutions on -1 branch of Lambert
idneg = find(omega < 1);
romega(idneg) = A*(B(idneg) + lambertw(-1, C(idneg)));
romega = real(romega);

