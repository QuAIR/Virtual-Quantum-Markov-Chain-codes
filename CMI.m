%% Comditional Mutual Information
function cmi = CMI(rhoABC)
d = size(rhoABC);
da = d(1)^(1/3);
db = da;
dc = da;

rhoAB = PartialTrace(rhoABC, 3, [da db dc]);
rhoBC = PartialTrace(rhoABC, 1, [da db dc]);
rhoB = PartialTrace(rhoAB, 1, [da db]);

cmi = Entropy(rhoAB) + Entropy(rhoBC) - Entropy(rhoABC) - Entropy(rhoB);
end