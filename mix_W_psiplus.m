clear;

%% Required

% CVX: http://www.cvxr.com
% QETLAB: http://www.qetlab.com

%% prepare states
da = 2;
db = 2;
dc = 2;

v0 = [1;0];
v1 = [0;1];

v01 = kron(v0, v1);
v10 = kron(v1, v0);
v = (v01 + v10)/sqrt(2);

Psiplus = v * v';
Psiplus = Psiplus/trace(Psiplus);

rAB =  kron(Psiplus, eye(2)/2);

v_w = 1/sqrt(3)*(Tensor(v0,v0,v1) + Tensor(v0,v1,v0) + Tensor(v1,v0,v0));  % W state

cost_list = [];
cmi_list = [];
prob = 0:0.02:1;

for i = 1:length(prob)

p = prob(i);
rhoABC = v_w*v_w';
rhoABC = (1-p)*rhoABC + p*rAB;
rhoABC = rhoABC/trace(rhoABC);

cmi_rho = CMI(rhoABC);

rhoAB = PartialTrace(rhoABC, 3, [da db dc]);


%% HPTP recover rho 

JI = MaxEntangled(2,0,1)*MaxEntangled(2,0,1)'*eye(4)*da;  % identity map
cvx_begin sdp quiet
cvx_solver sedumi
cvx_precision best
    variable JD1(db*db*dc, db*db*dc) hermitian
    variable JD2(db*db*dc, db*db*dc) hermitian
    variable p1
    variable p2
    JD = JD1 - JD2;
    cost = p1 + p2;

    JN = PermuteSystems(kron(JI, JD), [1 3 2 4], [da da db db*dc]);
    sigABC = ApplyMap(rhoAB, JN);
  
    minimize cost
    subject to
        JD1 >= 0; 
        PartialTrace(JD1, 2, [db db*dc]) == p1 * eye(db);
        JD2 >= 0; 
        PartialTrace(JD2, 2, [db db*dc]) == p2 * eye(db);
        sigABC == rhoABC; 

cvx_end

cost_list(i) = log2(cost);
cmi_list(i) = cmi_rho;

end
yyaxis left
scatter(prob, cost_list, 'x')
ylabel('$\nu(\rho_{ABC})$', 'interpreter','latex')

yyaxis right
scatter(prob, cmi_list)
ylabel('CMI', 'interpreter','latex')
xlabel('$p$','interpreter','latex')




