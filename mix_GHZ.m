clear;

%% Required

% CVX: http://www.cvxr.com
% QETLAB: http://www.qetlab.com

%% case study
da = 2;
db = 2;
dc = 2;
v_w = [0 1 1 0 1 0 0 0]';
v_ghz = [1 0 0 0 0 0 0 1]';
rhoGHZ = v_ghz * v_ghz';
rhoGHZ = rhoGHZ / trace(rhoGHZ);

rhoW = v_w * v_w';
rhoW = rhoW / trace(rhoW);

p = 0.0:0.01:1.0;
cost_store = zeros(1, numel(p));
for j=1:numel(p)

    rhoABC = (1-p(j))*rhoW + p(j)*eye(2^3) / 2^3;
    
    JI = MaxEntangled(2,0,1)*MaxEntangled(2,0,1)'*eye(4)*da; % identity map
    
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
    sigABC = ApplyMap(PartialTrace(rhoABC, 3, [da, db, dc]), JN);
    
    minimize cost
    subject to
        JD1 >= 0; 
        PartialTrace(JD1, 2, [db db*dc]) == p1 * eye(db);
        JD2 >= 0; 
        PartialTrace(JD2, 2, [db db*dc]) == p2 * eye(db);
        sigABC == rhoABC; 
    cvx_end
    
    cost

    cost_store(j) = log2(cost);
end



