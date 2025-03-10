clear;

%% Required

% CVX: http://www.cvxr.com
% QETLAB: http://www.qetlab.com

%% setup

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

fCPTP = zeros(1, numel(p));
fHPTP = zeros(1, numel(p));
for j=1:numel(p)

    rhoABC = (1-p(j))*rhoGHZ + p(j)*eye(2^3) / 2^3;
    
    JI = MaxEntangled(2,0,1)*MaxEntangled(2,0,1)'*eye(4)*da; % identity map
    
    %% CPTP
    cvx_begin sdp quiet
    cvx_solver SDPT3
    cvx_precision best
    variable JD(db*db*dc, db*db*dc) hermitian
    
    JN = PermuteSystems(kron(JI, JD), [1 3 2 4], [da da db db*dc]);
    sigABC = ApplyMap(PartialTrace(rhoABC, 3, [da, db, dc]), JN);
    
    dcptp = SchattenNorm(sigABC - rhoABC, 1);
    minimize dcptp
    subject to
        JD >= 0; 
        PartialTrace(JD, 2, [db db*dc]) == eye(db);
    cvx_end
    
    fCPTP(j) = dcptp;


    %% HPTP
    cvx_begin sdp quiet
    cvx_solver SDPT3
    cvx_precision best
    variable JD(db*db*dc, db*db*dc) hermitian
    
    JN = PermuteSystems(kron(JI, JD), [1 3 2 4], [da da db db*dc]);
    sigABC = ApplyMap(PartialTrace(rhoABC, 3, [da, db, dc]), JN);
    
    dhptp = SchattenNorm(sigABC - rhoABC, 1);
    minimize dhptp
    subject to
        PartialTrace(JD, 2, [db db*dc]) == eye(db);
    cvx_end
    
    fHPTP(j) = dhptp;

end

%% plot
c1 = [20/255 54/255 95/255];
c2 = [118/255 162/255 185/255];
c3 = [191/255 217/255 229/255];
c4 = [214/255 79/255 56/255];

grid on
hold on

plot(p, fCPTP, '-', 'LineWidth',  1.5, 'Color', c1);
plot(p, fHPTP, '-', 'LineWidth',  1.5, 'Color', c4);


leg = legend('CPTP', 'HPTP');
set(leg,'Interpreter','latex','FontSize',16,'Location','northwest');
xlabel('Noise rate $p$', 'Interpreter','latex', 'FontSize', 20, 'FontName','Times New Roman')
ylabel('$\varepsilon$', 'Interpreter','latex','FontSize',20, 'FontName','Times New Roman')

x_fill = [p, fliplr(p)];
inBetween = [fCPTP, fliplr(fHPTP)];

fill(x_fill, inBetween, c2, 'FaceAlpha', 0.15, 'LineStyle', 'none', 'Marker','none');