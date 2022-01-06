function [adjSUP,L] = SUP_adjustment_diag(SUP, B1_sub, Mask_sub, B1_ref, Mask_ref)

% pulse adjustment routine
% Adjustment matrix is forced here to be diagonal
% author : Vincent Gras
% contact : vincent.gras@cea.fr

Nc = size(B1_ref, 4);

X = map2vect(Mask_sub>0 & Mask_ref>0, B1_ref);
Y = map2vect(Mask_sub>0 & Mask_ref>0, B1_sub);
fun = @(X,Y,D) norm(Y - X*diag(D), 'fro')^2 ;

x0 = [ones(Nc,1); zeros(Nc, 1)];
A = [];
B = [];
Aeq = [];
Beq = [];
lb = [];
ub = [];
nonlcon = [];

options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'final', 'MaxFunctionEvaluations', 3e+05, 'MaxIterations', 2000 );
sol = fmincon(@(s) fun(X,Y,s(1:Nc)+1i*s(1+Nc:2*Nc)),x0,A,B,Aeq,Beq,lb,ub,nonlcon,options);

L = diag(sol(1:Nc)+1i*sol(1+Nc:2*Nc));
adjSUP = SUP.new();
adjSUP.dRF = L\SUP.dRF;


