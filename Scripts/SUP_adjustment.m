function [adjSUP,L] = SUP_adjustment(SUP, B1_sub, Mask_sub, B1_ref, Mask_ref)

% pulse adjustment routine
% author : Vincent Gras
% contact : vincent.gras@cea.fr

Nc = size(B1_ref, 4);

X = map2vect(Mask_sub>0 & Mask_ref>0, B1_ref);
Y = map2vect(Mask_sub>0 & Mask_ref>0, B1_sub);
fun = @(X,Y,L) norm(Y - X*L, 'fro')^2 ;

x0 = [reshape(eye(Nc), Nc*Nc, 1); zeros(Nc*Nc, 1)];
A = [];
B = [];
Aeq = [];
Beq = [];
lb = [];
ub = [];
nonlcon = [];

options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'final', 'MaxFunctionEvaluations', 3e+05, 'MaxIterations', 2000 );
sol = fmincon(@(s) fun(X,Y,reshape(s(1:Nc*Nc)+1i*s(1+Nc*Nc:2*Nc*Nc),Nc,Nc)),x0,A,B,Aeq,Beq,lb,ub,nonlcon,options);

L = reshape(sol(1:Nc*Nc)+1i*sol(1+Nc*Nc:2*Nc*Nc),Nc,Nc);
adjSUP = SUP.new();
adjSUP.dRF = L\SUP.dRF;


