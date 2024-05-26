function [Res]=Pade(delay, order)
s = poly(0,'s');
Res = (-delay*s + 2*order)^order / (delay*s +2*order)^order;
endfunction
n = 6;
T0 = 0.93;
K = 0.65;
Ti = 4.8;
T = 0;
Td = Ti/4;
Ts = Td/8;
h = T/100;
S = poly(0, 's');
Wobj = Pade(T, 5)/(1+S*T0)^n
W1 = (1 + 1/(Ti*S) + Td*S/(1 + Ts*S))*K*Wobj
W = W1/(1 + W1)
sl = syslin('c', W);
[A B C D] = abcd(sl);
E = -eye(A);
H = lyap(A, E, 'c');
l = spec(H)
// проверяем устойчивость
if l > 0 then
k = norm(H, 2);
else
k = %inf
end
printf("%4.4f\n", k);
