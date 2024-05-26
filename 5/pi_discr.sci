function [Res]=Pade(delay, order)
s = poly(0,'s');
Res = (-delay*s + 2*order)^order / (delay*s +2*order)^order;
endfunction
n = 6;
T0 = 0.93;
K = 0.45;
Ti = 9;
T = 1.5;
Td = 0;
Ts = Td/8;
h = T/100;
S = poly(0, 's');
Wobj = Pade(T, 5)/(1+S*T0)^n
W1 = (1 + 1/(Ti*S) + Td*S/(1 + Ts*S))*K*Wobj
W = W1/(1 + W1)
sl = syslin('c', W);
dicrMat = dscr(sl, h);
Ad = dicrMat.A;
Ed = -eye(Ad);
Hd = lyap(Ad, Ed, 'd');
ld = spec(Hd)
if ld > 0 then
kd = norm(Hd, 2);
else
kd = %inf
end
printf("%4.4f\n", kd);
