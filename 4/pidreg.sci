n = 6;
T0 = 0.93;
K = 0.531;
Ti = 4.933;
T = 1;
Td = Ti/4;
Ts = Td/8;
h = 1.5;
function [Res]=Pade(delay, order)
 s = poly(0, 's');
 Res = (-delay*s + 2*order)^order / (delay*s + 2*order)^order;
endfunction
S = poly(0, 's');

Wobj = Pade(T, n)/(1+S*T0)^n
W1 = (1 + 1/(Ti*S) + Td*S/(1+Ts*S))*K*Wobj
W = W1/(1 + W1);
sl = syslin('c', W);
dicrMat = dscr(sl, h);
t = [0:h:100];
v = zeros(dicrMat.B);
u = ones(t);
x = zeros(u);

for i=1:length(u)
  v = dicrMat.A * v + dicrMat.B * u(i);
  x(i) = dicrMat.C * v + dicrMat.D * u(i);
end
plot(t, x, 'blue');
y0 = [
 evstr(csvRead('./pidreg_T1_1.csv', ",", [], "string"));
]';
y = y0;
err = sum((x-y')*(x'-y))/(length(t));
disp("Err = ", err);
title("pidreg")
plot(t,y,'red');
