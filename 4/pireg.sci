n = 6;
T0 = 0.93;
K = 0.531;
Ti = 4.933;
T = 1;
Td = Ti;
Ts = Td;
h = T/10;
S = poly(0, 's');

Wobj = 2*(1-T*S + (T*S)^2/2 - (T*S)^3/6)/(1+S*T0)^n;
W1 = (1 + 1/(Ti*S))*K*Wobj;

W = W1/(1 + W1);
sl = syslin('c', W);
dicrMat = dscr(sl, h);
t = [0:h:100];
v = zeros(dicrMat.B);
u = ones(t);
x = zeros(u);

for i=1:length(u)
  v = dicrMat.A * v + dicrMat.B;
  x(i) = dicrMat.C * v + dicrMat.D;
end
plot(t, x, 'blue');
y0 = [
 evstr(csvRead('./pireg_T1_3.csv', ",", [], "string"));
]';
y = y0;
err = sum((x-y')*(x'-y))/(length(t));
disp("Err = ", err);

t_exact = [0:1/10:100];
y_exact = [
 evstr(csvRead('./pireg_T1.csv', ",", [], "string"));
]';
title("pireg_h")
plot(t_exact,y_exact,'red');
