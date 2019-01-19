clear, clc;
%dane wej�ciowe
f = @(x) -4*x;
g = @(x) exp(2) * (exp(4) - 1 )^-1 * (exp(2*x) - exp(-2*x)) + x;
a = 0;
b = 1;
ax = @(x) 1;
bx = @(x) 0;
cx = @(x) -4;
ua = 0;
ub = 2;
n = 750;
h = (b-a)/(n+1)
x = linspace((a+h),(b-h),n);
%obliczenia
a1 = (-2 .* ax(x) + cx(x) .*h^2) .* diag(eye(n));
A1 = diag(a1);
a2 = (ax(x(2:n)) -1/2 .* bx(x(2:n)) .*h)' .* diag(eye(n-1));
A2 = diag(a2, -1);
a3 = (ax(x(1:n-1)) +1/2 .* bx(x(1:n-1)) .*h)' .* diag(eye(n-1));
A3 = diag(a3, 1);
A = A1  + A2 + A3;
F = (h^2 * f(x));
F(1) = F(1) - ua * (ax(1) - 1/2 * bx(1)*h);
F(n) = F(n) - ub * (ax(n) + 1/2 * bx(n)*h);
U = linsolve(A,F');
U = [ua U' ub];
X = [a x b];
%wykres
plot(X, U, X, g(X), 'ro');
legend('Metoda Analityczna','Metoda Numeryczna',0);
%error
E = max(abs(g(X) - U))