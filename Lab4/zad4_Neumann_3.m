clear, clc;
%dane wej�ciowe
f = @(x) exp(x) .* (-x .^2 + x + 2);
g = @(x) (x) .* exp(x);
a = 0;
b = 1;
ax = @(x) 1;
bx = @(x) -x;
cx = @(x) 1;
ua = 0;
ubf = 2 * exp(1);
c = b-a;
n = 1000;
h = c/(n+1)
x = linspace((a+h),(b-h),n);
%obliczenia
a1 = [(-2 * ax(x) + cx(x) * h^2).* diag(eye(n)); 3];
A1 = diag(a1);
a2 = (ax(x) +1/2 * bx(x) * h)' .* diag(eye(n));
A2 = diag(a2, 1);
a3 = [(ax(x(2:n)) -1/2 * bx(x(2:n)) * h)' .* diag(eye(n-1)); -4];
A3 = diag(a3, -1);
a4 = [zeros(n-2,1); 1];
A4 = diag(a4, -2);

A = A1 + A2 + A3 + A4;

F = h^2 * f(x);
F(1) = F(1) - ua;
F = [F 2*h*ubf];
U = linsolve(A,F');
U = [ua U'];
X = [a x b];
%wykres
plot(X, U, X, g(X), 'ro');
plot(X, U, X, g(X), 'ro');
legend('Metoda Analityczna','Metoda Numeryczna');
%error
E = max(abs(g(X) - U))