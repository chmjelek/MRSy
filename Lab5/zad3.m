clear, clc;
% zad 3
% Przykładowe dane
nx = 10; ny = 10; xa = 1; xb = 2; ya = 1; yb = 2;
f = @(x, y) x ./ y + y ./ x;
g = @(x, y) x .* y .* log(x .* y);
uxa = @(x) x .* log(x); %u1
uxb = @(x) x .* log(4 * x .^ 2); %u3
uya = @(y) y .* log(y); %u4
uyb = @(y) 2 * y .* log(2 * y); %u2

% Obliczenia
h = (xb-xa)/(nx+1);
k = (yb-ya)/(ny+1);
xM = linspace(xa + h, xb - h, nx);
yM  = linspace(ya+k, yb-k, ny);
T = -2*(h^2 + k^2) * eye(nx) + k^2 * diag(diag(eye(nx-1)),-1) + k^2 * diag(diag(eye(nx-1)),1);
I = eye(nx);
B = kron(eye(ny), T);
C = kron(diag(diag(eye(ny-1)),-1), h^2 * I);
D = kron(diag(diag(eye(ny-1)),1), h^2 * I);
A = B + C + D;

[X1 Y1] = meshgrid(xM,yM);
F = f(X1, Y1)' .* ones(nx, ny);
F = h^2 * k^2 .* F;
F(:, 1) = F(:, 1) - h^2 * uxa(xM)';
F(:, ny) = F(:, ny) - h^2 * uxb(xM)';
F(1, :) = F(1, :) - k^2 * uya(yM);
F(nx, :) = F(nx, :) - k^2 * uyb(yM);
F = reshape(F,  nx*ny, 1);
% Rozwiązanie
U = linsolve(A,F);
U = reshape(U,  nx, ny)';
[X,Y] = meshgrid(xa:h:xb, ya:k:yb);
U = [uya(yM)' .* diag(eye(ny)), U, uyb(yM)' .* diag(eye(ny))];
XM = [xa xM xb];
G = g(X,Y);
U = [uxa(XM)' .* diag(eye(nx+2)), U', uxb(XM)' .* diag(eye(nx+2))]';
% Błąd
Error = max(max(abs(G-U)));
% Wykres
subplot(1,2,1)
surf(X,Y,U)
subplot(1,2,2)
surf(X,Y,G)