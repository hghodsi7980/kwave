figure
tiledlayout(2, 3)

% DEMOSPL
nexttile
hold on
x = [0.05, 0.6; 0.45, 0.5; 0.2, 0.2];
y = bsplpcw(x);
z = bsplpcw(x, []);
a = bspl(x);
b = bspl(y);
c = bspl(z);
x = [x; x(1, :)];
plot(x(:, 1), x(:, 2), 'o', 'Color', 'k')
plot(a(:, 1), a(:, 2), 'Color', '#EDB120')
plot(b(:, 1), b(:, 2), 'Color', '#D95319')
plot(c(:, 1), c(:, 2), 'Color', '#0072BD')
xlim([-0.05, 0.5])
pbaspect([1, 1, 1])
daspect([1, 1, 1])
bplt
axis off

% DEMOINT
nexttile
hold on
A = [0, 0.2; 1.85, 0.85; -0.85, 0.85; 1, 0.2];
B = [0, 0.8; 1.85, 0.15; -0.85, 0.15; 1, 0.8];
a = bspl(A, 1e3);
b = bspl(B, 1e3);
plot(a(:, 1), a(:, 2))
plot(b(:, 1), b(:, 2))
[~, val] = bsplint(A, B);
plot(val(:, 1), val(:, 2), 'o', 'Color', '#EDB120')
[~, val] = bsplint(A);
plot(val(:, 1), val(:, 2), 'o', 'Color', '#EDB120')
[~, val] = bsplint(B);
plot(val(:, 1), val(:, 2), 'o', 'Color', '#EDB120')
ylim([0.2, 0.8])
pbaspect([1 1 1])
daspect([1 1 1])
bplt
axis off

% DEMOFIT
nexttile
hold on
X = [0.59, 0.62; 0.44, 0.47; 0.49, 0.39; 0.7, 0.4; 0.9, 0.37];
A = bsplget(X);
B = bsplget(X, 0);
C = bsplget(X, 1);
D = bsplget(X, 1e3);
E = bsplfit(X, size(X, 1) - 2, 1);
a = bspl(A);
b = bspl(B);
c = bspl(C);
d = bspl(D);
e = bspl(E);
plot(X(:, 1), X(:, 2), 'o', 'Color', 'k')
plot(a(:, 1), a(:, 2), 'Color', '#0072BD')
plot(b(:, 1), b(:, 2), 'Color', '#D95319')
plot(c(:, 1), c(:, 2), 'Color', '#EDB120')
plot(d(:, 1), d(:, 2), 'Color', '#7E2F8E')
plot(e(:, 1), e(:, 2), 'Color', '#77AC30')
pbaspect([1, 1, 1])
daspect([1, 1, 1])
bplt
axis off

% DEMOBOX
nexttile
hold on
x = [0.42, 0.58, 0.47; 0.07, 0.83, 0.25; 0.22, 0.81, 0.67; ...
     0.79, 0.87 ,0.62; 0.28, 0.01, 0.61; 0.40, 0.75, 0.24];
E = bspl(x);
F = bsplbox(x, [], 1);
D = bsplbox(x, [], 0);
plot3(E(:, 1), E(:, 2), E(:, 3), 'Color', '#0072BD')
plot3(F([1, 2, 3, 4, 1], 1), F([1, 2, 3, 4, 1], 2), ...
      F([1, 2, 3, 4, 1], 3), 'Color', '#D95319')
plot3(F([5, 6, 7, 8, 5], 1), F([5, 6, 7, 8, 5], 2), ...
      F([5, 6, 7, 8, 5], 3), 'Color', '#D95319')
plot3(F([1, 2, 6, 5, 1], 1), F([1, 2, 6, 5, 1], 2), ...
      F([1, 2, 6, 5, 1], 3), 'Color', '#D95319')
plot3(F([4, 3, 7, 8, 4], 1), F([4, 3, 7, 8, 4], 2), ...
      F([4, 3, 7, 8, 4], 3), 'Color', '#D95319')
plot3(D([1, 2, 3, 4, 1], 1), D([1, 2, 3, 4, 1], 2), ...
      D([1, 2, 3, 4, 1], 3), 'Color', '#EDB120')
plot3(D([5, 6, 7, 8, 5], 1), D([5, 6, 7, 8, 5], 2), ...
      D([5, 6, 7, 8, 5], 3), 'Color', '#EDB120')
plot3(D([1, 2, 6, 5, 1], 1), D([1, 2, 6, 5, 1], 2), ...
      D([1, 2, 6, 5, 1], 3), 'Color', '#EDB120')
plot3(D([4, 3, 7, 8, 4], 1), D([4, 3, 7, 8, 4], 2), ...
      D([4, 3, 7, 8, 4], 3), 'Color', '#EDB120')
view(3)
pbaspect([1, 1, 1])
daspect([1, 1, 1])
bplt
axis off

% DEMOCRV
nexttile
hold on
A = [0.7, 0.75; 0, 0.5; 0.2, 0.8; 0.1, 0.4; 1, 0; 0, 0.25];
a = bspl(A, 1e3);
k = bsplcrv(A, 1e3);
[~, n] = bsplvec(A, 1e3);
k = a + 0.1 * range(a) .* n .* k / max(abs(k));
[b, idx, val] = bsplcsp(A);
plot(a(:, 1), a(:, 2))
plot(k(:, 1), k(:, 2))
for i = 1 : length(val)
    csp = val(i) > abs(bsplcrv(A, idx(i) + 0.02 * (rand - 0.5)));
    if csp
        plot(b(i, 1), b(i, 2), 'o', 'Color', '#7E2F8E')
    else
        plot(b(i, 1), b(i, 2), 'o', 'Color', '#EDB120')
    end
end
pbaspect([1, 1, 1])
daspect([1, 1, 1])
bplt
axis off

% DEMOVEC
nexttile
hold on
A = [0.80, 0.50, 0.46; 0.9, 0.95, 0.27; 0.43, 0.05, 0; 0.14, 0.65, 0.1; ...
     0.62, 0.70, 0.35];
x = bspl(A, 101);
[t, n] = bsplvec(A, 101);
t = x + 0.05 * t * range(x(:));
n = x + 0.05 * n * range(x(:));
plot3(x(:, 1), x(:, 2), x(:, 3))
for i = 1 : 2 : 101
    plot3([x(i, 1), t(i, 1)], [x(i, 2), t(i, 2)], [x(i, 3), t(i, 3)], ...
         'Color', '#D95319')
    plot3([x(i, 1), n(i, 1)], [x(i, 2), n(i, 2)], [x(i, 3), n(i, 3)], ...
         'Color', '#EDB120')
    plot3(t(i, 1), t(i, 2), t(i, 3), '.', 'Color', '#D95319', 'MarkerSize', 4)
    plot3(n(i, 1), n(i, 2), n(i, 3), '.', 'Color', '#EDB120', 'MarkerSize', 4)
end
pbaspect([1, 1, 1])
daspect([1, 1, 1])
bplt
axis off


% Format figure
function bplt()
    h = findall(gca, 'Type', 'axes');
    set(h, 'box', 'on', 'Color', 'w')
    set(h, 'TickDir', 'in')
    set(h, 'TickLabelInterpreter', 'LaTeX')
end