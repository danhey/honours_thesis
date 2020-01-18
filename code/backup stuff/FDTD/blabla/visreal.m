function visreal(array2d, xrange, yrange)

%% Get the maximum magnitude before taking the real part.
% cmax = max(abs(array2d(:)));
cmax = max(real(array2d(:)));


%% Attach a row and column at the ends (needed from pcolor)
array2d = real(array2d);
[Nx, Ny] = size(array2d);
array2d = [array2d, array2d(:,1)];
array2d = [array2d; array2d(1,:)];

%% Create the matrices for x,y,C
xs = linspace(xrange(1), xrange(2), Nx+1);
ys = linspace(yrange(1), yrange(2), Ny+1);
[X, Y] = meshgrid(xs, ys);
C = permute(array2d, [2 1]);

%% Draw with pcolor().
h = pcolor(X, Y, C);

%% Aesthetics
set(h, 'EdgeColor', 'none');
set(gca, 'TickDir', 'out');
axis image;

caxis([-cmax, cmax]);
colormap('b2r')
colorbar;
