function latexImage(array2d, xrange, yrange, saver)

%% Get the maximum magnitude before taking the real part.
% cmax = max(abs(array2d(:)));
cmax = max(real(array2d(:)));


%% Attach a row and column at the ends.  (Not visualized, but required by pcolor().)
array2d = real(array2d);
[Nx, Ny] = size(array2d);
array2d = [array2d, array2d(:,1)];
array2d = [array2d; array2d(1,:)];

%% Create the matrices for x-locations (X), y-locations (Y), and color (C).
xs = linspace(xrange(1), xrange(2),13);
ys = linspace(yrange(1), yrange(2), Ny+1);
[X, Y] = meshgrid(xs, ys);
C = permute(array2d, [2 1]);

%% Draw with pcolor().
h = imagesc(C);

%% Make the figure look better.
set(gca, 'TickDir', 'out');
set(gca,'Ydir','Normal')
axis image;

xticklabels = xrange(1):2:xrange(2);
yticklabels = yrange(1):2:yrange(2);
xticks = linspace(1, Nx, numel(xticklabels));
yticks = linspace(1, Ny, numel(yticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)

xlabel('x (\mum)')
ylabel('y (\mum)')

caxis([-cmax, cmax]);
colormap('colormap_br')
colorbar;

if saver
    cleanfigure('targetResolution', 100);
    matlab2tikz('nameoffile2.tex','width','\figW','height','\figH');
end