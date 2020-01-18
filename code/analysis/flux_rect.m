function P = flux_rect(L0, Sx, Sy, ind_x, ind_y, dL)
%% Input Parameters
% Sx, Sy: 2D arrays of Poynting vectors in x- and y-directions
% ind_x, ind_y: x- and y-indices of rectangular region
% dL: [dx dy]

dL = dL*L0; 

%% Output Parameters
% P: net power flux going out of rectangular region
Pxn = -1*dL(2) * sum(Sx(ind_x(1), ind_y(1) : ind_y(2))); 

Pxp = 1*dL(2) * sum(Sx(ind_x(2)+1, ind_y(1) : ind_y(2))); 

Pyp = 1*dL(1) * sum(Sy(ind_x(1) : ind_x(2), ind_y(2)+1)); 
Pyn = -1*dL(1) * sum(Sy(ind_x(1) : ind_x(2), ind_y(1))); 

P = Pxn + Pxp + Pyn + Pyp; 

end

