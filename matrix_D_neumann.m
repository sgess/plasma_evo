function M = matrix_D_neumann(n,D,q)

% j refers to index of radial position, starting from j = 0
j = 0:(n-1);
j = j';

top_diag_inds = 2:(n-2); % equations involve j = 1 to j = n-2
cen_diag_inds = 2:(n-1); % equations involve j = 1 to j = n-1
bot_diag_inds = 3:(n-1); % equations involve j = 2 to j = n-1

s = q./(2*j);

top_diag = -s(top_diag_inds).*(j(top_diag_inds).*D(top_diag_inds) + j(top_diag_inds+1).*D(top_diag_inds+1));
cen_diag = 1+s(cen_diag_inds).*(2*j(cen_diag_inds).*D(cen_diag_inds) + j(cen_diag_inds+1).*D(cen_diag_inds+1) + j(cen_diag_inds-1).*D(cen_diag_inds-1));
bot_diag = -s(bot_diag_inds).*(j(bot_diag_inds).*D(bot_diag_inds) + j(bot_diag_inds-1).*D(bot_diag_inds-1));

M = diag(cen_diag)+diag(top_diag,1)+diag(bot_diag,-1);
M(1,1) = 1+s(2)*(D(2)+2*D(3)); % neumann, s(2) refers to D at j = 1
M(end,end) = 1+s(end-1)*(j(end-1)*D(end-1)+j(end-2)*D(end-2)); % neumann, s(2) refers j = 1

%bc = zeros(n-2,1);
%bc(end) = s(end-1)*(j(end-1)*D(end-1)+j(end)*D(end))*bc_val;