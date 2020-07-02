function [M,bc] = matrix_T_input(n,N,K,q,bc_val)

% j refers to index of radial position, starting from j = 0
j = 0:(n-1);
j = j';

top_diag_inds = 2:(n-2); % equations involve j = 1 to j = n-2
cen_diag_inds = 2:(n-1); % equations involve j = 1 to j = n-1
bot_diag_inds = 3:(n-1); % equations involve j = 2 to j = n-1

s = q./(3*N.*j);

top_diag = -s(top_diag_inds).*(j(top_diag_inds).*K(top_diag_inds) + j(top_diag_inds+1).*K(top_diag_inds+1));
cen_diag = 1+s(cen_diag_inds).*(2*j(cen_diag_inds).*K(cen_diag_inds) + j(cen_diag_inds+1).*K(cen_diag_inds+1) + j(cen_diag_inds-1).*K(cen_diag_inds-1));
bot_diag = -s(bot_diag_inds).*(j(bot_diag_inds).*K(bot_diag_inds) + j(bot_diag_inds-1).*K(bot_diag_inds-1));

M = diag(cen_diag)+diag(top_diag,1)+diag(bot_diag,-1);
M(1,1) = 1+s(2)*(K(2)+2*K(3)); % neumann, s(2) refers j = 1

bc = zeros(n-2,1);
bc(end) = s(end-1)*(j(end-1)*K(end-1)+j(end)*K(end))*bc_val;