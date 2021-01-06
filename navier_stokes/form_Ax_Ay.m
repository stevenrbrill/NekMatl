function [A_x,A_y] = form_Ax_Ay(N1,E,w1d,J,dpdx_dpdx,dpdy_dpdy,Re)

A_x = zeros(N1*N1*E);
A_y = zeros(N1*N1*E);
JRe = reshape(J.*1./Re,[N1*N1,E]);
for ie=1:E
    A_x(1+(ie-1)*N1*N1:ie*N1*N1,1+(ie-1)*N1*N1:ie*N1*N1) = reshape((w1d.*JRe(:,ie)')*dpdx_dpdx(:,:,ie),[N1*N1,N1*N1]);
    A_y(1+(ie-1)*N1*N1:ie*N1*N1,1+(ie-1)*N1*N1:ie*N1*N1) = reshape((w1d.*JRe(:,ie)')*dpdy_dpdy(:,:,ie),[N1*N1,N1*N1]);
end