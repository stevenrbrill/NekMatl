function [plot1_new,u_mean] = post_channel(N,Ex,Ey,w,X,Y,Ys,en_on,step,time,u,psi_xy,N_en_y,plot1)

u_mean = zeros(Ey*length(u(1,:,1)),1);
for i = 1:Ey
    for j = 1:Ex
        for k = 1:N+1
            u_mean((i-1)*(N+1)+k) = u_mean((i-1)*(N+1)+k) + sum(w.*u(:,k,(i-1)*Ex+j));
        end
    end
end
if en_on
    for i = 1:Ey
        if ((i <= N_en_y) || (i > Ey-N_en_y))
            u_mean((i-1)*(N+1)+1:(i)*(N+1)) = u_mean((i-1)*(N+1)+1:(i)*(N+1)) + (Ex*2)*psi_xy{1}(1,:,(i-1)*Ex+1)';
        end
    end
end
u_mean = u_mean/(2*Ex);

u_recon = u;
if en_on
    for iy = 1:Ey
        for ix = 1:Ex
            if ((iy <= N_en_y) || (iy > Ey-N_en_y))
                i = (iy-1)*Ex+ix;
                u_recon(:,:,i) = u(:,:,i) + psi_xy{1}(:,:,i);
            end
        end
    end
end

if plot1
    figure(1);
else
    figure(3);
end
plotit(u_recon,X,Y);
if plot1
    figure(2);
    plot1_new = 0;
else
    figure(4);
    plot1_new = 1;
end
for i = 1:Ey
    plot(u_mean,Ys,'b-o')
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    yticks(linspace(-1,1,Ey+1));
    hold on
end
hold off
%         xlim([0,1])
ylim([-1,1])
xlabel('u')
ylabel('y')
%         pause(.1);
umax = glmax(u_recon);
disp([time umax 1/(2*pi/(umax*time))])

end