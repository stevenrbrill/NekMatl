function plot1_new = post_channel(N,Ex,Ey,X,Y,Ys,en_on,time,u,psi_xy,plot1)

u_mean = zeros(Ey*length(u(1,:,1)),1);
for i = 1:Ey
    for j = 1:Ex
        for k = 1:N+1
            u_mean((i-1)*(N+1)+k) = u_mean((i-1)*(N+1)+k) + sum(u(:,k,(i-1)*Ex+j));
        end
    end
end
if en_on
    for i = 1:Ey
        u_mean((i-1)*(N+1)+1:(i)*(N+1)) = u_mean((i-1)*(N+1)+1:(i)*(N+1)) + (N+1)*Ex*psi_xy(1,:,(i-1)*Ex+1);
    end
end
u_mean = u_mean/((N+1)*Ex);

if en_on
    u_recon = u + psi_xy;
else
    u_recon = u;
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
    hold on
end
hold off
%         xlim([0,1])
ylim([-1,1])
xlabel('u')
ylabel('y')
%         pause(.1);
umax = glmax(u_recon);
[time umax 1/(2*pi/(umax*time))]


end