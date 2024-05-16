%close all
clear all
%clf reset

Nmax = 5;   % -Nmax to Nmax Fourier modes in y
Mmax = 5;   % -Mmax to Mmax Fourier modes in x
Ntot = 2*Nmax+1; 
Mtot = 2*Mmax+1;
kx = 0.5; % Values of kx, ky; do not change.
ky = 1;
nk = 10;  % number of kz increments; set to 1 for plotting individual modes
kzmin = 0.0; %Min value of kz
kzmax = 10.0; %Max value of kz
maxn = 5*Mtot*Ntot;
lambda = zeros(nk,1);
position = zeros(nk,1);

Re = 10000;
Pe = 10000;
Ri = 30;

%Amplitude of u tilde 
A = 0.5;
%Coefficient v0 
v0 = -0.9145;

% % This section of code plots the background flow if needed
% Set nk = 1 and use any value of k_z to use it. 
% nx = 100;  % Number of points in x, y directions for plotting
% ny = 50;
% ubackx = zeros(nx,ny);
% ubacky = zeros(nx,ny);
% streamfunc = zeros(nx,ny); % u_x = -dpsi/dy, u_y = dpsi/dx
% ubackamp = zeros(nx,ny);
% dx = 4*pi/nx;
% dy = 2*pi/ny;
% for j=1:ny
%     y(j) = j*dy;
%     for i=1:nx
%         x(i)= i*dx;
%         ubackx(i,j) = sin(y(j)) + 2*A*cos(ky*y(j))*cos(kx*x(i));
%         ubacky(i,j) = A*(v0*cos(kx*x(i))+sin(ky*y(j))*sin(kx*x(i)));
%         ubackamp(i,j) = sqrt(ubackx(i,j)^2+ubacky(i,j)^2);
%         streamfunc(i,j) = cos(y(j)) - 2*A*sin(ky*y(j))*cos(kx*x(i)) + 2*A*v0*sin(kx*x(i));
%     end
% end
% figure(1)
% imagesc(x,y,transpose(ubackx),[-1 1]);
% colorbar;
% hold on
% contour(x,y,transpose(streamfunc),'-k');
% c.LineWidth=3;
% xlabel('x','Interpreter','latex');
% ylabel('y','Interpreter','latex');
% %imagesc(transpose(ubackx));
% %imagesc(transpose(ubacky));
% set(gcf,'color','w')
% saveas(gcf,'Background0.2.pdf')


for j=1:nk
    %kz(1) = 0; % This is to use when nk = 1 
    kz(j) = kzmin+(j-1)*(kzmax-kzmin)/(nk-1)
    Amat = zeros(maxn,maxn);  % zero A matrix
    Bmat = zeros(maxn,maxn);  % zero B matrix
    realvalues = zeros(maxn,maxn);
    D = zeros(maxn,maxn);
   
    for n = (-Nmax):Nmax    % n scans y
        for m = (-Mmax):Mmax    % m scans x
            
            % Set to 0 to remove viscous/dissipative terms.
            term = (kx^2*m^2 +ky^2*n^2+kz(j)^2); 
            %term = 0.0;

            mind = Mmax + m + 1; %Necessary reshuffling because 
            nind = Nmax + n + 1; %Matlab does not do negative indices
            
            ind = indexing(mind,nind,Mmax,Nmax);

            %u-equation 
            if(n>-Nmax) % q_(n-1, m)
            Amat(ind.u(mind,nind),ind.u(mind,nind-1)) = complex(-m*kx/2, 0); %
            Amat(ind.u(mind,nind),ind.v(mind,nind-1)) = complex(-ky/2, 0); %
            end
            if(n<Nmax) % q_(n+1, m)
            Amat(ind.u(mind,nind),ind.u(mind,nind+1)) = complex(m*kx/2, 0); %
            Amat(ind.u(mind,nind),ind.v(mind,nind+1)) = complex(-ky/2, 0); %
            end
            if(m>-Mmax && n>-Nmax) % q_(n-1, m-1)
            Amat(ind.u(mind,nind),ind.u(mind-1,nind-1)) = complex(0, -m*kx*A/2 + (n-1)*ky*A/4); %
            Amat(ind.u(mind,nind),ind.v(mind-1,nind-1)) = complex(0,-A*ky/2); %
            end
            if(m>-Mmax && n<Nmax) % q_(n+1, m-1)
            Amat(ind.u(mind,nind),ind.u(mind-1,nind+1)) = complex(0, -m*kx*A/2 - (n+1)*ky*A/4); %
            Amat(ind.u(mind,nind),ind.v(mind-1,nind+1)) = complex(0,A*ky/2); %
            end
            if(m<Mmax && n>-Nmax) % q_(n-1, m+1)
            Amat(ind.u(mind,nind),ind.u(mind+1,nind-1)) = complex(0, -m*kx*A/2 - (n-1)*ky*A/4);%
            Amat(ind.u(mind,nind),ind.v(mind+1,nind-1)) = complex(0,-A*ky/2);%
            end
            if(m<Mmax && n<Nmax) % q_(n+1, m+1)           
            Amat(ind.u(mind,nind),ind.u(mind+1,nind+1)) = complex(0, -m*kx*A/2 + (n+1)*ky*A/4); %
            Amat(ind.u(mind,nind),ind.v(mind+1,nind+1)) = complex(0,A*ky/2); %
            end
            if(m>-Mmax) % q_(n, m-1)            
            Amat(ind.u(mind,nind),ind.u(mind-1,nind)) = complex(0,-n*v0*ky*A/2);%
            end
            if(m<Mmax) % q_(n, m+1)
            Amat(ind.u(mind,nind),ind.u(mind+1,nind)) = complex(0,-n*v0*ky*A/2);%
            end
            Amat(ind.u(mind,nind),ind.p(mind,nind))= -m*kx; %
            Amat(ind.u(mind,nind),ind.u(mind,nind))= -term/Re; %
            Bmat(ind.u(mind,nind),ind.u(mind,nind))= 1.0; %

            %v  equation   
            if(n>-Nmax) % q_(n-1, m)
            Amat(ind.v(mind,nind),ind.v(mind,nind-1)) = complex(-m*kx/2, 0); %
            end
            if(n<Nmax) % q_(n+1, m)
            Amat(ind.v(mind,nind),ind.v(mind,nind+1)) = complex(m*kx/2, 0); %
            end
            if(m>-Mmax && n>-Nmax) % q_(n-1, m-1)
            Amat(ind.v(mind,nind),ind.v(mind-1,nind-1)) = complex(0, -(m-1)*kx*A/2 + n*ky*A/4); %
            Amat(ind.v(mind,nind),ind.u(mind-1,nind-1)) = complex(0,A*kx/4); %
            end
            if(m>-Mmax && n<Nmax) % q_(n+1, m-1)
            Amat(ind.v(mind,nind),ind.v(mind-1,nind+1)) = complex(0, -(m-1)*kx*A/2 - n*ky*A/4); %
            Amat(ind.v(mind,nind),ind.u(mind-1,nind+1)) = complex(0,-A*kx/4); %
            end
            if(m<Mmax && n>-Nmax) % q_(n-1, m+1)
            Amat(ind.v(mind,nind),ind.v(mind+1,nind-1)) = complex(0, -(m+1)*kx*A/2 - n*ky*A/4); %
            Amat(ind.v(mind,nind),ind.u(mind+1,nind-1)) = complex(0,A*kx/4); %
            end
            if(m<Mmax && n<Nmax) % q_(n+1, m+1)
            Amat(ind.v(mind,nind),ind.v(mind+1,nind+1)) = complex(0, -(m+1)*kx*A/2 + n*ky*A/4); %
            Amat(ind.v(mind,nind),ind.u(mind+1,nind+1)) = complex(0,-A*kx/4); %
            end
            if(m>-Mmax) % q_(n, m-1)
            Amat(ind.v(mind,nind),ind.v(mind-1,nind)) = complex(0,-n*v0*ky*A/2); %
            Amat(ind.v(mind,nind),ind.u(mind-1,nind)) = complex(0,-v0*kx*A/2); %
            end
            if(m<Mmax) % q_(n, m+1)
            Amat(ind.v(mind,nind),ind.v(mind+1,nind)) = complex(0,-n*v0*ky*A/2); %
            Amat(ind.v(mind,nind),ind.u(mind+1,nind)) = complex(0,v0*kx*A/2); % 
            end
            Amat(ind.v(mind,nind),ind.p(mind,nind))= -n*ky; %
            Amat(ind.v(mind,nind),ind.v(mind,nind))= -term/Re; %
            Bmat(ind.v(mind,nind),ind.v(mind,nind))= 1.0; %


            %w  equation   
            if(n>-Nmax) % q_(n-1, m)
            Amat(ind.w(mind,nind),ind.w(mind,nind-1)) = complex(-m*kx/2,0); %
            end
            if(n<Nmax) % q_(n+1, m)
            Amat(ind.w(mind,nind),ind.w(mind,nind+1)) = complex(m*kx/2,0); %
            end
            if(m>-Mmax && n>-Nmax) % q_(n-1, m-1)
            Amat(ind.w(mind,nind),ind.w(mind-1,nind-1)) = complex(0, -(m-1)*kx*A/2 + (n-1)*ky*A/4); %
            end
            if(m>-Mmax && n<Nmax) % q_(n+1, m-1)
            Amat(ind.w(mind,nind),ind.w(mind-1,nind+1)) = complex(0, -(m-1)*kx*A/2 - (n+1)*ky*A/4); %
            end
            if(m<Mmax && n>-Nmax) % q_(n-1, m+1)
            Amat(ind.w(mind,nind),ind.w(mind+1,nind-1)) = complex(0, -(m+1)*kx*A/2 - (n-1)*ky*A/4); %
            end
            if(m<Mmax && n<Nmax) % q_(n+1, m+1)
            Amat(ind.w(mind,nind),ind.w(mind+1,nind+1)) = complex(0, -(m+1)*kx*A/2 + (n+1)*ky*A/4); %
            end
            if(m>-Mmax) % q_(n, m-1)
            Amat(ind.w(mind,nind),ind.w(mind-1,nind)) = complex(0,-n*v0*ky*A/2); %
            end
            if(m<Mmax) % q_(n, m+1)
            Amat(ind.w(mind,nind),ind.w(mind+1,nind)) = complex(0,-n*v0*ky*A/2); %
            end
            Amat(ind.w(mind,nind),ind.t(mind,nind))= sqrt(Ri); %Ri; %
            Amat(ind.w(mind,nind),ind.p(mind,nind))= -kz(j); %
            Amat(ind.w(mind,nind),ind.w(mind,nind))= -term/Re; %
            Bmat(ind.w(mind,nind),ind.w(mind,nind))= 1.0; %

            %Temperature equation   
            if(n>-Nmax) % q_(n-1, m)
            Amat(ind.t(mind,nind),ind.t(mind,nind-1)) = complex(-m*kx/2, 0); %
            end
            if(n<Nmax) % q_(n+1, m)
            Amat(ind.t(mind,nind),ind.t(mind,nind+1)) = complex(m*kx/2, 0); %
            end
            if(m>-Mmax && n>-Nmax) % q_(n-1, m-1)
            Amat(ind.t(mind,nind),ind.t(mind-1,nind-1)) = complex(0, -(m-1)*kx*A/2 + (n-1)*ky*A/4); %
            end
            if(m>-Mmax && n<Nmax) % q_(n+1, m-1)
            Amat(ind.t(mind,nind),ind.t(mind-1,nind+1)) = complex(0, -(m-1)*kx*A/2 - (n+1)*ky*A/4); %
            end
            if(m<Mmax && n>-Nmax) % q_(n-1, m+1)
            Amat(ind.t(mind,nind),ind.t(mind+1,nind-1)) = complex(0, -(m+1)*kx*A/2 - (n-1)*ky*A/4); %
            end
            if(m<Mmax && n<Nmax) % q_(n+1, m+1)
            Amat(ind.t(mind,nind),ind.t(mind+1,nind+1)) = complex(0, -(m+1)*kx*A/2 + (n+1)*ky*A/4); %
            end
            if(m>-Mmax) % q_(n-1, m)
            Amat(ind.t(mind,nind),ind.t(mind-1,nind)) = complex(0,-n*v0*ky*A/2); %
            end
            if(m<Mmax) % q_(n-1, m)
            Amat(ind.t(mind,nind),ind.t(mind+1,nind)) = complex(0,-n*v0*ky*A/2); %
            end
            Amat(ind.t(mind,nind),ind.w(mind,nind))= -1.0 %- sqrt(Ri); %-1.0;  %
            Amat(ind.t(mind,nind),ind.t(mind,nind))= -term/Pe; %
            Bmat(ind.t(mind,nind),ind.t(mind,nind))= 1.0; %
 
            %Cont equation
            Amat(ind.p(mind,nind),ind.u(mind,nind))= kx*m; %
            Amat(ind.p(mind,nind),ind.v(mind,nind))= ky*n; %
            Amat(ind.p(mind,nind),ind.w(mind,nind))= kz(j); %
            Bmat(ind.p(mind,nind),ind.p(mind,nind))= 0.0; %
        end
    end
    % figure(1)  % This plots Amat, which is useful for debugging.
    % Cmat = log(Amat.*conj(Amat)+0.1);
    % imagesc(Cmat);

    [V,D] = eig(Amat,Bmat,"qz");
    realvalues = real(diag(D));

    % This next bit helps ignore negative and Inf eigenvalues
    for q=1:maxn
        if(realvalues(q)<=0.)
         realvalues(q)= 0.;   
        end
        if(realvalues(q)>=1.e2)
           realvalues(q)= 0;
        end
    end
    % This saves the eigenvalue of correct mode 
    [lambda(j,1),position(j,1)] = max(realvalues); 
    % This saves its functional form; only use if needed, set nk = 1.
     % bestmode = V(:,position(j,1)); 
     ## uplot = makemode(bestmode,kz(j),Mmax,Nmax); % This is to plot it if needed.

end
% This saves lambda(kz) to file
fileID = fopen('lambda.dat','w');
fprintf(fileID,'# Re %d\n',Re);
fprintf(fileID,'# Pe %d\n',Pe);
fprintf(fileID,'# Ri %d\n',Ri);
fprintf(fileID,'# Nmax %d\n',Nmax);
fprintf(fileID,'# Mmax %d\n',Mmax);
fprintf(fileID,'# A %d\n',A);
fprintf(fileID,"# kz, lambda\n");
for i=1:nk
    fprintf(fileID,'%8.4f %8.4f\n',kz(i),lambda(i,1));
end
fclose(fileID);

% This plots lambda(kz)
figure(3)
plot(kz,lambda(:,1))
xlim([0 10]);
ylim([0 0.5]);
xlabel('$k_z$','Interpreter','latex')
ylabel('Re($\lambda$)','Interpreter','latex')
%title('Re($\lambda$) vs. $k_z$',...
%     'FontSize', 16, 'Interpreter','latex')
set(gcf,'color','w')
%grid minor
saveas(gcf,'Lamba.pdf')


function [ind] = indexing(mind,nind,Mmax,Nmax)

    Ntot = 2*Nmax+1;
    Mtot = 2*Mmax+1;

    ind.u = zeros(Mtot,Ntot);
    ind.v = zeros(Mtot,Ntot);
    ind.w = zeros(Mtot,Ntot);
    ind.t = zeros(Mtot,Ntot);
    ind.p = zeros(Mtot,Ntot);

    ind.u(mind,nind) = mind + Mtot*(nind-1);
    ind.v(mind,nind) = Mtot*Ntot+ind.u(mind,nind);
    ind.w(mind,nind) = Mtot*Ntot+ind.v(mind,nind);
    ind.t(mind,nind) = Mtot*Ntot+ind.w(mind,nind);
    ind.p(mind,nind) = Mtot*Ntot+ind.t(mind,nind);

    if(mind>1)
       ind.u(mind-1,nind) = ind.u(mind,nind)-1;
       ind.v(mind-1,nind) = ind.v(mind,nind)-1;
       ind.w(mind-1,nind) = ind.w(mind,nind)-1;
       ind.t(mind-1,nind) = ind.t(mind,nind)-1;
       ind.p(mind-1,nind) = ind.p(mind,nind)-1;
    end
    if(mind<Mtot)
       ind.u(mind+1,nind) = ind.u(mind,nind)+1;
       ind.v(mind+1,nind) = ind.v(mind,nind)+1;
       ind.w(mind+1,nind) = ind.w(mind,nind)+1;
       ind.t(mind+1,nind) = ind.t(mind,nind)+1;
       ind.p(mind+1,nind) = ind.p(mind,nind)+1;
    end
  if(nind>1)
       ind.u(mind,nind-1) = ind.u(mind,nind)-Mtot;
       ind.v(mind,nind-1) = ind.v(mind,nind)-Mtot;
       ind.w(mind,nind-1) = ind.w(mind,nind)-Mtot;
       ind.t(mind,nind-1) = ind.t(mind,nind)-Mtot;
       ind.p(mind,nind-1) = ind.p(mind,nind)-Mtot;
  end
    if(nind<Ntot)
       ind.u(mind,nind+1) = ind.u(mind,nind)+Mtot;
       ind.v(mind,nind+1) = ind.v(mind,nind)+Mtot;
       ind.w(mind,nind+1) = ind.w(mind,nind)+Mtot;
       ind.t(mind,nind+1) = ind.t(mind,nind)+Mtot;
       ind.p(mind,nind+1) = ind.p(mind,nind)+Mtot;
    end
    if(mind>1 && nind>1)
       ind.u(mind-1,nind-1) = ind.u(mind,nind)-Mtot-1;
       ind.v(mind-1,nind-1) = ind.v(mind,nind)-Mtot-1;
       ind.w(mind-1,nind-1) = ind.w(mind,nind)-Mtot-1;
       ind.t(mind-1,nind-1) = ind.t(mind,nind)-Mtot-1;
       ind.p(mind-1,nind-1) = ind.p(mind,nind)-Mtot-1;
    end       
    if(mind<Mtot && nind>1)
       ind.u(mind+1,nind-1) = ind.u(mind,nind)-Mtot+1;
       ind.v(mind+1,nind-1) = ind.v(mind,nind)-Mtot+1;
       ind.w(mind+1,nind-1) = ind.w(mind,nind)-Mtot+1;
       ind.t(mind+1,nind-1) = ind.t(mind,nind)-Mtot+1;
       ind.p(mind+1,nind-1) = ind.p(mind,nind)-Mtot+1;
    end       
    if(mind<Mtot && nind<Ntot)
       ind.u(mind+1,nind+1) = ind.u(mind,nind)+Mtot+1;
       ind.v(mind+1,nind+1) = ind.v(mind,nind)+Mtot+1;
       ind.w(mind+1,nind+1) = ind.w(mind,nind)+Mtot+1;
       ind.t(mind+1,nind+1) = ind.t(mind,nind)+Mtot+1;
       ind.p(mind+1,nind+1) = ind.p(mind,nind)+Mtot+1;
    end
    if(mind>1 && nind<Ntot)
       ind.u(mind-1,nind+1) = ind.u(mind,nind)+Mtot-1;
       ind.v(mind-1,nind+1) = ind.v(mind,nind)+Mtot-1;
       ind.w(mind-1,nind+1) = ind.w(mind,nind)+Mtot-1;
       ind.t(mind-1,nind+1) = ind.t(mind,nind)+Mtot-1;
       ind.p(mind-1,nind+1) = ind.p(mind,nind)+Mtot-1;
    end

    % ind.u
    % ind.v
    % ind.w
    % ind.p
    % ind.t
end


function [uphys] = makemode(vec,kz,Mmax,Nmax)

Ntot = 2*Nmax+1;
Mtot = 2*Mmax+1;
kx = 0.5;
ky = 1.;

nx = 100;
ny = 50;
nz = 50;

dx = 4*pi/nx;
dy = 2*pi/ny;
dz = 2*pi/nz;

uphys.x = zeros(nx,ny,nz);
uphys.y = zeros(nx,ny,nz);
uphys.z = zeros(nx,ny,nz);
uphys.vortx = zeros(nx,ny,nz);
uphys.vorty = zeros(nx,ny,nz);
uphys.vortz = zeros(nx,ny,nz);


for k=1:nz
   z =k*dz;
     for j=1:ny
        y = j*dy;
        for i=1:nx
            x = i*dx;
            
            sumx = 0.;
            sumy = 0.;
            sumz = 0.;
            sumvortx = 0.;
            sumvorty = 0.;
            sumvortz = 0.;
            for m=-Mmax:Mmax
                for n = -Nmax:Nmax
                    mind = Mmax + m + 1;  
                    nind = Nmax + n + 1;
                    umn = vec(mind + Mtot*(nind-1));
                    vmn = vec(Mtot*Ntot + mind + Mtot*(nind-1));
                    wmn = vec(2*Mtot*Ntot + mind + Mtot*(nind-1));
                    term = exp(complex(0,m*kx*x+n*ky*y));
                    sumx = sumx + umn*term;
                    sumy = sumy + vmn*term;
                    sumz = sumz + wmn*term;
                    sumvortx = sumvortx + (complex(0,n*ky)*wmn-complex(0,kz)*vmn)*term;
                    sumvorty = sumvorty + (complex(0,kz)*umn-complex(0,m*kx)*wmn)*term;
                    sumvortz = sumvortz + (complex(0,m*kx)*vmn-complex(0,n*ky)*umn)*term;
                end
            end
            term = exp(complex(0,kz*z));
            sumx = sumx*term;
            sumy = sumy*term;
            sumz = sumz*term;
            sumvortx = sumvortx*term;
            sumvorty = sumvorty*term;
            sumvortz = sumvortz*term;
            uphys.x(i,j,k) = real(sumx);
            uphys.y(i,j,k) = real(sumy);
            uphys.z(i,j,k) = real(sumz);
            uphys.vortx(i,j,k) = real(sumvortx);
            uphys.vorty(i,j,k) = real(sumvorty);
            uphys.vortz(i,j,k) = real(sumvortz);

        end
    end
end

end



  
