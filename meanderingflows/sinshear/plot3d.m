

close all;
clear;
clc;

    datfileID = fopen('nk.dat', 'r');
    nk = fscanf(datfileID, '%d');

    %kxfileID = fopen('kx.dat','r');
    kzfileID = fopen('kz.dat','r');
    lfileID = fopen('lambda.dat','r');

    kformat = '%f';
    ksize = [1 nk];

    lformat = '%f';
    lsize = [1 nk];


    %kx = fscanf(kxfileID, kformat, ksize);
    kz = fscanf(kzfileID, kformat, ksize);
    lambda = fscanf(lfileID, lformat, lsize);


    %[KX, KZ] = meshgrid(kx, kz);

    fig = figure(1);

    %surf(KX, KZ, lambda);
    plot(kz, lambda)
    xlabel("$k_z$", 'Interpreter', 'latex')
    ylabel("Re($\lambda$)", 'Interpreter', 'latex')
    title("Re($\lambda$) versus $k_z$", 'Interpreter', 'latex')
    %zlabel("Re($lambda$)", 'Interpreter', 'latex')

%   fig = figure(2);

%   contourf(KX, KZ, lambda)
%   xlabel("kx")
%   ylabel("kz")
%   title("Contour plot of Re($\lambda$), Ri = " + string(Ri), 'Interpreter', 'latex')
%   print("ContourPlotRi="+ string(Ri) + ".pdf", "-dpdf")


