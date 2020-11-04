%close all;
clearvars;
clc

%I=double(imread('zigzag_mask.png'));
%I=mean(I,3); %To 2D matrix
I=double(imread('bike.jpg'));
%I=double(imread('noisedCircles.tif'));
%I=double(imread('phantom17.bmp'));
%I=double(imread('phantom18.bmp'));
I=mean(I,3);
I=I-min(I(:));
I=I/max(I(:));

[ni, nj]=size(I);

%%


%Lenght and area parameters
%circles.png mu=1, mu=2, mu=10
%noisedCircles.tif mu=0.1
%phantom17 mu=1, mu=2, mu=10
%phantom18 mu=0.2 mu=0.5
%hola carola mu=1
nu=0;


mu=0.2;
lambda1=4;
lambda2=1;
reIni=500;


%%Parameters
%lambda1=10^-3; %Hola carola problem
%lambda2=10^-3; %Hola carola problem

epHeaviside=0.1;
%eta=0.01;
eta=1;
tol=0.1;
%dt=(10^-2)/mu; 
dt=(10^-1)/mu;
iterMax=100000;
%reIni=0; %Try both of them
%reIni=500;
[X, Y]=meshgrid(1:nj, 1:ni);

%%Initial phi
%phi_0=(-sqrt( ( X-round(ni/2)).^2 + (Y-round(nj/2)).^2)+50);
phi_0=sin((pi/50)*X).*sin((pi/50)*Y);

%%% This initialization allows a faster convergence for phantom 18
%phi_0=(-sqrt( ( X-round(ni/2)).^2 + (Y-round(nj/4)).^2)+50);
%Normalization of the initial phi to [-1 1]
%phi_0=phi_0-min(phi_0(:));
%phi_0=2*phi_0/max(phi_0(:));
%phi_0=phi_0-1;

%phi_0=I; %For the Hola carola problem
phi_0=phi_0-min(phi_0(:));
phi_0=2*phi_0/max(phi_0(:));
phi_0=phi_0-1;

%%

%%Explicit Gradient Descent
seg=sol_ChanVeseIpol_GDExp( I, phi_0, mu, nu, eta, lambda1, lambda2, tol, epHeaviside, dt, iterMax, reIni );

