function [ phi ] = sol_ChanVeseIpol_GDExp( I, phi_0, mu, nu, eta, lambda1, lambda2, tol, epHeaviside, dt, iterMax, reIni )
%Implementation of the Chan-Vese segmentation following the explicit
%gradient descent in the paper of Pascal Getreur "Chan-Vese Segmentation".
%It is the equation 19 from that paper

%I     : Gray color image to segment
%phi_0 : Initial phi
%mu    : mu lenght parameter (regularizer term)
%nu    : nu area parameter (regularizer term)
%eta   : epsilon for the total variation regularization
%lambda1, lambda2: data fidelity parameters
%tol   : tolerance for the sopping criterium
% epHeaviside: epsilon for the regularized heaviside. 
% dt     : time step
%iterMax : MAximum number of iterations
%reIni   : Iterations for reinitialization. 0 means no reinitializacion

[ni,nj]=size(I);
hi=1;
hj=1;


phi=phi_0;
dif=inf;
nIter=0;
while dif>tol && nIter<iterMax
    
    phi_old=phi;
    nIter=nIter+1;        
    
    % Calculate H(phi)
    H = (1/2)*(1+(2/pi)*atan((phi_old/epHeaviside)));
    
    %Fixed phi, Minimization with.respect.to c1 and c2 (constant estimation)
    c1 = sum(I.*H, 'all')/sum(H, 'all');
    c2 = (sum(I.*(1 - H), 'all'))/(sum(1 - H, 'all'));
    
    %Boundary conditions
    phi(1,:)   = phi(2,:);
    phi(end,:) = phi(end-1,:);

    phi(:,1)   = phi(:,2);
    phi(:,end) = phi(:,end-1);

    
    %Regularized Dirac's Delta computation
    delta_phi = G11_diracReg(phi_old, epHeaviside);   %notice delta_phi=H'(phi)	
    
    %derivatives estimation
    %i direction, forward finite differences
    phi_iFwd  = DiFwd(phi, hi);
    phi_iBwd  = DiBwd(phi, hj);
    
    %j direction, forward finitie differences
    phi_jFwd  = DjFwd(phi, hi);
    phi_jBwd  = DjBwd(phi, hj);
    
    %centered finite diferences
    phi_icent = (phi_iFwd+phi_iBwd)./2;
    phi_jcent = (phi_jFwd+phi_jBwd)./2;
    
    %A and B estimation (A y B from the Pascal Getreuer's IPOL paper "Chan
    %Vese segmentation
    A = (mu)./(sqrt(eta^2+phi_jFwd.^2+phi_icent.^2));
    B = (mu)./(sqrt(eta^2+phi_iFwd.^2+phi_jcent.^2));
    
    
    %%Equation 22, for inner points
    for i = 2:ni-1
        for j = 2:nj-1
            phi(i,j) = (phi_old(i,j) + dt * delta_phi(i,j) * ...
                ( A(i,j) * phi_old(i+1,j) + A(i-1, j) * phi(i-1,j) + B(i,j) * phi_old(i,j+1) ...
                  + B(i,j-1) * phi(i,j-1) ...
                 - nu - lambda1 * (I(i,j) - c1)^2 + lambda2 * (I(i,j) - c2)^2)) ...
                 / (1 + dt * delta_phi(i,j) * (A(i,j) + A(i-1,j) + B(i,j) + B(i,j-1)));
        end
    end   
 
    %Reinitialization of phi
    if reIni>0 && mod(nIter, reIni)==0
        indGT = phi >= 0;
        indLT = phi < 0;
        
        phi=double(bwdist(indLT) - bwdist(indGT));
        
        %Normalization [-1 1]
        nor = min(abs(min(phi(:))), max(phi(:)));
        phi=phi/nor;
    end
  
    %Diference. This stopping criterium has the problem that phi can
    %change, but not the zero level set, that it really is what we are
    %looking for.
    dif = mean(sum( (phi(:) - phi_old(:)).^2 ))
          
    %Plot the level sets surface
    subplot(1,2,1) 
        %The level set function
        surfc(1:ni,1:nj,phi,'EdgeColor','none');
        hold on
        %The zero level set over the surface
        contour(phi, [0, 0], "red"); %TODO 17: Line to complete
        hold off
        title('Phi Function');
    
    %Plot the curve evolution over the image
    subplot(1,2,2)
        imagesc(I);        
        colormap gray;
        hold on;
        contour(phi, [0, 0], "red") %TODO 18: Line to complete
        title('Image and zero level set of Phi')

        axis off;
        hold off
    drawnow;
    pause(.0001); 
end