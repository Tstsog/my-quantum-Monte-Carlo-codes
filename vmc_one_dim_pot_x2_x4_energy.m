% This Matlab code performs variational Monte-Carlo (VMC) simulation and computes 
% an average energy for a partcile in ground state of one-dimensional (one dim/1D)
% harmonic-type oscillator, using a Metropolis algorithm [1]. 
% A plot of energy depending on variational parameter is obtained. 
%
% Ref. [1] E. Curotto, "Stochastic Simulations of Clusters: Quantum Methods in Flat and Curved Spaces", CPC Press (2010).
%      
% Double-well potential: V(x) = 0.5*x^2 + 0.25*x^4; 
% A trail function: psi = exp(-beta*x*x), where beta is parameter, of which optimal value needs to be found.
% The local energy analytically found: En = (beta - 2.*beta.^2.*x.^2) + (0.5.*x^2 + 0.25*x^4)
%
% An atomic units are used in calculation. 
%
% Written by Tsogbayar Tsednee (PhD)
% Contact email: tsog215@gmail.com
%
% January 12, 2025 & University of North Dakota
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [] = vmc_one_dim_pot_x2_x4_energy
%
clc; clear; format long
%
dbeta = 0.1;
beta = 0.1:dbeta:2.000;
%
fileID_save_data_1 = fopen('vmc_one_dim_pot_x2_x4_energy.txt','w');
%
for ii = 1:length(beta)
    %
    [En0_ave, sigma_std, rejection_ratio] = one_dim_ho_mc(beta(ii));
    %
    output = [ii, ii*dbeta, En0_ave, sigma_std, rejection_ratio];
    %
    fprintf(fileID_save_data_1, '%4.4f \t %8.6f \t %8.12f \t %8.12f \t %8.4f\n', output);
    %
end
%
fclose(fileID_save_data_1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_data = fopen('vmc_one_dim_pot_x2_x4_energy.txt', 'r');               % 
read_data = textscan(read_data, '%f %f %f %f %f');
beta_vals = read_data{2};
En0_vals = read_data{3};
%

En0_exact_value = 0.621; % an exact value of ground state of one-dimensional 0.5*x^2 + 0.25*x^4 system

figure(1)
hold on
yline(En0_exact_value, 'r--', 'LineWidth', 1.8);
plot(beta_vals, En0_vals, 'bo-', 'LineWidth',1.8)
hold off
axis([0.1 2 0.5 2])
xlabel('$\beta$','interpreter','latex')
ylabel('$\langle E_{0} \rangle$','interpreter','latex')
set(gca,'FontSize',18)
box on 


%%%
return
end
%
function [En0_ave, sigma_std, rejection_ratio] = one_dim_ho_mc(beta)
%
n_moves = 10^5;
%
deltax = 3.500;
%
x_init = -5.00;
%
n_rej = 0.;
sm_en0 = 0.;
sm_en0_2 = 0.;
%
x = x_init;
%
for moves = 1:2*n_moves
    %
    xt = x + deltax * (rand(1) - 0.5);
    %
    [rho] = psi_trail(beta,x);    
    [rho_t] = psi_trail(beta,xt);
    %
    q = rho_t./rho;
    %
    if (rand(1) <= q)
        x = xt ;
    else
        n_rej = n_rej + 1;
    end
    %
    if (moves < n_moves)
        %
        sm_en0 = sm_en0 + (beta - 2.*beta.^2.*x.^2) + (0.5.*x^2 + 0.25*x^4);
        sm_en0_2 = sm_en0_2 + ((beta - 2.*beta.^2.*x.^2) + (0.5.*x^2 + 0.25*x^4)).^2;
        %
    end
    %
end
%
En0_ave = sm_en0./n_moves;                            % numerical value
En0_sq_ave = sm_en0_2./n_moves;
sigma_std = sqrt((En0_sq_ave - En0_ave.^2)./n_moves); % standard deviation
rejection_ratio = 100*(n_rej./(2.*n_moves));            % rejection in percent
%%%
return
end

%%%%
function [rho] = psi_trail(beta,x)
%
rho = exp(-2.*beta.*x.*x);
%%%
return
end

