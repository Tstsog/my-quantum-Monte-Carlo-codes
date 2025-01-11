% This Matlab code performs variational Monte-Carlo (VMC) simulation and computes 
% an average energy for a partcile in ground state of one-dimensional (one dim/1D)
% harmonic oscillator (ho), and probability distribution function (pdf) for it using a Metropolis algorithm [1]. 
% The PDF is compared with an exact probability function of ground state 1D HO. 
%
% Ref. [1] E. Curotto, "Stochastic Simulations of Clusters: Quantum Methods in Flat and Curved Spaces", CPC Press (2010).
%      
% Double-well potential: V(x) = -0.5; 
% A trail function: psi = exp(-beta*x*x), where beta is parameter, of which optimal value needs to be found.
% The local energy analytically found: En = (beta - 2*beta^2*x^2) + 0.5*x^2
%
% An atomic units are used in calculation. 
%
% Written by Tsogbayar Tsednee (PhD)
% Contact email: tsog215@gmail.com
%
% January 11, 2025 & University of North Dakota
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [] = vmc_one_dim_ho_pdf
%
clc; clear; format short
%
beta = 0.500;    % optimal value, which gives an accurate smallest value of ground state energy, En0 = 0.5
%
one_dim_ho_mc(beta);

%%%
return
end
%
function [] = one_dim_ho_mc(beta)
%
n_moves = 10^5;
%
deltax = 5.00;
%
x_init = -5.50;
%
n_rej = 0.;
sm_en0 = 0.;
sm_en0_2 = 0.;
%
x = x_init;
%
fileID_save_data_1 = fopen('vmc_one_dim_ho_pdf.txt','w');
%
for moves = 1:2*n_moves
    %
    xt = x + deltax * (rand(1) - 0.5);
    %
    [rho] = psi_trail(beta,x);    
    [rho_t] = psi_trail(beta,xt);
    %
    q = rho_t/rho;
    %
    if (rand(1) <= q)
        x = xt ;
    else
        n_rej = n_rej + 1;
    end
    %
    if (moves < n_moves)
        %
        output = [moves, x];
        %
        fprintf(fileID_save_data_1, '%4.4f \t %8.12f\n', output);
        %
        sm_en0 = sm_en0 + (beta - 2*beta^2*x^2) + 0.5*x^2;        % energy 
        sm_en0_2 = sm_en0_2 + ((beta - 2*beta^2*x^2) + 0.5*x^2)^2;
        %
    end
%    [sm_pot/nmoves, (nrej/(2*nmoves))*100];

end
%
fclose(fileID_save_data_1);
%
En0_ave = sm_en0./n_moves;                            % numerical value
En0_sq_ave = sm_en0_2./n_moves;
sigma_std = sqrt((En0_sq_ave - En0_ave.^2)./n_moves); % standard deviation
rejection_ratio = 100*(n_rej./(2.*n_moves));            % rejection in percent
%
[beta, En0_ave, sigma_std, rejection_ratio]
%[beta, En0_ave, sigma_std, rejection_ratio]
% 0.5000    0.5000    0.0000   56.6920

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_data = fopen('vmc_one_dim_ho_pdf.txt', 'r');               % 
read_data = textscan(read_data, '%f %f ');
%number_of_moves = read_data{1};
x_val = read_data{2};
%

%%% exact distribution 
xx = x_init:0.1:abs(x_init);
psi = (1/pi).^(1/4).*exp(-0.5.*xx.^2); % ground state wave function for one-dimensional harmonic oscillator 
%
nbins = 30;
figure(1)
hold on
plot(xx, psi.^2, 'b-', LineWidth=1.5)
histogram(x_val, nbins,'Normalization','pdf') % probability distribution function by the Monte-Carlo metropolis method 
xlabel('$x$','interpreter','latex')
ylabel('$\rho(x)$','interpreter','latex')
hold off
set(gca,'FontSize',18)
box on

%%%
return
end

%%%%
function [rho] = psi_trail(beta,x)
%
rho = exp(-2.*beta*x*x);
%%%
return
end

