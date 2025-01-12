% This Matlab code performs variational Monte-Carlo (VMC) simulation and computes 
% an average energy for a partcile in ground state of one-dimensional (one dim/1D)
% harmonic-type oscillator, and probability distribution function (pdf) for it using a Metropolis algorithm [1]. 
% The PDF is compared with an exact probability function of ground state 1D harmonic-type oscillator potential. 
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
function [] = vmc_one_dim_pot_x2_x4_pdf
%
clc; clear; format short
%
beta = 0.700;    % optimal value, which gives an accurate smallest value of ground state energy, En0 = 0.621
%
one_dim_ho_mc(beta);

%%%
return
end
%
function [] = one_dim_ho_mc(beta)
%
n_moves = 10^6;
%
deltax = 3.50;
%
x_init = -5.00;
%
n_rej = 0.;
sm_en0 = 0.;
sm_en0_2 = 0.;
%
x = x_init;
%
fileID_save_data_1 = fopen('vmc_one_dim_pot_x2_x4_pdf.txt','w');
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
        sm_en0 = sm_en0 + (beta - 2.*beta.^2.*x.^2) + (0.5.*x^2 + 0.25*x^4);        % energy 
        sm_en0_2 = sm_en0_2 + ((beta - 2.*beta.^2.*x.^2) + (0.5.*x^2 + 0.25*x^4)).^2;
        %
    end
    %
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
% 0.7000    0.6242    0.0002   49.8272

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_data = fopen('vmc_one_dim_pot_x2_x4_pdf.txt', 'r');               % 
read_data = textscan(read_data, '%f %f ');
%number_of_moves = read_data{1};
x_val = read_data{2};
%

%%% exact distribution 
data_rho = [
    2.4558    0.0000
    2.3486    0.0001
    2.2400    0.0002
    2.1301    0.0005
    2.0188    0.0011
    1.9064    0.0025
    1.7928    0.0052
    1.6782    0.0103
    1.5625    0.0191
    1.4460    0.0336
    1.3285    0.0558
    1.2103    0.0878
    1.0913    0.1312
    0.9717    0.1865
    0.8515    0.2526
    0.7308    0.3271
    0.6097    0.4054
    0.4882    0.4819
    0.3664    0.5503
    0.2444    0.6044
    0.1222    0.6392
         0    0.6512
   -0.1222    0.6392
   -0.2444    0.6044
   -0.3664    0.5503
   -0.4882    0.4819
   -0.6097    0.4054
   -0.7308    0.3271
   -0.8515    0.2526
   -0.9717    0.1865
   -1.0913    0.1312
   -1.2103    0.0878
   -1.3285    0.0558
   -1.4460    0.0336
   -1.5625    0.0191
   -1.6782    0.0103
   -1.7928    0.0052
   -1.9064    0.0025
   -2.0188    0.0011
   -2.1301    0.0005
   -2.2400    0.0002
   -2.3486    0.0001
   -2.4558    0.0000
];

%
nbins = 30;
figure(1)
hold on
plot(data_rho(:,1), data_rho(:,2), 'b-', LineWidth=1.5)
histogram(x_val, nbins,'Normalization','pdf') % probability distribution function by the Monte-Carlo metropolis method 
hold off
axis([-2 2 0 0.7])
xlabel('$x$','interpreter','latex')
ylabel('$\rho(x)$','interpreter','latex')
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

