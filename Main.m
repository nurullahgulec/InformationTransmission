%%
clear
%%
% Define the parameters
F.alpha = 0.65; % Example value for alpha
F.sigma = 3; % Elasticity of substitution (example value)
%mu = sigma/(sigma-1); % Example value for mu
F.sigma_eps = 2; % Standard deviation for the log-normal distribution of epsilon
N = 5e2; % Number of points in the meshgrid
F.sigma_phi=2; %Sigma
F.m_phi=0;
F.m_eps=4;
F.tolerance_epsilon = 1e-3; % Tolerance level for convergence
F.kappa=0.75;
F.theta=0.8;
%We need to have an exogenous seperation rate:
F.sep=0.04;
F.sigma_sm=1;
F.nu=1;
F.f_sm=0.001;


% Create a meshgrid for chi
chi_min = exp(F.m_phi - 3 * F.sigma_phi); % Minimum value of x
chi_max = exp(F.m_phi + 3 * F.sigma_phi); % Maximum value of x
chi_values = linspace(chi_min, chi_max, N);
%F.latent_SM=normpdf(chi_values', chi_values, F.sigma_sm)';
F.latent_SM=1/N*ones(N,N);

[m_1,Delta,Phi, Delta_H]=period1(F,chi_values);

%%
init.Delta=Delta;
init.Phi=Phi;
init.Delta_H=Delta_H;
init.m_1=m_1;
[m_2,Delta,Phi, D_SM]=period2(F,init,chi_values);

[m_2_nosm,~,~]=period2_nosm(F,init,chi_values);

%%
m_fark=m_2-m_2_nosm;
totalconnectiondif=sum(m_fark,"all")

%%



