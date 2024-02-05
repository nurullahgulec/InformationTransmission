function[m_2,Delta,Phi,D_SM,Delta_H]=period2(F,init,chi_values)
mu = F.sigma/(F.sigma-1); % Markup

dx = mean(diff(chi_values)); % Differential element, assuming uniform spacing
L=10;
% Define the log-normal distribution function Gx for chi
Gx = @(chi) lognpdf(chi, F.m_phi, F.sigma_phi);
Gx_values = Gx(chi_values); % Evaluate the log-normal PDF at each grid point
N=length(chi_values);

Phi = init.Phi; % Placeholder for actual Phi function values at grid points
Delta = init.Delta; % Placeholder for actual Delta function values at grid points
Delta_H = init.Delta_H; % Placeholder for Delta_H, to be computed based on your model
%Have a guess for D_SM
D_SM=[zeros(N/2,1); ones(N/2,1)];
%D_SM=[ ones(N/2,1); zeros(N/2,1)];
D_SM_prime=NaN(N,1);
R=1;

m_1=init.m_1;
m_2=init.m_1;
iter=1;
while R > F.tolerance_epsilon & iter<50
    
    %F.latent_SM(X_r, X_o) note that
    SM= (F.sep*F.latent_SM + (1-F.sep)*(F.latent_SM .*D_SM')).*Gx_values ./ (sum((1-F.sep)*F.latent_SM .* D_SM'+F.sep*F.latent_SM,2) .*Gx_values );%% think about the distribution better
 %WE PUT D_SM' instead of D_SM since we want to get rid of those origin
 %firms that does not change sales managers
    %%
    % BE VERY CAREFUL
    %I(x,xr)=*m1(x,xo_1)*sm(xo_1, xr)+ m1(x,xo_2)*sm(xo_2,x_r)+....
    
    F.I_SM=@(x, x_r) F.kappa .* m_1(x,x_r)+ F.theta .*sum(m_2(x,:) .* SM(x_r,:)' );
    I_NOSM=m_1;
    I_SM=F.I_SM(1:N,1:N);
    I = zeros(N, N); % Initialize the final NxN matrix
    
    % Logical indexing
    I(:, D_SM == 1) = I_SM(:, D_SM == 1);
    I(:, D_SM == 0) = I_NOSM(:, D_SM == 0);
    
    
    %%
    
    % Calculate linear indices for each D_SM decision
    %
    %%Delta=[Delta_NOSM;Delta_SM];
    %Phi=[Phi_NOSM; Phi_SM];
    %linear_indices = sub2ind(size(Phi), D_SM'+1, 1:N);
    % Extract elements using linear indices
    %Phi = Phi(linear_indices);
    %Delta=Delta(linear_indices);
    
    %%
    % Step 2: Compute the implied profit function for SM and NOSM from equation (23)
    
    F.computeProfit = @(chi_ind, chi_prime_ind) F.alpha^(F.sigma - 1) * mu^(-F.sigma) * (mu - 1) ...
        .* Delta(chi_ind)' * Delta_H  .* Phi(chi_prime_ind);
    % Step 3: Compute the information flow matrix from equation (30) and (29)
    % ... Code for information flow matrix here ...
    
    % Step 4: Compute the implied matching functions from equations (33) and (32)
    % 
    F.matchingFunction2= @(chi_ind, chi_prime_ind) ...
                logncdf(F.computeProfit(chi_ind, chi_prime_ind)+F.nu*I ,F.m_eps, F.sigma_eps);
    
    m_2=F.matchingFunction2(1:N,1:N);
    %% Step 5: Compute the implied network productivity and quality functions
    % Define the Phi_prime and Delta_prime functions using summation
    Phi_prime = @(chi_ind) chi_values(chi_ind).^(F.sigma-1) + F.alpha^(F.sigma-1) * mu^(1-F.sigma) * ...
        sum(F.matchingFunction2(chi_ind, 1:N) .* Phi .* Gx_values,2)' * dx;
    Delta_prime = @(chi_ind) mu^(-F.sigma) + mu^(-F.sigma) * F.alpha^(F.sigma-1) * ...
        sum(F.matchingFunction2(1:N, chi_ind)' .* Delta .* Gx_values,2) * dx;
    
    % Compute Phi_prime and Delta_prime for the entire grid
    Phi_prime_values = Phi_prime(1:N);
    Delta_prime_values = Delta_prime(1:N);
    %% Step 6: Compute the implied household demand shifter from equations (38) to (39) and then to (24)
    Pi_matrix=F.computeProfit(1:N, 1:N);
    
    % Calculate the cumulative probabilities using the normal CDF
    
    % Calculate the cumulative probabilities using the normal CDF
    cumulative_prob_below_mean_adjusted = normcdf((log(Pi_matrix+F.nu*I)-F.m_eps - F.sigma_eps^2) / F.sigma_eps);
    
    % Compute the conditional expectation
    eps_bar = exp(F.m_eps + F.sigma_eps^2 / 2) .* (cumulative_prob_below_mean_adjusted);
    
    
    
    %%
    L_f=sum(sum(eps_bar .* Gx_values,2)*dx.*Gx_values')*dx;
    L_p=sum(chi_values.^(F.sigma-1) .*Delta .*Gx_values )*dx;
    
    Delta_H_prime=(L-L_f-F.f_sm)/L_p;
    
    %%
    % Step 7: Compute the optimal decisions of firms from equations (34) and (35)
    %Very important simplification, only matching functions will change when we
    %calculate the ideal function:
    for x=1:N
        m_sm=m_2;
        m_nosm=m_2;
        m_sm(:,x)=logncdf(F.computeProfit(1:N, x)+F.nu*I_SM(:,x) ,F.m_eps, F.sigma_eps);
        m_nosm(:,x)=logncdf(F.computeProfit(1:N, x)+F.nu*I_NOSM(:,x) ,F.m_eps, F.sigma_eps);
        pi_direct_sm= (mu-1) .* mu^(-F.sigma) *Delta_H .* Phi(x)+  ...
        sum( m_sm(:,x) .* Pi_matrix(1:N, x) .* Gx_values',1) * dx;
        pi_direct_nosm= (mu-1) .* mu^(-F.sigma) *Delta_H .* Phi(x)+  ...
        sum( m_nosm(:,x) .* Pi_matrix(1:N, x) .* Gx_values',1) * dx;
        V_nosm=pi_direct_nosm; %%VERY CRITICAL DELTA_H_PRIME
        V_sm=pi_direct_sm-F.f_sm;
        D_SM_prime(x,1)=V_sm>V_nosm;
    end
    
    %%
    % Step 8: Compute the residuals and update the guesses if needed
    % Residuals for SM and NOSM
    R_Phi = max(abs(Phi_prime_values - Phi)); % Residual for Phi 
    R_Delta = max(abs(Delta_prime_values'*Delta_H_prime - Delta_H .* Delta)); % Residual for Delta 
    R_D = max(abs(D_SM-D_SM_prime)); % Residual for Delta 

    % Overall residual R
    R = max([R_Phi, R_Delta,R_D])
    
    % Update the guesses for network productivity and scaled quality functions
    if R > F.tolerance_epsilon
        Phi = Phi_prime_values;
        Delta = Delta_prime_values';
        Delta_H=Delta_H_prime;
        D_SM=D_SM_prime;
        % Update other variables as necessary
    end
iter=iter+1;
end



end
