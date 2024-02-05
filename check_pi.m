%%

%pi_direct from formula
pi_direct= @(chi_ind) (mu-1) .* mu^(-sigma) *Delta_H .* Phi(chi_ind)+  ...
    sum(F.matchingFunction(1:N, chi_ind) .* Pi_matrix(1:N, chi_ind) .* Gx_values',1) * dx;


pi_normal= @(chi_ind) (mu-1)*Delta(chi_ind)*Phi(chi_ind)*Delta_H_prime

pi_direct(58)
pi_normal(58)