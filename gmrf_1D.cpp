
#include <TMB.hpp>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  // Data
  DATA_VECTOR( Y_i );
  DATA_SPARSE_MATRIX( M0 );
  DATA_SPARSE_MATRIX( M2 );
  
  // Parameters
  PARAMETER(lambda);
  PARAMETER( ln_kappa );
  PARAMETER( ln_tau );

  // Random effects
  PARAMETER_VECTOR( omega_i ); //random process

  // Objective funcction
  int n_i = Y_i.size();
  vector<Type> nll_comp(2);
  nll_comp.setZero();
  Type kappa = exp(ln_kappa);
  Type tau = exp(ln_tau);
  Type range = 2/kappa; // range = sqrt(8*nu)/kappa


  SparseMatrix<Type> Q = pow(kappa,2) * M0 + M2; //from Lindgren et al. 2011 when alpha=1
  nll_comp(0) += SCALE( GMRF(Q), 1/tau)( omega_i );
  SIMULATE{
    GMRF(Q).simulate(omega_i);
    omega_i = omega_i / tau;
  }

  // Probability of data conditional on random effects
  vector<Type> eta(n_i);
  for( int i=0; i<n_i; i++){
    eta(i) = lambda + omega_i(i);
    nll_comp(1) -= dpois( Y_i(i), exp(eta(i)), true ); 
    SIMULATE{
      Y_i(i) = rpois( exp(eta(i)) );
    }
  }

  // Reporting
  Type nll = nll_comp.sum();
  SIMULATE{
    REPORT( omega_i );
    REPORT( Y_i );
  }
  REPORT( eta );
  REPORT( nll_comp );
  REPORT( nll );
  REPORT( kappa );
  REPORT( range );
  REPORT( Q );
  REPORT( tau );
  REPORT( omega_i );

  return nll;
}
