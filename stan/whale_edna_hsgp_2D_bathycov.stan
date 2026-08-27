// =============================================================================
// whale_edna_hsgp_2D_bathycov.stan
//
// Demonstration B: a 2-D HSGP over (X, Y) PLUS bottom depth as a PARAMETRIC
// covariate (centred natural spline), rather than a 3rd GP axis:
//   log(lambda_si) = mu_s + f_s(X, Y) + B_bathy_i . beta_bathy_s
// f_s is the dimension-general HSGP (here D1 = 2); B_bathy is the spline
// design passed as data and beta_bathy the estimated per-species coefficients.
// SF to US/Canada border, UTM Zone 10N.
//
// Latent field:
//   log(lambda_si)      = mu_s + f_s(X_i, Y_i)   // HSGP GP
//   log(lambda_edna_si) = log_lambda_si
//                       + log_zsample_effect_si
//                       + log_conv_factor[s]
//                       + log_vol_filtered                  // <-- v4.1 fix
//
// Observation models - qPCR hurdle on Pacific hake (s = 1) and
// zero-inflated Beta-Binomial metabarcoding on all three species.

//

// =============================================================================

functions {

  // -------------------------------------------------------------------------
  // Zero-inflated Beta-Binomial log-PMF
  // log_p_zero = log P(absent from bottle) = -lambda  (Poisson)
  // log_p_pos  = log(1 - exp(log_p_zero))
  // count == 0 : log_sum_exp(log_p_zero, log_p_pos + BB(0|n,a,b))
  // count  > 0 : log_p_pos + BB(count|n,a,b)
  // -------------------------------------------------------------------------
  real zi_beta_binomial_lpmf(int count, int reads,
                              real log_p_zero, real log_p_pos,
                              real alpha_bb, real beta_bb) {
    real log_bb = beta_binomial_lpmf(count | reads, alpha_bb, beta_bb);
    if (count == 0) {
      return log_sum_exp(log_p_zero, log_p_pos + log_bb);
    } else {
      return log_p_pos + log_bb;
    }
  }

  // -------------------------------------------------------------------------
  // Calculated Spectral Lambda for each basis and each dimension.
  // -------------------------------------------------------------------------

	vector lambda_nD(array[] real L, array[] int m, int D) {
		vector[D] lam;
		for(i in 1:D){
			lam[i] = ((m[i]*pi())/(2*L[i]))^2; }
				
		return lam;
	}

 // -------------------------------------------------------------------------
  // Spectral density of 2-D squared-exponential kernel
  // -------------------------------------------------------------------------

	real spd_2D(real alpha, real rho1, real rho2, real w1, real w2) {
		real S;
		S = alpha^2 * sqrt(2*pi())^2 * rho1*rho2 * exp(-0.5*(rho1^2*w1^2 + rho2^2*w2^2));
				
		return S;
	}
	real spd_nD(real alpha, row_vector rho, vector w, int D) {
		real S;
		S = alpha^2 * sqrt(2*pi())^D * prod(rho) * exp(-0.5*((rho .* rho) * (w .* w)));
				
		return S;
	}

	
  // -------------------------------------------------------------------------
  // N x M HSGP design matrix (2-D product of 1-D eigenfunctions)
  // -------------------------------------------------------------------------
  vector phi_2D(real L1, real L2, int m1, int m2, vector x1, vector x2) {
		vector[rows(x1)] fi;
		vector[rows(x1)] fi1;
		vector[rows(x1)] fi2;
		fi1 = 1/sqrt(L1)*sin(m1*pi()*(x1+L1)/(2*L1));
		fi2 = 1/sqrt(L2)*sin(m2*pi()*(x2+L2)/(2*L2));
		fi = fi1 .* fi2;
		return fi;
	}
	vector phi_nD(array[] real L, array[] int m, matrix x) {
		int c = cols(x);
		int r = rows(x);
		
		matrix[r,c] fi;
		vector[r] fi1;
		for (i in 1:c){
			fi[,i] = 1/sqrt(L[i])*sin(m[i]*pi()*(x[,i]+L[i])/(2*L[i]));
		}
		fi1 = fi[,1];
		for (i in 2:c){
			fi1 = fi1 .* fi[,i];
		}
		return fi1;
	}
}

// =============================================================================
data {

  // Dimensions
  int<lower=1> N;          // total samples (station × depth combinations)
  int<lower=1> S;          // species (3)
  int<lower=1> M;          // HSGP basis dimension = prod(HSGP_M)
  int<lower=1> D1;         // Number of dimensions for the primary GP.

  array[M,D1] int INDICES;	//indices of combinations of basis functions 

  // Species-specific conversion factor (log scale): log(copies per animal
  // per litre per km^2). Whales shed more eDNA per animal than hake, so
  // this is per-species.
  vector[S+1] log_conv_factor; // +1 is for the junk category

  // log(vol_filtered) - litres of seawater filtered. 
  real log_vol_filtered;

  // Coordinates normalised to [-1, 1], one column per GP dimension (D1).
  // vS1: col 1 = X km, 2 = Y km. 3-D demo adds col 3 = Z_bathy (m).
  matrix[N, D1] coords;

  int<lower=1> N_pred;          // total prediction locations
  matrix[N_pred, D1] pred_coords;

  // Original-unit half-ranges used for normalisation, one per dimension.
  vector<lower=0>[D1] coord_scale;

  // HSGP boundary extension factors and basis counts
  array[D1] real<lower=0> L_hsgp;
  array[D1] int<lower=1> m_hsgp;

  // Fixed water-column eDNA log-offsets (pre-computed in R): N × S
  matrix[N, S] log_zsample_effect;

  // Aliquot volume (microlitres); dilution fraction = vol_aliquot / 100
  real<lower=0> vol_aliquot;

  // ------ qPCR data (hake only, s=1), long form ------
  int<lower=1>                          N_qpcr_long;
  array[N_qpcr_long] int<lower=1, upper=N> qpcr_sample_idx;
  array[N_qpcr_long] int<lower=0, upper=1> qpcr_detect;
  vector[N_qpcr_long]                      qpcr_ct;      // 0 if not detected

  // ------ Metabarcoding data (all species), long form ------
  int<lower=1>                          N_mb_long;
  array[N_mb_long] int<lower=1, upper=N> mb_sample_idx;
  array[N_mb_long, S+1] int<lower=0>      mb_reads;
  array[N_mb_long]    int<lower=0>       mb_total;

  // 1 = zero-inflated BB, 0 = plain BB
  // int<lower=0, upper=1> use_zi;

  // Fixed qPCR standard-curve coefficients (pre-estimated calibration).
  // alpha_ct, beta_ct were already data in v4. v4.1: kappa, sigma_ct
  // join them as data - the entire standard-curve fit is treated as
  // known. sigma_ct is passed at an inflated value (~0.7) to absorb the
  // Ct ~ log(integer count) discrete-count noise that the continuous
  // Ct ~ Normal model can't otherwise account for.
  real alpha_ct;
  real beta_ct;
  real<lower=0, upper=1> kappa;
  real gamma0_ct;
  real gamma1_ct;
  real<lower=0> sigma0_ct;


  // Prior hyperparameters (passed as data for easy tuning).
  // gp_sigma uses Gamma(shape, rate) - mode (alpha-1)/beta, mean
  // alpha/beta. With shape=8, rate=4: mode 1.75, mean 2.0, sd 0.71;
  // P(sigma < 0.3) ~ 4e-5 so the low-sigma trap from v3.2 is
  // unreachable.
  real prior_mu_sp_mu;
  real prior_mu_sp_sig;
  real<lower=0> prior_gp_sigma_shape;
  real<lower=0> prior_gp_sigma_rate;
  // Shared transformed gp length scale parameters (gamma distributed)
  real<lower=0> prior_gp_raw_alpha;       
  real<lower=0> prior_gp_raw_beta;
  
  real prior_beta0_phi_mu;   real<lower=0> prior_beta0_phi_sig;
  real prior_gamma0_phi_mu;  real<lower=0> prior_gamma0_phi_sig;
  real prior_gamma1_phi_mu;  real<lower=0> prior_gamma1_phi_sig;

  // Parametric bottom-depth covariate: centred natural-spline design.
  // Enters log_lambda additively; the (X,Y) field stays a 2-D HSGP.
  int<lower=1> K_bathy;
  matrix[N, K_bathy]      B_bathy;
  matrix[N_pred, K_bathy] B_bathy_pred;
  real<lower=0> prior_beta_bathy_sig;
}

// =============================================================================
transformed data {

  // HSGP design matrix — built once, reused every iteration
  //matrix[N, M] PHI = hsgp_phi(coords, L_hsgp, m_hsgp);

  matrix[N,M] PHI;
  matrix[N_pred,M] PHI_pred;
  // Calculate PHI onces from the indexes, buffer, and  coordinate locations
  for(m in 1:M){
    PHI[,m] = phi_nD(L_hsgp, INDICES[m,], coords);
    PHI_pred[,m] = phi_nD(L_hsgp, INDICES[m,], pred_coords);
  }
  // Log aliquot dilution: log(vol_aliquot / 100)
  real log_vol_frac = log(vol_aliquot / 100.0);

}

// =============================================================================
parameters {

  // Metabarcoding overdispersion
  real          beta0_phi;
  real<lower=0> gamma0_phi;
  real<lower=0> gamma1_phi;
  
  // qPCR parameters
      // alpha_ct, beta_ct, kappa, sigma_ct are all fixed in the data block.

  // GP parameters
  vector[S+1]             mu_sp;      // log-density intercepts
  vector<lower=0>[S]    gp_sigma;   // GP marginal SD per species
  matrix<lower=0>[S, D1] gp_l_raw;  // GP length-scales: lx(km), ly(km) 
  array[S] vector[M]     z_beta;    // non-centred basis coefficients
  matrix[S, K_bathy]     beta_bathy; // bottom-depth spline coefficients
}

// =============================================================================
transformed parameters {

  matrix[N, S] log_lambda;       // log true animal density
  matrix[N, S+1] log_lambda_edna;  // log eDNA (copies)
  
  {// Local variable start
  matrix[N, S] f_s;  // log eDNA (copies)
    for (s in 1:S) {
      vector[M] diagSPD ;
    
      for(m in 1:M){
        diagSPD[m] =  sqrt(spd_nD(gp_sigma[s], gp_l_raw[s,], 
                            sqrt(lambda_nD(L_hsgp, INDICES[m,], D1)), D1)); 
      } // end M loop

      f_s[,s] = PHI * (diagSPD .* z_beta[s,]);
      

    log_lambda[, s]      = mu_sp[s] + f_s[,s] + B_bathy * beta_bathy[s]';
    log_lambda_edna[, s] = log_lambda[ ,s]
                              + log_zsample_effect[, s]
                              + log_conv_factor[s]
                              + log_vol_filtered;
    
    } // end S loop
  } // end local variable.
  // add latent variable for the junk for log_lambda_edna
  // NON-SPATIAL
  //log_RE_junk = log_RE_junk_raw * sigma_junk;

      log_lambda_edna[, S+1] = rep_vector(mu_sp[S+1]
                              // + log_RE_junk[i] 
                              + log_conv_factor[S+1]
                              + log_vol_filtered                              ,N);
}

// =============================================================================
model {

  // ------------------------------------------------------------------
  // Priors
  // ------------------------------------------------------------------
  mu_sp    ~ normal(prior_mu_sp_mu, prior_mu_sp_sig);
  gp_sigma ~ gamma(prior_gp_sigma_shape, prior_gp_sigma_rate);

  for (s in 1:S) {
    for(d in 1:D1){
      gp_l_raw[s, d] ~ gamma(prior_gp_raw_alpha, prior_gp_raw_beta);
    }
  }

  for (s in 1:S) {
      z_beta[s] ~ std_normal();
  }

  // kappa, sigma_ct are fixed in the data block (no priors for now).

  to_vector(beta_bathy) ~ normal(0, prior_beta_bathy_sig);

  beta0_phi  ~ normal(prior_beta0_phi_mu,  prior_beta0_phi_sig);
  gamma0_phi ~ normal(prior_gamma0_phi_mu, prior_gamma0_phi_sig);
  gamma1_phi ~ normal(prior_gamma1_phi_mu, prior_gamma1_phi_sig);

  // ------------------------------------------------------------------
  // Likelihood 1 — qPCR hurdle (hake only, s=1)
  // ------------------------------------------------------------------
    {
    for (r in 1:N_qpcr_long) {
      int  i            = qpcr_sample_idx[r];
      // log copies in aliquot = log eDNA at depth + log dilution fraction
      real log_lam      = log_lambda_edna[i, 1] +log_vol_frac;
      real lam          = exp(log_lam);   

      if (qpcr_detect[r] == 1) {
        real mu_ct        = alpha_ct - beta_ct * log_lam ;
        real sigma_ct     = pow(pow(sigma0_ct,2) + exp(2*(gamma0_ct + gamma1_ct*log_lam )),0.5) ;
          //  qpcr_like1 += normal_lpdf(qpcr_ct[r] | mu_ct, sigma_ct);
        
        target += normal_lpdf(qpcr_ct[r] | mu_ct, sigma_ct);
            //qpcr_like1 += log1m_exp(-kappa * lam);
        target += log1m_exp(-kappa * lam);
      } else{ // if qpcr_detect[r] == 0
            //qpcr_like1 += -(kappa * lam);
        target += -(kappa * lam);
      }
    }
  }

  // ------------------------------------------------------------------
  // Likelihood 2 — Metabarcoding ZI-Beta-Binomial (all species)
  // ------------------------------------------------------------------
  {
    for (r in 1:N_mb_long) {
      int i = mb_sample_idx[r];

      // eDNA concentrations at this sample
      vector[S+1] lam_edna_i;
      vector[S+1] lam_K_edna_i; // Conditional density on presence
      vector[S+1] log_lam_edna_i;
      vector[S+1] log_lam_K_edna_i; // Conditional density on presnce
      
      for (s in 1:(S+1)){ 
        log_lam_edna_i[s]     = log_lambda_edna[i, s] +log_vol_frac  ;
        lam_edna_i[s]         = exp(log_lam_edna_i[s]);
        log_lam_K_edna_i[s]   = log_lam_edna_i[s] - log1m_exp(-lam_edna_i[s]);
        lam_K_edna_i[s]       = exp(log_lam_K_edna_i[s]);
      }
      vector[S+1] log_nu;
      vector[S+1] log_pi_i;
      vector[S+1] pi_i;
      
      log_nu = log_lam_K_edna_i - log_lam_K_edna_i[S+1];
      log_pi_i = log_nu - log_sum_exp(log_nu);
      pi_i = exp(log_pi_i);
      
      real log_lam_sum  = log_sum_exp(log_lam_edna_i);
      
      for (s in 1:(S+1)){
        // Compositional proportions
          real log_phi_s  = beta0_phi  + log_lam_sum //+beta0_phi
                              + fmax(gamma0_phi - gamma1_phi * log_lam_K_edna_i[s], 0.0);
          real phi_s      = exp(log_phi_s);   // cap phi to avoid overflow

        real alpha_bb   =  exp(log_pi_i[s] + log_phi_s) ;
        real beta_bb    =  exp(log1m_exp(log_pi_i[s]) + log_phi_s) ;

          real log_p_zero = -lam_edna_i[s] ;
          real log_p_pos  = log1m_exp(log_p_zero);

          // mb_like +=  zi_beta_binomial_lpmf(
          //                 mb_reads[r, s] | mb_total[r],
          //                 log_p_zero, log_p_pos,
          //                 alpha_bb, beta_bb
          //             );
          target += zi_beta_binomial_lpmf(
                      mb_reads[r, s] | mb_total[r],
                      log_p_zero, log_p_pos,
                      alpha_bb, beta_bb
          );
        }
        // if(r==1){
        //   print("log_lam_edna(1:4): ", log_lambda_edna[1:4,] ) ;
        //   print("log_lam_edna_i (1:4):",log_lam_edna_i) ;
        // }
      }
    }
       //print("QPCR1: ",qpcr_like1);
      // print("QPCR2: ",qpcr_like2);
      // print("QPCR3: ",qpcr_like3);
    
       //print("MB: ",mb_like);
} // End MODEL BLOCK

// =============================================================================
generated quantities {
  // Back transformed gp_ length scales.
  matrix[S,D1] gp_l;
  for(d in 1:D1){
    gp_l[,d] = (gp_l_raw[,d] * coord_scale[d]);
  }

  //vector[N_qpcr_long] log_lik_qpcr;
  //vector[N_mb_long]   log_lik_mb;

  // array[N_qpcr_long] int  pp_qpcr_detect;
  // vector[N_qpcr_long]     pp_qpcr_ct;
  
  matrix[N_mb_long,S+1] pp_pi_hat;
  array[N_mb_long, S+1] int pp_mb_reads;

  matrix[N, S] lambda_hat;
  matrix[N, S+1] lambda_edna_hat;

    for (s in 1:S) {
      lambda_hat[, s]      = exp(log_lambda[, s]);
      lambda_edna_hat[, s] = exp(log_lambda_edna[, s]);
    }
    lambda_edna_hat[, S+1] = exp(log_lambda_edna[, S+1]);
  

  // // ------------------------------------------------------------------
  // // qPCR generated quantities
  // // ------------------------------------------------------------------
  // for (r in 1:N_qpcr_long) {
  //   int  i            = qpcr_sample_idx[r];
  //   real log_lam      = log_lambda_edna[i, 1] + log_vol_frac;
  //   real lam          = exp(fmin(log_lam, 15.0));
  //   real p_det        = fmax(fmin(1.0 - exp(-kappa * lam), 1.0 - 1e-9), 1e-9);
  // 
  //   log_lik_qpcr[r] = bernoulli_lpmf(qpcr_detect[r] | p_det);
  //   if (qpcr_detect[r] == 1) {
  //     real log_lam_safe = fmax(fmin(log_lam, 15.0), -10.0);
  //     real mu_ct        = alpha_ct - beta_ct * log_lam_safe;
  //     log_lik_qpcr[r]  += normal_lpdf(qpcr_ct[r] | mu_ct, fmax(sigma_ct, 1e-6));
  //   }
  // 
  //   pp_qpcr_detect[r] = bernoulli_rng(p_det);
  //   if (pp_qpcr_detect[r] == 1) {
  //     real log_lam_safe = fmax(fmin(log_lam, 15.0), -10.0);
  //     real mu_ct        = alpha_ct - beta_ct * log_lam_safe;
  //     pp_qpcr_ct[r]     = normal_rng(mu_ct, fmax(sigma_ct, 1e-6));
  //   } else {
  //     pp_qpcr_ct[r] = 0.0;
  //   }
  // }
  // 
  // // ------------------------------------------------------------------
  // // Metabarcoding generated quantities
  // // ------------------------------------------------------------------
   {
    for (r in 1:N_mb_long) {
      int i = mb_sample_idx[r];

      // eDNA concentrations at this sample
      vector[S+1] lam_edna_i;
      vector[S+1] lam_K_edna_i; // Conditional density on presence
      vector[S+1] log_lam_edna_i;
      vector[S+1] log_lam_K_edna_i; // Conditional density on presnce
      for (s in 1:(S+1)){ 
        log_lam_edna_i[s]     = log_lambda_edna[i, s] + log_vol_frac ;
        lam_edna_i[s]         = exp(log_lam_edna_i[s]);
        log_lam_K_edna_i[s]   = log_lam_edna_i[s] - log1m_exp(-lam_edna_i[s]);
        lam_K_edna_i[s]       = exp(log_lam_K_edna_i[s]);
      }
      vector[S+1] log_nu;
      vector[S+1] log_pi_i;
      vector[S+1] pi_i;
      
      log_nu = log_lam_K_edna_i - log_lam_K_edna_i[S+1];
      log_pi_i = log_nu - log_sum_exp(log_nu);
      pi_i = exp(log_pi_i);
      pp_pi_hat[r,] = to_row_vector(pi_i);
      
      real log_lam_sum  = log_sum_exp(log_lam_edna_i);
      
      for (s in 1:(S+1)){
        
        // Compositional proportions
          real log_phi_s  = beta0_phi  + log_lam_sum //+beta0_phi
                              + fmax(gamma0_phi - gamma1_phi * log_lam_K_edna_i[s], 0.0);
          real phi_s      = exp(log_phi_s);   // cap phi to avoid overflow

        // if(r==1){  print("pi_s ",r,": ",pi_s);
        //            print("log_lam_s ",r,": ",log_lambda_edna[i, ]);
        //           print("lam_s ",r,": ",lam_s);
        //            print("lam_K_s ",r,": ",lam_K_s);
        //           print("phi_s ",r,": ",phi_s);
        //           print("lam_sum ",r,": ",lam_sum);
        // }

        // FIX 2: clamp alpha_bb and beta_bb away from zero
        real alpha_bb   =  exp(log_pi_i[s] + log_phi_s) ;
        real beta_bb    =  exp(log1m_exp(log_pi_i[s]) + log_phi_s) ;

          // print(alpha_bb);
          // print(beta_bb);

          // real log_p_zero = -lam_edna_i[s] ;
          // real log_p_pos  = log1m_exp(log_p_zero);
          // 
          // log_lik_mb[r] = zi_beta_binomial_lpmf(
          //             mb_reads[r, s] | mb_total[r],
          //             log_p_zero, log_p_pos,
          //             alpha_bb, beta_bb
          // );
        
        real p_bb         = beta_rng(alpha_bb, beta_bb);
        pp_mb_reads[r, s] = binomial_rng(mb_total[r],
                            fmax(fmin(p_bb, 1.0 - 1e-9), 1e-9));
        }
      }
    }
    
  matrix[N_pred, S] log_lambda_pred;        // log true animal density
  matrix[N_pred, S+1] log_lambda_edna_pred;  // log eDNA (copies)

  for (s in 1:S) {
    {// Local variable start
      vector[M] diagSPD ;
      vector[N_pred] f_s_tmp; // GP - smooth effect
    
      for(m in 1:M){
        diagSPD[m] =  sqrt(spd_nD(gp_sigma[s], gp_l_raw[s,], 
                            sqrt(lambda_nD(L_hsgp, INDICES[m,], D1)), D1)); 
      } // end M loop

      f_s_tmp = PHI_pred * (diagSPD .* z_beta[s,]);
    
    log_lambda_pred[, s]      = mu_sp[s] + f_s_tmp + B_bathy_pred * beta_bathy[s]';
    log_lambda_edna_pred[, s] = log_lambda_pred[ ,s]
                              //+ log_zsample_effect[, s]
                              + log_conv_factor[s]
                              + log_vol_filtered;
    } // end local variable.  
  } // end S loop
  
  // add latent variable for the junk for log_lambda_edna
  // NON-SPATIAL
  //log_RE_junk = log_RE_junk_raw * sigma_junk;

      log_lambda_edna_pred[, S+1] = rep_vector(mu_sp[S+1]
                              // + log_RE_junk[i] 
                              + log_conv_factor[S+1]
                              + log_vol_filtered,N_pred);      
}
