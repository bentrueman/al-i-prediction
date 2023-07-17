// generated with brms 2.19.0
functions {
  /* compute the logm1 link
   * Args:
   *   p: a positive scalar
   * Returns:
   *   a scalar in (-Inf, Inf)
   */
  real logm1(real y) {
    return log(y - 1.0);
  }
  /* compute the logm1 link (vectorized)
   * Args:
   *   p: a positive vector
   * Returns:
   *   a vector in (-Inf, Inf)
   */
  vector logm1_vector(vector y) {
    return log(y - 1.0);
  }
  /* compute the inverse of the logm1 link
   * Args:
   *   y: a scalar in (-Inf, Inf)
   * Returns:
   *   a positive scalar
   */
  real expp1(real y) {
    return exp(y) + 1.0;
  }
  /* compute the inverse of the logm1 link (vectorized)
   * Args:
   *   y: a vector in (-Inf, Inf)
   * Returns:
   *   a positive vector
   */
  vector expp1_vector(vector y) {
    return exp(y) + 1.0;
  }
}
data {
  int<lower=1> N; // total number of observations
  int<lower=1> N_aliugl; // number of observations
  vector[N_aliugl] Y_aliugl; // response variable
  int<lower=0> Nmi_aliugl; // number of missings
  array[Nmi_aliugl] int<lower=1> Jmi_aliugl; // positions of missings
  int<lower=1> Ksp_aliugl; // number of special effects terms
  int<lower=1> N_aluminumdisugl; // number of observations
  vector[N_aluminumdisugl] Y_aluminumdisugl; // response variable
  int<lower=0> Nmi_aluminumdisugl; // number of missings
  array[Nmi_aluminumdisugl] int<lower=1> Jmi_aluminumdisugl; // positions of missings
  int<lower=1> Ksp_aluminumdisugl; // number of special effects terms
  int<lower=1> N_titaniumdisugl; // number of observations
  vector[N_titaniumdisugl] Y_titaniumdisugl; // response variable
  int<lower=0> Nmi_titaniumdisugl; // number of missings
  array[Nmi_titaniumdisugl] int<lower=1> Jmi_titaniumdisugl; // positions of missings
  int<lower=1> Ksp_titaniumdisugl; // number of special effects terms
  int<lower=1> N_avgcolourcolourunits; // number of observations
  vector[N_avgcolourcolourunits] Y_avgcolourcolourunits; // response variable
  int<lower=0> Nmi_avgcolourcolourunits; // number of missings
  array[Nmi_avgcolourcolourunits] int<lower=1> Jmi_avgcolourcolourunits; // positions of missings
  int<lower=1> Ksp_avgcolourcolourunits; // number of special effects terms
  int<lower=1> N_avgdocmgl; // number of observations
  vector[N_avgdocmgl] Y_avgdocmgl; // response variable
  int<lower=0> Nmi_avgdocmgl; // number of missings
  array[Nmi_avgdocmgl] int<lower=1> Jmi_avgdocmgl; // positions of missings
  int<lower=1> N_irondisugl; // number of observations
  vector[N_irondisugl] Y_irondisugl; // response variable
  int<lower=0> Nmi_irondisugl; // number of missings
  array[Nmi_irondisugl] int<lower=1> Jmi_irondisugl; // positions of missings
  int<lower=1> Ksp_irondisugl; // number of special effects terms
  int<lower=1> N_ceriumdisugl; // number of observations
  vector[N_ceriumdisugl] Y_ceriumdisugl; // response variable
  int<lower=0> Nmi_ceriumdisugl; // number of missings
  array[Nmi_ceriumdisugl] int<lower=1> Jmi_ceriumdisugl; // positions of missings
  int<lower=1> Ksp_ceriumdisugl; // number of special effects terms
  int<lower=1> N_phsonde; // number of observations
  vector[N_phsonde] Y_phsonde; // response variable
  int<lower=0> Nmi_phsonde; // number of missings
  array[Nmi_phsonde] int<lower=1> Jmi_phsonde; // positions of missings
  int<lower=1> Ksp_phsonde; // number of special effects terms
  int<lower=1> N_calciumdismgl; // number of observations
  vector[N_calciumdismgl] Y_calciumdismgl; // response variable
  int<lower=0> Nmi_calciumdismgl; // number of missings
  array[Nmi_calciumdismgl] int<lower=1> Jmi_calciumdismgl; // positions of missings
  int<lower=1> N_tempdegcsonde; // number of observations
  vector[N_tempdegcsonde] Y_tempdegcsonde; // response variable
  int<lower=0> Nmi_tempdegcsonde; // number of missings
  array[Nmi_tempdegcsonde] int<lower=1> Jmi_tempdegcsonde; // positions of missings
  // data for spline s(yday, bs = "cc")
  int nb_tempdegcsonde_1; // number of bases
  array[nb_tempdegcsonde_1] int knots_tempdegcsonde_1; // number of knots
  // basis function matrices
  matrix[N_tempdegcsonde, knots_tempdegcsonde_1[1]] Zs_tempdegcsonde_1_1;
  int<lower=1> N_sulfatemgl; // number of observations
  vector[N_sulfatemgl] Y_sulfatemgl; // response variable
  int<lower=0> Nmi_sulfatemgl; // number of missings
  array[Nmi_sulfatemgl] int<lower=1> Jmi_sulfatemgl; // positions of missings
  int<lower=1> Ksp_sulfatemgl; // number of special effects terms
  int<lower=1> N_alkalinitymgcaco3l; // number of observations
  vector[N_alkalinitymgcaco3l] Y_alkalinitymgcaco3l; // response variable
  int<lower=0> Nmi_alkalinitymgcaco3l; // number of missings
  array[Nmi_alkalinitymgcaco3l] int<lower=1> Jmi_alkalinitymgcaco3l; // positions of missings
  int<lower=1> Ksp_alkalinitymgcaco3l; // number of special effects terms
  int<lower=1> N_lithiumdisugl; // number of observations
  vector[N_lithiumdisugl] Y_lithiumdisugl; // response variable
  int<lower=0> Nmi_lithiumdisugl; // number of missings
  array[Nmi_lithiumdisugl] int<lower=1> Jmi_lithiumdisugl; // positions of missings
  int<lower=1> Ksp_lithiumdisugl; // number of special effects terms
  int<lower=1> N_fluorideugl; // number of observations
  vector[N_fluorideugl] Y_fluorideugl; // response variable
  int<lower=0> Nmi_fluorideugl; // number of missings
  array[Nmi_fluorideugl] int<lower=1> Jmi_fluorideugl; // positions of missings
  // data for group-level effects of ID 1
  int<lower=1> N_1; // number of grouping levels
  int<lower=1> M_1; // number of coefficients per level
  array[N_aliugl] int<lower=1> J_1_aliugl; // grouping indicator per observation
  // group-level predictor values
  vector[N_aliugl] Z_1_aliugl_1;
  int prior_only; // should the likelihood be ignored?
  int<lower=0> Ncens_aliugl; // number of left-censored
  array[Ncens_aliugl] int<lower=1> Jcens_aliugl; // positions of left-censored
  real U_aliugl; // left-censoring limit
  int<lower=0> Ncens_alkalinitymgcaco3l; // number of left-censored
  array[Ncens_alkalinitymgcaco3l] int<lower=1> Jcens_alkalinitymgcaco3l; // positions of left-censored
  real U_alkalinitymgcaco3l; // left-censoring limit
  int<lower=0> Ncens_calciumdismgl; // number of left-censored
  array[Ncens_calciumdismgl] int<lower=1> Jcens_calciumdismgl; // positions of left-censored
  real U_calciumdismgl; // left-censoring limit
  int<lower=0> Ncens_ceriumdisugl; // number of left-censored
  array[Ncens_ceriumdisugl] int<lower=1> Jcens_ceriumdisugl; // positions of left-censored
  real U_ceriumdisugl; // left-censoring limit
  int<lower=0> Ncens_fluorideugl; // number of left-censored
  array[Ncens_fluorideugl] int<lower=1> Jcens_fluorideugl; // positions of left-censored
  real U_fluorideugl; // left-censoring limit
  int<lower=0> Ncens_irondisugl; // number of left-censored
  array[Ncens_irondisugl] int<lower=1> Jcens_irondisugl; // positions of left-censored
  real U_irondisugl; // left-censoring limit
  int<lower=0> Ncens_lithiumdisugl; // number of left-censored
  array[Ncens_lithiumdisugl] int<lower=1> Jcens_lithiumdisugl; // positions of left-censored
  real U_lithiumdisugl; // left-censoring limit
  int<lower=0> Ncens_sulfatemgl; // number of left-censored
  array[Ncens_sulfatemgl] int<lower=1> Jcens_sulfatemgl; // positions of left-censored
  real U_sulfatemgl; // left-censoring limit
  int<lower=0> Ncens_titaniumdisugl; // number of left-censored
  array[Ncens_titaniumdisugl] int<lower=1> Jcens_titaniumdisugl; // positions of left-censored
  real U_titaniumdisugl; // left-censoring limit
}
transformed data {
  
}
parameters {
  vector[Nmi_aliugl] Ymi_aliugl; // estimated missings
  vector[Ksp_aliugl] bsp_aliugl; // special effects coefficients
  real<lower=0> sigma_aliugl; // dispersion parameter
  real<lower=1> nu_aliugl; // degrees of freedom or shape
  vector[Nmi_aluminumdisugl] Ymi_aluminumdisugl; // estimated missings
  vector[Ksp_aluminumdisugl] bsp_aluminumdisugl; // special effects coefficients
  real<lower=0> sigma_aluminumdisugl; // dispersion parameter
  real<lower=1> nu_aluminumdisugl; // degrees of freedom or shape
  vector[Nmi_titaniumdisugl] Ymi_titaniumdisugl; // estimated missings
  vector[Ksp_titaniumdisugl] bsp_titaniumdisugl; // special effects coefficients
  real<lower=0> sigma_titaniumdisugl; // dispersion parameter
  real<lower=1> nu_titaniumdisugl; // degrees of freedom or shape
  vector[Nmi_avgcolourcolourunits] Ymi_avgcolourcolourunits; // estimated missings
  vector[Ksp_avgcolourcolourunits] bsp_avgcolourcolourunits; // special effects coefficients
  real<lower=0> sigma_avgcolourcolourunits; // dispersion parameter
  real<lower=1> nu_avgcolourcolourunits; // degrees of freedom or shape
  vector[Nmi_avgdocmgl] Ymi_avgdocmgl; // estimated missings
  real Intercept_avgdocmgl; // temporary intercept for centered predictors
  real<lower=0> sigma_avgdocmgl; // dispersion parameter
  real<lower=1> nu_avgdocmgl; // degrees of freedom or shape
  vector[Nmi_irondisugl] Ymi_irondisugl; // estimated missings
  vector[Ksp_irondisugl] bsp_irondisugl; // special effects coefficients
  real<lower=0> sigma_irondisugl; // dispersion parameter
  real<lower=1> nu_irondisugl; // degrees of freedom or shape
  vector[Nmi_ceriumdisugl] Ymi_ceriumdisugl; // estimated missings
  vector[Ksp_ceriumdisugl] bsp_ceriumdisugl; // special effects coefficients
  real<lower=0> sigma_ceriumdisugl; // dispersion parameter
  real<lower=1> nu_ceriumdisugl; // degrees of freedom or shape
  vector[Nmi_phsonde] Ymi_phsonde; // estimated missings
  vector[Ksp_phsonde] bsp_phsonde; // special effects coefficients
  real<lower=0> sigma_phsonde; // dispersion parameter
  real<lower=1> nu_phsonde; // degrees of freedom or shape
  vector[Nmi_calciumdismgl] Ymi_calciumdismgl; // estimated missings
  real Intercept_calciumdismgl; // temporary intercept for centered predictors
  real<lower=0> sigma_calciumdismgl; // dispersion parameter
  real<lower=1> nu_calciumdismgl; // degrees of freedom or shape
  vector[Nmi_tempdegcsonde] Ymi_tempdegcsonde; // estimated missings
  // parameters for spline s(yday, bs = "cc")
  // standarized spline coefficients
  vector[knots_tempdegcsonde_1[1]] zs_tempdegcsonde_1_1;
  real<lower=0> sds_tempdegcsonde_1_1; // standard deviations of spline coefficients
  real<lower=0> sigma_tempdegcsonde; // dispersion parameter
  real<lower=1> nu_tempdegcsonde; // degrees of freedom or shape
  vector[Nmi_sulfatemgl] Ymi_sulfatemgl; // estimated missings
  vector[Ksp_sulfatemgl] bsp_sulfatemgl; // special effects coefficients
  real<lower=0> sigma_sulfatemgl; // dispersion parameter
  real<lower=1> nu_sulfatemgl; // degrees of freedom or shape
  vector[Nmi_alkalinitymgcaco3l] Ymi_alkalinitymgcaco3l; // estimated missings
  vector[Ksp_alkalinitymgcaco3l] bsp_alkalinitymgcaco3l; // special effects coefficients
  real<lower=0> sigma_alkalinitymgcaco3l; // dispersion parameter
  real<lower=1> nu_alkalinitymgcaco3l; // degrees of freedom or shape
  vector[Nmi_lithiumdisugl] Ymi_lithiumdisugl; // estimated missings
  vector[Ksp_lithiumdisugl] bsp_lithiumdisugl; // special effects coefficients
  real<lower=0> sigma_lithiumdisugl; // dispersion parameter
  real<lower=1> nu_lithiumdisugl; // degrees of freedom or shape
  vector[Nmi_fluorideugl] Ymi_fluorideugl; // estimated missings
  real Intercept_fluorideugl; // temporary intercept for centered predictors
  real<lower=0> sigma_fluorideugl; // dispersion parameter
  real<lower=1> nu_fluorideugl; // degrees of freedom or shape
  vector<lower=0>[M_1] sd_1; // group-level standard deviations
  array[M_1] vector[N_1] z_1; // standardized group-level effects
  real<lower=0> tau;
  vector<upper=U_aliugl>[Ncens_aliugl] Ycens_aliugl; // estimated left-censored
  vector<upper=U_alkalinitymgcaco3l>[Ncens_alkalinitymgcaco3l] Ycens_alkalinitymgcaco3l; // estimated left-censored
  vector<upper=U_calciumdismgl>[Ncens_calciumdismgl] Ycens_calciumdismgl; // estimated left-censored
  vector<upper=U_ceriumdisugl>[Ncens_ceriumdisugl] Ycens_ceriumdisugl; // estimated left-censored
  vector<upper=U_fluorideugl>[Ncens_fluorideugl] Ycens_fluorideugl; // estimated left-censored
  vector<upper=U_irondisugl>[Ncens_irondisugl] Ycens_irondisugl; // estimated left-censored
  vector<upper=U_lithiumdisugl>[Ncens_lithiumdisugl] Ycens_lithiumdisugl; // estimated left-censored
  vector<upper=U_sulfatemgl>[Ncens_sulfatemgl] Ycens_sulfatemgl; // estimated left-censored
  vector<upper=U_titaniumdisugl>[Ncens_titaniumdisugl] Ycens_titaniumdisugl; // estimated left-censored
}
transformed parameters {
  // actual spline coefficients
  vector[knots_tempdegcsonde_1[1]] s_tempdegcsonde_1_1;
  vector[N_1] r_1_aliugl_1; // actual group-level effects
  real lprior = 0; // prior contributions to the log posterior
  // compute actual spline coefficients
  s_tempdegcsonde_1_1 = sds_tempdegcsonde_1_1 * zs_tempdegcsonde_1_1;
  r_1_aliugl_1 = sd_1[1] * z_1[1];
  lprior += normal_lpdf(bsp_aliugl | 0, tau);
  lprior += student_t_lpdf(sigma_aliugl | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(nu_aliugl | 2, 0.1) - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += normal_lpdf(bsp_aluminumdisugl | 0, 1);
  lprior += student_t_lpdf(sigma_aluminumdisugl | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(nu_aluminumdisugl | 2, 0.1)
            - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += normal_lpdf(bsp_titaniumdisugl | 0, 1);
  lprior += student_t_lpdf(sigma_titaniumdisugl | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(nu_titaniumdisugl | 2, 0.1)
            - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += normal_lpdf(bsp_avgcolourcolourunits | 0, 1);
  lprior += student_t_lpdf(sigma_avgcolourcolourunits | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(nu_avgcolourcolourunits | 2, 0.1)
            - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += normal_lpdf(Intercept_avgdocmgl | 0, 0.1);
  lprior += student_t_lpdf(sigma_avgdocmgl | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(nu_avgdocmgl | 2, 0.1) - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += normal_lpdf(bsp_irondisugl | 0, 1);
  lprior += student_t_lpdf(sigma_irondisugl | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(nu_irondisugl | 2, 0.1) - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += normal_lpdf(bsp_ceriumdisugl | 0, 1);
  lprior += student_t_lpdf(sigma_ceriumdisugl | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(nu_ceriumdisugl | 2, 0.1)
            - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += normal_lpdf(bsp_phsonde | 0, 1);
  lprior += student_t_lpdf(sigma_phsonde | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(nu_phsonde | 2, 0.1) - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += normal_lpdf(Intercept_calciumdismgl | 0, 0.1);
  lprior += student_t_lpdf(sigma_calciumdismgl | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(nu_calciumdismgl | 2, 0.1)
            - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += student_t_lpdf(sds_tempdegcsonde_1_1 | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sigma_tempdegcsonde | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(nu_tempdegcsonde | 2, 0.1)
            - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += normal_lpdf(bsp_sulfatemgl | 0, 1);
  lprior += student_t_lpdf(sigma_sulfatemgl | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(nu_sulfatemgl | 2, 0.1) - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += normal_lpdf(bsp_alkalinitymgcaco3l | 0, 1);
  lprior += student_t_lpdf(sigma_alkalinitymgcaco3l | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(nu_alkalinitymgcaco3l | 2, 0.1)
            - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += normal_lpdf(bsp_lithiumdisugl | 0, 1);
  lprior += student_t_lpdf(sigma_lithiumdisugl | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(nu_lithiumdisugl | 2, 0.1)
            - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += normal_lpdf(Intercept_fluorideugl | 0, 0.1);
  lprior += student_t_lpdf(sigma_fluorideugl | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += gamma_lpdf(nu_fluorideugl | 2, 0.1) - 1 * gamma_lccdf(1 | 2, 0.1);
  lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // vector combining observed and missing responses
    vector[N_aliugl] Yl_aliugl = Y_aliugl;
    // vector combining observed and missing responses
    vector[N_aluminumdisugl] Yl_aluminumdisugl = Y_aluminumdisugl;
    // vector combining observed and missing responses
    vector[N_titaniumdisugl] Yl_titaniumdisugl = Y_titaniumdisugl;
    // vector combining observed and missing responses
    vector[N_avgcolourcolourunits] Yl_avgcolourcolourunits = Y_avgcolourcolourunits;
    // vector combining observed and missing responses
    vector[N_avgdocmgl] Yl_avgdocmgl = Y_avgdocmgl;
    // vector combining observed and missing responses
    vector[N_irondisugl] Yl_irondisugl = Y_irondisugl;
    // vector combining observed and missing responses
    vector[N_ceriumdisugl] Yl_ceriumdisugl = Y_ceriumdisugl;
    // vector combining observed and missing responses
    vector[N_phsonde] Yl_phsonde = Y_phsonde;
    // vector combining observed and missing responses
    vector[N_calciumdismgl] Yl_calciumdismgl = Y_calciumdismgl;
    // vector combining observed and missing responses
    vector[N_tempdegcsonde] Yl_tempdegcsonde = Y_tempdegcsonde;
    // vector combining observed and missing responses
    vector[N_sulfatemgl] Yl_sulfatemgl = Y_sulfatemgl;
    // vector combining observed and missing responses
    vector[N_alkalinitymgcaco3l] Yl_alkalinitymgcaco3l = Y_alkalinitymgcaco3l;
    // vector combining observed and missing responses
    vector[N_lithiumdisugl] Yl_lithiumdisugl = Y_lithiumdisugl;
    // vector combining observed and missing responses
    vector[N_fluorideugl] Yl_fluorideugl = Y_fluorideugl;
    // initialize linear predictor term
    vector[N_aliugl] mu_aliugl = rep_vector(0.0, N_aliugl);
    // initialize linear predictor term
    vector[N_aluminumdisugl] mu_aluminumdisugl = rep_vector(0.0,
                                                            N_aluminumdisugl);
    // initialize linear predictor term
    vector[N_titaniumdisugl] mu_titaniumdisugl = rep_vector(0.0,
                                                            N_titaniumdisugl);
    // initialize linear predictor term
    vector[N_avgcolourcolourunits] mu_avgcolourcolourunits = rep_vector(0.0,
                                                                    N_avgcolourcolourunits);
    // initialize linear predictor term
    vector[N_avgdocmgl] mu_avgdocmgl = rep_vector(0.0, N_avgdocmgl);
    // initialize linear predictor term
    vector[N_irondisugl] mu_irondisugl = rep_vector(0.0, N_irondisugl);
    // initialize linear predictor term
    vector[N_ceriumdisugl] mu_ceriumdisugl = rep_vector(0.0, N_ceriumdisugl);
    // initialize linear predictor term
    vector[N_phsonde] mu_phsonde = rep_vector(0.0, N_phsonde);
    // initialize linear predictor term
    vector[N_calciumdismgl] mu_calciumdismgl = rep_vector(0.0,
                                                          N_calciumdismgl);
    // initialize linear predictor term
    vector[N_tempdegcsonde] mu_tempdegcsonde = rep_vector(0.0,
                                                          N_tempdegcsonde);
    // initialize linear predictor term
    vector[N_sulfatemgl] mu_sulfatemgl = rep_vector(0.0, N_sulfatemgl);
    // initialize linear predictor term
    vector[N_alkalinitymgcaco3l] mu_alkalinitymgcaco3l = rep_vector(0.0,
                                                                    N_alkalinitymgcaco3l);
    // initialize linear predictor term
    vector[N_lithiumdisugl] mu_lithiumdisugl = rep_vector(0.0,
                                                          N_lithiumdisugl);
    // initialize linear predictor term
    vector[N_fluorideugl] mu_fluorideugl = rep_vector(0.0, N_fluorideugl);
    Yl_aliugl[Jmi_aliugl] = Ymi_aliugl;
    Yl_aluminumdisugl[Jmi_aluminumdisugl] = Ymi_aluminumdisugl;
    Yl_titaniumdisugl[Jmi_titaniumdisugl] = Ymi_titaniumdisugl;
    Yl_avgcolourcolourunits[Jmi_avgcolourcolourunits] = Ymi_avgcolourcolourunits;
    Yl_avgdocmgl[Jmi_avgdocmgl] = Ymi_avgdocmgl;
    Yl_irondisugl[Jmi_irondisugl] = Ymi_irondisugl;
    Yl_ceriumdisugl[Jmi_ceriumdisugl] = Ymi_ceriumdisugl;
    Yl_phsonde[Jmi_phsonde] = Ymi_phsonde;
    Yl_calciumdismgl[Jmi_calciumdismgl] = Ymi_calciumdismgl;
    Yl_tempdegcsonde[Jmi_tempdegcsonde] = Ymi_tempdegcsonde;
    Yl_sulfatemgl[Jmi_sulfatemgl] = Ymi_sulfatemgl;
    Yl_alkalinitymgcaco3l[Jmi_alkalinitymgcaco3l] = Ymi_alkalinitymgcaco3l;
    Yl_lithiumdisugl[Jmi_lithiumdisugl] = Ymi_lithiumdisugl;
    Yl_fluorideugl[Jmi_fluorideugl] = Ymi_fluorideugl;
    Yl_aliugl[Jcens_aliugl] = Ycens_aliugl; // add imputed left-censored values
    Yl_alkalinitymgcaco3l[Jcens_alkalinitymgcaco3l] = Ycens_alkalinitymgcaco3l; // add imputed left-censored values
    Yl_calciumdismgl[Jcens_calciumdismgl] = Ycens_calciumdismgl; // add imputed left-censored values
    Yl_ceriumdisugl[Jcens_ceriumdisugl] = Ycens_ceriumdisugl; // add imputed left-censored values
    Yl_fluorideugl[Jcens_fluorideugl] = Ycens_fluorideugl; // add imputed left-censored values
    Yl_irondisugl[Jcens_irondisugl] = Ycens_irondisugl; // add imputed left-censored values
    Yl_lithiumdisugl[Jcens_lithiumdisugl] = Ycens_lithiumdisugl; // add imputed left-censored values
    Yl_sulfatemgl[Jcens_sulfatemgl] = Ycens_sulfatemgl; // add imputed left-censored values
    Yl_titaniumdisugl[Jcens_titaniumdisugl] = Ycens_titaniumdisugl; // add imputed left-censored values
    mu_avgdocmgl += Intercept_avgdocmgl;
    mu_calciumdismgl += Intercept_calciumdismgl;
    mu_tempdegcsonde += Zs_tempdegcsonde_1_1 * s_tempdegcsonde_1_1;
    mu_fluorideugl += Intercept_fluorideugl;
    for (n in 1 : N_aliugl) {
      // add more terms to the linear predictor
      mu_aliugl[n] += bsp_aliugl[1] * Yl_aluminumdisugl[n]
                      + bsp_aliugl[2] * Yl_titaniumdisugl[n]
                      + bsp_aliugl[3] * Yl_avgdocmgl[n]
                      + bsp_aliugl[4] * Yl_avgcolourcolourunits[n]
                      + bsp_aliugl[5] * Yl_irondisugl[n]
                      + bsp_aliugl[6] * Yl_ceriumdisugl[n]
                      + bsp_aliugl[7] * Yl_phsonde[n]
                      + bsp_aliugl[8] * Yl_tempdegcsonde[n]
                      + bsp_aliugl[9] * Yl_sulfatemgl[n]
                      + bsp_aliugl[10] * Yl_alkalinitymgcaco3l[n]
                      + bsp_aliugl[11] * Yl_lithiumdisugl[n]
                      + bsp_aliugl[12] * Yl_calciumdismgl[n]
                      + bsp_aliugl[13] * Yl_fluorideugl[n]
                      + r_1_aliugl_1[J_1_aliugl[n]] * Z_1_aliugl_1[n];
    }
    for (n in 1 : N_aluminumdisugl) {
      // add more terms to the linear predictor
      mu_aluminumdisugl[n] += bsp_aluminumdisugl[1] * Yl_avgdocmgl[n];
    }
    for (n in 1 : N_titaniumdisugl) {
      // add more terms to the linear predictor
      mu_titaniumdisugl[n] += bsp_titaniumdisugl[1] * Yl_aluminumdisugl[n];
    }
    for (n in 1 : N_avgcolourcolourunits) {
      // add more terms to the linear predictor
      mu_avgcolourcolourunits[n] += bsp_avgcolourcolourunits[1]
                                    * Yl_avgdocmgl[n];
    }
    for (n in 1 : N_irondisugl) {
      // add more terms to the linear predictor
      mu_irondisugl[n] += bsp_irondisugl[1] * Yl_avgdocmgl[n];
    }
    for (n in 1 : N_ceriumdisugl) {
      // add more terms to the linear predictor
      mu_ceriumdisugl[n] += bsp_ceriumdisugl[1] * Yl_irondisugl[n];
    }
    for (n in 1 : N_phsonde) {
      // add more terms to the linear predictor
      mu_phsonde[n] += bsp_phsonde[1] * Yl_aluminumdisugl[n];
    }
    for (n in 1 : N_sulfatemgl) {
      // add more terms to the linear predictor
      mu_sulfatemgl[n] += bsp_sulfatemgl[1] * Yl_calciumdismgl[n];
    }
    for (n in 1 : N_alkalinitymgcaco3l) {
      // add more terms to the linear predictor
      mu_alkalinitymgcaco3l[n] += bsp_alkalinitymgcaco3l[1]
                                  * Yl_calciumdismgl[n];
    }
    for (n in 1 : N_lithiumdisugl) {
      // add more terms to the linear predictor
      mu_lithiumdisugl[n] += bsp_lithiumdisugl[1] * Yl_fluorideugl[n];
    }
    target += student_t_lpdf(Yl_aliugl | nu_aliugl, mu_aliugl, sigma_aliugl);
    target += student_t_lpdf(Yl_aluminumdisugl | nu_aluminumdisugl, mu_aluminumdisugl, sigma_aluminumdisugl);
    target += student_t_lpdf(Yl_titaniumdisugl | nu_titaniumdisugl, mu_titaniumdisugl, sigma_titaniumdisugl);
    target += student_t_lpdf(Yl_avgcolourcolourunits | nu_avgcolourcolourunits, mu_avgcolourcolourunits, sigma_avgcolourcolourunits);
    target += student_t_lpdf(Yl_avgdocmgl | nu_avgdocmgl, mu_avgdocmgl, sigma_avgdocmgl);
    target += student_t_lpdf(Yl_irondisugl | nu_irondisugl, mu_irondisugl, sigma_irondisugl);
    target += student_t_lpdf(Yl_ceriumdisugl | nu_ceriumdisugl, mu_ceriumdisugl, sigma_ceriumdisugl);
    target += student_t_lpdf(Yl_phsonde | nu_phsonde, mu_phsonde, sigma_phsonde);
    target += student_t_lpdf(Yl_calciumdismgl | nu_calciumdismgl, mu_calciumdismgl, sigma_calciumdismgl);
    target += student_t_lpdf(Yl_tempdegcsonde | nu_tempdegcsonde, mu_tempdegcsonde, sigma_tempdegcsonde);
    target += student_t_lpdf(Yl_sulfatemgl | nu_sulfatemgl, mu_sulfatemgl, sigma_sulfatemgl);
    target += student_t_lpdf(Yl_alkalinitymgcaco3l | nu_alkalinitymgcaco3l, mu_alkalinitymgcaco3l, sigma_alkalinitymgcaco3l);
    target += student_t_lpdf(Yl_lithiumdisugl | nu_lithiumdisugl, mu_lithiumdisugl, sigma_lithiumdisugl);
    target += student_t_lpdf(Yl_fluorideugl | nu_fluorideugl, mu_fluorideugl, sigma_fluorideugl);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zs_tempdegcsonde_1_1);
  target += std_normal_lpdf(z_1[1]);
  target += cauchy_lpdf(tau | 0, 1) - cauchy_lccdf(0 | 0, 1);
}
generated quantities {
  // actual population-level intercept
  real b_avgdocmgl_Intercept = Intercept_avgdocmgl;
  // actual population-level intercept
  real b_calciumdismgl_Intercept = Intercept_calciumdismgl;
  // actual population-level intercept
  real b_fluorideugl_Intercept = Intercept_fluorideugl;
}


