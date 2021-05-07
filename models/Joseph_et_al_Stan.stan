data {
int<lower = 1> m_imbcr;
int<lower = 1> time_bins;
int<lower = 1> dist_bins;
real<lower = 0> max_dist;
real<lower = 0, upper = max_dist> dist_bin_width;
vector[dist_bins] midpoints;
int<lower = 1> ndim;
int<lower = 1> n_observer;
int<lower = 1> n_transect;
int<lower = 1> n_cells;
int<lower = 1> n_edges;
int<lower=1, upper=n_cells> node1[n_edges];
int<lower=1, upper=n_cells> node2[n_edges];
// train data
int<lower = 0> n_train_detections;
int<lower = 1> n_train_counts;
int<lower = 0> n_train[n_train_counts];
matrix[n_train_counts, m_imbcr] X_train;
int<lower = 1, upper = n_cells> train_cell[n_train_counts];
int<lower = 1, upper = time_bins> t_train[n_train_detections];
int<lower = 1, upper = dist_bins> d_train[n_train_detections];
int<lower = 1, upper = n_train_counts> idx_train_det[n_train_detections];
int<lower = 1, upper = n_observer> observer_train[n_train_counts];
int<lower = 1, upper = n_transect> transect_train[n_train_counts];
// bbs data
int<lower = 1> nsite_bbs;
int<lower = 1, upper = n_cells> cell_bbs[nsite_bbs];
int<lower = 1> m_bbs;
matrix[nsite_bbs, m_bbs] X_bbs;
int<lower = 0> bbs_count[nsite_bbs];
// ebird data
int<lower = 1> nsite_ebd;
int<lower = 1, upper = n_cells> cell_ebd[nsite_ebd];
int<lower = 1> m_ebd;
matrix[nsite_ebd, m_ebd] X_ebd;
int<lower = 0> ebd_count[nsite_ebd];
}
transformed data {
vector[dist_bins] log_f = log(2 * midpoints * dist_bin_width / max_dist^2);
vector[dist_bins] midpts_sq_divby_2;
real<lower = 0> eps = 1e-8;
for (i in 1:dist_bins)
midpts_sq_divby_2[i] = midpoints[i]^2 / 2;
}


parameters {
real<lower=0, upper = 1> a;
vector<lower = 0>[ndim] alpha;
vector[m_imbcr] beta_imbcr;
matrix[n_cells, ndim] u;
cholesky_factor_corr[ndim] L;
real<lower = 0> alpha_d;
real<lower = 0> beta_d;
vector<lower=0>[n_observer] sig_sq;
vector[m_bbs] beta_bbs;
real<lower = 0> phi_bbs;
vector[m_ebd] beta_ebd;
real<lower = 0> phi_ebd;
}


transformed parameters {
matrix[ndim, n_cells] Theta;
vector[n_train_counts] log_p_detect;
vector[n_train_counts] log_lambda;
vector[time_bins] pi_a_c;
matrix[n_train_counts, dist_bins] pi_d_c;
real log_a = log(a);
real log1m_a = log1m(a);
vector[time_bins] log_pi_a;
// coregionalization component (add sd)
Theta = diag_pre_multiply(alpha, L) * u';
log_lambda = X_train * beta_imbcr
+ Theta[1, train_cell]';
for (t in 1:time_bins)
log_pi_a[t] = log_a + (t-1) * log1m_a;
pi_a_c = softmax(log_pi_a);
{
// local scope for temp variables
vector[dist_bins] log_g;
vector[dist_bins] log_pi_d;
for (i in 1:n_train_counts) {
log_g = - midpts_sq_divby_2 / sig_sq[observer_train[i]];
log_pi_d = log_g + log_f;
pi_d_c[i, ] = softmax(log_pi_d)';
log_p_detect[i] = log_sum_exp(log_pi_d) + log_sum_exp(log_pi_a);
}
}
}
model {
beta_imbcr ~ std_normal();
alpha ~ std_normal();
L ~ lkj_corr_cholesky(5);
alpha_d ~ std_normal();
beta_d ~ std_normal();
sig_sq ~ gamma(alpha_d, beta_d);
// icar
for (i in 1:ndim) {
target += -0.5 * dot_self(u[node1, i] - u[node2, i]);
// soft sum-to-zero constraint on phi)
sum(u[, i]) ~ normal(0, 0.01 * n_cells);
}
n_train ~ poisson_log(log_lambda + log_p_detect);
for (i in 1:n_train_detections) {
t_train[i] ~ categorical(pi_a_c);
d_train[i] ~ categorical(to_vector(pi_d_c[idx_train_det[i], ]));
}
beta_bbs ~ std_normal();
phi_bbs ~ std_normal();
bbs_count ~ neg_binomial_2_log(X_bbs * beta_bbs + Theta[2, cell_bbs]', phi_bbs);
beta_ebd ~ std_normal();
phi_ebd ~ std_normal();
ebd_count ~ neg_binomial_2_log(X_ebd * beta_ebd + Theta[3, cell_ebd]', phi_ebd);
}
