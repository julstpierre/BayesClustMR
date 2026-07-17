data {
  int<lower=1> J;  
  int<lower=1> K;  
  vector[K] rho;
  vector[K] phi;
  vector[J] theta_estimates;  
  vector<lower=0>[J] sigma_estimates;  
  int<lower=0, upper=1> num_null;
  int<lower=0, upper=1> num_junk;
  real<lower=1> deg_freedom;
  matrix[J, K + num_null + num_junk] etas;
}


parameters {
  ordered[K] cluster_center;
  //vector[K] cluster_center;
  vector<lower=0>[K + num_null + num_junk] alpha;
  array[J] simplex[K + num_null + num_junk] pi;
}


model {
  int num_of_clusters = K + num_null + num_junk;
  
  cluster_center ~ normal(rho, phi);
  //cluster_center ~ normal(0, 5);
  
  //alpha ~ gamma(ga, kappa);
  //alpha ~ normal(0, 2);
  //alpha ~ exponential(0.02);
  //alpha ~ gamma(2, 1);
  
  alpha ~ lognormal(0, 0.5);

  
  if (num_of_clusters > 1) {
    for (j in 1:J) {
      vector[num_of_clusters] diri_param = to_vector(etas[j]);
      pi[j] ~ dirichlet(alpha .* diri_param);
      //pi[j] ~ dirichlet(alpha * diri_param);
    }
  }
  
  for (j in 1:J) {
    vector[num_of_clusters] likelihood = log(pi[j]);
    real abs_sigma_hat = abs(sigma_estimates[j]);
    for (num in 1:num_of_clusters) {
      if (num_null == 1){  // There is a null cluster
        if (num == 1){
          likelihood[num] += normal_lpdf(theta_estimates[j] | 0, abs_sigma_hat);
        }
        else{
          if (num_junk == 1){  // There is a junk cluster
            if (num < num_of_clusters){
                likelihood[num] += normal_lpdf(theta_estimates[j] | cluster_center[num - 1], abs_sigma_hat);
            }
            else{
              likelihood[num] += student_t_lpdf(theta_estimates[j] | deg_freedom, 0, 1);
            }
          }
          else { // There isn’t a junk cluster
            likelihood[num] += normal_lpdf(theta_estimates[j] | cluster_center[num - 1], abs_sigma_hat);
          }
        }
      }
      else{ // There isn't a null cluster
        if (num_junk == 1){ // There is a junk cluster
          if (num < num_of_clusters){
              likelihood[num] += normal_lpdf(theta_estimates[j] | cluster_center[num], abs_sigma_hat);
          }
          else{ // There is a junk cluster
            likelihood[num] += student_t_lpdf(theta_estimates[j] | deg_freedom, 0, 1);
          }
        }
          
        else {// There isn’t a junk cluster
          likelihood[num] += normal_lpdf(theta_estimates[j] | cluster_center[num], abs_sigma_hat);
        }
      }
    }
  target += log_sum_exp(likelihood);
  }
}

generated quantities {
  int num_of_clusters = K + num_null + num_junk;
  vector[J] log_lik;
  array[J] simplex[K + num_null + num_junk] post_prob;

  for (j in 1:J) {
    vector[num_of_clusters] likelihood = rep_vector(0, num_of_clusters);
    real abs_sigma_hat = abs(sigma_estimates[j]);

    for (num in 1:num_of_clusters) {
      likelihood[num] = log(pi[j][num]);

      if (num_null == 1){
        if (num == 1){
          likelihood[num] += normal_lpdf(theta_estimates[j] | 0, abs_sigma_hat);
        } else {
          if (num_junk == 1){
            if (num < num_of_clusters){
              likelihood[num] += normal_lpdf(theta_estimates[j] | cluster_center[num - 1], abs_sigma_hat);
            } else {
              likelihood[num] += student_t_lpdf(theta_estimates[j] | deg_freedom, 0, 1);
            }
          } else {
            likelihood[num] += normal_lpdf(theta_estimates[j] | cluster_center[num - 1], abs_sigma_hat);
          }
        }
      } else {
        if (num_junk == 1){
          if (num < num_of_clusters){
            likelihood[num] += normal_lpdf(theta_estimates[j] | cluster_center[num], abs_sigma_hat);
          } else {
            likelihood[num] += student_t_lpdf(theta_estimates[j] | deg_freedom, 0, 1);
          }
        } else {
          likelihood[num] += normal_lpdf(theta_estimates[j] | cluster_center[num], abs_sigma_hat);
        }
      }
    }

    log_lik[j] = log_sum_exp(likelihood);
    post_prob[j] = softmax(likelihood);
  }
}
