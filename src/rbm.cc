/*
Copyright 2016 WaizungTaam.  All rights reserved.

License:       Apache License 2.0
Email:         waizungtaam@gmail.com
Creation time: 2016-07-15
Last modified: 2016-07-15

*/

#include "rbm.h"
#include "./math/matrix.h"
#include "./math/utils.h"

RBM::RBM(int num_vis, int num_hid) {
  weight = Matrix(num_vis, num_hid, "uniform", -1.0, 1.0);
  w_bias_vis = Matrix(1, num_vis, "uniform", -1.0, 1.0);
  w_bias_hid = Matrix(1, num_hid, "uniform", -1.0, 1.0);
}
Matrix RBM::reconstruct(const Matrix & mat_vis) {
  return prop_hv(prop_vh(mat_vis));
}
void RBM::train(const Matrix & mat_vis, 
                double learning_rate, int num_steps) {
  contrastive_divergence(mat_vis, learning_rate, num_steps);
}
Matrix RBM::prop_vh(const Matrix & mat_vis) {
  return nn::logistic(mat_vis * weight + w_bias_hid);
}
Matrix RBM::prop_hv(const Matrix & mat_hid) {
  return nn::logistic(mat_hid * weight.T() + w_bias_vis);
}
Matrix RBM::sample_vh(const Matrix & mat_vis) {
  return nn::binomial_sample(1, prop_vh(mat_vis));
}
Matrix RBM::sample_hv(const Matrix & mat_hid) {
  return nn::binomial_sample(1, prop_hv(mat_hid));
}
Matrix RBM::gibbs_vhv(const Matrix & mat_vis) {
  return sample_hv(sample_vh(mat_vis));
}
Matrix RBM::gibbs_hvh(const Matrix & mat_hid) {
  return sample_vh(sample_hv(mat_hid));
}
void RBM::contrastive_divergence(const Matrix & mat_vis, 
                                 double learning_rate, 
                                 int num_steps) {
  int idx_step;
  Matrix prob_hid_begin = prop_vh(mat_vis);
  Matrix samp_hid = sample_vh(mat_vis);
  for (idx_step = 1; idx_step < num_steps - 1; ++idx_step) {
    samp_hid = gibbs_hvh(samp_hid);
  }
  Matrix samp_vis_end = sample_hv(samp_hid);
  Matrix prob_hid_end = prop_vh(samp_vis_end);

  Matrix bias = Matrix(mat_vis.shape()[0], 1, 1.0);

  Matrix delta_weight = mat_vis.T() * prob_hid_begin - 
                        samp_vis_end.T() * prob_hid_end;
  Matrix delta_w_bias_vis = bias.T() * (mat_vis - samp_vis_end);
  Matrix delta_w_bias_hid = bias.T() * (prob_hid_begin - prob_hid_end);

  weight += learning_rate * delta_weight;
  w_bias_vis += learning_rate * delta_w_bias_vis;
  w_bias_hid += learning_rate * delta_w_bias_hid;
} 