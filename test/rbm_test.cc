/*
Copyright 2016 WaizungTaam.  All rights reserved.

License:       Apache License 2.0
Email:         waizungtaam@gmail.com
Creation time: 2016-07-15
Last modified: 2016-07-15

*/

#include <iostream>
#include <iomanip>
#include "../src/rbm.h"
#include "../src/math/matrix.h"
#include "../src/math/utils.h"

int main() {
  int num_epochs = 1000, idx_epoch, num_vis = 8, 
  num_hid_1 = 6, num_steps_1 = 1, 
  num_hid_2 = 2, num_steps_2 = 1,
  num_hid_3 = 6, num_steps_3 = 8;
  double lr = 1e-4;
  Matrix data_train(100, num_vis, "binomial", 1.0, 0.5),
         data_test(4, num_vis, "binomial", 1.0, 0.5);
  RBM machine_1(num_vis, num_hid_1),
      machine_2(num_vis, num_hid_2),
      machine_3(num_vis, num_hid_3);
  for (idx_epoch = 0; idx_epoch < num_epochs; ++idx_epoch) {
    machine_1.train(data_train, lr, num_steps_1);
    machine_2.train(data_train, lr, num_steps_2);
    machine_3.train(data_train, lr, num_steps_3);
    std::cout << idx_epoch << "\t";
    std::cout << std::setprecision(8) << std::setw(10) << std::left << std::setfill('0') 
              << nn::pow(data_train - machine_1.reconstruct(data_train), 2).sum()
                 / data_train.shape()[0] / data_train.shape()[1] << "\t";
    std::cout << std::setprecision(8) << std::setw(10) << std::left << std::setfill('0') 
              << nn::pow(data_train - machine_2.reconstruct(data_train), 2).sum()
                 / data_train.shape()[0] / data_train.shape()[1] << "\t";
    std::cout << std::setprecision(8) << std::setw(10) << std::left << std::setfill('0') 
              << nn::pow(data_train - machine_3.reconstruct(data_train), 2).sum()
                 / data_train.shape()[0] / data_train.shape()[1] << "\n";                 
  }
  std::cout << data_test << std::endl;
  std::cout << machine_1.reconstruct(data_test) << std::endl;
  std::cout << machine_2.reconstruct(data_test) << std::endl;
  std::cout << machine_3.reconstruct(data_test) << std::endl;
  return 0;
}