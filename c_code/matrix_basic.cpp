// SPDX-FileCopyrightText: 2021 Benjamin Brock
//
// SPDX-License-Identifier: BSD-3-Clause

#include <bcl/bcl.hpp>
#include <bcl/containers/DMatrix.hpp>

float multiply_by_two(float x) {
  return 2 * x;
}

int main(int argc, char** argv) {
  BCL::init();

  BCL::DMatrix<float> a({8, 8});
  BCL::DMatrix<float> b({8, 8});

  a = 1;
  b = 1;

  BCL::print("A matrix:\n");
  a.print();

  a.apply_inplace(multiply_by_two);

  BCL::print("A matrix after applying \"multiply_by_two\":\n");
  a.print();

  b.apply_inplace([](float x) { return x + 12; });

  auto c = a.dot(b);

  BCL::print("Result of multiply A and B:\n");
  c.print();

  BCL::finalize();
  return 0;
}
