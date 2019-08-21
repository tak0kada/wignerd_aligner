#include "soft20_api.hpp"
#include "read_coef.hpp"
#include <iostream>
#include <cmath>

int main()
{
    std::vector<std::vector<double>> sig0 = read_coef("./cell00001.coef");
    std::vector<std::vector<double>> sig1 = read_coef("./cell00001.coef");
    soft20::Aligner aligner(std::sqrt(sig0[0].size()));
    const auto [dist, rot] = aligner.dist_rot(sig0, sig1);
    std::cout << rot[0] << " " << rot[1] << " " << rot[2] << std::endl;
}
