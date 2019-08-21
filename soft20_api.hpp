extern "C"
{
#include "soft20_wrap.h"
}
#include <vector>
#include <array>
#include <utility>
#include <cassert>
#include <iostream>

namespace soft20
{

struct Aligner
{
private:
    workspace w;
    int bwIn;

public:
    Aligner(const int bwIn)
    : bwIn(bwIn)
    {
        reserve(&this->w, bwIn, bwIn);
    }
    ~Aligner()
    {
        clear(&(this->w));
    }

    double distance(const std::vector<std::vector<double>>& sig0, const std::vector<std::vector<double>>& sig1)
    {
        assert(sig0[0].size() == sig1[0].size());
        assert(sig0[0].size() == static_cast<std::size_t>(bwIn * bwIn));
        double alpha, beta, gamma;
        double dist;
        calc_rot_dist_c(sig0[0].size(), bwIn,
                sig0[0].data(), sig1[0].data(),
                sig0[1].data(), sig1[1].data(),
                sig0[2].data(), sig1[2].data(),
                &alpha, &beta, &gamma, &dist, &this->w);
            // (sig0.size(), bwOut, sig1.size(), sig1.data(), bwOut, corr.data(), &(this->w));
        return dist;
    }

    std::pair<double, std::array<double, 3>> dist_rot(
            const std::vector<std::vector<double>>& sig0, const std::vector<std::vector<double>>& sig1)
    {
        assert(sig0[0].size() == sig1[0].size());
        assert(sig0[0].size() == static_cast<std::size_t>(bwIn * bwIn));
        double alpha, beta, gamma;
        double dist;
        calc_rot_dist_c(sig0[0].size(), bwIn,
                sig0[0].data(), sig1[0].data(),
                sig0[1].data(), sig1[1].data(),
                sig0[2].data(), sig1[2].data(),
                &alpha, &beta, &gamma, &dist, &this->w);
        return {dist, {alpha, beta, gamma}};
    }
};

} // end of namespace soft20
