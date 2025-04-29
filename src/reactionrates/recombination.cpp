// SPDX-License-Identifier: MIT

#include <pdi.h>

#include "geometry.hpp"
#include "recombination.hpp"

RecombinationRate::RecombinationRate(double const n_0, double const T_0, double const K_cx_0)
    : m_slope_coefficient("slope_coefficient")
    , m_intercept_coefficient("intercept_coefficient")
    , m_n_0(n_0)
    , m_T_0(T_0)
    , m_K_cx_0(K_cx_0)
{
    double slope_coefficient_alloc[s_coefficients_size]
            = {0.03849324670136343, -0.03713312941661376};
    double intercept_coefficient_alloc[s_coefficients_size]
            = {-5.117013348529228, -1.1497231886316353};

    DKokkosView_h<s_coefficients_size> slope_coefficient_host(slope_coefficient_alloc);
    DKokkosView_h<s_coefficients_size> intercept_coefficient_host(intercept_coefficient_alloc);
    Kokkos::deep_copy(m_slope_coefficient, slope_coefficient_host);
    Kokkos::deep_copy(m_intercept_coefficient, intercept_coefficient_host);

    PDI_multi_expose(
            "r_rate_coeff_pol_expose",
            "recombination_slope_coefficients",
            slope_coefficient_host.data(),
            PDI_OUT,
            "recombination_intercept_coefficients",
            intercept_coefficient_host.data(),
            PDI_OUT,
            NULL);
}

RecombinationRate::RecombinationRate(double const norm_coeff_rate)
    : RecombinationRate(1e20, 10, 1 / norm_coeff_rate)
{
}

DFieldSpX RecombinationRate::operator()(
        DFieldSpX rate,
        DConstFieldSpX density,
        DConstFieldSpX temperature) const
{
    DKokkosView<s_coefficients_size> slope_coeff_proxy(m_slope_coefficient);
    DKokkosView<s_coefficients_size> intercept_coeff_proxy(m_intercept_coefficient);
    double T_0_proxy(m_T_0);
    double n_0_proxy(m_n_0);
    double k_cx_0_proxy(m_K_cx_0);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(rate),
            KOKKOS_LAMBDA(IdxSpX const ifspx) {
                IdxX idx(ifspx);
                // we convert the temperature from a normalisation to T_0 to a normalisation
                // on T=10eV, on which the log approximation has been made
                double temp = temperature(ielec(), idx) * T_0_proxy / 10.;
                double temperature_log10 = Kokkos::log10(temp);
                // same with n_0
                double dens = density(ielec(), idx) * n_0_proxy / 1e20;
                double density_log10 = Kokkos::log10(dens);

                double rate_log10 = 0.0;
                for (int i = 0; i < slope_coeff_proxy.size(); ++i) {
                    double polynomial_cs
                            = slope_coeff_proxy(i) * density_log10 + intercept_coeff_proxy(i);
                    rate_log10 += polynomial_cs * Kokkos::pow(temperature_log10, i);
                }
                // normalise the result to K_cx_0
                rate(ifspx) = Kokkos::pow(10, rate_log10) / k_cx_0_proxy;
            });
    return rate;
}
