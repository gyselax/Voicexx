// SPDX-License-Identifier: MIT

#include <pdi.h>

#include "charge_exchange.hpp"
#include "geometry.hpp"
#include "species_info.hpp"

ChargeExchangeRate::ChargeExchangeRate(double const n_0, double const T_0)
    : m_cx_coefficients("chargexchange_coefficients")
    , m_n_0(n_0)
    , m_T_0(T_0)
    , m_K_cx_0(1.)
{
    double cx_coefficients_alloc[s_coefficients_size]
            = {0.25125030648153424,
               0.37349850950319474,
               -0.009966677207985306,
               0.017922185310747015,
               -0.010699854255929037};
    DKokkosView_h<s_coefficients_size> cx_coefficients_host(cx_coefficients_alloc);
    Kokkos::deep_copy(m_cx_coefficients, cx_coefficients_host);

    // we construct the K_cx_0 which depends on n_0 and T_0
    // this K_cx_0 is actually computed as in the operator() of this class,
    // see the comments below for more details
    double temperature_log10 = Kokkos::log10(T_0 / 10);
    double rate_log10 = 0.0;
    for (int i = 0; i < cx_coefficients_host.extent(0); ++i) {
        rate_log10 += cx_coefficients_host[i] * Kokkos::pow(temperature_log10, i);
    }
    m_K_cx_0 = pow(10, rate_log10);

    PDI_multi_expose(
            "cx_rate_coeff_pol_expose",
            "charge_exchange_coefficients",
            cx_coefficients_host.data(),
            PDI_OUT,
            NULL);
}

ChargeExchangeRate::ChargeExchangeRate(double const n_0, double const T_0, double const K_cx_0)
    : ChargeExchangeRate(n_0, T_0)
{
    m_K_cx_0 = K_cx_0;
}

ChargeExchangeRate::ChargeExchangeRate(double const norm_coeff_rate)
    : ChargeExchangeRate(1e20, 10, 1 / norm_coeff_rate)
{
}

double ChargeExchangeRate::get_Kcx0() const
{
    return m_K_cx_0;
}

DFieldSpX ChargeExchangeRate::operator()(
        DFieldSpX rate,
        DConstFieldSpX density,
        DConstFieldSpX temperature) const
{
    DKokkosView<s_coefficients_size> cx_coeff_proxy(m_cx_coefficients);
    double k_cx_0_proxy(m_K_cx_0);
    double T_0_proxy(m_T_0);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(rate),
            KOKKOS_LAMBDA(IdxSpX const ifspx) {
                // we convert the temperature from a normalisation to T_0 to a normalisation
                // on T=10eV, on which the polynomial fitting (on log(n),log(T)) has been made
                double temp = temperature(ielec(), ddc::select<GridX>(ifspx)) * T_0_proxy / 10.;
                double temperature_log10 = Kokkos::log10(temp);

                double rate_log10 = 0.0;
                for (int i = 0; i < cx_coeff_proxy.size(); ++i) {
                    rate_log10 += cx_coeff_proxy(i) * Kokkos::pow(temperature_log10, i);
                }
                // normalise the result to K_cx_0
                rate(ifspx) = Kokkos::pow(10, rate_log10) / k_cx_0_proxy;
            });
    return rate;
}
