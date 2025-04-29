// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "ireactionrate.hpp"
#include "view.hpp"

/**
 * @brief A class that describes a ionisation reaction rate.
 */
class IonisationRate : public IReactionRate
{
private:
    /**
     * @brief Stores the coefficients depending on electron density for calculating the ionisation reaction rate polynomial as a function of electron temperature.
     * n_0, T_0 and K_cx_0 are the normalisation choices. They are used to transform normalised
     * temperatures and densities to real ones because the reaction rates depend on the real ones.
     */
    static constexpr std::size_t s_coefficients_size = 6;
    DKokkosView<s_coefficients_size> m_slope_coefficient;
    DKokkosView<s_coefficients_size> m_intercept_coefficient;

    double const m_n_0;
    double const m_T_0;
    double const m_K_cx_0; // normalised to 1e-14 for numerical reasons

public:
    /**
     * @brief Creates an instance of the ConstantIonisationRate class.
     * A polynomial of reaction rate vs. ion temperature is fitted in log-log space.
     * The polynomial fits data from the OPEN-ADAS database for Hydrogen.
     *
     * @param[in] n_0 The normalisation density
     * @param[in] T_0 The normalisation temperature
     * @param[in] K_cx_0 The normalisation charge_exchange rate.
     * It is always given by the K_cx_0 of a ChargeExchangeRate.
     */
    explicit IonisationRate(double n_0, double T_0, double K_cx_0);

    /**
     * @brief Creates an instance of the ConstantIonisationRate class.
     * a normalisation coefficient.
     * This constructor is used to conserve retro-compatibility. The normalisation density
     * is chosen to be 1e20, the normalisation temperature is 10eV and the K_cx_0 is
     * 1e-14 divided by this normalisation coefficient.
     *
     * @param[in] norm_coeff_rate The normalisation coefficient
     */
    explicit IonisationRate(double const norm_coeff_rate);

    /**
     * @brief Compute the ionisation reaction rate depending on density and temperature.
     *
     * @param[out] rate The ionisation reaction rates.
     * @param[in] density The normalised plasma density at which the reaction rate is computed.
     * @param[in] temperature The normalised plasma temperature at which the reaction rate is computed.
     * @return The ionisation reaction rate.
     */
    DFieldSpX operator()(DFieldSpX rate, DConstFieldSpX density, DConstFieldSpX temperature)
            const override;
};
