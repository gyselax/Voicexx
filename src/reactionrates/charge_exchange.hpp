// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "ireactionrate.hpp"
#include "view.hpp"

/**
 * @brief A class that describes a charge-exchange reaction rate.
 */
class ChargeExchangeRate : public IReactionRate
{
private:
    /*
     * @brief Stores the linear coefficients to compute the polynomial for the charge-exchange reaction rate vs. ion temperature.
     * n_0, T_0 and K_cx_0 are the normalisation choices. They are used to transform normalised
     * temperatures and densities to real ones because the reaction rates depend on the real ones.
     * As usual T is wrongly called the temperature. It is actually the energy.
    */
    static constexpr std::size_t s_coefficients_size = 5;
    DKokkosView<s_coefficients_size> m_cx_coefficients;

    double const m_n_0;
    double const m_T_0;
    double m_K_cx_0; // this K_cx_0 is normalised to 1e-14 for numerical reasons

public:
    /**
     * @brief Creates an instance of the ConstantChargeExchangeRate class.
     * A polynomial of reaction rate vs. ion temperature is fitted in log-log space.
     * The polynomial fits data from the OPEN-ADAS database for Hydrogen.
     * m_K_cx_0 is computed by taking a typical charge-exchange reaction rate
     * for the given density and normalisation. This is equivalent to use the operator()
     * of this class with n_0 and T_0.
     *
     * @param[in] n_0 The normalisation density
     * @param[in] T_0 The normalisation temperature
     */
    ChargeExchangeRate(double n_0, double T_0);

    /**
     * @brief Creates an instance of the ConstantChargeExchangeRate class by specifying
     * a charge-exchange reaction rate normalisation.
     * This allow to change the K_cx_0 for given n_0 and T_0
     * This method should not be used, it is only for testing!
     *
     * @param[in] n_0 The normalisation density
     * @param[in] T_0 The normalisation temperature
     * @param[in] K_cx_0 The normalisation charge-exchange reaction rate
     */
    ChargeExchangeRate(double n_0, double T_0, double K_cx_0);

    /**
     * @brief Creates an instance of the ConstantChargeExchangeRate class by specifying
     * a normalisation coefficient.
     * This constructor is used to conserve retro-compatibility. The normalisation density
     * is chosen to be 1e20, the normalisation temperature is 10eV and the K_cx_0 is
     * 1e-14 divided by this normalisation coefficient.
     *
     * @param[in] norm_coeff_rate The normalisation coefficient
     */
    ChargeExchangeRate(double const norm_coeff_rate);


    /**
     * @brief Return the m_K_cx_0. This is used to initialised the ionisation and recombination
     * reaction rate to ensure consistency
     */
    [[nodiscard]] double get_Kcx0() const;

    /**
     * @brief Compute the charge-exchange reaction rate depending on density and temperature.
     *
     * @param[out] rate The charge-exchange reaction rates.
     * @param[in] density The normalised plasma density at which the reaction rate is computed. Although inputted, charge-exchange reaction rate is not affected by density.
     * @param[in] temperature The plasma normalised temperature at which the reaction rate is computed.
     * @return The charge-exchange reaction rate.
     */
    DFieldSpX operator()(DFieldSpX rate, DConstFieldSpX density, DConstFieldSpX temperature)
            const override;
};
