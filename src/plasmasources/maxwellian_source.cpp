#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "mask_tanh.hpp"
#include "maxwellian_source.hpp"
#include "paraconfpp.hpp"
#include "species_info.hpp"

MaxwellianSource::MaxwellianSource(
        IdxRangeX const& idx_range_x,
        IdxRangeVx const& idx_range_vx,
        IdxRangeSp const& idx_range_species,
        double const extent,
        double const stiffness,
        double const amplitude,
        double const density,
        double const temperature_elec,
        double const temperature_ions)
    : m_amplitude(amplitude)
    , m_density(density)
    , m_temperature_elec(temperature_elec)
    , m_temperature_ions(temperature_ions)
    , m_spatial_extent(idx_range_x)
    , m_velocity_shape(IdxRangeSpVx(idx_range_species, idx_range_vx))
{
    host_t<DFieldMemX> spatial_extent_host
            = mask_tanh(idx_range_x, extent, stiffness, MaskType::Normal, true);
    ddc::parallel_deepcopy(get_field(m_spatial_extent), get_const_field(spatial_extent_host));

    // compute the source velocity profile (maxwellian profile)
    double const coeff_elec(1.0 / std::sqrt(2 * M_PI * m_temperature_elec));
    double const coeff_ions(1.0 / std::sqrt(2 * M_PI * m_temperature_ions));
    IdxRangeSpVx idx_range_spvx(idx_range_species, idx_range_vx);
    host_t<DFieldMemSpVx> velocity_shape_host(idx_range_spvx);
    ddc::host_for_each(idx_range_spvx, [&](IdxSpVx const ispvx) {
        IdxVx ivx(ispvx);
        IdxSp isp(ispvx);
        CoordVx const coordvx = ddc::coordinate(ivx);
        double const coordvx_sq = coordvx * coordvx;
        double source;
        if (isp == ielec()) {
            source = coeff_elec * std::exp(-coordvx_sq / (2 * temperature_elec));
        } else {
            source = coeff_ions * std::exp(-coordvx_sq / (2 * temperature_ions));
        }
        velocity_shape_host(ispvx) = density * source;
    });
    ddc::parallel_deepcopy(get_field(m_velocity_shape), velocity_shape_host);
    ddc::expose_to_pdi("maxwellian_source_extent", extent);
    ddc::expose_to_pdi("maxwellian_source_stiffness", stiffness);
    ddc::expose_to_pdi("maxwellian_source_amplitude", m_amplitude);
    ddc::expose_to_pdi("maxwellian_source_density", m_density);
    ddc::expose_to_pdi("maxwellian_source_temperature_elec", m_temperature_elec);
    ddc::expose_to_pdi("maxwellian_source_temperature_ions", m_temperature_ions);
    ddc::expose_to_pdi("maxwellian_source_velocity_shape", velocity_shape_host);
    ddc::expose_to_pdi("maxwellian_source_spatial_extent", spatial_extent_host);
}

DFieldSpXVx MaxwellianSource::operator()(DFieldSpXVx const allfdistribu, double const dt) const
{
    Kokkos::Profiling::pushRegion("MaxwellianSource");
    DConstFieldSpVx velocity_shape = get_const_field(m_velocity_shape);

    DConstFieldX spatial_extent = get_const_field(m_spatial_extent);

    double const amplitude = m_amplitude;

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(allfdistribu),
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                double const df(
                        amplitude * spatial_extent(ddc::select<GridX>(ispxvx))
                        * velocity_shape(ddc::select<Species, GridVx>(ispxvx)) * dt);
                allfdistribu(ispxvx) += df;
            });

    Kokkos::Profiling::popRegion();
    return allfdistribu;
}

MaxwellianSource maxwellian_source::init_from_input(
        IdxRangeSpXVx grid_idx_range,
        PC_tree_t const& yaml_input_file)
{
    IdxRangeX mesh_x(grid_idx_range);
    IdxRangeVx mesh_vx(grid_idx_range);
    IdxRangeSp mesh_sp(grid_idx_range);
    return MaxwellianSource(
            mesh_x,
            mesh_vx,
            mesh_sp,
            PCpp_double(yaml_input_file, ".PlasmaSource.extent"),
            PCpp_double(yaml_input_file, ".PlasmaSource.stiffness"),
            PCpp_double(yaml_input_file, ".PlasmaSource.amplitude"),
            PCpp_double(yaml_input_file, ".PlasmaSource.density"),
            PCpp_double(yaml_input_file, ".PlasmaSource.temperature_elec"),
            PCpp_double(yaml_input_file, ".PlasmaSource.temperature_ions"));
}
