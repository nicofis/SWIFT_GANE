/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/* This file's header */
#include "feedback.h"

/* Local includes. */
#include "../GANE/enrichment.h"
#include "../GANE/yield_tables.h"
#include "../GANE/energy_tables.h"
#include "hydro_properties.h"
#include "inline.h"
#include "random.h"
#include "timers.h"

/**
 * @brief Computes the fraction of the available SNII/HMXB energy to
 * inject for a given event [GANE].
 *
 * Note that the fraction can be > 1.
 *
 * This allows a choice of different f_th scaling functions.
 *
 * @param sp The #spart.
 * @param props The properties of the feedback model.
 */
double gane_feedback_energy_fraction(const struct spart *sp,
                                      const struct feedback_props *props) {
  /* Model parameter */
  const double f_E = props->f_E;
  return f_E;
}

/**
 * @brief Compute the properties of the SNII/HMXB feedback energy injection. [GANE]
 * Only does something if the particle reached the SNII or HMXB active age during this time
 * step.
 *
 * @param sp The star particle.
 * @param star_age Age of star at the beginning of the step in Gyrs.
 * @param dt Length of time-step in Gyrs.
 * @param num_gas_ngbs Total (integer) number of gas neighbours within the
 * star's kernel.
 * @param ngb_gas_mass Total un-weighted mass in the star's kernel (internal
 * units)
 * @param feedback_props The properties of the feedback model.
 */
INLINE static void compute_stellar_feedback(
    struct spart *sp, const double star_age, const double dt, const int ngb_gas_N,
    const float ngb_gas_mass, const struct feedback_props *feedback_props) {

  /* Properties of the model (all in internal units) */
  /* Get metallicity (metal mass fraction) at birth time of the star */
  const double Z_birth =
      chemistry_get_star_total_metal_mass_fraction_for_feedback(sp);
  
  double delta_E_SNII = 0.;
  double delta_E_HMXB = 0.;

  /* Compute SNII feedback, if corresponds */
  if (feedback_props->with_SNII_feedback) {
    const double E_start = interpolate_feedback_energy(star_age, Z_birth, 
                           feedback_props, gane_source_SNII);
    const double E_end = interpolate_feedback_energy(star_age + dt, Z_birth, 
                           feedback_props, gane_source_SNII);
    delta_E_SNII = (E_end - E_start) * sp->mass_init *
                    feedback_props->feedback_energy_to_internal;
  }

  /* Compute HMXB feedback, if corresponds */
  if (feedback_props->with_HMXB_feedback) {
    const double E_start = interpolate_feedback_energy(star_age, Z_birth, 
                           feedback_props, gane_source_HMXB);
    const double E_end = interpolate_feedback_energy(star_age + dt, Z_birth, 
                           feedback_props, gane_source_HMXB);
    delta_E_HMXB = (E_end - E_start) * sp->mass_init *
                    feedback_props->feedback_energy_to_internal;
  }
  
  const double delta_E_SNII_tot = delta_E_SNII + sp->feedback_data.delta_E_SNII_reservoir;
  const double delta_E_HMXB_tot = delta_E_HMXB + sp->feedback_data.delta_E_HMXB_reservoir;
  const double delta_E_tot = delta_E_SNII_tot + delta_E_HMXB_tot;
  const double f_E =
      gane_feedback_energy_fraction(sp, feedback_props);
  
  /* Update the time since the last SNII/HMXB feedback event */
  sp->feedback_data.time_since_last_feedback_event += dt;

  /* Abort if there is no SNII/HMXB energy available this step */
  if (f_E * delta_E_tot <= 0.) return;

  /* Number of SNII/HMXB events for this stellar particle */
  int number_of_feedback_events = 0;
  const char *injection_type = feedback_props->injection_type;

  if (strcmp(injection_type, "Uniform") == 0) {
    /* If we are doing uniform injection, we just inject the energy into all 
    neighbouring gas particles. */
    number_of_feedback_events = ngb_gas_N;
  }
  else if (strcmp(injection_type, "Stochastic") == 0) {
    /* If we are doing stochastic injection, we inject the energy into a 
    random subset of neighbouring gas particle equal to the number of SNII rays. */
    number_of_feedback_events = eagle_SNII_feedback_num_of_rays;
  }
  else {
    error("Invalid feedback injection type!");
    return;
  }

  const double delta_u = f_E * delta_E_tot * ngb_gas_N / (number_of_feedback_events * ngb_gas_mass);

  /* Decide whether to apply feedback based on minimum timestep or temperature, if greater than zero.
  If the star is old enough, we also apply feedback. */
  const bool do_feedback = 
  ((sp->feedback_data.time_since_last_feedback_event >= feedback_props->min_timestep_Gyr) 
  && (delta_u / feedback_props->temp_to_u_factor >= feedback_props->min_delta_T_K)) 
  || (star_age + dt >= 0.05);

  if (!do_feedback) {
    /* If we are not doing feedback, we store the energy in the reservoir for 
    future timesteps and abort. */
    sp->feedback_data.delta_E_SNII_reservoir += delta_E_SNII;
    sp->feedback_data.delta_E_HMXB_reservoir += delta_E_HMXB;
    return;
  }

  /* Current total f_E for this star */
  double star_f_E = sp->f_E * sp->number_of_feedback_events;

  /* New total */
  star_f_E = (star_f_E + f_E) / (sp->number_of_feedback_events + 1.);

  /* Store all of this in the star for delivery onto the gas and recording */
  sp->f_E = star_f_E;
  sp->number_of_feedback_events++;
  sp->feedback_data.to_distribute.delta_u = delta_u;
  sp->feedback_data.to_distribute.num_of_thermal_energy_inj =
      number_of_feedback_events;

  /* Track cumulative energy ejected by this star via SNII feedback [tracking] */
  sp->cumulative_SNII_energy_ejected += f_E * delta_E_SNII_tot;

  /* Track cumulative energy ejected by this star via HMXB feedback [tracking] */
  sp->cumulative_HMXB_energy_ejected += f_E * delta_E_HMXB_tot;

  /* Reset the reservoir since we have now injected the energy */
  sp->feedback_data.delta_E_SNII_reservoir = 0.f;
  sp->feedback_data.delta_E_HMXB_reservoir = 0.f;
  sp->feedback_data.time_since_last_feedback_event = 0.f;
}

/**
 * @brief calculates stellar mass in spart that died over the timestep, calls
 * functions to calculate feedback due to SNIa, SNII and AGB
 *
 * @param feedback_props feedback_props data structure
 * @param phys_const The physical constants in internal units.
 * @param cosmo The cosmological model.
 * @param sp spart that we're evolving
 * @param us unit_system data structure
 * @param age age of spart at beginning of step
 * @param dt length of current timestep
 * @param ti_begin The current integer time (for random number hashing).
 */
void compute_stellar_evolution(const struct feedback_props *feedback_props,
                               const struct phys_const *phys_const,
                               const struct cosmology *cosmo, struct spart *sp,
                               const struct unit_system *us, const double age,
                               const double dt, const integertime_t ti_begin) {

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (age < 0.f) error("Negative age for a star.");

  if (sp->count_since_last_enrichment != 0 && engine_current_step > 0)
    error("Computing feedback on a star that should not");
#endif

  /* Convert dt and stellar age from internal units to Gyr. */
  const double Gyr_inv = 1. / (phys_const->const_year * 1e9);
  const double dt_Gyr = dt * Gyr_inv;
  const double star_age_Gyr = age * Gyr_inv;

  /* Get the birth mass of the star */
  const double M_init = sp->mass_init;

  /* Get the total metallicity (metal mass fraction) at birth time and impose a
   * minimum */
  const double Z =
      max(chemistry_get_star_total_metal_mass_fraction_for_feedback(sp),
          exp10(log10_min_metallicity));

  /* Get the individual abundances (mass fractions at birth time) */
  const float *const abundances =
      chemistry_get_star_metal_mass_fraction_for_feedback(sp);

  /* Properties collected in the stellar density loop. */
  const int ngb_Number = sp->feedback_data.to_collect.ngb_N;
  const float ngb_gas_mass = sp->feedback_data.to_collect.ngb_mass;

  /* Check if there are neighbours, otherwise exit */
  if (ngb_gas_mass == 0.f || sp->density.wcount * pow_dimension(sp->h) < 1e-4) {
    feedback_reset_feedback(sp, feedback_props);
    return;
  }

  /* Update the enrichment weights */
  const float enrichment_weight_inv =
      sp->feedback_data.to_collect.enrichment_weight_inv;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_collect.enrichment_weight_inv < 0.)
    error("Negative inverse weight!");
#endif

  /* Now we start filling the data structure for information to apply to the
   * particles. Do _NOT_ read from the to_collect substructure any more. */

  /* Zero all the output fields */
  feedback_reset_feedback(sp, feedback_props);

  /* Update the weights used for distribution */
  const float enrichment_weight =
      (enrichment_weight_inv != 0.f) ? 1.f / enrichment_weight_inv : 0.f;
  sp->feedback_data.to_distribute.enrichment_weight = enrichment_weight;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_distribute.enrichment_weight < 0.)
    error("Negative weight!");
#endif

  /* Calculate mass of stars that has died from the star's birth up to the
   * beginning and end of timestep */
  const double max_dying_mass_Msun =
      dying_mass_msun(star_age_Gyr, Z, feedback_props);
  const double min_dying_mass_Msun =
      dying_mass_msun(star_age_Gyr + dt_Gyr, Z, feedback_props);

#ifdef SWIFT_DEBUG_CHECKS
  /* Sanity check. Worth investigating if necessary as functions for evaluating
   * mass of stars dying might be strictly decreasing.  */
  if (min_dying_mass_Msun > max_dying_mass_Msun)
    error("min dying mass is greater than max dying mass");
#endif

  /* Compute properties of the SNII/HMXB feedback model [GANE] */
  /* Check first if the star hasn't died yet */
  if ((feedback_props->with_SNII_feedback || feedback_props->with_HMXB_feedback) 
      && star_age_Gyr < 0.05) {
    compute_stellar_feedback(sp, star_age_Gyr, dt_Gyr, ngb_Number, 
                          ngb_gas_mass, feedback_props);
  }
  
  /* Integration interval is zero - this can happen if minimum and maximum
   * dying masses are above imf_max_mass_Msun. Return without doing any
   * enrichment. */
  if (min_dying_mass_Msun == max_dying_mass_Msun) return;

  /* Life is better in log */
  const double log10_max_dying_mass_Msun = log10(max_dying_mass_Msun);
  const double log10_min_dying_mass_Msun = log10(min_dying_mass_Msun);

  /* Compute elements, energy and momentum to distribute from the
   *  three channels SNIa, SNII, AGB */
  if (feedback_props->with_SNIa_enrichment) {
    evolve_SNIa(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun, M_init, Z,
                feedback_props, star_age_Gyr, dt_Gyr, &sp->feedback_data);
  }
  if (feedback_props->with_SNII_enrichment) {
    evolve_SNII(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun, M_init, Z,
                abundances, feedback_props, &sp->feedback_data);
  }
  if (feedback_props->with_AGB_enrichment) {
    evolve_AGB(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun, M_init, Z,
               abundances, feedback_props, &sp->feedback_data);
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_distribute.mass != 0.f)
    error("Injected mass will be lost");
#endif

  /* Compute the total mass to distribute (H + He  metals) */
  sp->feedback_data.to_distribute.mass =
      sp->feedback_data.to_distribute.total_metal_mass +
      sp->feedback_data.to_distribute.metal_mass[chemistry_element_H] +
      sp->feedback_data.to_distribute.metal_mass[chemistry_element_He];

  /* Compute energy change due to kinetic energy of ejectas */
  sp->feedback_data.to_distribute.energy +=
      sp->feedback_data.to_distribute.mass *
      feedback_props->AGB_ejecta_specific_kinetic_energy;

  /* Compute energy change due to kinetic energy of the star */
  sp->feedback_data.to_distribute.energy +=
      sp->feedback_data.to_distribute.mass * 0.5f *
      (sp->v[0] * sp->v[0] + sp->v[1] * sp->v[1] + sp->v[2] * sp->v[2]) *
      cosmo->a2_inv;

  TIMER_TOC(timer_do_star_evol);
}

/**
 * @brief Initialize the global properties of the feedback scheme.
 *
 * @param fp The #feedback_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
void feedback_props_init(struct feedback_props *fp,
                         const struct phys_const *phys_const,
                         const struct unit_system *us,
                         struct swift_params *params,
                         const struct hydro_props *hydro_props,
                         const struct cosmology *cosmo) {

  /* Main operation modes ------------------------------------------------- */

  fp->with_SNIa_feedback =
      parser_get_param_int(params, "GANEFeedback:use_SNIa_feedback");
  
  /* [GANE] */
  fp->with_SNII_feedback =
      parser_get_param_int(params, "GANEFeedback:use_SNII_feedback");

  fp->with_HMXB_feedback =
      parser_get_param_int(params, "GANEFeedback:use_HMXB_feedback");

  fp->with_AGB_enrichment =
      parser_get_param_int(params, "GANEFeedback:use_AGB_enrichment");

  fp->with_SNII_enrichment =
      parser_get_param_int(params, "GANEFeedback:use_SNII_enrichment");

  fp->with_SNIa_enrichment =
      parser_get_param_int(params, "GANEFeedback:use_SNIa_enrichment");

  if (fp->with_SNIa_feedback && !fp->with_SNIa_enrichment) {
    error("Cannot run with SNIa feedback without SNIa enrichment.");
  }

  /* Properties of the IMF model ------------------------------------------ */

  /* Minimal and maximal mass considered */
  fp->imf_max_mass_msun =
      parser_get_param_double(params, "GANEFeedback:IMF_max_mass_Msun");
  fp->imf_min_mass_msun =
      parser_get_param_double(params, "GANEFeedback:IMF_min_mass_Msun");

  /* Check that it makes sense. */
  if (fp->imf_max_mass_msun < fp->imf_min_mass_msun) {
    error("Can't have the max IMF mass smaller than the min IMF mass!");
  }

  fp->log10_imf_max_mass_msun = log10(fp->imf_max_mass_msun);
  fp->log10_imf_min_mass_msun = log10(fp->imf_min_mass_msun);

  /* Properties of the SNII energy feedback model ------------------------- */

  char model[64];
  parser_get_param_string(params, "GANEFeedback:SNII_feedback_model", model);
  if (strcmp(model, "Random") == 0)
    fp->feedback_model = SNII_random_ngb_model;
  else if (strcmp(model, "Isotropic") == 0)
    fp->feedback_model = SNII_isotropic_model;
  else if (strcmp(model, "MinimumDistance") == 0)
    fp->feedback_model = SNII_minimum_distance_model;
  else if (strcmp(model, "MinimumDensity") == 0)
    fp->feedback_model = SNII_minimum_density_model;
  else
    error(
        "The SNII feedback model must be either 'Random', 'MinimumDistance', "
        "'MinimumDensity' or 'Isotropic', not %s",
        model);

  /* Stellar mass limits for SNII [~feedback~ enrichment] [GANE] */
  const double SNII_min_mass_msun =
      parser_get_param_double(params, "GANEFeedback:SNII_min_mass_Msun");
  const double SNII_max_mass_msun =
      parser_get_param_double(params, "GANEFeedback:SNII_max_mass_Msun");

  /* Check that it makes sense. */
  if (SNII_max_mass_msun < SNII_min_mass_msun) {
    error("Can't have the max SNII mass smaller than the min SNII mass!");
  }

  fp->SNII_min_mass_msun = SNII_min_mass_msun;
  fp->SNII_max_mass_msun = SNII_max_mass_msun;
  fp->log10_SNII_min_mass_msun = log10(SNII_min_mass_msun);
  fp->log10_SNII_max_mass_msun = log10(SNII_max_mass_msun);

  /* Properties of the SNII/HMXB feedback model [GANE] --------------------------- */
  fp->f_E =
      parser_get_param_double(params, "GANEFeedback:SNII_HMXB_energy_fraction");
  parser_get_param_string(params, "GANEFeedback:SNII_HMXB_injection_type",
                        fp->injection_type);
  fp->min_timestep_Gyr =
      parser_get_opt_param_double(params, "GANEFeedback:SNII_HMXB_min_timestep_Gyr", 0.0);
  fp->min_delta_T_K =
      parser_get_opt_param_double(params, "GANEFeedback:SNII_HMXB_min_delta_T_K", 0.0);

  /* Properties of the SNII enrichment model -------------------------------- */

  /* Set factors for each element adjusting SNII yield */
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    char buffer[50];
    sprintf(buffer, "GANEFeedback:SNII_yield_factor_%s",
            chemistry_get_element_name((enum chemistry_element)elem));

    fp->SNII_yield_factor[elem] =
        parser_get_opt_param_float(params, buffer, 1.f);
  }

  /* Properties of the SNIa enrichment model -------------------------------- */

  fp->SNIa_DTD_delay_Gyr =
      parser_get_param_double(params, "GANEFeedback:SNIa_DTD_delay_Gyr");

  char temp[32] = {0};
  parser_get_param_string(params, "GANEFeedback:SNIa_DTD", temp);

  if (strcmp(temp, "Exponential") == 0) {

    fp->SNIa_DTD = eagle_feedback_SNIa_DTD_exponential;

    /* Read SNIa exponential DTD model parameters */
    fp->SNIa_DTD_exp_norm = parser_get_param_float(
        params, "GANEFeedback:SNIa_DTD_exp_norm_p_Msun");
    fp->SNIa_DTD_exp_timescale_Gyr = parser_get_param_float(
        params, "GANEFeedback:SNIa_DTD_exp_timescale_Gyr");
    fp->SNIa_DTD_exp_timescale_Gyr_inv = 1.f / fp->SNIa_DTD_exp_timescale_Gyr;

  } else if (strcmp(temp, "PowerLaw") == 0) {

    fp->SNIa_DTD = eagle_feedback_SNIa_DTD_power_law;

    /* Read SNIa power-law DTD model parameters */
    fp->SNIa_DTD_power_law_norm = parser_get_param_float(
        params, "GANEFeedback:SNIa_DTD_power_law_norm_p_Msun");

    /* Renormalize everything such that the integral converges to
       'SNIa_DTD_power_law_norm' over 13.6 Gyr. */
    fp->SNIa_DTD_power_law_norm /= log(13.6 / fp->SNIa_DTD_delay_Gyr);

  } else {
    error("Invalid SNIa DTD model: '%s'", temp);
  }

  /* Energy released by supernova type Ia */
  fp->E_SNIa_cgs =
      parser_get_param_double(params, "GANEFeedback:SNIa_energy_erg");
  fp->E_SNIa =
      fp->E_SNIa_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Properties of the SNIa enrichment model -------------------------------- */

  /* Read AGB ejecta velocity */
  const float ejecta_velocity_km_p_s = parser_get_param_float(
      params, "GANEFeedback:AGB_ejecta_velocity_km_p_s");

  /* Convert to internal units */
  const float ejecta_velocity_cgs = ejecta_velocity_km_p_s * 1e5;
  const float ejecta_velocity =
      ejecta_velocity_cgs / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);

  /* Convert to specific thermal energy */
  fp->AGB_ejecta_specific_kinetic_energy =
      0.5f * ejecta_velocity * ejecta_velocity;

  /* Properties of the enrichment down-sampling ----------------------------- */

  fp->stellar_evolution_age_cut =
      parser_get_param_double(params,
                              "GANEFeedback:stellar_evolution_age_cut_Gyr") *
      phys_const->const_year * 1e9;

  fp->stellar_evolution_sampling_rate = parser_get_param_double(
      params, "GANEFeedback:stellar_evolution_sampling_rate");

  if (fp->stellar_evolution_sampling_rate < 1 ||
      fp->stellar_evolution_sampling_rate >= (1 << (8 * sizeof(char) - 1)))
    error("Stellar evolution sampling rate too large. Must be >0 and <%d",
          (1 << (8 * sizeof(char) - 1)));

  /* Gather common conversion factors --------------------------------------- */

  /* Calculate internal mass to solar mass conversion factor */
  const double Msun_cgs = phys_const->const_solar_mass *
                          units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  const double unit_mass_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  fp->mass_to_solar_mass = unit_mass_cgs / Msun_cgs;
  fp->solar_mass_to_mass = 1. / fp->mass_to_solar_mass;

  /* Calculate temperature to internal energy conversion factor (all internal
   * units) */
  const double k_B = phys_const->const_boltzmann_k;
  const double m_p = phys_const->const_proton_mass;
  const double mu = hydro_props->mu_ionised;
  fp->temp_to_u_factor = k_B / (mu * hydro_gamma_minus_one * m_p);

  /* Calculate conversion factor from rho to n_H
   * Note this assumes primoridal abundance */
  const double X_H = hydro_props->hydrogen_mass_fraction;
  fp->rho_to_n_cgs =
      (X_H / m_p) * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  /* Initialise the IMF ------------------------------------------------- */

  init_imf(fp);

  /* Initialise the yields ---------------------------------------------- */

  /* Read yield table filepath  */
  parser_get_param_string(params, "GANEFeedback:yield_filename",
                          fp->yield_table_path);

  /* Allocate yield tables  */
  allocate_yield_tables(fp);

  /* Read the tables  */
  read_yield_tables(fp);

  /* Set yield_mass_bins array */
  const float imf_log10_mass_bin_size =
      (fp->log10_imf_max_mass_msun - fp->log10_imf_min_mass_msun) /
      (eagle_feedback_N_imf_bins - 1);

  for (int i = 0; i < eagle_feedback_N_imf_bins; i++)
    fp->yield_mass_bins[i] =
        imf_log10_mass_bin_size * i + fp->log10_imf_min_mass_msun;

  /* Resample yields from mass bins used in tables to mass bins used in IMF  */
  compute_yields(fp);

  /* Resample ejecta contribution to enrichment from mass bins used in tables to
   * mass bins used in IMF  */
  compute_ejecta(fp);

  /* Initialise the SNII/HMXB tables, if corresponds [GANE] -------------------------------------- */

  if (fp->with_SNII_feedback) {
    /* Read feedback energy table filepath  */
    parser_get_param_string(params, "GANEFeedback:SNII_filename",
                            fp->SNII_table_path);

    /* Allocate feedback energy tables */
    allocate_energy_tables(fp, gane_source_SNII);

    /* Read the tables  */
    read_energy_tables(fp, gane_source_SNII);
  }

  if (fp->with_HMXB_feedback) {
    /* Read feedback energy table filepath  */
    parser_get_param_string(params, "GANEFeedback:HMXB_filename",
                            fp->HMXB_table_path);

    /* Allocate feedback energy tables */
    allocate_energy_tables(fp, gane_source_HMXB);

    /* Read the tables  */
    read_energy_tables(fp, gane_source_HMXB);
  }

  /* Table energies are in erg / Msun */
  fp->feedback_energy_to_internal = fp->mass_to_solar_mass /
                                units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);
}

/**
 * @brief Clean-up the memory allocated for the feedback routines
 *
 * We simply free all the arrays.
 *
 * @param fp the feedback data structure.
 */
void feedback_clean(struct feedback_props *fp) {

  swift_free("imf-tables", fp->imf);
  swift_free("imf-tables", fp->imf_mass_bin);
  swift_free("imf-tables", fp->imf_mass_bin_log10);
  swift_free("feedback-tables", fp->yields_SNIa);
  swift_free("feedback-tables", fp->yield_SNIa_IMF_resampled);
  swift_free("feedback-tables", fp->yield_AGB.mass);
  swift_free("feedback-tables", fp->yield_AGB.metallicity);
  swift_free("feedback-tables", fp->yield_AGB.yield);
  swift_free("feedback-tables", fp->yield_AGB.yield_IMF_resampled);
  swift_free("feedback-tables", fp->yield_AGB.ejecta);
  swift_free("feedback-tables", fp->yield_AGB.ejecta_IMF_resampled);
  swift_free("feedback-tables", fp->yield_AGB.total_metals);
  swift_free("feedback-tables", fp->yield_AGB.total_metals_IMF_resampled);
  swift_free("feedback-tables", fp->yield_SNII.mass);
  swift_free("feedback-tables", fp->yield_SNII.metallicity);
  swift_free("feedback-tables", fp->yield_SNII.yield);
  swift_free("feedback-tables", fp->yield_SNII.yield_IMF_resampled);
  swift_free("feedback-tables", fp->yield_SNII.ejecta);
  swift_free("feedback-tables", fp->yield_SNII.ejecta_IMF_resampled);
  swift_free("feedback-tables", fp->yield_SNII.total_metals);
  swift_free("feedback-tables", fp->yield_SNII.total_metals_IMF_resampled);
  if (fp->with_SNII_feedback) {
    swift_free("feedback-tables", fp->SNII_energies.age); // [GANE]
    swift_free("feedback-tables", fp->SNII_energies.metallicity); // [GANE]
  }
  if (fp->with_HMXB_feedback) {
    swift_free("feedback-tables", fp->HMXB_energies.age); // [GANE]
    swift_free("feedback-tables", fp->HMXB_energies.metallicity); // [GANE]
  }
  swift_free("feedback-tables", fp->lifetimes.mass);
  swift_free("feedback-tables", fp->lifetimes.metallicity);
  swift_free("feedback-tables", fp->yield_mass_bins);
  for (int i = 0; i < eagle_feedback_lifetime_N_metals; i++) {
    free(fp->lifetimes.dyingtime[i]);
  }
  free(fp->lifetimes.dyingtime);
  for (int i = 0; i < eagle_feedback_SNIa_N_elements; i++) {
    free(fp->SNIa_element_names[i]);
  }
  free(fp->SNIa_element_names);
  for (int i = 0; i < eagle_feedback_SNII_N_elements; i++) {
    free(fp->SNII_element_names[i]);
  }
  free(fp->SNII_element_names);
  for (int i = 0; i < eagle_feedback_AGB_N_elements; i++) {
    free(fp->AGB_element_names[i]);
  }
  free(fp->AGB_element_names);
  /* [GANE] */
  if (fp->with_SNII_feedback) {
    for (int i = 0; i < gane_feedback_N_metals; i++) {
      free(fp->SNII_energies.energy[i]);
    }
    free(fp->SNII_energies.energy);
  }
  if (fp->with_HMXB_feedback) {
    for (int i = 0; i < gane_feedback_N_metals; i++) {
      free(fp->HMXB_energies.energy[i]);
    }
    free(fp->HMXB_energies.energy);
  }
}

/**
 * @brief Write a feedback struct to the given FILE as a stream of bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_dump(const struct feedback_props *feedback, FILE *stream) {

  /* To make sure everything is restored correctly, we zero all the pointers to
     tables. If they are not restored correctly, we would crash after restart on
     the first call to the feedback routines. Helps debugging. */
  struct feedback_props feedback_copy = *feedback;

  /* zero AGB and SNII table pointers */
  zero_yield_table_pointers(&feedback_copy.yield_AGB);
  zero_yield_table_pointers(&feedback_copy.yield_SNII);

  /* zero SNIa table pointers */
  feedback_copy.yield_SNIa_IMF_resampled = NULL;
  feedback_copy.yields_SNIa = NULL;
  feedback_copy.yield_SNIa_total_metals_IMF_resampled = 0;

  /* zero element name tables */
  feedback_copy.SNIa_element_names = NULL;
  feedback_copy.SNII_element_names = NULL;
  feedback_copy.AGB_element_names = NULL;

  /* zero SNII energy table pointers [GANE] */
  feedback_copy.SNII_energies.age = NULL;
  feedback_copy.SNII_energies.metallicity = NULL;
  feedback_copy.SNII_energies.energy = NULL;

  /* zero HMXB energy table pointers [GANE] */
  feedback_copy.HMXB_energies.age = NULL;
  feedback_copy.HMXB_energies.metallicity = NULL;
  feedback_copy.HMXB_energies.energy = NULL;

  /* zero mass bins table */
  feedback_copy.yield_mass_bins = NULL;

  /* zero lifetime tracks */
  feedback_copy.lifetimes.mass = NULL;
  feedback_copy.lifetimes.metallicity = NULL;
  feedback_copy.lifetimes.dyingtime = NULL;

  /* zero IMF tables */
  feedback_copy.imf = NULL;
  feedback_copy.imf_mass_bin = NULL;
  feedback_copy.imf_mass_bin_log10 = NULL;

  restart_write_blocks((void *)&feedback_copy, sizeof(struct feedback_props), 1,
                       stream, "feedback", "feedback function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * Read the structure from the stream and restore the feedback tables by
 * re-reading them.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_restore(struct feedback_props *feedback, FILE *stream) {
  restart_read_blocks((void *)feedback, sizeof(struct feedback_props), 1,
                      stream, NULL, "feedback function");

  if (strlen(feedback->yield_table_path) != 0)
    feedback_restore_tables(feedback);

  /* [GANE] */
  if (feedback->with_SNII_feedback &&
    strlen(feedback->SNII_table_path) != 0)
    feedback_restore_energy_tables(feedback, gane_source_SNII);
  
  if (feedback->with_HMXB_feedback &&
    strlen(feedback->HMXB_table_path) != 0)
    feedback_restore_energy_tables(feedback, gane_source_HMXB);
}
