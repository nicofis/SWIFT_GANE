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
#ifndef SWIFT_GANE_FEEDBACK_PROPERTIES_THERMAL_H
#define SWIFT_GANE_FEEDBACK_PROPERTIES_THERMAL_H

/* Config parameters. */
#include <config.h>

/* Local includes. */
#include "chemistry.h"
#include "hydro_properties.h"

/**
 * @brief Modes of energy injection for SNII feedback
 */
enum SNII_feedback_models {
  SNII_random_ngb_model,       /*< Random neighbour model for SNII feedback */
  SNII_isotropic_model,        /*< Isotropic model of SNII feedback */
  SNII_minimum_distance_model, /*< Minimum-distance model of SNII feedback */
  SNII_minimum_density_model   /*< Minimum-density model of SNII feedback */
};

/**
 * @brief Stores AGB and SNII yield tables
 */
struct yield_table {

  /*! Yield table mass bins */
  double *mass;

  /*! Yield table metallicity (metal mass fractions) bins */
  double *metallicity;

  /*! Array to store yield table (individual metals produced by the star)
     resampled by IMF mass bins */
  double *yield_IMF_resampled;

  /*! Array to store yield table being read in */
  double *yield;

  /*! Array to store table of ejecta (metals alredy in the stars that are
    ejected) resampled by IMF mass bins */
  double *ejecta_IMF_resampled;

  /*! Array to store table of ejecta being read in */
  double *ejecta;

  /*! Array to store table of total mass released ( metals produced by the star)
    resampled by IMF mass bins */
  double *total_metals_IMF_resampled;

  /* Array to store table of total mass released being read in */
  double *total_metals;
};

/**
 * @brief Stores tables to determine stellar lifetimes. Used for calculation of
 * IMF
 */
struct lifetime_table {

  /* table of masses */
  double *mass;

  /* table of metallicities */
  double *metallicity;

  /* table of lifetimes depending on mass an metallicity */
  double **dyingtime;
};

/**
 * @brief Functional form of the SNIa delay time distribution.
 */
enum eagle_feedback_SNIa_DTD {

  /*! Power-law with slope -1 */
  eagle_feedback_SNIa_DTD_power_law = 1,

  /*! Exponential model (EAGLE default) */
  eagle_feedback_SNIa_DTD_exponential = 2
};

/**
 * @brief Stores tables to determine SNII/HMXB feedback energies [GANE].
 */
struct energy_table {

  /* table of metallicities */
  double *metallicity;

  /* table of ages */
  double *age;

  /* table of energies depending on ages and metallicity */
  double **energy;
};

/**
 * @brief Properties of the GANE feedback model.
 */
struct feedback_props {

  /* ------------ Main operation modes ------------- */

  /*! Are we doing AGB enrichment? */
  int with_AGB_enrichment;

  /*! Are we doing SNII enrichment? */
  int with_SNII_enrichment;

  /*! Are we doing SNIa enrichment? */
  int with_SNIa_enrichment;

  /*! Are we doing SNIa feedback? */
  int with_SNIa_feedback;

  /*! Are we doing SNII/HMXB feedback? [GANE] */
  int with_SNII_HMXB_feedback;

  /* ------------ Yield tables    ----------------- */

  /* Yield tables for AGB and SNII  */
  struct yield_table yield_AGB;
  struct yield_table yield_SNII;

  /* Arrays of yield tables for SNIa */
  double *yield_SNIa_IMF_resampled;
  double yield_SNIa_total_metals_IMF_resampled;
  double *yields_SNIa;

  /* Arrays for names of elements being tracked for each enrichment channel */
  char **SNIa_element_names;
  char **SNII_element_names;
  char **AGB_element_names;

  /* Array of mass bins for yield calculations */
  double *yield_mass_bins;

  /* Location of yield tables */
  char yield_table_path[200];

  /* ------------- Lifetime tracks   --------------- */

  /* Table of lifetime values */
  struct lifetime_table lifetimes;

  /* ------------- SNII/HMXB feedback energy table [GANE] --------------- */

  /* Table of SNII/HMXB feedback energy values */
  struct energy_table feedback_energies;

  /* Location of HMXB tables */
  char energy_table_path[200];

  /* ------------- SNII parameters    --------------- */

  /* Array of adjustment factors for SNII  */
  float SNII_yield_factor[chemistry_element_count];

  /* ------------- SNIa parameters    --------------- */

  /* What delay time distribution are we using? */
  enum eagle_feedback_SNIa_DTD SNIa_DTD;

  /*! Normalisation of the SNIa DTD in the exponential model */
  float SNIa_DTD_exp_norm;

  /*! Time-scale of the SNIa decay function in the exponential model in
   * Giga-years */
  float SNIa_DTD_exp_timescale_Gyr;

  /*! Inverse of time-scale of the SNIa decay function in the exponential model
   * in Giga-years */
  float SNIa_DTD_exp_timescale_Gyr_inv;

  /*! Normalisation of the SNIa DTD in the power-law model */
  float SNIa_DTD_power_law_norm;

  /*! Stellar age below which no SNIa explode in Giga-years */
  float SNIa_DTD_delay_Gyr;

  /*! Energy released by one supernova type II in cgs units */
  double E_SNIa_cgs;

  /*! Energy released by one supernova type II in internal units */
  float E_SNIa;

  /* ------------- AGB parameters    ---------------- */

  /*! Specific kinetic energy injected from AGB ejectas (in internal units). */
  float AGB_ejecta_specific_kinetic_energy;

  /* ------------- SNII/HMXB parameters [GANE] ---------------- */

  /*! Energy fraction for SNII/HMXB feedback */
  float f_E;

  /*! Type of feedback energy injection */
  char injection_type[20];

  /*! Minimum timestep for SNII/HMXB feedback */
  float min_timestep_Gyr;

  /*! Minimum temperature change for SNII/HMXB feedback */
  float min_delta_T_K;

  /* ------------- Conversion factors --------------- */

  /*! Conversion factor from internal mass unit to solar mass */
  double mass_to_solar_mass;

  /*! Conversion factor from internal mass unit to solar mass */
  double solar_mass_to_mass;

  /*! Conversion factor from density in internal units to Hydrogen number
   * density in cgs */
  double rho_to_n_cgs;

  /*! Conversion factor from temperature to internal energy */
  float temp_to_u_factor;

  /*! Conversion from (table energy) x (stellar mass in internal units) to
   *  internal energy units. Table energies are in erg / Msun [GANE] */
  double feedback_energy_to_internal;

  /* ------------- Parameters for IMF --------------- */

  /*! Array to store calculated IMF */
  double *imf;

  /*! Arrays to store IMF mass bins */
  double *imf_mass_bin;

  /*! Arrays to store IMF mass bins (log10)*/
  double *imf_mass_bin_log10;

  /*! Minimal stellar mass considered by the IMF (in solar masses) */
  double imf_min_mass_msun;

  /*! Maximal stellar mass considered by the IMF (in solar masses) */
  double imf_max_mass_msun;

  /*! Log 10 of the minimal stellar mass considered by the IMF (in solar masses)
   */
  double log10_imf_min_mass_msun;

  /*! Log 10 of the maximal stellar mass considered by the IMF (in solar masses)
   */
  double log10_imf_max_mass_msun;

  /* ------------ SNe feedback properties ------------ */

  /*! SNII feedback model: random, isotropic or minimum distance */
  enum SNII_feedback_models feedback_model;

  /*! Minimal stellar mass considered for SNII feedback (in solar masses) */
  double SNII_min_mass_msun;

  /*! Maximal stellar mass considered for SNII feedback (in solar masses) */
  double SNII_max_mass_msun;

  /*! Log 10 of the minimal stellar mass considered for SNII feedback (in solar
   * masses) */
  double log10_SNII_min_mass_msun;

  /*! Log 10 of the maximal stellar mass considered for SNII feedback (in solar
   * masses) */
  double log10_SNII_max_mass_msun;

  /* ------------ Enrichment sampling properties ------------ */

  /*! Star age above which the enrichment will be downsampled (in internal
   * units) */
  double stellar_evolution_age_cut;

  /*! Number of time-steps in-between two enrichment events */
  int stellar_evolution_sampling_rate;
};

void feedback_props_init(struct feedback_props *fp,
                         const struct phys_const *phys_const,
                         const struct unit_system *us,
                         struct swift_params *params,
                         const struct hydro_props *hydro_props,
                         const struct cosmology *cosmo);

#endif /* SWIFT_GANE_FEEDBACK_PROPERTIES_THERMAL_H */
