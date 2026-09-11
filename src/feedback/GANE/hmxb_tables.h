/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_GANE_FEEDBACK_HMXB_TABLES_H
#define SWIFT_GANE_FEEDBACK_HMXB_TABLES_H

/* Local includes. */
#include "inline.h"

/*! Number of metallicity bins considered for the HMXB feedback */
#define gane_feedback_HMXB_N_metals ?

/*! Number of age bins considered for the HMXB feedback */
#define gane_feedback_HMXB_N_ages ?

/**
 * @brief reads HMXB energy tables and stores them in stars_props data
 * struct
 *
 * @param feedback_props the #feedback_props data struct to read the table into.
 */
INLINE static void read_HMXB_tables(struct feedback_props *feedback_props) {

#ifdef HAVE_HDF5

  /* filenames to read HDF5 files */
  char fname[256];

  hid_t file_id, dataset;
  herr_t status;

  /* Open HMXB tables for reading */
  sprintf(fname, "%s/HMXB.hdf5", feedback_props->HMXB_table_path);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) error("unable to open file %s\n", fname);

  /* read array of metallicities */
  dataset = H5Dopen(file_id, "Metallicities", H5P_DEFAULT);
  if (dataset < 0) error("Error opening HMXB Metallicities dataset");

  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   feedback_props->HMXB_energies.metallicity);
  if (status < 0) error("error reading HMXB metallicities");

  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* read array of ages */
  dataset = H5Dopen(file_id, "Ages", H5P_DEFAULT);
  if (dataset < 0) error("Error opening HMXB Ages dataset");

  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   feedback_props->HMXB_energies.age);
  if (status < 0) error("error reading HMXB ages");

  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* read HMXB energy tables */

  /* Temporary array with identical static dimensions, exactly like lifetimes */
  double temp_energy[gane_feedback_HMXB_N_ages]
                         [gane_feedback_HMXB_N_metals];

  dataset = H5Dopen(file_id, "Energy", H5P_DEFAULT);
  if (dataset < 0) error("error opening HMXB Energy dataset");

  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   temp_energy);
  if (status < 0) error("error reading HMXB Energy table");

  status = H5Dclose(dataset);
  if (status < 0) error("error closing Energy dataset");

  /* Copy from temporary array to final destination struct */
  for (int i = 0; i < gane_feedback_HMXB_N_ages; i++) {
    for (int j = 0; j < gane_feedback_HMXB_N_metals; j++) {
      feedback_props->HMXB_energies.energy[i][j] = temp_energy[i][j];
    }
  }

  /* Close file */
  status = H5Fclose(file_id);
  if (status < 0) error("error closing HMXB file %s", fname);

#endif
}

/**
 * @brief allocates space for the HMXB tables
 *
 * @param feedback_props the #feedback_props data struct to store the tables in
 */
INLINE static void allocate_HMXB_tables(
    struct feedback_props *feedback_props) {

  /* Allocate array for HMXB metallicity bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->HMXB_energies.metallicity,
                     SWIFT_STRUCT_ALIGNMENT,
                     gane_feedback_HMXB_N_metals * sizeof(double)) != 0) {
    error("Failed to allocate HMXB metallicity array");
  }

  /* Allocate array for HMXB age bins */
  if (swift_memalign("feedback-tables",
                     (void **)&feedback_props->HMXB_energies.age,
                     SWIFT_STRUCT_ALIGNMENT,
                     gane_feedback_HMXB_N_ages * sizeof(double)) != 0) {
    error("Failed to allocate HMXB age array");
  }

  /* Allocate 2D energy table: rows (ages) and columns (metals) */
  feedback_props->HMXB_energies.energy =
      (double **)malloc(gane_feedback_HMXB_N_ages * sizeof(double *));
  if (feedback_props->HMXB_energies.energy == NULL) {
    error("Failed to allocate HMXB energy row pointers");
  }

  for (int i = 0; i < gane_feedback_HMXB_N_ages; i++) {
    if (swift_memalign("feedback-tables",
                       (void **)&feedback_props->HMXB_energies.energy[i],
                       SWIFT_STRUCT_ALIGNMENT,
                       gane_feedback_HMXB_N_metals * sizeof(double)) != 0) {
      error("Failed to allocate HMXB energy column array");
    }
  }
}

/**
 * @brief Interpolating function for HMXB tables. Returns interpolated energy
 * given the #spart's age and metallicity. (REVISAR)
 * 
 * @param star_age age of the #spart.
 * @param Z_birth #spart's metallicity at birth.
 * @param feedback_props the properties of the feedback model.
 */
INLINE static double interpolate_HMXB_energy(
  const double star_age, const double Z_birth,
  const struct feedback_props *feedback_props) {
  
  /* Get table */
  const struct HMXB_table *HMXB_table = &feedback_props->HMXB_energies;
  
  /* Find closest cells in HMXB table to (star_age, Z_birth) */
  
  /* Find age cells */
  int i;
  for (i=0; i < gane_feedback_HMXB_N_ages - 1 && HMXB_table->age[i] <= star_age; i++) {
    continue;
  }
  const int age_index = i - 1;
  const double age_1 = HMXB_table->age[age_index];
  const double age_2 = HMXB_table->age[age_index + 1];

  /* Find metallicity cells */
  int j;
  for (j=0; j < gane_feedback_HMXB_N_metals - 1 && HMXB_table->metallicity[j] <= Z_birth; j++) {
    continue;
  }
  const int Z_index = j - 1;
  const double Z_1 = HMXB_table->metallicity[Z_index];
  const double Z_2 = HMXB_table->metallicity[Z_index + 1];
  
  /* Normalize cell's age and metallicity locations */
  const float d_age = (float)((star_age - age_1) / (age_2 - age_1));
  const float d_Z = (float)((Z_birth - Z_1) / (Z_2 - Z_1));

  /* Interpolate to obtain energy */
  const double E = interpolate_2d(HMXB_table->energy, age_index, Z_index, d_age, d_Z);

  return E;
  }

#endif /* SWIFT_GANE_FEEDBACK_HMXB_TABLES_H */
