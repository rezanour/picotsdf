#pragma once

//
// ultra minimal and compact tsdf implementation
// to simplify math, assumes that the axis-aligned volume starts at origin and extends to size
// clients must transform inputs & outputs accordingly!
//

typedef struct vector3_t {
  float x, y, z;
} vector3;

typedef struct vector3i_t {
  long x, y, z;
} vector3i;

typedef struct tsdf_cell_t {
  float distance;
  float weight;
} tsdf_cell;

typedef struct tsdf_t {
  // number of cells per axis
  vector3i resolution;
  // total size of the volume
  vector3 size;
  // size of a single cell
  vector3 cell_size;
  // the cell data for the volume
  tsdf_cell* cells;
} tsdf;

tsdf* tsdf_create(const vector3* size, const vector3i* resolution);
void tsdf_destroy(tsdf* tsdf);

// sample_z is the camera space z coordinate (ie. depth buffer value) of the surface sample
void tsdf_add_sample_z(tsdf* tsdf, const vector3* ray_start, const vector3* ray_dir, const vector3* forward, float sample_z);

// sample_dist is the euclidean distance from the camera position to the surface sample
void tsdf_add_sample_euclidean(tsdf* tsdf, const vector3* ray_start, const vector3* ray_dir, float sample_dist);

// NOTE not implemented yet!
int tsdf_sample(tsdf* tsdf, const vector3* ray_start, const vector3* ray_dir, vector3* out_sample);
