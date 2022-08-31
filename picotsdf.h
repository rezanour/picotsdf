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
  // size of a single cell (inverse)
  vector3 inv_cell_size;
  // the cell data for the volume
  tsdf_cell* cells;

  long pitch;
  long slice;
} tsdf;

tsdf* tsdf_create(const vector3* size, const vector3i* resolution);
void tsdf_destroy(tsdf* tsdf);

// sample_z is camera z value (ie. depth buffer), not euclidean distance
void tsdf_addsample(tsdf* tsdf, const vector3* ray_start, const vector3* ray_dir, float sample_z);

int tsdf_sample(tsdf* tsdf, const vector3* ray_start, const vector3* ray_dir, vector3* out_sample);
