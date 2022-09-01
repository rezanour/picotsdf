#include "picotsdf.h"

#include <malloc.h>
#include <math.h>
#include <memory.h>
#include <minmax.h>

inline void vector3_mul(const vector3* v0, const vector3* v1, vector3* out) {
  out->x = v0->x * v1->x;
  out->y = v0->y * v1->y;
  out->z = v0->z * v1->z;
}

inline void vector3_div(const vector3* v0, const vector3* v1, vector3* out) {
  out->x = v0->x / v1->x;
  out->y = v0->y / v1->y;
  out->z = v0->z / v1->z;
}

inline void vector3_divi(const vector3* v0, const vector3i* v1, vector3* out) {
  out->x = v0->x / v1->x;
  out->y = v0->y / v1->y;
  out->z = v0->z / v1->z;
}

tsdf* tsdf_create(const vector3* size, const vector3i* resolution) {
  const long total_cells = resolution->x * resolution->y * resolution->z;
  tsdf* result = (tsdf*)malloc(sizeof(tsdf));
  if (result) {
    result->cells = (tsdf_cell*)malloc(total_cells * sizeof(tsdf_cell));
    if (result->cells) {
      memset(result->cells, 0, total_cells * sizeof(tsdf_cell));
      result->resolution = *resolution;
      result->size = *size;
      vector3_divi(size, resolution, &result->cell_size);
      return result;
    }
    free(result);
  }
  return 0;
}

void tsdf_destroy(tsdf* tsdf) {
  if (tsdf) {
    if (tsdf->cells)
      free(tsdf->cells);
    free(tsdf);
  }
}

typedef struct dda_params_t {
  vector3 p;
  vector3 delta_dist;
  vector3 side_dist;
  vector3 cell_center;
  vector3 cell_size;
  vector3i cell;
  vector3i step;
  vector3i resolution;
  tsdf_cell* cur_cell;
  long pitch;
  long slice;
  int hit;
} dda_params;

int test_ray_volume(const vector3* ray_start, const vector3* ray_dir, const vector3i* resolution, vector3* out_point) {
  // quick rejects
  if (ray_start->x < 0 && ray_dir->x <= 0) return 0;
  if (ray_start->y < 0 && ray_dir->y <= 0) return 0;
  if (ray_start->z < 0 && ray_dir->z <= 0) return 0;
  if (ray_start->x >= (float)resolution->x && ray_dir->x >= 0) return 0;
  if (ray_start->y >= (float)resolution->y && ray_dir->y >= 0) return 0;
  if (ray_start->z >= (float)resolution->z && ray_dir->z >= 0) return 0;

  if (ray_start->x < 0) {
    float scale = -ray_start->x / ray_dir->x;
    out_point->x = ray_start->x + ray_dir->x * scale;
    out_point->y = ray_start->y + ray_dir->y * scale;
    out_point->z = ray_start->z + ray_dir->z * scale;
    if (out_point->y < 0 || out_point->y >= (float)resolution->y) return 0;
    if (out_point->z < 0 || out_point->z >= (float)resolution->z) return 0;
  }
  else {
    float scale = (ray_start->x - resolution->x) / ray_dir->x;
    out_point->x = ray_start->x + ray_dir->x * scale;
    out_point->y = ray_start->y + ray_dir->y * scale;
    out_point->z = ray_start->z + ray_dir->z * scale;
    if (out_point->y < 0 || out_point->y >= (float)resolution->y) return 0;
    if (out_point->z < 0 || out_point->z >= (float)resolution->z) return 0;
  }
  return 1;
}

int dda_init(tsdf* tsdf, const vector3* ray_start, const vector3* ray_dir, dda_params* out_dda) {
  // translate start into floating point cell index space (each cell 1 unit apart)
  vector3_div(ray_start, &tsdf->cell_size, &out_dda->p);

  out_dda->cell.x = (int)out_dda->p.x;
  out_dda->cell.y = (int)out_dda->p.y;
  out_dda->cell.z = (int)out_dda->p.z;

  // if the ray is starting outside of the volume, do a ray test to see if we even hit the volume.
  // If we do, return the new ray start position on the boundary of the volume that we should use.
  if (out_dda->cell.x < 0 || out_dda->cell.x >= tsdf->resolution.x ||
    out_dda->cell.y < 0 || out_dda->cell.y >= tsdf->resolution.y ||
    out_dda->cell.z < 0 || out_dda->cell.z >= tsdf->resolution.z) {
    if (!test_ray_volume(ray_start, ray_dir, &tsdf->resolution, &out_dda->p)) {
      return 0;
    }
    out_dda->cell.x = (int)out_dda->p.x;
    out_dda->cell.y = (int)out_dda->p.y;
    out_dda->cell.z = (int)out_dda->p.z;
  }

  out_dda->delta_dist.x = (float)fabs(1.f / ray_dir->x);
  out_dda->delta_dist.y = (float)fabs(1.f / ray_dir->y);
  out_dda->delta_dist.z = (float)fabs(1.f / ray_dir->z);

  // init
  if (ray_dir->x < 0.f) {
    out_dda->step.x = -1;
    out_dda->side_dist.x = (out_dda->p.x - out_dda->cell.x) * out_dda->delta_dist.x;
  }
  else {
    out_dda->step.x = 1;
    out_dda->side_dist.x = (out_dda->cell.x + 1.0f - out_dda->p.x) * out_dda->delta_dist.x;
  }
  if (ray_dir->y < 0.f) {
    out_dda->step.y = -1;
    out_dda->side_dist.y = (out_dda->p.y - out_dda->cell.y) * out_dda->delta_dist.y;
  }
  else {
    out_dda->step.y = 1;
    out_dda->side_dist.y = (out_dda->cell.y + 1.0f - out_dda->p.y) * out_dda->delta_dist.y;
  }
  if (ray_dir->z < 0.f) {
    out_dda->step.z = -1;
    out_dda->side_dist.z = (out_dda->p.z - out_dda->cell.z) * out_dda->delta_dist.z;
  }
  else {
    out_dda->step.z = 1;
    out_dda->side_dist.z = (out_dda->cell.z + 1.0f - out_dda->p.z) * out_dda->delta_dist.z;
  }

  out_dda->hit = 0;
  out_dda->resolution = tsdf->resolution;
  out_dda->pitch = out_dda->resolution.x;
  out_dda->slice = out_dda->resolution.x * out_dda->resolution.y;
  out_dda->cur_cell = &tsdf->cells[out_dda->cell.z * out_dda->slice + out_dda->cell.y * out_dda->pitch + out_dda->cell.x];
  out_dda->cell_center.x = (out_dda->cell.x + 0.5f) * tsdf->cell_size.x;
  out_dda->cell_center.y = (out_dda->cell.y + 0.5f) * tsdf->cell_size.y;
  out_dda->cell_center.z = (out_dda->cell.z + 0.5f) * tsdf->cell_size.z;
  out_dda->cell_size = tsdf->cell_size;

  return 1;
}

inline int dda_step(dda_params* dda) {
  if (dda->side_dist.x < dda->side_dist.y && dda->side_dist.x < dda->side_dist.z) {
    dda->side_dist.x += dda->delta_dist.x;
    dda->cell.x += dda->step.x;
    dda->cell_center.x += dda->step.x * dda->cell_size.x;
    dda->cur_cell += dda->step.x;
    if (dda->cell.x < 0 || dda->cell.x >= dda->resolution.x) {
      return 0;
    }
  }
  else if (dda->side_dist.y < dda->side_dist.x && dda->side_dist.y < dda->side_dist.z) {
    dda->side_dist.y += dda->delta_dist.y;
    dda->cell.y += dda->step.y;
    dda->cell_center.y += dda->step.y * dda->cell_size.y;
    dda->cur_cell += dda->step.y * dda->pitch;
    if (dda->cell.y < 0 || dda->cell.y >= dda->resolution.y) {
      return 0;
    }
  }
  else {
    dda->side_dist.z += dda->delta_dist.z;
    dda->cell.z += dda->step.z;
    dda->cell_center.z += dda->step.z * dda->cell_size.z;
    dda->cur_cell += dda->step.z * dda->slice;
    if (dda->cell.z < 0 || dda->cell.z >= dda->resolution.z) {
      return 0;
    }
  }
  return 1;
}

void tsdf_add_sample_z(tsdf* tsdf, const vector3* ray_start, const vector3* ray_dir, const vector3* forward, float sample_z) {
  dda_params dda;
  if (!dda_init(tsdf, ray_start, ray_dir, &dda)) {
    return; // doesn't intersect volume
  }

  float dist_x, dist_y, dist_z, dist, z;
  while (!dda.hit) {
    if (!dda_step(&dda)) {
      return; // exited the volume
    }
    dist_x = (dda.cell_center.x - ray_start->x);
    dist_y = (dda.cell_center.y - ray_start->y);
    dist_z = (dda.cell_center.z - ray_start->z);
    z = (dist_x * forward->x) + (dist_y * forward->y) + (dist_z * forward->z);
    dist = sample_z - z;
    dda.cur_cell->distance = (dda.cur_cell->distance * dda.cur_cell->weight + dist) / (dda.cur_cell->weight + 1);
    dda.cur_cell->distance = min(dda.cur_cell->distance, 1.0f);
    dda.cur_cell->distance = max(dda.cur_cell->distance, -1.0f);
    dda.cur_cell->weight += 1.0f;
  }
}

void tsdf_add_sample_euclidean(tsdf* tsdf, const vector3* ray_start, const vector3* ray_dir, float sample_dist) {
  dda_params dda;
  if (!dda_init(tsdf, ray_start, ray_dir, &dda)) {
    return; // doesn't intersect volume
  }

  float dist_x, dist_y, dist_z, dist;
  while (!dda.hit) {
    if (!dda_step(&dda)) {
      return; // exited the volume
    }
    dist_x = (dda.cell_center.x - ray_start->x);
    dist_y = (dda.cell_center.y - ray_start->y);
    dist_z = (dda.cell_center.z - ray_start->z);
    dist = sample_dist - sqrtf(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);
    dda.cur_cell->distance = (dda.cur_cell->distance * dda.cur_cell->weight + dist) / (dda.cur_cell->weight + 1);
    dda.cur_cell->distance = min(dda.cur_cell->distance, 1.0f);
    dda.cur_cell->distance = max(dda.cur_cell->distance, -1.0f);
    dda.cur_cell->weight += 1.0f;
  }
}

int tsdf_sample(tsdf* tsdf, const vector3* ray_start, const vector3* ray_dir, vector3* out_sample) {
  dda_params dda;
  if (!dda_init(tsdf, ray_start, ray_dir, &dda)) {
    return 0; // doesn't intersect volume
  }

  while (!dda.hit) {
    if (!dda_step(&dda)) {
      return 0; // exited the volume
    }
    if (dda.cur_cell->distance < 1.0f && dda.cur_cell->weight > 0.0f) {
      dda.hit = 1;
    }
  }

  return dda.hit;
}
