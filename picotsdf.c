#include "picotsdf.h"

#include <malloc.h>
#include <math.h>
#include <memory.h>

static void vector3_mul(const vector3* v0, const vector3* v1, vector3* out) {
  out->x = v0->x * v1->x;
  out->y = v0->y * v1->y;
  out->z = v0->z * v1->z;
}

static void vector3_div(const vector3* v0, const vector3* v1, vector3* out) {
  out->x = v0->x / v1->x;
  out->y = v0->y / v1->y;
  out->z = v0->z / v1->z;
}

static void vector3i_div(const vector3i* v0, const vector3* v1, vector3* out) {
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
      vector3i_div(resolution, size, &result->inv_cell_size);
      result->pitch = resolution->x;
      result->slice = resolution->x * resolution->y;
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

#define DDA_SETUP() \
  /* TODO: does ray ever intersect volume ? early out if not (if it's worth testing) */ \
\
/* translate start into floating point cell index */ \
vector3 p; \
vector3_mul(ray_start, &tsdf->inv_cell_size, &p); \
\
vector3i cell = { (int)p.x, (int)p.y, (int)p.z }; \
\
/* TODO: check if out of bounds */ \
\
/* TODO: need to do abs() here, but have to emulate it as SSE does't have */ \
vector3 delta_dist = { (float)fabs(1.f / ray_dir->x), (float)fabs(1.f / ray_dir->y), (float)fabs(1.f / ray_dir->z) }; \
\
vector3 side_dist = { 0.f }; \
vector3i step = { 0 }; \
\
/* init */ \
if (ray_dir->x < 0.f) { \
  step.x = -1; \
  side_dist.x = (p.x - cell.x) * delta_dist.x; \
} \
else { \
  step.x = 1; \
  side_dist.x = (cell.x + 1.0f - p.x) * delta_dist.x; \
} \
if (ray_dir->y < 0.f) { \
  step.y = -1; \
  side_dist.y = (p.y - cell.y) * delta_dist.y; \
} \
else { \
  step.y = 1; \
  side_dist.y = (cell.y + 1.0f - p.y) * delta_dist.y; \
} \
if (ray_dir->z < 0.f) { \
  step.z = -1; \
  side_dist.z = (p.z - cell.z) * delta_dist.z; \
} \
else { \
  step.z = 1; \
  side_dist.z = (cell.z + 1.0f - p.z) * delta_dist.z; \
} \


#define DDA_BEGIN() \
  /* DDA */ \
int hit = 0; \
tsdf_cell* c = &tsdf->cells[cell.z * tsdf->slice + cell.y * tsdf->pitch + cell.x]; \
while (!hit) { \
  if (side_dist.x < side_dist.y && side_dist.x < side_dist.z) { \
    side_dist.x += delta_dist.x; \
    cell.x += step.x; \
    c += step.x; \
    if (cell.x < 0 || cell.x >= tsdf->resolution.x) { \
      break; \
    } \
  } \
  else if (side_dist.y < side_dist.x && side_dist.y < side_dist.z) { \
    side_dist.y += delta_dist.y; \
    cell.y += step.y; \
    c += step.y * tsdf->pitch; \
    if (cell.y < 0 || cell.y >= tsdf->resolution.y) { \
      break; \
    } \
  } \
  else { \
    side_dist.z += delta_dist.z; \
    cell.z += step.z; \
    c += step.z * tsdf->slice; \
    if (cell.z < 0 || cell.z >= tsdf->resolution.z) { \
      break; \
    } \
  } \

#define DDA_END() }



// sample_z is camera z value (ie. depth buffer), not euclidean distance
void tsdf_addsample(tsdf* tsdf, const vector3* ray_start, const vector3* ray_dir, float sample_z) {
  DDA_SETUP();

  DDA_BEGIN();
    c->distance = 0.0f;
    c->weight = 1.0f;
  DDA_END();
}

int tsdf_sample(tsdf* tsdf, const vector3* ray_start, const vector3* ray_dir, vector3* out_sample) {
  DDA_SETUP();

  DDA_BEGIN();
    if (c->distance < 1.f && c->weight > 0.f) {
      //*out_sample = { 1.f, 1.f, 1.f };
      hit = 1;
    }
  DDA_END();

  return hit;
}
