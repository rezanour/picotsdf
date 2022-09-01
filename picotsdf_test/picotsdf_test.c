#include <stdio.h>
#include <time.h>

#include <picotsdf.h>

int main(int argc, char* argv[]) {
  vector3 size = { 1.f, 1.f, 0.5f };
  vector3i resolution = { 10, 10, 5 };

  tsdf* tsdf = tsdf_create(&size, &resolution);

  vector3 ray_start = { 0.21f, 0.18f, 0.f };
  vector3 ray_dir = { 1.f, 0.f, 0.f };

  long i;
  clock_t start, stop;
  start = clock();

  for (i = 0; i < 640 * 480; ++i) {
    ray_dir.y += 0.001f;
    ray_dir.z += 0.001f;
    tsdf_add_sample_euclidean(tsdf, &ray_start, &ray_dir, 0.4f);
  }

  stop = clock();
  printf("elapsed time: %ums", (stop - start) * 1000 / CLOCKS_PER_SEC);

  //vector3 sample;
  //ray_start.y = 3;
  //if (tsdf_sample(tsdf, &ray_start, &ray_dir, &sample)) {
  //  printf("hit the sample.\n");
  //}

  tsdf_destroy(tsdf);

  return 0;
}
