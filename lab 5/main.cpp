#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>

#define M_PI 3.14159265358979323846

#include <algorithm>

class Vector
{
public:
  explicit Vector(double x = 0, double y = 0, double z = 0)
  {
    data[0] = x;
    data[1] = y;
    data[2] = z;
  }
  double norm2() const
  {
    return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
  }
  double norm() const
  {
    return sqrt(norm2());
  }
  void normalize()
  {
    double n = norm();
    data[0] /= n;
    data[1] /= n;
    data[2] /= n;
  }
  double operator[](int i) const { return data[i]; };
  double &operator[](int i) { return data[i]; };
  double data[3];
};

Vector operator+(const Vector &a, const Vector &b)
{
  return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector &a, const Vector &b)
{
  return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector &b)
{
  return Vector(a * b[0], a * b[1], a * b[2]);
}
Vector operator*(const Vector &a, const double b)
{
  return Vector(a[0] * b, a[1] * b, a[2] * b);
}
Vector operator/(const Vector &a, const double b)
{
  return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector &a, const Vector &b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector &a, const Vector &b)
{
  return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
Vector gamma_correction(const Vector &a)
{
  return Vector(std::min<double>(255, pow(a[0], 1 / 2.2)), std::min<double>(255, pow(a[1], 1 / 2.2)), std::min<double>(255, pow(a[2], 1 / 2.2)));
}
Vector operator*(const Vector &a, const Vector &b)
{
  return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

Vector random_direction()
{
  double r1 = ((double)rand() / (RAND_MAX));
  double r2 = ((double)rand() / (RAND_MAX));
  double x = cos(2 * M_PI * r1) / sqrt(r2 * (1 - r2));
  double y = sin(2 * M_PI * r1) / sqrt(r2 * (1 - r2));
  double z = 1 - 2 * r2;
  return Vector(x, y, z);
}

int main()
{
  int W, H, C;
  int model_W, model_H, model_C;
  int nbiter = 2;

  // stbi_set_flip_vertically_on_load(true);
  unsigned char *image = stbi_load("input.jpg",
                                   &W,
                                   &H,
                                   &C,
                                   STBI_rgb);

  unsigned char *model = stbi_load("redim.jpg",
                                   &model_W,
                                   &model_H,
                                   &model_C,
                                   STBI_rgb);
  

  int n = W * H;
  std::vector<std::pair<int, int>> proj_I(n);
  std::vector<std::pair<int, int>> proj_M(n);
  for (int iter = 0; iter < nbiter; iter++)
  {
    Vector v = random_direction();
    for (int i = 0; i < n; i++)
    {
      unsigned char *current_image_pixel_index = image + C * i;
      unsigned char *current_model_pixel_index = model + model_C * i;
      Vector input_pixel(*current_image_pixel_index, *(current_image_pixel_index+1), *(current_image_pixel_index+2));
      Vector model_pixel(*current_model_pixel_index, *(current_model_pixel_index+1), *(current_model_pixel_index+2));
      int projI1 = dot(input_pixel, v);
      int projM1 = dot(model_pixel, v);
      proj_I[i] = std::pair<int, int>(projI1, i);
      proj_M[i] = std::pair<int, int>(projM1, i);
    }
    std::sort(proj_I.begin(), proj_I.end());
    std::sort(proj_M.begin(), proj_M.end());
    for (size_t i = 0; i < n; i++)
    {
      int change_index = proj_I[i].second;
      Vector change_vector = v * (proj_M[i].first - proj_I[i].first);
      unsigned char *current_index = image + change_index*C;
      *current_index = *current_index + change_vector[0];
      *(current_index+1) = *(current_index+1) + change_vector[1];
      *(current_index+2) = *(current_index+2) + change_vector[2];
    }
  }
  
  stbi_write_png("image.png", W, H, 3, &image[0], 0);
  return 0;
}

