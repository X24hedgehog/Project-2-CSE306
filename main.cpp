#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#include <math.h>
#include <cmath>
#include <ctime>
#include <iostream>
#define M_PI 3.14159265358979323846

#include <string>
#include <stdio.h>
#include <algorithm>
#include <random>
#include <omp.h>
#include "lbfgs.h"
// #include <stb_image_write.h>
#include <sstream>

std::random_device my_random_device;
std::default_random_engine my_random_engine(my_random_device());
std::uniform_real_distribution<float> my_distribution(0, 1);
std::uniform_real_distribution<float> my_weight_distribution(0.8, 1);

double generate()
{
  return my_distribution(my_random_engine);
}

double generate_weight()
{
  return my_weight_distribution(my_random_engine);
}

class Vector
{
public:
  explicit Vector(double x = 0, double y = 0)
  {
    data[0] = x;
    data[1] = y;
  }
  double norm2() const
  {
    return data[0] * data[0] + data[1] * data[1];
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
  }
  void print_position()
  {
    std::cout << "x: " << data[0] << "  " << "y: " << data[1] << std::endl;
  }
  double operator[](int i) const { return data[i]; };
  double &operator[](int i) { return data[i]; };
  double data[3];
};

Vector operator+(const Vector &a, const Vector &b)
{
  return Vector(a[0] + b[0], a[1] + b[1]);
}
Vector operator-(const Vector &a, const Vector &b)
{
  return Vector(a[0] - b[0], a[1] - b[1]);
}
Vector operator*(const double a, const Vector &b)
{
  return Vector(a * b[0], a * b[1]);
}
Vector operator*(const Vector &a, const double b)
{
  return Vector(a[0] * b, a[1] * b);
}
Vector operator/(const Vector &a, const double b)
{
  return Vector(a[0] / b, a[1] / b);
}
double dot(const Vector &a, const Vector &b)
{
  return a[0] * b[0] + a[1] * b[1];
}
double cross(const Vector &a, const Vector &b)
{
  return (a[0] * b[1] - a[1] * b[0]);
}
Vector operator*(const Vector &a, const Vector &b)
{
  return Vector(a[0] * b[0], a[1] * b[1]);
}
double norm(const Vector &a)
{
  return sqrt(pow(a[0], 2) + pow(a[1], 2));
}
Vector normalization(const Vector &a)
{
  return a / norm(a);
}

Vector check_intersect(Vector A, Vector B, std::pair<Vector, Vector> l)
{
  Vector u = l.first;
  Vector v = l.second;
  Vector N = Vector(v[1] - u[1], u[0] - v[0]);
  double t = dot(u - A, N) / dot(B - A, N);
  if (t < 0 || t > 1)
  {
    return Vector(INFINITY, INFINITY);
    // No intersection is equivalent to intersection at infinity
  }
  return A + t * (B - A);
};

Vector check_intersect_voronoi(Vector A, Vector B, std::pair<Vector, Vector> l, double w1, double w2)
{
  Vector u = l.first;
  Vector v = l.second;
  Vector M = (u + v) / 2 + (w1 - w2) / (2 * pow(norm(v - u), 2)) * (v - u);
  Vector N = u - v;
  double t = dot(M - A, N) / dot(B - A, N);
  if (t < 0 || t > 1)
  {
    return Vector(INFINITY, INFINITY);
    // No intersection is equivalent to intersection at infinity
  }
  return A + t * (B - A);
};

bool check_inside(Vector P, std::pair<Vector, Vector> edge)
{
  Vector u = edge.first;
  Vector v = edge.second;
  Vector N = Vector(v[1] - u[1], u[0] - v[0]);
  if (dot(P - u, N) <= 0)
    return true;
  return false;
}

bool check_inside_voronoi(Vector X, std::pair<Vector, Vector> edge, double w1, double w2)
{
  Vector u = edge.first;
  Vector v = edge.second;
  Vector M = (u + v) / 2 + (w1 - w2) / (2 * pow(norm(v - u), 2)) * (v - u);
  if (dot(X - M, v - u) < 0)
  {
    return true;
  }
  return false;
}

// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon
{
public:
  std::vector<Vector> vertices;
  Polygon(){};
  Polygon(std::vector<Vector> &my_vertices)
  {
    vertices = my_vertices;
  }

  Vector get_vertice(size_t index) const
  {
    return vertices[index];
  }

  double get_area() const
  {
    double area = 0;
    if (vertices.size() > 0)
    {
      Vector o = vertices[0];
      for (int i = 1; i + 1 < vertices.size(); i++)
      {
        Vector a = vertices[i] - o;
        Vector b = vertices[i + 1] - o;
        area += abs(cross(a, b));
      }
    }

    return area / 2;
  }

  void print_vertices()
  {
    std::cout << "Start printing vertices of a polygon" << std::endl;
    for (int i = 0; i < vertices.size(); i++)
    {
      std::cout << "Coord x: " << vertices[i][0] << "    " << "Coord y: " << vertices[i][1] << std::endl;
    }
    std::cout << " " << std::endl;
  }

  void clip_edge_voronoi(std::pair<Vector, Vector> edge, double w1, double w2)
  {
    Polygon outPolygon = Polygon();
    for (int i = 0; i < vertices.size(); i++)
    {
      Vector curVertex = vertices[i];
      Vector prevVertex = vertices[(i > 0) ? (i - 1) : vertices.size() - 1];

      Vector intersection = check_intersect_voronoi(prevVertex, curVertex, edge, w1, w2);
      if (check_inside_voronoi(curVertex, edge, w1, w2))
      {
        if (!check_inside_voronoi(prevVertex, edge, w1, w2))
        {
          outPolygon.vertices.push_back(intersection);
        }
        outPolygon.vertices.push_back(curVertex);
      }
      else if (check_inside_voronoi(prevVertex, edge, w1, w2))
      {
        outPolygon.vertices.push_back(intersection);
      }
    }
    this->vertices = outPolygon.vertices;
  }

  void clip_polygon(Polygon c)
  {
    Polygon outPolygon;
    for (int i = 0; i < c.vertices.size(); i++)
    {
      int prev_index;
      if (i != 0)
      {
        prev_index = i - 1;
      }
      else
      {
        prev_index = c.vertices.size() - 1;
      }
      std::pair<Vector, Vector> clipEdge(c.vertices[prev_index], c.vertices[i]);
      outPolygon = Polygon();
      // std::cout << "Clip polygon edge:" << std::endl;
      // std::cout << "Point 1: " << c.vertices[prev_index][0] << " and " << c.vertices[prev_index][1] << std::endl;
      // std::cout << "Point 2: " << c.vertices[i][0] << " and " << c.vertices[i][1] << std::endl;
      // c.vertices[prev_index].print_position();
      // c.vertices[i].print_position();
      // std::cout << " " << std::endl;
      // std::cout << "New subject polygon: " << std::endl;
      // this->print_vertices();
      for (int j = 0; j < this->vertices.size(); j++)
      {
        Vector curVertex = vertices[j];
        Vector prevVertex = vertices[(j > 0) ? (j - 1) : vertices.size() - 1];
        // std::cout << "Subject polygon edge:" << std::endl;
        // prevVertex.print_position();
        // curVertex.print_position();
        Vector intersection = check_intersect(prevVertex, curVertex, clipEdge);
        // std::cout << "Intersect: " << std::endl;
        // intersection.print_position();
        if (check_inside(curVertex, clipEdge))
        {
          // std::cout << "Current vertext " << curVertex[0] << " " << curVertex[1] << " inside clip edge" << std::endl;
          if (!check_inside(prevVertex, clipEdge))
          {
            // std::cout << "Previous vertext " << prevVertex[0] << " " << prevVertex[1] << " not inside clip edge" << std::endl;
            outPolygon.vertices.push_back(intersection);
            // std::cout << "outPolygon add intersection " << intersection[0] << " " << intersection[1] << std::endl;
          }

          outPolygon.vertices.push_back(curVertex);
          // std::cout << "outPolygon add current vertext " << curVertex[0] << " " << curVertex[1] << std::endl;
        }
        else if (check_inside(prevVertex, clipEdge))
        {
          // std::cout << "Current vertext " << curVertex[0] << " " << curVertex[1] << " not inside clip edge" << std::endl;
          // std::cout << "Previous vertext " << prevVertex[0] << " " << prevVertex[1] << " inside clip edge" << std::endl;
          outPolygon.vertices.push_back(intersection);
          // std::cout << "outPolygon add intersection " << intersection[0] << " " << intersection[1] << std::endl;
        }
      }
      this->vertices = outPolygon.vertices;
    }
    this->vertices = outPolygon.vertices;
  }

  Vector compute_centroid()
  {
    Vector centroid;
    double area = this->get_area();
    double pos_x = 0;
    double pos_y = 0;
    for (int i = 0; i < vertices.size() - 1; i++)
    {
      double xi = vertices[i][0];
      double xip1 = vertices[i + 1][0];
      double yi = vertices[i][1];
      double yip1 = vertices[i + 1][1];
      pos_x += (xi + xip1) * (xi * yip1 - xip1 * yi);
      pos_y += (yi + yip1) * (xi * yip1 - xip1 * yi);
    }
    pos_x /= (6 * area);
    pos_y /= (6 * area);
    centroid[0] = pos_x;
    centroid[1] = pos_y;
    return centroid;
  }
};

int sgn(double num)
{
  if (num < 0)
  {
    return -1;
  }
  return 1;
}

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, const std::vector<Vector> &vertices, std::string filename, std::string fillcol = "none")
{
  FILE *f = fopen(filename.c_str(), "w+");
  fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
  for (int i = 0; i < polygons.size(); i++)
  {
    fprintf(f, "<g>\n");
    fprintf(f, "<polygon points = \"");
    for (int j = 0; j < polygons[i].vertices.size(); j++)
    {
      fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
    }
    fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
    fprintf(f, "</g>\n");
  }
  for (int i = 0; i < vertices.size(); i++)
  {
    fprintf(f, "<g>\n");
    fprintf(f, "<circle cx= \"");
    fprintf(f, "%3.3f\" cy=\"%3.3f", vertices[i][0] * 1000, 1000 - vertices[i][1] * 1000);
    fprintf(f, "\" r=\"3\"\nfill = \"red\" stroke = \"red\"/>\n", fillcol.c_str());
    fprintf(f, "</g>\n");
  }
  fprintf(f, "</svg>\n");
  fclose(f);
}

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
// void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none")
// {
//   FILE *f = fopen(filename.c_str(), "w+");
//   fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
//   for (int i = 0; i < polygons.size(); i++)
//   {
//     fprintf(f, "<g>\n");
//     fprintf(f, "<polygon points = \"");
//     for (int j = 0; j < polygons[i].vertices.size(); j++)
//     {
//       fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
//     }
//     fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
//     fprintf(f, "</g>\n");
//   }
//   fprintf(f, "</svg>\n");
//   fclose(f);
// }

// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes)
{
  FILE *f;
  if (frameid == 0)
  {
    f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    fprintf(f, "<g>\n");
  }
  else
  {
    f = fopen(filename.c_str(), "a+");
  }
  fprintf(f, "<g>\n");
  for (int i = 0; i < polygons.size(); i++)
  {
    fprintf(f, "<polygon points = \"");
    for (int j = 0; j < polygons[i].vertices.size(); j++)
    {
      fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
    }
    fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
  }
  fprintf(f, "<animate\n");
  fprintf(f, "    id = \"frame%u\"\n", frameid);
  fprintf(f, "    attributeName = \"display\"\n");
  fprintf(f, "    values = \"");
  for (int j = 0; j < nbframes; j++)
  {
    if (frameid == j)
    {
      fprintf(f, "inline");
    }
    else
    {
      fprintf(f, "none");
    }
    fprintf(f, ";");
  }
  fprintf(f, "none\"\n    keyTimes = \"");
  for (int j = 0; j < nbframes; j++)
  {
    fprintf(f, "%2.3f", j / (double)(nbframes));
    fprintf(f, ";");
  }
  fprintf(f, "1\"\n   dur = \"5s\"\n");
  fprintf(f, "    begin = \"0s\"\n");
  fprintf(f, "    repeatCount = \"indefinite\"/>\n");
  fprintf(f, "</g>\n");
  if (frameid == nbframes - 1)
  {
    fprintf(f, "</g>\n");
    fprintf(f, "</svg>\n");
  }
  fclose(f);
}

std::vector<Polygon> compute_voronoi(std::vector<Vector> points, const double *weights)
{
  std::vector<Polygon> polygons_list;
  std::vector<Vector> square_vertices = {Vector(1, 0), Vector(1, 1), Vector(0, 1), Vector(0, 0)};
  Polygon subject_polygon = Polygon(square_vertices);

  for (int i = 0; i < points.size(); i++)
  {
    std::vector<Vector> my_vertices = {Vector(1, 0), Vector(1, 1), Vector(0, 1), Vector(0, 0)};
    Polygon subject_polygon = Polygon(my_vertices);
    for (int j = 0; j < points.size(); j++)
    {
      if (j != i)
      {
        Vector Pi = points[i];
        Vector Pj = points[j];
        std::pair<Vector, Vector> current_edge;
        current_edge.first = Pi;
        current_edge.second = Pj;
        subject_polygon.clip_edge_voronoi(current_edge, weights[i], weights[j]);
      }
    }
    polygons_list.push_back(subject_polygon);
  }
  return polygons_list;
}

struct data
{
  std::vector<Vector> *points;
  double *lambdas;
};

static lbfgsfloatval_t
evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step)
{
  data *params = (data *)instance;
  std::vector<Vector> p = *(params->points);
  double *l = params->lambdas;
  std::vector<double> s;
  std::vector<Polygon> polygons = compute_voronoi(p, x);

  for (int i = 0; i < polygons.size(); i++)
  {
    s.push_back(polygons[i].get_area());
  }

  lbfgsfloatval_t fx = 0.0;
  double integral_part;
  double polygon_integral;
  for (int i = 0; i < n; i++)
  {
    integral_part = 0;
    std::vector<Vector> current_vertices = polygons[i].vertices;
    int num_vertices = current_vertices.size();
    // std::cout << "Number of vertices: " << num_vertices << std::endl;
    // std::cout << " " << std::endl;
    // std::cout << " " << std::endl;
    // for (int k = 1; k < num_vertices + 1; k++)
    // {
    //   double xkm1 = current_vertices[k - 1][0];
    //   double ykm1 = current_vertices[k - 1][1];
    //   double xk;
    //   double yk;
    //   if (k == num_vertices)
    //   {
    //     xk = current_vertices[0][0];
    //     yk = current_vertices[0][1];
    //   }
    //   else
    //   {
    //     xk = current_vertices[k][0];
    //     yk = current_vertices[k][1];
    //   }
    //   integral_part += (xkm1 * yk - xk * ykm1) * (pow(xkm1, 2) + xkm1 * xk + pow(xk, 2) +
    //                                               pow(ykm1, 2) + ykm1 * yk + pow(yk, 2) - 4 * (p[i][0] * (xkm1 + xk) + p[i][1] * (ykm1 + yk)) +
    //                                               6 * pow(norm(p[i]), 2));
    // }
    // integral_part /= 12;
    if (num_vertices > 2)
    {
      Vector o = current_vertices[0];
      polygon_integral = 0;
      for (int j = 1; j + 1 < num_vertices; j++)
      {
        Vector a = current_vertices[j];
        Vector b = current_vertices[j + 1];
        std::vector<Vector> triangle_vertices;
        triangle_vertices.push_back(o);
        triangle_vertices.push_back(a);
        triangle_vertices.push_back(b);
        double current_sum = 0;
        for (int v1 = 0; v1 < 3; v1++)
        {
          for (int v2 = 0; v2 < v1 + 1; v2++)
          {
            current_sum += dot(triangle_vertices[v1] - p[i], triangle_vertices[v2] - p[i]);
          }
        }
        double current_area = abs(cross(b - o, a - o)) / 2;
        current_sum *= current_area / 6;
        polygon_integral += current_sum;
      }
    }
    else
    {
      polygon_integral = 0;
    }

    fx += polygon_integral - x[i] * (s[i] - l[i]);
    g[i] = s[i] - l[i];
  }
  return -1 * fx;
}

Polygon generate_disk(Vector center, double radius, int num_edges = 20)
{
  std::vector<Vector> vertices;
  for (int i = 0; i < num_edges; i++)
  {
    Vector current_vertext;
    current_vertext[0] = center[0] + radius * cos(2 * i * M_PI / num_edges);
    current_vertext[1] = center[1] + radius * sin(2 * i * M_PI / num_edges);
    vertices.push_back(current_vertext);
  }
  Polygon out_polygon = Polygon(vertices);

  // out_polygon.print_vertices();
  return out_polygon;
}

std::vector<Polygon> intersect_disks(std::vector<Vector> points, double *weights, int num_edges = 100)
{
  std::vector<Polygon> polygons = compute_voronoi(points, weights);
  std::vector<double> radii;
  for (int r = 0; r < points.size(); r++)
  {
    radii.push_back(sqrt(weights[r] - weights[points.size()]));
  }
  for (int i = 0; i < polygons.size(); i++)
  {
    Polygon clip_polygon = generate_disk(points[i], radii[i]);
    polygons[i].clip_polygon(clip_polygon);
  }
  return polygons;
}

struct fluid_data
{
  std::vector<Vector> *points;
  double dfv;
};

static lbfgsfloatval_t fluid_evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step)
{
  fluid_data *params = (fluid_data *)instance;
  std::vector<Vector> p = *(params->points);
  double desired_fluid_volumn = params->dfv;
  std::vector<double> s;
  std::vector<Polygon> polygons = compute_voronoi(p, x);

  for (int i = 0; i < polygons.size(); i++)
  {
    Polygon clip_polygon = generate_disk(p[i], sqrt(x[i] - x[n - 1]));
    polygons[i].clip_polygon(clip_polygon);
    s.push_back(polygons[i].get_area());
  }

  lbfgsfloatval_t fx = 0.0;
  double integral_part;
  double polygon_integral;
  for (int i = 0; i < n - 1; i++)
  {
    integral_part = 0;
    std::vector<Vector> current_vertices = polygons[i].vertices;
    int num_vertices = current_vertices.size();
    if (num_vertices > 2)
    {
      Vector o = current_vertices[0];
      polygon_integral = 0;
      for (int j = 1; j + 1 < num_vertices; j++)
      {
        Vector a = current_vertices[j];
        Vector b = current_vertices[j + 1];
        std::vector<Vector> triangle_vertices;
        triangle_vertices.push_back(o);
        triangle_vertices.push_back(a);
        triangle_vertices.push_back(b);
        double current_sum = 0;
        for (int v1 = 0; v1 < 3; v1++)
        {
          for (int v2 = 0; v2 < v1 + 1; v2++)
          {
            current_sum += dot(triangle_vertices[v1] - p[i], triangle_vertices[v2] - p[i]);
          }
        }
        double current_area = abs(cross(b - o, a - o)) / 2;
        current_sum *= current_area / 6;
        polygon_integral += current_sum;
      }
    }
    else
    {
      polygon_integral = 0;
    }
    fx += polygon_integral + desired_fluid_volumn / (n - 1) * x[i] - s[i] * x[i];
    g[i] = s[i] - desired_fluid_volumn / (n - 1);
  }
  double total_area = 0;
  for (int polygon = 0; polygon < n - 1; polygon++)
  {
    total_area += polygons[polygon].get_area();
  }
  double estimated_air_volumn = 1 - total_area;
  fx += x[n - 1] * (1 - desired_fluid_volumn - estimated_air_volumn);
  g[n - 1] = estimated_air_volumn - (1 - desired_fluid_volumn);
  return -1 * fx;
}

static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls)
{
  printf("Iteration %d:\n", k);
  printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
  printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
  printf("\n");
  if (k > 200)
  {
    return 1;
  }
  return 0;
}

std::vector<Polygon> gallouet_merigot_scheme(std::vector<Vector> &positions,
                                             std::vector<Vector> &velocities,
                                             const std::vector<double> &masses,
                                             double desired_fluid_volume)
{
  double eps = 0.004;
  double dt = 0.002;
  Vector g;
  g[0] = 0.0;
  g[1] = -9.8;

  fluid_data *params = new fluid_data();
  params->points = &positions;
  params->dfv = desired_fluid_volume;
  const int n_points = positions.size();

  double *weights = new double[n_points + 1];
  for (int i = 0; i < positions.size(); i++)
  {
    double w = generate_weight();
    weights[i] = w;
  }
  weights[n_points] = 0.6;
  int ret = lbfgs(positions.size() + 1, weights, NULL, fluid_evaluate, NULL, (void *)params, nullptr);
  std::cout << ret << std::endl;

  std::vector<Polygon> polygons = intersect_disks(positions, weights);
  for (int i = 0; i < positions.size(); i++)
  {
    Vector f_spring = (polygons[i].compute_centroid() - positions[i]) / (pow(eps, 2));
    Vector f = f_spring + masses[i] * g;
    velocities[i] = velocities[i] + dt * f / masses[i];
    positions[i] = positions[i] + dt * velocities[i];
    for (int axis = 0; axis < 2; axis++)
    {
      if (positions[i][axis] > 1)
      {
        positions[i][axis] = 1;
        velocities[i][axis] = 0;
      }
      else if (positions[i][axis] < 0)
      {
        positions[i][axis] = 0;
        velocities[i][axis] = 0;
      }
    }
  }
  std::vector<Polygon> new_polygon = intersect_disks(positions, weights);
  return new_polygon;
}

void generate_voronoi() {
  const int n_points = 10;
  std::vector<Vector> points;
  for (int i = 0; i < n_points; i++)
  {
    double pos1 = generate();
    double pos2 = generate();
    points.push_back(Vector(pos1, pos2));
  }

  double weights[n_points];
  for (int i = 0; i < n_points; i++)
  {
    double w = generate();
    weights[i] = w / 4;
  }

  double balance_weights[n_points];
  for (int i = 0; i < n_points; i++)
  {
    double w = 1.0;
    balance_weights[i] = w;
  }

  std::vector<Polygon> polygons = compute_voronoi(points, balance_weights);
  std::vector<Polygon> weighted_polygons = compute_voronoi(points, weights);

  save_svg(polygons, points, "test_balance_image.svg");
  save_svg(weighted_polygons, points, "test_weighted_image.svg");
}

void generate_optimized_voronoi() {
  const int n_points = 200;
  std::vector<Vector> points;
  for (int i = 0; i < n_points; i++)
  {
    double pos1 = generate();
    double pos2 = generate();
    points.push_back(Vector(pos1, pos2));
  }

  double weights[n_points];
  for (int i = 0; i < n_points; i++)
  {
    double w = generate();
    weights[i] = w / 4;
  }

  double lambda_arr[n_points];
  const Vector C = Vector(0.5, 0.5);
  double sum_lambdas = 0;
  for (int i = 0; i < n_points; i++)
  {
    lambda_arr[i] = exp(-pow(norm(points[i] - C), 2) / 0.02);
    sum_lambdas += lambda_arr[i];
  }
  for (int i = 0; i < n_points; i++)
  {
    lambda_arr[i] /= sum_lambdas;
  }

  data *params = new data();
  params->points = &points;
  params->lambdas = lambda_arr;
  int ret = lbfgs(n_points, weights, NULL, evaluate, NULL, (void *)params, nullptr);
  std::cout << ret << std::endl;
  std::vector<Polygon> lbfsg_weighted_polygons = compute_voronoi(points, weights);


  save_svg(lbfsg_weighted_polygons, points, "test_lbfsg_weighted_image.svg");
}

void fluid_simulation()
{
  const int n_points = 20;
  int nb_frame = 150;
  std::vector<Vector> points;
  for (int i = 0; i < n_points; i++)
  {
    double pos1 = generate();
    double pos2 = generate();
    points.push_back(Vector(pos1, pos2));
  }

  std::vector<Vector> vels(n_points);
  std::vector<double> masses(n_points, 200);

  std::vector<Polygon> animated_polygons;

  for (int i = 0; i < nb_frame; i++)
  {
    animated_polygons = gallouet_merigot_scheme(points, vels, masses, 0.5);
    save_svg_animated(animated_polygons, "test_fluid_simulation.svg", i, nb_frame);
  }
}

void test_clip_disk(int num_edge)
{
  std::vector<Vector> square_vertices = {Vector(1, 0), Vector(1, 1), Vector(0, 1), Vector(0, 0)};
  Polygon subject_polygon = Polygon(square_vertices);
  Polygon disk = generate_disk(Vector(0.5, 0.5), 0.7, num_edge);
  subject_polygon.clip_polygon(disk);
  std::vector<Polygon> pols;
  std::vector<Vector> ps;
  pols.push_back(subject_polygon);
  // pols.push_back(disk);
  ps.push_back(Vector(0.5, 0.5));
  save_svg(pols, ps, "disk_intersect_square.svg");
}

void test_intersect_disk()
{
  const int n_p = 30;
  std::vector<Vector> points;
  for (int i = 0; i < n_p; i++)
  {
    double pos1 = generate();
    double pos2 = generate();
    points.push_back(Vector(pos1, pos2));
  }
  double weights[n_p + 1];
  for (int i = 0; i < n_p + 1; i++)
  {
    if (i < n_p)
    {
      weights[i] = 0.3;
    }
    else
    {
      weights[i] = 0.05;
    }
  }
  std::vector<Polygon> p = intersect_disks(points, weights, 6);
  std::vector<Polygon> p1 = compute_voronoi(points, weights);
  save_svg(p, points, "test_intersect_disk.svg");
  save_svg(p1, points, "test_voronoi_disk.svg");
}

void test_fluid_evaluate()
{
  const int n_points = 10;
  double desired_fluid_volume = 0.1;
  std::vector<Vector> positions;
  for (int i = 0; i < n_points; i++)
  {
    double pos1 = generate();
    double pos2 = generate();
    positions.push_back(Vector(pos1, pos2));
  }

  double weights[n_points + 1];
  for (int i = 0; i < n_points; i++)
  {
    double w = generate_weight();
    weights[i] = w;
    std::cout << "Weight of fluid " << i << "th: " << weights[i] << std::endl;
  }
  weights[n_points] = 0.6;
  std::cout << "Weight of air: " << weights[n_points] << std::endl;

  fluid_data *params = new fluid_data();
  params->points = &positions;
  params->dfv = desired_fluid_volume;

  int ret = lbfgs(n_points + 1, weights, NULL, fluid_evaluate, progress, (void *)params, nullptr);
  std::cout << "Status of lbfgs: " << ret << std::endl;
  std::vector<Polygon> output_polygons = intersect_disks(positions, weights, 6);
  save_svg(output_polygons, positions, "test_fluid.svg");
}

// Lab 9: Tutte-embedding
std::vector<Vector> tutte_embedding(std::vector<Vector> vertices, std::vector<Vector> boundaries, std::vector<std::vector<Vector>> adj, int num_iter, std::vector<int> interior_indexes, std::vector<int> boundary_indexes) {
  double s = 0;
  for (int i=0; i<boundaries.size(); i++) {
    s += norm(boundaries[i+1] - boundaries[i]);
  }
  double cs = 0;
  for (int i = 0; i< boundaries.size(); i++) {
    double theta = 2*M_PI*cs/s;
    vertices[i] = Vector(cos(theta), sin(theta));
    cs += norm(boundaries[i+1] - boundaries[i]);
  }
  for (int i=0; i<num_iter; i++) {
    std::vector<Vector> current_vertices = vertices;
    for (int v=0; v < interior_indexes.size(); v++) {
      int current_index = interior_indexes[v];
      current_vertices[current_index] = Vector(0, 0);
      int num_neighbor = adj[current_index].size();
      for (int n = 0; n < num_neighbor; n++) {
        current_vertices[current_index] = current_vertices[current_index] + adj[current_index][n];
      }
      current_vertices[current_index] = current_vertices[current_index] / num_neighbor;
    }
    for (int v=0; v < boundary_indexes.size(); v++) {
      current_vertices[boundary_indexes[v]] = vertices[boundary_indexes[v]];
    }
    vertices = current_vertices;
  }
  return vertices;
}


int main(int argc, char **argv)
{
  // The 3 following functions are created to test helping functions

  // test_intersect_disk();
  // test_intersect_disk();
  // test_fluid_evaluate();

  // The 3 following functions are created to run the 3 labs, in order: lab 6, lab 7, lab 8

  // generate_voronoi();
  // generate_optimized_voronoi();
  fluid_simulation();
  
  return 0;
}