#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#include <math.h>
#include <ctime>
#include <iostream>
#define M_PI 3.14159265358979323846

#include <string>
#include <stdio.h>
#include <algorithm>
#include <random>
#include <omp.h>
#include "lbfgs.h"

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
  // std::cout << M[0] << " " << M[1] << std::endl;
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
  // std::cout << M[0] << " " << M[1] << std::endl;
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
    for (int i = 0; i < vertices.size(); i++)
    {
      std::cout << vertices[i][0] << "    " << vertices[i][1] << std::endl;
    }
    std::cout << "End of a polygon" << std::endl;
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
        // std::cout << "curVertex inside" << std::endl;
        if (!check_inside_voronoi(prevVertex, edge, w1, w2))
        {
          // std::cout << "prevVertex not inside" << std::endl;
          outPolygon.vertices.push_back(intersection);
        }
        outPolygon.vertices.push_back(curVertex);
      }
      else if (check_inside_voronoi(prevVertex, edge, w1, w2))
      {
        // std::cout << "curVertext not inside but prevVertext inside" << std::endl;
        outPolygon.vertices.push_back(intersection);
      }
      else
      {
        // std::cout << "Both not inside" << std::endl;
      }
    }
    this->vertices = outPolygon.vertices;
    // this->print_vertices();
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
      std::pair<Vector, Vector> clipEdge(c.vertices[i], c.vertices[prev_index]);
      outPolygon = Polygon();
      for (int i = 0; i < vertices.size(); i++)
      {
        Vector curVertex = vertices[i];
        Vector prevVertex = vertices[(i > 0) ? (i - 1) : vertices.size() - 1];
        Vector intersection = check_intersect(prevVertex, curVertex, clipEdge);
        if (check_inside(curVertex, clipEdge))
        {
          if (!check_inside(prevVertex, clipEdge))
            outPolygon.vertices.push_back(intersection);
          outPolygon.vertices.push_back(curVertex);
        }
        else if (check_inside(prevVertex, clipEdge))
          outPolygon.vertices.push_back(intersection);
      }
    }
    *this = outPolygon;
  }
};

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
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none")
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
  fprintf(f, "</svg>\n");
  fclose(f);
}

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
  // polygons_list.push_back(subject_polygon);
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

static lbfgsfloatval_t evaluate(
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
  // std::cout << "Number of polygons is: " << polygons.size() << std::endl;
  // std::cout << " " << std::endl;
  // for (int c = 0; c < polygons.size(); c++)
  // {
  //   polygons[c].print_vertices();
  // }
  // std::cout << " " << std::endl;
  // std::cout << " " << std::endl;

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
    // polygons[i].print_vertices();
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
      // std::cout << "Root vertice: " << o[0] << " " << o[1] << std::endl;
      // std::cout << " " << std::endl;
      polygon_integral = 0;
      for (int j = 1; j + 1 < num_vertices; j++)
      {
        Vector a = current_vertices[j];
        Vector b = current_vertices[j + 1];
        std::vector<Vector> triangle_vertices;
        triangle_vertices.push_back(o);
        triangle_vertices.push_back(a);
        triangle_vertices.push_back(b);
        // std::cout << "Vertices of triangle:" << std::endl;
        // for (int ver = 0; ver < 3; ver++)
        // {
        //   std::cout << triangle_vertices[ver][0] << " " << triangle_vertices[ver][1] << std::endl;
        // }
        // std::cout << " " << std::endl;
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
    // std::cout << "Real area: " << s[i] << ", Wanted area: " << l[i] << std::endl;
    // std::cout << " " << std::endl;
  }
  // std::cout << fx << std::endl;
  // std::cout << "another iteration over" << std::endl;
  // std::cout << "another iteration over" << std::endl;
  // std::cout << "another iteration over" << std::endl;
  // std::cout << "another iteration over" << std::endl;
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

std::random_device my_random_device;
std::default_random_engine my_random_engine(my_random_device());
std::uniform_real_distribution<float> my_distribution(0, 1);

double generate()
{
  return my_distribution(my_random_engine);
}

void generate_diagram()
{
  const int n_points = 2000;
  std::vector<Vector> points;
  for (int i = 0; i < n_points; i++)
  {
    double pos1 = generate();
    double pos2 = generate();
    // std::cout << pos1 << "    " << pos2 << std::endl;
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

  std::vector<Polygon> polygons = compute_voronoi(points, balance_weights);
  std::vector<Polygon> weighted_polygons = compute_voronoi(points, weights);

  data *params = new data();
  params->points = &points;
  params->lambdas = lambda_arr;
  int ret = lbfgs(n_points, weights, NULL, evaluate, progress, (void *)params, nullptr);
  std::cout << ret << std::endl;
  std::vector<Polygon> lbfsg_weighted_polygons = compute_voronoi(points, weights);

  save_svg(polygons, points, "balance_image.svg");
  save_svg(weighted_polygons, points, "weighted_image.svg");
  save_svg(lbfsg_weighted_polygons, points, "lbfsg_weighted_image.svg");
}

int main(int argc, char **argv)
{
  generate_diagram();
  return 0;
}