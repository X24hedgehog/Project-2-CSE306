// std::vector<Polygon> Laguerre_polygons = compute_voronoi(p, x);

// for (int i = 0; i < Laguerre_polygons.size(); i++)
// {
//   Polygon clip_polygon = generate_polygon(p[i], sqrt(x[i] - x[n-1]));
//   Laguerre_polygons[i].clip_polygon(clip_polygon);
//   s.push_back(Laguerre_polygons[i].get_area());
// }

//   Polygon generate_polygon(Vector center, double radius, int num_edges)
// {
//   std::vector<Vector> vertices;
//   for (int i = 0; i < num_edges; i++)
//   {
//     Vector current_vertext;
//     current_vertext[0] = radius * cos(2 * i * M_PI / num_edges);
//     current_vertext[1] = radius * sin(2 * i * M_PI / num_edges);
//     vertices.push_back(current_vertext);
//   }
//   Polygon out_polygon = Polygon(vertices);
//   return out_polygon;
// }