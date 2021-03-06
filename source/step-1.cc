/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * based on deal.II step-1
 */


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <cmath>
#include <fstream>
#include <iostream>

using namespace dealii;



std::tuple<int, int, int>
grid_parameters(const Triangulation<2> &tria)
{
  return std::make_tuple(tria.n_levels(),
                         tria.n_cells(),
                         tria.n_active_cells());
}


void
first_grid()
{
  Triangulation<2> triangulation;

  GridGenerator::hyper_cube(triangulation);
  std::cout << "Number of original vertices:" << triangulation.n_vertices()
            << std::endl;

  triangulation.refine_global(4);

  std::cout << "Number of original vertices after 4 refinements:"
            << triangulation.n_vertices() << std::endl;

  {
    std::ofstream out("grid-1.svg");
    GridOut       grid_out;
    grid_out.write_svg(triangulation, out);
    std::cout << "Grid written to grid-1.svg" << std::endl;
  }

  {
    std::ofstream out("grid-1.vtk");
    GridOut       grid_out;
    grid_out.write_vtk(triangulation, out);
    std::cout << "Grid written to grid-1.vtk" << std::endl;
  }
  auto params = grid_parameters(triangulation);
  std::cout << std::get<0>(params) << " " << std::get<1>(params) << "  "
            << std::get<2>(params) << std::endl;
}



void
second_grid()
{
  Triangulation<2> triangulation;

  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 10);

  triangulation.reset_all_manifolds();
  for (unsigned int step = 0; step < 5; ++step)
    {
      std::ofstream out("grid-2-" + std::to_string(step) + ".vtk");
      GridOut       grid_out;
      grid_out.write_vtk(triangulation, out);

      for (auto &cell : triangulation.active_cell_iterators())
        {
          for (const auto v : cell->vertex_indices())
            {
              const double distance_from_center =
                center.distance(cell->vertex(v));

              if (std::fabs(distance_from_center - inner_radius) <=
                  1e-6 * inner_radius)
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }

      triangulation.execute_coarsening_and_refinement();
    }
  {
    std::ofstream out("grid-2.svg");
    GridOut       grid_out;
    grid_out.write_svg(triangulation, out);
    std::cout << "Grid written to grid-2.svg" << std::endl;
  }
  auto params = grid_parameters(triangulation);
  std::cout << std::get<0>(params) << " " << std::get<1>(params) << "  "
            << std::get<2>(params) << std::endl;
}

void
third_grid()
{
  Triangulation<2> triangulation;
  GridGenerator::hyper_L(triangulation);
  std::cout << "Number of original vertices:" << triangulation.n_vertices()
            << std::endl;

  triangulation.refine_global(1);

  std::cout << "Number of original vertices after 1 refinement:"
            << triangulation.n_vertices() << std::endl;


  {
    std::ofstream out("grid-3.vtk");
    GridOut       grid_out;
    grid_out.write_vtk(triangulation, out);
    std::cout << "Grid written to grid-3.vtk" << std::endl;
  }

  const Point<2> corner(0, 0);
  for (unsigned int step = 0; step < 5; ++step)
    {
      std::ofstream out("grid-3-" + std::to_string(step) + ".vtk");
      GridOut       grid_out;
      grid_out.write_vtk(triangulation, out);

      for (auto &cell : triangulation.active_cell_iterators())
        {
          for (const auto v : cell->vertex_indices())
            {
              const double distance_from_center =
                corner.distance(cell->vertex(v));

              if (std::fabs(distance_from_center) <= 1. / 3.)
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }

      triangulation.execute_coarsening_and_refinement();
    }

  auto params = grid_parameters(triangulation);
  std::cout << std::get<0>(params) << " " << std::get<1>(params) << "  "
            << std::get<2>(params) << std::endl;
}


void
circle_grid()
{
  Triangulation<2> triangulation;
  GridGenerator::hyper_ball<2>(triangulation);
  // triangulation.set_all_manifold_ids(0);
  triangulation.set_manifold(0, SphericalManifold<2>());
  const Point<2> mesh_center;
  for (const auto &cell : triangulation.active_cell_iterators())
    if (mesh_center.distance(cell->center()) > cell->diameter() / 10)
      cell->set_all_manifold_ids(0);
  triangulation.refine_global(2);
  std::ofstream out("circle.vtk");
  GridOut       grid_out;
  grid_out.write_vtk(triangulation, out);
}


void
torus_grid()
{
  Triangulation<2, 3> triangulation;
  GridGenerator::torus<2, 3>(triangulation, 4, 1);
  triangulation.refine_global(2);
  std::ofstream out("torus.vtk");
  GridOut       grid_out;
  grid_out.write_vtk(triangulation, out);
}

int
main()
{
  first_grid();
  second_grid();
  third_grid();
  circle_grid();
  torus_grid();
}
