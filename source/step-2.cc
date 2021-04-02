/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * based on deal.II step-2
 */


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

#include <fstream>

using namespace dealii;


void
make_grid(Triangulation<2> &triangulation)
{
  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 5);

  static const SphericalManifold<2> manifold_description(center);
  triangulation.set_all_manifold_ids(0);
  triangulation.set_manifold(0, manifold_description);

  for (unsigned int step = 0; step < 3; ++step)
    {
      Triangulation<2>::active_cell_iterator cell =
                                               triangulation.begin_active(),
                                             endc = triangulation.end();

      for (; cell != endc; ++cell)
        for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v)
          {
            const double distance_from_center =
              center.distance(cell->vertex(v));

            if (std::fabs(distance_from_center - inner_radius) < 1e-10)
              {
                cell->set_refine_flag();
                break;
              }
          }

      triangulation.execute_coarsening_and_refinement();
    }
}

void
make_square_grid(Triangulation<2> &triangulation)
{
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);
}



SparsityPattern
distribute_dofs(DoFHandler<2> &dof_handler)
{
  static const FE_Q<2> finite_element(1);
  dof_handler.distribute_dofs(finite_element);

  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());

  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);

  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);

  std::ofstream out("sparsity_pattern1.svg");
  sparsity_pattern.print_svg(out);

  return sparsity_pattern;
}



SparsityPattern
renumber_dofs(DoFHandler<2> &dof_handler)
{
  DoFRenumbering::Cuthill_McKee(dof_handler);

  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);

  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);

  std::ofstream out("sparsity_pattern2.svg");
  sparsity_pattern.print_svg(out);

  return sparsity_pattern;
}

void
row_lenths(const SparsityPattern &sparsity_pattern)
{
  for (unsigned int i = 0; i < sparsity_pattern.n_rows(); i++)
    {
      std::cout << "Row" << i << "- - row length "
                << sparsity_pattern.row_length(i) << std::endl;
    }
}

std::tuple<int, int, double, double>
compute_pattern_statistics(const SparsityPattern &sparsity_pattern)
{
  double avarage_per_row = 0.;
  for (unsigned int i = 0; i < sparsity_pattern.n_rows(); i++)
    {
      avarage_per_row += sparsity_pattern.row_length(i);
    }

  double fill_ratio = avarage_per_row / (double)(sparsity_pattern.n_rows() *
                                                 sparsity_pattern.n_cols());
  avarage_per_row /= (double)sparsity_pattern.n_rows();
  return std::make_tuple(sparsity_pattern.n_rows(),
                         sparsity_pattern.bandwidth(),
                         avarage_per_row,
                         fill_ratio);
}


int
main()
{
  Triangulation<2> triangulation;
  // make_grid(triangulation);
  make_square_grid(triangulation);

  DoFHandler<2> dof_handler(triangulation);

  auto sparcity_pattern  = distribute_dofs(dof_handler);
  auto sparcity_pattern2 = renumber_dofs(dof_handler);

  row_lenths(sparcity_pattern2);

  std::cout << sparcity_pattern.row_length(41) << std::endl;

  auto patern_statistics = compute_pattern_statistics(sparcity_pattern2);
  std::cout << std::get<0>(patern_statistics) << " "
            << std::get<1>(patern_statistics) << "  "
            << std::get<2>(patern_statistics) << "  "
            << std::get<3>(patern_statistics) << std::endl;
}
