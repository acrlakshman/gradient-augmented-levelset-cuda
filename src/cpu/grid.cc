/*
 * Copyright (c) 2019, Lakshman Anumolu, Raunak Bardia
 * All rights reserved.
 *
 * This file is part of gradient-augmented-levelset-cuda project whose
 * distribution is governed by the BSD 2-Clause License contained in the
 * accompanying LICENSE.txt file.
 */

#include "grid.h"

template<typename T, int DIM>
GALS::CPU::Grid<T, DIM>::Grid(int nx, int ny, int nz) : m_dimension(DIM), m_nx(nx), m_ny(ny), m_nz(nz) {
    m_grid.resize((m_nx + 2) * (m_ny + 2) * (m_nz + 2));
}

template<typename T, int DIM>
GALS::CPU::Grid<T, DIM>::Grid(int nx, int ny) : Grid(nx, ny, 1)
{}

template<typename T, int DIM>
GALS::CPU::Grid<T, DIM>::Grid(int nx) : Grid(nx, 1, 1)
{}

template<typename T, int DIM>
GALS::CPU::Grid<T, DIM>::~Grid()
{}

template<typename T, int DIM>
std::vector<T> GALS::CPU::Grid<T, DIM>::x(int i, int j, int k)
{
   std::vector<T> pos(m_dimension);

   for (int c = 0; c < m_dimension; ++c) pos[i] = m_grid[((k * m_nx * m_ny) + (j * m_nx) + i) * m_dimension + c];

   return pos;
}

template<typename T, int DIM>
const int GALS::CPU::Grid<T, DIM>::dimension() const
{
   return m_dimension;
}

template<typename T, int DIM>
const int GALS::CPU::Grid<T, DIM>::size() const
{
   return m_grid.size();
}

template<typename T, int DIM>
const std::vector<int> GALS::CPU::Grid<T, DIM>::getNumCells() const
{
   return std::vector<int>{m_nx, m_ny, m_nz};
}

template<typename T, int DIM>
const int GALS::CPU::Grid<T, DIM>::getPadding() const
{
   return m_pad;
}

template<typename T, int DIM>
void GALS::CPU::Grid<T, DIM>::setPadding(const int pad)
{
   m_pad = pad;
}

template<typename T, int DIM>
void GALS::CPU::Grid<T, DIM>::generate(T x_min, T x_max, T y_min, T y_max, T z_min, T z_max) {
   // TODO
}

template class GALS::CPU::Grid<double, 1>;
template class GALS::CPU::Grid<double, 2>;
template class GALS::CPU::Grid<double, 3>;
