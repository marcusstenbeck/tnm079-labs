/*************************************************************************************************
*
* Modeling and animation (TNM079) 2007
* Code base for lab assignments. Copyright:
*   Gunnar Johansson (gunnar.johansson@itn.liu.se)
*   Ken Museth (ken.museth@itn.liu.se)
*   Michael Bang Nielsen (bang@daimi.au.dk)
*   Ola Nilsson (ola.nilsson@itn.liu.se)
*   Andreas Sˆderstrˆm (andreas.soderstrom@itn.liu.se)
*
*************************************************************************************************/
#include "LevelSetOperator.h"


/*! Computes the squares of the partial derivatives in x, y, z using the
 * Godunov method
 *
 * \f[
 * \left( \dfrac{\partial\phi}{\partial x} \right)^2 \approx \begin{cases}
 * \text{max}\left[\text{max}(\phi_x^-,0)^2,       \text{min}(\phi_x^+,0)^2\right] & F > 0 \\
 * \text{max}\left[\text{min}(\phi_x^-,0)^2,       \text{max}(\phi_x^+,0)^2\right] & F < 0 \\
 * \end{cases}
 * \f]
 *
 * \param[in] i grid x coordinate
 * \param[in] j grid y coordinate
 * \param[in] k grid z coordinate
 * \param[in] a speed function
 * \param[out] ddx2 (dphi/dx)^2
 * \param[out] ddy2 (dphi/dy)^2
 * \param[out] ddz2 (dphi/dy)^2
 */
void LevelSetOperator::Godunov(unsigned int i, unsigned int j, unsigned int k, float a,
                               float & ddx2, float & ddy2, float & ddz2)
{
  float ddxm = mLS->DiffXm(i,j,k);
  float ddxp = mLS->DiffXp(i,j,k);
  float ddym = mLS->DiffYm(i,j,k);
  float ddyp = mLS->DiffYp(i,j,k);
  float ddzm = mLS->DiffZm(i,j,k);
  float ddzp = mLS->DiffZp(i,j,k);

  if (a > 0) {
    ddx2 = std::max( std::pow(std::max(ddxm,0.0f),2), std::pow(std::min(ddxp,0.0f),2) );
    ddy2 = std::max( std::pow(std::max(ddym,0.0f),2), std::pow(std::min(ddyp,0.0f),2) );
    ddz2 = std::max( std::pow(std::max(ddzm,0.0f),2), std::pow(std::min(ddzp,0.0f),2) );
  }
  else {
    ddx2 = std::max( std::pow(std::min(ddxm,0.0f),2), std::pow(std::max(ddxp,0.0f),2) );
    ddy2 = std::max( std::pow(std::min(ddym,0.0f),2), std::pow(std::max(ddyp,0.0f),2) );
    ddz2 = std::max( std::pow(std::min(ddzm,0.0f),2), std::pow(std::max(ddzp,0.0f),2) );
  }
}


void LevelSetOperator::IntegrateEuler(float dt)
{
  // Create grid used to store next time step
  LevelSetGrid grid = GetGrid();

  // Iterate over grid and compute the grid values for the next timestep
  LevelSetGrid::Iterator iter = GetGrid().BeginNarrowBand();
  LevelSetGrid::Iterator iend = GetGrid().EndNarrowBand();
  while (iter != iend) {
    unsigned int i = iter.GetI();
    unsigned int j = iter.GetJ();
    unsigned int k = iter.GetK();

    // Compute rate of change
    float ddt = Evaluate(i,j,k);

    // Compute the next time step and store it in the grid
    grid.SetValue(i,j,k, GetGrid().GetValue(i,j,k) + ddt*dt);

    iter++;
  }

  // Update the grid with the next time step
  GetGrid() = grid;
}


void LevelSetOperator::IntegrateRungeKutta(float dt)
{
  // Advance the solution one time step (dt) using the Runge-Kutta scheme
  // Hint: This scheme is a sequence of Euler steps..
}

