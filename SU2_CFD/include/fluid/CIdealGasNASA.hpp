/*!
 * \file CIdealGasNASA.hpp
 * \brief Defines the ideal gas with NASA 9-coefficient polynomial model.
 * \author E. L. Gevard
 * \version 8.6.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "CFluidModel.hpp"

using namespace std;

/*!
 * \class CIdealGasNASA
 * \brief Child class for defining the ideal gas with NASA 9-coefficient polynomial model
 * \author: E. L. Gevard
 */

class CIdealGasNASA final : public CFluidModel{
 public:
  CIdealGasNASA(su2double val_gas_constant, su2double val_operating_pressure, su2double val_Temperature_Ref);

  /*!
   * \brief Set the NASA temperature polynomial coefficients for variable Cp.
   * \param[in] config - configuration container for the problem.
   * \param[in] val_Temperature_Ref - reference temperature value
  */
  void SetCpModel(const CConfig* config, su2double val_Temperature_Ref) override;

  /*!
   * \brief Set the Dimensionless State using Temperature.
   * \param[in] t - Temperature value at the point.
   */
  void SetTDState_T(su2double t, const su2double *val_scalars = nullptr) override;

  /*!
   * \brief Set the Dimensionless State using enthalpy.
   * \param[in] val_enthalpy - Enthalpy value at the point.
   */
  void SetTDState_h(su2double val_enthalpy, const su2double* val_scalars = nullptr);
  
 private:
  su2double Gas_Constant{0.0}; /*!< \brief Specific Gas Constant. */
  su2double Gamma{0.0};        /*!< \brief Ratio of specific heats. */
  su2double Std_Ref_Temp_ND{0.0}; /*!< \brief Nondimensional standard reference temperature for enthalpy. */
  array<su2double, N_NASA_COEFFS> coeffs_; /*!< \brief Polynomial coefficients for heat capacity as a function of temperature. */
  su2double Temperature_Min{0};   /*!< \brief Minimum temperature value for the polynomial model. */
  su2double Temperature_Max{5e3};   /*!< \brief Maximum temperature value for the polynomial model. */
  su2double Enthalpy_Ref{0}; /*!< \brief Nondimensional standard reference enthalpy*/

  /*!
   * \brief Compute dimensionless Cp and Enthalpy given a dimensionless temperature
   * \param[in] t - static temperature
   * \param[out] cp - reference to specific heat Cp
   * \param[out] h - reference to specific enthalpy h
  */
  void Compute_CpH(su2double t, su2double& cp, su2double& h);
};
