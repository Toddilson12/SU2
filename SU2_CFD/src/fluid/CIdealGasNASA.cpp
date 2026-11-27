/*!
 * \file CIdealGasNASA.cpp
 * \brief Source of the ideal gas with NASA 9-coefficient polynomial model.
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

#include "../../include/fluid/CIdealGasNASA.hpp"

CIdealGasNASA::CIdealGasNASA(su2double val_gas_constant, su2double val_operating_pressure, su2double val_Temperature_Ref)
    : CFluidModel()
{
    /* In the incompressible ideal gas model, the thermodynamic pressure is decoupled
    from the governing equations and held constant. The density is therefore only a
    function of temperature variations. We also use a molecular weight (g/mol) and the
    universal gas constant to compute the specific gas constant for the fluid. The gas
    is incompressible, so Cp = Cv (gamma = 1). */

    Gas_Constant = val_gas_constant;
    Pressure = val_operating_pressure;
    Gamma = 1.0;
    Std_Ref_Temp_ND = val_Temperature_Ref;
}

void CIdealGasNASA::SetCpModel(const CConfig* config, su2double val_Temperature_Ref)
{
    const auto& CpPolyOptions = config->GetCpPolynomialParsedOptions();
    for (int i = 0; i < CpPolyOptions.n_coeffs; ++i) {
        /*--- Note that we use cp = dh/dt here. To transform to mass Cp, each coefficient has to multiplied by 
        non-dim. gas constant. ---*/
        coeffs_[i] = config->GetCp_PolyCoeffND(i) * config->GetGas_ConstantND(); 
    }
    /*--- Non-dim. temperatures ---*/
    Temperature_Min = CpPolyOptions.temperature_limits[0] / val_Temperature_Ref;
    Temperature_Max = CpPolyOptions.temperature_limits[1] / val_Temperature_Ref;

    /*--- Precompute reference enthalpy for sensible enthalpy (incompressible) ---*/
    Compute_CpH(Std_Ref_Temp_ND, Cp, Enthalpy_Ref);
}


void CIdealGasNASA::SetTDState_T(su2double t, const su2double *val_scalars)
{
    /* The EoS only depends upon temperature. */
    if (t < Temperature_Min || t > Temperature_Max)
      cout << "Warning: temperature outside valid range for NASA polynomial.";

    Temperature = t;
    Density = Pressure / (Temperature * Gas_Constant);

    Compute_CpH(Temperature, Cp, Enthalpy);
    Enthalpy -= Enthalpy_Ref;

    Cv = Cp / Gamma;
}


void CIdealGasNASA::SetTDState_h(su2double val_enthalpy, const su2double* val_scalars) 
{
    Enthalpy = val_enthalpy;
    /*--- convergence criterion for temperature in [K], high accuracy needed for restarts. ---*/
    const su2double toll = 1e-5;
    su2double temp_iter = Std_Ref_Temp_ND;
    su2double Cp_iter = 0.0;
    su2double Enthalpy_iter = 0.0;
    su2double delta_temp_iter = 1e10;
    su2double delta_enthalpy_iter;
    const int counter_limit = 30;

    int counter = 0;

    /*--- Compute temperature given enthalpy using Newton-Raphson. ---*/
    while ((abs(delta_temp_iter) > toll) && (counter++ < counter_limit)) {
      /* Evaluate the new Cp and enthalpy from the coefficients and temperature. */
      
      /*--- Compute sensible enthalpy ---*/
      Compute_CpH(temp_iter, Cp_iter, Enthalpy_iter);
      Enthalpy_iter -= Enthalpy_Ref;

      delta_enthalpy_iter = Enthalpy - Enthalpy_iter;
      
      /*--- Calculate clamped temperature if it is outside range. Ensures delta_temp_iter
            is calculated based on the clamped temperature and prevents breaking early due to toll ---*/
      const su2double next_temp = max(min(temp_iter + delta_enthalpy_iter / Cp_iter, Temperature_Max), Temperature_Min);
      delta_temp_iter = next_temp - temp_iter;

      temp_iter = next_temp;
    }
    
    if (temp_iter == Temperature_Min || temp_iter == Temperature_Max) {
      cout << "Warning: temperature outside valid range for NASA polynomials during Newton-Raphson." << endl;
    }

    if (counter == counter_limit) {
      cout << "Warning: Newton-Raphson exceeds max. iterations in temperature computation." << endl;
    }

    Temperature = temp_iter;
    Compute_CpH(Temperature, Cp, Enthalpy_iter); /* use Enthalpy_iter as Enthalpy is already prescribed */

    Density = Pressure / (Temperature * Gas_Constant);
    Cv = Cp / Gamma;
}
  

void CIdealGasNASA::Compute_CpH(su2double t, su2double& cp, su2double& h)
{
    const auto& A = coeffs_;
    const su2double invT = 1.0 / t;

    /*--- T^-2 and T^-1 terms -> A1, A2 ---*/
    cp = A[0] * invT * invT + A[1] * invT;

    h = -A[0] * invT + A[1] * log(t) + A[7]; /*--- A[7] = b1 = enthalpy integration constant ---*/

    /*--- T^0 to T^4 terms -> A3 to A7 ---*/
    su2double T_i = 1.0;

    /*--- Only computes a reference enthalpy, without the integration constant ---*/
    for (unsigned short i = 2; i < 7; i++) {
        cp += A[i] * T_i;

        h += (A[i] * T_i * t) / (i - 1);

        T_i *= t;
    }
}