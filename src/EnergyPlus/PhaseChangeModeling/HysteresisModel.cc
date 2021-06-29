// EnergyPlus, Copyright (c) 1996-2020, The Board of Trustees of the University of Illinois,
// The Regents of the University of California, through Lawrence Berkeley National Laboratory
// (subject to receipt of any required approvals from the U.S. Dept. of Energy), Oak Ridge
// National Laboratory, managed by UT-Battelle, Alliance for Sustainable Energy, LLC, and other
// contributors. All rights reserved.
//
// NOTICE: This Software was developed under funding from the U.S. Department of Energy and the
// U.S. Government consequently retains certain rights. As such, the U.S. Government has been
// granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable,
// worldwide license in the Software to reproduce, distribute copies to the public, prepare
// derivative works, and perform publicly and display publicly, and to permit others to do so.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
// (1) Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//
// (2) Redistributions in binary form must reproduce the above copyright notice, this list of
//     conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//
// (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory,
//     the University of Illinois, U.S. Dept. of Energy nor the names of its contributors may be
//     used to endorse or promote products derived from this software without specific prior
//     written permission.
//
// (4) Use of EnergyPlus(TM) Name. If Licensee (i) distributes the software in stand-alone form
//     without changes from the version obtained under this License, or (ii) Licensee makes a
//     reference solely to the software portion of its product, Licensee must refer to the
//     software as "EnergyPlus version X" software, where "X" is the version number Licensee
//     obtained under this License and may not use a different name for the software. Except as
//     specifically required in this Section (4), Licensee shall not use in a company name, a
//     product name, in advertising, publicity, or other promotional activities any name, trade
//     name, trademark, logo, or other designation of "EnergyPlus", "E+", "e+" or confusingly
//     similar designation, without the U.S. Department of Energy's prior written consent.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
// AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include <ObjexxFCL/Array1D.hh>

#include <EnergyPlus/DataHeatBalance.hh>
#include <EnergyPlus/DataIPShortCuts.hh>
#include <EnergyPlus/EnergyPlus.hh>
#include <EnergyPlus/InputProcessing/InputProcessor.hh>
#include <EnergyPlus/PhaseChangeModeling/HysteresisModel.hh>
#include <EnergyPlus/UtilityRoutines.hh>

#include <cmath>
#pragma warning(disable : 4756)

namespace EnergyPlus {

namespace HysteresisPhaseChange {

    bool getHysteresisModels(true);
    int numHysteresisModels = 0;
    std::vector<HysteresisPhaseChange> hysteresisPhaseChangeModels;

    HysteresisPhaseChange *HysteresisPhaseChange::factory(const std::string &objectName)
    {
        if (getHysteresisModels) {
            readAllHysteresisModels();
            getHysteresisModels = false;
        }
        for (auto &hm : hysteresisPhaseChangeModels) {
            if (hm.name == objectName) {
                return &hm;
            }
        }
        // because of the passive linking between materials and material property objects,
        // we don't know ahead of time for sure whether we will have a material property
        // so we can't return fatal here if it isn't found, just leave it null
        return nullptr;
    }

    Real64 HysteresisPhaseChange::getPhaseFraction(Real64 T, Real64 Tc, Real64 tau1, Real64 tau2,bool ForM)
    {
        Real64 const bigDiff(50); // To avoid float-point underflow
        
        // Looks up the enthalpy on the characteristic curve defined by the parameters Tc, tau1, and tau2,
        // and the position on that curve defined by T.
        Real64 eta1 = std::abs(T - Tc) / tau1;
        Real64 eta2 = std::abs(T - Tc) / tau2;
        Real64 eta;
        if (eta1 > bigDiff) {
            eta1 = bigDiff;
        }
        if (eta2 > bigDiff) {
            eta2 = bigDiff;
        }
        eta1 = (1.0/ 2.0) * exp(-2 * eta1);
        eta2 = (1.0/ 2.0) * exp(-2 * eta2);            
        if (T <= Tc) {
            eta = eta1;   // Bug
        } else {
            eta = 1 - eta2;
        }
        /*
        if (ForM) {
            if (T >= Tc + tau2 - 1.4) {
                eta = 1;
            } else if (T <= Tc + tau2 - 1.5) {
                eta = eta;
            } else {
                eta1 = this->getPhaseFraction(Tc + tau2 - 1.5, Tc, tau1, tau2, ForM);
                eta2 = this->getPhaseFraction(Tc + tau2 - 1.4, Tc, tau1, tau2, ForM);
                eta = (T - Tc - tau2 + 1.5) / 0.1 * (eta2 - eta1) + eta1;
            }
        }
        */
        return eta;
    }


    Real64 HysteresisPhaseChange::getCurrentSpecificHeat(
        Real64 prevTempTD, Real64 updatedTempTDT, Real64 prevEnth,Real64 &updatedEnth, Real64 prevPhaseFraction,Real64 &updatedPhaseFraction, int prevPhaseChangeState, int &phaseChangeState)
    {
        // Main public facing function; returns the current specific heat based on input properties, and current and previous conditions.
        // In a future version, this could be compartmentalized to track all states and histories, but it would require some further modification to
        // the HBFDManager
        Real64 Tc;   // assigned later
        Real64 Tau1; // assigned later
        Real64 Tau2; // assigned later
        Real64 Cp;
        Real64 phaseChangeDeltaT = updatedTempTDT - prevTempTD; 
        Real64 enth_PrevTemp; //assigned later

        Real64 prevPhaseFraction_Freezing; // assigned later
        Real64 UpdatedPhaseFraction_Freezing;// assigned later

        Real64 prevPhaseFraction_Melting; //
        Real64 UpdatedPhaseFraction_Melting;//
        Real64 const smalldiff(1.e-3);
        bool ForM;
        // determine phase fraction and curve characteristics based on delta T direction, updated temp, and previous state
        if (phaseChangeDeltaT <= 0){ //Decreasing
            ForM = 1;
            Tc = this->peakTempFreezing;   // Freezing curve
            Tau1 = this->deltaTempFreezingLow;
            Tau2 = this->deltaTempFreezingHigh;
            
            prevPhaseFraction_Freezing = this->getPhaseFraction(prevTempTD, Tc, Tau1, Tau2, ForM);
            if (std::abs(prevPhaseFraction_Freezing)<smalldiff){ // previous state is solid.Otherwise, the denominator will be "zero!"
                updatedPhaseFraction = 0;    
            }else{
                UpdatedPhaseFraction_Freezing = this->getPhaseFraction(updatedTempTDT, Tc, Tau1, Tau2,ForM);
                if (std::abs(UpdatedPhaseFraction_Freezing) < smalldiff) {
                    UpdatedPhaseFraction_Freezing = 0;
                }
                updatedPhaseFraction = prevPhaseFraction/prevPhaseFraction_Freezing*UpdatedPhaseFraction_Freezing;
            }
            phaseChangeState = 0;  //always 0
        }else{ //Increasing
            ForM = 0;
            Tc = this->peakTempMelting;   // Melting curve
            Tau1 = this->deltaTempMeltingLow;
            Tau2 = this->deltaTempMeltingHigh;

            prevPhaseFraction_Melting = this->getPhaseFraction(prevTempTD, Tc, Tau1, Tau2,ForM);
            if (std::abs(1-prevPhaseFraction_Melting)<smalldiff){ //previous state is liquid!
                updatedPhaseFraction = 1;  
            }else{
                UpdatedPhaseFraction_Melting = this->getPhaseFraction(updatedTempTDT, Tc, Tau1, Tau2,ForM);
                updatedPhaseFraction = 1 - (1-prevPhaseFraction)/(1-prevPhaseFraction_Melting)*(1-UpdatedPhaseFraction_Melting);
            }
            phaseChangeState = 0;  //always 0
        }
        // then calculate the specific heat and return it
        Cp = this->specHeat(prevTempTD, updatedTempTDT, prevPhaseFraction, updatedPhaseFraction);

        this->CpOld = Cp;

        updatedEnth = prevEnth + Cp*(updatedTempTDT-prevTempTD);
        if (updatedEnth > 120000) {
            int pause = 1; // ForDebug
        }
        return Cp;
    }

    Real64 HysteresisPhaseChange::specHeat(Real64 temperaturePrev,
                                           Real64 temperatureCurrent,
                                           Real64 prevPhaseFraction,
                                           Real64 updatedPhaseFraction)
    {

        //	Tc                  ! Critical (Melting/Freezing) Temperature of PCM
        //	Tau1                ! Width of Melting Zone low
        //	Tau2                ! Width of Melting Zone high
        //	EnthalpyOld         ! Previous Timestep Nodal Enthalpy
        //	EnthalpyNew         ! Current Timestep Nodal Enthalpy
        Real64 const smalldiff(1.e-8);
        Real64 T = temperatureCurrent;

        Real64 Cp;
        
        
        if (std::abs(T - temperaturePrev) < smalldiff) {
            if ((std::abs(prevPhaseFraction - updatedPhaseFraction) < smalldiff) | (T == temperaturePrev)) {
                Cp = this->CpOld;    
            } else {
                Cp = this->specificHeatLiquid * updatedPhaseFraction + (1 - updatedPhaseFraction) * this->specificHeatSolid +
                     (updatedPhaseFraction - prevPhaseFraction) / (temperatureCurrent - temperaturePrev) * this->totalLatentHeat;
            }
            
            //Cp = this->specificHeatLiquid * updatedPhaseFraction +(1-updatedPhaseFraction)*this->specificHeatSolid;
        }else{
            Cp = this->specificHeatLiquid * updatedPhaseFraction +(1-updatedPhaseFraction)*this->specificHeatSolid +
                (updatedPhaseFraction - prevPhaseFraction)/(temperatureCurrent-temperaturePrev)*this->totalLatentHeat;
        }
        return Cp;
    }

    Real64 HysteresisPhaseChange::getConductivity(Real64 T)
    {
        if (T < this->peakTempMelting) {
            return this->fullySolidThermalConductivity;
        } else if (T > this->peakTempFreezing) {
            return this->fullyLiquidThermalConductivity;
        } else {
            return (this->fullySolidThermalConductivity + this->fullyLiquidThermalConductivity) / 2.0;
        }
    }

    Real64 HysteresisPhaseChange::getDensity(Real64 T)
    {
        if (T < this->peakTempMelting) {
            return this->fullySolidDensity;
        } else if (T > this->peakTempFreezing) {
            return this->fullyLiquidDensity;
        } else {
            return (this->fullySolidDensity + this->fullyLiquidDensity) / 2.0;
        }
    }

    void readAllHysteresisModels()
    {

        // convenience variables
        DataIPShortCuts::cCurrentModuleObject = "MaterialProperty:PhaseChangeHysteresis";
        numHysteresisModels = inputProcessor->getNumObjectsFound(DataIPShortCuts::cCurrentModuleObject);

        // loop over all hysteresis input instances, if zero, this will simply not do anything
        for (int hmNum = 1; hmNum <= numHysteresisModels; ++hmNum) {

            // just a few vars to pass in and out to GetObjectItem
            int ioStatus;
            int numAlphas;
            int numNumbers;

            // get the input data and store it in the Shortcuts structures
            inputProcessor->getObjectItem(DataIPShortCuts::cCurrentModuleObject,
                                          hmNum,
                                          DataIPShortCuts::cAlphaArgs,
                                          numAlphas,
                                          DataIPShortCuts::rNumericArgs,
                                          numNumbers,
                                          ioStatus,
                                          DataIPShortCuts::lNumericFieldBlanks,
                                          DataIPShortCuts::lAlphaFieldBlanks,
                                          DataIPShortCuts::cAlphaFieldNames,
                                          DataIPShortCuts::cNumericFieldNames);

            // the input processor validates the numeric inputs based on the IDD definition
            // still validate the name to make sure there aren't any duplicates or blanks
            // blanks are easy: fatal if blank
            if (DataIPShortCuts::lAlphaFieldBlanks[0]) {
                ShowFatalError("Invalid input for " + DataIPShortCuts::cCurrentModuleObject + " object: Name cannot be blank");
            }

            // we just need to loop over the existing vector elements to check for duplicates since we haven't add this one yet
            for (auto &existingHysteresisModel : hysteresisPhaseChangeModels) {
                if (DataIPShortCuts::cAlphaArgs(1) == existingHysteresisModel.name) {
                    ShowFatalError("Invalid input for " + DataIPShortCuts::cCurrentModuleObject +
                                   " object: Duplicate name found: " + existingHysteresisModel.name);
                }
            }

            // now build out a new hysteresis instance and add it to the vector
            HysteresisPhaseChange thisHM;
            thisHM.name = DataIPShortCuts::cAlphaArgs(1);
            thisHM.totalLatentHeat = DataIPShortCuts::rNumericArgs(1);
            thisHM.fullyLiquidThermalConductivity = DataIPShortCuts::rNumericArgs(2);
            thisHM.fullyLiquidDensity = DataIPShortCuts::rNumericArgs(3);
            thisHM.specificHeatLiquid = DataIPShortCuts::rNumericArgs(4);
            thisHM.deltaTempMeltingHigh = DataIPShortCuts::rNumericArgs(5);
            thisHM.peakTempMelting = DataIPShortCuts::rNumericArgs(6);
            thisHM.deltaTempMeltingLow = DataIPShortCuts::rNumericArgs(7);
            thisHM.fullySolidThermalConductivity = DataIPShortCuts::rNumericArgs(8);
            thisHM.fullySolidDensity = DataIPShortCuts::rNumericArgs(9);
            thisHM.specificHeatSolid = DataIPShortCuts::rNumericArgs(10);
            thisHM.deltaTempFreezingHigh = DataIPShortCuts::rNumericArgs(11);
            thisHM.peakTempFreezing = DataIPShortCuts::rNumericArgs(12);
            thisHM.deltaTempFreezingLow = DataIPShortCuts::rNumericArgs(13);
            thisHM.specHeatTransition = (thisHM.specificHeatSolid + thisHM.specificHeatLiquid) / 2.0;
            thisHM.CpOld = thisHM.specificHeatSolid;
            hysteresisPhaseChangeModels.push_back(thisHM);
        }
    }

    void clear_state()
    {
        numHysteresisModels = 0;
        getHysteresisModels = true;
        hysteresisPhaseChangeModels.clear();
    }

} // namespace HysteresisPhaseChange

} // namespace EnergyPlus
