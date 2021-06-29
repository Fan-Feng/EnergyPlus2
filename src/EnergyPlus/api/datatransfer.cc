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

#include <cmath>

#include <EnergyPlus/api/datatransfer.h>
#include <EnergyPlus/DataEnvironment.hh>
#include <EnergyPlus/DataHeatBalance.hh>
#include <EnergyPlus/DataHVACGlobals.hh>
#include <EnergyPlus/UtilityRoutines.hh>
#include <EnergyPlus/DataGlobals.hh>
#include <EnergyPlus/DataRuntimeLanguage.hh>
#include <EnergyPlus/OutputProcessor.hh>
#include <EnergyPlus/PluginManager.hh>

const char * listAllAPIDataCSV() {
    std::string output;
    output.append("**ACTUATORS**\n");
    for (auto const & availActuator : EnergyPlus::DataRuntimeLanguage::EMSActuatorAvailable) {
        if (availActuator.ComponentTypeName.empty() && availActuator.UniqueIDName.empty() && availActuator.ControlTypeName.empty()) {
            break;
        }
        output.append("Actuator,");
        output.append(availActuator.ComponentTypeName).append(",");
        output.append(availActuator.ControlTypeName).append(",");
        output.append(availActuator.UniqueIDName).append(";\n");
    }
    output.append("**INTERNAL_VARIABLES**\n");
    for (auto const & availVariable : EnergyPlus::DataRuntimeLanguage::EMSInternalVarsAvailable) {
        if (availVariable.DataTypeName.empty() && availVariable.UniqueIDName.empty()) {
            break;
        }
        output.append("InternalVariable,");
        output.append(availVariable.DataTypeName).append(",");
        output.append(availVariable.UniqueIDName).append("\n");
    }
    output.append("**PLUGIN_GLOBAL_VARIABLES**\n");
    for (auto const & gVarName : EnergyPlus::PluginManagement::globalVariableNames) {
        output.append("PluginGlobalVariable,");
        output.append(gVarName).append("\n");
    }
    output.append("**TRENDS**\n");
    for (auto const & trend : EnergyPlus::PluginManagement::trends) {
        output.append("PluginTrendVariable,");
        output.append(trend.name).append("\n");
    }
    output.append("**METERS**\n");
    for (auto const & meter : EnergyPlus::OutputProcessor::EnergyMeters) {
        if (meter.Name.empty()) {
            break;
        }
        output.append("OutputMeter,");
        output.append(meter.Name).append("\n");
    }
    output.append("**VARIABLES**\n");
    for (auto const & variable : EnergyPlus::OutputProcessor::RVariableTypes) {
        if (variable.VarNameOnly.empty() && variable.KeyNameOnlyUC.empty()) {
            break;
        }
        output.append("OutputVariable,");
        output.append(variable.VarNameOnly).append(",");
        output.append(variable.KeyNameOnlyUC).append("\n");
    }
    // add output vars and meters
    return output.c_str();
}

int apiDataFullyReady() {
    if (EnergyPlus::PluginManagement::fullyReady) {
        return 0;
    }
    return 1;
}


void requestVariable(const char* type, const char* key) {
    // allow specifying a request for an output variable, so that E+ does not have to keep all of them in memory
    // should be called before energyplus is run!
    // note that the variable request array is cleared during clear_state, so if you run multiple E+ runs, these must be requested again each time.
    EnergyPlus::OutputProcessor::APIOutputVariableRequest request;
    request.varName = type;
    request.varKey = key;
    EnergyPlus::OutputProcessor::apiVarRequests.push_back(request);
}

int getVariableHandle(const char* type, const char* key) {
    // Variables are accessed through a single integer ID, but there are multiple internal types: real and integer.
    // I am going to make the integer handle span all both types, by carefully defining the handle.
    // basically, the handles are contiguous, with:
    //  - index 1 being the first real variable handle
    //  - index N being the highest real variable handle
    //  - index N+1 being the first integer variable handle
    //  - index N+M being the highest integer variable handle
    // In this function, it is as simple as looping over both types and continuing to increment
    // the handle carefully.  In the getValue function it is just a matter of checking array sizes.
    std::string const typeUC = EnergyPlus::UtilityRoutines::MakeUPPERCase(type);
    std::string const keyUC = EnergyPlus::UtilityRoutines::MakeUPPERCase(key);
    int handle = -1; // initialize to -1 as a flag
    if (EnergyPlus::OutputProcessor::RVariableTypes.allocated()) {
        handle = 0; // initialize to 0 to get a 1 based Array1D index
        for (int i = 1; i <= EnergyPlus::OutputProcessor::NumOfRVariable; i++) {
            auto &availOutputVar = EnergyPlus::OutputProcessor::RVariableTypes(i);
            handle++;
            if (typeUC == availOutputVar.VarNameOnlyUC && keyUC == availOutputVar.KeyNameOnlyUC) {
                return handle;
            }
        }
    }
    if (EnergyPlus::OutputProcessor::IVariableTypes.allocated()) {
        // now, if real variables *were* searched, we need to pick up the handle where it left off, otherwise initialize it to zero
        if (handle == -1) {
            // real variables were not searched, init to zero
            handle = 0;
        } else {
            // real variables where searched, let it just continue where it left off
        }
        for (int i = 1; i <= EnergyPlus::OutputProcessor::NumOfIVariable; i++) {
            auto &availOutputVar = EnergyPlus::OutputProcessor::IVariableTypes(i);
            handle++;
            if (typeUC == availOutputVar.VarNameOnlyUC && keyUC == availOutputVar.KeyNameOnlyUC) {
                return handle;
            }
        }
    }
    return -1; // return -1 if it wasn't found
}

Real64 getVariableValue(const int handle) {
    // this function works in conjunction with the plan set up in getVariableHandle
    // basically, the handles are contiguous, with:
    //  - index 1 being the first real variable handle
    //  - index N being the highest real variable handle
    //  - index N+1 being the first integer variable handle
    //  - index N+M being the highest integer variable handle
    if (handle > 0 && handle <= EnergyPlus::OutputProcessor::NumOfRVariable) {
        auto &thisOutputVar = EnergyPlus::OutputProcessor::RVariableTypes(handle);
        return *thisOutputVar.VarPtr.Which;
    } else if (handle <= EnergyPlus::OutputProcessor::NumOfRVariable + EnergyPlus::OutputProcessor::NumOfIVariable) {
        int thisHandle = handle - EnergyPlus::OutputProcessor::NumOfRVariable;
        auto &thisOutputVar = EnergyPlus::OutputProcessor::IVariableTypes(thisHandle);
        return (Real64)*thisOutputVar.VarPtr.Which;
    } else {
        if (EnergyPlus::DataGlobals::eplusRunningViaAPI) {
            std::cout << "ERROR: Variable handle out of range in getVariableValue" << std::endl;
            throw std::runtime_error("Variable handle out of range in getVariableValue");
        } else {
            // must be running from python plugin, need to fatal out once the plugin is done
            // throw an error, set the fatal flag, and then return zero
            EnergyPlus::ShowSevereError("Data Exchange API: Index error in getVariableValue; received handle: " + std::to_string(handle));
            EnergyPlus::ShowContinueError("The getVariableValue function will return 0 for now to allow the plugin to finish, then EnergyPlus will abort");
            EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
            return 0;
        }
    }
}


int getMeterHandle(const char* meterName) {
    std::string const meterNameUC = EnergyPlus::UtilityRoutines::MakeUPPERCase(meterName);
    auto i = EnergyPlus::GetMeterIndex(meterNameUC);
    if (i == 0) {
        // inside E+, zero is meaningful, but through the API, I want to use negative one as a signal of a bad lookup
        return -1;
    } else {
        return i;
    }
}

Real64 getMeterValue(int handle) {
    if (handle >= 1 && handle <= (int)EnergyPlus::OutputProcessor::EnergyMeters.size()) {
        return EnergyPlus::GetCurrentMeterValue(handle);
    } else {
        if (EnergyPlus::DataGlobals::eplusRunningViaAPI) {
            std::cout << "ERROR: Meter handle out of range in getMeterValue" << std::endl;
            throw std::runtime_error("Meter handle out of range in getMeterValue");
        } else {
            // must be running from python plugin, need to fatal out once the plugin is done
            // throw an error, set the fatal flag, and then return zero
            EnergyPlus::ShowSevereError("Data Exchange API: Index error in getMeterValue; received handle: " + std::to_string(handle));
            EnergyPlus::ShowContinueError("The getMeterValue function will return 0 for now to allow the plugin to finish, then EnergyPlus will abort");
            EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
            return 0;
        }
    }
}


int getActuatorHandle(const char* componentType, const char* controlType, const char* uniqueKey) {
    int handle;
    handle = 0;
    std::string const typeUC = EnergyPlus::UtilityRoutines::MakeUPPERCase(componentType);
    std::string const keyUC = EnergyPlus::UtilityRoutines::MakeUPPERCase(uniqueKey);
    std::string const controlUC = EnergyPlus::UtilityRoutines::MakeUPPERCase(controlType);
    for (auto const & availActuator : EnergyPlus::DataRuntimeLanguage::EMSActuatorAvailable) { // is this iterating over cushion space?
        handle++;
        std::string const actuatorTypeUC = EnergyPlus::UtilityRoutines::MakeUPPERCase(availActuator.ComponentTypeName);
        std::string const actuatorIDUC = EnergyPlus::UtilityRoutines::MakeUPPERCase(availActuator.UniqueIDName);
        std::string const actuatorControlUC = EnergyPlus::UtilityRoutines::MakeUPPERCase(availActuator.ControlTypeName);
        if (typeUC == actuatorTypeUC && keyUC == actuatorIDUC && controlUC == actuatorControlUC) {
            return handle;
        }
    }
    return -1;
}

void resetActuator(int handle) {
    // resets the actuator so that E+ will use the internally calculated value again
    if (handle >= 1 && handle <= (int)EnergyPlus::DataRuntimeLanguage::EMSActuatorAvailable.size()) {
        auto & theActuator(EnergyPlus::DataRuntimeLanguage::EMSActuatorAvailable(handle));
        *theActuator.Actuated = false;
    } else {
        if (EnergyPlus::DataGlobals::eplusRunningViaAPI) {
            std::cout << "ERROR: Actuator handle out of range in resetActuator" << std::endl;
            throw std::runtime_error("Actuator handle out of range in resetActuator");
        } else {
            // must be running from python plugin, need to fatal out once the plugin is done
            // throw an error, set the fatal flag, and then return
            EnergyPlus::ShowSevereError("Data Exchange API: index error in resetActuator; received handle: " + std::to_string(handle));
            EnergyPlus::ShowContinueError("The resetActuator function will return to allow the plugin to finish, then EnergyPlus will abort");
            EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
        }
    }
}

void setActuatorValue(const int handle, const Real64 value) {
    if (handle >= 1 && handle <= (int)EnergyPlus::DataRuntimeLanguage::EMSActuatorAvailable.size()) {
        auto & theActuator(EnergyPlus::DataRuntimeLanguage::EMSActuatorAvailable(handle));
        if (theActuator.RealValue) {
            *theActuator.RealValue = value;
        } else if (theActuator.IntValue) {
            *theActuator.IntValue = (int)std::lround(value);
        } else {
            // follow protocol from EMS manager, where 1.0 is true, 0.0 is false, and anything else is also false
            *theActuator.LogValue = value > 0.99999 && value < 1.00001; // allow small tolerance while passing between languages and types
        }
        *theActuator.Actuated = true;
    } else {
        if (EnergyPlus::DataGlobals::eplusRunningViaAPI) {
            std::cout << "ERROR: Actuator handle out of range in setActuatorValue" << std::endl;
            throw std::runtime_error("Actuator handle out of range in setActuatorValue");
        } else {
            // must be running from python plugin, need to fatal out once the plugin is done
            // throw an error, set the fatal flag, and then return
            EnergyPlus::ShowSevereError("Data Exchange API: index error in setActuatorValue; received handle: " + std::to_string(handle));
            EnergyPlus::ShowContinueError("The setActuatorValue function will return to allow the plugin to finish, then EnergyPlus will abort");
            EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
        }
    }
}

Real64 getActuatorValue(const int handle) {
    if (handle >= 1 && handle <= (int)EnergyPlus::DataRuntimeLanguage::EMSActuatorAvailable.size()) {
        auto &theActuator(EnergyPlus::DataRuntimeLanguage::EMSActuatorAvailable(handle));
        if (theActuator.RealValue) {
            return *theActuator.RealValue;
        } else if (theActuator.IntValue) {
            return (float)*theActuator.IntValue;
        } else {
            // follow protocol from EMS manager, where 1.0 is true, 0.0 is false, and anything else is also false
            if (*theActuator.LogValue) {
                return 1;
            } else {
                return 0;
            }
        }
    } else {
        if (EnergyPlus::DataGlobals::eplusRunningViaAPI) {
            std::cout << "ERROR: Actuator handle out of range in getActuatorValue" << std::endl;
            throw std::runtime_error("Actuator handle out of range in getActuatorValue");
        } else {
            // must be running from python plugin, need to fatal out once the plugin is done
            // throw an error, set the fatal flag, and then return 0
            EnergyPlus::ShowSevereError("Data Exchange API: index error in getActuatorValue; received handle: " + std::to_string(handle));
            EnergyPlus::ShowContinueError("The getActuatorValue function will return 0 for now to allow the plugin to finish, then EnergyPlus will abort");
            EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
            return 0;
        }
    }
}

int getInternalVariableHandle(const char* type, const char* key) {
    int handle;
    handle = 0;
    std::string const typeUC = EnergyPlus::UtilityRoutines::MakeUPPERCase(type);
    std::string const keyUC = EnergyPlus::UtilityRoutines::MakeUPPERCase(key);
    for (auto const & availVariable : EnergyPlus::DataRuntimeLanguage::EMSInternalVarsAvailable) {
        handle++;
        std::string const variableTypeUC = EnergyPlus::UtilityRoutines::MakeUPPERCase(availVariable.DataTypeName);
        std::string const variableIDUC = EnergyPlus::UtilityRoutines::MakeUPPERCase(availVariable.UniqueIDName);
        if (typeUC == variableTypeUC && keyUC == variableIDUC) {
            return handle;
        }
    }
    return -1;
}

Real64 getInternalVariableValue(int handle) {
    if (handle >= 1 && handle <= (int)EnergyPlus::DataRuntimeLanguage::EMSInternalVarsAvailable.size()) {
        auto thisVar = EnergyPlus::DataRuntimeLanguage::EMSInternalVarsAvailable(handle);
        if (thisVar.PntrVarTypeUsed == EnergyPlus::DataRuntimeLanguage::PntrReal) {
            return *thisVar.RealValue;
        } else if (thisVar.PntrVarTypeUsed == EnergyPlus::DataRuntimeLanguage::PntrInteger) {
            return (Real64)(*thisVar.IntValue);
        } else {
            // Doesn't look like this struct actually has a logical member type, so uh, throw here?
            std::cout << "ERROR: Invalid internal variable type here, developer issue." << std::endl;
            throw std::runtime_error("Invalid internal variable type here, developer issue.");
        }
    } else {
        if (EnergyPlus::DataGlobals::eplusRunningViaAPI) {
            std::cout << "ERROR: Internal variable handle out of range in getInternalVariableValue" << std::endl;
            throw std::runtime_error("Internal variable handle out of range in getInternalVariableValue");
        } else {
            // must be running from python plugin, need to fatal out once the plugin is done
            // throw an error, set the fatal flag, and then return 0
            EnergyPlus::ShowSevereError("Data Exchange API: index error in getInternalVariableValue; received handle: " + std::to_string(handle));
            EnergyPlus::ShowContinueError("The getInternalVariableValue function will return 0 for now to allow the plugin to finish, then EnergyPlus will abort");
            EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
            return 0;
        }
    }
}


int getPluginGlobalVariableHandle(const char* name) {
    return EnergyPlus::PluginManagement::pluginManager->getGlobalVariableHandle(name);
}

Real64 getPluginGlobalVariableValue(int handle) {
    if (handle < 0 || handle > EnergyPlus::PluginManagement::pluginManager->maxGlobalVariableIndex) {
        // need to fatal out once the plugin is done
        // throw an error, set the fatal flag, and then return 0
        EnergyPlus::ShowSevereError("Data Exchange API: Problem -- index error in getPluginGlobalVariableValue; received handle: " + std::to_string(handle));
        EnergyPlus::ShowContinueError("The getPluginGlobalVariableValue function will return 0 for now to allow the plugin to finish, then EnergyPlus will abort");
        EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
        return 0;
    }
    return EnergyPlus::PluginManagement::pluginManager->getGlobalVariableValue(handle);
}

void setPluginGlobalVariableValue(int handle, Real64 value) {
    if (handle < 0 || handle > EnergyPlus::PluginManagement::pluginManager->maxGlobalVariableIndex) {
        // need to fatal out once the plugin is done
        // throw an error, set the fatal flag, and then return
        EnergyPlus::ShowSevereError("Data Exchange API: Problem -- index error in setPluginGlobalVariableValue; received handle: " + std::to_string(handle));
        EnergyPlus::ShowContinueError("The getPluginGlobalVariableValue function will return to allow the plugin to finish, then EnergyPlus will abort");
        EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
    }
    EnergyPlus::PluginManagement::pluginManager->setGlobalVariableValue(handle, value);
}

int getPluginTrendVariableHandle(const char* name) {
    return EnergyPlus::PluginManagement::pluginManager->getTrendVariableHandle(name);
}

Real64 getPluginTrendVariableValue(int handle, int timeIndex) {
    if (handle < 0 || handle > EnergyPlus::PluginManagement::pluginManager->maxTrendVariableIndex) {
        // need to fatal out once the plugin is done
        // throw an error, set the fatal flag, and then return
        EnergyPlus::ShowSevereError("Data Exchange API: Problem -- index error in getPluginTrendVariableValue; received handle: " + std::to_string(handle));
        EnergyPlus::ShowContinueError("The getPluginTrendVariableValue function will return 0 to allow the plugin to finish, then EnergyPlus will abort");
        EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
        return 0;
    }
    if (timeIndex < 1 || timeIndex > ((int)EnergyPlus::PluginManagement::pluginManager->getTrendVariableHistorySize(handle))) {
        // need to fatal out once the plugin is done
        // throw an error, set the fatal flag, and then return
        EnergyPlus::ShowSevereError("Data Exchange API: Problem -- trend history count argument out of range in getPluginTrendVariableValue; received value: " + std::to_string(timeIndex));
        EnergyPlus::ShowContinueError("The getPluginTrendVariableValue function will return 0 to allow the plugin to finish, then EnergyPlus will abort");
        EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
        return 0;
    }
    return EnergyPlus::PluginManagement::pluginManager->getTrendVariableValue(handle, timeIndex);
}

Real64 getPluginTrendVariableAverage(int handle, int count) {
    if (handle < 0 || handle > EnergyPlus::PluginManagement::pluginManager->maxTrendVariableIndex) {
        // need to fatal out once the plugin is done
        // throw an error, set the fatal flag, and then return
        EnergyPlus::ShowSevereError("Data Exchange API: Problem -- index error in getPluginTrendVariableAverage; received handle: " + std::to_string(handle));
        EnergyPlus::ShowContinueError("The getPluginTrendVariableAverage function will return 0 to allow the plugin to finish, then EnergyPlus will abort");
        EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
        return 0;
    }
    if (count < 2 || count > ((int)EnergyPlus::PluginManagement::pluginManager->getTrendVariableHistorySize(handle))) {
        // need to fatal out once the plugin is done
        // throw an error, set the fatal flag, and then return
        EnergyPlus::ShowSevereError("Data Exchange API: Problem -- trend history count argument out of range in getPluginTrendVariableAverage; received value: " + std::to_string(count));
        EnergyPlus::ShowContinueError("The getPluginTrendVariableAverage function will return 0 to allow the plugin to finish, then EnergyPlus will abort");
        EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
        return 0;
    }
    return EnergyPlus::PluginManagement::pluginManager->getTrendVariableAverage(handle, count);
}

Real64 getPluginTrendVariableMin(int handle, int count) {
    if (handle < 0 || handle > EnergyPlus::PluginManagement::pluginManager->maxTrendVariableIndex) {
        // need to fatal out once the plugin is done
        // throw an error, set the fatal flag, and then return
        EnergyPlus::ShowSevereError("Data Exchange API: Problem -- index error in getPluginTrendVariableMin; received handle: " + std::to_string(handle));
        EnergyPlus::ShowContinueError("The getPluginTrendVariableMin function will return 0 to allow the plugin to finish, then EnergyPlus will abort");
        EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
        return 0;
    }
    if (count < 2 || count > ((int)EnergyPlus::PluginManagement::pluginManager->getTrendVariableHistorySize(handle))) {
        // need to fatal out once the plugin is done
        // throw an error, set the fatal flag, and then return
        EnergyPlus::ShowSevereError("Data Exchange API: Problem -- trend history count argument out of range in getPluginTrendVariableMin; received value: " + std::to_string(count));
        EnergyPlus::ShowContinueError("The getPluginTrendVariableMin function will return 0 to allow the plugin to finish, then EnergyPlus will abort");
        EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
        return 0;
    }
    return EnergyPlus::PluginManagement::pluginManager->getTrendVariableMin(handle, count);
}

Real64 getPluginTrendVariableMax(int handle, int count) {
    if (handle < 0 || handle > EnergyPlus::PluginManagement::pluginManager->maxTrendVariableIndex) {
        // need to fatal out once the plugin is done
        // throw an error, set the fatal flag, and then return
        EnergyPlus::ShowSevereError("Data Exchange API: Problem -- index error in getPluginTrendVariableMax; received handle: " + std::to_string(handle));
        EnergyPlus::ShowContinueError("The getPluginTrendVariableMax function will return 0 to allow the plugin to finish, then EnergyPlus will abort");
        EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
        return 0;
    }
    if (count < 2 || count > ((int)EnergyPlus::PluginManagement::pluginManager->getTrendVariableHistorySize(handle))) {
        // need to fatal out once the plugin is done
        // throw an error, set the fatal flag, and then return
        EnergyPlus::ShowSevereError("Data Exchange API: Problem -- trend history count argument out of range in getPluginTrendVariableMax; received value: " + std::to_string(count));
        EnergyPlus::ShowContinueError("The getPluginTrendVariableMax function will return 0 to allow the plugin to finish, then EnergyPlus will abort");
        EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
        return 0;
    }
    return EnergyPlus::PluginManagement::pluginManager->getTrendVariableMax(handle, count);
}

Real64 getPluginTrendVariableSum(int handle, int count) {
    if (handle < 0 || handle > EnergyPlus::PluginManagement::pluginManager->maxTrendVariableIndex) {
        // need to fatal out once the plugin is done
        // throw an error, set the fatal flag, and then return
        EnergyPlus::ShowSevereError("Data Exchange API: Problem -- index error in getPluginTrendVariableSum; received handle: " + std::to_string(handle));
        EnergyPlus::ShowContinueError("The getPluginTrendVariableSum function will return 0 to allow the plugin to finish, then EnergyPlus will abort");
        EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
        return 0;
    }
    if (count < 2 || count > ((int)EnergyPlus::PluginManagement::pluginManager->getTrendVariableHistorySize(handle))) {
        // need to fatal out once the plugin is done
        // throw an error, set the fatal flag, and then return
        EnergyPlus::ShowSevereError("Data Exchange API: Problem -- trend history count argument out of range in getPluginTrendVariableSum; received value: " + std::to_string(count));
        EnergyPlus::ShowContinueError("The getPluginTrendVariableSum function will return 0 to allow the plugin to finish, then EnergyPlus will abort");
        EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
        return 0;
    }
    return EnergyPlus::PluginManagement::pluginManager->getTrendVariableSum(handle, count);
}

Real64 getPluginTrendVariableDirection(int handle, int count) {
    if (handle < 0 || handle > EnergyPlus::PluginManagement::pluginManager->maxTrendVariableIndex) {
        // need to fatal out once the plugin is done
        // throw an error, set the fatal flag, and then return
        EnergyPlus::ShowSevereError("Data Exchange API: Problem -- index error in getPluginTrendVariableDirection; received handle: " + std::to_string(handle));
        EnergyPlus::ShowContinueError("The getPluginTrendVariableDirection function will return 0 to allow the plugin to finish, then EnergyPlus will abort");
        EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
        return 0;
    }
    if (count < 2 || count > ((int)EnergyPlus::PluginManagement::pluginManager->getTrendVariableHistorySize(handle))) {
        // need to fatal out once the plugin is done
        // throw an error, set the fatal flag, and then return
        EnergyPlus::ShowSevereError("Data Exchange API: Problem -- trend history count argument out of range in getPluginTrendVariableDirection; received value: " + std::to_string(count));
        EnergyPlus::ShowContinueError("The getPluginTrendVariableDirection function will return 0 to allow the plugin to finish, then EnergyPlus will abort");
        EnergyPlus::PluginManagement::shouldIssueFatalAfterPluginCompletes = true;
        return 0;
    }
    return EnergyPlus::PluginManagement::pluginManager->getTrendVariableDirection(handle, count);
}


int year() {
    return EnergyPlus::DataEnvironment::Year;
}

int month() {
    return EnergyPlus::DataEnvironment::Month;
}

int dayOfMonth() {
    return EnergyPlus::DataEnvironment::DayOfMonth;
}

int dayOfWeek() {
    return EnergyPlus::DataEnvironment::DayOfWeek;
}

int dayOfYear() {
    return EnergyPlus::DataEnvironment::DayOfYear;
}

int daylightSavingsTimeIndicator() {
    return EnergyPlus::DataEnvironment::DSTIndicator;
}

int hour() {
    return EnergyPlus::DataGlobals::HourOfDay - 1; // no, just stay on 0..23+ DSTadjust ! offset by 1 and daylight savings time
}

Real64 currentTime() {
    if (EnergyPlus::DataHVACGlobals::TimeStepSys < EnergyPlus::DataGlobals::TimeStepZone) {
        // CurrentTime is for end of zone timestep, need to move back to beginning of current zone timestep, then add on system time elapsed already plus current system timestep
        return EnergyPlus::DataGlobals::CurrentTime - EnergyPlus::DataGlobals::TimeStepZone + EnergyPlus::DataHVACGlobals::SysTimeElapsed + EnergyPlus::DataHVACGlobals::TimeStepSys;
    } else {
        return EnergyPlus::DataGlobals::CurrentTime;
    }
}

int minutes() {
    // the -1 is to push us to the right minute, but this should be handled cautiously because if we are inside the HVAC iteration loop,
    // currentTime() returns a floating point fractional hour, so truncation could put this a few seconds from the expected minute.
    Real64 currentTimeVal = currentTime();
    Real64 fractionalHoursIntoTheDay = currentTimeVal - double(EnergyPlus::DataGlobals::HourOfDay - 1);
    Real64 fractionalMinutesIntoTheDay = fractionalHoursIntoTheDay * 60.0;
    return (int)(fractionalMinutesIntoTheDay);
}

int holidayIndex() {
    return EnergyPlus::DataEnvironment::HolidayIndex;
}

int sunIsUp() { // maintain response convention from previous (EMS) implementation
    if (EnergyPlus::DataEnvironment::SunIsUp) {
        return 1;
    } else {
        return 0;
    }
}

int isRaining() {
    if (EnergyPlus::DataEnvironment::IsRain) {
        return 1;
    } else {
        return 0;
    }
}

int warmupFlag() {
    if (EnergyPlus::DataGlobals::WarmupFlag) {
        return 1;
    } else {
        return 0;
    }
}

Real64 systemTimeStep() {
    return EnergyPlus::DataHVACGlobals::TimeStepSys;
}

int currentEnvironmentNum() {
    return EnergyPlus::DataEnvironment::CurEnvirNum;
}

int kindOfSim() {
    return EnergyPlus::DataGlobals::KindOfSim;
}

int getConstructionHandle(const char* constructionName) {
    int handle;
    handle = 0;
    std::string const nameUC = EnergyPlus::UtilityRoutines::MakeUPPERCase(constructionName);
    for (auto const & construct : EnergyPlus::DataHeatBalance::Construct) {
        handle++;
        std::string const thisNameUC = EnergyPlus::UtilityRoutines::MakeUPPERCase(construct.Name);
        if (nameUC == thisNameUC) {
            return handle;
        }
    }
    return 0;
}
