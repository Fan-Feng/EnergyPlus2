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

#ifndef EnergyPlusAPIRuntime_h_INCLUDED
#define EnergyPlusAPIRuntime_h_INCLUDED

#include <EnergyPlus/api/TypeDefs.h>
#include <EnergyPlus/api/EnergyPlusAPI.h>

#ifdef __cplusplus

#include <functional>

ENERGYPLUSLIB_API void callbackBeginNewEnvironment(std::function<void ()> f);
ENERGYPLUSLIB_API void callbackAfterNewEnvironmentWarmupComplete(std::function<void ()> f);
ENERGYPLUSLIB_API void callbackBeginZoneTimeStepBeforeInitHeatBalance(std::function<void ()> f);
ENERGYPLUSLIB_API void callbackBeginZoneTimeStepAfterInitHeatBalance(std::function<void ()> f);
ENERGYPLUSLIB_API void callbackBeginTimeStepBeforePredictor(std::function<void ()> f);
ENERGYPLUSLIB_API void callbackAfterPredictorBeforeHVACManagers(std::function<void ()> f);
ENERGYPLUSLIB_API void callbackAfterPredictorAfterHVACManagers(std::function<void ()> f);
ENERGYPLUSLIB_API void callbackInsideSystemIterationLoop(std::function<void ()> f);
ENERGYPLUSLIB_API void callbackEndOfZoneTimeStepBeforeZoneReporting(std::function<void ()> f);
ENERGYPLUSLIB_API void callbackEndOfZoneTimeStepAfterZoneReporting(std::function<void ()> f);
ENERGYPLUSLIB_API void callbackEndOfSystemTimeStepBeforeHVACReporting(std::function<void ()> f);
ENERGYPLUSLIB_API void callbackEndOfSystemTimeStepAfterHVACReporting(std::function<void ()> f);
ENERGYPLUSLIB_API void callbackEndOfZoneSizing(std::function<void ()> f);
ENERGYPLUSLIB_API void callbackEndOfSystemSizing(std::function<void ()> f);
ENERGYPLUSLIB_API void callbackEndOfAfterComponentGetInput(std::function<void ()> f);
ENERGYPLUSLIB_API void callbackUnitarySystemSizing(std::function<void ()> f);
//ENERGYPLUSLIB_API void callbackUserDefinedComponentModel(std::function<void ()> f);

extern "C" {

#endif  // __cplusplus

/// \file runtime.h
/// \brief This "runtime" API category provides access to link into a running simulation by providing callback functions.
/// \details While an EnergyPlus simulation is running, registered callback functions are "called back" at user-specified
///          points in the simulation.  When a function is called back, it can then leverage other APIs to do program
///          manipulation.  For example, it could leverage the functional API to look up a fluid property value, then
///          use the data exchange API to assign a new control property on a component using an actuator.
///          Once the callback functions have been registered, the client should then call the main `energyplus` function
///          to initiate a new run.  If multiple runs are to be performed in series, the `cClearAllStates` function can be
///          used to "reset" the simulation.  There are also methods to allow issuing warnings through internal EnergyPlus
///          methods and to callback with updates on simulation progress and output messages.
/// \see energyplus
/// \see cClearAllStates

/// \brief Clears the simulation state of EnergyPlus
/// \details After a simulation is complete, if a second is to be run using the same memory space, the simulation state
///          must be cleared, or unexpected errors will occur.
/// \remark This function will also clear any callback functions, so callback functions must be registered again.
ENERGYPLUSLIB_API void cClearAllStates();

/// \brief Runs an EnergyPlus simulation
/// \details This function launches an EnergyPlus simulation using the given arguments.  The arguments are identical to
///          the command line for EnergyPlus(.exe).  Prior to calling this function, it is expected that most workflows
///          will register callback functions to be called back from EnergyPlus at user-specified points in the simulation.
///          If no callbacks are registered, EnergyPlus will proceed as in a regular command line fashion and return
///          when complete.  An example usage of this API endpoint is as follows:
/// \code{.c}
///   void progressCallback(int const i) {
///     printf("Updated progress: %d", i);
///   }
///   int main(int argc, const char * argv[]) {
///     registerProgressCallback(progressCallback);
///     energyplus(argc, argv);
///   }
/// \endcode
/// \remark If this function is called multiple times in the same thread, the client must call `cClearAllStates` in
///         between each run, or unexpected errors will occur.
/// \remark Prior to calling this function, it is also expected that most workflows will want to use the data exchange API
///         function `requestVariable` to pre-request Output:Variables to make them accessible in the data exchange functions.
/// \see cClearAllStates
/// \see requestVariable
ENERGYPLUSLIB_API int energyplus(int argc, const char *argv[]);

/// \brief Asks EnergyPlus to issue a warning message to the error file.
/// \details During an EnergyPlus simulation, if certain conditions arise, it may be useful to alert the user using
///          this function, which will issue a warning note in the standard error file and continue the simulation.
/// \param[in] message The warning message to be issued during the simulation.
/// \remark The usefulness of this function during API applications is questionable, however it is highly valuable
///         during Python Plugin applications.  In these applications, the user is still interfacing with EnergyPlus
///         as a black-box program and will rely on output messages primarily through EnergyPlus output files.  Using
///         this function to issue a message using standard EnergyPlus techniques will make the process familiar.
ENERGYPLUSLIB_API void issueWarning(const char * message); // Issue the warning text to the err file
/// \brief Asks EnergyPlus to issue a severe message to the error file.
/// \details During an EnergyPlus simulation, if certain conditions arise, it may be useful to alert the user using
///          this function, which will issue a severe error note in the standard error file and continue the simulation.
///          Severe errors should lead to program termination, so this should be followed with aborting the simulation in
///          one of these ways:
///          - If running as a Python Plugin, the function should return 1 to inform EnergyPlus to throw a fatal.
///          - If running as an API, the client should abort however is appropriate.  If needed, `cClearAllStates` can be
///            called to reinitialize the program.
///
/// \param[in] message The severe error message to be issued during the simulation.
/// \remark The usefulness of this function during API applications is questionable, however it is highly valuable
///         during Python Plugin applications.  In these applications, the user is still interfacing with EnergyPlus
///         as a black-box program and will rely on output messages primarily through EnergyPlus output files.  Using
///         this function to issue a message using standard EnergyPlus techniques will make the process familiar.
ENERGYPLUSLIB_API void issueSevere(const char * message); // Issue the severe text to the err file
/// \brief Asks EnergyPlus to issue a plain text message to the error file.
/// \details During an EnergyPlus simulation, if certain conditions arise, it may be useful to send information to the
///          user with this function, either to provide standard information, or supplemental information to a previously
///          sent warning or error.
/// \param[in] message The text message to be issued during the simulation.
/// \remark The usefulness of this function during API applications is questionable, however it is highly valuable
///         during Python Plugin applications.  In these applications, the user is still interfacing with EnergyPlus
///         as a black-box program and will rely on output messages primarily through EnergyPlus output files.  Using
///         this function to issue a message using standard EnergyPlus techniques will make the process familiar.
ENERGYPLUSLIB_API void issueText(const char * message); // Issue additional supporting text to the err file


/// \brief Register a callback function to receive updates on simulation progress
/// \details During an EnergyPlus simulation, the progress of the simulation will move from zero to one hundred percent.
///          This function can be used to get updates during the simulation that could be used to better provide
///          feedback to the user, perhaps in the form of a progress bar.
/// \param[in] f The function to be called back during the simulation.  The function should expect one constant integer
///              argument, which will be the percent complete for the current simulation run.
/// \remark This function has limited usefulness in Python Plugin applications, but is highly valuable for being able
///         to report progress during an API workflow simulation.
ENERGYPLUSLIB_API void registerProgressCallback(void (*f)(int const));
/// \brief Register a callback function to receive standard output messages being sent from the simulation.
/// \details During an EnergyPlus simulation, there are a number of output messages being sent to standard output on the
///          command line.  This function can be used to get copies of those messages during the simulation which could
///          be captured and better presented to the user during an API workflow.
/// \param[in] f The function to be called back during the simulation.  The function should expect one constant char *
///              argument, which will be the message being sent to standard output.
/// \remark This function has limited usefulness in Python Plugin applications, but is highly valuable for being able
///         to report messages during an API workflow simulation.
ENERGYPLUSLIB_API void registerStdOutCallback(void (*f)(const char * message));


/// \brief Register a callback function to be called at the beginning of each new environment.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
ENERGYPLUSLIB_API void callbackBeginNewEnvironment(void (*f)());
/// \brief Register a callback function to be called at each new environment once warmup is complete.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
ENERGYPLUSLIB_API void callbackAfterNewEnvironmentWarmupComplete(void (*f)());
/// \brief Register a callback function to be called at the beginning of each zone time step before the heat balance is initialized.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
ENERGYPLUSLIB_API void callbackBeginZoneTimeStepBeforeInitHeatBalance(void (*f)());
/// \brief Register a callback function to be called at the beginning of each zone time step after the heat balance is initialized.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
ENERGYPLUSLIB_API void callbackBeginZoneTimeStepAfterInitHeatBalance(void (*f)());
/// \brief Register a callback function to be called at each zone time step just before the predictor step.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
ENERGYPLUSLIB_API void callbackBeginTimeStepBeforePredictor(void (*f)());
/// \brief Register a callback function to be called at each zone time step before HVAC managers are called.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
ENERGYPLUSLIB_API void callbackAfterPredictorBeforeHVACManagers(void (*f)());
/// \brief Register a callback function to be called at each zone time step after HVAC managers are called.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
ENERGYPLUSLIB_API void callbackAfterPredictorAfterHVACManagers(void (*f)());
/// \brief Register a callback function to be called inside the HVAC system iteration loop.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
ENERGYPLUSLIB_API void callbackInsideSystemIterationLoop(void (*f)());
/// \brief Register a callback function to be called at the end of each zone time step but before zones have reported data.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
ENERGYPLUSLIB_API void callbackEndOfZoneTimeStepBeforeZoneReporting(void (*f)());
/// \brief Register a callback function to be called at the end of each zone time step after zones have reported data.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
ENERGYPLUSLIB_API void callbackEndOfZoneTimeStepAfterZoneReporting(void (*f)());
/// \brief Register a callback function to be called at the end of each HVAC system time step but before HVAC data has been reported.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
ENERGYPLUSLIB_API void callbackEndOfSystemTimeStepBeforeHVACReporting(void (*f)());
/// \brief Register a callback function to be called at the end of each HVAC system time step after HVAC data has been reported.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
ENERGYPLUSLIB_API void callbackEndOfSystemTimeStepAfterHVACReporting(void (*f)());
/// \brief Register a callback function to be called at the end of the zone sizing process.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
ENERGYPLUSLIB_API void callbackEndOfZoneSizing(void (*f)());
/// \brief Register a callback function to be called at the end of the HVAC system sizing process.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
ENERGYPLUSLIB_API void callbackEndOfSystemSizing(void (*f)());
/// \brief Register a callback function to be called after some specific components have processed their input.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
/// \remark Currently this calling point is called after each instance of a DX coil, fan, and unitary system object.
ENERGYPLUSLIB_API void callbackEndOfAfterComponentGetInput(void (*f)());
/// \brief Register a callback function to be called at unitary system sizing.
/// \details During an EnergyPlus simulation, a number of predetermined calling points have been established at which
///          any registered callback functions are "called back".  This API function allows a client to register a function
///          with no arguments to be called at this specific calling point.  From inside this function, the client can
///          leverage other API categories to look up property values or exchange data with the simulation as needed.
/// \param[in] f The function to be called back at this specific calling point in the simulation.  The function expects no arguments.
/// \remark This function is only allowed during API simulations.  For Python Plugin applications, the client will
///         create a custom Python class and override specific functions to be called at equivalent points in the simulation.
ENERGYPLUSLIB_API void callbackUnitarySystemSizing(void (*f)());
// The user defined component model won't actually call out to this API endpoint -- it is coupled with a specific
// plugin instance in the code.
// ENERGYPLUSLIB_API void callbackUserDefinedComponentModel(void (*f)());

#ifdef __cplusplus
}
#endif // __cplusplug

#endif // EnergyPlusAPIRuntime_h_INCLUDED
