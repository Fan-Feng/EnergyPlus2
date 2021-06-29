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

#ifndef ConvectionCoefficients_hh_INCLUDED
#define ConvectionCoefficients_hh_INCLUDED

// ObjexxFCL Headers
#include <ObjexxFCL/Array1D.hh>
#include <ObjexxFCL/Array1S.hh>
#include <ObjexxFCL/Optional.hh>

// EnergyPlus Headers
#include <EnergyPlus/DataGlobals.hh>
#include <EnergyPlus/DataVectorTypes.hh>
#include <EnergyPlus/EnergyPlus.hh>

namespace EnergyPlus {
    class OutputFiles;

namespace ConvectionCoefficients {

    // Using/Aliasing
    using DataVectorTypes::Vector;

    // Data
    // MODULE PARAMETER DEFINITIONS:
    extern Real64 const AdaptiveHcInsideLowLimit;  // W/m2-K
    extern Real64 const AdaptiveHcOutsideLowLimit; // W/m2-K

    extern Real64 const OneThird;   // 1/3 in highest precision
    extern Real64 const OneFourth;  // 1/4 in highest precision
    extern Real64 const OneFifth;   // 1/5 in highest precision
    extern Real64 const OneSixth;   // 1/6 in highest precision
    extern Real64 const FourFifths; // 4/5 in highest precision

    // Coefficients that modify the convection coeff based on surface roughness
    extern Array1D<Real64> const RoughnessMultiplier;

    // parameters for identifying more specific hc model equations, inside face
    extern int const HcInt_UserValue;
    extern int const HcInt_UserSchedule;
    extern int const HcInt_UserCurve;
    extern int const HcInt_ASHRAEVerticalWall;
    extern int const HcInt_WaltonUnstableHorizontalOrTilt;
    extern int const HcInt_WaltonStableHorizontalOrTilt;
    extern int const HcInt_FisherPedersenCeilDiffuserFloor;
    extern int const HcInt_FisherPedersenCeilDiffuserCeiling;
    extern int const HcInt_FisherPedersenCeilDiffuserWalls;
    extern int const HcInt_AlamdariHammondStableHorizontal;
    extern int const HcInt_AlamdariHammondVerticalWall;
    extern int const HcInt_AlamdariHammondUnstableHorizontal;
    extern int const HcInt_KhalifaEq3WallAwayFromHeat;
    extern int const HcInt_KhalifaEq4CeilingAwayFromHeat;
    extern int const HcInt_KhalifaEq5WallNearHeat;
    extern int const HcInt_KhalifaEq6NonHeatedWalls;
    extern int const HcInt_KhalifaEq7Ceiling;
    extern int const HcInt_AwbiHattonHeatedFloor;
    extern int const HcInt_AwbiHattonHeatedWall;
    extern int const HcInt_BeausoleilMorrisonMixedAssistingWall;
    extern int const HcInt_BeausoleilMorrisonMixedOppossingWall;
    extern int const HcInt_BeausoleilMorrisonMixedStableCeiling;
    extern int const HcInt_BeausoleilMorrisonMixedUnstableCeiling;
    extern int const HcInt_BeausoleilMorrisonMixedStableFloor;
    extern int const HcInt_BeausoleilMorrisonMixedUnstableFloor;
    extern int const HcInt_FohannoPolidoriVerticalWall;
    extern int const HcInt_KaradagChilledCeiling;
    extern int const HcInt_ISO15099Windows;
    extern int const HcInt_GoldsteinNovoselacCeilingDiffuserWindow;
    extern int const HcInt_GoldsteinNovoselacCeilingDiffuserWalls;
    extern int const HcInt_GoldsteinNovoselacCeilingDiffuserFloor;

    // parameters for identifying more specific hc model equations, outside face
    extern int const HcExt_None; // none is allowed because Hn and Hf are split
    extern int const HcExt_UserValue;
    extern int const HcExt_UserSchedule;
    extern int const HcExt_UserCurve;
    extern int const HcExt_ASHRAESimpleCombined;
    extern int const HcExt_NaturalASHRAEVerticalWall;
    extern int const HcExt_NaturalWaltonUnstableHorizontalOrTilt;
    extern int const HcExt_NaturalWaltonStableHorizontalOrTilt;
    extern int const HcExt_SparrowWindward;
    extern int const HcExt_SparrowLeeward;
    extern int const HcExt_MoWiTTWindward;
    extern int const HcExt_MoWiTTLeeward;
    extern int const HcExt_DOE2Windward;
    extern int const HcExt_DOE2Leeward;
    extern int const HcExt_NusseltJurges;
    extern int const HcExt_McAdams;
    extern int const HcExt_Mitchell;
    extern int const HcExt_ClearRoof;
    extern int const HcExt_BlockenWindward;
    extern int const HcExt_EmmelVertical;
    extern int const HcExt_EmmelRoof;
    extern int const HcExt_AlamdariHammondVerticalWall;
    extern int const HcExt_FohannoPolidoriVerticalWall;
    extern int const HcExt_ISO15099Windows;
    extern int const HcExt_AlamdariHammondStableHorizontal;
    extern int const HcExt_AlamdariHammondUnstableHorizontal;

    // parameters, by zone, for flow regimes for adaptive convection on inside face
    extern int const InConvFlowRegime_A1; // In-floor heating or in-ceiling cooling
    extern int const InConvFlowRegime_A2; // In-wall heating
    extern int const InConvFlowRegime_A3; // no HVAC system, all bouyancy
    extern int const InConvFlowRegime_B;  // Convective heater in zone
    extern int const InConvFlowRegime_C;  // central mechanical air
    extern int const InConvFlowRegime_D;  // zone mechanical air
    extern int const InConvFlowRegime_E;  // mixed. mechancial air and bouyancy

    // params for reference temperature type
    extern int const RefTempMeanAirTemp;
    extern int const RefTempAdjacentAirTemp;
    extern int const RefTempSupplyAirTemp;

    // params for wind speed type
    extern int const RefWindWeatherFile;
    extern int const RefWindAtZ;
    extern int const RefWindParallComp;
    extern int const RefWindParallCompAtZ;

    extern bool GetUserSuppliedConvectionCoeffs; // Get user input first call for Init

    extern Real64 CubeRootOfOverallBuildingVolume; // building meta data. cube root of the volume of all the zones
    extern Real64 RoofLongAxisOutwardAzimuth;      // roof surfaces meta data. outward normal azimuth for longest roof edge

    struct HcInsideFaceUserCurveStruct
    {
        // Members
        std::string Name; // user's name for object
        int ReferenceTempType;
        int HcFnTempDiffCurveNum;
        int HcFnTempDiffDivHeightCurveNum;
        int HcFnACHCurveNum;
        int HcFnACHDivPerimLengthCurveNum;

        // Default Constructor
        HcInsideFaceUserCurveStruct()
            : ReferenceTempType(0), HcFnTempDiffCurveNum(0), HcFnTempDiffDivHeightCurveNum(0), HcFnACHCurveNum(0), HcFnACHDivPerimLengthCurveNum(0)
        {
        }
    };

    struct HcOutsideFaceUserCurveStruct
    {
        // Members
        std::string Name;
        int ReferenceTempType;
        bool SuppressRainChange;
        int WindSpeedType;
        int HfFnWindSpeedCurveNum;
        int HnFnTempDiffCurveNum;
        int HnFnTempDiffDivHeightCurveNum;

        // Default Constructor
        HcOutsideFaceUserCurveStruct()
            : ReferenceTempType(0), SuppressRainChange(false), WindSpeedType(0), HfFnWindSpeedCurveNum(0), HnFnTempDiffCurveNum(0),
              HnFnTempDiffDivHeightCurveNum(0)
        {
        }
    };

    struct InsideFaceAdaptiveConvAlgoStruct
    {
        // Members
        bool EnteredByUser;
        std::string Name;
        int SimpleBouyVertWallEqNum; // InConvClass_A3_VertWalls
        int SimpleBouyVertWallUserCurveNum;
        int SimpleBouyStableHorizEqNum; // InConvClass_A3_StableHoriz
        int SimpleBouyStableHorizUserCurveNum;
        int SimpleBouyUnstableHorizEqNum; // InConvClass_A3_UnstableHoriz
        int SimpleBouyUnstableHorizUserCurveNum;
        int SimpleBouyStableTiltedEqNum; // InConvClass_A3_StableTilted
        int SimpleBouyStableTiltedUserCurveNum;
        int SimpleBouyUnstableTiltedEqNum; // InConvClass_A3_UnstableTilted
        int SimpleBouyUnstableTiltedUserCurveNum;
        int SimpleBouyWindowsEqNum; // InConvClass_A3_Windows
        int SimpleBouyWindowsUserCurveNum;
        int FloorHeatCeilingCoolVertWallEqNum; // InConvClass_A1_VertWalls
        int FloorHeatCeilingCoolVertWallUserCurveNum;
        int FloorHeatCeilingCoolStableHorizEqNum; // InConvClass_A1_StableHoriz
        int FloorHeatCeilingCoolStableHorizUserCurveNum;
        int FloorHeatCeilingCoolUnstableHorizEqNum; // InConvClass_A1_UntableHoriz
        int FloorHeatCeilingCoolUnstableHorizUserCurveNum;
        int FloorHeatCeilingCoolHeatedFloorEqNum; // InConvClass_A1_HeatedFloor
        int FloorHeatCeilingCoolHeatedFloorUserCurveNum;
        int FloorHeatCeilingCoolChilledCeilingEqNum; // InConvClass_A1_ChilledCeil
        int FloorHeatCeilingCoolChilledCeilingUserCurveNum;
        int FloorHeatCeilingCoolStableTiltedEqNum; // InConvClass_A1_StableTilted
        int FloorHeatCeilingCoolStableTiltedUserCurveNum;
        int FloorHeatCeilingCoolUnstableTiltedEqNum; // InConvClass_A1_UnstableTilted
        int FloorHeatCeilingCoolUnstableTiltedUserCurveNum;
        int FloorHeatCeilingCoolWindowsEqNum; // InConvClass_A1_Windows
        int FloorHeatCeilingCoolWindowsUserCurveNum;
        int WallPanelHeatVertWallEqNum; // InConvClass_A2_VertWallsNonHeated
        int WallPanelHeatVertWallUserCurveNum;
        int WallPanelHeatHeatedWallEqNum; // InConvClass_A2_HeatedVerticalWall
        int WallPanelHeatHeatedWallUserCurveNum;
        int WallPanelHeatStableHorizEqNum; // InConvClass_A2_StableHoriz
        int WallPanelHeatStableHorizUserCurveNum;
        int WallPanelHeatUnstableHorizEqNum; // InConvClass_A2_UnstableHoriz
        int WallPanelHeatUnstableHorizUserCurveNum;
        int WallPanelHeatStableTiltedEqNum; // InConvClass_A2_StableTilted
        int WallPanelHeatStableTiltedUserCurveNum;
        int WallPanelHeatUnstableTiltedEqNum; // InConvClass_A2_UnstableTilted
        int WallPanelHeatUnstableTiltedUserCurveNum;
        int WallPanelHeatWindowsEqNum; // InConvClass_A2_Windows
        int WallPanelHeatWindowsUserCurveNum;
        int ConvectiveHeatVertWallEqNum;
        int ConvectiveHeatVertWallUserCurveNum;
        int ConvectiveHeatVertWallNearHeaterEqNum;
        int ConvectiveHeatVertWallNearHeaterUserCurveNum;
        int ConvectiveHeatStableHorizEqNum;
        int ConvectiveHeatStableHorizUserCurveNum;
        int ConvectiveHeatUnstableHorizEqNum;
        int ConvectiveHeatUnstableHorizUserCurveNum;
        int ConvectiveHeatStableTiltedEqNum;
        int ConvectiveHeatStableTiltedUserCurveNum;
        int ConvectiveHeatUnstableTiltedEqNum;
        int ConvectiveHeatUnstableTiltedUserCurveNum;
        int ConvectiveHeatWindowsEqNum;
        int ConvectiveHeatWindowsUserCurveNum;
        int CentralAirWallEqNum;
        int CentralAirWallUserCurveNum;
        int CentralAirCeilingEqNum;
        int CentralAirCeilingUserCurveNum;
        int CentralAirFloorEqNum;
        int CentralAirFloorUserCurveNum;
        int CentralAirWindowsEqNum;
        int CentralAirWindowsUserCurveNum;
        int ZoneFanCircVertWallEqNum;
        int ZoneFanCircVertWallUserCurveNum;
        int ZoneFanCircStableHorizEqNum;
        int ZoneFanCircStableHorizUserCurveNum;
        int ZoneFanCircUnstableHorizEqNum;
        int ZoneFanCircUnstableHorizUserCurveNum;
        int ZoneFanCircStableTiltedEqNum;
        int ZoneFanCircStableTiltedUserCurveNum;
        int ZoneFanCircUnstableTiltedEqNum;
        int ZoneFanCircUnstableTiltedUserCurveNum;
        int ZoneFanCircWindowsEqNum;
        int ZoneFanCircWindowsUserCurveNum;
        int MixedBouyAssistingFlowWallEqNum;
        int MixedBouyAssistingFlowWallUserCurveNum;
        int MixedBouyOppossingFlowWallEqNum;
        int MixedBouyOppossingFlowWallUserCurveNum;
        int MixedStableFloorEqNum;
        int MixedStableFloorUserCurveNum;
        int MixedUnstableFloorEqNum;
        int MixedUnstableFloorUserCurveNum;
        int MixedStableCeilingEqNum;
        int MixedStableCeilingUserCurveNum;
        int MixedUnstableCeilingEqNum;
        int MixedUnstableCeilingUserCurveNum;
        int MixedWindowsEqNum;
        int MixedWindowsUserCurveNum;

        // Default Constructor
        InsideFaceAdaptiveConvAlgoStruct()
            : EnteredByUser(false), SimpleBouyVertWallEqNum(HcInt_FohannoPolidoriVerticalWall), SimpleBouyVertWallUserCurveNum(0),
              SimpleBouyStableHorizEqNum(HcInt_AlamdariHammondStableHorizontal), SimpleBouyStableHorizUserCurveNum(0),
              SimpleBouyUnstableHorizEqNum(HcInt_AlamdariHammondUnstableHorizontal), SimpleBouyUnstableHorizUserCurveNum(0),
              SimpleBouyStableTiltedEqNum(HcInt_WaltonStableHorizontalOrTilt), SimpleBouyStableTiltedUserCurveNum(0),
              SimpleBouyUnstableTiltedEqNum(HcInt_WaltonUnstableHorizontalOrTilt), SimpleBouyUnstableTiltedUserCurveNum(0),
              SimpleBouyWindowsEqNum(HcInt_ISO15099Windows), SimpleBouyWindowsUserCurveNum(0),
              FloorHeatCeilingCoolVertWallEqNum(HcInt_KhalifaEq3WallAwayFromHeat), FloorHeatCeilingCoolVertWallUserCurveNum(0),
              FloorHeatCeilingCoolStableHorizEqNum(HcInt_AlamdariHammondStableHorizontal), FloorHeatCeilingCoolStableHorizUserCurveNum(0),
              FloorHeatCeilingCoolUnstableHorizEqNum(HcInt_KhalifaEq4CeilingAwayFromHeat), FloorHeatCeilingCoolUnstableHorizUserCurveNum(0),
              FloorHeatCeilingCoolHeatedFloorEqNum(HcInt_AwbiHattonHeatedFloor), FloorHeatCeilingCoolHeatedFloorUserCurveNum(0),
              FloorHeatCeilingCoolChilledCeilingEqNum(HcInt_KaradagChilledCeiling), FloorHeatCeilingCoolChilledCeilingUserCurveNum(0),
              FloorHeatCeilingCoolStableTiltedEqNum(HcInt_WaltonStableHorizontalOrTilt), FloorHeatCeilingCoolStableTiltedUserCurveNum(0),
              FloorHeatCeilingCoolUnstableTiltedEqNum(HcInt_WaltonUnstableHorizontalOrTilt), FloorHeatCeilingCoolUnstableTiltedUserCurveNum(0),
              FloorHeatCeilingCoolWindowsEqNum(HcInt_ISO15099Windows), FloorHeatCeilingCoolWindowsUserCurveNum(0),
              WallPanelHeatVertWallEqNum(HcInt_KhalifaEq6NonHeatedWalls), WallPanelHeatVertWallUserCurveNum(0),
              WallPanelHeatHeatedWallEqNum(HcInt_AwbiHattonHeatedWall), WallPanelHeatHeatedWallUserCurveNum(0),
              WallPanelHeatStableHorizEqNum(HcInt_AlamdariHammondStableHorizontal), WallPanelHeatStableHorizUserCurveNum(0),
              WallPanelHeatUnstableHorizEqNum(HcInt_KhalifaEq7Ceiling), WallPanelHeatUnstableHorizUserCurveNum(0),
              WallPanelHeatStableTiltedEqNum(HcInt_WaltonStableHorizontalOrTilt), WallPanelHeatStableTiltedUserCurveNum(0),
              WallPanelHeatUnstableTiltedEqNum(HcInt_WaltonUnstableHorizontalOrTilt), WallPanelHeatUnstableTiltedUserCurveNum(0),
              WallPanelHeatWindowsEqNum(HcInt_ISO15099Windows), WallPanelHeatWindowsUserCurveNum(0),
              ConvectiveHeatVertWallEqNum(HcInt_FohannoPolidoriVerticalWall), ConvectiveHeatVertWallUserCurveNum(0),
              ConvectiveHeatVertWallNearHeaterEqNum(HcInt_KhalifaEq5WallNearHeat), ConvectiveHeatVertWallNearHeaterUserCurveNum(0),
              ConvectiveHeatStableHorizEqNum(HcInt_AlamdariHammondStableHorizontal), ConvectiveHeatStableHorizUserCurveNum(0),
              ConvectiveHeatUnstableHorizEqNum(HcInt_KhalifaEq7Ceiling), ConvectiveHeatUnstableHorizUserCurveNum(0),
              ConvectiveHeatStableTiltedEqNum(HcInt_WaltonStableHorizontalOrTilt), ConvectiveHeatStableTiltedUserCurveNum(0),
              ConvectiveHeatUnstableTiltedEqNum(HcInt_WaltonUnstableHorizontalOrTilt), ConvectiveHeatUnstableTiltedUserCurveNum(0),
              ConvectiveHeatWindowsEqNum(HcInt_ISO15099Windows), ConvectiveHeatWindowsUserCurveNum(0),
              CentralAirWallEqNum(HcInt_GoldsteinNovoselacCeilingDiffuserWalls), CentralAirWallUserCurveNum(0),
              CentralAirCeilingEqNum(HcInt_FisherPedersenCeilDiffuserCeiling), CentralAirCeilingUserCurveNum(0),
              CentralAirFloorEqNum(HcInt_GoldsteinNovoselacCeilingDiffuserFloor), CentralAirFloorUserCurveNum(0),
              CentralAirWindowsEqNum(HcInt_GoldsteinNovoselacCeilingDiffuserWindow), CentralAirWindowsUserCurveNum(0),
              ZoneFanCircVertWallEqNum(HcInt_KhalifaEq3WallAwayFromHeat), ZoneFanCircVertWallUserCurveNum(0),
              ZoneFanCircStableHorizEqNum(HcInt_AlamdariHammondStableHorizontal), ZoneFanCircStableHorizUserCurveNum(0),
              ZoneFanCircUnstableHorizEqNum(HcInt_KhalifaEq4CeilingAwayFromHeat), ZoneFanCircUnstableHorizUserCurveNum(0),
              ZoneFanCircStableTiltedEqNum(HcInt_WaltonStableHorizontalOrTilt), ZoneFanCircStableTiltedUserCurveNum(0),
              ZoneFanCircUnstableTiltedEqNum(HcInt_WaltonUnstableHorizontalOrTilt), ZoneFanCircUnstableTiltedUserCurveNum(0),
              ZoneFanCircWindowsEqNum(HcInt_ISO15099Windows), ZoneFanCircWindowsUserCurveNum(0),
              MixedBouyAssistingFlowWallEqNum(HcInt_BeausoleilMorrisonMixedAssistingWall), MixedBouyAssistingFlowWallUserCurveNum(0),
              MixedBouyOppossingFlowWallEqNum(HcInt_BeausoleilMorrisonMixedOppossingWall), MixedBouyOppossingFlowWallUserCurveNum(0),
              MixedStableFloorEqNum(HcInt_BeausoleilMorrisonMixedStableFloor), MixedStableFloorUserCurveNum(0),
              MixedUnstableFloorEqNum(HcInt_BeausoleilMorrisonMixedUnstableFloor), MixedUnstableFloorUserCurveNum(0),
              MixedStableCeilingEqNum(HcInt_BeausoleilMorrisonMixedStableCeiling), MixedStableCeilingUserCurveNum(0),
              MixedUnstableCeilingEqNum(HcInt_BeausoleilMorrisonMixedUnstableCeiling), MixedUnstableCeilingUserCurveNum(0),
              MixedWindowsEqNum(HcInt_GoldsteinNovoselacCeilingDiffuserWindow), MixedWindowsUserCurveNum(0)
        {
        }
    };

    struct OutsideFaceAdpativeConvAlgoStruct
    {
        // Members
        bool EnteredByUser;
        std::string Name;
        bool SuppressRainChange;
        int HWindWallWindwardEqNum;
        int HWindWallWindwardUserCurveNum;
        int HWindWallLeewardEqNum;
        int HWindWallLeewardUserCurveNum;
        int HWindHorizRoofEqNum;
        int HWindHorizRoofUserCurveNum;
        int HNatVertWallEqNum;
        int HNatVertWallUserCurveNum;
        int HNatStableHorizEqNum;
        int HNatStableHorizUserCurveNum;
        int HNatUnstableHorizEqNum;
        int HNatUstableHorizUserCurveNum;

        // Default Constructor
        OutsideFaceAdpativeConvAlgoStruct()
            : EnteredByUser(false), SuppressRainChange(false), HWindWallWindwardEqNum(HcExt_BlockenWindward), HWindWallWindwardUserCurveNum(0),
              HWindWallLeewardEqNum(HcExt_EmmelVertical), HWindWallLeewardUserCurveNum(0), HWindHorizRoofEqNum(HcExt_ClearRoof),
              HWindHorizRoofUserCurveNum(0), HNatVertWallEqNum(HcExt_NaturalASHRAEVerticalWall), HNatVertWallUserCurveNum(0),
              HNatStableHorizEqNum(HcExt_NaturalWaltonStableHorizontalOrTilt), HNatStableHorizUserCurveNum(0),
              HNatUnstableHorizEqNum(HcExt_NaturalWaltonUnstableHorizontalOrTilt), HNatUstableHorizUserCurveNum(0)
        {
        }
    };

    struct BoundingBoxVertStruct
    {
        // Members
        int SurfNum;
        int VertNum;
        Vector Vertex;

        // Default Constructor
        BoundingBoxVertStruct() : SurfNum(0), VertNum(0), Vertex(0.0, 0.0, 0.0)
        {
        }
    };

    struct RoofGeoCharactisticsStruct
    {
        // Members
        BoundingBoxVertStruct XdYdZd; // 1 low x, low y, low z
        BoundingBoxVertStruct XdYdZu; // 2 low x, low y, hi z
        BoundingBoxVertStruct XdYuZd; // 3 low x, hi y, low z
        BoundingBoxVertStruct XdYuZu; // 4 low x, hi y, hi z
        BoundingBoxVertStruct XuYdZd; // 5 hi x, low y, low z
        BoundingBoxVertStruct XuYuZd; // 6 hi x, hi y, low z
        BoundingBoxVertStruct XuYdZu; // 7 hi x, low y, hi z
        BoundingBoxVertStruct XuYuZu; // 8 hi x, hi y, hi z
        Array1D<Vector> BoundSurf;    // long edge of roof group bounding surface
        Real64 Area;
        Real64 Perimeter;
        Real64 Height;

        // Default Constructor
        RoofGeoCharactisticsStruct() : BoundSurf(4), Area(0.0), Perimeter(0.0), Height(0.0)
        {
        }
    };

    // Object Data
    extern InsideFaceAdaptiveConvAlgoStruct InsideFaceAdaptiveConvectionAlgo; // stores rules for Hc model equations
    extern OutsideFaceAdpativeConvAlgoStruct OutsideFaceAdaptiveConvectionAlgo;
    extern Array1D<HcInsideFaceUserCurveStruct> HcInsideUserCurve;
    extern Array1D<HcOutsideFaceUserCurveStruct> HcOutsideUserCurve;
    extern RoofGeoCharactisticsStruct RoofGeo;

    // Functions

    void InitInteriorConvectionCoeffs(const Array1D<Real64> &SurfaceTemperatures, // Temperature of surfaces for evaluation of HcIn
                                      Optional_int_const ZoneToResimulate = _    // if passed in, then only calculate surfaces that have this zone
    );

    void InitExteriorConvectionCoeff(int SurfNum,      // Surface number (in Surface derived type)
                                     Real64 HMovInsul, // Equivalent convection coefficient of movable insulation
                                     int Roughness,    // Roughness index (1-6), see DataHeatBalance parameters
                                     Real64 AbsExt,    // Exterior thermal absorptance
                                     Real64 TempExt,   // Exterior surface temperature (C)
                                     Real64 &HExt,           // Convection coefficient to exterior air
                                     Real64 &HSky,           // "Convection" coefficient to sky temperature
                                     Real64 &HGround,        // "Convection" coefficient to ground temperature
                                     Real64 &HAir            // Radiation to Air Component
    );

    Real64 CalcHfExteriorSparrow(Real64 SurfWindSpeed, // Local wind speed at height of the heat transfer surface (m/s)
                                 Real64 GrossArea,     // Gross surface area {m2}
                                 Real64 Perimeter,     // Surface perimeter length {m}
                                 Real64 CosTilt,       // Cosine of the Surface Tilt Angle
                                 Real64 Azimuth,       // Facing angle (degrees) of the surface outward normal
                                 int Roughness,        // Surface roughness index (6=very smooth, 5=smooth, 4=medium smooth,
                                 Real64 WindDirection  // Wind (compass) direction (degrees)
    );

    bool Windward(Real64 CosTilt,      // Cosine of the surface tilt angle
                  Real64 Azimuth,      // or Facing, Direction the surface outward normal faces (degrees)
                  Real64 WindDirection // Wind direction measured clockwise from geographhic North
    );

    void GetUserConvectionCoefficients();

    void ApplyConvectionValue(std::string const &SurfaceTypes, std::string const &ConvectionType, int Value);

    Real64 CalcASHRAESimpExtConvectCoeff(int Roughness,       // Integer index for roughness, relates to parameter array indices
                                         Real64 SurfWindSpeed // Current wind speed, m/s
    );

    Real64 CalcASHRAESimpleIntConvCoeff(Real64 Tsurf, Real64 Tamb, Real64 cosTilt);

    void CalcASHRAESimpleIntConvCoeff(int SurfNum,                  // surface number for which coefficients are being calculated
                                      Real64 SurfaceTemperature,    // Temperature of surface for evaluation of HcIn
                                      Real64 ZoneMeanAirTemperature // Mean Air Temperature of Zone
    );

    Real64 CalcASHRAETARPNatural(Real64 Tsurf, Real64 Tamb, Real64 cosTilt);

    void CalcASHRAEDetailedIntConvCoeff(int SurfNum,                  // surface number for which coefficients are being calculated
                                        Real64 SurfaceTemperature,    // Temperature of surface for evaluation of HcIn
                                        Real64 ZoneMeanAirTemperature // Mean Air Temperature of Zone
    );

    void CalcDetailedHcInForDVModel(int SurfNum,                          // surface number for which coefficients are being calculated
                                    const Array1D<Real64> &SurfaceTemperatures, // Temperature of surfaces for evaluation of HcIn
                                    Array1D<Real64> &HcIn,                      // Interior Convection Coeff Array
                                    Optional<Array1S<Real64> const> Vhc = _     // Velocity array for forced convection coeff calculation
    );

    Real64 CalcZoneSupplyAirTemp(int ZoneNum);

    Real64 CalcZoneSystemVolFlowRate(int ZoneNum);

    Real64 CalcZoneSystemACH(int ZoneNum);

    Real64 CalcCeilingDiffuserACH(int ZoneNum);

    Real64 CalcCeilingDiffuserIntConvCoeff(Real64 ACH,  // [1/hr] air system air change rate
                                           Real64 Tsurf,
                                           Real64 Tair,
                                           Real64 cosTilt,
                                           Real64 humRat,
                                           Real64 height,
                                           bool isWindow=false);

    void CalcCeilingDiffuserIntConvCoeff(int ZoneNum, const Array1D<Real64> &SurfaceTemperatures); // zone number for which coefficients are being calculated

    // CalcCeilingDiffuserInletCorr should replace CalcCeilingDiffuser (above), if ZoneTempPredictorCorrector can
    // ever be made to work correctly with the inlet air temperature.

    void CalcCeilingDiffuserInletCorr(int ZoneNum,                        // Zone number
                                      Array1D<Real64> &SurfaceTemperatures // For CalcASHRAEDetailed, if called
    );

    void CalcTrombeWallIntConvCoeff(int ZoneNum,                        // Zone number for which coefficients are being calculated
                                    const Array1D<Real64> &SurfaceTemperatures // Temperature of surfaces for evaluation of HcIn
    );

    void CalcNusselt(int SurfNum, // Surface number
                     Real64 asp,  // Aspect ratio: window height to gap width
                     Real64 tso,  // Temperature of gap surface closest to outside (K)
                     Real64 tsi,  // Temperature of gap surface closest to zone (K)
                     Real64 gr,   // Gap gas Grashof number
                     Real64 pr,   // Gap gas Prandtl number
                     Real64 &gnu        // Gap gas Nusselt number
    );

    Real64 SetExtConvectionCoeff(int SurfNum); // Surface Number

    Real64 SetIntConvectionCoeff(int SurfNum); // Surface Number

    Real64 CalcISO15099WindowIntConvCoeff(Real64 SurfaceTemperature, // Temperature of surface for evaluation of HcIn
                                          Real64 AirTemperature,     // Mean Air Temperature of Zone (or adjacent air temperature)
                                          Real64 AirHumRat,          // air humidity ratio
                                          Real64 Height,             // window cavity height [m]
                                          Real64 TiltDeg,                  // glazing tilt in degrees
                                          Real64 sineTilt            // sine of glazing tilt
    );

    void CalcISO15099WindowIntConvCoeff(int SurfNum,               // surface number for which coefficients are being calculated
                                        Real64 SurfaceTemperature, // Temperature of surface for evaluation of HcIn
                                        Real64 AirTemperature      // Mean Air Temperature of Zone (or adjacent air temperature)
    );

    void SetupAdaptiveConvectionStaticMetaData(EnergyPlus::OutputFiles &outputFiles);

    void SetupAdaptiveConvectionRadiantSurfaceData();

    void ManageInsideAdaptiveConvectionAlgo(int SurfNum); // surface number for which coefficients are being calculated

    void ManageOutsideAdaptiveConvectionAlgo(int SurfNum, // surface number for which coefficients are being calculated
                                             Real64 &Hc         // result for Hc Outside face, becomes HExt.
    );

    void EvaluateIntHcModels(int SurfNum,
                             int ConvModelEquationNum,
                             Real64 &Hc // calculated Hc value
    );

    void EvaluateExtHcModels(int SurfNum, int NaturalConvModelEqNum, int ForcedConvModelEqNum, Real64 &Hc);

    void DynamicExtConvSurfaceClassification(int SurfNum); // surface number

    void MapExtConvClassificationToHcModels(int SurfNum); // surface number

    void DynamicIntConvSurfaceClassification(int SurfNum); // surface number

    void MapIntConvClassificationToHcModels(int SurfNum); // surface pointer index

    void CalcUserDefinedInsideHcModel(int SurfNum, int UserCurveNum, Real64 &Hc);

    void CalcUserDefinedOutsideHcModel(int SurfNum, int UserCurveNum, Real64 &H);

    //** Begin catalog of Hc equation functions. **** !*************************************************

    Real64 CalcASHRAEVerticalWall(Real64 DeltaTemp); // [C] temperature difference between surface and air

    Real64 CalcWaltonUnstableHorizontalOrTilt(Real64 DeltaTemp, // [C] temperature difference between surface and air
                                              Real64 CosineTilt // Cosine of tilt angle
    );

    Real64 CalcWaltonStableHorizontalOrTilt(Real64 DeltaTemp, // [C] temperature difference between surface and air
                                            Real64 CosineTilt // Cosine of tilt angle
    );

    Real64 CalcFisherPedersenCeilDiffuserFloor(Real64 ACH,  // [1/hr] air system air change rate
                                               Real64 Tsurf,
                                               Real64 Tair,
                                               Real64 cosTilt,
                                               Real64 humRat,
                                               Real64 height,
                                               bool isWindow=false);

    Real64 CalcFisherPedersenCeilDiffuserCeiling(Real64 ACH,  // [1/hr] air system air change rate
                                                 Real64 Tsurf,
                                                 Real64 Tair,
                                                 Real64 cosTilt,
                                                 Real64 humRat,
                                                 Real64 height,
                                                 bool isWindow=false);

    Real64 CalcFisherPedersenCeilDiffuserWalls(Real64 ACH,  // [1/hr] air system air change rate
                                               Real64 Tsurf,
                                               Real64 Tair,
                                               Real64 cosTilt,
                                               Real64 humRat,
                                               Real64 height,
                                               bool isWindow=false);

    Real64 CalcFisherPedersenCeilDiffuserNatConv(Real64 Hforced,
                                                 Real64 ACH,
                                                 Real64 Tsurf,
                                                 Real64 Tair,
                                                 Real64 cosTilt,
                                                 Real64 humRat,
                                                 Real64 height,
                                                 bool isWindow);

    Real64 CalcAlamdariHammondUnstableHorizontal(Real64 DeltaTemp,         // [C] temperature difference between surface and air
                                                   Real64 HydraulicDiameter  // [m] characteristic size, = (4 * area) / perimeter
    );

    Real64 CalcAlamdariHammondUnstableHorizontal(Real64 DeltaTemp,         // [C] temperature difference between surface and air
                                                 Real64 HydraulicDiameter, // [m] characteristic size, = (4 * area) / perimeter
                                                 int SurfNum               // for messages
    );

    Real64 CalcAlamdariHammondStableHorizontal(Real64 DeltaTemp,         // [C] temperature difference between surface and air
                                               Real64 HydraulicDiameter    // [m] characteristic size, = (4 * area) / perimeter
    );

    Real64 CalcAlamdariHammondStableHorizontal(Real64 DeltaTemp,         // [C] temperature difference between surface and air
                                               Real64 HydraulicDiameter, // [m] characteristic size, = (4 * area) / perimeter
                                               int SurfNum               // for messages
    );

    Real64 CalcAlamdariHammondVerticalWall(Real64 DeltaTemp, // [C] temperature difference between surface and air
                                             Real64 Height     // [m] characteristic size, = zone height
    );

    Real64 CalcAlamdariHammondVerticalWall(Real64 DeltaTemp, // [C] temperature difference between surface and air
                                           Real64 Height,    // [m] characteristic size, = zone height
                                           int SurfNum       // for messages
    );

    Real64 CalcKhalifaEq3WallAwayFromHeat(Real64 DeltaTemp); // [C] temperature difference between surface and air

    Real64 CalcKhalifaEq4CeilingAwayFromHeat(Real64 DeltaTemp); // [C] temperature difference between surface and air

    Real64 CalcKhalifaEq5WallsNearHeat(Real64 DeltaTemp); // [C] temperature difference between surface and air

    Real64 CalcKhalifaEq6NonHeatedWalls(Real64 DeltaTemp); // [C] temperature difference between surface and air

    Real64 CalcKhalifaEq7Ceiling(Real64 DeltaTemp); // [C] temperature difference between surface and air

    Real64 CalcAwbiHattonHeatedFloor(Real64 DeltaTemp,        // [C] temperature difference between surface and air
                                     Real64 HydraulicDiameter // [m] characteristic size, = (4 * area) / perimeter
    );

    Real64 CalcAwbiHattonHeatedWall(Real64 DeltaTemp,        // [C] temperature difference between surface and air
                                    Real64 HydraulicDiameter // [m] characteristic size, = (4 * area) / perimeter
    );

    Real64 CalcBeausoleilMorrisonMixedAssistedWall(Real64 const &DeltaTemp,     // [C] temperature difference between surface and air
                                                   Real64 const &Height,        // [m] characteristic size
                                                   Real64 const &SurfTemp,      // [C] surface temperature
                                                   Real64 const &SupplyAirTemp, // [C] temperature of supply air into zone
                                                   Real64 const &AirChangeRate  // [ACH] [1/hour] supply air ACH for zone
    );

    Real64 CalcBeausoleilMorrisonMixedAssistedWall(Real64 const &DeltaTemp,     // [C] temperature difference between surface and air
                                                   Real64 const &Height,        // [m] characteristic size
                                                   Real64 const &SurfTemp,      // [C] surface temperature
                                                   int ZoneNum           // index of zone for messaging
    );

    Real64 CalcBeausoleilMorrisonMixedOpposingWall(Real64 const &DeltaTemp,     // [C] temperature difference between surface and air
                                                   Real64 const &Height,        // [m] characteristic size
                                                   Real64 const &SurfTemp,      // [C] surface temperature
                                                   Real64 const &SupplyAirTemp, // [C] temperature of supply air into zone
                                                   Real64 const &AirChangeRate // [ACH] [1/hour] supply air ACH for zone
    );

    Real64 CalcBeausoleilMorrisonMixedOpposingWall(Real64 const &DeltaTemp,     // [C] temperature difference between surface and air
                                                   Real64 const &Height,        // [m] characteristic size
                                                   Real64 const &SurfTemp,      // [C] surface temperature
                                                   int ZoneNum           // index of zone for messaging
    );

    Real64 CalcBeausoleilMorrisonMixedStableFloor(Real64 const &DeltaTemp,         // [C] temperature difference between surface and air
                                                  Real64 const &HydraulicDiameter, // [m] characteristic size, = (4 * area) / perimeter
                                                  Real64 const &SurfTemp,          // [C] surface temperature
                                                  Real64 const &SupplyAirTemp,     // [C] temperature of supply air into zone
                                                  Real64 const &AirChangeRate      // [ACH] [1/hour] supply air ACH for zone
    );

    Real64 CalcBeausoleilMorrisonMixedStableFloor(Real64 const &DeltaTemp,         // [C] temperature difference between surface and air
                                                  Real64 const &HydraulicDiameter, // [m] characteristic size, = (4 * area) / perimeter
                                                  Real64 const &SurfTemp,          // [C] surface temperature
                                                  int ZoneNum               // index of zone for messaging
    );

    Real64 CalcBeausoleilMorrisonMixedUnstableFloor(Real64 const &DeltaTemp,         // [C] temperature difference between surface and air
                                                    Real64 const &HydraulicDiameter, // [m] characteristic size, = (4 * area) / perimeter
                                                    Real64 const &SurfTemp,          // [C] surface temperature
                                                    Real64 const &SupplyAirTemp,     // [C] temperature of supply air into zone
                                                    Real64 const &AirChangeRate      // [ACH] [1/hour] supply air ACH for zone
    );

    Real64 CalcBeausoleilMorrisonMixedUnstableFloor(Real64 const &DeltaTemp,         // [C] temperature difference between surface and air
                                                    Real64 const &HydraulicDiameter, // [m] characteristic size, = (4 * area) / perimeter
                                                    Real64 const &SurfTemp,          // [C] surface temperature
                                                    int ZoneNum               // index of zone for messaging
    );

    Real64 CalcBeausoleilMorrisonMixedStableCeiling(Real64 const &DeltaTemp,         // [C] temperature difference between surface and air
                                                    Real64 const &HydraulicDiameter, // [m] characteristic size, = (4 * area) / perimeter
                                                    Real64 const &SurfTemp,          // [C] surface temperature
                                                    Real64 const &SupplyAirTemp,     // [C] temperature of supply air into zone
                                                    Real64 const &AirChangeRate      // [ACH] [1/hour] supply air ACH for zone
    );

    Real64 CalcBeausoleilMorrisonMixedStableCeiling(Real64 const &DeltaTemp,         // [C] temperature difference between surface and air
                                                    Real64 const &HydraulicDiameter, // [m] characteristic size, = (4 * area) / perimeter
                                                    Real64 const &SurfTemp,          // [C] surface temperature
                                                    int ZoneNum               // index of zone for messaging
    );

    Real64 CalcBeausoleilMorrisonMixedUnstableCeiling(Real64 const &DeltaTemp,         // [C] temperature difference between surface and air
                                                      Real64 const &HydraulicDiameter, // [m] characteristic size, = (4 * area) / perimeter
                                                      Real64 const &SurfTemp,          // [C] surface temperature
                                                      Real64 const &SupplyAirTemp,     // [C] temperature of supply air into zone
                                                      Real64 const &AirChangeRate     // [ACH] [1/hour] supply air ACH for zone
    );

    Real64 CalcBeausoleilMorrisonMixedUnstableCeiling(Real64 const &DeltaTemp,         // [C] temperature difference between surface and air
                                                      Real64 const &HydraulicDiameter, // [m] characteristic size, = (4 * area) / perimeter
                                                      Real64 const &SurfTemp,          // [C] surface temperature
                                                      int ZoneNum               // index of zone for messaging
    );

    Real64 CalcFohannoPolidoriVerticalWall(Real64 DeltaTemp, // [C] temperature difference between surface and air
                                           Real64 Height,    // [m] characteristic size, height of zone
                                           Real64 SurfTemp,  // [C] surface temperature
                                           Real64 QdotConv   // [W/m2] heat flux rate for rayleigh #
    );

    Real64 CalcFohannoPolidoriVerticalWall(Real64 DeltaTemp, // [C] temperature difference between surface and air
                                           Real64 Height,    // [m] characteristic size, height of zone
                                           Real64 SurfTemp,  // [C] surface temperature
                                           Real64 QdotConv,  // [W/m2] heat flux rate for rayleigh #
                                           int SurfNum       // for messages
    );

    Real64 CalcKaradagChilledCeiling(Real64 DeltaTemp); // [C] temperature difference between surface and air

    Real64 CalcGoldsteinNovoselacCeilingDiffuserWindow(Real64 AirSystemFlowRate,  // [m3/s] air system flow rate
                                                       Real64 ZoneExtPerimLength, // [m] length of zone perimeter with exterior walls
                                                       Real64 WindWallRatio,      // [ ] fraction of window area to wall area for zone
                                                       int WindowLocationType     // index for location types
    );

    Real64 CalcGoldsteinNovoselacCeilingDiffuserWindow(Real64 ZoneExtPerimLength, // [m] length of zone perimeter with exterior walls
                                                       Real64 WindWallRatio,      // [ ] fraction of window area to wall area for zone
                                                       int WindowLocationType,    // index for location types
                                                       int ZoneNum                // for messages
    );

    Real64 CalcGoldsteinNovoselacCeilingDiffuserWall(Real64 AirSystemFlowRate,  // [m3/s] air system flow rate
                                                     Real64 ZoneExtPerimLength, // [m] length of zone perimeter with exterior walls
                                                     int WindowLocationType    // index for location types
    );

    Real64 CalcGoldsteinNovoselacCeilingDiffuserWall(Real64 ZoneExtPerimLength, // [m] length of zone perimeter with exterior walls
                                                     int WindowLocationType,    // index for location types
                                                     int ZoneNum                // for messages
    );

    Real64 CalcGoldsteinNovoselacCeilingDiffuserFloor(Real64 AirSystemFlowRate,  // [m3/s] air system flow rate
                                                      Real64 ZoneExtPerimLength // [m] length of zone perimeter with exterior walls
    );

    Real64 CalcGoldsteinNovoselacCeilingDiffuserFloor(Real64 ZoneExtPerimLength, // [m] length of zone perimeter with exterior walls
                                                      int ZoneNum                // for messages
    );

    Real64 CalcSparrowWindward(int RoughnessIndex, Real64 FacePerimeter, Real64 FaceArea, Real64 WindAtZ);

    Real64 CalcSparrowLeeward(int RoughnessIndex, Real64 FacePerimeter, Real64 FaceArea, Real64 WindAtZ);

    Real64 CalcSparrowWindward(int RoughnessIndex, Real64 FacePerimeter, Real64 FaceArea, Real64 WindAtZ, int SurfNum);

    Real64 CalcSparrowLeeward(int RoughnessIndex, Real64 FacePerimeter, Real64 FaceArea, Real64 WindAtZ, int SurfNum);

    Real64 CalcMoWITTNatural(Real64 DeltaTemp);

    Real64 CalcMoWITTForcedWindward(Real64 WindAtZ);

    Real64 CalcMoWITTForcedLeeward(Real64 WindAtZ);

    Real64 CalcMoWITTWindward(Real64 DeltaTemp, Real64 WindAtZ);

    Real64 CalcMoWITTLeeward(Real64 DeltaTemp, Real64 WindAtZ);

    Real64 CalcDOE2Forced(Real64 SurfaceTemp, Real64 AirTemp, Real64 CosineTilt, Real64 HfSmooth, int RoughnessIndex);

    Real64 CalcDOE2Windward(Real64 SurfaceTemp, Real64 AirTemp, Real64 CosineTilt, Real64 WindAtZ, int RoughnessIndex);

    Real64 CalcDOE2Leeward(Real64 SurfaceTemp, Real64 AirTemp, Real64 CosineTilt, Real64 WindAtZ, int RoughnessIndex);

    Real64 CalcNusseltJurges(Real64 WindAtZ);

    Real64 CalcMcAdams(Real64 WindAtZ);

    Real64 CalcMitchell(Real64 WindAtZ, Real64 LengthScale);

    Real64 CalcMitchell(Real64 WindAtZ, Real64 LengthScale, int SurfNum);

    Real64 CalcBlockenWindward(Real64 WindAt10m,
                               Real64 WindDir,    // Wind direction measured clockwise from geographhic North
                               Real64 SurfAzimuth // or Facing, Direction the surface outward normal faces (degrees)
    );

    Real64 CalcEmmelVertical(Real64 WindAt10m,
                             Real64 WindDir,     // Wind direction measured clockwise from geographhic North
                             Real64 SurfAzimuth, // or Facing, Direction the surface outward normal faces (degrees)
                             int SurfNum);

    Real64 CalcEmmelRoof(Real64 WindAt10m,
                         Real64 WindDir,                // Wind direction measured clockwise from geographhic North
                         Real64 LongAxisOutwardAzimuth, // or Facing, Direction the surface outward normal faces (degrees)
                         int SurfNum);

    Real64 CalcClearRoof(Real64 AirTemp,
                         Real64 WindAtZ,
                         Real64 WindDirect, // Wind direction measured clockwise from geographhic North
                         Real64 RoofArea,
                         Real64 RoofPerimeter,
                         int RoughnessIndex);

    Real64 CalcClearRoof(int SurfNum,
                         Real64 SurfTemp,
                         Real64 AirTemp,
                         Real64 WindAtZ,
                         Real64 WindDirect, // Wind direction measured clockwise from geographhic North
                         Real64 RoofArea,
                         Real64 RoofPerimeter);

    void clear_state();

} // namespace ConvectionCoefficients

} // namespace EnergyPlus

#endif
