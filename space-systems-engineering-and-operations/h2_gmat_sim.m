%General Mission Analysis Tool(GMAT) Script
%Created: 2026-04-09 00:04:01


%----------------------------------------
%---------- Spacecraft
%----------------------------------------

Create Spacecraft SMOS;
SMOS.DateFormat = UTCGregorian;
SMOS.Epoch = '20 Mar 2000 11:59:28.000';
SMOS.CoordinateSystem = EarthMJ2000Eq;
SMOS.DisplayStateType = Keplerian;
SMOS.SMA = 7060.999999999999;
SMOS.ECC = 0.0116;
SMOS.INC = 98.45;
SMOS.RAAN = 270;
SMOS.AOP = 360;
SMOS.TA = 360;
SMOS.DryMass = 630;
SMOS.Cd = 2.7;
SMOS.Cr = 1.5;
SMOS.DragArea = 20;
SMOS.SRPArea = 20;
SMOS.SPADDragScaleFactor = 1;
SMOS.SPADSRPScaleFactor = 1;
SMOS.AtmosDensityScaleFactor = 1;
SMOS.ExtendedMassPropertiesModel = 'None';
SMOS.NAIFId = -10000001;
SMOS.NAIFIdReferenceFrame = -9000001;
SMOS.OrbitColor = Red;
SMOS.TargetColor = Teal;
SMOS.OrbitErrorCovariance = [ 1e+70 0 0 0 0 0 ; 0 1e+70 0 0 0 0 ; 0 0 1e+70 0 0 0 ; 0 0 0 1e+70 0 0 ; 0 0 0 0 1e+70 0 ; 0 0 0 0 0 1e+70 ];
SMOS.CdSigma = 1e+70;
SMOS.CrSigma = 1e+70;
SMOS.Id = 'SatId';
SMOS.Attitude = CoordinateSystemFixed;
SMOS.SPADSRPInterpolationMethod = Bilinear;
SMOS.SPADSRPScaleFactorSigma = 1e+70;
SMOS.SPADDragInterpolationMethod = Bilinear;
SMOS.SPADDragScaleFactorSigma = 1e+70;
SMOS.AtmosDensityScaleFactorSigma = 1e+70;
SMOS.ModelFile = 'aura.3ds';
SMOS.ModelOffsetX = 0;
SMOS.ModelOffsetY = 0;
SMOS.ModelOffsetZ = 0;
SMOS.ModelRotationX = 0;
SMOS.ModelRotationY = 0;
SMOS.ModelRotationZ = 0;
SMOS.ModelScale = 1;
SMOS.AttitudeDisplayStateType = 'Quaternion';
SMOS.AttitudeRateDisplayStateType = 'AngularVelocity';
SMOS.AttitudeCoordinateSystem = EarthMJ2000Eq;
SMOS.EulerAngleSequence = '321';








%----------------------------------------
%---------- ForceModels
%----------------------------------------

Create ForceModel SMOS_PROPAGATOR_ForceModel;
SMOS_PROPAGATOR_ForceModel.CentralBody = Earth;
SMOS_PROPAGATOR_ForceModel.PrimaryBodies = {Earth};
SMOS_PROPAGATOR_ForceModel.SRP = Off;
SMOS_PROPAGATOR_ForceModel.RelativisticCorrection = Off;
SMOS_PROPAGATOR_ForceModel.ErrorControl = RSSStep;
SMOS_PROPAGATOR_ForceModel.GravityField.Earth.Degree = 4;
SMOS_PROPAGATOR_ForceModel.GravityField.Earth.Order = 4;
SMOS_PROPAGATOR_ForceModel.GravityField.Earth.StmLimit = 100;
SMOS_PROPAGATOR_ForceModel.GravityField.Earth.PotentialFile = 'EGM96.cof';
SMOS_PROPAGATOR_ForceModel.GravityField.Earth.TideModel = 'None';
SMOS_PROPAGATOR_ForceModel.Drag.AtmosphereModel = MSISE90;
SMOS_PROPAGATOR_ForceModel.Drag.HistoricWeatherSource = 'ConstantFluxAndGeoMag';
SMOS_PROPAGATOR_ForceModel.Drag.PredictedWeatherSource = 'ConstantFluxAndGeoMag';
SMOS_PROPAGATOR_ForceModel.Drag.CSSISpaceWeatherFile = 'SpaceWeather-All-v1.2.txt';
SMOS_PROPAGATOR_ForceModel.Drag.SchattenFile = 'SchattenPredict.txt';
SMOS_PROPAGATOR_ForceModel.Drag.F107 = 150;
SMOS_PROPAGATOR_ForceModel.Drag.F107A = 150;
SMOS_PROPAGATOR_ForceModel.Drag.MagneticIndex = 3;
SMOS_PROPAGATOR_ForceModel.Drag.SchattenErrorModel = 'Nominal';
SMOS_PROPAGATOR_ForceModel.Drag.SchattenTimingModel = 'NominalCycle';
SMOS_PROPAGATOR_ForceModel.Drag.DragModel = 'Spherical';

%----------------------------------------
%---------- Propagators
%----------------------------------------

Create Propagator SMOS_PROPAGATOR;
SMOS_PROPAGATOR.FM = SMOS_PROPAGATOR_ForceModel;
SMOS_PROPAGATOR.Type = RungeKutta89;
SMOS_PROPAGATOR.InitialStepSize = 60;
SMOS_PROPAGATOR.Accuracy = 1e-09;
SMOS_PROPAGATOR.MinStep = 0.001;
SMOS_PROPAGATOR.MaxStep = 2700;
SMOS_PROPAGATOR.MaxStepAttempts = 50;
SMOS_PROPAGATOR.StopIfAccuracyIsViolated = true;

%----------------------------------------
%---------- Subscribers
%----------------------------------------

Create ReportFile SIM_OUT;
SIM_OUT.SolverIterations = Current;
SIM_OUT.UpperLeft = [ 0 0 ];
SIM_OUT.Size = [ 0.9984189723320158 0.9548611111111112 ];
SIM_OUT.RelativeZOrder = 53;
SIM_OUT.Maximized = true;
SIM_OUT.Filename = '/home/nico/Documents/POLIMI/YEAR1/SEM2/SSEO/COURSEWORK/HOMEWORK2/gmat-sims/container/RESULTS/long_sim.txt';
SIM_OUT.Precision = 16;
SIM_OUT.Add = {SMOS.A1Gregorian, SMOS.ElapsedSecs, SMOS.ElapsedDays, SMOS.Earth.SMA, SMOS.Earth.RadApo, SMOS.Earth.RadPer};
SIM_OUT.WriteHeaders = true;
SIM_OUT.LeftJustify = On;
SIM_OUT.ZeroFill = Off;
SIM_OUT.FixedWidth = true;
SIM_OUT.Delimiter = ' ';
SIM_OUT.ColumnWidth = 23;
SIM_OUT.WriteReport = true;
SIM_OUT.AppendToExistingFile = false;

%----------------------------------------
%---------- User Objects
%----------------------------------------

Create OpenFramesView CoordinateSystemView1;
CoordinateSystemView1.ViewFrame = CoordinateSystem;
CoordinateSystemView1.ViewTrajectory = Off;
CoordinateSystemView1.InertialFrame = Off;
CoordinateSystemView1.SetDefaultLocation = Off;
CoordinateSystemView1.SetCurrentLocation = Off;
CoordinateSystemView1.FOVy = 45;

Create OpenFramesView EarthView1;
EarthView1.ViewFrame = Earth;
EarthView1.ViewTrajectory = Off;
EarthView1.InertialFrame = Off;
EarthView1.SetDefaultLocation = Off;
EarthView1.SetCurrentLocation = Off;
EarthView1.FOVy = 45;

Create OpenFramesView SMOSView1;
SMOSView1.ViewFrame = SMOS;
SMOSView1.ViewTrajectory = Off;
SMOSView1.InertialFrame = Off;
SMOSView1.SetDefaultLocation = Off;
SMOSView1.SetCurrentLocation = Off;
SMOSView1.FOVy = 45;


%----------------------------------------
%---------- Mission Sequence
%----------------------------------------

BeginMissionSequence;
Propagate SMOS_PROPAGATOR(SMOS) {SMOS.ElapsedDays = 2922};
