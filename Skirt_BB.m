close all
% Script to create a pattern block for a size 12 skirt using the Beazley and Bond 
% method of drafting.
% Software implementation copyright Adrian Harwood 2017

%% General Measurements
% Following will be part of a skirt/trouser measurement data class and read
% from sizestream text file.
b_Waist                     = 71.8;
p_Hip                       = 97.7;
q_UpperHip                  = 88.8;
r_Thigh                     = 50.0;
s_Knee                      = 40.0;
t_Ankle                     = 20.0;
u_HipLevel                  = 30.0;
v_KneeLength                = 70.0;
w_BackLength                = 120.0;
x_AnkleLength               = 105.0;
y_OutsideLegLength          = 115.0;
z_FrontLength               = 110.0;
zz_InsideLegLength          = 80.0;

%% Skirt Measurements
% For now add in the standard size 12 measurements, although in future will
% be inferred from the measurement data class members.
a_Waist                     = 70.0;
b_UpperHip                  = 90.0;
c_Hip                       = 96.0;
d_CentreBack                = 60.0;
e_SideSeam                  = 61.0;
f_CentreFront               = 60.0;

%% Easement
% Add ease to all measurements (specific to sized 12 skirt)
a_Waist = a_Waist + 4.0;
b_UpperHip = b_UpperHip + 4.0;
c_Hip = c_Hip + 4.0;

%% Arbitrary Measurements
% Some of the following can be inferred from body scan information but for
% now assume that these follow the empirically driven values.
Arb_WaistLevel = 1.0;       % Ensures the waistline drops by 1cm to allow 
                            % it to curve round the body. This can be 
                            % informed from the body scan.
Arb_UpperHipLevel = 10.0;   % Generic assumption that can in future be
                            % informed from the body scan.
Arb_HipLevel = 20.0;        % Generic assumption that can in future be
                            % informed from the body scan.

% Waist suppression process required calulcation of a front and back dart
% by dividing up the circumference of the waist. For now we assume a fixed
% percentage is assigned to each although this could be adjusted in future.
Arb_BackDartPercent = 0.35;
Arb_FrontDartPercent = 0.20;
Arb_SideSeamPercent = 0.45;

% Dart length is arbitrary but can be inferred from body scan data.
Arb_BackDartLength = 14.0;
Arb_FrontDartLength = 8.0;

%% Assembly
% Points that make up the shape are listed considering the bottom left
% corner of the space to be the origin.

% Create component representing half back of skirt
backBlock1 = Block();

% Add all the fixed points to the block than coincide with the basic
% rectangle
backBlock1.addKeypoint([d_CentreBack, 0]);
backBlock1.addKeypoint([d_CentreBack, c_Hip / 4]);
backBlock1.addKeypoint([Arb_HipLevel, c_Hip / 4]);

% Compute the waistline suppression by finding the difference between 
% the waist measurement and half the hip measurement and then divide by 4
% for a quarter distance.
Int_WaistSupp = (c_Hip - a_Waist) / 4;

% Add point for waist
backBlock1.addKeypointNextTo([Arb_WaistLevel, 0], [d_CentreBack, 0], 'before');

% Add point for suppressed side seam at waist
% Can be computed using side seam percentage of total suppression required
Int_SuppressedSS = (c_Hip / 4) - Arb_SideSeamPercent * Int_WaistSupp;
backBlock1.addKeypointNextTo([0 Int_SuppressedSS], [Arb_WaistLevel, 0], 'before');

% Add the Upper Hip Level point
% Can be computed from the difference between Hip and Upper Hip
Int_SuppressedUpHip = (c_Hip / 4) - (c_Hip - b_UpperHip) / 4;
backBlock1.addKeypointNextTo([Arb_UpperHipLevel, Int_SuppressedUpHip], ...
    [0 Int_SuppressedSS], 'before');

% Add curve between waist point and upper hip point
backBlock1.addCircularCurve([0 Int_SuppressedSS], ...
    [Arb_UpperHipLevel, Int_SuppressedUpHip], 0.5);

% Add back dart
dartEdges = backBlock1.addDart([0 Int_SuppressedSS], [Arb_WaistLevel, 0], ...
    0.5, Arb_BackDartPercent * Int_WaistSupp, Arb_BackDartLength);

% Add curves either side of dart

% TODO: Bug in the gradient directions I need to work through here
backBlock1.addRightAngleCurve([0 Int_SuppressedSS], dartEdges(3,:),...
    [0 Int_SuppressedSS] - [Arb_UpperHipLevel, Int_SuppressedUpHip],...
    dartEdges(2,:) - dartEdges(3,:));

backBlock1.addRightAngleCurve(dartEdges(1,:), [Arb_WaistLevel, 0], ...
    dartEdges(2,:) - dartEdges(1,:), ...
    [d_CentreBack, 0] - [Arb_WaistLevel, 0]);

% Plot block
backBlock1.plotBlock();


