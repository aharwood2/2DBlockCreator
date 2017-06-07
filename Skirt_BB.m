close all
clear
% Script to create a pattern block for a size 12 skirt using the Beazley
% and Bond method. This software is a proof of concept and should be
% straightforward to extend to other patterns or method as required.
% Software implementation copyright Adrian Harwood 2017.
% The University of Manchester, UK.

%% General Measurements
% Following will be part of a skirt/trouser measurement data class and read
% from a body scanner text file.
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
% be inferred from the body scan data.
a_Waist                     = 70.0;
b_UpperHip                  = 90.0;
c_Hip                       = 96.0;
d_CentreBack                = 60.0;
e_SideSeam                  = 61.0;
f_CentreFront               = 60.0;

%% Easement
% Add ease to all measurements (specific to size 12 skirt in this case).
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

% Dart placement is also arbitrary and is specified as a percentage of
% quarter waist as measured from the start point of the waist (using strict
% connectivity order)
Arb_BackDartPlacement = 0.5;
Arb_FrontDartPlacement = 1 / 3;

%% Assembly
% Points that make up the shape are listed in a strict anti-clokwise order
% to maintain correct connectivity for plotting. The bottom left corner of
% the space to be the origin.

% Create component representing half back of skirt folded in half.
backBlock = Block();

% Add all the fixed points to the block that coincide with the basic
% rectangle. These points do not move throughout the drafting process.
backBlock.addKeypoint([d_CentreBack, 0]);
backBlock.addKeypoint([d_CentreBack, c_Hip / 4]);
backBlock.addKeypoint([Arb_HipLevel, c_Hip / 4]);

% Compute the waistline suppression by finding the difference between 
% the waist measurement and half the hip measurement and then divide by 4
% for a quarter distance.
Int_WaistSupp = (c_Hip - a_Waist) / 4;

% Add point for waist line drop.
backBlock.addKeypointNextTo(...
    [Arb_WaistLevel, 0], ...
    [d_CentreBack, 0], 'before');

% Add point for suppressed side seam at waist.
% Can be computed using side seam percentage of total suppression required.
Int_SuppressedSS = (c_Hip / 4) - Arb_SideSeamPercent * Int_WaistSupp;
backBlock.addKeypointNextTo(...
    [0 Int_SuppressedSS], ...
    [Arb_WaistLevel, 0], 'before');

% Add the suppressed Upper Hip Level point.
% Can be computed from the difference between Hip and Upper Hip
Int_SuppressedUpHip = (c_Hip / 4) - (c_Hip - b_UpperHip) / 4;
backBlock.addKeypointNextTo(...
    [Arb_UpperHipLevel, Int_SuppressedUpHip], ...
    [0 Int_SuppressedSS], 'before');

% Add curve between waist point and upper hip point as per BB method.
% Assume for now, in the absence of vary form curve that this is a curve
% defined by a circle.
backBlock.addCircularCurve(...
    [Arb_UpperHipLevel, Int_SuppressedUpHip], ...
    [0 Int_SuppressedSS], 0.5, true);

% Trace off block
frontBlock = Block(backBlock);

% Add back dart.
dartEdges = backBlock.addDart(...
    [0 Int_SuppressedSS], [Arb_WaistLevel, 0], ...
    Arb_BackDartPlacement, ...
    Arb_BackDartPercent * Int_WaistSupp, Arb_BackDartLength);

% Add curves either side of dart ensuring the curve intersects the joining
% edges at a right angle.
backBlock.addRightAngleCurve(...
    [0 Int_SuppressedSS], dartEdges(1,:),...
    [Arb_UpperHipLevel, Int_SuppressedUpHip] - [0 Int_SuppressedSS],...
    dartEdges(2,:) - dartEdges(1,:),...
    [true, true]);

backBlock.addRightAngleCurve(...
    dartEdges(3,:), [Arb_WaistLevel, 0], ...
    dartEdges(2,:) - dartEdges(3,:), ...
    [d_CentreBack, 0] - [Arb_WaistLevel, 0],...
    [true, true]);

% Add front dart
dartEdges = frontBlock.addDart(...
    [0 Int_SuppressedSS], [Arb_WaistLevel, 0], ...
    Arb_FrontDartPlacement, ...
    Arb_FrontDartPercent * Int_WaistSupp, Arb_FrontDartLength);

% Add curves
frontBlock.addRightAngleCurve(...
    [0 Int_SuppressedSS], dartEdges(1,:),...
    [Arb_UpperHipLevel, Int_SuppressedUpHip] - [0 Int_SuppressedSS],...
    dartEdges(2,:) - dartEdges(1,:), [true, true]);

frontBlock.addRightAngleCurve(...
    dartEdges(3,:), [Arb_WaistLevel, 0], ...
    dartEdges(2,:) - dartEdges(3,:), ...
    [d_CentreBack, 0] - [Arb_WaistLevel, 0], [true, true]);

% Plot blocks
figure
subplot(2,1,1), backBlock.plotBlock();
title('Back Block')
subplot(2,1,2), frontBlock.plotBlock();
title('Front Block')