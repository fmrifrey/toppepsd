%% Function to generate a double-refocussed hyperbolic-tangent velocity selective module
%
% [B1, GLabel, GCont, T] = genDRHT(T, bSection)
%
% in:
%      T            - struct of gradient and RF parameters 
%      bSection     - sections of the module to generate: 'excite',
%                     'refocus', 'VSgrad', 'combine', or 'all'
%      bvelCompCont - flag to use velocity-compensated control
%
% out:
%      B1     - complex B1+ waveform (units of B1max)
%      GLabel - label module gradient waveform  (units of Gmax)
%      GCont  - control module gradient waveform (units of Gmax)
%      T      - updated struct of gradient and RF parameters
%
% T parameter descriptions:
%      B1max  - B1+ maximum amplitdue (units)
%      Gmax   - maximum gradient amplitude (units/cm)
%      RFUP   - RF update time (µs)
%      GUP    - gradient update time (µs)
%      units  - B1+ and gradient units ('G', 'T', 'Hz')
%      f      - gradient flat top time
%      r      - gradient rise time
%      ta[n]  - start of trapezoid attack time for nth gradient
%      td[n]  - end of trapezoid decay time for nth gradient
%      RFe    - duration of excitation pulse
%      RFe_2  - isodelay of excitation pulse
%      RFe2   - start time of last excitation pulse
%      RFr    - duration of refocussing pulse
%      RFr1   - start time of 1st refocussing pulse
%      RFr2   - start time of 2nd refocussing pulse
%      RFrpad - time between the HS pulses
%      All timings in ms
%
% Written by Joseph G. Woods, University of Oxford, April 2022

function [B1, GLabel, GCont, T] = genDRHT(T, bSection, bvelCompCont)

import psdutils.asl.*

if ~exist('T'       ,'var') || isempty(T);        error('T must be specified!'); end
if ~exist('bSection','var') || isempty(bSection); bSection = 'all'; end

% Initialise outputs in case they are not set
B1     = [];
GLabel = [];
GCont  = [];

RUP_GRD_ms = @(A) round(ceil(round(round(A,12)*1e3/T.GUP,9))*T.GUP*1e-3, 3);

% Set the gradient polarities
T.polLabel = [ 1,-1, 1,-1];
if bvelCompCont; T.polCont = [ 1, 1, 1, 1];
else;            T.polCont = [ 0, 0, 0, 0]; end

%% Generate the hard excitation pulses

if strcmpi(bSection,'excite') || strcmpi(bSection,'all')

    FA     = 90;  % Flip angle (degrees)
    phase1 = 0;   % Phase of flip down
    phase2 = 180; % Phase of flip up
    
    % Generate the hard pulses
    T.RFe = 1e-3 * ceil(1e6*FA*pi/180/(T.gamrad*T.B1max)/T.GUP)*T.GUP; % ms
    T.B1_excite1 = genhard(FA, phase1, [], T.B1max, T.RFUP, T.units);
    T.B1_excite2 = genhard(FA, phase2, [], T.B1max, T.RFUP, T.units);
    
    % Approximate isodelay (for off-resonance robustness)
    vsisd = 1.3125;
    T.RFe_2 = T.RFe - RUP_GRD_ms(T.RFe/2/vsisd);

    % No gradients on
    T.grad_excite1 = zeros(size(T.B1_excite1));
    T.grad_excite2 = zeros(size(T.B1_excite2));

end

%% Generate the HT refocussing pulses
% Uses hyperbolic-tanget pulses (see Qin and van Zijl. MRM 2016 https://doi.org/10.1002/mrm.26010)

if strcmpi(bSection,'refocus') || strcmpi(bSection,'all')

    % Set the BIR parameters (Guo and Wong, MRM 2012. http://doi.wiley.com/10.1002/mrm.24145)
    wmax = 42520.0;   % max frequency sweep (hz)
    zeta = 43.58;     % (s^-1)
    tkap = 69.65;     % tan of kapp
    phase = [-90,90]; % Initial phase of refocussing pulses (180° phase between them)
    
    [rho, theta, ~] = genBIR(wmax, zeta, tkap, RUP_GRD_ms(T.RFr/2), T.RFUP);

    T.B1_refocus(:,1) = T.B1max * rho.birmid .* exp(1i * theta.birmid) * exp(1i * phase(1)*pi/180);
    T.B1_refocus(:,2) = T.B1max * rho.birmid .* exp(1i * theta.birmid) * exp(1i * phase(2)*pi/180);
    
    % Generate the 90° phase increments (Liu et al. MRM 2021. https://doi.org/10.1002/mrm.28622)
    T.B1_refocus2 = T.B1_refocus * exp(1i *  90*pi/180);
    T.B1_refocus3 = T.B1_refocus * exp(1i * 180*pi/180);
    T.B1_refocus4 = T.B1_refocus * exp(1i * 270*pi/180);
    
end

%% Generate the velocity encoding gradients

if strcmpi(bSection,'VSgrad') || strcmpi(bSection,'all')
    
    gLabel = genVSGrad(T, T.polLabel);
    gCont  = genVSGrad(T, T.polCont );
    
    % Increase resolution of waveform to match the RF pulses
    T.gLabel = increaseres(gLabel, round(T.GUP/T.RFUP));
    T.gCont  = increaseres(gCont , round(T.GUP/T.RFUP));
    
end
    
%% Combine the whole VS module
    
if strcmpi(bSection,'combine') || strcmpi(bSection,'all')

    % Gaps for VS gradients
    gap1 = zeros(round(T.RFr1*1e3/T.RFUP), 1);                % Gap during G1
    gap2 = zeros(round((T.RFr2-T.RFr1-T.RFr)*1e3/T.RFUP), 1); % Gap during G2 and G3
    gap3 = zeros(round((T.RFe2-T.RFr2-T.RFr)*1e3/T.RFUP), 1); % Gap during G4
    
    %        [     excite;   G1 ;      refocus      ;G2+G3;      refocus      ;  G4 ;      excite ];
    B1     = [T.B1_excite1; gap1; T.B1_refocus(:,1) ; gap2; T.B1_refocus(:,2) ; gap3; T.B1_excite2];

    % Dynamic phase cycling (Liu et al. MRM 2021. https://doi.org/10.1002/mrm.28622)
    T.B1(:,1) = B1;
    T.B1(:,2) = [T.B1_excite1; gap1; T.B1_refocus2(:,1); gap2; T.B1_refocus2(:,2); gap3; T.B1_excite2]; % 2nd phase
    T.B1(:,3) = [T.B1_excite1; gap1; T.B1_refocus3(:,1); gap2; T.B1_refocus3(:,2); gap3; T.B1_excite2]; % 3rd phase
    T.B1(:,4) = [T.B1_excite1; gap1; T.B1_refocus4(:,1); gap2; T.B1_refocus4(:,2); gap3; T.B1_excite2]; % 4th phase
    
    GLabel = [T.grad_excite1; T.gLabel; T.grad_excite2];
    GCont  = [T.grad_excite1; T.gCont ; T.grad_excite2];
    
end
 
end
