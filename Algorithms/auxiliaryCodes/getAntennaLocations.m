%% Copyright(C) 2018
%  Developed By: Muhammet Emin Yanik
%  Advisor: Prof. Murat Torlak
%  The University of Texas at Dallas
%  Department of Electrical and Computer Engineering

%  This work was supported by the Semiconductor Research Corporation (SRC) task 2712.029
%  through The University of Texas at Dallas' Texas Analog Center of Excellence (TxACE).

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice
%%

function [rxAntPos,txAntPos,virtualChPos,distAntennas] = getAntennaLocations(radarType,isFigure)
% Radar Types:
% 'IWR1443'
% 'Simulation'

if (nargin < 2)
    isFigure = 0; % Plot the antenna locations
end
if (nargin < 1)
    radarType = 'IWR1443'; % Default is IWR1443
end

%% Speed of light
%-------------------------------------------------------------------------%
c = physconst('lightspeed');

switch radarType
    case 'IWR1443' % It is AWR1443
        
        fC = 79e9; % center frequency
        lambda = c/fC;
        
        dTxRx = 5e-3;
        % Rx Antenna 1 is the reference
        % Coordinates: [x y z], x-Horizontal, y-Vertical, z-Depth
        rxAntPos = [0 0             0;...
                    0 lambda/2      0;...
                    0 lambda        0;...
                    0 3*lambda/2    0];
        
        txAntPos = [0           3*lambda/2+dTxRx            0;...
                    -lambda/2   3*lambda/2+dTxRx+lambda     0;...
                    0           3*lambda/2+dTxRx+2*lambda   0];
    
        
    case 'Simulation'
        inch2mm = 25.4;
        lambda = inch2mm * 153.56e-3 * 1e-3;
        
        %% Uniform Cascaded 12Tx and 16Rx
        yAxisRxBlock = (0:7).'*lambda/2;
        rxAntPos = [zeros(16,1), [yAxisRxBlock ; yAxisRxBlock+48*lambda], zeros(16,1)];
        yAxisTx = (0:11).'*4*lambda;
        
        dTx_Y = 50e-3;
        txAntPos = [dTx_Y*ones(12,1), 15/4*lambda+yAxisTx, zeros(12,1)];
    
    otherwise
        error('Please enter a correct radar type.')
end


[nRx,~] = size(rxAntPos);
[nTx,~] = size(txAntPos);

txT = reshape(txAntPos,nTx,1,[]);
rxT = reshape(rxAntPos,1,nRx,[]);
virtualChPos = (txT+rxT)/2;
virtualChPos = reshape(permute(virtualChPos,[2,1,3]),[],3);

distAntennas = txT-rxT;
distAntennas = reshape(permute(distAntennas,[2,1,3]),[],3);

%% Plot antenna locations
if isFigure
    figure('OuterPosition',[695 166 670 712]);
    plot(rxAntPos(:,1),rxAntPos(:,2),'s')
    hold on
    plot(txAntPos(:,1),txAntPos(:,2),'o')
    plot((max(virtualChPos(:,1)) + min(virtualChPos(:,1)))/2,(max(virtualChPos(:,2))+min(virtualChPos(:,2)))/2,'x')
    xlabel('x-axis (meters)')
    ylabel('y-axis (meters)')
    title('MIMO Array')
    
    figure('OuterPosition',[695 166 670 712]);
    plot(virtualChPos(:,1),virtualChPos(:,2),'s')
    hold on
    plot((max(virtualChPos(:,1))+min(virtualChPos(:,1)))/2,(max(virtualChPos(:,2))+min(virtualChPos(:,2)))/2,'x')
    xlabel('x-axis (meters)')
    ylabel('y-axis (meters)')
    title('Virtual Array')
end