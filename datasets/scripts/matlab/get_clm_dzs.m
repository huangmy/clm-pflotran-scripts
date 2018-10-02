% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Returns soil layer thickness for CLM vertical soil layers.
%
% Gautam Bisht (gbisht@lbl.gov)
% 01-06-2014
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function dzs = get_clm_dzs (nlayers)

dzsoi_15lyrs = [
       18
       28
       45
       75
      124
      204
      336
      554
      913
     1506
     2483
     4093
     6748
    11126
    13851]; % [mm]

dzsoi_30lyrs = [
       18
       28
       45
       75
      124
      177
      200
      200
      200
      200
      200
      200
      200
      200
      200
      200
      200
      200
      200
      200
      227
      336
      554
      913
     1506
     2483
     4093
     6748
    11126
    13851
    ]; % [mm]

switch nlayers
    case 10
        dzs = dzsoi_15lyrs(1:10);
     case 30
         dzs = dzsoi_30lyrs;
    otherwise
        error(['Unknown nlayers: ' num2str(nlayers)]);
end

% [mm] --> [m]
dzs = dzs'/1000;  
