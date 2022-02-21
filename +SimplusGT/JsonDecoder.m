%
% Decodes json stored in filename and returns the Matlab struct that it
% creates.
%
% Author(s) Rob Oldaker

function [inputData] = JsonDecoder(fileName)
    FID = fopen(fileName);
  	raw = fread(FID,inf); % Reading the contents
  	str = char(raw'); % Transformation
    fclose(FID);
    
    % need to unquote these so they are valid Matlab numeric literals
    str = strrep(str,'"NaN"','NaN');
    str = strrep(str,'"-Inf"','-Inf');
    str = strrep(str,'"Inf"','Inf');
    inputData = jsondecode(str);
end