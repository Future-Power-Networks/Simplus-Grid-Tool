%
% Saves a Matlab struct as a json file
%
% Author(s): Rob Oldaker
%
% Modified by Yitong Li:
% The option ['Pretty Print',true] in jsonencode is only available for
% matlab 2021a or later. Therefore, this setting is removed and replaced by
% a tedious manual method.
%
function SaveAsJsonToFile(data,fn)
    
    % Check the matlab version
    matlabVersion = version('-release');
    matlabVersion = matlabVersion(1:(end-1));
    matlabVersion= str2double(matlabVersion);
    
    % Call jsonencode
    if matlabVersion>=2021
        json = jsonencode(data,'PrettyPrint',true,'ConvertInfAndNaN',false);
    else
        json = jsonencode(data,'ConvertInfAndNaN',false);
    end
    
    % Make the json file cleaner
    json = strrep(json,'{','\n{\n');
    json = strrep(json,'}','\n}\n');
    json = strrep(json,'[','\n[\n');
    json = strrep(json,']','\n]\n');
    
    % NaN and Inf are *not* valid json literals so put them in quotes to make it valid json
    json = strrep(json,'-Inf','"-Inf"');
    json = strrep(json,'NaN','"NaN"');
    json = strrep(json,'Inf','"Inf"');
    %
    fid = fopen(fn,'wt');
    fprintf(fid, json, + '\n');
    fclose(fid);
end