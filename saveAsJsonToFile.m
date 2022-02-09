function saveAsJsonToFile(data,fn)
    json = jsonencode(data,'PrettyPrint',true,'ConvertInfAndNaN',false);
    % NaN and Inf are not valid json literals so put them in quotes to make it valid json
    json = strrep(json,'-Inf','"-Inf"');
    json = strrep(json,'NaN','"NaN"');
    json = strrep(json,'Inf','"Inf"');
    %
    fid = fopen(fn,'wt');
    fprintf(fid, json, + '\n');
    fclose(fid);
end