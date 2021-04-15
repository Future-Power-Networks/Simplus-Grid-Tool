% This function specifies that the library should appear in the Library
% Browser and be cached in the browser repository.

% Author(s): Yitong Li

function blkStruct = slblocks

% Name of the library file saved in the same path
Browser.Library = 'SimplexPS';

% Library name shown in the Simulink library browser
Browser.Name = 'Simplex Power System';

blkStruct.Browser = Browser;

end