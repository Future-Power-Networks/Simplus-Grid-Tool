% This function specifies that the library should appear in the Library
% Browser and be cached in the browser repository.

% Author(s): Yitong Li

function blkStruct = slblocks

% Name of the library file saved in the same path
Browser.Library = 'SimplusGT';

% Library name shown in the Simulink library browser
Browser.Name = 'Simplus Grid Tool';

blkStruct.Browser = Browser;

end