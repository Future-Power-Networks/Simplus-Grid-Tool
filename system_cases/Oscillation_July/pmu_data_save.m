% Structure containing all PMU data
% (e.g., out.pmu1, out.pmu2, etc.)
% Make sure 'out' is loaded/available in your workspace

% List of your timeseries field names
pmuNames = {'pmu1', 'pmu2', 'pmu3', 'pmu6', 'pmu8', 'pmu11', 'pmu12', 'pmu13'};

% Define your custom column names
columnNames = { ...
    'Time', ...
    'Frequency_a_Hz', ...
    'Frequency_b_Hz', ...
    'Frequency_c_Hz', ...
    'Va_mag', ...
    'Va_ang', ...
    'Vb_mag', ...
    'Vb_ang', ...
    'Vc_mag', ...
    'Vc_ang', ...
    'Ia_mag', ...
    'Ia_ang', ...
    'Ib_mag', ...
    'Ib_ang', ...
    'Ic_mag', ...
    'Ic_ang', ...
    'P_3ph', ...
    'Q_3ph' ...
};

% Loop over each PMU
for i = 1:length(pmuNames)
    pmuFieldName = pmuNames{i};
    
    % Access timeseries inside the structure
    pmuData = out.(pmuFieldName);
    
    % Extract time and data
    time = pmuData.Time;    % Column vector
    data = pmuData.Data;    % Matrix (Nx17)
    
    % Combine time and data
    combinedData = [time, data]; % (Nx18)
    
    % Create a table with your custom headers
    T = array2table(combinedData, 'VariableNames', columnNames);
    
    % Prepare the CSV filename
    csvFileName = sprintf('%s.csv', pmuFieldName); % e.g., pmu1.csv
    
    % Open file manually to write the note first
    % fid = fopen(csvFileName, 'w');
    % if fid == -1
    %     error('Cannot open file: %s', csvFileName);
    % end
    
    % Write the sample rate note
    %fprintf(fid, 'Sample Rate 100 Hz\n');
    %fclose(fid);
    
    % Now append the table *with column headers*
    writetable(T, csvFileName, 'WriteMode', 'append', 'WriteVariableNames', true);
    
    fprintf('Saved %s with sample rate note and column names to %s\n', pmuFieldName, csvFileName);
    %fclose(fid);
end
