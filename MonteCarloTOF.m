% Clear workspace and command window
clear, clc
format long

% Set the loop variable
loop = true;

% Ensures that the menu reopens when finished with subprogram
while loop

% Opens dialogue box    
answer = questdlg( ...
    'Please select an option', ...
	'SIMION Monte Carlo Simulation', ...
	'Create Data','Format Data','Show Data', ...
    'Create Data' ...
);

% Handle response
switch answer
    
    % Create data in .ion file for inserting into SIMION
    case 'Create Data'
        
        % Input the number of ions
        numIons = inputdlg( ...
            'Number of Ions', ...
            'Generate ION File', ...
            [1 50], ...
            {'100'} ...
        );

        % Input the time of birth range
        tob = inputdlg( ...
            {'Time of Birth Min', 'Time of Birth Max'}, ...
            'Generate ION File', ...
            [1 50; 1 50], ...
            [{'0'}, {'0'}] ...
        );

        % Input the mass range
        mass = inputdlg({'Mass Min', 'Mass Max'}, ...
            'Generate ION File', ...
            [1 50; 1 50], ...
            [{'100'}, {'100'}] ...
        );
        
        % Input the charge range
        charge = inputdlg( ...
            {'Charge Min', 'Charge Max'}, ...
            'Generate ION File', ...
            [1 50; 1 50], ...
            [{'1'}, {'1'}] ...
        );
       % Input the x position range
        x = inputdlg( ...
            {'x Min', 'x Max'}, ...
            'Generate ION File', ...
            [1 50; 1 50], ...
            [{'3.5'}, {'3.5'}] ...
        ); 

        % Input the y position range
        y = inputdlg( ...
            {'y Min', 'y Max'}, ...
            'Generate ION File', ...
            [1 50; 1 50], ...
            [{'0'}, {'0'}] ...
        );

        % Input the z position range
        z = inputdlg( ...
            {'z Min', 'z Max'}, ...
            'Generate ION File', ...
            [1 50; 1 50], ...
            [{'0'}, {'0'}] ...
        ); 

        % Input the azimuthal angle range
        azimuth = inputdlg( ...
            {'Azimuth Min', 'Azimuth Max'}, ...
            'Generate ION File', ...
            [1 50; 1 50], ...
            [{'0'}, {'0'}] ...
        ); 

        % Input the elevation angle range
        elevation = inputdlg( ...
            {'Elevation Min', 'Elevation Max'}, ...
            'Generate ION File', ...
            [1 50; 1 50], ...
            [{'0'}, {'0'}] ...
        ); 

        % Input the kinetic energy range
        ke = inputdlg( ...
            {'Kinetic Energy Min', 'Kinetic Energy Max'}, ...
            'Generate ION File', ...
            [1 50; 1 50], ...
            [{'0'}, {'100'}] ...
        ); 

        % Input the charge weighting factor range
        cwf = inputdlg( ...
            {'Charge Weighting Factor Min', 'Charge Weighting Factor Max'}, ...
            'Generate ION File', ...
            [1 50; 1 50], ...
            [{'1'}, {'1'}] ...
        ); 

        % Input the colour range
        colour = inputdlg( ...
            {'Colour Min', 'Colour Max'}, ...
            'Generate ION File', ...
            [1 50; 1 50], ...
            [{'1'}, {'1'}] ...
        );

        % Randomly generate the requested number of ions within the ranges
        % of the parameters
        tob       = str2double(cell2mat(tob(1,1)))       + (str2double(cell2mat(tob(2,1)))       - str2double(cell2mat(tob(1,1))))       .* rand(str2double(cell2mat(numIons(1,1))),1);
        mass      = str2double(cell2mat(mass(1,1)))      + (str2double(cell2mat(mass(2,1)))      - str2double(cell2mat(mass(1,1))))      .* rand(str2double(cell2mat(numIons(1,1))),1);
        charge    = str2double(cell2mat(charge(1,1)))    + (str2double(cell2mat(charge(2,1)))    - str2double(cell2mat(charge(1,1))))    .* rand(str2double(cell2mat(numIons(1,1))),1);
        x         = str2double(cell2mat(x(1,1)))         + (str2double(cell2mat(x(2,1)))         - str2double(cell2mat(x(1,1))))         .* rand(str2double(cell2mat(numIons(1,1))),1);
        y         = str2double(cell2mat(y(1,1)))         + (str2double(cell2mat(y(2,1)))         - str2double(cell2mat(y(1,1))))         .* rand(str2double(cell2mat(numIons(1,1))),1);
        z         = str2double(cell2mat(z(1,1)))         + (str2double(cell2mat(z(2,1)))         - str2double(cell2mat(z(1,1))))         .* rand(str2double(cell2mat(numIons(1,1))),1);
        azimuth   = str2double(cell2mat(azimuth(1,1)))   + (str2double(cell2mat(azimuth(2,1)))   - str2double(cell2mat(azimuth(1,1))))   .* rand(str2double(cell2mat(numIons(1,1))),1);
        elevation = str2double(cell2mat(elevation(1,1))) + (str2double(cell2mat(elevation(2,1))) - str2double(cell2mat(elevation(1,1)))) .* rand(str2double(cell2mat(numIons(1,1))), 1);
        ke        = str2double(cell2mat(ke(1,1)))        + (str2double(cell2mat(ke(2,1)))        - str2double(cell2mat(ke(1,1))))        .* rand(str2double(cell2mat(numIons(1,1))),1);
        cwf       = str2double(cell2mat(cwf(1,1)))       + (str2double(cell2mat(cwf(2,1)))       - str2double(cell2mat(cwf(1,1))))       .* rand(str2double(cell2mat(numIons(1,1))),1);
        colour    = str2double(cell2mat(colour(1,1)))    + (str2double(cell2mat(colour(2,1)))    - str2double(cell2mat(colour(1,1))))    .* rand(str2double(cell2mat(numIons(1,1))),1);

        % OVERWRITE
        azimuth   = 360 * rand(str2double(cell2mat(numIons(1,1))), 1);
        elevation = acosd(2 * rand(str2double(cell2mat(numIons(1,1))), 1) - 1) - 90;
        
        KEx = -0.02585 * log(rand(str2double(cell2mat(numIons(1,1))), 1));
        KEy = -0.02585 * log(rand(str2double(cell2mat(numIons(1,1))), 1));
        KEz = -0.02585 * log(rand(str2double(cell2mat(numIons(1,1))), 1));
        
        ke = KEx + KEy + KEz;
        
%         vx = 9822.69385*(2*KEx./mass).^(0.5);
%         vy = 9822.69385*(2*KEy./mass).^(0.5);
%         vz = 9822.69385*(2*KEz./mass).^(0.5);
%         
%         rx = 1-2*floor(2*rand(str2double(cell2mat(numIons(1,1))),1));
%         ry = 1-2*floor(2*rand(str2double(cell2mat(numIons(1,1))),1));
%         rz = 1-2*floor(2*rand(str2double(cell2mat(numIons(1,1))),1));
%         
%         vx = rx.*vx;
%         vy = ry.*vy;
%         vz = rz.*vz;
        
%         vr2d = (vx.^2+vy.^2).^(0.5);
%         vr = (vx.^2+vy.^2+vz.^2).^(0.5);
%         
%         azimuth = atand(vy./vx);
%         elevation = cosd(90-atand(vz./vr2d));
        
%         [azimuth,elevation,v] = cart2sph(vx,vy,vz);
%         
%         h1 = figure;
%         h1 = histogram(180*azimuth./pi,10);
%         h2 = figure;
%         h2 = histogram(cos(elevation+pi/2),10);
%         h3 = figure;
%         h3 = histogram(v,10);
%         waitfor(h1);
%         waitfor(h2);
%         waitfor(h3);
        
        % Write the data to a ION file
        fileContents = [tob, mass, charge, x, y, z, azimuth, elevation, ke, cwf, colour];
        csvwrite('..\simion\matlab.ion', fileContents)
    
    % Process the output file from SIMION
    case 'Format Data'
        
        % Select the output file
        [file, path] = uigetfile('..\simion\*.txt');
        
        % Set the table column names
        opts = detectImportOptions(fullfile(path, file));
        opts.VariableNames = {'TOF', 'a', 'X', 'b', 'Vx', 'c', 'd', 'e'};
        
        % Read the file into a table
        T = readtable(fullfile(path, file), opts);
        
        % Remove the unused table columns
        T = removevars(T, {'a', 'b', 'c', 'd', 'e'});
        
        % Convert the table into an array
        T = table2array(T);
        
        % Find the start position of the data from the last run
        dataStart = find(strcmp(T(:, 2), 'Next')', 1, 'last') + 1;
        
        % Find the end position of the data from the last run
        dataEnd = size(T, 1) - 1;
        
        % Create a new array with only the data from the last run
        T = T(dataStart:dataEnd, :);
        
        % Create a new array with twice the number of columns as T
        % C = zeros(6, dataEnd-dataStart);
        C = ['', '', '', '', '', '';];
        
        % Loop through every second row in the table
        for ii = 1:2:size(T, 1)
            
            % Combine the current row and next row of T into one row of C
            C = [C; T(ii, :), T(ii+1, :)];
            % C((ii+1)/2) = [T(ii,:), T(ii+1,:)];
            
        end
        
        % Remove the unused array columns
        C(:,6) = [];
        C(:,5) = [];
        C(:,1) = [];
        
        % Extract everything after the left bracket in each cell
        C = extractAfter(C, "(");
        
        % Convert each cell from a string to a double
        C = cellfun(@str2double, C);
        
        % Input the simulation parameters
        parameters = inputdlg( ...
            {'Potential Difference of the Acceleration Region', ...
                'Distance of the Acceleration Region'}, ...
            'Format Data', ...
            [1 50; 1 50], ...
            [{'2000'}, {'10'}] ...
        );
        
        % Set the parameters
        mass = 100;
        charge = 1;
        parameters = cellfun(@str2double, parameters);
        
        % Calculate the acceleration experience by the ions
        a = (((charge * 1.6E-19) * parameters(1)) / ((mass * 1.66E-27) * (parameters(2) * 1E-3))) * 1E-9;
        
        % Calculate the distance from the midpoint of the acceleration
        % region to the start of the drift region
        s0 = parameters(2) / 2;
        
        % Calculate the beta parameter
        beta = (C(:, 2).^2 + 2 * a .* ((s0 + 1) - C(:, 1))) / (2 * a * s0);
        
        % Calculate the turnaround time
        turnaroundTime = (C(:, 3) + C(:, 2) ./ a) * 1000;
        
        % Create an array of the data to be saved to the file
        outputData = [beta, turnaroundTime];
        
        % Input the output file name
        fileName = inputdlg( ...
            {'Enter File Name'}, ...
            'Format Data', ...
            [1 50], ...
            {''} ...
        );
        
        % Set the output file name
        outputFile = strcat('..\MATLAB\data\', fileName{1});
    
        % Save the data to the selected file
        csvwrite(outputFile, outputData)
        
        
        % Data not saved in long format?
        % Data needed: mass,charge,potential difference,distance between
        % plates
        
        
%         % Plot beta against turnaround time
%         s = scatter(beta, turnaroundTime, 2, 'k','filled');
%         
        % Calculate the number of ions used in the simulation
        numIons = (dataEnd - dataStart + 1) / 2;
%         
%         % Set the title and axes labels of the plot
%         title(sprintf('Monte Carlo Simulation of %.f ions', numIons))
%         xlabel('Beta')
%         ylabel('Time of Flight / nanoseconds')
%         
%         % Wait until the figure has been closed
%         waitfor(s)
        
        
        % Set the value for the nanoseconds per channel
        nspc = 1;
        
        % Convert the time of flight for each ion into nanoseconds and sort
        % in ascending order
        time = sort(C(:, 3) * 1000);
        
        % Calculate the number of bars needed for the histogram
        bars = round((time(end) - time(1)) / nspc);
        
        % Store the time of flight histogram data to an array
        % n -> frequency
        % x -> time of flight
        [n, x] = hist(time, bars);
        % n_normalised = n/numel(C(:,3))/(x(2)-x(1)); %normalise to unit area
        
        % Initialise the variables for the plot
        graphStart = 475;
        graphEnd = 500;
        t = graphStart:nspc:graphEnd;
        y = zeros(1, length(t));
        
        % Create a vector using the initialised variables
        vector = [t', y'];
        
        % Create a vector using the histogram data
        vector2 = [x', n'];
        
        % Calculate the start and end positions of the histogram data with
        % respect to the initialised variables
        startPos = floor(vector2(1, 1) / nspc) - graphStart;
            %endPos = ceil(vector2(end,1)/nspc)+1;
        endPos = startPos + size(vector2, 1) + 1;
        
        % Replace the initialised data with the histogram data
        %    vector(startPos:endPos,:) = [vector2(1,1)-nspc, 0; vector2; vector2(end,1)+nspc, 0];
        
        % Create an array of the data to be saved to the file
        %    outputData = [vector(:,1), vector(:,2)];
        outputData = [vector2(1, 1) - 2 * nspc, 0; vector2(1, 1) - nspc, 0; vector2; vector2(end, 1) + nspc, 0; vector2(end, 1) + 2 * nspc, 0];
        
        % Set the output file name
        outputFile = strcat('..\MATLAB\data\histogram\', fileName{1});
    
        % Save the data to the selected file
        csvwrite(outputFile, outputData)
        
        % Plot the time of flight histogram curve
        %    h = plot(vector(:,1), vector(:,2), 'r');
        h = plot(outputData(:, 1), outputData(:, 2), 'r');
        
        % Set the title and axes labels of the plot
        title(sprintf('Monte Carlo Simulation of %.f ions', numIons))
        xlabel('Time of Flight / nanoseconds')
        ylabel('Frequency')
        
        % Wait until the figure has been closed
        waitfor(h)
    
    % Display the data from processed output files
    case 'Show Data'
        
        % Select the output file
        [file,path] = uigetfile('..\MATLAB\data\histogram\*.csv');

        % Set the table column names
        opts = detectImportOptions(fullfile(path, file));
        opts.VariableNames = {'tof', 'freq'};
        
        % Read the file into a table
        T = readtable(fullfile(path, file), opts);
        
        % Plot the time of flight histogram curve
        h = plot(T.tof, smooth(T.freq), 'r');
        
%         halfFreq = (max(T.freq)/2);
%         
%         hold on
%         plot(1:3000, halfFreq*ones(1,3000))
        
        dim = [.55 .3 .3 .55];
        str = sprintf('FWHM frequencies = %.f', halfFreq);
        annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');
        
        % Set the title and axes labels of the plot
        title(sprintf('Monte Carlo Simulation of %.f ions', height(T)))
        xlabel('Time of Flight / nanoseconds')
        ylabel('Frequency')
        
        % Wait until the figure has been closed
        waitfor(h)
    
    % Window closed without selecting any option    
    case ''
        
        % Exit the menu loop
        loop = false;
        
end

end