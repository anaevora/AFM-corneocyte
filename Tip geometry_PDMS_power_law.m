% This algorithm considers PDMS AFM data to obtain the parameters 'c' and
% 'n' to describe the geometry of the AFM probe. Young's modulus (E*) was
% obtained using Micromanipulation technique.
% files are after contact point determination with JPK software and in the
% form: "Vertical Deflection" "Segment Time" "Tip-Sample Separation"

%Note power law F=bh^m where 'b' and 'm' are the coefficients found with this script and used to characterise probe geometry.
%These were used to calculate AFM tip geometry parameters as explained in the manuscript (please see main text and SI for detailed derivation aof equations).
%The load index, m, of the power law was calculated 
%at different contact depths, hc, based on a polynomial fit of PDMS loading curves. 
%The plane strain elastic modulus, E^*, of the PDMS specimen required to derive the tip geometry was independently calculated 
%from micromanipulation experiments based on the indentation of 3 different regions. 
%The derived parameters c and n define the local geometry of the AFM tip at different selected contact depths. 

% here information about the data/experiment need to be entered
input_folder = 'C:\Users\anasf\Desktop\Work\AFM\txt'; % where are the data files with CP fitted
% where are the files going to be saved?
output_folder = 'C:\Users\anasf\Desktop\Work\AFM\PDMS'; % name folder
mkdir(output_folder);   % create folder
% what is the working folder for Matlab?
working_folder = 'C:\Users\anasf\Desktop\Work\OliverPharr_Method';

% 1_ open folder and list files
data_folder = cd (input_folder);
D = dir('*.txt');	% make a file list (D) of the (.txt) data in data_folder
[~,index] = sortrows({D.date}.'); D = D(index); clear index     % order data by acquisition time
D_cell = struct2cell(D); D_cell_filename = D_cell(1,:)';	% create cell array of strings with file-names

% 2_ output arrays initialisation
young = zeros(size(D_cell_filename,1),2);      % Young's modulus

% 3_ FOR cycle which opens one file at the time and perform post-processing steps
for i = 1:size(D_cell_filename,1) 
    
    % 3a_ open file
    cd (input_folder);
    myfilename = D_cell_filename{i};
    fileID = fopen(myfilename);
    C = textscan(fileID, '%f%f%f%f', 'CommentStyle', '#');	% raw files contain 4 columns
    mydata = cell2mat(C);	% save data of file(i) into matrix mydata
    fclose(fileID);
    cd (working_folder)
    
    % 3b_ save data from file into arrays
    height = mydata(:,1)*10^9;	% cantilever height [m]
    force = mydata(:,2)*10^9;	% vertical deflection [N]
    segment = mydata(:,3);      % time for extend/retract [s]
    
    segment_start = zeros(5,1);
    jj = 1;
    for ii = 1:length(segment)-1
        if segment(ii)-segment(ii+1) > 0.1
            segment_start(jj,1) = (ii+1);	% index of [segment] change from extend to retract
            jj = jj+1;
        end
    end
    
    % extend (E) data
    force_E = force(1:segment_start(1)-1);
    height_E = height(1:segment_start(1)-1);
    segment_E = segment(1:segment_start(1)-1);
    
    
    % retract (R) data
    force_R = force(segment_start(1):end);
    height_R = height(segment_start(1):end);
    segment_R = segment(segment_start(1):end);
    
    Fmin = min(force_R); %find maximum of force
    fmin1 = find(force_R == Fmin);
    indentation = height_R(fmin1);
  
    
    %Approach curve
    approach = find (height_E < 0);
    height_a = abs(height_E(approach));    % change sign - positive indentation
    force_a = force_E(approach);
    
   
    height_a1 = find(height_a < 250);
    X = height_a(height_a1);
    Y = force_a(height_a1); 
    
    %% Power law fit to PDMS data
    
    % Now we have noisy training data that we can send to fitnlm().
    % Plot the noisy initial data.
    plot(X, Y, 'b*', 'LineWidth', 2, 'MarkerSize', 1);
    grid on;
    % Convert X and Y into a table, which is the form fitnlm() likes the input data to be in.
    tbl = table(X', Y');
    % Define the model as Y = a * (x .^ b) + c
    % Note how this "x" of modelfun is related to big X and big Y.
    % x((:, 1) is actually X and x(:, 2) is actually Y - the first and second columns of the table.
    modelfun = @(b,x) b(1).*x(:, 1).^b(2);  

    beta0 = [1 1]; % Guess values to start with.  Just make your best guess.
    % Now the next line is where the actual model computation is done.
    mdl = fitnlm(X, Y, modelfun, beta0);
    % Extract the coefficient values from the the model object.
    % The actual coefficients are in the "Estimate" column of the "Coefficients" table that's part of the mode.
    coefficients = mdl.Coefficients{:, 'Estimate'};
    % Now we have the coefficients and we can plot y for ANY x, not just the training set.
    % Let's make up a bunch of x (50) from the min to the max.
    xFitted = linspace(min(X), max(X), 50);
    % Create smoothed/regressed data using the model:
    yFitted = coefficients(1)* xFitted.^ coefficients(2); %F=bh^m (power law to define curve behaviour)

    % Now we're done and we can plot the smooth model as a red line going through the noisy blue markers.
    hold on;
    plot(xFitted, yFitted, 'r*-', 'LineWidth', 3);
    grid on;
    title('Power Law Regression with fitnlm()');
    xlabel('X (m)');
    ylabel('Y (N)');
    legendHandle = legend('Data', 'Fitted', 'Location', 'north');
    legendHandle.FontSize = 14;

    message = sprintf('Coefficients for Y = k * X ^ m:\n  k = %8.5f\n  m = %8.5f\n ',...
    coefficients(1), coefficients(2));
    text(2, -8, message, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold', 'Interpreter', 'none');
    % Set up figure properties:
    % Enlarge figure to full screen.
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    % Get rid of tool bar and pulldown menus that are along top of figure.
    % set(gcf, 'Toolbar', 'none', 'Menu', 'none');
    % Give a name to the title bar.
    set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off')
    
    b = coefficients(1);
    m = coefficients(2);
   

%Save in output arrays 
    young(i,1) = b; 
    young(i,2) = m;
    
end

% SAVE
cd(output_folder);
filename1 = 'PDMS.xlsx';
xlswrite(filename1,young);
