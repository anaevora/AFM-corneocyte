% This algorithm uses the Oliver-Pharr method to analyse the retract curves
% and calculate the stiffness of the sample as the slope of the tangent at
% maximum applied force. Viscoelastic properties are also measured using force-relaxation curves.

% It takes the .txt files saved from JPK software as input and 
% give as output the Young's modulus for each file (saved as .xslx file).
% Contact point (CP) determination with JPK software and in the
% form: "Vertical Deflection", "segment time" and "Tip-Sample Separation"
%Note: "Segment time" is necessary to tell difference in txt file between
%loading and unloading curves.

% here information about the data/experiment need to be entered
input_folder = 'C:\Users\anasf\Desktop\Work\AFM\txt'; % where are the data files with CP fitted

v = 0.4; %poisson ratio of keratin

% where are the files going to be saved?
output_folder = 'C:\Users\anasf\Desktop\Work\AFM\viscoelastic'; % name folder
% what is the working folder for Matlab?
working_folder = 'C:\Users\anasf\Desktop\Work\OliverPharr_Method';

% 1_ open folder and list files
data_folder = cd (input_folder);
D = dir('*.txt');	% make a file list (D) of the (.txt) data in data_folder
[~,index] = sortrows({D.date}.'); D = D(index); clear index     % order data by acquisition time
D_cell = struct2cell(D); D_cell_filename = D_cell(1,:)';	% create cell array of strings with file-names

% 2_ output arrays initialisation
young = zeros(size(D_cell_filename,1),6);      % Young's modulus

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
    height = mydata(:,1)*(10^9);	% cantilever height [nm]
    force = mydata(:,2)*(10^9);	% vertical deflection [nN]
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
    
    % Obtain unloading curve with zero at minimum value
    Fmin = min(force_R); %find mminimum of retract force
    fmin1 = find(force_R == Fmin); %find correspondent index on array
    indentation = height_R(fmin1); %find correspondent value of height array
   
    hind = find(height_R < indentation);
    height_1 = abs(height_R(hind));    % change sign - positive indentation
    force_1 = force_R(hind);
    
    ind = find(force_1 > 0);
    height_ind = height_1(ind);
    force_ind = force_1(ind);


    % consider indentation region only
    
    indent = find (height_E < 0);
    height_comp1 = abs(height_E(indent));    % change sign - positive indentation
    force_comp1 = force_E(indent);
    time_1 = segment(indent);
    time_ind = time_1 - time_1(1);
    
    hmax_load1 = max(height_comp1); %maximum indentation depth
    rising_time = max(time_ind) - time_ind(1); %rising time from contact point to hmax
    speed = hmax_load1/rising_time; % speed of compression in um
    
 
    %% Polynomial fit to unloading curve

    % Find tangent to Fmax (maximum force of force vs indentation curve)
    Fmax = find(max(force_ind)); %find maximum of force
    fmaxI = force_ind(Fmax);
    fmaxE = max(force_E);
    
    n = 2; %polynomial of degree 2
    [p,~,mu] = polyfit(height_ind,force_ind,n); % degree 2 polynomial fit, centered and scaled
   
    % Plot the fit with the data
    xfit = linspace(min(height_ind),max(height_ind),500);
    yfit = polyval(p,xfit,[],mu);
    %plot(xfit,yfit,'r') %to plot polynomial fit - check if all is ok
    
    x0 = height_ind(Fmax);
    f = @(height_ind) polyval(p,height_ind,[],mu);
    df = derivest(f,x0); % Here is the derivative
    
    % Calculate the end points of a tangent line
    xt = height_ind; %[x0+0.5 x0+0.5];
    yt = f(x0) + (xt-x0)*df;
    %plot(xt,yt) %to plot tangent - confirm it is all ok)
    %xline(x0)   %vertical line meaning maximum indentation depth
    
    %Calculate R^2 (goodness of polyvat fit)
    y_coeff = p(1)*height_ind.^2+p(2)*height_ind+p(3);
    fit_coeff = [height_ind force_ind y_coeff];
    R_matrix = corrcoef(fit_coeff);
    R_fit_height_ind_F = R_matrix(2,1); 
    R_fit_polyfit_F = R_matrix(3,2);

    %% Calculate biomechanical parameters 

    %Calculate slope of tangent -> stiffness of sample at Fmax
    coefficients = polyfit(xt, yt, 1);
    S = abs(coefficients(1));
    
    %Calculate contact depth band area of contaxct based on AFM tip
    %geometry (parameters calculated from PDMS calibration)

    hmax_unload= max(height_ind);

    %Parameters of cantilevers obtained based on PDMS nanoindentation experiments

	%Cantilever x parameters
 	%Note: should be updated for the cantilever used! Parameters are calculated using an excel file 
  	%(can be found as SI info: "Tip geometry parameters")

    k30 = 4.130775;
    n30 = 9.459987;
    geo30 = 0.831133;
    c30 = 1.42283E-14;

    k40 = 3.914392;
    n40 = 9.241828;
    geo40 =	0.825094;
    c40 = 5.28943E-20;

    k50 = 3.769292;
    n50 = 8.531067;
    geo50 = 0.820819;
    c50 = 1.76948E-18;

    k60 = 4.0507;
    n60 = 9.93427;
    geo60 = 0.828946;
    c60 = 2.81845E-21;

    k80 = 3.933467;
    n80 = 9.337248;
    geo80 = 0.825642;
    c80 = 2.51162E-20;

    k100 = 4.10022;
    n100 = 10.19109;
    geo100 =0.830304;
    c100 = 3.37467E-22;

    k120 = 4.262689;
    n120 =	11.05694;
    geo120 = 0.83463;
    c120 = 4.22925E-24;

    k150 = 4.092838;
    n150 = 10.15255;
    geo150 = 0.830102;
    c150 = 3.62795E-22;

    k200 = 3.811338;
    n200 = 8.734277;
    geo200 = 0.822077;
    c200 = 6.04279E-19;

    k250 = 3.70194;
    n250 = 8.210221;
    geo250 = 0.818769;
    c250 = 8.57891E-18;


    %calculate area of contact based on maximum indentation depth to
    %account for complex AFM tip geometry

    if hmax_unload<30
           
        hs = geo30*(fmaxI/S);
        hc = hmax_unload-hs; % contact depth based on contact radius
        Area = pi*((hc/c30)^(2/n30)); % area of contact based on contact depth and parameter c (from PDMS fit)
        constant = pi*((hmax_unload/c30)^(2/n30));

    elseif hmax_unload > 30 && hmax_unload < 40

        hs = geo40*(fmaxI/S);
        hc = hmax_unload-hs; % contact depth based on contact radius
        Area = pi*((hc/c40)^(2/n40)); % area of contact based on contact depth and parameter c (from PDMS fit)
        constant = pi*((hmax_unload/c40)^(2/n40));

    elseif hmax_unload > 40 && hmax_unload < 50
        
        hs = geo40*(fmaxI/S);
        hc = hmax_unload-hs; % contact depth based on contact radius
        Area = pi*((hc/c50)^(2/n50)); % area of contact based on contact depth and parameter c (from PDMS fit)
        constant = pi*((hmax_unload/c50)^(2/n50));

    elseif hmax_unload > 50 && hmax_unload < 60
        
        hs = geo60*(fmaxI/S);
        hc = hmax_unload-hs; % contact depth based on contact radius
        Area = pi*((hc/c60)^(2/n60)); % area of contact based on contact depth and parameter c (from PDMS fit)
        constant = pi*((hmax_unload/c60)^(2/n60));

    elseif hmax_unload > 60 && hmax_unload < 80
        
        hs = geo80*(fmaxI/S);
        hc = hmax_unload-hs; % contact depth based on contact radius
        Area = pi*((hc/c80)^(2/n80)); % area of contact based on contact depth and parameter c (from PDMS fit)
        constant = pi*((hmax_unload/c80)^(2/n80));

   elseif hmax_unload > 80 && hmax_unload < 100
        
        hs = geo100*(fmaxI/S);
        hc = hmax_unload-hs; % contact depth based on contact radius
        Area = pi*((hc/c100)^(2/n100)); % area of contact based on contact depth and parameter c (from PDMS fit)
        constant = pi*((hmax_unload/c100)^(2/n100));

    elseif hmax_unload > 100 && hmax_unload < 120
        
        hs = geo120*(fmaxI/S);
        hc = hmax_unload-hs; % contact depth based on contact radius
        Area = pi*((hc/c120)^(2/n120)); % area of contact based on contact depth and parameter c (from PDMS fit)
        constant = pi*((hmax_unload/c120)^(2/n120));

    elseif hmax_unload > 120 && hmax_unload < 150
        
        hs = geo150*(fmaxI/S);
        hc = hmax_unload-hs; % contact depth based on contact radius
        Area = pi*((hc/c150)^(2/n150)); % area of contact based on contact depth and parameter c (from PDMS fit)
        constant = pi*((hmax_unload/c150)^(2/n150));

     
    elseif hmax_unload > 150 && hmax_unload < 200
        
        hs = geo200*(fmaxI/S);
        hc = hmax_unload-hs; % contact depth based on contact radius
        Area = pi*((hc/c200)^(2/n200)); % area of contact based on contact depth and parameter c (from PDMS fit)
        constant = pi*((hmax_unload/c200)^(2/n200));

    elseif hmax_unload > 200 && hmax_unload < 600
        
        hs = geo250*(fmaxI/S);
        hc = hmax_unload-hs; % contact depth based on contact radius
        Area = pi*((hc/c250)^(2/n250)); % area of contact based on contact depth and parameter c (from PDMS fit)
        constant = pi*((hmax_unload/c250)^(2/n250));

    end


    %Calculate Young's modulus based on area of contact and stiffness at
    %fmax
    E = (sqrt(pi)/2)*(1-v^2)*(S/sqrt(Area))*(10^3); % in MPa

    
    % plot original height_E, force_E and new graph based on the contactpoint
    %figure('OuterPosition',[460 250 550 450])
    %plot(height_ind, force_ind, 'c-')
    %hold on
    %plot(xfit,yfit,'r')
    %hold on
    %plot(xt,yt)
    %legend('Data','Oliver-Phaar', 'Tangent to Fmax', 'Location','northwest')
    %xlabel('Indentation depth (nm)')
    %ylabel('Force (nN)')
    
   
    %Save in output arrays 
    young(i,1) = E; %Young's modulus in MPa
    young(i,2) = S; %stiffness of cell (N/m)
    young(i,3) = hmax_unload; %maximum contact depth (nm);
    young(i,4) = hc; %real contact depth;
    young(i,5) = fmaxI; %maximum applied force at loading
    young(i,6) = constant; %constant to be use in the calculus of Cs
    
    
end


%%
%Enter information of about data/experiment for stress relaxation curve
input_folder = 'C:\Users\anasf\Desktop\WORK\AFM\pause'; % where are the data files
% where are the files going to be saved?
output_folder = 'C:\Users\anasf\Desktop\WORK\AFM\viscoelastic'; % name folder
mkdir(output_folder);   % create folder
% what is the working folder for Matlab?
working_folder = 'C:\Users\anasf\Desktop\WORK\Viscoplastic_Model';

% 1_ open folder and list files
data_folder = cd (input_folder);
D = dir('*.txt');	% make a file list (D) of the (.txt) data in data_folder
[~,index] = sortrows({D.date}.'); D = D(index); clear index     % order data by acquisition time
D_cell = struct2cell(D); D_cell_filename = D_cell(1,:)';	% create cell array of strings with file-names

% 4_ output arrays initialisation
prony = zeros(size(D_cell_filename,1),10);  %Prony series fitting parameters
data = zeros(size(D_cell_filename,1),10);


% 5_ FOR cycle which opens one file at the time and perform post-processing steps
for i = 1:size(D_cell_filename,1) 
    
    % 5a_ open file
    cd (input_folder);
    myfilename = D_cell_filename{i};
    fileID = fopen(myfilename);
    C = textscan(fileID, '%f%f%f%f', 'CommentStyle', '#');	% raw files contain 4 columns
    mydata = cell2mat(C);	% save data of file(i) into matrix mydata
    fclose(fileID);
    cd (working_folder)
    
    % 5b_ save data from file into arrays
    force_pause = mydata(:,2)*(10^9);	% vertical deflection [nN]
    time_pause = mydata(:,3);       % time [s]
    segment = mydata(:,4);      % time for extend/retract [s]
    displacement_1 = mydata(:,1)*(10^9);	% displacement [nm]
        
%6_Fit prony series to F vs time data (relaxation data)
    tempo = time_pause;
    t = time_pause - max(time_1);
    y = force_pause;

    %t_1 = time_pause - max(time_1);
    %t_2 = find (t_1 <2);
    %t = t_1(t_2);
    %y = force_pause(t_2);

    F = @(x,t) x(1)+x(2)*exp(-t/x(3)) + x(4)*exp(-t/x(5));
    %x0 = [1251.32 135.82 0.0141 179.97 0.70];
    %x0 =[1477.56 100.619 0.05 150.601	1.30]
    %x0 = [2054.10 500.94 0.01 130.66 0.2];

    x0 = [1634.46 2559 0.1 800 1.5];


    lb = [0, 0, 0, 0, 0];
    xunc = lsqcurvefit(F, x0, t, y, lb);
    tlist = linspace(min(t), max(t));   % Plot Finer Resolution
    figure(1)
    plot(t,y,'bo')
    title('Data points')
    hold on
    plot(tlist, F(xunc,tlist), '-r', 'LineWidth', 2) 
    hold off
   
    %5_Coefficients of fitting prony series
    
    B0 = xunc(1);
    B1 = xunc(2);
    B2 = xunc(4);
    tau1 = xunc(3);
    tau2 = xunc(5);
    
    %Save array of coefficients
    
    prony(i,1) = B0;
    prony(i,2) = B1;
    prony(i,3) = B2;
    prony(i,4) = tau1;
    prony(i,5) = tau2;
        
   
    for index = 1:size(D_cell_filename)
       
          b_0 = prony(i,1);
          b_1 = prony(i,2);
          b_2 = prony(i,3);
          constant = young(i, 6);

    %Calculate Cs used to calculate hardness values      
        C0 = b_0/constant;
        C1 = b_1/constant;
        C2 = b_2/constant;
    
    %Calculate hardness values      
        H0 = (C0+C1+C2)*(10^3); % in MPa
        H_eq = (C0)*(10^3); %in MPa

    prony(i,6) = C0;
    prony(i,7) = C1;
    prony(i,8) = C2;
    prony(i,9) = H0;
    prony(i,10) = H_eq;
    end

      
    %Save array of coefficients
    
    data(:,1) = prony(:,1); %B0
    data(:,2) = prony(:,2); %B1
    data(:,3) = prony(:,3); %B2
    data(:,4) = prony(:,6); %C0
    data(:,5) = prony(:,7); %C1
    data(:,6) = prony(:,8); %C2
    data(:,7) = prony(:,4); %tau1
    data(:,8) = prony(:,5); %tau2
    data(:,9) = prony(:,9); %H0
    data(:,10) = prony(:,10); %Heq
end
   
% SAVE
cd(output_folder);
filename1 = 'Hardness.xlsx';
xlswrite(filename1,data)
cd(output_folder);
filename1 = 'YoungsModulus.xlsx';
xlswrite(filename1,young)
cd(output_folder);

