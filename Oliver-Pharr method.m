% This algorithm uses the Oliver-Pharr method to analyse the retract curves
% and calculate the stiffness of the sample as the slope of the tangent at maxium applied force.
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
    height_1 = abs(height_R(hind)); % change sign - positive indentation
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

  %Here must input cantilever parameters based on PDMS nanoindentation
  %Cantilever 2 parameters:

    k30 = 3.826475017;	
    n30 = 8.807985907;
    geo30 = 0.8222;
    c30 = 4.36032E-13;

    k40 = 3.505255143;
    n40 = 7.306158483;
    geo40 =	0.8225;
    c40 = 7.33195E-11;

    k50 = 3.283507921;
    n50 =6.345618053;
    geo50 =	0.8125;
    c50 = 2.01979E-09;

    k60 = 3.094249691;	
    n60 = 5.574961454;
    geo60 = 0.805;
    c60 = 2.94556E-08;

    k80 = 2.835055552;
    n80 = 4.592806994;
    geo80 = 0.7982;
    c80 = 9.29815E-07;

    k100 = 2.694876159;
    n100 = 4.096830291;
    geo100 = 0.7882;
    c100 = 5.55235E-06;

    k120 = 2.62324218;
    n120 =	3.852891776;
    geo120 =0.7824;
    c120 =	1.35293E-05;

    k150 = 2.530843396;
    n150 =3.547724283;
    geo150 = 0.7754;
    c150 = 4.20776E-05;

    k200 = 2.330683953;
    n200 = 2.923170849;
    geo200 = 0.7663;
    c200 =	0.000453102;

    k250 = 2.141581572;
    n250 =	2.378775377;
    geo250 =0.7571;
    c250 =0.003843771;   

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

% SAVE
cd(output_folder);
filename1 = 'YoungsModulus.xlsx';
xlswrite(filename1,young)
cd(output_folder);
