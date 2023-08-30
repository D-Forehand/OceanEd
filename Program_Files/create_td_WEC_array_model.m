%% CODE DESCRIPTION AND AUTHORS
%
% *************************** Code Description ****************************
% 
% This code derives hydrodynamic time-domain models of arrays of rigid-body
% wave energy converters. It models all the hydrodynamic interactions (i.e.
% it takes into account all the radiated and diffracted waves) between all
% the converters in the array. It takes WAMIT output data as input. It can
% accept input from multiple WAMIT runs for the same array but for
% different sets of angular wave frequencies. At least one of the WAMIT
% runs should include the zero period/infinite frequency case. The code
% reads in the ".out" and ".1" files from the first WAMIT run and it reads
% in the ".1" files from all the other runs. The reason for this is that
% when this code plots the radiation impedance functions against frequency,
% you might see that the frequency resolution is not high enough to capture
% all the variations in these functions. As a result, this code allows you
% to rerun WAMIT for just the missing intermediate frequencies and then
% load the data from the two WAMIT runs, rather than having to rerun WAMIT
% for the complete set of finely spaced frequencies.
%
% The heart of this code is the derivation of the single radiation state-
% space model which approximates the complete set of radiation convolution
% integrals. It is difficult to derive a numerically stable radiation
% state-space model and this code is designed to attempt to do this. More
% details on the theory behind, and implementation of, this code can be
% found in [1].
%
% The output of this code is a Matlab “.mat” file which contains the
% radiation state-space system (sys), the inverse of the reduced mass
% matrix (inv_reduced_mass_matrix), the reduced damping matrix
% (reduced_damping_matrix) and the reduced stiffness matrix
% (reduced_stiffness_matrix). These quantities can then be used in Simulink
% models like the “time_domain_array_sinusoidal_input” model to perform
% time-domain simulations. The matrices are “reduced” because they are only
% for the free (i.e. not the locked) degrees of freedom.
%
% This code uses the function "get_TF.m", which uses Matlab's "invfreqs.m"
% function to find transfer functions which approximate the radiation
% impedance functions.
%
% ******************************** Authors ********************************
%
% Dr David Forehand
% The University of Edinburgh
%
% With help from the late
% Dr Andy McCabe
% Lancaster University
%
% With acknowledgements to:
% Dr Anup Nambiar
% Dr Aristides Kiprakis
% Prof Robin Wallace
% All from The University of Edinburgh
% For their help with the wider wave-to-wire (i.e. from the waves all the
% way to the electricity network) numerical model.
%
% This version - 01/12/2018
%
% *************************************************************************

%% CLEAR THE COMMAND WINDOW, THE WORKSPACE, AND CLOSE ANY OPEN FIGURES

clc % Clear command window
clearvars % Clear variables from memory
clear global % Clear all global variables
close all % Closes any currently open figures

%% INITIALISING CONSTANTS AND PARAMETERS

two_pi=2*pi;

significance_index = 0.01;  % Parameter determining whether a transfer
                            % function is derived.

% Time-domain model fit parameters:

% The following two parameters are used by the Matlab "invfreqs" function:
iter    = 10;       % Maximum number of iterations in numerical search. 
tol     = 0.01;     % (default value) numerical search stops when the norm 
                    % of the gradient is less than 'tol'.
                    
maxN    = 16;       % Highest approximating transfer function order
                    % considered.

%% INPUT THE NUMBER OF WAMIT RUNS

string{01}='Give the number of WAMIT runs  performed  for';
string{02}='the  present  problem.   The  only difference';
string{03}='between these WAMIT runs  must be  the set of';
string{04}='angular wave frequencies analysed:';
string{05}=sprintf('%s\n%s\n%s\n%s\n\n',string{01},string{02}, ...
    string{03},string{04});

num_WAMIT_runs=input(string{05});

clear ('string')

%% INPUT THE NAMES OF THE WAMIT FRC FILES

string{01}='Give  the  filename  of  the  WAMIT  FRC file';
string{02}='number ';
string{03}=', without the ".frc" extension:';
string{04}='WAMIT  filenames   without   their  extension';
string{05}='should  be  less  than  16  characters  long.';
string{06}='Please re-enter the filename:';
string{07}=sprintf('%s\n%s\n%s\n\n',string{04},string{05},string{06});

% Initialising the frc_filename cell array:
frc_filename=cell(num_WAMIT_runs); 

for file_count=1:num_WAMIT_runs
    string{08}=sprintf('\n%s\n%s%d%s\n\n',string{01},string{02}, ...
        file_count,string{03});
    frc_filename{file_count}=input(string{08},'s');
    while length(frc_filename{file_count})>16
        disp(' ');
        frc_filename{file_count}=input(string{07},'s');
    end
end

clear ('file_count','string')

%% INPUT THE OPTIONAL FILE EXTENSION FOR THE OUTPUT FILE

string{01}='Give  the optional extension  for  the output'; 
string{02}='data (.mat) file.  This file will contain the';
string{03}='inverse  of  the  reduced  mass  matrix,  the';
string{04}='reduced damping  and stiffness matrices,  and';
string{05}='the  state-space  system  for  the  radiation';
string{06}='terms.   The name  of  this file  will  be  a';
string{07}='concatenation     of     the    frc_filename,';
string{08}='"_MassDampStiffMatsRadSS", this optional file';
string{09}='extension and ".mat".';
string{10}='(This optional extension  should be less than';
string{11}='16 characters long):';
string{12}=sprintf('\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n', ...
    string{01},string{02},string{03},string{04}, ...
    string{05},string{06},string{07},string{08}, ...
    string{09},string{10},string{11});
string{13}='This optional extension  should be  less than';
string{14}='16 characters long.  Please re-enter it:';
string{15}=sprintf('\n%s\n%s\n\n',string{13},string{14});
optional_extension=input(string{12},'s');
while length(optional_extension)>16
   optional_extension=input(string{15},'s');
end

clear ('string')

%% EXTRACT g, L AND rho FROM THE FIRST WAMIT ".OUT" FILE

fid=fopen(strcat(frc_filename{1},'.out'),'rt');%open the ".out" file
% On a Windows system, the 'rt' opens the file for reading in text mode
% instead of binary mode (the default).
if fid < 0
   error(['Could not open ',frc_filename{1},'.out']);
end

s=fgetl(fid);%read first line
k=strfind(s,'Gravity');%searches s

while isempty(k)==1
    s=fgetl(fid);%reads next line
    k=strfind(s,'Gravity');%searches s
end
g_and_L_data=cell2mat(textscan(s,'%*s %f %*s %*s %f'));
% The asterisks in the above line result in the associated fields being
% ignored (skipped).  See the "Skipping Fields" section of the "textscan"
% page in Matlab's help.
g=g_and_L_data(1);
L=g_and_L_data(2);

s=fgetl(fid);%reads next line
k=strfind(s,'infinite');%searches s
if isempty(k)==1
    rho=cell2mat(textscan(s,'%*s %*s %*f %*s %*s %f'));
else
    rho=cell2mat(textscan(s,'%*s %*s %*s %*s %*s %f'));    
end

rho_g_L2=rho*g*(L^2);
rho_g_L3=rho*g*(L^3);
rho_g_L4=rho*g*(L^4);

fclose(fid);%close the ".out" file

clear('ans','g_and_L_data','fid','k','s')

%% EXTRACT THE NUMBER OF BODIES
% This section extracts the number of bodies from the first WAMIT ".out"
% file.

fid=fopen(strcat(frc_filename{1},'.out'),'rt');%open the ".out" file
if fid < 0
   error(['Could not open ',frc_filename{1},'.out']);
end

s=fgetl(fid);%read first line

num_bodies=0;%set the number of bodies to zero

while ~feof(fid)%while not at the end of the file
    p=strfind(s,'Body number: N=');%searches s
    
    while (isempty(p)==1) && (~feof(fid))%while search empty and not at the
                                         %end of the file
        s=fgetl(fid);%reads next line
        p=strfind(s,'Body number: N=');%searches s
    end
    
    % The above "inner" while loop will exit if we're on a line with a
    % the string 'Body number: N=' OR we have reached the end
    % of the file.
    
    if feof(fid), break, end%if at end of file exit "outer" while loop
    
    num_bodies=num_bodies+1;
    
    s=fgetl(fid);%reads next line
end

% If there is only one body there will be no "Body number: N=" in the 
% ".out" file.  Hence the above code will return num_bodies=0.  The
% following "if" statement deals with that situation.  
if num_bodies==0
    num_bodies=1;
end

fclose(fid);%close the ".out" file

clear('ans','fid','p','s')

%% INITIALISING THE BODY MATRICES
% This section initialises:

% (1) The matrices that will contain the hydrostatic and gravitational
%     restoring coefficients for each boby:
C_H_G=zeros(6,6,num_bodies);

% (2) The body and external mass matrices for each body:
M_B_E=zeros(6,6,num_bodies);

% (3) The external damping matrices for each body:
B_E=zeros(6,6,num_bodies);

% (4) The external stiffness matrices for each body:
C_E=zeros(6,6,num_bodies);

%% EXTRACT THE HYDROSTATIC AND GRAVITATIONAL RESTORING COEFFICIENTS
% This section extracts the hydrostatic and gravitational restoring
% coefficients for each boby from the first WAMIT ".out" file.

fid=fopen(strcat(frc_filename{1},'.out'),'rt');%open the ".out" file
if fid < 0
   error(['Could not open ',frc_filename{1},'.out']);
end

s=fgetl(fid);%read first line

N_body=0;%set the body count to zero

while ~feof(fid)%while not at the end of the file
    p=strfind(s,'Hydrostatic and gravitational');%searches s
    
    while (isempty(p)==1) && (~feof(fid))%while search empty and not at the
                                         %end of the file
        s=fgetl(fid);%reads next line
        p=strfind(s,'Hydrostatic and gravitational');%searches s
    end
    
    % The above "inner" while loop will exit if we're on a line with a
    % the string 'Hydrostatic and gravitational' OR we have reached the end
    % of the file.
    
    if feof(fid), break, end%if at end of file exit "outer" while loop
    
    N_body=N_body+1;
    
    s=fgetl(fid);
    C_H_G(3,3:5,N_body)=cell2mat(textscan(s,'%*s %f %f %f'));
    s=fgetl(fid);
    C_H_G(4,4:6,N_body)=cell2mat(textscan(s,'%*s %f %f %f'));
    s=fgetl(fid);
    C_H_G(5,5:6,N_body)=cell2mat(textscan(s,'%*s %f %f'));
    
    %Dimensionalise the hydrostatic and gravitational restoring
    %coefficients:
    C_H_G(3,3,N_body)=rho_g_L2*C_H_G(3,3,N_body);
    C_H_G(3,4:5,N_body)=rho_g_L3*C_H_G(3,4:5,N_body);
    C_H_G(4,4:6,N_body)=rho_g_L4*C_H_G(4,4:6,N_body);
    C_H_G(5,5:6,N_body)=rho_g_L4*C_H_G(5,5:6,N_body);
    
    %There are 3 additional hydrostatic restoring coefficients which can be
    %non-zero (see the text below the equations for the 8 hydrostatic and
    %gravitational restoring coefficients on page 4-2 of [2]):
    C_H_G(4,3,N_body)=C_H_G(3,4,N_body);
    C_H_G(5,3,N_body)=C_H_G(3,5,N_body);
    C_H_G(5,4,N_body)=C_H_G(4,5,N_body);
    
end

fclose(fid);%close the ".out" file

clear('N_body','ans','fid','p','s')

%% EXTRACT THE MASS MATRICES FROM THE WAMIT ".OUT" FILE
% This section extracts the body and external mass matrices for each body
% from the first WAMIT ".out" file.  Note these matrices are already 
% dimensional.

fid=fopen(strcat(frc_filename{1},'.out'),'rt');%open the ".out" file
if fid < 0
   error(['Could not open ',frc_filename{1},'.out']);
end

s=fgetl(fid);%read first line

N_body=0;%set the body count to zero

while ~feof(fid)%while not at the end of the file
    p=strfind(s,'external mass matrix:');%searches s
    
    while (isempty(p)==1) && (~feof(fid))%while search empty and not at the
                                         %end of the file
        s=fgetl(fid);%reads next line
        p=strfind(s,'external mass matrix:');%searches s
    end
    
    if feof(fid), break, end%if at end of file exit while loop
    
    N_body=N_body+1;
    
    s=fgetl(fid);
    M_B_E(1,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
    s=fgetl(fid);
    M_B_E(2,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
    s=fgetl(fid);
    M_B_E(3,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
    s=fgetl(fid);
    M_B_E(4,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
    s=fgetl(fid);
    M_B_E(5,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
    s=fgetl(fid);
    M_B_E(6,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
    
end

fclose(fid);%close the ".out" file

clear('N_body','ans','fid','p','s')

%% EXTRACT THE EXTERNAL DAMPING MATRICES FROM THE WAMIT ".OUT" FILE
% This section extracts the external damping matrices for each body from
% the first WAMIT ".out" file.  Note these matrices are already
% dimensional.

fid=fopen(strcat(frc_filename{1},'.out'),'rt');%open the ".out" file
if fid < 0
   error(['Could not open ',frc_filename{1},'.out']);
end

s=fgetl(fid);%read first line

N_body=0;%set the body count to zero

while ~feof(fid)%while not at the end of the file
    p1=strfind(s,'Body number: N=');%searches s
    p2=strfind(s,'external damping matrix');%searches s
    
    while (isempty(p1)==1) && (isempty(p2)==1) && (~feof(fid))%while search
                                                              %empty and
                                                              %not at the
                                                              %end of the
                                                              %file
        s=fgetl(fid);%reads next line
        p1=strfind(s,'Body number: N=');%searches s
        p2=strfind(s,'external damping matrix');%searches s
    end
    
    if feof(fid), break, end%if at end of file exit while loop
    
    if isempty(p1)==0
        N_body=N_body+1;
    end
    
    if isempty(p2)==0
        s=fgetl(fid);
        B_E(1,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
        s=fgetl(fid);
        B_E(2,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
        s=fgetl(fid);
        B_E(3,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
        s=fgetl(fid);
        B_E(4,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
        s=fgetl(fid);
        B_E(5,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
        s=fgetl(fid);
        B_E(6,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
    end
    
    s=fgetl(fid);%reads next line
end

fclose(fid);%close the ".out" file

clear('N_body','ans','fid','p1','p2','s')

%% EXTRACT THE EXTERNAL STIFFNESS MATRICES FROM THE WAMIT ".OUT" FILE
% This section extracts the external stiffness matrices for each body from
% the first WAMIT ".out" file.  Note these matrices are already
% dimensional.

fid=fopen(strcat(frc_filename{1},'.out'),'rt');%open the ".out" file
if fid < 0
   error(['Could not open ',frc_filename{1},'.out']);
end

s=fgetl(fid);%read first line

N_body=0;%set the body count to zero

while ~feof(fid)%while not at the end of the file
    p1=strfind(s,'Body number: N=');%searches s
    p2=strfind(s,'external stiffness matrix');%searches s
    
    while (isempty(p1)==1) && (isempty(p2)==1) && (~feof(fid))%while search
                                                              %empty and
                                                              %not at the
                                                              %end of the
                                                              %file
        s=fgetl(fid);%reads next line
        p1=strfind(s,'Body number: N=');%searches s
        p2=strfind(s,'external stiffness matrix');%searches s
    end
    
    if feof(fid), break, end%if at end of file exit while loop
    
    if isempty(p1)==0
        N_body=N_body+1;
    end
    
    if isempty(p2)==0
        s=fgetl(fid);
        C_E(1,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
        s=fgetl(fid);
        C_E(2,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
        s=fgetl(fid);
        C_E(3,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
        s=fgetl(fid);
        C_E(4,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
        s=fgetl(fid);
        C_E(5,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
        s=fgetl(fid);
        C_E(6,:,N_body)=cell2mat(textscan(s,'%f %f %f %f %f %f'));
    end
    
    s=fgetl(fid);%reads next line
end

fclose(fid);%close the ".out" file

clear('N_body','ans','fid','p1','p2','s')

%% LOAD ADDED MASS AND ADDED DAMPING DATA
% In addition to a number of finite positive frequencies that the WAMIT run
% will analyse, we are assuming it will also analyse the infinite frequency
% (zero period) case.  This is obtained by placing 0.0 in the PER array in
% the WAMIT POT input file (see page 4-9 of [2]).  If this is the case, the
% length of the rows in the ".1" output file will vary between 4 and 5
% columns.  An "infinite frequency (zero period)" row will just have 4
% columns because the added damping at infinite frequency is zero and is,
% consequently, not given in this row.  A "finite frequency (finite,
% non-zero period)" row will have 5 columns because both the added mass and
% added damping are given.
%
% Due to these different row lengths in the ".1" WAMIT output file, it
% cannot be opened with the command:
%
% hydro_data=load(strcat(frc_filename,'.1'));
%
% Instead the code below reads the ".1" WAMIT output file and creates two
% matrices.  The inf_data matrix contains the data corresponding to the
% infinite frequency case and the non_inf_data matrix contains the data
% corresponding to the finite frequency case.

% This first section counts how many rows there should be in inf_data and
% non_inf_data so that these two arrays can be initialised:

non_inf_count=1;%initialising row index for non_inf_data matrix
max_inf_count=1;

for file_count=1:num_WAMIT_runs
    
    fid=fopen(strcat(frc_filename{file_count},'.1'),'rt');%open the ".1"
                                                          %file
    if fid < 0
        error(['Could not open ',frc_filename{file_count},'.1']);
    end

    s=fgetl(fid);%read first line

    inf_count=1;%initialising row index for inf_data matrix
    
    while s~=-1%while not at the end of the file
        
        temp_s=s;
        
        %Counting the number of columns in row "s":
        row_count=1;
        [token,remain]=strtok(temp_s);
        while isempty(remain)==0
            temp_s=remain;
            [token,remain]=strtok(temp_s);
            row_count=row_count+1;
        end

        if row_count==4
            inf_count=inf_count+1;
        elseif row_count==5
            non_inf_count=non_inf_count+1;
        else
            error(['Wrong number of rows in',frc_filename{file_count}, ...
                  '.1 file']);
        end

        s=fgetl(fid);%read the next line
    
    end

    fclose(fid);%close the ".1" file

    if inf_count>1
        max_inf_count=inf_count;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This second section initialises the inf_data and non_inf_data arrays:

inf_data=zeros(max_inf_count-1,4);
non_inf_data=zeros(non_inf_count-1,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This third section populates the inf_data and non_inf_data arrays:

non_inf_count=1;%initialising row index for non_inf_data matrix

for file_count=1:num_WAMIT_runs

    fid=fopen(strcat(frc_filename{file_count},'.1'),'rt');%open the ".1"
                                                          %file
    if fid < 0
    error(['Could not open ',frc_filename{file_count},'.1']);
    end

    s=fgetl(fid);%read first line

    inf_count=1;%initialising row index for inf_data matrix

    while s~=-1%while not at the end of the file
        
        temp_s=s;
        
        %Counting the number of columns in row "s":
        row_count=1;
        [token,remain]=strtok(temp_s);
        while isempty(remain)==0
            temp_s=remain;
            [token,remain]=strtok(temp_s);
            row_count=row_count+1;
        end

        if row_count==4
            inf_data(inf_count,:)=cell2mat(textscan(s,'%f %f %f %f'));
            inf_count=inf_count+1;
        elseif row_count==5
            non_inf_data(non_inf_count,:)= ...
                                    cell2mat(textscan(s,'%f %f %f %f %f'));
            non_inf_count=non_inf_count+1;
        else
            error(['Wrong number of rows in',frc_filename{file_count}, ...
                  '.1 file']);
        end

        s=fgetl(fid);%read the next line
    end

    fclose(fid);%close the ".1" file

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear('ans','fid','file_count','inf_count','max_inf_count', ...
      'non_inf_count','remain','row_count','s','temp_s','token')

%% CALCULATING THE NUMBER OF DIFFERENT I,J MODES
% By the end of this cell, this value will be stored in the variable
% ij_mode_count.

first_mode_i=non_inf_data(1,2);
first_mode_j=non_inf_data(1,3);
ij_mode_count=1;
mode_i=non_inf_data(2,2);
mode_j=non_inf_data(2,3);
while (mode_i ~= first_mode_i) || (mode_j ~= first_mode_j) % Note the use 
                                                           % of || for "OR"
                                                           % here.
    ij_mode_count=ij_mode_count+1;
    mode_i=non_inf_data(ij_mode_count+1,2);
    mode_j=non_inf_data(ij_mode_count+1,3);
end
% To display the number of different i,j modes on the screen, uncomment the
% next line.
% disp(['Number of i,j modes = ' num2str(ij_mode_count)])

if ij_mode_count ~= size(inf_data,1)
   error('The number of rows of inf_data does not equal ij_mode_count');
end

clear('first_mode_i','first_mode_j','mode_i','mode_j')

%% CONVERTING THE FIRST COLUMN OF NON_INF_DATA FROM PERIOD TO FREQUENCY

non_inf_data(:,1)=two_pi./non_inf_data(:,1);

%% ODTAINING INDEX_1 WHICH CONTAINS THE ORDER OF THE FREQUENCIES IN
% NON_INF_DATA

[temp,index_1]=sort(non_inf_data(1:ij_mode_count:end,1));

clear('temp')

%% CALCULATING THE NUMBER OF FREQUENCIES

num_freqs=size(non_inf_data,1)/ij_mode_count;

if num_freqs ~= length(index_1)
   error('The number of frequencies does not equal the lenght of index_1');
end

%% CALCULATING INDEX_2 WHICH WILL BE USED TO SORT NON_INF_DATA INTO
% ASCENDING FREQUENCY ORDER

for freq_count=1:num_freqs
    index_2((freq_count-1)*ij_mode_count+1:freq_count*ij_mode_count)= ...
        (index_1(freq_count)-1)*ij_mode_count+1: ...
        index_1(freq_count)*ij_mode_count;
end

clear('freq_count','index_1')

%% SORTING NON_INF_DATA INTO ASCENDING FREQUENCY ORDER

non_inf_data=non_inf_data(index_2,:);

clear('index_2')

%% LOOP TO DIMENSIONALISE COLUMN 4 OF INF_DATA FOR FUTURE USE:

for ij_mode=1:ij_mode_count
    i_mode=non_inf_data(ij_mode,2);
    j_mode=non_inf_data(ij_mode,3);
    
    % CALCULATING k (for the definition of k see the top of page 3-4 in
    % [2]). For the example vector y below, see the difference between
    % mod(y,6) and y-(ceil(y/6)*6)+6
    %
    %                 y = [1  2  3  4  5  6  7  8  9 10 11 12 13 14 15]
    %          mod(y,6) = [1  2  3  4  5  0  1  2  3  4  5  0  1  2  3]
    % y-(ceil(y/6)*6)+6 = [1  2  3  4  5  6  1  2  3  4  5  6  1  2  3] 
    
    i_between_1_and_6=i_mode-ceil(i_mode/6)*6+6;
    j_between_1_and_6=j_mode-ceil(j_mode/6)*6+6;
    k=3+floor(i_between_1_and_6/3.5)+floor(j_between_1_and_6/3.5);
    rho_times_L_to_the_k=rho*(L^k);
    
     % Dimensionalising column 4 of inf_data for future use:
    inf_data(ij_mode,4)=rho_times_L_to_the_k*inf_data(ij_mode,4);

end

clear('ij_mode','i_mode','j_mode','i_between_1_and_6', ...
      'j_between_1_and_6','k','rho_times_L_to_the_k')

%% CALCULATING THE NUMBER OF DOFS AND THE VARIOUS MODES PRESENT ETC
% Suppose we have only two modes present: heave of body 1 (i.e. mode 3) and
% heave of body 2 (i.e. mode 9).  Then 'mode' will be the row vector [3 9]
% and 'mode_index' will be the "inverse" of this.  That is, mode_index will
% be the row vector [0 0 1 0 0 0 0 0 2], i.e. mode_index(3) = 1 and
% mode_index(9) = 2.

% Calculating the number of degrees of freedom (number of modes):
num_dofs=sqrt(size(inf_data,1));

% Calculating the 'mode' row vector:
mode=(inf_data(1:num_dofs:end,2))';

% Initialising the 'mode_index' row vector:
mode_index=zeros(1,max(mode));

% Calculating the 'mode_index' row vector:
mode_index(mode)=1:num_dofs;

%% CREATING A STRUCTURE FOR THE PERIOD, FREQUENCY, ADDED MASS AND DAMPING

freaknumber = 1;
for n1 = 1:size(non_inf_data,1)
    if (n1>1) && (non_inf_data(n1,1)~=non_inf_data((n1-1),1))
        freaknumber = freaknumber + 1;
    end
    row = non_inf_data(n1,2);
    col = non_inf_data(n1,3);
    
    row_between_1_and_6=row-ceil(row/6)*6+6;
    col_between_1_and_6=col-ceil(col/6)*6+6;
    k=3+floor(row_between_1_and_6/3.5)+floor(col_between_1_and_6/3.5);
    
    row=mode_index(row);
    col=mode_index(col);
    
    ABdata(row,col).wavefrq(freaknumber,1) = non_inf_data(n1,1);

    ABdata(row,col).addedmass(freaknumber,1) = ...
                                              rho*(L^k)*non_inf_data(n1,4);
    ABdata(row,col).damping(freaknumber,1) = ...
       rho*(L^k)*non_inf_data(n1,5)*ABdata(row,col).wavefrq(freaknumber,1);
end

clear('freaknumber','n1','row','col','row_between_1_and_6', ...
      'col_between_1_and_6','k')

%% CREATING A 2D MATRIX WITH THE INFINITE-FREQUENCY ADDED-MASS DATA

Ainfs=zeros(num_dofs,num_dofs);

for n1 = 1:size(inf_data,1)
    row = mode_index(inf_data(n1,2));
    col = mode_index(inf_data(n1,3));
    
    Ainfs(row,col) = inf_data(n1,4);
end

clear('n1','row','col')

%% ADDING INFORMATION ABOUT THE REAL AND IMAGINARY PARTS OF THE RADIATION
% IMPEDANCE FUNCTIONS TO THE ABDATA STRUCTURE

for i_mode = 1:num_dofs
    for j_mode = 1:num_dofs
        ABdata(i_mode,j_mode).imaginarybit = ...
            (ABdata(i_mode,j_mode).addedmass ...
            - Ainfs(i_mode,j_mode)).*ABdata(i_mode,j_mode).wavefrq;
        ABdata(i_mode,j_mode).realbit      = ...
            ABdata(i_mode,j_mode).damping;
    end
end

clear('i_mode','j_mode')

%% PLOTTING THE REAL AND IMAGINARY PARTS OF THE RADIATION IMPEDANCE FUNCS

string{01}='Do you  want to see  plots of  the  real  and';
string{02}='imaginary parts  of  the  radiation impedance';
string{03}='functions?';
string{04}='Type "y" for yes or "n" for no:';
string{05}=sprintf('\n%s\n%s\n%s\n%s\n\n',string{01},string{02}, ...
    string{03},string{04});
string{06}='Press return to plot the next function:';
string{07}=sprintf('%s',string{06});
string{08}='Last function plotted, press return to exit:';
string{09}=sprintf('%s',string{08});

max_omega=max(ABdata(1,1).wavefrq);

plot_or_not=input(string{05},'s');
if (strcmp(plot_or_not,'y')) || (strcmp(plot_or_not,'"y"') || ...
        strcmp(plot_or_not,'Y')) || (strcmp(plot_or_not,'"Y"'))
    for i_mode=1:num_dofs
        for j_mode=i_mode:num_dofs
            
            i_between_1_and_6=mode(i_mode)-ceil(mode(i_mode)/6)*6+6;
            j_between_1_and_6=mode(j_mode)-ceil(mode(j_mode)/6)*6+6;
            k=3+floor(i_between_1_and_6/3.5)+floor(j_between_1_and_6/3.5);
            
            if i_mode == j_mode
                plot([0; ABdata(i_mode,j_mode).wavefrq], ...
                     [0; ABdata(i_mode,j_mode).realbit])
                title({'Real Part of the Radiation Impedance Function K_{ij}(\omega),'; ...
                       ['with i=' num2str(mode(i_mode)) ' and j=' num2str(mode(j_mode)) ', against \omega']});
                axis([0 max_omega -1 1]); axis('auto y');
                xlabel('\omega (rad/s)');
                if k == 3
                    ylabel('Re(K_{ij}(\omega)) (kg/s)');
                elseif k == 4
                    ylabel('Re(K_{ij}(\omega)) (kg m/s)');
                elseif k == 5
                    ylabel('Re(K_{ij}(\omega)) (kg m^{2}/s)');
                end
                next_plot=input(string{07},'s');
                plot([0; ABdata(i_mode,j_mode).wavefrq], ...
                     [0; ABdata(i_mode,j_mode).imaginarybit])
                title({'Imaginary Part of the Radiation Impedance Function K_{ij}(\omega),'; ...
                       ['with i=' num2str(mode(i_mode)) ' and j=' num2str(mode(j_mode)) ', against \omega']});
                axis([0 max_omega -1 1]); axis('auto y');
                xlabel('\omega (rad/s)');
                if k == 3
                    ylabel('Im(K_{ij}(\omega)) (kg/s)');
                elseif k == 4
                    ylabel('Im(K_{ij}(\omega)) (kg m/s)');
                elseif k == 5
                    ylabel('Im(K_{ij}(\omega)) (kg m^{2}/s)');
                end
                if i_mode == num_dofs
                    next_plot=input(string{09},'s');
                else
                    next_plot=input(string{07},'s');
                end
            else
                plot([0; ABdata(i_mode,j_mode).wavefrq], ...
                     [0; ABdata(i_mode,j_mode).realbit], ...
                     [0; ABdata(j_mode,i_mode).wavefrq], ...
                     [0; ABdata(j_mode,i_mode).realbit])
                title({'Real Part of the Radiation Impedance Function K_{ij}(\omega) and K_{ji}(\omega),'; ...
                       ['with i=' num2str(mode(i_mode)) ' and j=' num2str(mode(j_mode)) ', against \omega']});
                axis([0 max_omega -1 1]); axis('auto y');
                xlabel('\omega (rad/s)');
                if k == 3
                    ylabel('Re(K_{ij}(\omega)) and Re(K_{ji}(\omega)) (kg/s)');
                elseif k == 4
                    ylabel('Re(K_{ij}(\omega)) and Re(K_{ji}(\omega)) (kg m/s)');
                elseif k == 5
                    ylabel('Re(K_{ij}(\omega)) and Re(K_{ji}(\omega)) (kg m^{2}/s)');
                end
                legend(['i=' num2str(mode(i_mode)) ' and j=' num2str(mode(j_mode))], ...
                       ['i=' num2str(mode(j_mode)) ' and j=' num2str(mode(i_mode))]);
                next_plot=input(string{07},'s');
                plot([0; ABdata(i_mode,j_mode).wavefrq], ...
                     [0; ABdata(i_mode,j_mode).imaginarybit], ...
                     [0; ABdata(j_mode,i_mode).wavefrq], ...
                     [0; ABdata(j_mode,i_mode).imaginarybit])
                title({'Imaginary Part of the Radiation Impedance Function K_{ij}(\omega) and K_{ji}(\omega),'; ...
                       ['with i=' num2str(mode(i_mode)) ' and j=' num2str(mode(j_mode)) ', against \omega']});
                axis([0 max_omega -1 1]); axis('auto y');
                xlabel('\omega (rad/s)');
                if k == 3
                    ylabel('Im(K_{ij}(\omega)) and Im(K_{ji}(\omega)) (kg/s)');
                elseif k == 4
                    ylabel('Im(K_{ij}(\omega)) and Im(K_{ji}(\omega)) (kg m/s)');
                elseif k == 5
                    ylabel('Im(K_{ij}(\omega)) and Im(K_{ji}(\omega)) (kg m^{2}/s)');
                end
                legend(['i=' num2str(mode(i_mode)) ' and j=' num2str(mode(j_mode))], ...
                       ['i=' num2str(mode(j_mode)) ' and j=' num2str(mode(i_mode))]);
                next_plot=input(string{07},'s');
            end
        end
    end
end

clear('string','plot_or_not','i_mode','j_mode','i_between_1_and_6', ...
      'j_between_1_and_6','k','next_plot')

%% CREATE 'SIGNIFICANCE MATRIX' - 1 IF TF IS DERIVED, 0 IF TF NOT DERIVED
% Note that modes are in blocks of 3 (translational/rotational/.... etc.)
% blocknumber = ceil(size(ABdata,1)/3);
blocknumber = ceil(max(mode)/3);

% Mean magnitude of data
meanmag=zeros(3*blocknumber,3*blocknumber);
for i_mode = 1:num_dofs
    for j_mode = 1:num_dofs
        meanmag(mode(i_mode),mode(j_mode))             = ...
            mean(sqrt((ABdata(i_mode,j_mode).imaginarybit.^2) ...
            + (ABdata(i_mode,j_mode).realbit.^2)));
    end
end

clear('i_mode','j_mode')

% Check if mean magnitude of data is more than the product of maximum
% mean magnitude in block and significance_index....
Asubmat1=zeros(3,3);
sigmat=zeros(size(meanmag));
for n2 = 1:blocknumber
    for n3 = 1:blocknumber
        rows = (1:3) + ((n2-1)*3);
        cols = (1:3) + ((n3-1)*3);
        Asubmat = meanmag(rows,cols);
        Asmax   = max(max(Asubmat));
        
        for n4 = 1:3
            for n5 = 1:3
                if Asmax == 0
                    Asubmat1(n4,n5) = 0;
                elseif Asubmat(n4,n5)/Asmax < significance_index
                    Asubmat1(n4,n5) = 0;
                else
                    Asubmat1(n4,n5) = 1;
                end
            end
        end
        sigmat(rows,cols) = Asubmat1;        
    end
end

clear('Asubmat1','n2','n3','rows','cols','Asubmat','Asmax','n4','n5', ...
    'blocknumber','meanmag')

%% FIND ALL THE APPROXIMATING TRANSFER FUNCTIONS FOR ALL THE RADIATION
% IMPEDANCE FUNCTIONS

% This is a parameter of Andy's get_TF function which is not used in this
% code:
fitdatalength = 1;

% **** TRANSFER FUNCTIONS ARE ONLY FITTED UP TO THIS OMEGA (FREQUENCY) ****
% ********************** THE USER CAN EDIT THIS PART **********************
% 
% UNCOMMENT THE NEXT LINE TO FIT THE TRANSFER FUNCTIONS TO ALL FREQUENCIES: 
% last_omega_index = num_freqs;
%
% ALTERNATIVELY,  COMMENT OUT THE ABOVE LINE OF CODE AND UNCOMMENT AND EDIT
% THE NEXT LINE TO FIT TRANSFER FUNCTIONS JUST UP TO A CERTAIN FREQUENCY:
last_omega_index = 66;
%
% *************************************************************************

datascaler_in = 1; % This just scales the input radiation impedance
                   % functions
datascaler_out = 1/datascaler_in; % This just scales the approximating
                                  % transfer functions back again
fit_errormax = 0.01; % This is the maximum relative root-mean-square error
                     % between the transfer function and radiation
                     % impedance function that it is approximating.

% Create a cell array which will contain the numerator coefficients of all
% the approximating transfer functions:                     
tfnums=cell(num_dofs,num_dofs);

% Create a cell array which will contain the denominator coefficients of
% all the approximating transfer functions:
tfdens=cell(num_dofs,num_dofs);

% Create cell arrays which will contain the zeros and poles of all the
% approximating transfer functions:
zs=cell(num_dofs,num_dofs);
ps=cell(num_dofs,num_dofs);

% Create a numeric array to store the gains of all the approximating
% transfer functions:
ks=zeros(num_dofs,num_dofs);
            
not_all_transfer_funcs_found = 0;
            
for i_mode = 1:num_dofs
                
    if not_all_transfer_funcs_found
        break
    end
                
    for j_mode = 1:num_dofs
                    
        clc
        fprintf('\n');
        fprintf('Last Omega Index: \n\n');
        fprintf('%4d \n\n',last_omega_index);
        fprintf('RADIATION TRANSFER FUNCTIONS: \n\n');
        fprintf('mode = %2d %2d \n\n',mode(i_mode),mode(j_mode));
                    
        transfer_func_found = 0;
        
        % If transfer function is to be derived:
        if sigmat(mode(i_mode),mode(j_mode)) ~= 0
                        
            modeflag=[i_mode j_mode 0];
            
            % For all approximating transfer function numerator orders up
            % to maxN and all denominator orders up to the numerator order
            % minus one, search until a transfer function is found which is
            % stable and which has a relative root-mean-square error to the
            % radiation impedance function of less than fit_errormax:
            for a_order = 2:maxN % denominator order
                            
                if transfer_func_found
                    break
                end
                            
                for b_order = 1:a_order-1 % numerator order
                    
                    % Construct the radiation impedance function:
                    h = ABdata(i_mode,j_mode).realbit(1:last_omega_index) ...
                        + (1i*ABdata(i_mode,j_mode).imaginarybit(1:last_omega_index));
                    
                    % Scale the radiation impedance function:
                    h = datascaler_in*h;
                                
                    % Compute the approximating transfer function with
                    % numerator order b_order and denominator order
                    % a_order. Just jump (continue) to next transfer
                    % function order if no transfer function is found:
                    try
                        [TFnum,TFden,fitdata,no_TFs] = ...
                            get_TF(modeflag,fitdatalength,h,ABdata(i_mode,j_mode).wavefrq(1:last_omega_index),b_order,a_order,iter,tol);
                    catch ME
                        continue
                    end
                    
                    % Apply the opposite scale to the approximating
                    % transfer function:
                    TFnum = datascaler_out*TFnum;
                    
                    % Check to see if approxinating transfer function is
                    % stable and 
                    if (max(real(roots(TFden)))<=-1e-3) && ...
                            (fitdata(3)<=fit_errormax)
                        transfer_func_found = 1;
                        break
                    end
                                
                end

            end

        else
            % Set transfer function to zero if it was not to be derived
            % according to sigmat:
            TFnum   = 0;
            TFden   = 1;
            transfer_func_found = 1;

        end
                    
        if ~transfer_func_found
            not_all_transfer_funcs_found = 1;
            break
        end
        
        tfnums(i_mode,j_mode) = {TFnum};
        tfdens(i_mode,j_mode) = {TFden};
                    
        % Convert transfer function to zero-pole-gain form:
        [z0,p0,k0] = tf2zp(TFnum,TFden);
        zs(i_mode,j_mode) = {z0};   % Zeros and poles into the cell arrays
        ps(i_mode,j_mode) = {p0};
        ks(i_mode,j_mode) = k0;     % Gain data into its numeric array
        
        clear TFnum TFden fitdata z0 p0 k0
    end
end

%% CONVERT THE SYSTEM OF ZERO-POLE-GAIN FORM TRANSFER FUNCTIONS TO A SINGLE
% RADIATION STATE-SPACE MODEL

% If all the transfer functions were found:
if ~not_all_transfer_funcs_found
    % Convert all the transfer functions in zero-pole-gain format into one
    % zero-pole-gain model:
    H = zpk(zs,ps,ks);
                
    % Then convert to a single raditation state-space model:
    % The "try" construct is to catch any errors if, say,
    % memory runs out during this process. 
    try
        sys = ss(H);
    catch ME
        return % This "return" will hopefully exit the containing "if"
               % statement without executing the remaining code inside it.
               % I've tested example code and it seems to work.
    end
                
    % Extract the system matrix for simulink block:
    Arad = get(sys,'a');
                
    % Test for stability:
    eigenvalues = eig(Arad);
    max_real_e_value = max(real(eigenvalues));
    fprintf('\n\n maximum real part of eigenvalues: %0.3f \n\n',max_real_e_value);
                
    % Nice plot:
    plot(real(eigenvalues),imag(eigenvalues),'ok')
    xlabel('REAL')
    ylabel('IMAGINARY')
    title('EIGENVALUES')
    grid
                
    % Size and condition number of system matrix:
    size_Arad = size(Arad,1);
    cond_num = cond(Arad);                
end

%% PLOTTING THE REAL AND IMAGINARY PARTS OF BOTH K_ij(omega) AND THE
%  FREQUENCY RESPONSE H(omega) OF ITS APPROXIMATING TRANSFER FUNCTION
%  AGAINST OMEGA.

string{01}='Do you  want to see  plots of  the  real  and';
string{02}='imaginary  parts   of   both   the  radiation';
string{03}='impedance   functions   and   the   frequency';
string{04}='response   of  their  approximating  transfer';
string{05}='functions?';
string{06}='Type "y" for yes or "n" for no:';
string{07}=sprintf('\n%s\n%s\n%s\n%s\n%s\n%s\n\n',string{01}, ...
    string{02},string{03},string{04},string{05},string{06});
string{08}='Press return to plot the next function:';
string{09}=sprintf('%s',string{08});
string{10}='Last function plotted, press return to exit:';
string{11}=sprintf('%s',string{10});

max_omega=max(ABdata(1,1).wavefrq);

plot_or_not=input(string{07},'s');

if (strcmp(plot_or_not,'y')) || (strcmp(plot_or_not,'"y"') || ...
        strcmp(plot_or_not,'Y')) || (strcmp(plot_or_not,'"Y"'))
    for i_mode=1:num_dofs
        for j_mode=1:num_dofs
            
            i_between_1_and_6=mode(i_mode)-ceil(mode(i_mode)/6)*6+6;
            j_between_1_and_6=mode(j_mode)-ceil(mode(j_mode)/6)*6+6;
            k=3+floor(i_between_1_and_6/3.5)+floor(j_between_1_and_6/3.5);
            
            tf_sys = tf(tfnums(i_mode,j_mode),tfdens(i_mode,j_mode));
            f1 = freqresp(tf_sys,ABdata(i_mode,j_mode).wavefrq);
            fr = squeeze(f1);
            
            % For the two lines below see the Matlab help for "invfreqs":
            num_order = length(tfnums{i_mode,j_mode})-1;
            den_order = length(tfdens{i_mode,j_mode})-1;
                
            plot([0; ABdata(i_mode,j_mode).wavefrq], ...
                 [0; ABdata(i_mode,j_mode).realbit], ...
                 [0; ABdata(i_mode,j_mode).wavefrq], ...
                 [0; real(fr)])
            title({'Real Part of both K_{ij}(\omega) and the Frequency Response H(\omega) of its', ...
                   ['Approximating Transfer Function against \omega for ' ...
                    'i=' num2str(mode(i_mode)) ' and j=' num2str(mode(j_mode))]});
            axis([0 max_omega -1 1]); axis('auto y');
            xlabel('\omega (rad/s)');
            if k == 3
                ylabel('Real Part of Freq Response (kg/s)');
            elseif k == 4
                ylabel('Real Part of Freq Response (kg m/s)');
            elseif k == 5
                ylabel('Real Part of Freq Response (kg m^{2}/s)');
            end
            legend('Re(K_{ij}(\omega))', ...
                   'Re(H(\omega))');
            text(0.05,0.94, ...
                 ['Numer order = ' num2str(num_order)], ...
                 'Units','normalized', ...
                 'BackgroundColor',[1 1 0]);
            text(0.05,0.88, ...
                 ['Denom order = ' num2str(den_order)], ...
                 'Units','normalized', ...
                 'BackgroundColor',[1 1 0]);
                 
            next_plot=input(string{09},'s');
                
            plot([0; ABdata(i_mode,j_mode).wavefrq], ...
                 [0; ABdata(i_mode,j_mode).imaginarybit], ...
                 [0; ABdata(i_mode,j_mode).wavefrq], ...
                 [0; imag(fr)])
            title({'Imaginary Part of both K_{ij}(\omega) and the Frequency Response H(\omega) of its', ...
                   ['Approximating Transfer Function against \omega for ' ...
                    'i=' num2str(mode(i_mode)) ' and j=' num2str(mode(j_mode))]});                
            axis([0 max_omega -1 1]); axis('auto y');
            xlabel('\omega (rad/s)');
            if k == 3
                ylabel('Imag Part of Freq Response (kg/s)');
            elseif k == 4
                ylabel('Imag Part of Freq Response (kg m/s)');
            elseif k == 5
                ylabel('Imag Part of Freq Response (kg m^{2}/s)');
            end
            legend('Im(K_{ij}(\omega))', ...
                   'Im(H(\omega))');
            text(0.05,0.94, ...
                 ['Numer order = ' num2str(num_order)], ...
                 'Units','normalized', ...
                 'BackgroundColor',[1 1 0]);
            text(0.05,0.88, ...
                 ['Denom order = ' num2str(den_order)], ...
                 'Units','normalized', ...
                 'BackgroundColor',[1 1 0]);
               
            if (i_mode == num_dofs && j_mode == num_dofs)
                next_plot=input(string{11},'s');
            else
                next_plot=input(string{09},'s');
            end

        end
    end
end

clear('string','plot_or_not','i_mode','j_mode','i_between_1_and_6', ...
      'j_between_1_and_6','k','next_plot')

%% CONSTRUCTING THE REDUCED MASS MATRIX
% Here we will use the following constructs:
%                 y = [1  2  3  4  5  6  7  8  9 10 11 12 13 14 15]
%         ceil(y/6) = [1  1  1  1  1  1  2  2  2  2  2  2  3  3  3]
% y-(ceil(y/6)*6)+6 = [1  2  3  4  5  6  1  2  3  4  5  6  1  2  3]

reduced_mass_matrix=zeros(num_dofs,num_dofs);%Preallocating the array
for i_index=1:num_dofs
    i_body=ceil(mode(i_index)/6);
    i_body_mode=mode(i_index)-i_body*6+6;
    for j_index=1:num_dofs
        j_body=ceil(mode(j_index)/6);
        j_body_mode=mode(j_index)-j_body*6+6;
        if j_body==i_body
            reduced_mass_matrix(i_index,j_index)= ...
                M_B_E(i_body_mode,j_body_mode,i_body) ...
                +inf_data(num_dofs*(i_index-1)+j_index,4);
        else
            reduced_mass_matrix(i_index,j_index)= ...
                inf_data(num_dofs*(i_index-1)+j_index,4);
        end
    end
end

clear('i_index','i_body','i_body_mode','j_index','j_body','j_body_mode')

%% CONSTRUCTING THE INVERSE OF THE REDUCED MASS MATRIX

inv_reduced_mass_matrix=inv(reduced_mass_matrix);

clear('reduced_mass_matrix')

%% CONSTRUCTING THE REDUCED DAMPING MATRIX
% Here we will use the following constructs:
%                 y = [1  2  3  4  5  6  7  8  9 10 11 12 13 14 15]
%         ceil(y/6) = [1  1  1  1  1  1  2  2  2  2  2  2  3  3  3]
% y-(ceil(y/6)*6)+6 = [1  2  3  4  5  6  1  2  3  4  5  6  1  2  3]

reduced_damping_matrix=zeros(num_dofs,num_dofs);%Preallocating the array
for i_index=1:num_dofs
    i_body=ceil(mode(i_index)/6);
    i_body_mode=mode(i_index)-i_body*6+6;
    for j_index=1:num_dofs
        j_body=ceil(mode(j_index)/6);
        j_body_mode=mode(j_index)-j_body*6+6;
        if j_body==i_body
            reduced_damping_matrix(i_index,j_index)= ...
                B_E(i_body_mode,j_body_mode,i_body);
        end
    end
end

clear('i_index','i_body','i_body_mode','j_index','j_body','j_body_mode')

%% CONSTRUCTING THE REDUCED STIFFNESS MATRIX
% Here we will use the following constructs:
%                 y = [1  2  3  4  5  6  7  8  9 10 11 12 13 14 15]
%         ceil(y/6) = [1  1  1  1  1  1  2  2  2  2  2  2  3  3  3]
% y-(ceil(y/6)*6)+6 = [1  2  3  4  5  6  1  2  3  4  5  6  1  2  3]

reduced_stiffness_matrix=zeros(num_dofs,num_dofs);%Preallocating the array
    for i_index=1:num_dofs
        i_body=ceil(mode(i_index)/6);
        i_body_mode=mode(i_index)-i_body*6+6;
        for j_index=1:num_dofs
            j_body=ceil(mode(j_index)/6);
            j_body_mode=mode(j_index)-j_body*6+6;
            if j_body==i_body
                reduced_stiffness_matrix(i_index,j_index)= ...
                    C_H_G(i_body_mode,j_body_mode,i_body) ...
                    +C_E(i_body_mode,j_body_mode,i_body);
            end
        end
    end

%% SAVING THE DATA FOR THE TIME-DOMAIN HYDRODYNAMIC ARRAY MODEL

Matrices_and_SS_file=strcat(frc_filename{1},'_MassDampStiffMatsRadSS', ...
    optional_extension,'.mat');
save(Matrices_and_SS_file,'inv_reduced_mass_matrix', ...
                          'reduced_damping_matrix', ...
                          'reduced_stiffness_matrix', ...
                          'sys', ...
                          'cond_num', ...
                          'size_Arad', ...
                          'max_real_e_value', ...
                          'max_omega', ...
                          'fit_errormax', ...
                          'datascaler_in');

%% REFERENCES
% [1] Forehand, D.I.M., Kiprakis, A.E., Nambiar, A.J. and Wallace, A.R.
%     (2016) “A fully coupled wave-to-wire model of an array of wave energy
%     converters”. IEEE Transactions on Sustainable Energy. Vol. 7, No. 1,
%     pp. 118-128.
% [2] The WAMIT User Manual (version 7.0)
