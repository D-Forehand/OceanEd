%% CODE DESCRIPTION AND AUTHORS
%
% *************************** Code Description ****************************
% 
% This code checks that the hydrodynamic time-domain model is working
% correctly (by comparing against the frequency-domain results) for all the
% incident wave frequencies (omegas) and all the incident wave directions
% (betas). For an incident sinusoidal wave of a particular frequency
% (omega) and wave direction (beta), the time-domain body responses should
% match the frequency-domain body responses (in amplitude and phase) once
% the initial transients in the time-domain motions have died away. Thus
% the time-domain simulations will need to be run for long enough for the
% initial transients to die away and for steady-state oscillatory behaviour
% to be reached. The main part of this code runs the time-domain model for
% all the omegas and all the betas.  This means that the time-domain model
% will be run for (num_omega)x(num_beta) simulations, were num_omega is the
% number of omegas and num_beta is the number of betas. The last section of
% this code plots the Response Amplitude Operator (RAO) and error surfaces.
% For a particular mode of motion, the top plot is the magnitude of the RAO
% of that mode plotted against all omegas and all betas. The middle plot is
% the relative error in amplitude for that mode between the frequency and
% time-domain results, plotted against all omegas and all betas. The bottom
% plot is the error in phase (in degrees) for that mode between the
% frequency and time-domain results, plotted against all omegas and all
% betas.
%
% This code uses the Simulink model
% “time_domain_array_sinusoidal_input.slx”, which is used to run the
% time-domain model for an incident sinusoidal wave. The code also uses the
% function “least_squares_sine_fit.m”, which fits a sine wave to the
% time-domain responses (which are assumed to have reached steady-state
% oscillatory conditions).
%
% ******************************** Author *********************************
%
% Dr David Forehand
% The University of Edinburgh
%
% This version - 01/12/2018
%
% *************************************************************************


%% CLEAR THE COMMAND WINDOW, THE WORKSPACE, AND CLOSE ANY OPEN FIGURES

clc % Clear command window
clearvars % Clear variables from memory
clear global % Clear all global variables
close all % Closes any currently open figures

%% TURN ALL WARNINGS OFF

warning off all

%% INITIALISING CONSTANTS

two_pi=2*pi;
pi_over_180=pi/180;

%% INPUT THE NAME OF THE WAMIT FRC FILE

string{01}='Give  the  filename  of  the  WAMIT FRC file,';
string{02}='without the ".frc" extension: ';
string{03}=sprintf('%s\n%s\n\n',string{01},string{02});
string{04}='WAMIT  filenames   without   their  extension';
string{05}='should  be  less  than  16  characters  long.';
string{06}='Please re-enter the filename: ';
string{07}=sprintf('%s\n%s\n%s\n\n',string{04},string{05},string{06});
frc_filename=input(string{03}, 's');
% The  's'  here  returns  the  entered  string  as a text
% variable rather  than  as a variable name  or  numerical
% value.
while length(frc_filename)>16
    disp(' ');
    frc_filename=input(string{07}, 's');
end
disp(' ');

clear('string')

%% INPUT THE FILE EXTENSION OF THE STATE-SPACE MATRICES INPUT FILE

string{01}='Give the optional extension that was added to';
string{02}='the input data (.mat) file which contains the';
string{03}='inverse  of  the  reduced  mass  matrix,  the';
string{04}='reduced damping  and stiffness matrices,  and';
string{05}='the  state-space  system  for  the  radiation';
string{06}='terms.  The  filename  of  this  file  was  a';
string{07}='concatenation     of     the    frc_filename,';
string{08}='"_MassDampStiffMatsRadSS", this optional file';
string{09}='extension and ".mat",  where the frc_filename';
string{10}='is  the filename of the WAMIT FRC file  input';
string{11}='above.';
string{12}='(This optional extension  should be less than';
string{13}='16 characters long):';
string{14}=sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n', ...
    string{01},string{02},string{03},string{04}, ...
    string{05},string{06},string{07},string{08}, ...
    string{09},string{10},string{11},string{12}, ...
    string{13});
string{15}='This optional extension  should be  less than';
string{16}='16 characters long.  Please re-enter it:';
string{17}=sprintf('%s\n%s\n\n',string{15},string{16});
optional_extension=input(string{14}, 's');
while length(optional_extension)>16
    disp(' ');
    optional_extension=input(string{17}, 's');
end
disp(' ');

clear('string')

%% INPUT WHETHER HASKIND RELATIONS OR DIFFRACTION POTENTIAL IS USED
% Note  if NBETAH is set to a value greater than zero  in
% the main Force Control File, then the excitation forces
% and body motions would come from the Haskind relations.
% That is, in particular, IOPTN(3) would be automatically
% set to 0  and  I think you would not be allowed  to set
% IOPTN(4)  to  +2  or  -2  as  this  would  generate  an
% inconsistency error.  See pages 4-12 to 4-17 of [1] for
% more details.

string{01}='Are  the excitation forces  and  body motions';
string{02}='calculated using the Haskind relations or the';
string{03}='diffraction potential?';
string{04}='Type 1 for the Haskind relations or 2 for the';
string{05}='diffraction potential:';
string{06}=sprintf('%s\n%s\n%s\n%s\n%s\n\n', ...
    string{01},string{02},string{03},string{04},string{05});
string{07}='Please re-enter your selection  - either type';
string{08}='1 for  the Haskind relations  or  2  for  the';
string{09}='diffraction potential:';
string{10}=sprintf('%s\n%s\n%s\n\n', ...
    string{07},string{08},string{09});
excitation_force_string=input(string{06}, 's');
while (~strcmp(excitation_force_string,'1')) && (~strcmp(excitation_force_string,'2'))
    % The function strcmp compares strings.
    disp(' ');
    excitation_force_string=input(string{10}, 's');
end
excitation_force_index=str2num(excitation_force_string);
disp(' ');

clear('string','excitation_force_string')

%% LOAD THE DATA FOR THE TIME-DOMAIN HYDRODYNAMIC ARRAY MODEL

load(strcat(frc_filename,'_MassDampStiffMatsRadSS', ...
    optional_extension,'.mat'));

clear('optional_extension')

%% LOAD THE EXCITATION FORCES

string{01}='The WAMIT file ';
string{02}='which contains the excitation forces from the';
string{03}='Haskind relations does not exist.';
string{04}=sprintf('%s%s\n%s\n%s',string{01}, ...
    strcat(frc_filename,'.2'),string{02},string{03});
string{05}='diffraction potential does not exist.';
string{06}=sprintf('%s%s\n%s\n%s',string{01}, ...
    strcat(frc_filename,'.3'),string{02},string{05});
switch excitation_force_index
    case{1}
        if exist(strcat(frc_filename,'.2')) ~= 2
            error(string{04});
        else
            excite_forces=load(strcat(frc_filename,'.2'));
        end
    case{2}
        if exist(strcat(frc_filename,'.3')) ~= 2
            error(string{06});
        else
            excite_forces=load(strcat(frc_filename,'.3'));
        end
end

clear('string','excitation_force_index')

%% LOAD THE BODY MOTIONS

string{01}='The WAMIT file ';
string{02}='which  contains  the  body motions  does  not';
string{03}='exist.';
string{04}=sprintf('%s%s\n%s\n%s',string{01}, ...
    strcat(frc_filename,'.4'),string{02},string{03});
if exist(strcat(frc_filename,'.4')) ~= 2
    error(string{04});
else
    body_motions=load(strcat(frc_filename,'.4'));
end
   
clear('string')

%% INPUT THE WAVE AMPLITUDE 
% This  is  not really necessary  for this code  but it is
% added here for completeness.

string{01}='Give the wave amplitude in metres:';
string{02}=sprintf('%s\n\n',string{01});
amplitude=input(string{02});

clear('string')

%% CALCULATING THE NUMBER OF DIFFERENT MODES

num_modes=1;
first_mode=excite_forces(1,3);
next_mode=excite_forces(2,3);
while next_mode ~= first_mode
    num_modes=num_modes+1;
    next_mode=excite_forces(num_modes+1,3);
end
mode=excite_forces(1:num_modes,3);

clear('first_mode','next_mode')

%% CALCULATING THE NUMBER OF BETAS (WAVE HEADING ANGLES)

num_betas=1;
first_beta=excite_forces(1,2);
next_beta=excite_forces(num_modes+1,2);
while next_beta ~= first_beta
    num_betas=num_betas+1;
    next_beta=excite_forces(num_betas*num_modes+1,2);
end
beta=excite_forces(1:num_modes:num_betas*num_modes,2);

clear('first_beta','next_beta')

%% EXTRACTING THE FREQUENCY DATA

omega=two_pi./excite_forces(1:num_betas*num_modes:end,1);
num_omegas=length(omega);

%% EXTRACT g, L AND rho FROM THE WAMIT ".OUT" FILE

fid=fopen(strcat(frc_filename,'.out'),'rt');%open the ".out" file
% On a Windows system, the 'rt' opens the file for reading in text mode
% instead of binary mode (the default).
if fid < 0
   error(['Could not open ',frc_filename,'.out']);
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
rho=cell2mat(textscan(s,'%*s %*s %*f %*s %*s %f'));

% Calculating     the    multipliers     necessary     for
% dimensionalising  the excitation forces   (see  the last
% equation in Section 4.3 on page 4-3 of [1]):
rho_g_A_L_squared=rho*g*amplitude*(L^2);
rho_g_A_L_cubed=rho*g*amplitude*(L^3);

fclose(fid);%close the ".out" file

clear('ans','g_and_L_data','fid','k','s')

%% EXAMINATION OF BODY MOTION TIME-SERIES

string{01}='EXAMINATION OF BODY MOTIONS';
string{02}='Do you want to examine  the body motion time-';
string{03}='series  for  one  frequency  (omega)  of  the';
string{04}='incident wave and one wave heading (beta)?';
string{05}='Type "y" for yes or "n" for no:';
string{06}=sprintf('\n%s\n%s\n%s\n%s\n%s\n\n',string{01} ...
    ,string{02},string{03},string{04},string{05});
string{07}='You can choose  an omega  from  the following';
string{08}='list:';
string{09}=sprintf('\n%s\n%s',string{07},string{08});
string{10}='Enter the omega you want from the above list:';
string{11}=sprintf('\n%s\n\n',string{10});
string{12}='That omega is not in the list above,   please';
string{13}='choose again:';
string{14}=sprintf('\n%s\n%s\n\n',string{12},string{13});
string{15}='You  can  choose  a beta  from  the following';
string{16}='list:';
string{17}=sprintf('\n%s\n%s',string{15},string{16});
string{18}='Enter the beta you want  from the above list:';
string{19}=sprintf('\n%s\n\n',string{18});
string{20}='That beta is not in the list above,    please';
string{21}='choose again:';
string{22}=sprintf('\n%s\n%s\n\n',string{20},string{21});
string{23}='Enter the finishing time for the simulation,';
string{24}='this should be a real number greater than 0:';
string{25}=sprintf('\n%s\n%s\n\n',string{23},string{24});

program_test=input(string{06}, 's');
while (strcmp(program_test,'y')) || (strcmp(program_test,'"y"') || ...
        strcmp(program_test,'Y')) || (strcmp(program_test,'"Y"'))
    
    disp(string{09});
    disp(omega');
    if num_omegas>1
        picked_omega=input(string{11});
        while ~isequal(size(picked_omega),[1 1]) || ...    %not a scalar
                ~isreal(picked_omega) || ...               %not real
                min(abs(omega-picked_omega)) > 1.0e-4      %not in the list
            picked_omega=input(string{14});
        end
        [temp_omega_error,omega_num]=min(abs(omega-picked_omega));
    else
        omega_num=1;        
    end

    disp(string{17});
    disp(beta');
    if num_betas>1
        picked_beta=input(string{19});
        while ~isequal(size(picked_beta),[1 1]) || ...     %not a scalar
                ~isreal(picked_beta) || ...                %not real
                min(abs(beta-picked_beta)) > 1.0e-4        %not in the list
            picked_beta=input(string{22});
        end
        [temp_beta_error,beta_num]=min(abs(beta-picked_beta));
    else
        beta_num=1;
    end
       
    string{26}='The omega and beta you have chosen are:';
    string{27}=['omega = ' num2str(omega(omega_num)) ...
        ' and beta = ' num2str(beta(beta_num))];
    string{28}=sprintf('\n%s\n%s',string{26},string{27});
    disp(string{28});
    
    t_final=input(string{25});
    while ~isequal(size(t_final),[1 1]) || ...              %not a scalar
                ~isreal(t_final) || ...                     %not real
                t_final < 0
            t_final=input(string{25});
    end
    
    offset=num_betas*num_modes*(omega_num-1)+num_modes*(beta_num-1);

    % Calculating the vector of  the dimensional magnitdes
    % of the excitation forces for the selected omega  and
    % the selected beta:
    Mod_X_i=zeros(1,num_modes);%creating a row vector of zeros
    for mode_count=1:num_modes
        mode_between_1_and_6=mode(mode_count) ...
            -ceil(mode(mode_count))*6+6;
        if mode_between_1_and_6<3.5
            Mod_X_i(mode_count)= ...
                rho_g_A_L_squared*excite_forces(offset+mode_count,4);
        else
            Mod_X_i(mode_count)= ...
                rho_g_A_L_cubed*excite_forces(offset+mode_count,4);
        end
    end

    % Calculating the vector of phases (in radians) of the
    % excitation forces  for  the selected omega  and  the
    % selected beta:
    Pha_X_i=pi_over_180* ...
        excite_forces(offset+1:offset+num_modes,5);
    
    excitation_frequency=omega(omega_num);
    
    [test_t,test_x,test_y]= ...
    sim('time_domain_array_sinusoidal_input',t_final);

    % Calculating the vector of  the dimensional magnitdes
    % of  the  WAMIT  calculated  body  motions   for  the
    % selected omega and the selected beta:
    Mod_Xi_i=zeros(1,num_modes);%creating a row vector of zeros
    for mode_count=1:num_modes
        mode_between_1_and_6=mode(mode_count) ...
            -ceil(mode(mode_count))*6+6;
        if mode_between_1_and_6<3.5
            Mod_Xi_i(mode_count)= ...
                amplitude*body_motions(offset+mode_count,4);
        else
            Mod_Xi_i(mode_count)= ...
                (amplitude/L)*body_motions(offset+mode_count,4);
        end
    end
    
    % Calculating the vector of phases (in radians) of the
    % WAMIT calculated body motions for the selected omega
    % and the selected beta:
    Pha_Xi_i=pi_over_180* ...
        body_motions(offset+1:offset+num_modes,5);
    
    wamit_y=zeros(length(test_t),num_modes);
    for mode_count=1:num_modes
        wamit_y(:,mode_count)=Mod_Xi_i(mode_count)* ...
            sin(excitation_frequency*test_t+Pha_Xi_i(mode_count));
    end
    
    string{29}='Do you want to see a plot of the response  of';
    string{30}='one of the modes?';
    string{31}='Type "y" for yes or "n" for no:';
    string{32}=sprintf('\n%s\n%s\n%s\n\n',string{29} ...
        ,string{30},string{31});
    string{33}='Enter the mode number you want plotted:';
    string{34}=sprintf('%s\n\n',string{33});
    string{35}='Please re-enter a valid mode number:';
    string{36}=sprintf('\n%s\n\n',string{35});
    string{37}='Do you want to see the plot for another mode?';
    string{38}=sprintf('\n%s\n%s\n\n',string{37},string{31});
    plot_or_not=input(string{32}, 's');
    while (strcmp(plot_or_not,'y')) || (strcmp(plot_or_not,'"y"') || ...
            strcmp(plot_or_not,'Y')) || (strcmp(plot_or_not,'"Y"'))
        disp(' ');
        disp(['There are ' num2str(num_modes) ' modes in total.']);
        disp('These modes are:');
        disp(num2str(mode','%6.0d'));
        mode_num=input(string{34});
        while ~isequal(size(mode_num),[1 1]) || ...      %not a scalar
                ~ismember(mode_num,mode)                 %not a mode number
            mode_num=input(string{36});
        end
        [unused_var,mode_index]=ismember(mode_num,mode);
        plot(test_t,test_y(:,mode_index), ...
            test_t,wamit_y(:,mode_index));
        title(['Response for Mode ' num2str(mode_num) ...
            ',  Omega = ' num2str(excitation_frequency) ...
            ' rad/s  and  Beta = ' num2str(beta(beta_num)) ' degrees']);
        xlabel('Time (s)');
        ylabel('Displacement (m)');
        legend('Time-domain solution','Frequency-domain solution (WAMIT)');
        plot_or_not=input(string{38}, 's');
    end
    string{39}='Do you want to examine  the body motion time-';
    string{40}='series  for another frequency (omega)  of the';
    string{41}='incident  wave   and   another  wave  heading';
    string{42}='(beta)?';
    string{43}='Type "y" for yes or "n" for no:';
    string{44}=sprintf('\n%s\n%s\n%s\n%s\n%s\n\n',string{39} ...
        ,string{40},string{41},string{42},string{43});
    program_test=input(string{44}, 's');
end

clear('beta_num','excitation_frequency','Mod_X_i','Mod_Xi_i', ...
      'mode_between_1_and_6','mode_count','mode_index','mode_num', ...
      'offset','omega_num','Pha_X_i','Pha_Xi_i','picked_beta', ...
      'picked_omega','plot_or_not','program_test','string','t_final', ...
      'temp_beta_error','temp_omega_error','test_t','test_x','test_y', ...
      'unused_var','wamit_y')

%% INPUT OF TIME WHEN STEADY STATE IS REACHED FOR ALL OMEGA AND BETA

string{01}='Give the time at which you think steady state';
string{02}='will have definitely been reached for all the';
string{03}='wave frequencies (omega)   and  wave headings';
string{04}='(beta):';
string{05}=sprintf('\n%s\n%s\n%s\n%s\n\n',string{01},string{02}, ...
    string{03},string{04});
string{06}='This  should be  a number  greater than zero.';
string{07}='Please re-enter it:';
string{08}=sprintf('\n%s\n%s\n\n',string{06},string{07});

t_steady_state=input(string{05});
while ~isequal(size(t_steady_state),[1 1]) || ...             %not a scalar
        ~isreal(t_steady_state) || ...                            %not real
        t_steady_state < 0
    t_steady_state=input(string{08});
end

clear('string')

%% INPUT OF THE SIMULATION END TIME

string{01}='Give the simulation end time.  This should be';
string{02}='sufficiently  bigger   than   the  previously';
string{03}='inputted time  in order for the least squares';
string{04}='sine fit  to  work properly.    However,   it';
string{05}='shouldn''t be too big,  otherwise  the set of';
string{06}='simulations will take too long to run:';
string{07}=sprintf('\n%s\n%s\n%s\n%s\n%s\n%s\n\n',string{01}, ...
    string{02},string{03},string{04},string{05},string{06});
string{08}=['This should be a number greater than ' ...
    num2str(t_steady_state)];
string{09}='Please re-enter it:';
string{10}=sprintf('\n%s\n%s\n\n',string{08},string{09});

t_end=input(string{07});
while ~isequal(size(t_end),[1 1]) || ...                      %not a scalar
        ~isreal(t_end) || ...                                     %not real
        t_end < t_steady_state
    t_end=input(string{10});
end

clear('string')

%% MAIN LOOP TO RUN THE TIME-DOMAIN MODEL FOR ALL WAVE FREQUENCIES (OMEGAS)
% AND ALL WAVE HEADING ANGLES (BETAS) 

% c = cell(m, n) or c = cell([m, n]) creates an m-by-n cell array of empty
% matrices. Arguments m and n must be scalars:
relative_amp_error=cell(1,num_modes);
phase_error=cell(1,num_modes);
RAO=cell(1,num_modes);

for mode_count=1:num_modes
    relative_amp_error{mode_count}=zeros(num_betas,num_omegas);
    phase_error{mode_count}=zeros(num_betas,num_omegas);
    RAO{mode_count}=zeros(num_betas,num_omegas);
end
    
for beta_num=1:num_betas
    for omega_num=1:num_omegas
        
        clc
        fprintf('There are %4d betas  (wave headings) in total and \n', ...
                num_betas);
        fprintf('there are %4d omegas (wave frequencies) in total.\n\n', ...
                num_omegas);
        fprintf('Now performing time-domain simulations for all \n'); ...
        fprintf('betas and all omegas: \n\n');    
        fprintf('beta_num, omega_num = %4d %4d \n\n',beta_num,omega_num);
        
        offset=num_betas*num_modes*(omega_num-1)+num_modes*(beta_num-1);
        
        % Calculating the vector of  the dimensional magnitdes
        % of the excitation forces for the selected omega  and
        % the selected beta:
        Mod_X_i=zeros(1,num_modes);%creating a row vector of zeros
        for mode_count=1:num_modes
            mode_between_1_and_6=mode(mode_count) ...
                -ceil(mode(mode_count))*6+6;
            if mode_between_1_and_6<3.5
                Mod_X_i(mode_count)= ...
                    rho_g_A_L_squared*excite_forces(offset+mode_count,4);
            else
                Mod_X_i(mode_count)= ...
                    rho_g_A_L_cubed*excite_forces(offset+mode_count,4);
            end
        end
        
        % Calculating the vector of phases (in radians) of the
        % excitation forces  for  the selected omega  and  the
        % selected beta:
        Pha_X_i=pi_over_180* ...
            excite_forces(offset+1:offset+num_modes,5);
        
        % Calculating the vector of  the dimensional magnitdes
        % of  the  WAMIT  calculated  body  motions   for  the
        % selected omega and the selected beta:
        Mod_Xi_i=zeros(1,num_modes);%creating a row vector of zeros
        for mode_count=1:num_modes
            mode_between_1_and_6=mode(mode_count) ...
                -ceil(mode(mode_count))*6+6;
            if mode_between_1_and_6<3.5
                Mod_Xi_i(mode_count)= ...
                    amplitude*body_motions(offset+mode_count,4);
            else
                Mod_Xi_i(mode_count)= ...
                    (amplitude/L)*body_motions(offset+mode_count,4);
            end
        end
        
        for mode_count=1:num_modes
            RAO{mode_count}(beta_num,omega_num)= ...
                body_motions(offset+mode_count,4);
        end

        % Calculating the vector of phases (in radians) of the
        % WAMIT calculated body motions for the selected omega
        % and the selected beta:
        Pha_Xi_i=pi_over_180* ...
            body_motions(offset+1:offset+num_modes,5);
        
        excitation_frequency=omega(omega_num);
        
        [t_sim,temp_states,y_sim]= ...
            sim('time_domain_array_sinusoidal_input',t_end);
        
        % Finding the index  of the first element  of the time
        % column   vector  t_sim   which   is   greater   than
        % t_steady_state.
        start_index=find(t_sim>t_steady_state,1,'first');
        
        for mode_count=1:num_modes
            [td_amp,td_phase]= ...
                least_squares_sine_fit(excitation_frequency, ...
                t_sim(start_index:end)', ...
                y_sim(start_index:end,mode_count)');
            fd_amp=Mod_Xi_i(mode_count);
            fd_phase=Pha_Xi_i(mode_count);
            relative_amp_error{mode_count}(beta_num,omega_num)= ...
                (td_amp-fd_amp)/fd_amp;
            phase_error{mode_count}(beta_num,omega_num)= ...
                td_phase-fd_phase;
            if phase_error{mode_count}(beta_num,omega_num)>two_pi
                phase_error{mode_count}(beta_num,omega_num)= ...
                    phase_error{mode_count}(beta_num,omega_num)-two_pi;
            elseif phase_error{mode_count}(beta_num,omega_num)<-two_pi
                phase_error{mode_count}(beta_num,omega_num)= ...
                    phase_error{mode_count}(beta_num,omega_num)+two_pi;
            end
        end
        
    end
end

fprintf('Finished all the time-domain simulations. \n');

clear('beta_num','excitation_frequency','fd_amp','fd_phase','Mod_X_i', ...
      'Mod_Xi_i','mode_between_1_and_6','mode_count','offset', ...
      'omega_num','Pha_X_i','Pha_Xi_i','start_index','t_sim','td_amp', ...
      'td_phase','temp_states','y_sim')

% Save all the workspace variables to a file in case we don't want to rerun
% all the above code (which can take a long time) and we just want to
% experiment with the code below: 
save('error_surf_data')

%% SECTION TO PLOT THE RAO AND ERROR SURFACES

% Load all the workspace variables from a file when we don't want to rerun
% all the above code (which can take a long time) and we just want to
% experiment with the code below: 
clearvars
load('error_surf_data')

for mode_count=1:num_modes
    for beta_num=1:num_betas
        for omega_num=1:num_omegas
            if abs(abs(phase_error{mode_count}(beta_num,omega_num))-two_pi) < 0.3
                phase_error{mode_count}(beta_num,omega_num) = ...
                    phase_error{mode_count}(beta_num,omega_num) ...
                    - sign(phase_error{mode_count}(beta_num,omega_num))*two_pi;
            end
        end
    end 
end

% The  following  code  calculates  the  maximum absolute
% relative error in amplitude,  the rms relative error in
% amplitude,  the maximum absolute error in phase and the
% rms error in phase  for each  of  the modes.    For the
% phase errors,  it produces these  in both  radians  and
% degrees. 
rms_relative_amp_error=zeros(1,num_modes);
rms_phase_error=zeros(1,num_modes);
max_abs_relative_amp_error=zeros(1,num_modes);
max_abs_phase_error=zeros(1,num_modes);
num_points=num_betas*num_omegas;
for mode_count=1:num_modes
    % Note sum(sum(A)) sums all the elements of a 2D array (i.e. matrix) A.
    % See the help on "sum".
    rms_relative_amp_error(mode_count)= ...
        sqrt(sum(sum(relative_amp_error{mode_count}.^2))/num_points);
    rms_phase_error(mode_count)= ...
        sqrt(sum(sum(phase_error{mode_count}.^2))/num_points);
    % Note max(max(A)) returns the largest element of a 2D array (i.e.
    % matrix) A.  See the help on "max".
    max_abs_relative_amp_error(mode_count)= ...
        max(max(abs(relative_amp_error{mode_count})));
    max_abs_phase_error(mode_count)= ...
        max(max(abs(phase_error{mode_count})));
end
rms_phase_error_degrees=rms_phase_error/pi_over_180;
max_abs_phase_error_degrees=max_abs_phase_error/pi_over_180;

% The following code  displays all the errors  calculated
% above  in the command window.    It does this  for each
% mode.
string{01}='The maximum absolute relative error in amplitude is:';
string{02}='The rms relative error in amplitude is:';
string{03}='The maximum absolute error in phase is:';
string{04}='The rms error in phase is:';
string{05}='--------------------------';
for mode_count=1:num_modes
    string{06}=['For MODE ' num2str(mode(mode_count))];
    string{07}=num2str(max_abs_relative_amp_error(mode_count),'%9.4e');
    string{08}=num2str(rms_relative_amp_error(mode_count),'%9.4e');
    string{09}=[num2str(max_abs_phase_error(mode_count),'%9.4e') ...
        ' radians = ' ...
        num2str(max_abs_phase_error_degrees(mode_count),'%9.4e') ...
        ' degrees'];
    string{10}=[num2str(rms_phase_error(mode_count),'%9.4e') ...
        ' radians = ' ...
        num2str(rms_phase_error_degrees(mode_count),'%9.4e') ...
        ' degrees'];
    string{11}=sprintf('\n%s\n\n%s\n%s\n\n%s\n%s\n\n%s\n%s\n\n%s\n%s\n\n%s', ...
        string{06},string{01},string{07},string{02},string{08}, ...
        string{03},string{09},string{04},string{10},string{05});
    disp(string{11});
end

% NEW FOR LOOP TO CONVERT THE PHASE ERROR TO DEGREES
for mode_count=1:num_modes
    phase_error{mode_count} = (180.0/pi)*phase_error{mode_count};
end

string{22}='Do you want to see the error plots for one of';
string{23}='the modes?';
string{24}='Type "y" for yes or "n" for no:';
string{25}=sprintf('\n%s\n%s\n%s\n\n',string{22} ...
    ,string{23},string{24});
string{26}='Enter the mode number you want plotted:';
string{27}=sprintf('%s\n\n',string{26});
string{28}='Please re-enter a valid mode number:';
string{29}=sprintf('\n%s\n\n',string{28});
string{30}='Do you want to see another mode?';
string{31}=sprintf('\n%s\n%s\n\n',string{30},string{24});
plot_or_not=input(string{25}, 's');
while (strcmp(plot_or_not,'y')) || (strcmp(plot_or_not,'"y"') || ...
        strcmp(plot_or_not,'Y')) || (strcmp(plot_or_not,'"Y"'))
    disp(' ');
    disp(['There are ' num2str(num_modes) ' modes in total.']);
    disp('These modes are:');
    disp(num2str(mode','%6.0d'));
    mode_num=input(string{27});
    while ~isequal(size(mode_num),[1 1]) || ...          %not a scalar
            ~ismember(mode_num,mode)                     %not a mode number
        mode_num=input(string{29});
    end
    [unused_var,mode_index]=ismember(mode_num,mode);
    
    if num_betas==1
        
        subplot(3,1,1);
        plot(omega,RAO{mode_index});
        title(['Modulus (Magnitude) of RAO for Mode ' num2str(mode_num)]);
        xlabel('Omega (radians/s)');
        ylabel('Modulus of RAO');
        
        subplot(3,1,2);
        plot(omega,relative_amp_error{mode_index});
        title({['Relative Error in Amplitude for Mode ' ...
            num2str(mode_num)]; ...
            ['(max abs rel error in amp = ' ...
            num2str(max_abs_relative_amp_error(mode_index),'%9.4e') ...
            ', rms rel error in amp = ' ...
            num2str(rms_relative_amp_error(mode_index),'%9.4e') ...
            ')']});
        xlabel('Omega (radians/s)');
        ylabel('Relative Error in Amplitude');
        
        subplot(3,1,3);
        plot(omega,phase_error{mode_index});
        title({['Error in Phase for Mode ' num2str(mode_num)]; ...
            ['(max abs error in phase = ' ...
            num2str(max_abs_phase_error_degrees(mode_index),'%9.4e') ...
            ' degrees, rms error in phase = ' ...
            num2str(rms_phase_error_degrees(mode_index),'%9.4e') ...
            ' degrees)']});
        xlabel('Omega (radians/s)');
        ylabel('Error in Phase (degrees)');
        
    elseif num_omegas~=1
        [omega_mesh,beta_mesh]=meshgrid(omega,beta);
        
        subplot(3,1,1);
        surfl(omega_mesh,beta_mesh,RAO{mode_index});
        shading interp
        colormap copper
        title(['Modulus (Magnitude) of RAO for Mode ' num2str(mode_num)]);
        xlabel('Omega (radians/s)');
        xlim([0 max(omega)]);
        ylabel('Beta (degrees)');
        ylim([min(beta) max(beta)]);
        set(gca,'YTick',[min(beta) (3*min(beta)+max(beta))/4 (min(beta)+max(beta))/2 (min(beta)+3*max(beta))/4 max(beta)]);
        zlabel('Modulus of RAO');
        
        subplot(3,1,2);
        surfl(omega_mesh,beta_mesh,relative_amp_error{mode_index});
        shading interp
        colormap copper
        title({['Relative Error in Amplitude for Mode ' ...
            num2str(mode_num)]; ...
            ['(max abs rel error in amp = ' ...
            num2str(max_abs_relative_amp_error(mode_index),'%9.4e') ...
            ', rms rel error in amp = ' ...
            num2str(rms_relative_amp_error(mode_index),'%9.4e') ...
            ')']});
        xlabel('Omega (radians/s)');
        xlim([0 max(omega)]);
        ylabel('Beta (degrees)');
        ylim([min(beta) max(beta)]);
        set(gca,'YTick',[min(beta) (3*min(beta)+max(beta))/4 (min(beta)+max(beta))/2 (min(beta)+3*max(beta))/4 max(beta)]);
        zlabel('Relative Error in Amplitude');
        
        subplot(3,1,3);
        surfl(omega_mesh,beta_mesh,phase_error{mode_index});
        shading interp
        colormap copper
        title({['Error in Phase for Mode ' num2str(mode_num)]; ...
            ['(max abs error in phase = ' ...
            num2str(max_abs_phase_error_degrees(mode_index),'%9.4e') ...
            ' degrees, rms error in phase = ' ...
            num2str(rms_phase_error_degrees(mode_index),'%9.4e') ...
            ' degrees)']});
        xlabel('Omega (radians/s)');
        xlim([0 max(omega)]);
        ylabel('Beta (degrees)');
        ylim([min(beta) max(beta)]);
        set(gca,'YTick',[min(beta) (3*min(beta)+max(beta))/4 (min(beta)+max(beta))/2 (min(beta)+3*max(beta))/4 max(beta)]);
        zlabel('Error in Phase (degrees)');
        
    end
    
    plot_or_not=input(string{31}, 's');
    
end

clear('string')

%% REFERENCES
% [1] The WAMIT User Manual (version 7.0)