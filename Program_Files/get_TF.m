function [TFnum,TFden,fitdata,no_TFs] = get_TF(modeflag,fitdatalength,h,omegas,b_order_in,a_order_in,iter,tol)
% function [TFnum,TFden,fitdata,no_TFs] = get_TF(modeflag,fitdatalength,h,omegas,b_order_in,a_order_in,iter,tol)
% OR -----
% function [TFnum,TFden,fitdata,no_TFs] = get_TF(modeflag,fitdatalength,h,omegas,b_order_in,a_order_in)
%
% Returns numerator and denominator coefficients of transfer function 
% with frequency response closest to that input to the function, according
% to the selection criterion (at the moment least relative 
% root-mean-square error between both real and imaginary parts).
%
%**************************************************************************
% INPUT ARGUMENTS:
%
% MANDATORY:
%
%         modeflag ........ row vector of indices of transfer function mode.
%         fitdatalength ... number of TFs to be included in 'fitdata' output.
%         h ............... complex frequency response vector.
%         omegas .......... vector of angular frequencies corresponding to 'h'.
%         b_order_in ...... >0: order of TF numerator to be derived.
%                    ...... <0: negative of lowest TF order to be derived.   
%         a_order_in ...... >0: order of TF denominator to be derived.
%                    ...... <-2: negative of highest TF order to be derived. 
%
% OPTIONAL (inclusion of these arguments forces the derivation of a stable
%           transfer function):
%
%         iter ............ maximum number of iterations in numerical search.
%         tol ............. norm of the gradient for numerical search stop.
%
%**************************************************************************
% OUTPUT ARGUMENTS:
%
%         TFnum ..... vector of coefficients of selected transfer function numerator.
%         TFden ..... vector of coefficients of selected transfer function denominator.
%         fitdata ... matrix of selection criteria values with associated
%                     denominator and numerator orders. 
%         no_TFs .... number of transfer functions derived and assessed.
%
%**************************************************************************

%**************************************************************************
% Andy McCabe
% Lancaster University
%
% Original form - 04-05-2006
% Revised form  - 06-11-2007
% 
%**************************************************************************

%**************************************************************************
% Set minimum number of TFs to assess
no_TFs = 1;

% Mean of frequency response magnitude
meanmag = mean(abs(h));

% Test for inclusion of 'iter' and 'tol'
if nargin > 7
    %**************************************************************************
    % STABLE (RADIATION) TRANSFER FUNCTION(S)            
    % Test for multiple TF assessment
    if a_order_in <= -2
        maxorder = -a_order_in;     % set highest denominator order
        minorder = -b_order_in;     % set lowest denominator order
        no_TFs   = sum([(minorder-1):(maxorder-1)]); % total number of TFs

        % Initialize r.m.s. error matrix
        errs = [0 0 0];

        if ~modeflag(3)
            clc
            fprintf('\n');
            fprintf('RADIATION TRANSFER FUNCTION \n\n');
            fprintf('mode   = %2d  %2d \n\n',modeflag(1),modeflag(2));
        end
        
        % Loop for multiple TF assessment
        for a_order = minorder:maxorder
            for b_order = 1:a_order-1

                % Something to watch while it loops
                if modeflag(3)
                    clc
                    fprintf('\n');
                    fprintf('RADIATION TRANSFER FUNCTION \n\n');
                    fprintf('mode   = %2d  %2d \n\n',modeflag(1),modeflag(2));
                    fprintf('orders = %2d  %2d \n\n',a_order,b_order);
                end

                % Stop a lot of warning messages flashing up for no useful purpose
                warning off all

                % Inverse Frequency-Response Function
                [num,den] = invfreqs(h,omegas,b_order,a_order,[],iter,tol);

                % Turn the warnings back on again for the rest of the program
                warning on all

                % Test for stability (belt-and-braces approach) and
                % calculate frequency response error
                if max(real(roots(den)))<=-1e-3
                    sys = tf(num,den);
                    f1 = freqresp(sys,omegas);
                    fr = squeeze(f1);

                    errors  = [[real(h) - real(fr)];[imag(h) - imag(fr)]];
                    sqerrs  = errors.^2;
                    rmserrs = sqrt(mean(sqerrs))/meanmag;

                    % store result
                    errs = [errs; a_order, b_order, rmserrs];
                end
            end
        end
        % Remove initial zero row
        errs1 = errs(2:size(errs,1),:);

        % Sort for best fit
        [errs2, indr] = sort(errs1(:,3)); 

        % List results to required number in order of best fit
        if size(errs1,1) < fitdatalength
            fitdatalength = length(errs2);
        end
        fitdata = [errs1(indr(1:fitdatalength),1),errs1(indr(1:fitdatalength),2),errs2(1:fitdatalength)];

        % Set numerator and denominator orders for best-fit TF
        a_order_in = errs1(indr(1),1);
        b_order_in = errs1(indr(1),2);

        % Stop a lot of warning messages flashing up for no useful purpose
        warning off all

        % TF numerator and denominator coefficient vectors
        [TFnum,TFden] = invfreqs(h,omegas,b_order_in,a_order_in,[],iter,tol);

        % Turn the warnings back on again for the rest of the program
        warning on all

    % Otherwise single TF assessment
    else
        % Stop a lot of warning messages flashing up for no useful purpose
        warning off all

        % TF numerator and denominator coefficient vectors
        [TFnum,TFden] = invfreqs(h,omegas,b_order_in,a_order_in,[],iter,tol);

        % Turn the warnings back on again for the rest of the program
        warning on all

        % calculate frequency response error
        sys = tf(TFnum,TFden);
        f1 = freqresp(sys,omegas);
        fr = squeeze(f1);

        errors  = [[real(h) - real(fr)];[imag(h) - imag(fr)]];
        sqerrs  = errors.^2;
        rmserrs = sqrt(mean(sqerrs))/meanmag;

        % store result
        fitdata = [a_order_in, b_order_in, rmserrs];
    end

    %************************************************************************** 

% exclusion of 'iter' and 'tol'
else
    %**************************************************************************
    % POSSIBLY UNSTABLE (WAVE EXCITATION) TRANSFER FUNCTION(S)            
    % Test for multiple TF assessment
    if a_order_in <= -2
        maxorder = -a_order_in;     % set highest denominator order
        minorder = -b_order_in;     % set lowest denominator order
        no_TFs   = sum([(minorder-1):(maxorder-1)]); % total number of TFs

        % Initialize r.m.s. error matrix
        errs = [0 0 0];

        if ~modeflag(3)
            clc
            fprintf('\n');
            fprintf('WAVE EXCITATION TRANSFER FUNCTION \n\n');
            fprintf('mode   = %2d  \n\n',modeflag(1));
        end

        % Loop for multiple TF assessment
        for a_order = minorder:maxorder
            for b_order = 1:a_order-1
                
                % Something to watch while it loops
                if modeflag(3)
                    clc
                    fprintf('\n');
                    fprintf('WAVE EXCITATION TRANSFER FUNCTION \n\n');
                    fprintf('mode   = %2d  \n\n',modeflag(1));
                    fprintf('orders = %d  %d \n\n',a_order,b_order);
                end

                % Stop a lot of warning messages flashing up for no useful purpose
                warning off all

                % Inverse Frequency-Response Function
                [num,den] = invfreqs(h,omegas,b_order,a_order);

                % Turn the warnings back on again for the rest of the program
                warning on all

                % calculate frequency response error
                sys = tf(num,den);
                f1 = freqresp(sys,omegas);
                fr = squeeze(f1);

                errors  = [[real(h) - real(fr)];[imag(h) - imag(fr)]];
                sqerrs  = errors.^2;
                rmserrs = sqrt(mean(sqerrs))/meanmag;

                % store result
                errs = [errs; a_order, b_order, rmserrs];
            end
        end
        % Remove initial zero row
        errs1 = errs(2:size(errs,1),:);

        % Sort for best fit
        [errs2, indr] = sort(errs1(:,3)); 

        % List results to required number in order of best fit
        if size(errs1,1) < fitdatalength
            fitdatalength = length(errs2);
        end
        fitdata = [errs1(indr(1:fitdatalength),1),errs1(indr(1:fitdatalength),2),errs2(1:fitdatalength)];

        % Set numerator and denominator orders for best-fit TF
        a_order_in = errs1(indr(1),1);
        b_order_in = errs1(indr(1),2);

        % Stop a lot of warning messages flashing up for no useful purpose
        warning off all
        
        % TF numerator and denominator coefficient vectors
        [TFnum,TFden] = invfreqs(h,omegas,b_order_in,a_order_in);
        
        % Turn the warnings back on again for the rest of the program
        warning on all

    % Otherwise single TF assessment
    else
        % Stop a lot of warning messages flashing up for no useful purpose
        warning off all
        
        % TF numerator and denominator coefficient vectors
        [TFnum,TFden] = invfreqs(h,omegas,b_order_in,a_order_in);

        % Turn the warnings back on again for the rest of the program
        warning on all

        % calculate frequency response error
        sys = tf(TFnum,TFden);
        f1 = freqresp(sys,omegas);
        fr = squeeze(f1);

        errors  = [[real(h) - real(fr)];[imag(h) - imag(fr)]];
        sqerrs  = errors.^2;
        rmserrs = sqrt(mean(sqerrs))/meanmag;

        % store result
        fitdata = [a_order_in, b_order_in, rmserrs];
    end
end

%************************************************************************** 
    
%************************************************************************** 
% END   ******************************************************************* 
%************************************************************************** 

