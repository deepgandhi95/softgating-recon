function[f] = run_soft_recon(Patient,PfileName,functype,MAGorPHS,const_width,sigma,cal)

    if nargin == 6
        cal = false;
    end

    p = read_gehdr14_NH(PfileName);
    gate_file = which_gating(p.rdb.user4);

    % Read in the gating file data
    gating  = PCVIPRgatingExport(gate_file,p.rdb.da_xres,0);

    % Get the temporal ordering
    [~,si] = sort(gating.time,'ascend');

    % Make a copy of the gating data with the respiratory waveform set to a
    % linear ramp
    outputfinal = zeros(numel(gating.ecg),4);
    outputfinal(:,1) = gating.ecg;
    outputfinal(:,2) = 1:size(outputfinal,1);
    outputfinal(:,3) = gating.time;
    outputfinal(:,4) = gating.prep;

    % Write out the new gating file
    fid = fopen('external_gating_file','w');
    fwrite(fid,outputfinal,'int32','b');
    fclose(fid);

    % Load in the analog binning.
    load(MAGorPHS); %#ok<*LOAD>
    if exist('SoftGating_Phs','var')
        %myPhs = SoftGating_Phs;
        SoftGating_Phs(SoftGating_Phs == 8) = Inf; %#ok<*NODEF>
        %SoftGating_Phs = SoftGating_Phs;
    elseif exist('SoftGating_Mag','var')
        %myPhs = SoftGating_Mag;
        SoftGating_Mag(SoftGating_Mag == 8) = Inf;
        SoftGating_Phs = SoftGating_Mag;
    end
    % Loop through all the time frames and reconstruct
    nframes = 8;
    
    mkdir recon         % create directory for recon results

    if cal==false
        recon_cmd = sprintf('/home/gan5eh/local/bin/pcvipr_recon_binary -f %s -pils -dat_plus_dicom -export_cal -external_gating_weights_name RespWeight_soft.dat -external_gating_weights -resp_gate_signal bellows -resp_gate thresh -resp_gate_efficiency .99999 -external_timestamps -external_timestamps_name external_gating_file',PfileName);
    else
        recon_cmd = sprintf('/home/gan5eh/local/bin/pcvipr_recon_binary -f %s -pils -dat_plus_dicom -fcal %s -export_cal -external_gating_weights_name RespWeight_soft.dat -external_gating_weights -resp_gate_signal bellows -resp_gate thresh -resp_gate_efficiency .99999 -external_timestamps -external_timestamps_name external_gating_file',PfileName,cal);
    end

    for fr = 1:nframes 
        switch functype
            case 'Inverse'
                myPhs = mod(SoftGating_Phs + 4 + fr - 1, 8);  % Shift the binning to the frame we want to reconstruct is between 0-1
                f = ones(size(myPhs,1),1);
                f( myPhs >= (4-const_width/2) & myPhs <= (4+const_width/2)) = 1; %Replace 0.5 by 4  % set the middle of the bin to a hard 0.5
                f( myPhs < (4-const_width/2)) = (1./(abs(myPhs( myPhs < (4-const_width/2)) - (4-const_width/2)) + sigma)/(1/sigma));
                f( myPhs > (4+const_width/2)) = (1./(abs(myPhs( myPhs > (4+const_width/2)) - (4+const_width/2)) + sigma)/(1/sigma)); 
                f = f-min(f);
                f = f.*(1/max(f));   
                f(isnan(f)) = 0;                            % Remove any infinities
                w = f;
                w(si) = f;                                  % Reorder
                fid = fopen('RespWeight_soft.dat','w');     % Write to file
                fwrite(fid,w,'single');
                fclose(fid);
                
                
            case 'Exponential'
                myPhs = mod(SoftGating_Phs + 0.5 - (fr - 1), 8);
                f = zeros(size(myPhs));
                dn = abs(myPhs-0.5);
                dn(dn>4) = 8 - dn(dn>4);
                for i = 1:numel(myPhs)
                    if dn(i) > const_width 
                        f(i) = exp(-1.*sigma.*(dn(i) - const_width));
                    else
                        f(i) = 1;
                    end
                end
                f = f-min(f);
                f = f.*(1/max(f));
                w = f;
                w(si) = f;                                  % Reorder
                fid = fopen('RespWeight_soft.dat','w');     % Write to file
                fwrite(fid,w,'single');
                fclose(fid);
        
            case 'Linear'
                myPhs = mod(SoftGating_Phs + 0.5 - (fr - 1), 8);
                f = zeros(size(myPhs));
                dn = abs(myPhs-0.5);
                dn(dn>4) = 8 - dn(dn>4);
                x = dn;
               
                A = [1 const_width; 1 sigma];
                B = [1; 0];
                X = A\B;

                c = X(1);
                m = X(2);
                y = zeros(size(x));
                for i = 1:numel(x)
                    if x(i) <= const_width
                        y(i) = 1;
                    elseif x(i) > const_width
                        y(i) = m.*(x(i))+c;
                        if y(i)<0
                            y(i) = 0;
                        end
                    end
                end
                
                w = y;
                w(si) = y;                                  % Reorder
                fid = fopen('RespWeight_soft.dat','w');     % Write to file
                fwrite(fid,w,'single');
                fclose(fid);
                
            case 'Hardgating'
                myPhs = mod(SoftGating_Phs + 0.5 - (fr - 1), 8);
                f = zeros(size(myPhs));
                dn = abs(myPhs-0.5);
                dn(dn>4) = 8 - dn(dn>4);
                x = dn;
                
                A = [1 const_width; 1 sigma];
                B = [1; 0];
                X = A\B;

                c = X(1);
                m = X(2);
                y = zeros(size(x));
                for i = 1:numel(x)
                    if x(i) <= const_width
                        y(i) = 1;
                    elseif x(i) > const_width
                        y(i) = 0;
                    end
                end
                w = y;
                w(si) = y;                                  % Reorder
                fid = fopen('RespWeight_soft.dat','w');     % Write to file
                fwrite(fid,w,'single');
                fclose(fid);
                
            case 'Sigmoid'
            myPhs = mod(SoftGating_Phs + 0.5 - (fr - 1), 8);
            f = zeros(size(myPhs));
            dn = abs(myPhs-0.5);
            dn(dn>4) = 8 - dn(dn>4);
            for i = 1:numel(myPhs)
                if dn(i) > const_width 
                    f(i) = exp(-1.*sigma.*(dn(i) - const_width)^2);
                else
                    f(i) = 1;
                end
            end
            f = f-min(f);
            f = f.*(1/max(f));
            w = f;
            w(si) = f;                                  % Reorder
            fid = fopen('RespWeight_soft.dat','w');     % Write to file
            fwrite(fid,w,'single');
            fclose(fid);
            
            case 'Inverse_Exp'
            myPhs = mod(SoftGating_Phs + 0.5 - (fr - 1), 8);
            f = zeros(size(myPhs));
            dn = abs(myPhs-0.5);
            dn(dn>4) = 8 - dn(dn>4);
            alpha = 2;
            for i = 1:numel(myPhs)
                if dn(i) > const_width 
                    f(i) = (-1.*exp(alpha.*((dn(i) - const_width))));
                else
                    f(i) = 1;
                end
            end
            f = f-min(f);
            f = f.*(1/max(f));
            w = f;
            v = (w - sigma)/(1-sigma);
            v(v<=0)=0;
            v(si) = f;                                  % Reorder
            fid = fopen('RespWeight_soft.dat','w');     % Write to file
            fwrite(fid,v,'single');
            fclose(fid);
        end

        % Run the recon for this frame
        system(recon_cmd);

        % Move the image into the recon folder
        mv_cmd = sprintf('mv X_000_000.dat recon/%s_%s_CW-%0.2f_Sigma-%0.2f_X_000_0%02d.dat',Patient,functype,const_width,sigma,fr-1); 
        system(mv_cmd); 
        formatSpec = '%s_%s_CW-%0.2f_Sigma-%0.2f_Frame%02d';
        str = sprintf(formatSpec,Patient,functype,const_width,sigma,fr-1);
        mkdir(str);
        movefile('*.dcm',str);
        %movefile('*.txt',str);
    end
    
