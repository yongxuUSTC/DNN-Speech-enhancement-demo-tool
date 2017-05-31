  %% ************************************************************************
        % writeHTK - convert data to htk format -> by anonymous Sue (2001)
        %**************************************************************************
        function retcode = writeHTK_new(filename, htkdata, nFrames, sampPeriod, SampSize, ParamKind, byte_order)
            % Write an HTK format file.
            %
            % Input parameters:
            %    filename		HTK data file
            %    htkdata      HTK data read: an m x n matrix with
            %                    m = no. of channels
            %                    n = no. of frames
            %  The following are from the HTK header (see HTK manual):
            %    nFrames      no. of frames (samples)
            %    sampPeriod   sample period (in 100 ns units?)
            %    SampSize     sample size
            %    ParamKind    parameter kind code
            %
            %    byteorder    'be' for big-endian (typical for Unix) (default)
            %                 'le' for little-endian (typical for MSWindows)
            %
            % Output parameters:
            %    retcode      0 if successful
            
            % Written by Sue 17/12/01
            
            retcode=-1;	% initialise in case of error
            if nargin < 6
                fprintf('Usage: %s(filename, htkdata, nFrames, sampPeriod, SampSize, ParamKind [, byte_order])', mfilename);
            end;
            
            % Default to big-endian (HTK format)
            if nargin < 7
                byte_order = 'be';
            end;
            
            fid = fopen (filename, 'w', sprintf('ieee-%s', byte_order));
            if fid < 1
                fprintf('%s: can''t open output file %s\n', mfilename, filename);
                return
            end
            
            % Write header
            fwrite (fid, nFrames, 'int32'); %nSamples in HTK
            fwrite (fid, sampPeriod, 'int32');
            fwrite (fid, SampSize, 'int16');
            fwrite (fid, ParamKind, 'int16');
            
            % Write actual data
            fwrite(fid, htkdata', 'float32');
            
            fclose(fid);
            
            retcode=0;
        end% ------ OF WRITEHTK
        