classdef pTXObj < handle
    
    % base class for pTX pulses
    % author : Vincent Gras
    % contact : vincent.gras@cea.fr

    properties
        
        Nc = 1                       % nb of tx channels
%        FA = 0                      % flip angle deg
        %refB1 = 0;                   % a reference B1 in Tesla (0 if undefined), typically the B1 at the center in CP mode for 1V applied
        refAmplitude = 1             % reference amplitude
        refAmplitudeUnit = 'W^1/2'   % Unit for the ref amplitude 'Volt' or  W^1/2 ....
        gyr = pTXUtils.gyromagneticRatio('1h');                      % Gyromagnetic ratio Hz/Tesla
        nominalb0 = 7                % Nominal B0 Tesla
        %gradAxes = eye(3)            % Direction for B0 gradients in projection on X,Y,Z
        %master_ = 0;                 % Master index
        flipped_rfphase_ = false;    % sign of RF phase
        subject_orientation_ = 'HFS'; % assume subject orientation = head first supine !
        XYZ = 'LPH';
    end
    
    methods
        
        function o = pTXObj (varargin) %Nchannels, MagneticField)
            
            
            if (nargin  == 0)
                return;
            end
            
            i = 1;
                
            if (isa(varargin{i}, 'pTXObj'))
                
                f = properties('pTXObj');
                
                for n = 1:numel(f)
                    o.(f{n}) = varargin{i}.(f{n});
                end
                
                i = i+1;
                
            end
            
            o.setProperties(varargin(i:end));
            
%             while (i < nargin)                
%                 assert(ischar(varargin{i}), 'Bad input #%d, expecting character array', i);
%                 assert(i< nargin, 'property %s is not associated with a value !', varargin{i});
%                 o.setProperty(varargin{i:i+1});
%                 i = i + 2;                
%             end   
%             
%             
%         
%             remainarg = [];
%             par = struct();
%             
%             if (isstruct(varargin{1}))
%                 
%                 assert(nargin == 1, 'bad input');
%                 remainarg = [];
%                 par = varargin{1};
%                 
%                 
%             elseif (ischar(varargin{1}))
%                 
%                 if (nargin > 1)
%                     [par,remainarg] = readConstrInput(varargin);
%                 end
%             elseif (isstruct(varargin{1}))
%                 par = varargin{1};
%             elseif (isa(varargin{1}, 'pTXObj'))
%                 
%                 f = properties('pTXObj');
%                 
%                 pTXObj.setProperties(o, varargin{1}, f);
% %                 for i = 1:numel(f)
% %                     o.(f{i}) = varargin{1}.(f{i});
% %                 end
%                 
%                 [par,remainarg] = readConstrInput(varargin(2:end));
%                 
%             else
%                 error('invalid input');
%             end
%             
%             
%             if (~isempty(remainarg))
%                 disp(remainarg)
%                 error('Some args have not been recognized (see above)');
%                 
%             end
%             
%             pTXObj.setProperties(o, par);
%             
        end
        
        function setProperties(o, l)
           
            i = 1;
            while (i <= numel(l))
                
                par = l{i};
                assert(ischar(par), 'Bad input #%d, expecting character array', i);
                assert(~isempty(par), 'Bad input #%d, expecting non-empty character array', i);
                
                switch par(1)
                    case '?'
                        o.setProperty(par(2:end), true);
                    case '!'
                        o.setProperty(par(2:end), false);
                    otherwise   
                        %assert(i+1<=numel(l), 'Input %d : Missing value for %s', i, par);
                        
                        if (i+1<=numel(l))
                            o.setProperty(par, l{i+1});
                            i = i+1;
                        end
                end
                i = i + 1;                
            end   

        end
        
        function success = setProperty(o, name, val)
            
            success = true;
            
            switch lower(name)
                
                case 'nc'
                    o.setNc(val);
                case 'nominalb0'
                    o.nominalb0(:) = val;
                case 'gyr'
                    o.gyr(:) = val;
                case 'nucleus'
                    o.gyr(:) = pTXUtils.gyromagneticRatio(val);
                case 'refamplitude'
                    assert(iscell(val)&& numel(val) == 2, 'bad value for property ''refamplitude'', expecting {value, unit} pair');
                    o.refAmplitude(:) = val{1};
                    o.refAmplitudeUnit = val{2};
%                 case 'fa'
%                     if (isnumeric(val))
%                         o.FA(:) = val;
%                     else
%                         assert(iscell(val) && numel(val) == 2, 'bad value for property ''fa'', expecting {value, unit} pair');
%                         
%                         if (strcmpi(val{2}, 'rad'))
%                             val{1} = val{1} * 180/pi;
%                         end
%                         o.FA(:) = val{1}; 
%                     end
                case 'xyz'
                    o.XYZ(:) = val;
                otherwise
                    success = false;
            end

            
        end
        
        function g = absgyr(o)
            g = abs(o.gyr);
        end
        
        function g = neggyr(o)
            g = -(o.gyr); 
        end
        
        function o = copy(o, ori)
            
            o.setProperties(ori);
            
        end
        
        function str = tostruct(o)
            
            f = properties(o);
            
            str = struct();
            
            for i = 1:numel(f)
                str.(f{i}) = o.(f{i});
            end
            
            
            
            
        end
        
        function p = setNc(p, nch)
            
            assert(nch>0, 'number of channels must be > 0')
            p.Nc(:) = nch;
            
        end
        
        function a = getRefAmplitude(o,unit)
           a = o.refAmplitude;
           
           if (nargin < 2)
               return;
           end
           
           assert(strcmpi(o.refAmplitudeUnit, unit), 'not automatic conversion (yet) !!');
            
            
        end
        
        function o = setRefAmplitude(o, A, unit)
            o.refAmplitude = A;
            if (nargin >2)
                o.refAmplitudeUnit = unit;
            end
        end
        
        function val = X(o)
            val = getAxis(o.XYZ(1));
        end
        
        function val = Y(o)
            val = getAxis(o.XYZ(2));
        end
        
        function val = Z(o)
            val = getAxis(o.XYZ(3));
        end
        
        function U = sRotMatUnity(o)
            
            % sRotMat unity means : XYZ = PRH
            % Compute U such that
            % P = U11 o.X + U21 o.Y + U31 o.Z
            % R = U12 o.X + U22 o.Y + U33 o.Z
            
            U = zeros(3);
            
            for i = 1:3
                switch (upper(o.XYZ(i)))
                    case 'A'
                        U(i,1) = -1;
                    case 'P'
                        U(i,1) = 1;
                    case 'R'
                        U(i,2) = 1;
                    case 'L'
                        U(i,2) = -1;
                    case 'H'
                        U(i,3) = 1;
                    case 'F'
                        U(i,3) = -1;
                    otherwise
                        error('unrecognized direction %s', upper(o.XYZ(i)));
                end
            end
            
        end
        
    end
    
    
    methods (Static)
        
    end
    
end

function val = getAxis(val)

switch (val)
    case 'A'
        val = 'P->A';
    case 'P'
        val = 'A->P';
    case 'R'
        val = 'L->R';
    case 'L'
        val = 'R->L';
    case 'H'
        val = 'F->H';
    case 'F'
        val = 'H->F';
    otherwise
        error('unrecognized direction %s', val);
end

end