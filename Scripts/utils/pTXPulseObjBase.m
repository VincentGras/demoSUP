classdef pTXPulseObjBase < pTXObj
    
    % Matlab class for handling pTX pulses
    % author : Vincent Gras
    % contact : vincent.gras@cea.fr

    properties
        
        
        FA = 0                       % flip angle deg
        N = 0;
        dt = []
        dRF = []        % in s
        dk = []         % in rad/s defined as
        minRasterTime = 1e-6;        % Min raster time
        rasterTime = inf;
        
    end
    
    methods
        
        function p = pTXPulseObjBase(varargin) % Nch, Npoints, dt, varargin)
            
            p@pTXObj(varargin{1:min(nargin,1)});
            
            if (nargin  == 0)
                return;
            end
            
            i = 1;
            
            if (isa(varargin{i}, 'pTXObj'))
                
                % copy constructor 
                
                if (isa(varargin{i}, 'pTXPulseObjBase'))
                    
                      assert(varargin{i}.Nc == p.Nc, 'Inconsistent Nc !');
                    
                      p.minRasterTime = varargin{i}.minRasterTime;
                      p.rasterTime = varargin{i}.rasterTime;
                      p.setN(varargin{i}.N);
                      p.dt(:) = varargin{i}.dt;
                      p.dRF(:) = varargin{i}.dRF;
                      p.dk(:) = varargin{i}.dk;
                      p.FA(:) =  varargin{i}.FA;
                    
%                     f = properties(varargin{i});
%                     f = setdiff(f, properties('pTXObj'));
%                     f = intersect(f, properties(p));
%                     
%                     for n = 1:numel(f)
%                         p.(f{n}) = varargin{i}.(f{n});
%                     end
                end
                
                i = i+1;
                
            end
            
            
            p.setProperties(varargin(i:end));

%             while (i < nargin)
%                 
%                 assert(ischar(varargin{i}), 'Bad input #%d, expecting character array', i);
%                 assert(i< nargin, 'property %s is not associated with a value !', varargin{i});
%                 assert(p.setProperty(varargin{i:i+1}), 'unknown property ''%s''', varargin{i});
%                 i = i + 2;
%                 
%             end
            
            %             p@pTXObj(varargin{1:min(nargin,1)});
            %
            %             if (nargin < 1)
            %                 return;
            %             end
            %
            %             if (isstruct(varargin{1}))
            %
            %                 assert(nargin == 1, 'bad input');
            %                 remainarg = [];
            %                 par = varargin{1};
            %
            %
            %             elseif (ischar(varargin{1}) && nargin == 1)
            %
            %                 remainarg = [];
            %                 par = struct();
            %
            %             elseif (ischar(varargin{1}) && nargin > 1)
            %
            %                 [par,remainarg] = readConstrInput(varargin);
            %
            %             elseif (isa(varargin{1}, 'pTXObj'))
            %
            %                 f = properties(varargin{1});
            %                 f = setdiff(f, properties('pTXObj'));
            %                 f = intersect(f, properties(p));
            %
            %                 pTXPulseObjBase.setProperties(p, varargin{1}, f);
            % %                 for i = 1:numel(f)
            % %                     p.(f{i}) = varargin{1}.(f{i});
            % %                 end
            %
            %
            %                 [par,remainarg] = readConstrInput(varargin(2:end));
            %
            %
            %             else
            %
            %                 error('bad input');
            %
            %             end
            %
            %
            %             if (~isempty(remainarg))
            %
            %                 disp(remainarg)
            %                 error('Some args have not been recognized (see above)');
            %
            %             end
            %
            %             pTXPulseObjBase.setProperties(p, par);
            %
        end
        
        function success = setProperty(p, name, val)
            
            success = p.setProperty@pTXObj(name, val);
            
            if (success)
                return;
            end
            
            success = true;
            
            switch lower(name)
                case 'fa'
                    if (isnumeric(val))
                        p.FA(:) = val;
                    else
                        assert(iscell(val) && numel(val) == 2, 'bad value for property ''fa'', expecting {value, unit} pair');
                        
                        if (strcmpi(val{2}, 'rad'))
                            val{1} = val{1} * 180/pi;
                        end
                        p.FA(:) = val{1};
                    end
                case 'n'
                    
                    p.setN(val);
                    
                case 'drf'
                    
                    p.dRF(:) = val;
                    
                case 'dk'
                    
                    p.dk(:) = val;
                    
                case 'dt'
                    
                    p.dt(:) = val;
                    
                case 'rastertime'
                    
                    p.rasterTime(:) =  val;
                    
                otherwise
                    success = false;
                    
            end
            
            
        end
        
        function p = reset(p)
            
            p.dt = zeros(1, p.N);
            p.dRF = zeros(p.Nc, p.N);
            p.dk = zeros(3, p.N);
            
        end
        
        function p = setN(p, n)
            
            assert(isscalar(n) && rem(n, 1)==0, 'Bad parameter N');
            if (n ~= p.N)
                p.N = n;
                p.reset();
            end
            
        end
        
        function p = setNc(p, nc)
            
            assert(isscalar(nc) && rem(nc, 1)==0, 'Bad parameter Nc');
            if (nc ~= p.Nc)
                p.Nc = nc;
                p.dRF = zeros(p.Nc, p.N);
            end
        end
        
        function p = setNcN(p, nc, n)
            
            assert(isscalar(nc) && rem(nc, 1)==0, 'Bad parameter Nc');
            assert(isscalar(n) && rem(n, 1)==0, 'Bad parameter N');
            if (nc ~= p.Nc || n ~= p.N)
                p.Nc = nc;
                p.N = n;
                p.reset();
            end
            
        end
        
        function q = new(p)
            
            q = pTXPulseObjBase(p);
            
        end
        
        function q = sampled(p, varargin)
            
            [gw,rfw,Dt] = p.getWaveform(varargin{:});
            
            q = pTXPulseObjBase(p);
            q.rasterTime = Dt(1);
            q.setN(numel(Dt));
            q.dt(:) = Dt;
            q.dRF(:,:) = (rfw/q.refAmplitude) .* q.dt;
            q.dk(:,:) = (gw * (2*pi*q.neggyr)) .* q.dt;
            
            
        end
        
        function o = newpTXObj(p)
            
            o = pTXObj();
            f = properties('pTXObj');
            
            for i = 1:numel(f)
                o.(f{i}) = p.(f{i});
            end
            
        end
        
        function str = tostruct(p)
            
            str = tostruct@pTXObj(p);
            str.N = p.N;
            str.dt = p.dt;
            str.dRF = p.dRF;
            str.dk = p.dk;
            
        end
        
        function  rasterTime = getMaxRasterTime(p)
            
            epsi = p.minRasterTime; % by default 0.1 us
            
            rasterTime = round([p.dt] / epsi);
            
            for i = 2:numel(rasterTime)
                rasterTime(1) = gcd(rasterTime(1), rasterTime(i));
            end
            
            if (isempty(rasterTime))
                return;
            end
            
            rasterTime = rasterTime(1)*epsi;
            
        end
        
        function p = interp(p, dtnew)
            
            [dk_, dtnew] = pTXUtils.interp1dsignal(p.dt, p.dk, dtnew);
            dRF_ = pTXUtils.interp1dsignal(p.dt, p.dRF, dtnew);
            
            p.setN(numel(dtnew));
            p.dt(:) = dtnew;
            p.dRF(:,:) = dRF_;
            p.dk(:,:) = dk_;
            
        end
        
        function t = getTimeVector(p, arg) % cumsum(dt) - 0.5 * dt
            
            if (nargin > 1)
                assert(isscalar(arg) && arg >= 0 && arg <= 1, 'expected scalar as input (in [0,1])');
                %    t = t - arg * p.dt;
            else
                arg = 0.5;
            end
            
            t = cumsum(p.dt) - (1-arg) * p.dt;
            
        end
        
        function fa = getFA(o, unit)
            
            
            if (nargin < 2)
                unit = 'deg';
            end
            
            switch lower(unit)
                case 'deg'
                    fa = o.FA;
                case 'rad'
                    fa = o.FA*pi/180;
                otherwise
                    error('unknown unit ''%s''', unit);
            end
        end

        function t = getRFDur(p)
            i = p.rfsupport();
            t = sum(p.dt(i(1):i(end)-1));
        end

        
        function t = getBaseDur(p)
            t = sum(p.dt);
        end
        
        function t = getDur(p)
            t = p.getBaseDur();
        end
        
        function [Etot, E] = getEnergy(p)
            
            i = any(abs(p.dRF)>0,1);
            E = p.refAmplitude^2 * bsxfun(@rdivide, abs(p.dRF(:,i)).^2, p.dt(i));
            E = sum(E, 2);
            Etot = sum(E);
            
        end
        
        function [E, Etot] = getEnergyPerChannel(p)
            
            [Etot, E] = p.getEnergy();
            
        end
        
        function [P,Pmax] = getMaxRFPowerPerChannel(p)
            
            i = any(abs(p.dRF)>0,1);
            P = p.refAmplitude^2 * max(abs(bsxfun(@rdivide, p.dRF(:,i), p.dt(i))).^2, [], 2);
            Pmax = max(P);
            
        end
        
        function SED = getSED(p, VOP)
            
            if (ismatrix(VOP))
                VOP = reshape(VOP, size(VOP, 1), size(VOP, 1), size(VOP, 2) / size(VOP, 1));
            end
            
            % VOP dimensions = Nc x Nc x NVOP
            assert (size(VOP,1) == p.Nc && size(VOP,2) == p.Nc, 'Bad VOP dimensions, expecting %d x %d x NVOP', p.Nc, p.Nc);
            
            NVOP = size(VOP, 3);
            SED = zeros(NVOP, 1);
            
            %             for i = 1:NVOP
            %                 for j = 1:p.N
            %                     RFc = p.dRF(:,j);
            %                     if (any(RFc))
            %                         SED(i) = SED(i) + real(RFc' * VOP(:,:,i) * RFc) / p.dt(j);
            %                     end
            %                 end
            %             end
            %
            %             SED_test = zeros(NVOP, 1);
            
            for j = 1:p.N
                RFc = p.dRF(:,j);
                if (any(RFc))
                    E = real(pagemtimes(RFc, 'ctranspose', pagemtimes(VOP, RFc), 'none')) / p.dt(j);
                    SED = SED + E(:);
                end
            end
            
            
            SED = SED * p.refAmplitude^2;
            
        end
        
        function [gw, rfw, Dt] = getWaveformBase(p, Dt)
            
            if (nargin < 2)
                Dt = p.rasterTime;
            end
            
            p.checkRasterTime(Dt);
                        
            if ((~isscalar(Dt) && numel(p.dt) ~= numel(Dt)) || any(abs(p.dt - Dt)>eps))
                p = p.new();
                p.interp(Dt);                 
            end
            
            gw = (p.dk ./ p.dt) / (2*pi*p.neggyr());
            rfw = (p.dRF ./ p.dt) .* p.refAmplitude;
            
            Dt = p.dt;
        end

        
        function [gw, rfw, Dt] = getWaveform(p, Dt)
            
            if (nargin < 2)
                Dt = p.rasterTime;
            end
            
            [gw, rfw, Dt] = p.getWaveformBase(Dt);
        
        end
        
        function p = multiplex(p, ch)
            
            p.dRF = p.dRF(ch, :);
            p.Nc = numel(ch);
            
        end
        
        function p = calibrate(p, b1, w)
            
            tx = mean(p.dRF, 2);
            
            assert(any(tx), 'tx configuration == 0 !!');
            
            if (p.getFA('rad') == 0)
                warning ('zero flip angle')
            end
            
            if (nargin < 3)
                w = 1;
            end
            
            b1Avg = mean(abs(sum(b1 .* tx,1)).* w) / mean(w);
            p.dRF = p.dRF * (p.getFA('rad') / b1Avg / p.refAmplitude); 
            
        end

        
        function rt = defaultRasterTime(p)
            rt = p.rasterTime;
        end
        
        function visualizeMP(p, rasterTime_s)
            
            if (isempty(p.dRF))
                warning('empty pulse');
                return;
            end
            
            if (nargin < 2)
                rasterTime_s = p.defaultRasterTime();
            end
            
            [kw, rfw, deltat] = p.getWaveform(rasterTime_s);
            
            t = [0, cumsum(deltat)];
            
            %figure();
            clf;
            
            %             if (numel(t)<= 1)
            %                 warning ('only time point --> nothing relevant to display');
            %                 return;
            %             end
            
            subplot(3,1,1);
            stairs(t*1e3, 1e3*cat(2, kw, kw(:, end))');
            xlabel('t (ms)')
            ylabel('G (mT/m)')
            hold on;
            %ylim([-20,20])
            subplot(3,1,2);
            stairs(t, cat(2, abs(rfw), abs(rfw(:, end)))');
            xlabel('t (ms)')
            ylabel(sprintf('abs(RF) (%s)', p.refAmplitudeUnit))
            subplot(3,1,3);
            stairs(t, cat(2, angle(rfw), angle(rfw(:, end)))');
            set(gca, 'ylim', [-pi,pi]);
            xlabel('t (ms)')
            ylabel(sprintf('angle(RF) (%s)', 'radian'));
            
            
        end

        
        function visualize(p, rasterTime_s)
            
            if (isempty(p.dRF))
                warning('empty pulse');
                return;
            end
            
            if (nargin < 2)
                rasterTime_s = p.defaultRasterTime();
            end
            
            [kw, rfw, deltat] = p.getWaveform(rasterTime_s);
            
            t = [0, cumsum(deltat)];
            
            %figure();
            clf;
            
            %             if (numel(t)<= 1)
            %                 warning ('only time point --> nothing relevant to display');
            %                 return;
            %             end
            
            subplot(3,1,1);
            stairs(t*1e3, 1e3*cat(2, kw, kw(:, end))');
            xlabel('t (ms)')
            ylabel('G (mT/m)')
            hold on;
            %ylim([-20,20])
            subplot(3,1,2);
            stairs(t, cat(2, real(rfw), real(rfw(:, end)))');
            xlabel('t (ms)')
            ylabel(sprintf('Re(RF) (%s)', p.refAmplitudeUnit))
            subplot(3,1,3);
            stairs(t, cat(2, imag(rfw), imag(rfw(:, end)))');
            xlabel('t (ms)')
            ylabel(sprintf('Im(RF) (%s)', p.refAmplitudeUnit));
            
            
        end
        
        function p = backProjectK(p, T)
            
            p.dk = T * p.dk;
            
        end
        
        function Q = blochsimQfullHistory(p, positions_m, b1, b0)
            
            % Computes propagator
            % position_m (3xNpos)  : voxel positions            [m]
            % b1         (NcxNpos) : B1 field in                [rad/s/<reafAmplitudeUnit>]
            % b0         (1xNpos)  : static field offset in Hz  [rad/s]
            %
            % Note : make sure that the [p.refAmplitude . b1] == [rad/s]
            
            
            if (size(positions_m, 1) ~= 3 && size(positions_m, 2) == 3)
                warning ('transposing positions_m');
                positions_m = positions_m';
            end
            
            if (size(b1, 1) ~= p.Nc && size(b1, 2) == p.Nc)
                warning ('transposing b1');
                b1 = b1.';
            end
            
            if (isscalar(b0))
                b0 = repmat(b0, 1, size(b1, 2));
            end
            
            assert(ndims(b0) <= 2 && numel(b0) == max(size(b0)), 'b0 must be a line or row vector');
            
            assert(size(positions_m, 1) == 3, 'positions_m must be 3xNpos');
            assert(size(b1, 1) == p.Nc, 'b1 must be %dxNpos', p.Nc);
            assert(size(b1, 2) == size(positions_m, 2), 'positions_m says Npos=%d but size(b1,2) = %d', size(positions_m,2), size(b1, 2));
            assert(numel(b0) == size(positions_m, 2), 'positions_m says Npos=%d but size(b1,2) = %d', size(positions_m,2),numel(b0));
            
            
            %Q = pTXUtils.RunBlochSimFullHistory(positions_m, b1, p.b0Sign_() * b0, p.dRF * p.refAmplitude, p.dk, p.dt);
            Q = pTXUtils.RunBlochSimFullHistory(positions_m, b1, b0, p.dRF * p.refAmplitude, p.dk, p.dt);
            Q = reshape(Q, 4, numel(p.dt), numel(positions_m)/3);
            
            
        end
        
        function Q = blochsimQ(p, positions_m, b1, b0)
            
            % Computes propagator
            % position_m (3xNpos)  : voxel positions            [m]
            % b1         (NcxNpos) : B1 field in                [rad/s/<reafAmplitudeUnit>]
            % b0         (1xNpos)  : static field offset in Hz  [rad/s]
            %
            % Note : make sure that the [p.refAmplitude . b1] == [rad/s]
            
            
            if (size(positions_m, 1) ~= 3 && size(positions_m, 2) == 3)
                warning ('transposing positions_m');
                positions_m = positions_m';
            end
            
            if (size(b1, 1) ~= p.Nc && size(b1, 2) == p.Nc)
                warning ('transposing b1');
                b1 = b1.';
            end
            
            if (isscalar(b0))
                b0 = repmat(b0, 1, size(b1, 2));
            end
            
            assert(ndims(b0) <= 2 && numel(b0) == max(size(b0)), 'b0 must be a line or row vector');
            
            assert(size(positions_m, 1) == 3, 'positions_m must be 3xNpos');
            assert(size(b1, 1) == p.Nc, 'b1 must be %dxNpos', p.Nc);
            assert(size(b1, 2) == size(positions_m, 2), 'positions_m says Npos=%d but size(b1,2) = %d', size(positions_m,2), size(b1, 2));
            assert(numel(b0) == size(positions_m, 2), 'positions_m says Npos=%d but size(b1,2) = %d', size(positions_m,2),numel(b0));
            
            %Q = pTXUtils.RunBlochsim(positions_m, b1, p.b0Sign_() * b0, p.dRF * p.refAmplitude, p.dk, p.dt);
            Q = pTXUtils.RunBlochsim(positions_m, b1, b0, p.dRF * p.refAmplitude, p.dk, p.dt);
            
        end
        
        function q = clip(p)
            
            q = p.new();
            Ampl_ratio = abs(q.dRF./q.dt);
            Ampl_ratio(Ampl_ratio>1) = 1;
            q.dRF = Ampl_ratio.*q.dt.*exp(1i*angle(q.dRF));
            
        end
        
        function FA = blochsim(p, positions_m, b1, b0)
            
            % Computes flip angle
            % position_m (3xNpos)  : voxel positions            [m]
            % b1         (NcxNpos) : B1 field in                [rad/s/<reafAmplitudeUnit>]
            % b0         (1xNpos)  : static field offset in Hz  [rad/s]
            %
            % Note : make sure that the [p.refAmplitude . b1] == [rad/s]
            
            Q = p.blochsimQ( positions_m, b1, b0);
            %R = qtn__RotMatrix(Q);
            %FA = acos(R(3,3,:));
            %FA = reshape(FA, 1, numel(FA));
            FA = pTXUtils.qtn__FA(Q);
        end
        
        function p = setReferenceVoltage(p, refV)
            p.dRF = p.dRF*p.refAmplitude/refV;
            p.refAmplitude = refV;
        end
        
        function p = setRasterTime(p, rasterT)
            
            %p.checkRasterTime(rasterT);
            
            p.rasterTime(:) = rasterT;
            
        end
        
        function p = checkRasterTimeRF(p, rt)
            
            if (nargin < 2)
                rt = p.rasterTime;
            end
            
            assert(ndims(rt) <= 2 && size(rt,1) == 1, 'raster time must be 1xn');
            assert(all(rt < inf), 'raster time was not set (current value is infty)');
            assert(all(rt >= p.minRasterTime), 'raster time cannot be set < %.2f us (minRasterTime)', p.minRasterTime * 1e6);
            
        end
        
        function p = checkRasterTime(p, varargin)
            
            p.checkRasterTimeRF(varargin{:});
            
        end
        
        function p = setFA(p, newFAdeg, unit)
            
            if (nargin < 2)
                return;
            end
            
            if (nargin > 2)
                
                if (strcmpi(unit, 'rad'))
                    newFAdeg = newFAdeg * 180/pi;
                end
                
            end
            
            
            if (p.FA>0)
                p.dRF(:) = p.dRF(:) * (newFAdeg/p.FA);
                %p = p.setRFCoeff (p.getRFCoeff() * (newFAdeg/p.FA));
            end
            
            p.FA = newFAdeg;
            
        end
        
        function o = setRefAmplitude(o, A, unit)
            
            if (nargin < 3)
                unit = o.refAmplitudeUnit;
            end
            
            VrefOld = o.getRefAmplitude(unit);
            o.setRefAmplitude@pTXObj(A, unit);
            VrefNew = o.getRefAmplitude(unit);

            o.dRF = o.dRF * VrefOld / VrefNew;
            
        end

        function p = split(p, n)
           
            if (p.N<=0)
                return;
            end
            
            p.N = n*p.N;
            p.dRF = kron(p.dRF/n, ones(1,n));
            p.dk = kron(p.dk/n, ones(1, n));
            p.dt = kron (p.dt/n, ones(1, n));
            
        end
        
        function p = doubleangle(p)
            
            dt_ = p.dt;
            dk_ = p.dk;
            dRF_ = p.dRF;
            
            p.setN(2*p.N);
            
            p.FA = 2 * p.FA;
            
            p.dt(:) = [dt_, fliplr(dt_)];
            p.dk(:,:) = [dk_, -fliplr(dk_)];
            p.dRF(:,:) = [dRF_, fliplr(dRF_)];
            
        end
        
        function p = inverse(p)
            
            p.dRF = -fliplr(p.dRF);
            p.dk = -fliplr(p.dk);
            
        end
        
        function p = scaletime(p, s)
            
            p.dt(:) = p.dt * s;
            
        end
        
        
        function p = phaseshift(p, ph)
            p.dRF = p.dRF * exp(1i*ph);            
        end
        
        function p = compose(p, q, varargin)
            
            assert(strcmp(p.refAmplitudeUnit, p.refAmplitudeUnit), 'p and q do not have the same ref ampl unit');
            assert(p.Nc == q.Nc, 'p and q do not have the same number of channels');
            %            assert(norm(p.deltakFinal + q.deltakInitial) <= 2*eps, 'cannot compose those two pulses');
            
            refAmplp = p.refAmplitude;
            refAmplq = q.refAmplitude;
            
            p.setFA(varargin{:}); 
            p.N = p.N + q.N;
            p.dt = [p.dt,q.dt];
            p.refAmplitude = max(refAmplp,refAmplq);
            p.dRF = [p.dRF * refAmplp, q.dRF * refAmplq]/p.refAmplitude;
            p.dk = [p.dk, q.dk];
            
        end
        
        
        function k = getKspaceTraj(p, arg)
            
            if (nargin < 2)
                arg = 0.5;
            end
            
            assert(isscalar(arg) && arg >= 0 && arg <= 1, 'bad input');
            
            k = cumsum(p.dk, 2) - (1-arg) * p.dk;
            k = k(:,end) - k; % !!! 9/7/21
            
        end
        
        function FA = stasim(p, pos, b1, b0, cmplx)
            
            if (nargin < 5)
                cmplx = false;
            end
            
            FA = p.flipAngleOper(pos, b0);
            
            if (cmplx)
                FA = sum(b1.*FA, 1);
            else
                FA = abs(sum(b1.*FA, 1));
            end
        end
        
        function F = flipAngleOper(p, pos, b0)
            
            % Compute first "flip angle" operator : FAop : b1 -> FA
            
            % verify dimensions of b1 and pos
            
            if (size(pos, 1) ~= 3 && size(pos, 2) == 3)
                warning ('transposing pos');
                pos = pos';
            end
            
            Npos = numel(pos)/3;
            
            % check for b0
            
            b0exists = (nargin > 2) && ~isempty(b0);
            
            if (b0exists && ~isscalar(b0))
                assert(numel(b0) == Npos, 'Bad b0 input, expected 1x%d array', Npos);
                assert(size(b0,2) == Npos, 'Bad b0 input, expected 1x%d array', 1, Npos);
            end
            
            deltat = p.dt';
            %t = cumsum(deltat)-0.5*deltat; % t = p.getTimeVecto(p.N)';
            %T = sum(deltat);
            t = p.getTimeVector()';
            T = p.getDur();            
            k = p.getKspaceTraj()';
            %u = cat(2, p.dRF * p.refAmplitude, zeros(p.Nc, 1)); % Nc x N+1
            u = p.dRF * p.refAmplitude;
            u = permute(u, [2, 3, 1]); % N x 1 x Nc
            
            F = pos(1,:) .* k(:,1) + pos(2,:) .* k(:,2) + pos(3,:) .* k(:,3);
            
            if (b0exists)
                %F = F + b0 .* (p.b0Sign_()*(T - t)); % N x Npos x 1    (T-t ---> t-T)
                F = F + b0 .* (T - t); % N x Npos x 1 % !!! 9/7/21
            end
            
            F = ((1i)*u) .* exp(-1i*F); % N x Npos x Nc % !!! 9/7/21
            
            if (b0exists)
                F = F .* pTXUtils.sinc(0.5 * b0 .* deltat);
            end
            
            F = sum(F, 1);   %sum(F, 1); % 1 x Npos x Nc
            F = permute(F, [3, 2, 1]); % Nc x Npos
            
        end
        
        function f = exportToIniFileVB17(p, filename, varargin)
            
            par.TimeStep = 1e-5;
            par.MaxAbsRF = 0;
            par.Asymmetry = 0.5;
            par.InitialPhase = 0;
            par.Comment = 'created by pTXPulseObjBase.exportToIniFileVB17';
            par.RotMat = eye(3);
            par.SliceLocation = [];
            
            par = pTXUtils.readopt(par, varargin{:});
            
            [gw, rfw] = p.getWaveform(par.TimeStep);
            
            par.NominalFlipAngle = p.FA;
            
            master_ = 0;
            rfw = circshift(rfw, -(master_ - 1), 1); % ensure that channelIndex(1) == 0 [Nc]
            
            rfw = conj(rfw); % !!! 9/7/21
            
            f = pTXUtils.makeIniFileVB17(filename, rfw, 1e3 * (p.sRotMatUnity()*par.RotMat)' * gw, par);
            
        end
        
        
        function [s, S] = rfsupport(p)
            
            S = sum(abs(p.dRF), 1)>0;
            s = pTXPulseObjBase.supp(S);
            
            
        end
        
        function [s,S] = gradientsupport(p)
            
            S = sum(abs(p.dk), 1)>0;
            s = pTXPulseObjBase.supp(S);
            
        end
        
        function s = support(p)
            
            s = pTXPulseObjBase.supp(((sum(abs(p.dk), 1))>0) | (sum(abs(p.dRF), 1)>0));
            
        end
        
        function q = cut2supp(q, s)
            
            if (nargin < 2)
                s = q.support(); % 0 if support is void
            end
            
            if (s(1)>1 || s(end) <= q.N)
                j = s(1):s(end)-1;
                q.N = numel(j);
                q.dRF = q.dRF(:,j);
                q.dk = q.dk(:,j);
                q.dt = q.dt(j);
            elseif (s(1) == 0)
                q.setN(0);
            end
            
        end
        
        function [q] = cut2rfsupp(q)
            
            q.cut2supp(q.rfsupport());
            
        end
        
        
        function q = zpadding(q, pre, post)
            
            dtpre = [];
            dtpost = [];
            RFpre = zeros(q.Nc,0);
            RFpost = zeros(q.Nc, 0);
            dkpre = zeros(3,0);
            dkpost = zeros(3,0);
            
            if (pre>0)
                q.N = q.N + 1;
                dtpre = pre;
                RFpre = zeros(q.Nc,1);
                dkpre = zeros(3,1);
            end
            if (post>0)
                q.N = q.N + 1;
                dtpost = post;
                RFpost = zeros(q.Nc,1);
                dkpost = zeros(3,1);
            end
            
            
            q.dt = [dtpre, q.dt, dtpost];
            q.dRF = cat(2, RFpre, q.dRF, RFpost);
            q.dk = cat(2, dkpre, q.dk, dkpost);
            
            
            
        end
        
        function rt = getRasterTime(q)
            rt = unique(q.dt);
        end
        
        
    end
    
    methods (Static)
        
        function p = HardPulse(dur_s, refAmp, unit)
            
            p = pTXPulseObjBase('Nc', 1, 'N', 1, 'dt', dur_s, 'refamplitude', {refAmp, unit});
            p.dRF(1,:) = dur_s;
            
        end
        
        function p = HyperbolicSecant(dur_s, dt, BW_Hz, beta)
            
            if (nargin < 4)
                beta = 6;
            end
            DeltaOmega = 2*pi*BW_Hz;
            mu = DeltaOmega*dur_s/(2*pi)/beta;
            if (mu+eps < 2)
                warning ('mu = %f should be higher than 2 for an adiabatic passage', mu);
            end
            
            p = pTXPulseObjBase('Nc', 1, 'N', round(dur_s/dt), 'dt', dt);
            T = p.getDur();
            t = cumsum(p.dt) - 0.5 * p.dt;
            t = t - mean(t);
            u = 2*t/T;
            a  = sech(beta * u);
            ph = (0.5*DeltaOmega*dt) * cumtrapz(tanh(-beta * u));
            ph = ph - median(ph);
            p.dRF(1,:) = dt * a .* exp(-1i*ph);
            p.rasterTime = min(dt);
            p.refAmplitude = DeltaOmega;
            p.refAmplitudeUnit = 'rad/s';
            
            
        end
        
        function [p, info] = importFromIniFileVB17(f, varargin)
            
            par.master = 0;
            %par.flip_rfphase = true;
            par.gyr = (pTXUtils.gyromagneticRatio('1h'));
            %par.gradAxes = eye(3);
            %par.trim = '';
            
            par = pTXUtils.readopt(par, varargin{:});
            
            if (~ischar(f))
                
                p(numel(f),1) = pTXPulseObjBase;
                info = cell(size(p));
                
                for i = 1:numel(p)
                    [p(i), info{i}] = pTXPulseObjBase.importFromIniFileVB17(fullfile(f(i).folder, f(i).name), par);
                end
                
                %try
                info = cat(1, info{:});
                %catch err
                %    err.getReport();
                %end
                
                return;
            end
            
            
            [info, RFShape, GradientShape, channelIndex] = pTXUtils.readIniFileVB17(f);
            assert(all(diff(channelIndex)==1) && channelIndex(1) == 0, 'bad channel index, please check ini file %s', f);
            
            if (isfield(info, 'SampleTime'))
                dt = sscanf(info.SampleTime, '%f'); % us
                dt = dt * 1e-6;
            elseif (isfield(info, 'TimeStep'))
                dt = sscanf(info.TimeStep, '%f'); % us
                dt = dt * 1e-6;
            else
                warning ('Cannot read TimeStep in header file. Assuming 10 us');
                dt = 1e-5;
            end
            
            p = pTXPulseObjBase('Nc', size(RFShape, 2), 'N', size(RFShape, 1), 'dt', dt, 'rastertime', dt);
            p.gyr = par.gyr;
            p.dRF(:,:) = dt * RFShape.';
            
            p.dRF(:,:) = circshift(p.dRF, par.master - find(channelIndex==0), 1);
            p.dk(:,:) = (dt * 2*pi*p.neggyr()*1e-3) * GradientShape.'; % gyromagnetic ratio is < 0, dk anf GradientSahep have opposite sign
            
            
            if (isfield(info, 'RotMat'))
                
                % define gradAxes
                
                R = sscanf(info.RotMat, '%f;%f;%f,%f;%f;%f,%f;%f;%f');
                assert(numel(R) == 9, 'bad object RotMat');
                R = reshape(R, 3, 3);
            else
                
                warning ('assuming RotMat = eye(3)');
                R = eye(3);
                
            end
            
            p.backProjectK(p.sRotMatUnity() * R);
            
            p.FA = sscanf(info.NominalFlipAngle, '%f');
            p.refAmplitude = sscanf(info.MaxAbsRF, '%f');
            p.refAmplitudeUnit = 'Volt';
                        
        end
        
        function index = supp(s)
            
            assert (ndims(s) <= 2 && (size(s,1) == 1 || size(s,2) == 1), 'expecting line or column');
            s = reshape(s, 1, numel(s));
            
            if (~islogical(s))
                s = abs(s)>0;
            end
            
            index = 0;
            
            if (~any(s))
                return;
            end
            
            index = 1+find(diff(s));
            
            if (s(1))
                
                index = [1, index];
                
            end
            
            
            if (s(end))
                index = [index, numel(s)+1];
            end
            
            
        end
        
        function p = Compose(p1,p2, varargin)
            
            p = p1.new();
            p.compose(p2, varargin{:});
            
            
        end
        
        function q = Inverse(p)
            q = p.new();
            q.inverse();
        end
        
        function q = PhaseShift(p, ph)
            q = p.new();
            q.phaseshift(ph);
        end
        
    end
    
end
