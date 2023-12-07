% EMD.M 
%
% G. Rilling, July 2002
%
% computes EMD (Empirical Mode Decomposition) according to:
%
% N. E. Huang et al., "The empirical mode decomposition and the 
% Hilbert spectrum for non-linear and non stationary time series analysis,"  
% Proc. Royal Soc. London A, Vol. 454, pp. 903-995, 1998
%
% with variations reported in:
%
% G. Rilling, P. Flandrin and P. Gonçalvès
% "On Empirical Mode Decomposition and its algorithms"
% IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing
% NSIP-03, Grado (I), June 2003
%
% stopping criterion for sifting : 
%   at each point : mean amplitude < threshold2*envelope amplitude
%   &
%   mean of boolean array ((mean amplitude)/(envelope amplitude) > threshold) < tolerance
%   &
%   |#zeros-#extrema|<=1
%
% inputs:  - x : analysed signal (line vector)
%          - t (optional) : sampling times (line vector) (default : 1:length(x))
%          - stop (optional) : threshold, threshold2 and tolerance (optional)
%                              for sifting stopping criterion 
%                              default : [0.05,0.5,0.05]
%          - tst (optional) : if equals to 1 shows sifting steps with pause
%                             if equals to 2 no pause
%
% outputs: - imf : intrinsic mode functions (last line = residual)
%          - ort : index of orthogonality
%          - nbits : number of iterations for each mode
%
% calls:   - extr (finds extrema and zero-crossings)
%          - io : computes the index of orthogonality

function [imf,ort,nbits] = func_emd2(x,t,stop,tst)

% default for stopping
defstop = [0.05,0.5,0.05];

if(nargin==1)
  t = 1:length(x);
  stop = defstop;
  tst = 0;
end

if(nargin==2)
  stop = defstop;
  tst = 0;
end

if (nargin==3)
  tst=0;
end

S = size(x);
if ((S(1) > 1) & (S(2) > 1)) | (length(S) > 2)
  error('x must have only one row or one column')
end

if S(1) > 1
  x = x';
end

S = size(t);
if ((S(1) > 1) & (S(2) > 1)) | (length(S) > 2)
  error('t must have only one row or one column')
end

if S(1) > 1
  t = t';
end

if (length(t)~=length(x))
  error('x and t must have the same length')
end

S = size(stop);
if ((S(1) > 1) & (S(2) > 1)) | (S(1) > 3) | (S(2) > 3) | (length(S) > 2)
  error('stop must have only one row or one column of max three elements')
end

if S(1) > 1
  stop = stop';
  S = size(stop);
end

if S(2) < 3
  stop(3)=defstop(3);
end

if S(2) < 2
  stop(2)=defstop(2);
end

sd = stop(1);
sd2 = stop(2);
tol = stop(3);

if tst
  figure
end

% maximum number of iterations
MAXITERATIONS=2000;

% maximum number of symmetrized points for interpolations
NBSYM = 2;

lx = length(x);

sdt(lx) = 0;
sdt = sdt+sd;
sd2t(lx) = 0;
sd2t = sd2t+sd2;

% maximum number of extrema and zero-crossings in residual
ner = lx;
nzr = lx;

r = x;
imf = [];
k = 1;

% iterations counter for extraction of 1 mode
nbit=0;

% total iterations counter
NbIt=0;

while ner > 2
        
  % current mode
  m = r;
  
  % mode at previous iteration
  mp = m;
  
  sx = sd+1;
  
  % tests if enough extrema to proceed
  test = 0;
  
  [indmin,indmax,indzer] = func_extr(m);
  lm=length(indmin);
  lM=length(indmax);
  nem=lm + lM;
  nzm=length(indzer);
  
  j=1;
  
  % sifting loop
  while ( mean(sx > sd) > tol | any(sx > sd2) | (abs(nzm-nem)>1)) & (test == 0) & nbit<MAXITERATIONS
    
    if(nbit>MAXITERATIONS/5 & mod(nbit,floor(MAXITERATIONS/10))==0)
      disp(['mode ',int2str(k),' nombre d iterations : ',int2str(nbit)])
      disp(['stop parameter mean value : ',num2str(s)])
    end
   
    % boundary conditions for interpolations :
        
     if indmax(1) < indmin(1)
      if m(1) > m(indmin(1))
        lmax = fliplr(indmax(2:min(end,NBSYM+1)));
        lmin = fliplr(indmin(1:min(end,NBSYM)));
        lsym = indmax(1);
      else
        lmax = fliplr(indmax(1:min(end,NBSYM)));
        lmin = [fliplr(indmin(1:min(end,NBSYM-1))),1];
        lsym = 1;
      end
    else

      if m(1) < m(indmax(1))
        lmax = fliplr(indmax(1:min(end,NBSYM)));
        lmin = fliplr(indmin(2:min(end,NBSYM+1)));
        lsym = indmin(1);
      else
        lmax = [fliplr(indmax(1:min(end,NBSYM-1))),1];
        lmin = fliplr(indmin(1:min(end,NBSYM)));
        lsym = 1;
      end
    end
    
    if indmax(end) < indmin(end)
      if m(end) < m(indmax(end))
        rmax = fliplr(indmax(max(end-NBSYM+1,1):end));
        rmin = fliplr(indmin(max(end-NBSYM,1):end-1));
        rsym = indmin(end);
      else
        rmax = [lx,fliplr(indmax(max(end-NBSYM+2,1):end))];
        rmin = fliplr(indmin(max(end-NBSYM+1,1):end));
        rsym = lx;
      end
    else
      if m(end) > m(indmin(end))
        rmax = fliplr(indmax(max(end-NBSYM,1):end-1));
        rmin = fliplr(indmin(max(end-NBSYM+1,1):end));
        rsym = indmax(end);
      else
        rmax = fliplr(indmax(max(end-NBSYM+1,1):end));
        rmin = [lx,fliplr(indmin(max(end-NBSYM+2,1):end))];
        rsym = lx;
      end
    end
    
    tlmin = 2*t(lsym)-t(lmin);
    tlmax = 2*t(lsym)-t(lmax);
    trmin = 2*t(rsym)-t(rmin);
    trmax = 2*t(rsym)-t(rmax);
    
    % in case symmetrized parts do not extend enough
    if tlmin(1) > t(1) | tlmax(1) > t(1)
      if lsym == indmax(1)
        lmax = fliplr(indmax(1:min(end,NBSYM)));
      else
        lmin = fliplr(indmin(1:min(end,NBSYM)));
      end
      if lsym == 1
        error('bug')
      end
      lsym = 1;
      tlmin = 2*t(lsym)-t(lmin);
      tlmax = 2*t(lsym)-t(lmax);
    end   
    
    if trmin(end) < t(lx) | trmax(end) < t(lx)
      if rsym == indmax(end)
        rmax = fliplr(indmax(max(end-NBSYM+1,1):end));
      else
        rmin = fliplr(indmin(max(end-NBSYM+1,1):end));
      end
      if rsym == lx
        error('bug')
      end
      rsym = lx;
      trmin = 2*t(rsym)-t(rmin);
      trmax = 2*t(rsym)-t(rmax);
    end 
          
    mlmax =m(lmax); 
    mlmin =m(lmin);
    mrmax =m(rmax); 
    mrmin =m(rmin);
     
    % definition of envelopes from interpolation
        
    envmax = interp1([tlmax t(indmax) trmax],[mlmax m(indmax) mrmax],t,'spline');
    envmin = interp1([tlmin t(indmin) trmin],[mlmin m(indmin) mrmin],t,'spline');

    envmoy = (envmax + envmin)/2;

    m = m - envmoy;
   
    [indmin,indmax,indzer] = func_extr(m);
    lm=length(indmin);
    lM=length(indmax);
    nem = lm + lM;
    nzm = length(indzer);

    % evaluation of mean zero
    sx=2*(abs(envmoy))./(abs(envmax-envmin));
    s = mean(sx);
    
    % display
        
     if tst
      subplot(4,1,1)
      plot(t,mp);hold on;
      plot(t,envmax,'--k');plot(t,envmin,'--k');plot(t,envmoy,'r');

      title(['IMF ',int2str(k),';   iteration ',int2str(nbit),' before sifting']);
      set(gca,'XTick',[])
      hold  off

      subplot(4,1,2)
      plot(t,sx)
      hold on
      plot(t,sdt,'--r')
      plot(t,sd2t,':k')
      title('stop parameter')
      set(gca,'XTick',[])
      hold off

      subplot(4,1,3)
      plot(t,m)
      title(['IMF ',int2str(k),';   iteration ',int2str(nbit),' after sifting']);
      set(gca,'XTick',[])

      subplot(4,1,4);
      plot(t,r-m)
      title('residue');
      disp(['stop parameter mean value : ',num2str(s)])
      if tst == 2
        pause(0.01)
      else
        pause
      end
      
    end
 
    % end loop : stops if not enough extrema
    if nem < 3
      test = 1;
    end

    mp = m;
    nbit=nbit+1;
    NbIt=NbIt+1;

    if(nbit==(MAXITERATIONS-1))
      warning(['forced stop of sifting : too many iterations... mode ',int2str(k),'. stop parameter mean value : ',num2str(s)])
    end
  
  end
  imf(k,:) = m;
  if tst
    disp(['mode ',int2str(k),' enregistre'])
  end
  nbits(k) = nbit;
  k = k+1;
  r = r - m;
  [indmin,indmax,indzer] = func_extr(r);
  ner = length(indmin) + length(indmax);
  nzr = length(indzer);
  nbit=1;

  if (max(r) - min(r)) < (1e-10)*(max(x) - min(x))
    if ner > 2
      warning('forced stop of EMD : too small amplitude')
    else

      disp('forced stop of EMD : too small amplitude')
    end
    break
  end
  
end

imf(k,:) = r;

ort = func_io(x,imf);

if tst
  close
end
