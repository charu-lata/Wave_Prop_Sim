function [hdr,seis] = filter_butter(hdr,seis_unfilt,fmini,fmax,dboct,tf_taper,flag,sampint)

%%%  Filters traces for high pass, low pass, or band pass using a
%%%%  butterworth filter
%
% fmini = low frequency  corner of filter in Hz
% fmax  = high frequency corner of filter in Hz
% dboct = roll off of filter in dB/octave
%
% If flag > 0 use a  causal filter  (good for active source data')
% If flag < 0 use an acausal filter (good for broadband data esp. surface waves)
%
% Give only fmini for high pass
% Give only fmax for low pass
% Give both for bandpass
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if isempty(hdr)==0
sampint = [hdr.sampint];
sampint=sampint(1); %Assuming all the same.
end

if (fmini>0 || fmax>0)

    %----Detrend the records
    seis = detrend(seis_unfilt);
       
    %----Window the records for filter application
    [m,n] = size(seis);
    
    if(tf_taper)
        %Note: for some reason, the window works best when it goes from -1
        %to 1, rather than 0 to 1 as is the usual method.
        %taper = 2*blackmanharris(m)-1;
        taper_fraction = .5;
        taper = tukeywin(m,taper_fraction);
        %taper = blackmanharris(m);
        taper = repmat(taper(:),[1 n]);
        seis = seis.*taper;
    end

    if fmini
        wmin=max(0,min(1,abs(fmini*sampint*2)));
    end
    if fmax
        wmax=max(0,min(1,abs(fmax*sampint*2)));
    end
    if fmini && ~fmax %High pass
        db=dboct;
        if dboct==0; db=24; end
        o=buttord(wmin,wmin/2,3,abs(db));
        if flag>0
            [b,a]=butter(o,wmin,'high');
            seis=filter(b,a,seis);
        else
            [b,a]=butter(ceil(o/2),wmin,'high');
            seis=filtfilt(b,a,seis);
        end
    elseif fmax && ~fmini %Low pass
        db=dboct;
        if dboct==0; db=24; end
        o=buttord(wmax,wmax/2,3,abs(db));
        if flag>0
            [b,a]=butter(o,wmax);
            seis=filter(b,a,seis);
        else
            [b,a]=butter(ceil(o/2),wmax);
            seis=filtfilt(b,a,seis);
        end
    elseif fmax && fmini %band pass
        db=dboct;
        if dboct==0; db=24; end

        o=buttord(wmin,wmin/2,3,abs(db));
        %Another way:
        %o=buttord([.8*wmin 1.2*wmax],[wmin wmax],3,abs(dboct))

        if flag>0
            [b,a]=butter(o,[wmin wmax]);
            seis=filter(b,a,seis);
        else
            [b,a]=butter(ceil(o/2),[wmin wmax]);
            seis=filtfilt(b,a,seis);
        end
    end
    if isempty(hdr)==0
    for i=1:size(hdr,2)
        hdr(i).hiCutFreq   = fmini;
        hdr(i).hiCutSlope  = dboct;
        hdr(i).lowCutFreq  = fmax;
        hdr(i).lowCutSlope = dboct;
    end
    end

else
    seis=seis_unfilt;
end
mmed = median(seis);
mmed = repmat(mmed,[m 1]);
seis = seis - mmed;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('    Butterworth Filter Applied\n\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%