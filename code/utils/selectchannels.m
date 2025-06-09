function [chidx, enumoffset] = selectchannels(monkey,datapath,filename)

shortfolder = 'shortedchan/';

switch monkey
    case 'Sa'
        PFCarray = 'PFCLR';
        switch PFCarray
            case 'PFCL'
                chans       = 1:96;   % recording array
                shortarray  = 'shortchmaskrec';
            case 'PFCR'
                chans       = 97:192; % uStim array
                shortarray  = 'shortchmaskstim';
            case 'PFCLR'
                chans1      = 1:96;   % recording array
                chans2      = 97:192; % uStim array
                shortarray1 = 'shortchmaskrec';
                shortarray2 = 'shortchmaskstim';
        end  
        % Shorted channels
        shortfile        = [filename '_shorted.mat'];
        if strcmp(PFCarray,'PFCLR')
            shortchmask1 = load([datapath shortfolder shortfile],shortarray1);
            chidx1       = chans1(~shortchmask1.(shortarray1){1});
            shortchmask2 = load([datapath shortfolder shortfile],shortarray2);
            chidx2       = chans2(~shortchmask2.(shortarray2){1});
            chidx        = [chidx1 chidx2];
        else
            shortchmask = load([datapath shortfolder shortfile],shortarray);
            chidx = chans(~shortchmask.(shortarray){1});
        end
        % adjust numbering for electrodes on different banks
        enumoffset = 256;
        
    case 'Wa'
        PFCarray = 'PFCL';
        switch PFCarray
            case 'PFCL'
                chans = 1:96; % uStim array
        end      
        % shorted channels
        shortfile   = [filename '_shorted.mat'];
        load([datapath shortfolder shortfile],'channelsKeep');
        chidx       = channelsKeep;
        % adjust numbering for electrodes on different banks
        enumoffset = 0;
end

end