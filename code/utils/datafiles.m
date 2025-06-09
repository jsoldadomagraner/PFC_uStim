function [filenames,numuStim,numang,stimorder] = datafiles(monkey)

switch monkey
    case 'Wa'
        filenames = {
            'Wa220801_s549'
            'Wa220802_s550'
            'Wa220803_s551'
            'Wa220804_s552'
            'Wa220805_s553'
            'Wa220808_s554'
            'Wa220809_s555'
            'Wa220810_s556'
            'Wa220811_s557'
            'Wa220812_s558'};
        numang    = 4;
        numuStim  = 4;
        stimorder = 1:4;
    case 'Sa'
        filenames = {'Sa210311_s224'};
        numang    = 4;
        numuStim  = 5;
        stimorder = [5,1,4,2,3]; %control is not first but last
end

end