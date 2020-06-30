function make
    % make.m
    % runs the Makefile to build the DAS3 mex function
    % should be configured depending on each computer's location of the make tool

    computer = getenv('COMPUTERNAME');

    if strcmp(computer, 'LRI-102855')       % Ton's computer
        !"C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.25.28610\bin\Hostx64\x64\nmake" /nologo
    elseif strcmp(computer, '...')

    else
       error('In make.m, please add MAKE setup for your computer: %s', computer); 
    end
end 