function savefig(dir, name)
% This function allows you to quickly save all currently open figures with
% a custom filename for each in multiple formats.  To use the function
% simply call savefigs with no arguments, then follow the prompts
%
% Upon execution this function will one-by-one bring each currently open
% figure to the foreground.  Then it will supply a text prompt in the main
% console window asking you for a filename.  It will save that figure to
% that filename in the .fig, .emf, .png, and .eps formats.  
%
% The formats that it saves in can be changed by commenting out or adding
% lines below.
%
% Copyright 2010 Matthew Guidry 
% matt.guidry ATT gmail DOTT com  (Email reformatted for anti-spam)

hfigs = get(0, 'children');                          %Get list of figures

names = {'timevssol', 'runtime', 'acctopo', 'accloss', 'accdup', 'au_test', 'allrfloglkl', 'rfsh', 'rfloglik', 'rfbycost', 'rfall', 'rf', 'reconall', 'recon'};
for m = 1:length(hfigs)
    figure(hfigs(m) ) %Bring Figure to foreground
    filename=fullfile(dir, strcat(names{m}, name));
    intext = strcat('Enter filename: (Current = ', filename, ' ) 1 to validate and 0 to skip)\n');
    whos intext
    an = input(intext, 's');%Prompt user
    if strcmp(an, '0')                        %Skip figure when user types 0
        continue
    elseif strcmp(an, '1')   
        eval(['print -depsc2 ' filename])   %Enhanced Postscript (Level 2 color) (Best for LaTeX documents)
    else
        eval(['print -depsc2 ' an])   %Enhanced Postscript (Level 2 color) (Best for LaTeX documents)      
    end
    
end