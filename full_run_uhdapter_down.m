function [ ] = full_run_uhdapter_down(outputdir,expt)


if nargin < 1 || isempty(outputdir), outputdir = 'C:\Users\Public\Documents\experiments\uhdapter_down'; end
if nargin < 2, expt = []; end

choice = input('You are running uhdapter DOWN. Was this your intent? (y/n): ', 's');

if strcmpi(choice,'y')

run_uhdapter_audapter_down(outputdir,[],0)
  
end