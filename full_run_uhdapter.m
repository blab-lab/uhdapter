function [ ] = full_run_uhdapter(outputdir,expt)


if nargin < 1 || isempty(outputdir), outputdir = 'C:\Users\Public\Documents\experiments\stress'; end
if nargin < 2, expt = []; end

run_uhdapter_audapter(outputdir,[],0)
