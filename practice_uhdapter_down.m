function [ ] = practice_uhdapter_down(outputdir,expt)


if nargin < 1 || isempty(outputdir), outputdir = 'C:\Users\Public\Documents\experiments\uhdapter_down'; end
if nargin < 2, expt = []; end

run_uhdapter_audapter_down(outputdir,[],1)
