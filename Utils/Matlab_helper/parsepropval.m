function prop = parsepropval(prop,varargin)
%parsepropval: Parse property/value pairs and return a structure.
%  Manages property/value pairs like MathWorks Handle Graphics functions.
%  This means that in addition to passing in Property name strings and
%  Values, you can also include structures with appropriately named fields.
%  The full, formal property names are defined by the defaults structure
%  (first input argument). This is followed by any number of property/value
%  pairs and/or structures with property names for the field names.  The
%  property name matching is case-insensitive and needs only be
%  unambiguous.
%
%  For example,
%
%    params.FileName = 'Untitled';
%    params.FileType = 'text';
%    params.Data = [1 2 3];
%    s.dat = [4 5 6];
%    parsepropval(params,'filenam','mydata.txt',s,'filety','binary')
%
%  returns a structure with the same field names as params, filled in
%  according to the property/value pairs and structures passed in.
%    ans = 
%        FileName: 'mydata.txt'
%        FileType: 'binary'
%        Data: [4 5 6]
%
%  The inputs are processed from left to right so if any property is
%  specified more than once the latest value is retained.
%
%  An error is generated if property names are ambiguous.  Values can be
%  any MATLAB variable.
%
%  Typical use is in a function with a variable number of input arguments.
%  For example,
%
%  function myfun(varargin)
%    properties.prop1 = [];
%    properties.prop2 = 'default';
%    properties = parsepropval(properties,varargin{:});
% 
%  New functionalities: 
%   'create': Create property if not exists (usefull for structs) [default]
%   'ignore': Ignore property if not exists (usefull for class objects)

% Version: 1.0, 13 January 2009
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})
%
%  File: parsepropval.m
%  Directory: 7_ftools/ftools/v12/utilities
%  Second author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. June 08. (2019b)
%
% 
% 'create':      strong match for fieldnames, create if not exist
% 'createloose': loose  match for fieldnames, create if not exist
% 'ignore':      strong match for fieldnames, ignore if not exist
% 'ignoreloose': loose  match for fieldnames, ignore if not exist

% POLCZ:
if ischar(prop)
    fun = prop;
    prop = varargin{1};
    varargin = varargin(2:end);
else
    fun = 'createloose';
end    

% Process inputs and set prop fields.
properties = fieldnames(prop);
arg_index = 1;
while arg_index <= length(varargin)
	arg = varargin{arg_index};
	if ischar(arg)
        [good,prop_index,properties] = match_property(arg,properties,fun);
		if good, prop.(properties{prop_index}) = varargin{arg_index + 1}; end
		arg_index = arg_index + 2;
	elseif isstruct(arg)
		arg_fn = fieldnames(arg);
		for i = 1:length(arg_fn)
			[good,prop_index,properties] = match_property(arg_fn{i},properties,fun);
			if good, prop.(properties{prop_index}) = arg.(arg_fn{i}); end
		end
		arg_index = arg_index + 1;
	else
		error(['Properties must be specified by property/value pairs',...
			' or structures.'])
	end
end


function [good,prop_index,properties] = match_property(arg,properties,fun)
good = 1;
if strcmp(fun,'create') || strcmp(fun,'ignore')
    prop_index = find(strcmp(arg,properties));
elseif strcmp(fun,'ignoreloose') || strcmp(fun,'createloose')
    prop_index = find(strcmpi(arg,properties));
    if isempty(prop_index)
        prop_index = find(strncmpi(arg,properties,length(arg)));
    end
end
if length(prop_index) ~= 1
    if strncmp(fun,'ignore',6) % ignore | ignoreloose
        good = 0;
    elseif strncmp(fun,'create',6) % create | createloose
        properties = [properties ; {arg}];
        prop_index = numel(properties);
    end
end
