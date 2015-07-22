function slash = GetSlash

if isunix || ismac
    slash = '/';
    
elseif ispc
    slash = '\';
        
else 
    error('unrecognised/unsupported operating system');
end

