function [atom_str] = getAtomStr(atom)
    
    if strcmp(atom, 'Na')
        atom_str = '^{23}Na';
    elseif strcmp(atom, 'Rb')
        atom_str = '^{87}Rb';
    else
        atom_str = atom;
    end
    
end