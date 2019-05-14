
%Pavithra
%April 2019
%dot_bracket_modify.m

%%%%%%
%Saves bpp and reactivity matrix to external .txt files named 'bpp.txt' and 'reaccs_only.txt'
%If prefer different names, modify function below.

%********Requires bpp matrix and reactivity matrix from HiTRACE pipeline, or any pipeline that delivers similar boostrapping values and chemical reactivity data in matrix format!******



%Example Run:
%dot_bracket_modify(bpp_1M7, d_1M7)
%%%%%%


function dot_bracket_modify(bpp, reactivity)
    
    fileID = fopen('reaccs_only.txt','w');
    for i = [1:length(reactivity)]
        if reactivity(i) >0
            fprintf(fileID,'+%.6f ',reactivity(i));
        end
        if reactivity(i)<0
            fprintf(fileID,'%.6f ',reactivity(i));
        end
    end
    
    fileID = fopen('bpp.txt','w');
    for i = [1:length(bpp)]
        for j = [1:length(bpp)]
            fprintf(fileID,'%.2f ',bpp(i,j));
        end
    fprintf(fileID,'\n');
    end
