clear all;
fpath_b = 'D:\\WinSCP\\PDE\\Matrix\\biharmonic.txt';
fpath_l = 'D:\\WinSCP\\PDE\\Matrix\\laplacian.txt';
fpath_s = 'D:\\WinSCP\\PDE\\Matrix\\setting.txt';
fpath_ls = 'D:\\WinSCP\\PDE\\Matrix\\ls.txt';

% read biharmonic matrix
fid = fopen(fpath_b,'r');
[A,Num] = fscanf(fid,'%f');
N = Num/3;
Row = zeros(N,1);
Col = zeros(N,1);
Element = zeros(N,1);
for i = 1:1:Num
    if mod(i,3) == 1
        Row(fix((i+2)/3)) = A(i);
    elseif mod(i,3) == 2
        Col(fix((i+1)/3)) = A(i);
    else
        Element(fix(i/3)) = A(i);
    end
end
clear A;
fclose(fid);
Row(Row==0) = [];
Col(Col==0) = [];
Element = Element(1:length(Row));
MatA = sparse(Row,Col,Element);

%read laplacian matrix
fid = fopen(fpath_l,'r');
[A,Num] = fscanf(fid,'%f');
N = Num/3;
Row = zeros(N,1);
Col = zeros(N,1);
Element = zeros(N,1);
for i = 1:1:Num
    if mod(i,3) == 1
        Row(fix((i+2)/3)) = A(i);
    elseif mod(i,3) == 2
        Col(fix((i+1)/3)) = A(i);
    else
        Element(fix(i/3)) = A(i);
    end
end
clear A;
fclose(fid);
Row(Row==0) = [];
Col(Col==0) = [];
Element = Element(1:length(Row));
MatB = sparse(Row,Col,Element);
clear Row;
clear Col;
clear Element;

fid = fopen(fpath_s,'r');
[A,N] = fscanf(fid,'%f');
L = A(1);
H = A(2);
Nx = A(3);
Nz = A(4);
fclose(fid);

fid = fopen(fpath_ls,'r');
[A,N] = fscanf(fid,'%f');
ls1 = zeros(Nx,1);
ls2 = zeros(Nx,1);
for i = 1:1:N
    if mod(i,2)==1
        ls1(fix((i+1)/2)) = A(i);
    else
        ls2(fix(i/2)) = A(i);
    end
end
fclose(fid);
        
        