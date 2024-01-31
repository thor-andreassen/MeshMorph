function rotated_TransMat=rotateTransMat(input_TransMat,desired_rotation)
% This function takes in a single Transformation matrix, determines the
% orientation of the values, and con convert it to a different ordered
% Transformation matrix with the following patterns chosen by changing the
% number for desired_rotation. NOTE 1-4 are with the translation vector in
% a column format, and 5-8 are with the translation vector in a row format.
% If a non valid number is chosen, or ther is no number chosen, it will go
% to format 2 as the standard.

% 1. [0,1;R,T]
% 2. [1,0;T,R];
% 3. [R,T;0,1];
% 4. [T,R;1,0];
% 5. [1,T;0,R];
% 6. [0,R;1,T];
% 7. [R,0;T,1];
% 8. [T,1;R,0];
if nargin<1
    error('Not Enough Inputs into rotateTransMat');
elseif nargin==1
    desired_rotation=2;
elseif nargin>2
    error('Too Many Inputs into rotateTransMat');
end

if input_TransMat(1,1)==0 && input_TransMat(1,2)==0 && input_TransMat(1,3)==0
    rotmat=input_TransMat(2:4,1:3);
    transvec_col=input_TransMat(2:4,4);
elseif input_TransMat(1,2)==0 && input_TransMat(1,3)==0 && input_TransMat(1,4)==0
    rotmat=input_TransMat(2:4,2:4);
    transvec_col=input_TransMat(2:4,1);
elseif input_TransMat(4,1)==0 && input_TransMat(4,2)==0 && input_TransMat(4,3)==0
    rotmat=input_TransMat(1:3,1:3);
    transvec_col=input_TransMat(1:3,4);
elseif input_TransMat(4,2)==0 && input_TransMat(4,3)==0 && input_TransMat(4,4)==0
    rotmat=input_TransMat(1:3,2:4);
    transvec_col=input_TransMat(1:3,1);
elseif input_TransMat(2,1)==0 && input_TransMat(3,1)==0 && input_TransMat(4,1)==0
    rotmat=input_TransMat(2:4,2:4)';
    transvec_col=input_TransMat(1,2:4)';
elseif input_TransMat(1,1)==0 && input_TransMat(2,1)==0 && input_TransMat(3,1)==0
    rotmat=input_TransMat(1:3,2:4)';
    transvec_col=input_TransMat(4,2:4)';
elseif input_TransMat(1,4)==0 && input_TransMat(2,4)==0 && input_TransMat(3,4)==0
    rotmat=input_TransMat(1:3,1:3)';
    transvec_col=input_TransMat(4,1:3)';
elseif input_TransMat(2,4)==0 && input_TransMat(3,4)==0 && input_TransMat(4,4)==0
    rotmat=input_TransMat(2:4,1:3)';
    transvec_col=input_TransMat(1,1:3)';
else
    error('The input Transformation Matrix is not a valid Transformation Matrix');
end

rotated_TransMat=zeros(4);
if desired_rotation==1
    rotated_TransMat(2:4,1:3)=rotmat;
    rotated_TransMat(2:4,4)=transvec_col;
    rotated_TransMat(1,4)=1;
elseif desired_rotation==2
    rotated_TransMat(2:4,2:4)=rotmat;
    rotated_TransMat(2:4,1)=transvec_col;
    rotated_TransMat(1,1)=1;
elseif desired_rotation==3
    rotated_TransMat(1:3,1:3)=rotmat;
    rotated_TransMat(1:3,4)=transvec_col;
    rotated_TransMat(4,4)=1;
elseif desired_rotation==4
    rotated_TransMat(1:3,2:4)=rotmat;
    rotated_TransMat(1:3,1)=transvec_col;
    rotated_TransMat(4,1)=1;
elseif desired_rotation==5
    rotated_TransMat(2:4,2:4)=rotmat';
    rotated_TransMat(1,2:4)=transvec_col';
    rotated_TransMat(1,1)=1;
elseif desired_rotation==6
    rotated_TransMat(1:3,2:4)=rotmat';
    rotated_TransMat(4,2:4)=transvec_col';
    rotated_TransMat(4,1)=1;
elseif desired_rotation==7
    rotated_TransMat(1:3,1:3)=rotmat';
    rotated_TransMat(4,1:3)=transvec_col';
    rotated_TransMat(4,4)=1;
elseif desired_rotation==8
    rotated_TransMat(2:4,1:3)=rotmat';
    rotated_TransMat(1,1:3)=transvec_col';
    rotated_TransMat(1,4)=1;
else
    rotated_TransMat(2:4,2:4)=rotmat;
    rotated_TransMat(2:4,1)=transvec_col;
    rotated_TransMat(1,1)=1;
end