function B=correlation2d(A,B,mode);
% fast correlation
% A:pxq
% B:p2xq2xk
% mode: 'full' | 'symmetric' | 'valid'
% output = correlation bw A and B(:,:,i) for all i
%
% switch according to nnz(A) that can revert back to conv2:
% conv2 faster when nnz(A)<numel(A)*0.1 (cf A=imageEdges)
% correlation2d faster when nnz(A)>numel(A)*0.1
% when A is pxqxKA and B is p2xq2, a for loop on KA is executed
% Timothee Cour, 29-Aug-2006 07:49:15


[pA,qA,KA] = size(A);
[pB,qB,K]=size(B);
if KA>1
    if K>1
        error('?');
    else
        for i=KA:-1:1
            C(:,:,i)=correlation2d(A(:,:,i),B,mode);            
        end
        B=C;
        return;
    end
end

thresSwitch=0.1;
if nnz(A)<numel(A)*thresSwitch
    K=size(B,3);
    A=double(A);
    for i=K:-1:1
        B2(:,:,i)=conv2(A,fliplr(flipud(B(:,:,i))),mode);
    end
    B=B2;
    return;
end

%TODO: check when A = all zero (shouldn't occur with thresSwitch)
%mode='full'|'same'|'valid'
if nargin <3
    mode='full';
end

B=B(pB:-1:1,qB:-1:1,:);

p2=pA+pB-1;
q2=qA+qB-1;
p2 = compute_good_fft_dim(p2);
q2 = compute_good_fft_dim(q2);
% disp(['factor(p2) = ' num2str(factor(p2)) , ' ; factor(q2) = ' num2str(factor(q2)) ]);

% fftw('planner','patient');
% fftw('planner','measure');%VOIR

A=fft2(A,p2,q2);
B=fft2(B,p2,q2);

% A=real(A);
% B=real(B);

% C=repmat(A,[1,1,K]).*B;
% C2 = mex_matrix_op_repmat_vector(B,A,'.*','col');
mex_matrix_op_repmat_vector(B,A,'.*','col');
clear A;
B=ifft2(B,'symmetric');

pC=pA+pB-1;
qC=qA+qB-1;

switch mode
    case 'full'
        B = B(1:pC,1:qC,:);
    case 'same'
        B=B(ceil((pB+1)/2):ceil((pB+1)/2)+pA-1,ceil((qB+1)/2):ceil((qB+1)/2)+qA-1,:);
%         B=B(ceil(pB/2):ceil(pB/2)+pA-1,ceil(qB/2):ceil(qB/2)+qA-1,:);
    case 'valid'
        B=B(pB:pC-pB+1,qB:qC-qB+1,:);
end
