% function [Matching,Cost] = Hungarian(Perf)
% % 
% % [MATCHING,COST] = Hungarian_New(WEIGHTS)
% %
% % A function for finding a minimum edge weight matching given a MxN Edge
% % weight matrix WEIGHTS using the Hungarian Algorithm.
% %
% % An edge weight of Inf indicates that the pair of vertices given by its
% % position have no adjacent edge.
% %
% % MATCHING return a MxN matrix with ones in the place of the matchings and
% % zeros elsewhere.
% % 
% % COST returns the cost of the minimum matching
% 
% % Written by: Alex Melin 30 June 2006
% 
% 
%  % Initialize Variables
%  Matching = zeros(size(Perf));
% 
% % Condense the Performance Matrix by removing any unconnected vertices to
% % increase the speed of the algorithm
% 
%   % Find the number in each column that are connected
%     num_y = sum(~isinf(Perf),1);
%   % Find the number in each row that are connected
%     num_x = sum(~isinf(Perf),2);
%     
%   % Find the columns(vertices) and rows(vertices) that are isolated
%     x_con = find(num_x~=0);
%     y_con = find(num_y~=0);
%     
%   % Assemble Condensed Performance Matrix
%     P_size = max(length(x_con),length(y_con));
%     P_cond = zeros(P_size);
%     P_cond(1:length(x_con),1:length(y_con)) = Perf(x_con,y_con);
%     if isempty(P_cond)
%       Cost = 0;
%       return
%     end
% 
%     % Ensure that a perfect matching exists
%       % Calculate a form of the Edge Matrix
%       Edge = P_cond;
%       Edge(P_cond~=Inf) = 0;
%       % Find the deficiency(CNUM) in the Edge Matrix
%       cnum = min_line_cover(Edge);
%     
%       % Project additional vertices and edges so that a perfect matching
%       % exists
%       Pmax = max(max(P_cond(P_cond~=Inf)));
%       P_size = length(P_cond)+cnum;
%       P_cond = ones(P_size)*Pmax;
%       P_cond(1:length(x_con),1:length(y_con)) = Perf(x_con,y_con);
%    
% %*************************************************
% % MAIN PROGRAM: CONTROLS WHICH STEP IS EXECUTED
% %*************************************************
%   exit_flag = 1;
%   stepnum = 1;
%   while exit_flag
%     switch stepnum
%       case 1
%         [P_cond,stepnum] = step1(P_cond);
%       case 2
%         [r_cov,c_cov,M,stepnum] = step2(P_cond);
%       case 3
%         [c_cov,stepnum] = step3(M,P_size);
%       case 4
%         [M,r_cov,c_cov,Z_r,Z_c,stepnum] = step4(P_cond,r_cov,c_cov,M);
%       case 5
%         [M,r_cov,c_cov,stepnum] = step5(M,Z_r,Z_c,r_cov,c_cov);
%       case 6
%         [P_cond,stepnum] = step6(P_cond,r_cov,c_cov);
%       case 7
%         exit_flag = 0;
%     end
%   end
% 
% % Remove all the virtual satellites and targets and uncondense the
% % Matching to the size of the original performance matrix.
% Matching(x_con,y_con) = M(1:length(x_con),1:length(y_con));
% Cost = sum(sum(Perf(Matching==1)));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   STEP 1: Find the smallest number of zeros in each row
% %           and subtract that minimum from its row
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [P_cond,stepnum] = step1(P_cond)
% 
%   P_size = length(P_cond);
%   
%   % Loop throught each row
%   for ii = 1:P_size
%     rmin = min(P_cond(ii,:));
%     P_cond(ii,:) = P_cond(ii,:)-rmin;
%   end
% 
%   stepnum = 2;
%   
% %**************************************************************************  
% %   STEP 2: Find a zero in P_cond. If there are no starred zeros in its
% %           column or row start the zero. Repeat for each zero
% %**************************************************************************
% 
% function [r_cov,c_cov,M,stepnum] = step2(P_cond)
% 
% % Define variables
%   P_size = length(P_cond);
%   r_cov = zeros(P_size,1);  % A vector that shows if a row is covered
%   c_cov = zeros(P_size,1);  % A vector that shows if a column is covered
%   M = zeros(P_size);        % A mask that shows if a position is starred or primed
%   
%   for ii = 1:P_size
%     for jj = 1:P_size
%       if P_cond(ii,jj) == 0 && r_cov(ii) == 0 && c_cov(jj) == 0
%         M(ii,jj) = 1;
%         r_cov(ii) = 1;
%         c_cov(jj) = 1;
%       end
%     end
%   end
%   
% % Re-initialize the cover vectors
%   r_cov = zeros(P_size,1);  % A vector that shows if a row is covered
%   c_cov = zeros(P_size,1);  % A vector that shows if a column is covered
%   stepnum = 3;
%   
% %**************************************************************************
% %   STEP 3: Cover each column with a starred zero. If all the columns are
% %           covered then the matching is maximum
% %**************************************************************************
% 
% function [c_cov,stepnum] = step3(M,P_size)
% 
%   c_cov = sum(M,1);
%   if sum(c_cov) == P_size
%     stepnum = 7;
%   else
%     stepnum = 4;
%   end
%   
% %**************************************************************************
% %   STEP 4: Find a noncovered zero and prime it.  If there is no starred
% %           zero in the row containing this primed zero, Go to Step 5.  
% %           Otherwise, cover this row and uncover the column containing 
% %           the starred zero. Continue in this manner until there are no 
% %           uncovered zeros left. Save the smallest uncovered value and 
% %           Go to Step 6.
% %**************************************************************************
% function [M,r_cov,c_cov,Z_r,Z_c,stepnum] = step4(P_cond,r_cov,c_cov,M)
% 
% P_size = length(P_cond);
% 
% zflag = 1;
% while zflag  
%     % Find the first uncovered zero
%       row = 0; col = 0; exit_flag = 1;
%       ii = 1; jj = 1;
%       while exit_flag
%           if P_cond(ii,jj) == 0 && r_cov(ii) == 0 && c_cov(jj) == 0
%             row = ii;
%             col = jj;
%             exit_flag = 0;
%           end      
%           jj = jj + 1;      
%           if jj > P_size; jj = 1; ii = ii+1; end      
%           if ii > P_size; exit_flag = 0; end      
%       end
% 
%     % If there are no uncovered zeros go to step 6
%       if row == 0
%         stepnum = 6;
%         zflag = 0;
%         Z_r = 0;
%         Z_c = 0;
%       else
%         % Prime the uncovered zero
%         M(row,col) = 2;
%         % If there is a starred zero in that row
%         % Cover the row and uncover the column containing the zero
%           if sum(find(M(row,:)==1)) ~= 0
%             r_cov(row) = 1;
%             zcol = find(M(row,:)==1);
%             c_cov(zcol) = 0;
%           else
%             stepnum = 5;
%             zflag = 0;
%             Z_r = row;
%             Z_c = col;
%           end            
%       end
% end
%   
% %**************************************************************************
% % STEP 5: Construct a series of alternating primed and starred zeros as
% %         follows.  Let Z0 represent the uncovered primed zero found in Step 4.
% %         Let Z1 denote the starred zero in the column of Z0 (if any). 
% %         Let Z2 denote the primed zero in the row of Z1 (there will always
% %         be one).  Continue until the series terminates at a primed zero
% %         that has no starred zero in its column.  Unstar each starred 
% %         zero of the series, star each primed zero of the series, erase 
% %         all primes and uncover every line in the matrix.  Return to Step 3.
% %**************************************************************************
% 
% function [M,r_cov,c_cov,stepnum] = step5(M,Z_r,Z_c,r_cov,c_cov)
% 
%   zflag = 1;
%   ii = 1;
%   while zflag 
%     % Find the index number of the starred zero in the column
%     rindex = find(M(:,Z_c(ii))==1);
%     if rindex > 0
%       % Save the starred zero
%       ii = ii+1;
%       % Save the row of the starred zero
%       Z_r(ii,1) = rindex;
%       % The column of the starred zero is the same as the column of the 
%       % primed zero
%       Z_c(ii,1) = Z_c(ii-1);
%     else
%       zflag = 0;
%     end
%     
%     % Continue if there is a starred zero in the column of the primed zero
%     if zflag == 1;
%       % Find the column of the primed zero in the last starred zeros row
%       cindex = find(M(Z_r(ii),:)==2);
%       ii = ii+1;
%       Z_r(ii,1) = Z_r(ii-1);
%       Z_c(ii,1) = cindex;    
%     end    
%   end
%   
%   % UNSTAR all the starred zeros in the path and STAR all primed zeros
%   for ii = 1:length(Z_r)
%     if M(Z_r(ii),Z_c(ii)) == 1
%       M(Z_r(ii),Z_c(ii)) = 0;
%     else
%       M(Z_r(ii),Z_c(ii)) = 1;
%     end
%   end
%   
%   % Clear the covers
%   r_cov = r_cov.*0;
%   c_cov = c_cov.*0;
%   
%   % Remove all the primes
%   M(M==2) = 0;
% 
% stepnum = 3;
% 
% % *************************************************************************
% % STEP 6: Add the minimum uncovered value to every element of each covered
% %         row, and subtract it from every element of each uncovered column.  
% %         Return to Step 4 without altering any stars, primes, or covered lines.
% %**************************************************************************
% 
% function [P_cond,stepnum] = step6(P_cond,r_cov,c_cov)
% a = find(r_cov == 0);
% b = find(c_cov == 0);
% minval = min(min(P_cond(a,b)));
% 
% P_cond(find(r_cov == 1),:) = P_cond(find(r_cov == 1),:) + minval;
% P_cond(:,find(c_cov == 0)) = P_cond(:,find(c_cov == 0)) - minval;
% 
% stepnum = 4;
% 
% function cnum = min_line_cover(Edge)
% 
%   % Step 2
%     [r_cov,c_cov,M,stepnum] = step2(Edge);
%   % Step 3
%     [c_cov,stepnum] = step3(M,length(Edge));
%   % Step 4
%     [M,r_cov,c_cov,Z_r,Z_c,stepnum] = step4(Edge,r_cov,c_cov,M);
%   % Calculate the deficiency
%     cnum = length(Edge)-sum(r_cov)-sum(c_cov);

function [C,T]=hungarian(A)
%HUNGARIAN Solve the Assignment problem using the Hungarian method.
%
%[C,T]=hungarian(A)
%A - a square cost matrix.
%C - the optimal assignment.
%T - the cost of the optimal assignment.
%s.t. T = trace(A(C,:)) is minimized over all possible assignments.

% Adapted from the FORTRAN IV code in Carpaneto and Toth, "Algorithm 548:
% Solution of the assignment problem [H]", ACM Transactions on
% Mathematical Software, 6(1):104-111, 1980.

% v1.0  96-06-14. Niclas Borlin, niclas@cs.umu.se.
%                 Department of Computing Science, Ume? University,
%                 Sweden. 
%                 All standard disclaimers apply.

% A substantial effort was put into this code. If you use it for a
% publication or otherwise, please include an acknowledgement or at least
% notify me by email. /Niclas

[m,n]=size(A);

if (m~=n)
    error('HUNGARIAN: Cost matrix must be square!');
end

% Save original cost matrix.
orig=A;

% Reduce matrix.
A=hminired(A);

% Do an initial assignment.
[A,C,U]=hminiass(A);

% Repeat while we have unassigned rows.
while (U(n+1))
    % Start with no path, no unchecked zeros, and no unexplored rows.
    LR=zeros(1,n);
    LC=zeros(1,n);
    CH=zeros(1,n);
    RH=[zeros(1,n) -1];
    
    % No labelled columns.
    SLC=[];
    
    % Start path in first unassigned row.
    r=U(n+1);
    % Mark row with end-of-path label.
    LR(r)=-1;
    % Insert row first in labelled row set.
    SLR=r;
    
    % Repeat until we manage to find an assignable zero.
    while (1)
        % If there are free zeros in row r
        if (A(r,n+1)~=0)
            % ...get column of first free zero.
            l=-A(r,n+1);
            
            % If there are more free zeros in row r and row r in not
            % yet marked as unexplored..
            if (A(r,l)~=0 & RH(r)==0)
                % Insert row r first in unexplored list.
                RH(r)=RH(n+1);
                RH(n+1)=r;
                
                % Mark in which column the next unexplored zero in this row
                % is.
                CH(r)=-A(r,l);
            end
        else
            % If all rows are explored..
            if (RH(n+1)<=0)
                % Reduce matrix.
                [A,CH,RH]=hmreduce(A,CH,RH,LC,LR,SLC,SLR);
            end
            
            % Re-start with first unexplored row.
            r=RH(n+1);
            % Get column of next free zero in row r.
            l=CH(r);
            % Advance "column of next free zero".
            CH(r)=-A(r,l);
            % If this zero is last in the list..
            if (A(r,l)==0)
                % ...remove row r from unexplored list.
                RH(n+1)=RH(r);
                RH(r)=0;
            end
        end
        
        % While the column l is labelled, i.e. in path.
        while (LC(l)~=0)
            % If row r is explored..
            if (RH(r)==0)
                % If all rows are explored..
                if (RH(n+1)<=0)
                    % Reduce cost matrix.
                    [A,CH,RH]=hmreduce(A,CH,RH,LC,LR,SLC,SLR);
                end
                
                % Re-start with first unexplored row.
                r=RH(n+1);
            end
            
            % Get column of next free zero in row r.
            l=CH(r);
            
            % Advance "column of next free zero".
            CH(r)=-A(r,l);
            
            % If this zero is last in list..
            if(A(r,l)==0)
                % ...remove row r from unexplored list.
                RH(n+1)=RH(r);
                RH(r)=0;
            end
        end
        
        % If the column found is unassigned..
        if (C(l)==0)
            % Flip all zeros along the path in LR,LC.
            [A,C,U]=hmflip(A,C,LC,LR,U,l,r);
            % ...and exit to continue with next unassigned row.
            break;
        else
            % ...else add zero to path.
            
            % Label column l with row r.
            LC(l)=r;
            
            % Add l to the set of labelled columns.
            SLC=[SLC l];
            
            % Continue with the row assigned to column l.
            r=C(l);
            
            % Label row r with column l.
            LR(r)=l;
            
            % Add r to the set of labelled rows.
            SLR=[SLR r];
        end
    end
end

% Calculate the total cost.
T=sum(orig(logical(sparse(C,1:size(orig,2),1))));


function A=hminired(A)
%HMINIRED Initial reduction of cost matrix for the Hungarian method.
%
%B=assredin(A)
%A - the unreduced cost matris.
%B - the reduced cost matrix with linked zeros in each row.

% v1.0  96-06-13. Niclas Borlin, niclas@cs.umu.se.

[m,n]=size(A);

% Subtract column-minimum values from each column.
colMin=min(A);
A=A-colMin(ones(n,1),:);

% Subtract row-minimum values from each row.
rowMin=min(A')';
A=A-rowMin(:,ones(1,n));

% Get positions of all zeros.
[i,j]=find(A==0);

% Extend A to give room for row zero list header column.
A(1,n+1)=0;
for k=1:n
    % Get all column in this row. 
    cols=j(k==i)';
    % Insert pointers in matrix.
    A(k,[n+1 cols])=[-cols 0];
end


function [A,C,U]=hminiass(A)
%HMINIASS Initial assignment of the Hungarian method.
%
%[B,C,U]=hminiass(A)
%A - the reduced cost matrix.
%B - the reduced cost matrix, with assigned zeros removed from lists.
%C - a vector. C(J)=I means row I is assigned to column J,
%              i.e. there is an assigned zero in position I,J.
%U - a vector with a linked list of unassigned rows.

% v1.0  96-06-14. Niclas Borlin, niclas@cs.umu.se.

[n,np1]=size(A);

% Initalize return vectors.
C=zeros(1,n);
U=zeros(1,n+1);

% Initialize last/next zero "pointers".
LZ=zeros(1,n);
NZ=zeros(1,n);

for i=1:n
    % Set j to first unassigned zero in row i.
	lj=n+1;
	j=-A(i,lj);

    % Repeat until we have no more zeros (j==0) or we find a zero
	% in an unassigned column (c(j)==0).
    
	while (C(j)~=0)
		% Advance lj and j in zero list.
		lj=j;
		j=-A(i,lj);
	
		% Stop if we hit end of list.
		if (j==0)
			break;
		end
	end

	if (j~=0)
		% We found a zero in an unassigned column.
		
		% Assign row i to column j.
		C(j)=i;
		
		% Remove A(i,j) from unassigned zero list.
		A(i,lj)=A(i,j);

		% Update next/last unassigned zero pointers.
		NZ(i)=-A(i,j);
		LZ(i)=lj;

		% Indicate A(i,j) is an assigned zero.
		A(i,j)=0;
	else
		% We found no zero in an unassigned column.

		% Check all zeros in this row.

		lj=n+1;
		j=-A(i,lj);
		
		% Check all zeros in this row for a suitable zero in another row.
		while (j~=0)
			% Check the in the row assigned to this column.
			r=C(j);
			
			% Pick up last/next pointers.
			lm=LZ(r);
			m=NZ(r);
			
			% Check all unchecked zeros in free list of this row.
			while (m~=0)
				% Stop if we find an unassigned column.
				if (C(m)==0)
					break;
				end
				
				% Advance one step in list.
				lm=m;
				m=-A(r,lm);
			end
			
			if (m==0)
				% We failed on row r. Continue with next zero on row i.
				lj=j;
				j=-A(i,lj);
			else
				% We found a zero in an unassigned column.
			
				% Replace zero at (r,m) in unassigned list with zero at (r,j)
				A(r,lm)=-j;
				A(r,j)=A(r,m);
			
				% Update last/next pointers in row r.
				NZ(r)=-A(r,m);
				LZ(r)=j;
			
				% Mark A(r,m) as an assigned zero in the matrix . . .
				A(r,m)=0;
			
				% ...and in the assignment vector.
				C(m)=r;
			
				% Remove A(i,j) from unassigned list.
				A(i,lj)=A(i,j);
			
				% Update last/next pointers in row r.
				NZ(i)=-A(i,j);
				LZ(i)=lj;
			
				% Mark A(r,m) as an assigned zero in the matrix . . .
				A(i,j)=0;
			
				% ...and in the assignment vector.
				C(j)=i;
				
				% Stop search.
				break;
			end
		end
	end
end

% Create vector with list of unassigned rows.

% Mark all rows have assignment.
r=zeros(1,n);
rows=C(C~=0);
r(rows)=rows;
empty=find(r==0);

% Create vector with linked list of unassigned rows.
U=zeros(1,n+1);
U([n+1 empty])=[empty 0];


function [A,C,U]=hmflip(A,C,LC,LR,U,l,r)
%HMFLIP Flip assignment state of all zeros along a path.
%
%[A,C,U]=hmflip(A,C,LC,LR,U,l,r)
%Input:
%A   - the cost matrix.
%C   - the assignment vector.
%LC  - the column label vector.
%LR  - the row label vector.
%U   - the 
%r,l - position of last zero in path.
%Output:
%A   - updated cost matrix.
%C   - updated assignment vector.
%U   - updated unassigned row list vector.

% v1.0  96-06-14. Niclas Borlin, niclas@cs.umu.se.

n=size(A,1);

while (1)
    % Move assignment in column l to row r.
    C(l)=r;
    
    % Find zero to be removed from zero list..
    
    % Find zero before this.
    m=find(A(r,:)==-l);
    
    % Link past this zero.
    A(r,m)=A(r,l);
    
    A(r,l)=0;
    
    % If this was the first zero of the path..
    if (LR(r)<0)
        ...remove row from unassigned row list and return.
        U(n+1)=U(r);
        U(r)=0;
        return;
    else
        
        % Move back in this row along the path and get column of next zero.
        l=LR(r);
        
        % Insert zero at (r,l) first in zero list.
        A(r,l)=A(r,n+1);
        A(r,n+1)=-l;
        
        % Continue back along the column to get row of next zero in path.
        r=LC(l);
    end
end


function [A,CH,RH]=hmreduce(A,CH,RH,LC,LR,SLC,SLR)
%HMREDUCE Reduce parts of cost matrix in the Hungerian method.
%
%[A,CH,RH]=hmreduce(A,CH,RH,LC,LR,SLC,SLR)
%Input:
%A   - Cost matrix.
%CH  - vector of column of 'next zeros' in each row.
%RH  - vector with list of unexplored rows.
%LC  - column labels.
%RC  - row labels.
%SLC - set of column labels.
%SLR - set of row labels.
%
%Output:
%A   - Reduced cost matrix.
%CH  - Updated vector of 'next zeros' in each row.
%RH  - Updated vector of unexplored rows.

% v1.0  96-06-14. Niclas Borlin, niclas@cs.umu.se.

n=size(A,1);

% Find which rows are covered, i.e. unlabelled.
coveredRows=LR==0;

% Find which columns are covered, i.e. labelled.
coveredCols=LC~=0;

r=find(~coveredRows);
c=find(~coveredCols);

% Get minimum of uncovered elements.
m=min(min(A(r,c)));

% Subtract minimum from all uncovered elements.
A(r,c)=A(r,c)-m;

% Check all uncovered columns..
for j=c
    % ...and uncovered rows in path order..
    for i=SLR
        % If this is a (new) zero..
        if (A(i,j)==0)
            % If the row is not in unexplored list..
            if (RH(i)==0)
                % ...insert it first in unexplored list.
                RH(i)=RH(n+1);
                RH(n+1)=i;
                % Mark this zero as "next free" in this row.
                CH(i)=j;
            end
            % Find last unassigned zero on row I.
            row=A(i,:);
            colsInList=-row(row<0);
            if (length(colsInList)==0)
                % No zeros in the list.
                l=n+1;
            else
                l=colsInList(row(colsInList)==0);
            end
            % Append this zero to end of list.
            A(i,l)=-j;
        end
    end
end

% Add minimum to all doubly covered elements.
r=find(coveredRows);
c=find(coveredCols);

% Take care of the zeros we will remove.
[i,j]=find(A(r,c)<=0);

i=r(i);
j=c(j);

for k=1:length(i)
    % Find zero before this in this row.
    lj=find(A(i(k),:)==-j(k));
    % Link past it.
    A(i(k),lj)=A(i(k),j(k));
    % Mark it as assigned.
    A(i(k),j(k))=0;
end

A(r,c)=A(r,c)+m;