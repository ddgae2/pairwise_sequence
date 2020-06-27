% pair-wise alignment amino acid/DNA and correlation analysis 
% copyright UCDAVIS 2015
% David Gae
% eg. align_pairwise('seq1',score,0.8,0,-0.8,0)
function align_pairwise = align_pairwise(seq1,score,a,b,c,d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% david gae 2/21/12 UC DAVIS. 
%correlation analysis of position- i,j amino acids in multiple sequence alignment. 
% for string use '%c'
%physiochemical values are based following papers.
%1. Volume (Chothia 1984).
%2. Annu. Rev. Biochem. (1984) 53 537-572.
%2. Average accessibility surface area (Janin et al. 1978)
%3. Accessible surface area in the standard state (Rose et al. 1985) Science (1975) 229 834-838


[value, amin]= textread(score, '%f %c', 20);
%[value, amin]= textread('janin.txt', '%f %c', 20)
%[value, amin]= textread('rose.txt', '%f %c', 20)
%convert column to row string
amino=transpose(amin);
maxamino= length(amino);

fileID1 = fopen('heatmap.csv','w');
%For nrow,ncol of sequences
seq1 = fopen('seq1','rt');          		    
%seq2 = fopen('seq2','rt');

%location of seq char
seq_a= char();
%seq_b= char();

%initalize
nrows=0;ncols=0;cna=0;
nrows1=0;ncols1=0;cna1=0;
cna=0;
%store seq_a 
rseq1 = fgetl(seq1);

	while (ischar(rseq1))
		
	cna = cna+1;

			seq_a(cna,1:length(rseq1)) = rseq1; %store string of eachline
			rseq1= fgetl(seq1);
			[nrows,ncols] =size(rseq1);
	end;
fclose(seq1);


%initialize
z=0;cnt=0;
fi = 0;
%fj = 0;
% determine the length of rows.
m = 0;
m = length(seq_a(1,:));


%sequence i determine sum(x^2) and sum(x)^2

num =0;   % counter
	for cnt = 1:length(seq_a(:,m))    % loop for number of rows
			num = num +1 ;            % count
					for i = 1:m       % per column/each residue 
			
							for z = 1:maxamino				% maxamino = length of 20 amino acid 

							if seq_a(num,i,:) == amino(z)	% match "seq_a mat" to amino acid value
								
								%use for sum[x*y] step
								fi(num,i) = value(z);		%store amino value of eachline into fi
								
								%use for correlation (x^2) and (x)^2
								fidoub(num,i)= value(z)*value(z);		%square of amino acid value
								
								%disp(fi);
								%use for sum[x*y] step
								fisum(i) = sum(fi(:,i));
								fidoubsum(i) = sum(fidoub(:,i));	    %sum(x^2)
								fisumdoub(i) = sum(fi(:,i))^2;          %sum(x)^2
								
							end
					end
			end
	end


%%%% determine for sum of expected value of (x*y) = fi(ii)*fj(jj)
%initialize

ii=0;
jj=0;

for ii = 1:m %  length of rows 

		pval1 = fi(:,ii);    % column value at position i

					for jj = 1:m   % length of rows
				
								if ii ~= jj && jj > ii					% count only  i !== j
									
									pval2 = fi(:,jj);					% column value at position j

									%	disp(ii); disp(jj);             %check positon i and j 
									%   disp(pval1); disp(pval2);       %check physiochemical values at position i and j

									for cnt = 1:length(pval1(:,1))      % count number of rows, position i
											
										meanxy(cnt) = pval1(cnt,1) * pval2(cnt,1);
															
									%	disp(meanxy(cnt))
									%	addxy(ii,jj) = sum(meanxy);

										end
										addxy(ii,jj) = sum(meanxy);		%meanxy in i,j

								%   disp(addxy);

							end							
	
					end
end


% Pearsons method of determination%
% *****unbiased covariance uses (N-1) of only the co-variance******
%  sum(xy) - sum(x)*sum(y)/N-1
%  sqrt ( (sum(x^2) - ((sum(x)^2)/N)) - sum((y^2) - ((sum(y)^2)/N)) )
%determination  covariance and correlation of position i and j
%determination cov using matlab function as comparsion (my self-check)
%i=0;
%j=0;
%	for i = 1:m 
%			for j = 2:m
%	covijmat = cov(fi(:,i),fj(:,j))      % using matlab function covariance
%	corijmat = corrcoef(fi(:,i),fj(:,j))  % using matlab function correlation
%	pause;
%			end
%	end 
%

%initialize
ii=0;
jj=0;

	for ii = 1:m
			for jj = 1:m
					if ii ~= jj && jj > ii
% disp(ii);disp(jj);

% *****unbiased covariance uses (N-1)******
%  sum(xy) - sum(x)*sum(y)/N-1
					covij(ii,jj)= addxy(ii,jj) - ((fisum(ii)*fisum(jj))/length(pval1(:,1))-1);

					%disp(addxy(ii,jj) - (fisum(ii)*fisum(jj))/length(pval1(:,1))-1);
					%pause;
					% matlab set double precisions, when multiplying floats.

% *****correlation values  (N-1)********
					corij1(ii,jj) = covij(ii,jj)/sqrt( (fidoubsum(ii)-(fisumdoub(ii)/length(pval1(:,1))) ) * (fidoubsum(jj)-(fisumdoub(jj)/length(pval1(:,1)))) );
					% double creates real and imaginary a+bi.correct only real values. 
					corij2(ii,jj) = real(corij1(ii,jj));
			

					end
			end
	end

% clean up corij-matrix of NaN=0 and Inf=1, Inf=-1 values, value1 = great than 1. 
% due to im-precision in floating point calculation. 
		NaN1 = find(isnan(corij2));
		corij2(NaN1)=zeros(size(NaN1));
				inf1 = find(isinf(corij2));
				corij2(inf1)=ones(size(inf1));
						value1= find(corij2>1);
						corij2(value1)=ones(size(value1));


%t-test- significance of bell-shape samples.
%initialize
	ii=0;
	jj=0;

			for ii = 1:m
				for jj = 1:m
					if ii ~= jj && jj > ii  %  i !=j

					% disp(ii);disp(jj);

					% student t-test and m-2 is length of aligned sequences -2 (eg.298-2)
					ttest(ii,jj) = sqrt(real(corij1(ii,jj).^2))*sqrt((m-2)/(1-real(corij1(ii,jj).^2)));
					trealtest(ii,jj) = real(ttest(ii,jj));					% t-test values.
			
					% disp(treal-test);
					% sigcorij3 most significant matrix array.
						%if trealtest(ii,jj) <= sqrt(real(corij1(ii,jj)).^2)
						%sigcorij3(ii,jj) = sqrt(real(corij1(ii,jj)).^2);						
						%else
						%sigcorij3(ii,jj)= 0;
						%end

					end
				end
			end



%%%correlation coefficient %%%%%%%%%%%%%%%%%%%
%%%conserved position analysis%%%%%%%%%%%%%%%%
%%% determination of (xi-meanx)
%%%( x-meanx)^2/numtotal-1 (variance values)

%initialize
		ii=0;
		cnt1=0;
		
		for ii = 1:m %  length of rows 
			cval1 = fi(:,ii);						% column value at position i
			csum1=0;								% starter counter and reset
			
				for cnt1 = 1:length(cval1(:,1))         % count number of rows, position i
						if cval1(cnt1,1) ~= 0.000;			% if column is zero, do not count. 
						csum1 = csum1 + 1;				%  mean counter only
						end

				end										% finish counter

					%disp(csum1);						% display counter total for mean.
		totalsum1 =sum(fi(:,ii));						% totalsum of ith position
		meanx(1:length(fi(:,ii)),ii)= totalsum1/csum1;  % mean of ith  within matrix (column via position).
					
		ctotal1(:,ii)=(csum1);							% store mean counter of column of i

%required for correlation coefficient of set with "real values"
%use of variance cal
					for cnt1 = 1:length(cval1(:,1))			% count number of rows, position i
								if fi(cnt1,ii) == 0.000		% value of row,col- position i
								meanx(cnt1,ii) = 0.000;		% meanx will be zero at position. 
								%disp(meanx);
								end
					end

Svx(:,ii) = (fi(:,ii)-meanx(:,ii));								% (x-meanx) value for covariance (x-meanx)*(y-meany)***
SSx(:,ii) = ((fi(:,ii)-meanx(:,ii)).^2);						% SS(x-meanx)^2 values
%SSx(:,ii) = sum((fi(:,ii)-meanx(:,ii)).^2);					% sum(x-meanx)^2 values
varx(ii) = sum((fi(:,ii)-meanx(:,ii)).^2)/(ctotal1(:,ii)-1);	% SS(x-meanx)^2/csum1-1 (variance values)

		end

%%%% determination of (yi-meany)
%%%( y-meany)^2/numtotal-1 (variance values)

%initialize
		jj=0;
		cnt2=0;

		for jj = 1:m %  length of rows 
			cval2 = fi(:,jj);						% column value at position j
			csum2=0;								% starter counter and reset

				for cnt2 = 1:length(cval2(:,1))				% count number of rows, position j
						if cval2(cnt2,1) ~= 0.000;				% if column is zero, do not count. 
						csum2 = csum2 + 1;					%  mean counter only
						end

				end											% finish counter

					%disp(csum2);							% display counter total for mean
		totalsum2 =sum(fi(:,jj));							%totalsum of jth position
		meany(1:length(fi(:,jj)),jj)= totalsum2/csum2;		%mean of jth position within matrix (column via position).

		ctotal2(:,jj)=(csum2);								% store mean counter of column of j

%required for correlation coefficient of set with "real values"
%use of variance cal
						for cnt2 = 1:length(cval2(:,1))			% count number of rows, position j
									if fi(cnt2,jj) == 0.000		% value of row,col- position j
									meany(cnt2,jj) = 0.000;		% meanx will be zero at position. 
									%disp(meany);
									end
						end
Svy(:,jj) =  (fi(:,jj) - meany(:,jj));								% (y-meany) value for covariance(x-meanx)*(y-meany)***
SSy(:,jj) =  ((fi(:,jj)-meany(:,jj)).^2);							% (y-meany)^2 values
%SSy(:,jj) =  sum((fi(:,jj)-meany(:,jj)).^2);						% sum(y-meany)^2 values
vary(jj) = sum((fi(:,jj)-meany(:,jj)).^2)/(ctotal2(:,jj)-1);		%( y-meany)^2/csum2-1 (variance values)
		
		end

%%%%values for Saddxy%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% X-meanX and Y-MeanY%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii=0;
jj=0;
cnt=0;

				for ii = 1:m %  length of rows 

%floating point precision issue, just round to the best round number with three significant digits. 

									Smeanx = round(Svx(:,ii)*1000)/1000;								% column (x-meanx) value at position i
									%Sdifx  =  round(SSx(:,ii)*1000)/1000;								% column (x-meanx)^2 value at position i, NOT USE.
										for jj = 1:m   % length of rows

												if ii ~= jj && jj > ii					% count only  i !== j

														Smeany = round(Svy(:,jj)*1000)/1000;			% column (y-meany) value at position j
														Sdify = round(SSy(:,jj)*1000)/1000;				% column (y-meany)^2 value at position j
														Sdify1 = round(SSy(:,ii)*1000)/1000;			% column (x-meanx)*2 value at postion i (fix)

															%disp(ii); disp(jj);						%check positon i and j 
														    %disp(Smeanx); disp(Smeany);				%check Smeanx and Smeany values at position i and j
															%pause();
													
														csum3=0;										%count number of rows with non-gaps for scale factor.

														for cnt = 1:length(Svy(:,1))					% count number of rows, position i
														
															Sxy(cnt) = Smeanx(cnt,1) * Smeany(cnt,1);   % using Smeanx and Smeany previous-(x-meanx)*(y-meany)

																if Sxy(cnt) ~= 0.000				    % if column is zero, do not count. 
																	csum3 = csum3 + 1;					%  counter of non-gaps for scaling. 
																end
																
																%disp(cnt);
																%disp(ii); disp(jj);
																%disp(Sxy(cnt));						%check answers
																%pause();

														%condition statement for aligned seq only variance.
																if Sxy(cnt) == 0.000
																	Sdify1(cnt,1)= 0;					%set (x-meanx)^2 = 0 in column values
																end

																if Sxy(cnt) == 0.000				
																	Sdify(cnt,1) = 0;					%set (y-meany)^2 = 0 in column values
																end
													    end
																									
											
											%%values for linear sqt fit values.
											%%below%%

											Saddxy(ii,jj) = sum(Sxy);						%Sumxy of i,j- values for covariance.
											Sdifx2(ii,jj) = sum(Sdify1);					%Sumx of x-meanx column
											Sdify2(ii,jj) = sum(Sdify);						%Sumy of y-meany column
											totalcsum3(ii,jj) = sum(csum3);					%Scale factor used in linear fit, R^2. 

											   end							

									  end

						end




%FOR conserved positions only where, variance of x and y is zero.
%least sqt fit is prefect.
%SR2 equal to 1.
ii=0;
jj=0;

		for ii = 1:m %  length of rows	
						for jj = 1:m   % length of rows
				
								if ii ~= jj && jj > ii					% count only  i !== j
								if (fi(1,ii)==fi(:,ii))					% print out if ROW value of fi = COLUMN value of fi 
								if (fi(1,jj)==fi(:,jj))					% print out if ROW value of fj = COLUMN value of fi
																		% logical values of TRUE =1 or FALSE = 0
										Saddxy(ii,jj) =1;				%NEW Sumxy of i,j- values for NEW covariance
										Sdifx2(ii,jj) =1;				%NEW Sumx of x-meanx column
										Sdify2(ii,jj) =1;				%NEW Sumy of y-meany column
						
								end							
								end
								end
						end
		end


%coefficient of determination R2
%%least sqaure fit coeff b= Saddxy/SSx
%%least square fit coeff b'' = Saddxy/SSy 
%% least squart fit coeff  r2= Saddxy^2/(SSx*SSy) 
%% r2 = b*b''

%initialize
cnt=0;
ii=0;
jj=0;
		for ii = 1:m %  length of rows 

				for jj = 1:m %  length of rows 
		
							  if ii ~= jj && jj > ii  %i !=j

								R2(ii,jj) = (Saddxy(ii,jj).^2)/((Sdifx2(ii,jj))*(Sdify2(ii,jj)));
	
						%condition 1
						%covariance is zero-- no relationship,hence least sqt fit is zero.
								
										if Saddxy(ii,jj) == 0.000
											R2(ii,jj) = 0.000;
											R2(ii,jj) = round(R2(ii,jj)*1000)/1000;
										%else SRadjust(ii,jj) = R2(ii,jj);
										end

						%condition 2, round down to integer.
										if R2(ii,jj) >= 1.00 
											SRadjust(ii,jj) = floor(R2(ii,jj));
										else  SRadjust(ii,jj) = R2(ii,jj); % value 1 or greater return to 1.
										end
							  end
				 end
		 end

%initialize of R^2
%my modification or cleanup to get close to crasp. 
cnt=0;
ii=0;
jj=0;
		for ii = 1:m %  length of rows 

			for jj = 1:m %  length of rows 

					if ii ~= jj && jj > ii  %i !=j


						%condition 3, scale factor for linear fit for values less than 7, better statistics and linear fit.  
						
							if totalcsum3(ii,jj) <= 7 && totalcsum3(ii,jj) > 1                %must for scaling
									SRadjust(ii,jj) = SRadjust(ii,jj)*(totalcsum3(ii,jj)-1);  %scaling SRadjust values. 
									else SRadjust(ii,jj) = SRadjust(ii,jj);
							end
						
						%condition 4,csum3 is small size, but correlation score is 0.75, round up to 1.
							if SRadjust(ii,jj) >= 0.75 && totalcsum3(ii,jj) <= 7 && totalcsum3(ii,jj) > 1
									SRadjust(ii,jj) = 1;
									else SRadjust(ii,jj) = SRadjust(ii,jj);  % do not use.
							end

					end
			 end
		end



%initialize
%how crasp does it. 
% position correlation analysis. 
cnt=0;
ii=0;
jj=0;
			for ii = 1:m %  length of rows 

				for jj = 1:m %  length of rows 

					if ii ~= jj && jj > ii  %i !=j	
						
								Svar(ii,jj) = (Saddxy(ii,jj));				                         %var(x,y) = 2cov(x,y)
							    SCRASP(ii,jj) = Svar(ii,jj)/sqrt(Sdifx2(ii,jj)*Sdify2(ii,jj));		 %correlation at positons x,y cov(x,y)/sqrt(varx+vary)
								
									%Svar(ii,jj) = ((varx(ii)+vary(jj))+(2*Saddxy(ii,jj)));				 %var(x,y) = var(x) +var(y) + 2cov(x,y), another statstics 

					%% scale values for correlation coefficient if less than 7. 
						if totalcsum3(ii,jj) <= 7 && totalcsum3(ii,jj) > 1   
								  Svar(ii,jj) = (Saddxy(ii,jj))*(totalcsum3(ii,jj)-1);				% scale values for values less than 7 better fit.
								  SCRASP(ii,jj) = Svar(ii,jj)/sqrt(Sdifx2(ii,jj)*Sdify2(ii,jj));
							else SCRASP(ii,jj) =SCRASP(ii,jj);
						end

							%SCRASP(ii,jj) = Svar(ii,jj)/(Sdifx2(ii,jj)+Sdify2(ii,jj));				%correlation at positons x,y cov(x,y)/sqrt(varx+vary)
							
						%condition2, covariance is zero-- no relationship,hence least sqt fit is zero.
							if Saddxy(ii,jj) == 0.000 && Saddxy(ii,jj) <= 0.01
							SCRASP(ii,jj) = 0.000;
							SCRASP(ii,jj) = round(SCRASP(ii,jj)*1000)/1000;
							end

						%condition3,csum3 is small size, but correlation score is 0.75, round up to 1.
							if SCRASP(ii,jj) >= 0.75 && totalcsum3(ii,jj) <= 7 && totalcsum3(ii,jj) > 1
							SCRASP(ii,jj) = 1;
							else SCRASP(ii,jj) = SCRASP(ii,jj);  % do not use.
							end
						%condition4,csum3 is small size, but correlation score is 0.75, round up to 1.
							if SCRASP(ii,jj) <= - 0.75 && totalcsum3(ii,jj) <= 7 && totalcsum3(ii,jj) > 1
							SCRASP(ii,jj) = -1;
							else SCRASP(ii,jj) = SCRASP(ii,jj);  % do not use.
							end

						  
					%clean up dont need this. 
					NaN1 = find(isnan(SCRASP));
					SCRASP(NaN1)=zeros(size(NaN1));
					inf1 = find(isinf(SCRASP));
					SCRASP(inf1)=zeros(size(inf1));
					
					end
				end
			end


SCRASP(a>SCRASP>b)=0;
SCRASP(c<SCRASP<d)=0;


fprintf(fileID1, [repmat('%d,', 1, size(SCRASP, 2)), '\n'], SCRASP');
        %Heatmap via MATLAB
        lowestValue = min(SCRASP(SCRASP(:)>0));
        highestValue = max(SCRASP(:));
        imagesc(SCRASP);
        cmap = jet(256);
        colormap(cmap);
        caxis(gca,[lowestValue-2/256, highestValue]);
        % lower than minimum colored as white:
        cmap(1,:)=[1,1,1];
        colormap(cmap)
        colorbar
fclose(fileID1);
