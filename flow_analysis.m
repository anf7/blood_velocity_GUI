function submatrix = flow_analysis(imstack,meanpixval,startstop,downsampleimg,executemask)

[z, r, c] = size(imstack);


steplimit = 4;
maxsteprad = 45;
stepradsize = 20;

ignore_areas_larger_than = 350;
ignore_continuous_areas_smaller_than = 1;
ignore_discontinous_areas_smaller_than = 3;
ignore_area_changes_greater_than = 4;
ignore_distance_changes_greater_than = 2.5;
min_continuous_step_difference = 3;
min_discontinuous_step_difference = 2;
max_angle_of_divergence = 30; % degrees
max_angle_of_divergence = max_angle_of_divergence*pi/180;


if downsampleimg
    maxsteprad = 25;
    stepradsize = 15;
    ignore_areas_larger_than = 200;
end



startr = startstop(1);
stopr = startstop(2);
startc = startstop(3);
stopc = startstop(4);



NFFT = 2^nextpow2(z);
zeroedstack = zeros(z,r,c);

for m = 1:z
    zeroedstack(m,:,:) = squeeze(imstack(m,:,:)) - meanpixval;
end

clear meanpixval
      

stdzeroedstack = zeros(r,c);
for n = 1:c
    for m = 1:r
        stdzeroedstack(m,n) = std(zeroedstack(:,m,n));
    end
end


reversezeroedstack = zeros(NFFT,r,c);
for m = 1:z
    reversezeroedstack(NFFT-m+1,:,:) = zeroedstack(m,:,:);
end

zeroedFT = fft(zeroedstack,NFFT,1);
reversezeroedFT = fft(reversezeroedstack,NFFT,1);

clear reversezeroedstack
clear zeroedstack
movement_matrix = zeros(r,c,3);


for n = startc:stopc
    
   
    for m = startr:stopr
        
        if ~executemask(m,n)
            continue
        end
        

        
        if sum(imstack(:,m,n)) == 0
            continue
        end

        
%% Continuous areas
        rindex = m;
        cindex = n;
        iimg = int8(127*ones(r,c));
        
        iimg(m,n) = 0;
        smask = iimg == 0;
        zeroedFTvec = zeroedFT(:,m,n);
        stdF = stdzeroedstack(m,n);
        
        [zsmask,zring,iimg] = dilatefill_mat_zerotest(iimg,smask,reversezeroedFT,zeroedFTvec,stdzeroedstack,stdF,NFFT,rindex,cindex,0,ignore_areas_larger_than);
        if isempty(zsmask)
            continue
        end
        totarea = sum(zsmask(:));
           
            
        sumsmask = totarea;
        oldsumsmask = sumsmask;
   
        [zrpos,zcpos] = centerofmass(zsmask);
        rpos = zrpos;
        cpos = zcpos;
        rupos = zrpos;
        cupos = zcpos;
             
        dist = 0;
        dilatedone = 0;
        stepnum = int8(0);
        while dilatedone == 0
            if stepnum ~= 0
                oldsumsmask = sumsmask;
                [smask,ring,iimg] = dilatefill_mat(iimg,smask,reversezeroedFT,zeroedFTvec,stdzeroedstack,stdF,NFFT,rindex,cindex,stepnum);
                sumsmask = sum(smask(:));
            else
                smask = zsmask;
                ring = zring;
                prior_rpos = zrpos;
                prior_cpos = zcpos;
            end
            
            if stepnum ~= 0 && (sumsmask > ignore_areas_larger_than || sumsmask < ignore_continuous_areas_smaller_than ||...
                    oldsumsmask > ignore_area_changes_greater_than*sumsmask ||...
                    oldsumsmask < sumsmask/ignore_area_changes_greater_than) 
                ustepnum = stepnum - 1;
                rupos = rpos;
                cupos = cpos;
                break
            end
            
 
            nextsmask = (ring & (iimg == (stepnum + 1)));
            
            if (isequal(nextsmask,false(r,c)) || stepnum == steplimit)
                ustepnum = stepnum;
                if stepnum ~= 0
                    [rupos,cupos] = centerofmass(smask);
                    dist = dist + ((rupos - prior_rpos)^2 + (cupos - prior_cpos)^2)^0.5;
                end
                totarea = totarea + sumsmask;
                dilatedone = 1;
            else
                [rindex,cindex] = find(nextsmask);
                if stepnum ~= 0
                    [rpos,cpos] = centerofmass(smask);    
                    dist = dist + ((rpos - prior_rpos)^2 + (cpos - prior_cpos)^2)^0.5;
                    prior_rpos = rpos;
                    prior_cpos = cpos;
                end
                totarea = totarea + sumsmask;
                smask = nextsmask;
                stepnum = stepnum + 1;
            end
        end
        
        
        dnstepnum = 0;
        rdpos = zrpos;
        cdpos = zcpos;
        if ~isempty(find((zring.*(iimg == -1)),1)) 
            smask = (zring.*(iimg == -1));
            [rindex,cindex] = find(smask);
            stepnum = int8(-1);
            prior_rpos = zrpos;
            prior_cpos = zcpos;
            rpos = prior_rpos;
            cpos = prior_cpos;
            oldsumsmask = sum(zsmask(:));
            
            dilatedone = 0;
            while dilatedone == 0
                if stepnum ~= -1
                    oldsumsmask = sumsmask;
                end

                [smask,ring,iimg] = dilatefill_mat(iimg,smask,reversezeroedFT,zeroedFTvec,stdzeroedstack,stdF,NFFT,rindex,cindex,stepnum);
                sumsmask = sum(smask(:));
                
                if sumsmask > ignore_areas_larger_than || sumsmask < ignore_continuous_areas_smaller_than ||...
                        oldsumsmask > ignore_area_changes_greater_than*sumsmask ||...
                        oldsumsmask < sumsmask/ignore_area_changes_greater_than
                    dnstepnum = stepnum + 1;
                    rdpos = rpos;
                    cdpos = cpos;
                    break
                end

                nextsmask = (ring & (iimg == (stepnum - 1)));
              
                if isequal(nextsmask,false(r,c)) || stepnum == -1*steplimit
                    dnstepnum = stepnum;
                    [rdpos,cdpos] = centerofmass(smask);
                    dist = dist + ((rdpos - prior_rpos)^2 + (cdpos - prior_cpos)^2)^0.5;
                    totarea = totarea + sumsmask;
                    dilatedone = 1;
                else
                    [rindex,cindex] = find(nextsmask);
                    [rpos,cpos] = centerofmass(smask);
                    dist = dist + ((rpos - prior_rpos)^2 + (cpos - prior_cpos)^2)^0.5;
                    prior_rpos = rpos;
                    prior_cpos = cpos;
                    totarea = totarea + sumsmask;
                    smask = nextsmask;
                    stepnum = stepnum - 1;
                end
            end
        end
        
        if ((ustepnum - dnstepnum >= min_continuous_step_difference && totarea >= 4) ||...
                (ustepnum - dnstepnum >= 2 && totarea >= 25)) &&...
                ((double(rdpos - rupos))^2 + (double(cdpos - cupos))^2)^0.5 > 0.75*double(ustepnum - dnstepnum)
            
            timeshift = double(ustepnum - dnstepnum);         
            rshift = double(rdpos - rupos);
            cshift = double(cdpos - cupos);
                
            velocity = (rshift^2 + cshift^2)^0.5/timeshift;
            movement_matrix(m,n,:) = [rshift, cshift, velocity];
            continue
        end
        
        if sum(zsmask(:)) < 6
            continue
        end
        
        
                
%%  Discontinuous areas
        
        rdpos = zrpos;
        cdpos = zcpos;
        rupos = zrpos;
        cupos = zcpos;
        ustepnum = 0;
        dnstepnum = 0;
        

            
        tryagain = true;
        trycount = 0;
        while tryagain == true && trycount < 2
        tryagain = false;
        dist = 0;
        neg1rpos = 0;
        neg1cpos = 0;
        for stepnumdir = -1:2:1
            stepnumrad = 5;
            stepnum = int8(0);
            smask = zsmask;
            prior_rpos = zrpos;
            prior_cpos = zcpos;
            
            rpos = zrpos;
            cpos = zcpos;

            while true
                if stepnum ~= 0
                    if ~isequal(smask,false(r,c))
                        [rindex,cindex] = find(smask);
                        [smask, ~, iimg] = dilatefill_mat(iimg,smask,reversezeroedFT,zeroedFTvec,stdzeroedstack,stdF,NFFT,rindex,cindex,stepnum);
                        sumsmask = sum(smask(:));                        
                        if (((sumsmask >= ignore_discontinous_areas_smaller_than && sumsmask < ignore_areas_larger_than) ||...
                                (stepnum == 1 || stepnum == -1)) && sumsmask < sum(priorsmask(:))*ignore_area_changes_greater_than &&...
                                sumsmask > sum(priorsmask(:))/ignore_area_changes_greater_than)

                            [rpos,cpos] = centerofmass(smask);
                            rvec = rpos - prior_rpos;
                            cvec = cpos - prior_cpos;
                            new_dist = ((rpos - prior_rpos)^2 + (cpos - prior_cpos)^2)^0.5;


                            if stepnum > 1 || stepnum < -1 || (stepnum == 1 && neg1rpos ~= 0 && neg1cpos ~= 0)
                                if stepnum == 1
                                    prior_rvec = zrpos - neg1rpos;
                                    prior_cvec = zcpos - neg1cpos;
                                    prior_dist = neg1dist;
                                end
                                angle = anglediff(rvec,cvec,prior_rvec,prior_cvec);

                                if (angle > max_angle_of_divergence) ||...
                                        (prior_dist/new_dist < 1/ignore_distance_changes_greater_than) ||...
                                        (prior_dist/new_dist > ignore_distance_changes_greater_than)
                                    iimg(smask == 1) = 126;
                                    tryagain = true;
                                    trycount = trycount + 1;
                                    break
                                end
                            end
                            if stepnum == -1
                                neg1dist = new_dist;
                                neg1area = sumsmask;
                            end
                            dist = dist + new_dist;
                            prior_dist = new_dist;
                            prior_rvec = rvec;
                            prior_cvec = cvec;                        
                            prior_rpos = rpos;
                            prior_cpos = cpos;
                            if stepnum == -1
                                neg1rpos = rpos;
                                neg1cpos = cpos;
                            end
                    
                        elseif (stepnumrad < maxsteprad) 
                            smask = priorsmask;
                            stepnumrad = stepnumrad + stepradsize;
                            stepnum = stepnum - stepnumdir;
                        else
                            break
                        end
                    else
                        if (stepnumrad < maxsteprad) 
                            smask = priorsmask;
                            stepnumrad = stepnumrad + stepradsize;
                            stepnum = stepnum - stepnumdir;
                        else
                            break
                        end
                    end
                end
                
                if stepnumdir == 1
                    ustepnum = stepnum;
                    rupos = rpos;
                    cupos = cpos;
                else
                    dnstepnum = stepnum;
                    rdpos = rpos;
                    cdpos = cpos;
                end

                [rindex,cindex] = find(smask);
                rr1 = max(1,(min(rindex)-stepnumrad));
                rr2 = min(r,(max(rindex)+stepnumrad));
                cc1 = max(1,(min(cindex)-stepnumrad));
                cc2 = min(c,(max(cindex)+stepnumrad));
                
                priorsmask = smask;

                iimg = process_iimg_region(iimg,cc1,cc2,rr1,rr2,reversezeroedFT,zeroedFTvec,stdzeroedstack,stdF,NFFT);

                stepnum = stepnum + stepnumdir;
                if (ustepnum - dnstepnum) > steplimit %(stepnum == steplimit) || (stepnum == -1*steplimit)
                    break
                end
                smask = false(r,c);            
                smask(rr1:rr2,cc1:cc2) = (iimg(rr1:rr2,cc1:cc2) == stepnum);
            end
        end
        end
        
        if ~exist('neg1area','var')
            neg1area = 0;
        end
        if (ustepnum - dnstepnum >= min_discontinuous_step_difference) || (neg1area > 4)
            timeshift = double(ustepnum - dnstepnum); 
            rshift = double(rdpos - rupos);
            cshift = double(cdpos - cupos);
            
            velocity = dist/timeshift;
            
            movement_matrix(m,n,:) = [rshift, cshift, velocity];
            continue
        end


        
    end
end


submatrix = movement_matrix(startr:stopr,startc:stopc,:);
submatrix(isnan(submatrix)) = 0;


