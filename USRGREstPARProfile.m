function depthz=USRGREstPARProfile(PAR0, Chlaprof, z, USRpercentage, z1)
% function [PAR,USR,GR,depthz]=USRGREstPARProfile(PAR0, Chlaprof, z, USRpercentage, z1)
% 
% PAR0=100,Chlaprof=Chl(:,i1,j1),z=depth,
% USRpercentage=USRpercentage2(i1,Day), z1=[1:400]
% Chlaprof and z are input Chl profile

if isempty(z1)
    z1(:,1)=1:400; % the standard depths from 1 to 400m
end
if isempty(PAR0)
    PAR0=100; % If without PAR0 when estimating z1%, PAR0 is set as 100.
end

[a,b]=size(Chlaprof);if a==1 Chlaprof=Chlaprof'; end
[a,b]=size(z);if a==1 z=z'; end 
data_chl=[z,Chlaprof];data_chl(isnan(Chlaprof),:)=[];
if ~isempty(data_chl)
    
    Chla1=interp1(data_chl(:,1),data_chl(:,2),z1,'linear'); %Interpolation of Chl to the standard depths
    iv=find(z1>max(data_chl(:,1)));
    if ~isempty(iv)
        Chla1(iv)=Chla1(iv(1)-1); %Extrapolation of Chl to the standard depths
    end
    clear data_chl
    
    Kd490=0.01660 + 0.07242.* Chla1.^0.68955; % MM01
    KdUSR=0.91.*Kd490.^0.89; % Lin16
    KdUSR(Kd490<=0.1)=0.0062+1.16.*Kd490(Kd490<=0.1)-0.00018./Kd490(Kd490<=0.1); % Lin16
    
    p1=0.1+0.79*nanmean(Kd490(1:10)); %sea-surface Kd490 for calculation of p1 and p2
    p2=0.21-0.23*nanmean(Kd490(1:10));
    KdGR=p1+p2.*exp(-0.082.*z1);
    GR=PAR0.*(1-USRpercentage).*exp((-1).*KdGR.*z1);
    
    USR=GR*nan;USR(1)=PAR0*USRpercentage*exp((-1)*KdUSR(1)*z1(1));
    for i=2:length(z1)
        USR(i)=USR(i-1).*exp((-1).*KdUSR(i).*(z1(i)-z1(i-1)));
    end
    
    PAR=USR + GR;
    depthz=interp1(PAR,z1,1,'linear');
else
    depthz=nan;PAR=nan;USR=nan;GR=nan;
end

end
