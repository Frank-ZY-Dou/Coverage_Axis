%Citation
%@inproceedings{dou2022coverage,
%  title={Coverage Axis: Inner Point Selection for 3D Shape Skeletonization},
%  author={Dou, Zhiyang and Lin, Cheng and Xu, Rui and Yang, Lei and Xin, Shiqing and Komura, Taku and Wang, Wenping},
%  booktitle={Computer Graphics Forum},
%  volume={41},
%  number={2},
%  pages={419--432},
%  year={2022},
%  organization={Wiley Online Library}
% }
% https://github.com/Frank-ZY-Dou/Coverage_Axis

clc,clear 
close all

real_name = '01Ants-27';
samplesnumber=  1500;
scale_setting = 0.03; % adding offset 
path='./inputs/';
output_path = './outputs/';
real_name2 = real_name
medial_path ='01Ants-27-raw.txt';
medial_path = [path,medial_path];
disp(['processing, offset: ',scale_setting]);
disp(medial_path);
[medial_x,medial_y,medial_z,medial_r]=textread(medial_path,'v %n %n %n %n'); 
medial_r = abs(medial_r);
[m,n] = size(medial_x)

medial_num = length(medial_x); 
init_r = medial_r;

medial_r  = medial_r + scale_setting;
disp(['./inputs/',real_name2,'.off']);
[V,F,UV,C,N] = readOFF(['./inputs/',[real_name2,'.off']]);
outputname = [output_path ,real_name2, '_original.obj',];
writeOBJ(outputname, V,F,UV,[],N,[]);
[N,I,B,r]  = random_points_on_mesh(V, F, samplesnumber, 'Color', 'blue', 'MaxIter', 1500);

  
filename = [output_path,real_name,'_','samples.pts'];
% 

fileID=fopen(filename,'w')
fprintf(fileID,'%f %f %f\n',N');
fclose(fileID);


V = N;
boundary_x = V(:,1);
boundary_y = V(:,2);
boundary_z = V(:,3);
boundary_num = length(boundary_x);
D = zeros(boundary_num,medial_num); % if_cover info




for i = 1:boundary_num
   for j  = 1:medial_num
        if((medial_x(j) - boundary_x(i))^2 + (medial_y(j) - boundary_y(i))^2 +...
                (medial_z(j) - boundary_z(i))^2<= medial_r(j)* medial_r(j))
            D(i,j) = 1;
        end    
    end
    
end
lambda = 1;
%first easy model
f =  ones(1,medial_num); %代表这是求和
A =  -lambda*D;
b =  -ones(boundary_num,1)*1;%here ,we fix p_i 
lb = zeros(medial_num,1);
ub = ones(medial_num,1);
iint = [1:medial_num];
tic;
[x,fval]=intlinprog(f,iint,A,b,[],[],lb,ub);
toc;
disp('min_number:');
disp(fval);
outputname = [output_path ,real_name2, '_original_', num2str(samplesnumber),'to',num2str(fval),'_scale',num2str(scale_setting),'.obj'];
writeOBJ(outputname, [medial_x,medial_y,medial_z],[],[],[],[],[])

fileID=fopen([output_path ,'Samples_',real_name2,'_',num2str(samplesnumber),'to',num2str(fval),'_scale',num2str(scale_setting),'.txt'],'w')
fprintf(fileID,'%f %f %f\n',N');
fclose(fileID);

GG = [medial_x,medial_y,medial_z,init_r];
GG2 = [medial_x,medial_y,medial_z,medial_r];
% save the ma with radius info.

filename = [output_path,'MA_init',real_name2,'_scale',num2str(scale_setting),'.txt'];

fileID=fopen(filename,'w')
fprintf(fileID,'v %f %f %f %f\n',GG(find(x),:)');
fclose(fileID);
  
filename = [output_path,'vis_MA_init',real_name2,'_scale',num2str(scale_setting),'.obj'];

fileID=fopen(filename,'w')
fprintf(fileID,'v %f %f %f\n',GG(find(x),1:3)');
fclose(fileID);

filename = [output_path,'MA_modified',real_name2,'_scale',num2str(scale_setting),'.txt'];

fileID=fopen(filename,'w')
fprintf(fileID,'v %f %f %f %f\n',GG2(find(x),:)');
fclose(fileID);

filename = [output_path,'MA_input_enlarged',real_name2,'_scale',num2str(scale_setting),'.txt'];

fileID=fopen(filename,'w')
fprintf(fileID,'v %f %f %f %f\n',[medial_x,medial_y,medial_z,medial_r]');
fclose(fileID);


filename = [output_path,'MA_input_initial',real_name2,'_scale',num2str(scale_setting),'.txt'];

fileID=fopen(filename,'w')
fprintf(fileID,'v %f %f %f %f\n',[medial_x,medial_y,medial_z,init_r]');
fclose(fileID);

disp('min_init_medial_r')
disp(min(init_r));
disp('solution info')
disp(sum(x))
disp(["radius:",r])
