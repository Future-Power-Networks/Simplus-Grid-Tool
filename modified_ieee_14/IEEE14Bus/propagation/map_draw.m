figure(1);
clf;
image((imread('bus14_orig.jpg')));

set(gca,'YDir','reverse');
axis off;
xgain=450/270;
ygain=361/216;

mp=cell(14);
mp{1,1}=[51,363];
mp{2,2}=[262,575];
mp{3,3}=[643,625];
mp{4,4}=[643,381];
mp{5,5}=[427,440];
mp{6,6}=[428,324];
mp{7,7}=[668,318];
mp{8,8}=[746,320];
mp{9,9}=[652,261];
mp{10,10}=[564,222];
mp{11,11}=[503,209];
mp{12,12}=[270,216];
mp{13,13}=[432,143];
mp{14,14}=[568,143];
mp{1,2} = [123,470];
mp{2,5} = [366,510];
mp{1,5} = [236,430];
mp{9,14} = [614,220];
mp{12,13} = [365,197];
mp{13,14} = [506,163];
mp{4,9} = [618,315];
mp{2,4} = [466,482];
mp{6,13} = [421,236];
mp{6,12} = [344,283];
mp{7,9} = [672,287];
mp{4,7} = [672,363];
mp{3,4} = [653,496];
mp{10,11} = [531,229];
mp{9,10} = [598,246];


size_gain=3.7;
mmode=1;
for i=1:14
    for j=i:14
        if MdSensResult(mmode).Layer2_mat(i,j)==0
        else
            radius=MdSensResult(mmode).Layer1_mat(i,j)*size_gain; % radius decided by layer-1
            if real(MdSensResult(mmode).Layer2_mat(i,j))>0 % positive real part: red
                mcolor=[0.7 0 0];
            elseif real(MdSensResult(mmode).Layer2_mat(i,j))<0 % negative real part: blue
                mcolor=[0 0.3 0.7];
            end
            if i==j
                mx=mp{i,i}(1)*xgain;
                my=mp{i,i}(2)*ygain;
                mpos=[mx-radius/2,my-radius/2,radius,radius];
            else
                if isempty(mp{i,j})
                    mx=(mp{i,i}(1)+mp{j,j}(1))/2*xgain;
                    my=(mp{i,i}(2)+mp{j,j}(2))/2*ygain;
                else
                    mx=mp{i,j}(1)*xgain;
                    my=mp{i,j}(2)*ygain;
                end
                mpos=[mx-radius/2,my-radius/2,radius,radius];
                %rectangle('position',mpos,'Curvature',[1 1],'edgecolor',[mcolor,0],'facecolor',[mcolor,0.2]);
            end
                rectangle('position',mpos,'Curvature',[1 1],'edgecolor',[mcolor,0],'facecolor',[mcolor,0.5]);

        end
    end
end
set(gca,'position',[0 0 1 1])
img=gcf;%getimage(gcf);
print(img,'-dpng','-r300','./IEEE14Bus/propagation/bus14_prop.png')


